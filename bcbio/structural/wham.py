"""Identify structural variants using association testing with WHAM.

https://github.com/jewmanchue/wham
"""
import collections
import os

import vcf

from bcbio import utils
from bcbio.distributed.multi import run_multicore, zeromq_aware_logging
from bcbio.distributed.transaction import file_transaction
from bcbio.heterogeneity import chromhacks
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import region
from bcbio.structural import shared
from bcbio.variation import bamprep, effects, vcfutils
from bcbio.provenance import do

def run(items, background=None):
    """Detect copy number variations from batched set of samples using WHAM.
    """
    if not background: background = []
    background_bams = []
    paired = vcfutils.get_paired_bams([x["align_bam"] for x in items], items)
    if paired:
        inputs = [paired.tumor_data]
        if paired.normal_bam:
            background = [paired.normal_data]
            background_bams = [paired.normal_bam]
    else:
        assert not background
        inputs, background = shared.find_case_control(items)
        background_bams = [x["align_bam"] for x in background]
    orig_vcf = _run_wham(inputs, background_bams)
    out = []
    for data in inputs:
        if "sv" not in data:
            data["sv"] = []
        sample_vcf = "%s-%s.vcf.gz" % (utils.splitext_plus(orig_vcf)[0], dd.get_sample_name(data))
        sample_vcf = vcfutils.select_sample(orig_vcf, dd.get_sample_name(data), sample_vcf, data["config"])
        if background:
            sample_vcf = filter_by_background(sample_vcf, orig_vcf, background, data)
        effects_vcf, _ = effects.add_to_vcf(sample_vcf, data, "snpeff")
        data["sv"].append({"variantcaller": "wham",
                           "vrn_file": effects_vcf or sample_vcf})
        out.append(data)
    return out

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "wham"))

def _run_wham(inputs, background_bams):
    """Run WHAM on a defined set of inputs and targets.
    """
    out_file = os.path.join(_sv_workdir(inputs[0]), "%s-wham.vcf.gz" % dd.get_sample_name(inputs[0]))
    if not utils.file_exists(out_file):
        with file_transaction(inputs[0], out_file) as tx_out_file:
            coords = chromhacks.autosomal_or_x_coords(dd.get_ref_file(inputs[0]))
            parallel = {"type": "local", "cores": dd.get_cores(inputs[0]), "progs": []}
            rs = run_multicore(_run_wham_coords,
                                [(inputs, background_bams, coord, out_file)
                                 for coord in coords],
                                inputs[0]["config"], parallel)
            rs = {coord: fname for (coord, fname) in rs}
            vcfutils.concat_variant_files([rs[c] for c in coords], tx_out_file, coords,
                                          dd.get_ref_file(inputs[0]), inputs[0]["config"])
    return out_file

@utils.map_wrap
@zeromq_aware_logging
def _run_wham_coords(inputs, background_bams, coords, final_file):
    """Run WHAM on a specific set of chromosome, start, end coordinates.
    """
    base, ext = utils.splitext_plus(final_file)
    raw_file = "%s-%s.vcf" % (base, region.to_safestr(coords))
    all_bams = ",".join([x["align_bam"] for x in inputs] + background_bams)
    if not utils.file_exists(raw_file):
        with file_transaction(inputs[0], raw_file) as tx_raw_file:
            cores = dd.get_cores(inputs[0])
            ref_file = dd.get_ref_file(inputs[0])
            coord_str = bamprep.region_to_gatk(coords)
            opts = "-k -m 30"
            cmd = ("WHAM-GRAPHENING {opts} -x {cores} -a {ref_file} -f {all_bams} -r {coord_str} "
                   "> {tx_raw_file}")
            do.run(cmd.format(**locals()), "Run WHAM: %s" % region.to_safestr(coords))
    merge_vcf = _run_wham_merge(raw_file, inputs[0])
    gt_vcf = _run_wham_genotype(merge_vcf, all_bams, coords, inputs[0])
    prep_vcf = vcfutils.sort_by_ref(gt_vcf, inputs[0])
    return [[coords, prep_vcf]]

def _run_wham_genotype(in_file, all_bams, coords, data):
    """Run genotyping on a prepped, merged VCF file.
    """
    out_file = "%s-wgts%s" % utils.splitext_plus(in_file)
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cores = dd.get_cores(data)
            ref_file = dd.get_ref_file(data)
            coord_str = bamprep.region_to_gatk(coords)
            cmd = ("WHAM-GRAPHENING -b {in_file} -x {cores} -a {ref_file} -f {all_bams} -r {coord_str} "
                   "> {tx_out_file}")
            do.run(cmd.format(**locals()), "Genotype WHAM: %s" % region.to_safestr(coords))
    return out_file

def _run_wham_merge(in_file, data):
    """Run WHAM merge functionality to combine closely spaced events.
    """
    out_file = "%s-merge%s" % utils.splitext_plus(in_file)
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = "mergeIndvs -f {in_file} > {tx_out_file}"
            do.run(cmd.format(**locals()), "Merge WHAM output")
    return out_file

def filter_by_background(in_vcf, full_vcf, background, data):
    """Filter SV calls also present in background samples.

    Skips filtering of inversions, which are not characterized differently
    between cases and controls in test datasets.
    """
    Filter = collections.namedtuple('Filter', ['id', 'desc'])
    back_filter = Filter(id='InBackground',
                         desc='Rejected due to presence in background sample')
    out_file = "%s-filter.vcf" % utils.splitext_plus(in_vcf)[0]
    if not utils.file_uptodate(out_file, in_vcf) and not utils.file_uptodate(out_file + ".vcf.gz", in_vcf):
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                reader = vcf.VCFReader(filename=in_vcf)
                reader.filters["InBackground"] = back_filter
                full_reader = vcf.VCFReader(filename=full_vcf)
                writer = vcf.VCFWriter(out_handle, template=reader)
                for out_rec, rec in zip(reader, full_reader):
                    rec_type = rec.genotype(dd.get_sample_name(data)).gt_type
                    if rec_type == 0 or any(rec_type == rec.genotype(dd.get_sample_name(x)).gt_type
                                            for x in background):
                        out_rec.add_filter("InBackground")
                    writer.write_record(out_rec)
    return vcfutils.bgzip_and_index(out_file, data["config"])
