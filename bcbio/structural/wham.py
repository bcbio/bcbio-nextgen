"""Identify structural variants using association testing with WHAM.

https://github.com/jewmanchue/wham
"""
import collections
import os

import vcf

from bcbio import utils
from bcbio.bam import ref
from bcbio.distributed.transaction import file_transaction
from bcbio.heterogeneity import chromhacks
from bcbio.pipeline import datadict as dd
from bcbio.structural import shared
from bcbio.variation import vcfutils
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
        final_vcf = shared.finalize_sv(orig_vcf, data, items)
        data["sv"].append({"variantcaller": "wham", "vrn_file": final_vcf})
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
            cores = dd.get_cores(inputs[0])
            ref_file = dd.get_ref_file(inputs[0])
            include_chroms = ",".join([c.name for c in ref.file_contigs(ref_file)
                                       if chromhacks.is_autosomal_or_x(c.name)])
            all_bams = ",".join([x["align_bam"] for x in inputs] + background_bams)
            cmd = ("whamg -x {cores} -a {ref_file} -f {all_bams} -c {include_chroms} "
                   "| bgzip -c > {tx_out_file}")
            do.run(cmd.format(**locals()), "WHAM SV caller: %s" % ", ".join(dd.get_sample_name(d) for d in inputs))
    return vcfutils.bgzip_and_index(out_file, inputs[0]["config"])

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
