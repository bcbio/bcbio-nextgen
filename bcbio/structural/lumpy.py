"""Structural variation detection for split and paired reads using lumpy.

Uses lumpyexpress for lumpy integration and samblaster for read preparation:
https://github.com/GregoryFaust/samblaster
https://github.com/arq5x/lumpy-sv
"""
import contextlib
import os
import sys
import shutil

import vcf

from bcbio import utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.structural import shared as sshared
from bcbio.variation import vcfutils, vfilter

# ## Lumpy main

def _run_lumpy(full_bams, sr_bams, disc_bams, work_dir, items):
    """Run lumpy-sv, using speedseq pipeline.
    """
    batch = sshared.get_cur_batch(items)
    ext = "-%s-svs" % batch if batch else "-svs"
    out_file = os.path.join(work_dir, "%s%s.vcf"
                            % (os.path.splitext(os.path.basename(items[0]["align_bam"]))[0], ext))
    sv_exclude_bed = sshared.prepare_exclude_file(items, out_file)
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            with tx_tmpdir(items[0]) as tmpdir:
                full_bams = ",".join(full_bams)
                sr_bams = ",".join(sr_bams)
                disc_bams = ",".join(disc_bams)
                exclude = "-x %s" % sv_exclude_bed if utils.file_exists(sv_exclude_bed) else ""
                ref_file = dd.get_ref_file(items[0])
                # use our bcbio python for runs within lumpyexpress
                curpython_dir = os.path.dirname(sys.executable)
                cmd = ("export PATH={curpython_dir}:$PATH && "
                       "lumpyexpress -v -B {full_bams} -S {sr_bams} -D {disc_bams} "
                       "{exclude} -T {tmpdir} -o {tx_out_file}")
                do.run(cmd.format(**locals()), "lumpyexpress", items[0])
    return vcfutils.sort_by_ref(out_file, items[0]), sv_exclude_bed

def _filter_by_support(in_file, data):
    """Filter call file based on supporting evidence, adding FILTER annotations to VCF.

    Filters based on the following criteria:
      - Minimum read support for the call (SU = total support)
      - Large calls need split read evidence.
    """
    rc_filter = ("FORMAT/SU < 4 || "
                 "(FORMAT/SR == 0 && ABS(SVLEN)>50000)")
    return vfilter.hard_w_expression(in_file, rc_filter, data, name="ReadCountSupport",
                                     limit_regions=None)

def _filter_by_background(base_samples, back_samples, gt_vcfs, data):
    """Filter base samples, marking any also present in the background.
    """
    filtname = "InBackground"
    filtdoc = "Variant also present in background samples with same genotype"
    for base_name in base_samples:
        orig_vcf = gt_vcfs[base_name]
        out_file = "%s-backfilter.vcf" % (utils.splitext_plus(orig_vcf)[0])
        if not utils.file_exists(out_file) and not utils.file_exists(out_file + ".gz"):
            with file_transaction(data, out_file) as tx_out_file:
                with utils.open_gzipsafe(orig_vcf) as in_handle:
                    with _vcf_readers([gt_vcfs[n] for n in back_samples]) as back_readers:
                        inp = vcf.Reader(in_handle, orig_vcf)
                        inp.filters[filtname] = vcf.parser._Filter(filtname, filtdoc)
                        with open(tx_out_file, "w") as out_handle:
                            outp = vcf.Writer(out_handle, inp)
                            for rec in inp:
                                back_recs = [r.next() for r in back_readers]
                                if _genotype_in_background(rec, back_recs):
                                    rec.add_filter(filtname)
                                outp.write_record(rec)
        gt_vcfs[base_name] = vcfutils.bgzip_and_index(out_file, data["config"])
    return gt_vcfs

def _genotype_in_background(rec, back_recs):
    """Check if the genotype in the record of interest is present in the background records.
    """
    def passes(rec):
        return not rec.FILTER or len(rec.FILTER) == 0
    return any([passes(brec) and passes(rec) and rec.samples[0].gt_alleles == brec.samples[0].gt_alleles
                for brec in back_recs])

@contextlib.contextmanager
def _vcf_readers(vcf_files):
    handles = []
    readers = []
    for vcf_file in vcf_files:
        in_handle = utils.open_gzipsafe(vcf_file)
        handles.append(in_handle)
        readers.append(vcf.Reader(in_handle, vcf_file))
    yield readers
    for handle in handles:
        handle.close()

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "lumpy"))

def run(items):
    """Perform detection of structural variations with lumpy, using bwa-mem alignment.
    """
    if not all(utils.get_in(data, ("config", "algorithm", "aligner")) in ["bwa", False, None] for data in items):
        raise ValueError("Require bwa-mem alignment input for lumpy structural variation detection")
    paired = vcfutils.get_paired_bams([x["align_bam"] for x in items], items)
    work_dir = _sv_workdir(paired.tumor_data if paired and paired.tumor_data else items[0])
    full_bams, sr_bams, disc_bams = [], [], []
    for data in items:
        dedup_bam, sr_bam, disc_bam = sshared.get_split_discordants(data, work_dir)
        full_bams.append(dedup_bam)
        sr_bams.append(sr_bam)
        disc_bams.append(disc_bam)
    lumpy_vcf, exclude_file = _run_lumpy(full_bams, sr_bams, disc_bams, work_dir, items)
    gt_vcfs = {}
    for data in items:
        sample = dd.get_sample_name(data)
        dedup_bam, sr_bam, _ = sshared.get_split_discordants(data, work_dir)
        sample_vcf = vcfutils.select_sample(lumpy_vcf, sample,
                                            utils.append_stem(lumpy_vcf, "-%s" % sample),
                                            data["config"])
        gt_vcf = _run_svtyper(sample_vcf, dedup_bam, sr_bam, data)
        gt_vcfs[dd.get_sample_name(data)] = _filter_by_support(gt_vcf, data)
    if paired and paired.normal_name:
        gt_vcfs = _filter_by_background([paired.tumor_name], [paired.normal_name], gt_vcfs, paired.tumor_data)
    out = []
    for data in items:
        if "sv" not in data:
            data["sv"] = []
        data["sv"].append({"variantcaller": "lumpy",
                           "vrn_file": gt_vcfs[dd.get_sample_name(data)],
                           "exclude_file": exclude_file})
        out.append(data)
    return out

def _run_svtyper(in_file, full_bam, sr_bam, data):
    """Genotype structural variant calls with SVtyper.
    """
    out_file = "%s-wgts.vcf.gz" % utils.splitext_plus(in_file)[0]
    if not utils.file_uptodate(out_file, in_file):
        with file_transaction(data, out_file) as tx_out_file:
            if not vcfutils.vcf_has_variants(in_file):
                shutil.copy(in_file, out_file)
            else:
                python = sys.executable
                svtyper = os.path.join(os.path.dirname(sys.executable), "svtyper")
                cmd = ("gunzip -c {in_file} | "
                       "{python} {svtyper} -B {full_bam} -S {sr_bam} | "
                       "bgzip -c > {tx_out_file}")
                do.run(cmd.format(**locals()), "SV genotyping with svtyper")
    return vcfutils.sort_by_ref(out_file, data)
