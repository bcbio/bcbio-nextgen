"""Structural variation detection for split and paired reads using lumpy.

Uses lumpyexpress for lumpy integration and samblaster for read preparation:
https://github.com/GregoryFaust/samblaster
https://github.com/arq5x/lumpy-sv
"""
import os
import sys

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
    return vcfutils.bgzip_and_index(out_file, items[0]["config"]), sv_exclude_bed

def _filter_by_support(in_file, data):
    """Filter call file based on supporting evidence, adding FILTER annotations to VCF.

    Filters based on the following criteria:
      - Minimum read support for the call.
    Other filters not currently applied due to being too restrictive:
      - Multiple forms of evidence in any sample (split and paired end)
    """
    rc_filter = "FORMAT/SU < 4"
    # approach_filter = "FORMAT/SR == 0 || FORMAT/PE == 0"
    return vfilter.hard_w_expression(in_file, rc_filter, data, name="ReadCountSupport")

def run(items):
    """Perform detection of structural variations with lumpy, using bwa-mem alignment.
    """
    if not all(utils.get_in(data, ("config", "algorithm", "aligner")) in ["bwa", False, None] for data in items):
        raise ValueError("Require bwa-mem alignment input for lumpy structural variation detection")
    work_dir = utils.safe_makedir(os.path.join(items[0]["dirs"]["work"], "structural", items[0]["name"][-1],
                                               "lumpy"))
    full_bams, sr_bams, disc_bams = [], [], []
    for data in items:
        dedup_bam, sr_bam, disc_bam = sshared.get_split_discordants(data, work_dir)
        full_bams.append(dedup_bam)
        sr_bams.append(sr_bam)
        disc_bams.append(disc_bam)
    lumpy_vcf, exclude_file = _run_lumpy(full_bams, sr_bams, disc_bams, work_dir, items)
    out = []
    for i, data in enumerate(items):
        if "sv" not in data:
            data["sv"] = []
        sample = dd.get_sample_name(data)
        sample_vcf = _filter_by_support(vcfutils.select_sample(lumpy_vcf, sample,
                                                               utils.append_stem(lumpy_vcf, "-%s" % sample),
                                                               data["config"]),
                                        data)
        data["sv"].append({"variantcaller": "lumpy",
                           "vrn_file": sample_vcf,
                           "exclude_file": exclude_file})
        out.append(data)
    return out
