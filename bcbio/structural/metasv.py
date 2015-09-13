"""Perform ensemble calling of structural variants using MetaSV.

https://github.com/chapmanb/metasv
http://dx.doi.org/10.1093/bioinformatics/btv204
"""
import os
import sys

from bcbio import utils
from bcbio.provenance import do
from bcbio.pipeline import datadict as dd
from bcbio.structural import shared
from bcbio.variation import vfilter

MIN_CALLERS = 2
SUPPORTED = set(["manta", "lumpy", "cnvkit", "wham"])

def run(calls, data):
    """Run MetaSV if we have enough supported callers, adding output to the set of calls.
    """
    work_dir = _sv_workdir(data)
    out_file = os.path.join(work_dir, "variants.vcf.gz")
    cmd = _get_cmd() + ["--sample", dd.get_sample_name(data), "--reference", dd.get_ref_file(data),
                        "--bam", dd.get_align_bam(data), "--outdir", work_dir]
    available_callers = 0
    for call in calls:
        if call["variantcaller"] in SUPPORTED:
            available_callers += 1
            cmd += ["--%s_vcf" % call["variantcaller"], call.get("vcf_file", call["vrn_file"])]
    if available_callers >= MIN_CALLERS:
        if not utils.file_exists(out_file):
            tx_work_dir = utils.safe_makedir(os.path.join(work_dir, "raw"))
            ins_stats = shared.calc_paired_insert_stats_save(dd.get_align_bam(data),
                                                            os.path.join(tx_work_dir, "insert-stats.yaml"))
            cmd += ["--workdir", tx_work_dir, "--num_threads", str(dd.get_num_cores(data))]
            cmd += ["--spades", utils.which("spades.py"), "--age", utils.which("age_align")]
            cmd += ["--assembly_max_tools=2", "--assembly_pad=500"]
            cmd += ["--boost_ins", "--isize_mean", ins_stats["mean"], "--isize_sd", ins_stats["std"]]
            do.run(cmd, "Combine variant calls with MetaSV")
        filters = ("(NUM_SVTOOLS = 1 && BA_NUM_GOOD_REC=0) || "
                   "(NUM_SVTOOLS = 1 && ABS(SVLEN)>10000) || "
                   "(NUM_SVTOOLS = 1 && ABS(SVLEN)<5000 && BA_FLANK_PERCENT>40) || "
                   "(ABS(SVLEN)<5000 && BA_NUM_GOOD_REC>1)")
        filter_file = vfilter.hard_w_expression(out_file, filters,
                                                data, name="ReassemblyStats", limit_regions=None)
        calls.append({"variantcaller": "metasv",
                      "vrn_file": filter_file})
    return calls

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "metasv"))

def _get_cmd():
    return [sys.executable, os.path.join(os.path.dirname(sys.executable), "run_metasv.py")]
