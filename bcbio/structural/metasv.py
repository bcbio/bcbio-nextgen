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

def run(items):
    """Run MetaSV if we have enough supported callers, adding output to the set of calls.
    """
    assert len(items) == 1, "Expect one input to MetaSV ensemble calling"
    data = items[0]
    work_dir = _sv_workdir(data)
    out_file = os.path.join(work_dir, "variants.vcf.gz")
    cmd = _get_cmd() + ["--sample", dd.get_sample_name(data), "--reference", dd.get_ref_file(data),
                        "--bam", dd.get_align_bam(data), "--outdir", work_dir]
    methods = []
    for call in data.get("sv", []):
        if call["variantcaller"] in SUPPORTED and call["variantcaller"] not in methods:
            methods.append(call["variantcaller"])
            cmd += ["--%s_vcf" % call["variantcaller"], call.get("vcf_file", call["vrn_file"])]
    if len(methods) >= MIN_CALLERS:
        if not utils.file_exists(out_file):
            tx_work_dir = utils.safe_makedir(os.path.join(work_dir, "raw"))
            ins_stats = shared.calc_paired_insert_stats_save(dd.get_align_bam(data),
                                                             os.path.join(tx_work_dir, "insert-stats.yaml"))
            cmd += ["--workdir", tx_work_dir, "--num_threads", str(dd.get_num_cores(data))]
            cmd += ["--spades", utils.which("spades.py"), "--age", utils.which("age_align")]
            cmd += ["--assembly_max_tools=1", "--assembly_pad=500"]
            cmd += ["--boost_ins", "--isize_mean", ins_stats["mean"], "--isize_sd", ins_stats["std"]]
            do.run(cmd, "Combine variant calls with MetaSV")
        filters = ("(NUM_SVTOOLS = 1 && ABS(SVLEN)>50000) || "
                   "(NUM_SVTOOLS = 1 && ABS(SVLEN)<4000 && BA_FLANK_PERCENT>80) || "
                   "(NUM_SVTOOLS = 1 && ABS(SVLEN)<4000 && BA_NUM_GOOD_REC=0) || "
                   "(ABS(SVLEN)<4000 && BA_NUM_GOOD_REC>2)")
        filter_file = vfilter.hard_w_expression(out_file, filters,
                                                data, name="ReassemblyStats", limit_regions=None)
        data["sv"].append({"variantcaller": "metasv",
                           "vrn_file": filter_file})
    return [data]

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "metasv"))

def _get_cmd():
    return [sys.executable, os.path.join(os.path.dirname(sys.executable), "run_metasv.py")]
