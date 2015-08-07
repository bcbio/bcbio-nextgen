"""Perform ensemble calling of structural variants using MetaSV.

https://github.com/chapmanb/metasv
http://dx.doi.org/10.1093/bioinformatics/btv204
"""
import os
import sys

from bcbio import utils
from bcbio.provenance import do
from bcbio.pipeline import datadict as dd

MIN_CALLERS = 2
SUPPORTED = set(["manta", "lumpy", "cnvkit"])

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
            cmd += ["--spades", utils.which("spades.py"), "--age", utils.which("age_align")]
            do.run(cmd, "Combine variant calls with MetaSV")
        calls.append({"variantcaller": "metasv",
                      "vrn_file": out_file})
    return calls

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "metasv"))

def _get_cmd():
    return [sys.executable, os.path.join(os.path.dirname(sys.executable), "run_metasv.py")]
