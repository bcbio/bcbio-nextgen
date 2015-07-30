"""Structural variant detection with the Manta caller from Illumina.

https://github.com/Illumina/manta
"""
import os
import sys

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.variation import vcfutils
from bcbio.provenance import do

def run(items):
    """Perform detection of structural variations with Manta.
    """
    paired = vcfutils.get_paired(items)
    work_dir = _sv_workdir(paired.tumor_data if paired else items[0])
    workflow_file = _prep_config(items, paired, work_dir)
    variant_file = _run_workflow(items, paired, workflow_file, work_dir)
    sample_file = _select_sample(items, paired, variant_file, work_dir)
    out = []
    for data in items:
        if "sv" not in data:
            data["sv"] = []
        data["sv"].append({"variantcaller": "manta",
                           "vrn_file": sample_file})
        out.append(data)
    return out

def _run_workflow(items, paired, workflow_file, work_dir):
    """Run manta analysis inside prepared workflow directory.
    """
    data = paired.tumor_data if paired else items[0]
    out_file = os.path.join(work_dir, "results", "variants",
                            "somaticSV.vcf.gz" if paired and paired.normal_bam else "diploidSV.vcf.gz")
    if not utils.file_exists(out_file):
        cmd = [sys.executable, workflow_file, "-m", "local", "-j", dd.get_num_cores(data),
               "--quiet"]
        do.run(cmd, "Run manta SV analysis")
    return out_file

def _select_sample(items, paired, variant_file, work_dir):
    """Fix VCF to have the correct sample name and select tumor samples from paired analyses.
    """
    sample_name = paired.tumor_name if paired else dd.get_sample_name(items[0])
    out_file = os.path.join(work_dir, "%s-%s.vcf.gz" % (utils.splitext_plus(os.path.basename(variant_file))[0],
                                                        sample_name))
    if not utils.file_uptodate(out_file, variant_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            cmd = "zcat {variant_file} | "
            if paired and paired.normal_bam:
                cmd += "sed 's/\tTUMOR/\t{sample_name}/' | bcftools view -s {sample_name}"
            else:
                cmd += "sed 's/\tSAMPLE/\t{sample_name}/' | bcftools view -s {sample_name}"
            cmd += " | bgzip -c > {tx_out_file}"
            do.run(cmd.format(**locals()), "Run manta SV analysis")
    return vcfutils.bgzip_and_index(out_file, items[0]["config"])

def _prep_config(items, paired, work_dir):
    """Run initial configuration, generating a run directory for Manta.
    """
    assert utils.which("configManta.py"), "Could not find installed configManta.py"
    out_file = os.path.join(work_dir, "runWorkflow.py")
    if not utils.file_exists(out_file):
        cmd = [sys.executable, utils.which("configManta.py")]
        if paired:
            if paired.normal_bam:
                cmd += ["--normalBam=%s" % paired.normal_bam, "--tumorBam=%s" % paired.tumor_bam]
            else:
                cmd += ["--normalBam=%s" % paired.tumor_bam]
        else:
            assert len(items) == 1, "Expect a single item if non-paired for manta: %s" % \
                ([dd.get_sample_name(d) for d in items])
            cmd += ["--normalBam=%s" % dd.get_align_bam(items[0])]
        data = paired.tumor_data if paired else items[0]
        cmd += ["--referenceFasta=%s" % dd.get_ref_file(data), "--runDir=%s" % work_dir]
        if dd.get_coverage_interval(data) not in ["genome"]:
            cmd += ["--exome"]
        do.run(cmd, "Configure manta SV analysis")
    return out_file

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "manta"))
