"""Structural variant detection with the Manta caller from Illumina.

https://github.com/Illumina/manta
"""
import collections
import glob
import os

from bcbio import utils
from bcbio.bam import ref
from bcbio.distributed.transaction import file_transaction
from bcbio.heterogeneity import chromhacks
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.variation import vcfutils
from bcbio.provenance import do, programs
from bcbio.structural import shared

def run(items):
    """Perform detection of structural variations with Manta.
    """
    paired = vcfutils.get_paired(items)
    data = paired.tumor_data if paired else items[0]
    work_dir = _sv_workdir(data)
    variant_file = _get_out_file(work_dir, paired)
    if not utils.file_exists(variant_file):
        with file_transaction(data, work_dir) as tx_work_dir:
            utils.safe_makedir(tx_work_dir)
            tx_workflow_file = _prep_config(items, paired, tx_work_dir)
            _run_workflow(items, paired, tx_workflow_file, tx_work_dir)
    assert utils.file_exists(variant_file), "Manta finished without output file %s" % variant_file
    variant_file = shared.annotate_with_depth(variant_file, items)
    out = []
    upload_counts = collections.defaultdict(int)
    for data in items:
        if "break-point-inspector" in dd.get_tools_on(data):
            if paired and paired.normal_bam and paired.tumor_name == dd.get_sample_name(data):
                variant_file = _run_break_point_inspector(data, variant_file, paired, work_dir)
        if "sv" not in data:
            data["sv"] = []
        final_vcf = shared.finalize_sv(variant_file, data, items)
        vc = {"variantcaller": "manta",
              "do_upload": upload_counts[final_vcf] == 0,  # only upload a single file per batch
              "vrn_file": final_vcf}
        evidence_bam = _get_evidence_bam(work_dir, data)
        if evidence_bam:
            vc["read_evidence"] = evidence_bam
        data["sv"].append(vc)
        upload_counts[final_vcf] += 1
        out.append(data)
    return out

def _run_break_point_inspector(data, variant_file, paired, work_dir):
    output_vcf = "%s-%s.vcf.gz" % (utils.splitext_plus(variant_file)[0], "bpi")
    stats_file = "%s-%s_stats.txt" % (utils.splitext_plus(variant_file)[0], "bpi")
    if not utils.file_exists(output_vcf):
        with file_transaction(data, output_vcf) as tx_output_vcf, file_transaction(data, stats_file) as tx_stats_file:
            cores = dd.get_num_cores(data)
            resources = config_utils.get_resources("break-point-inspector", data["config"])
            jvm_mem_opts = config_utils.adjust_opts(resources.get("jvm_opts", ["-Xms1000m", "-Xmx2000m"]),
                                                    {"algorithm": {"memory_adjust": {"magnitude": cores,
                                                                                     "direction": "increase"}}})
            jvm_tmp_arg = "-Djava.io.tmpdir=" + utils.safe_makedir(os.path.join(work_dir, "bpi_tmp"))
            cmd = ["break-point-inspector"] + jvm_mem_opts + [jvm_tmp_arg, "-vcf", variant_file]
            if paired:
                cmd += ["-ref", paired.normal_bam, "-tumor", paired.tumor_bam]
            cmd += ["-output_vcf", tx_output_vcf, ">", tx_stats_file]
            do.run(cmd, "Running Break Point Inspector for Manta SV calls")
    return output_vcf

def _get_out_file(work_dir, paired):
    """Retrieve manta output variant file, depending on analysis.
    """
    if paired:
        if paired.normal_bam:
            base_file = "somaticSV.vcf.gz"
        else:
            base_file = "tumorSV.vcf.gz"
    else:
        base_file = "diploidSV.vcf.gz"
    return os.path.join(work_dir, "results", "variants", base_file)

def _get_evidence_bam(work_dir, data):
    """Retrieve evidence BAM for the sample if it exists
    """
    evidence_bam = glob.glob(os.path.join(work_dir, "results", "evidence",
                                            "evidence_*.%s*.bam" % (dd.get_sample_name(data))))
    if evidence_bam:
        return evidence_bam[0]

def _run_workflow(items, paired, workflow_file, work_dir):
    """Run manta analysis inside prepared workflow directory.
    """
    utils.remove_safe(os.path.join(work_dir, "workspace"))
    data = paired.tumor_data if paired else items[0]
    cmd = [utils.get_program_python("configManta.py"), workflow_file, "-m", "local", "-j", dd.get_num_cores(data)]
    do.run(cmd, "Run manta SV analysis")
    utils.remove_safe(os.path.join(work_dir, "workspace"))

def _prep_config(items, paired, work_dir):
    """Run initial configuration, generating a run directory for Manta.
    """
    assert utils.which("configManta.py"), "Could not find installed configManta.py"
    out_file = os.path.join(work_dir, "runWorkflow.py")
    if not utils.file_exists(out_file) or _out_of_date(out_file):
        config_script = os.path.realpath(utils.which("configManta.py"))
        cmd = [utils.get_program_python("configManta.py"), config_script]
        if paired:
            if paired.normal_bam:
                cmd += ["--normalBam=%s" % paired.normal_bam, "--tumorBam=%s" % paired.tumor_bam]
            else:
                cmd += ["--tumorBam=%s" % paired.tumor_bam]
        else:
            cmd += ["--bam=%s" % dd.get_align_bam(data) for data in items]
        data = paired.tumor_data if paired else items[0]
        cmd += ["--referenceFasta=%s" % dd.get_ref_file(data), "--runDir=%s" % work_dir]
        if dd.get_coverage_interval(data) not in ["genome"]:
            cmd += ["--exome"]
        for region in _maybe_limit_chromosomes(data):
            cmd += ["--region", region]
        resources = config_utils.get_resources("manta", data["config"])
        if resources.get("options"):
            cmd += [str(x) for x in resources["options"]]
        # If we are removing polyX, avoid calling on small indels which require
        # excessively long runtimes on noisy WGS runs
        if "polyx" in dd.get_exclude_regions(data):
            cmd += ["--config", _prep_streamlined_config(config_script, work_dir)]
        do.run(cmd, "Configure manta SV analysis")
    return out_file

def _prep_streamlined_config(config_script, work_dir):
    """Create manta INI file without steps that potentially increase runtimes.

    This removes calling of small indels.
    """
    new_min_size = 100
    in_file = config_script + ".ini"
    out_file = os.path.join(work_dir, os.path.basename(in_file))
    with open(in_file) as in_handle:
        with open(out_file, "w") as out_handle:
            for line in in_handle:
                if line.startswith("minCandidateVariantSize"):
                    out_handle.write("minCandidateVariantSize = %s\n" % new_min_size)
                else:
                    out_handle.write(line)
    return out_file

def _maybe_limit_chromosomes(data):
    """Potentially limit chromosomes to avoid problematically named HLA contigs.

    HLAs have ':' characters in them which confuse downstream processing. If
    we have no problematic chromosomes we don't limit anything.
    """
    std_chroms = []
    prob_chroms = []
    noalt_calling = "noalt_calling" in dd.get_tools_on(data) or "altcontigs" in dd.get_exclude_regions(data)
    for contig in ref.file_contigs(dd.get_ref_file(data)):
        if contig.name.find(":") > 0 or (noalt_calling and not chromhacks.is_nonalt(contig.name)):
            prob_chroms.append(contig.name)
        else:
            std_chroms.append(contig.name)
    if len(prob_chroms) > 0:
        return std_chroms
    else:
        return []

def _sv_workdir(data):
    return os.path.join(
        data["dirs"]["work"], "structural", dd.get_sample_name(data), "manta")

def _out_of_date(rw_file):
    """Check if a run workflow file points to an older version of manta and needs a refresh.
    """
    with open(rw_file) as in_handle:
        for line in in_handle:
            if line.startswith("sys.path.append"):
                file_version = line.split("/lib/python")[0].split("Cellar/manta/")[-1]
                if file_version != programs.get_version_manifest("manta"):
                    return True
    return False
