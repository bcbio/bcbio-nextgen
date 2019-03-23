"""Quality control and summary metrics for next-gen alignments and analysis.
"""
import collections
import copy
import csv
import os

import yaml
from datetime import datetime
import pandas as pd
import glob

import toolz as tz

from bcbio import utils
from bcbio.cwl import cwlutils
from bcbio.log import logger
from bcbio.pipeline import config_utils, run_info
import bcbio.pipeline.datadict as dd
from bcbio.provenance import do
from bcbio.rnaseq import gtf
from bcbio.variation import damage, peddy, vcfutils, vcfanno

import six


# ## High level functions to generate summary

def qc_to_rec(samples):
    """CWL: Convert a set of input samples into records for parallelization.
    """
    samples = [utils.to_single_data(x) for x in samples]
    samples = cwlutils.assign_complex_to_samples(samples)
    to_analyze, extras = _split_samples_by_qc(samples)
    recs = cwlutils.samples_to_records([utils.to_single_data(x) for x in to_analyze + extras])
    return [[x] for x in recs]

def generate_parallel(samples, run_parallel):
    """Provide parallel preparation of summary information for alignment and variant calling.
    """
    to_analyze, extras = _split_samples_by_qc(samples)
    qced = run_parallel("pipeline_summary", to_analyze)
    samples = _combine_qc_samples(qced) + extras
    qsign_info = run_parallel("qsignature_summary", [samples])
    metadata_file = _merge_metadata([samples])
    summary_file = write_project_summary(samples, qsign_info)
    out = []
    for data in samples:
        if "summary" not in data[0]:
            data[0]["summary"] = {}
        data[0]["summary"]["project"] = summary_file
        data[0]["summary"]["metadata"] = metadata_file
        if qsign_info:
            data[0]["summary"]["mixup_check"] = qsign_info[0]["out_dir"]
        out.append(data)
    out = _add_researcher_summary(out, summary_file)
    # MultiQC must be run after all file outputs are set:
    return [[utils.to_single_data(d)] for d in run_parallel("multiqc_summary", [out])]

def pipeline_summary(data):
    """Provide summary information on processing sample.

    Handles standard and CWL (single QC output) cases.
    """
    data = utils.to_single_data(data)
    work_bam = dd.get_align_bam(data) or dd.get_work_bam(data)
    if not work_bam or not work_bam.endswith(".bam"):
        work_bam = None
    if dd.get_ref_file(data):
        if work_bam or (tz.get_in(["config", "algorithm", "kraken"], data)):  # kraken doesn't need bam
            logger.info("QC: %s %s" % (dd.get_sample_name(data), ", ".join(dd.get_algorithm_qc(data))))
            work_data = cwlutils.unpack_tarballs(utils.deepish_copy(data), data)
            data["summary"] = _run_qc_tools(work_bam, work_data)
            if (len(dd.get_algorithm_qc(data)) == 1 and "output_cwl_keys" in data):
                data["summary"]["qc"] = data["summary"]["qc"].get(dd.get_algorithm_qc(data)[0])
    return [[data]]

def get_qc_tools(data):
    """Retrieve a list of QC tools to use based on configuration and analysis type.

    Uses defaults if previously set.
    """
    if dd.get_algorithm_qc(data):
        return dd.get_algorithm_qc(data)
    analysis = data["analysis"].lower()
    to_run = []
    if tz.get_in(["config", "algorithm", "kraken"], data):
        to_run.append("kraken")
    if "fastqc" not in dd.get_tools_off(data):
        to_run.append("fastqc")
    if any([tool in dd.get_tools_on(data)
            for tool in ["qualimap", "qualimap_full"]]):
        to_run.append("qualimap")
    if analysis.startswith("rna-seq") or analysis == "smallrna-seq":
        if "qualimap" not in dd.get_tools_off(data):
            if gtf.is_qualimap_compatible(dd.get_gtf_file(data)):
                to_run.append("qualimap_rnaseq")
            else:
                logger.debug("GTF not compatible with Qualimap, skipping.")
    if analysis.startswith("chip-seq"):
        to_run.append("chipqc")
    if analysis.startswith("smallrna-seq"):
        to_run.append("small-rna")
        to_run.append("atropos")
    if "coverage_qc" not in dd.get_tools_off(data):
        to_run.append("samtools")
    if dd.has_variantcalls(data):
        if "coverage_qc" not in dd.get_tools_off(data):
            to_run += ["coverage", "picard"]
        to_run += ["qsignature", "variants"]
        if vcfanno.is_human(data):
            to_run += ["contamination", "peddy"]
        if vcfutils.get_paired_phenotype(data):
            to_run += ["viral"]
        if damage.should_filter([data]):
            to_run += ["damage"]
    if dd.get_umi_consensus(data):
        to_run += ["umi"]
    if tz.get_in(["config", "algorithm", "preseq"], data):
        to_run.append("preseq")
    to_run = [tool for tool in to_run if tool not in dd.get_tools_off(data)]
    to_run.sort()
    return to_run

def _run_qc_tools(bam_file, data):
    """Run a set of third party quality control tools, returning QC directory and metrics.

        :param bam_file: alignments in bam format
        :param data: dict with all configuration information

        :returns: dict with output of different tools
    """
    from bcbio.qc import (atropos, contamination, coverage, damage, fastqc, kraken,
                          qsignature, qualimap, samtools, picard, srna, umi, variant,
                          viral, preseq, chipseq)
    tools = {"fastqc": fastqc.run,
             "atropos": atropos.run,
             "small-rna": srna.run,
             "samtools": samtools.run,
             "qualimap": qualimap.run,
             "qualimap_rnaseq": qualimap.run_rnaseq,
             "qsignature": qsignature.run,
             "contamination": contamination.run,
             "coverage": coverage.run,
             "damage": damage.run,
             "variants": variant.run,
             "peddy": peddy.run_qc,
             "kraken": kraken.run,
             "picard": picard.run,
             "umi": umi.run,
             "viral": viral.run,
             "preseq": preseq.run,
             "chipqc": chipseq.run
             }
    qc_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "qc", data["description"]))
    metrics = {}
    qc_out = utils.deepish_copy(dd.get_summary_qc(data))
    for program_name in dd.get_algorithm_qc(data):
        if not bam_file and program_name != "kraken":  # kraken doesn't need bam
            continue
        if dd.get_phenotype(data) == "germline" and program_name != "variants":
            continue
        qc_fn = tools[program_name]
        cur_qc_dir = os.path.join(qc_dir, program_name)
        out = qc_fn(bam_file, data, cur_qc_dir)
        qc_files = None
        if out and isinstance(out, dict):
            # Check for metrics output, two cases:
            # 1. output with {"metrics"} and files ("base")
            if "metrics" in out:
                metrics.update(out.pop("metrics"))
            # 2. a dictionary of metrics
            elif "base" not in out:
                metrics.update(out)
            # Check for files only output
            if "base" in out:
                qc_files = out
        elif out and isinstance(out, six.string_types) and os.path.exists(out):
            qc_files = {"base": out, "secondary": []}
        if not qc_files:
            qc_files = _organize_qc_files(program_name, cur_qc_dir)
        if qc_files:
            qc_out[program_name] = qc_files

    metrics["Name"] = dd.get_sample_name(data)
    metrics["Quality format"] = dd.get_quality_format(data).lower()
    return {"qc": qc_out, "metrics": metrics}

def _organize_qc_files(program, qc_dir):
    """Organize outputs from quality control runs into a base file and secondary outputs.

    Provides compatibility with CWL output. Returns None if no files created during processing.
    """
    base_files = {"fastqc": "fastqc_report.html",
                  "qualimap_rnaseq": "qualimapReport.html",
                  "qualimap": "qualimapReport.html"}
    if os.path.exists(qc_dir):
        out_files = []
        for fname in [os.path.join(qc_dir, x) for x in os.listdir(qc_dir)]:
            if os.path.isfile(fname) and not fname.endswith(".bcbiotmp"):
                out_files.append(fname)
            elif os.path.isdir(fname) and not fname.endswith("tx"):
                for root, dirs, files in os.walk(fname):
                    for f in files:
                        if not f.endswith(".bcbiotmp"):
                            out_files.append(os.path.join(root, f))
        if len(out_files) > 0 and all([not f.endswith("-failed.log") for f in out_files]):
            if len(out_files) == 1:
                base = out_files[0]
                secondary = []
            else:
                base = None
                if program in base_files:
                    base_choices = [x for x in out_files if x.endswith("/%s" % base_files[program])]
                    if len(base_choices) == 1:
                        base = base_choices[0]
                if not base:
                    base = out_files[0]
                secondary = [x for x in out_files if x != base]
            return {"base": base, "secondary": secondary}

# ## Allow parallelization for separate QC runs

def _split_samples_by_qc(samples):
    """Split data into individual quality control steps for a run.
    """
    to_process = []
    extras = []
    for data in [utils.to_single_data(x) for x in samples]:
        qcs = dd.get_algorithm_qc(data)
        # kraken doesn't need bam
        if qcs and (dd.get_align_bam(data) or dd.get_work_bam(data) or
                    tz.get_in(["config", "algorithm", "kraken"], data)):
            for qc in qcs:
                add = copy.deepcopy(data)
                add["config"]["algorithm"]["qc"] = [qc]
                to_process.append([add])
        else:
            extras.append([data])
    return to_process, extras

def _combine_qc_samples(samples):
    """Combine split QC analyses into single samples based on BAM files.
    """
    by_bam = collections.defaultdict(list)
    for data in [utils.to_single_data(x) for x in samples]:
        batch = dd.get_batch(data) or dd.get_sample_name(data)
        if not isinstance(batch, (list, tuple)):
            batch = [batch]
        batch = tuple(batch)
        by_bam[(dd.get_align_bam(data) or dd.get_work_bam(data), batch)].append(data)
    out = []
    for data_group in by_bam.values():
        data = data_group[0]
        alg_qc = []
        qc = {}
        metrics = {}
        for d in data_group:
            qc.update(dd.get_summary_qc(d))
            metrics.update(dd.get_summary_metrics(d))
            alg_qc.extend(dd.get_algorithm_qc(d))
        data["config"]["algorithm"]["qc"] = alg_qc
        data["summary"]["qc"] = qc
        data["summary"]["metrics"] = metrics
        out.append([data])
    return out

# ## Generate project level QC summary for quickly assessing large projects

def write_project_summary(samples, qsign_info=None):
    """Write project summary information on the provided samples.
    write out dirs, genome resources,

    """
    work_dir = samples[0][0]["dirs"]["work"]
    out_file = os.path.join(work_dir, "project-summary.yaml")
    upload_dir = (os.path.join(work_dir, samples[0][0]["upload"]["dir"])
                  if "dir" in samples[0][0]["upload"] else "")
    date = str(datetime.now())
    prev_samples = _other_pipeline_samples(out_file, samples)
    with open(out_file, "w") as out_handle:
        yaml.safe_dump({"date": date}, out_handle,
                       default_flow_style=False, allow_unicode=False)
        if qsign_info:
            qsign_out = utils.deepish_copy(qsign_info[0])
            qsign_out.pop("out_dir", None)
            yaml.safe_dump({"qsignature": qsign_out}, out_handle, default_flow_style=False,
                           allow_unicode=False)
        yaml.safe_dump({"upload": upload_dir}, out_handle,
                       default_flow_style=False, allow_unicode=False)
        yaml.safe_dump({"bcbio_system": samples[0][0]["config"].get("bcbio_system", "")}, out_handle,
                       default_flow_style=False, allow_unicode=False)
        yaml.safe_dump({"samples": prev_samples + [_save_fields(sample[0]) for sample in samples]}, out_handle,
                       default_flow_style=False, allow_unicode=False)
    return out_file

def _merge_metadata(samples):
    """Merge all metadata into CSV file"""
    samples = list(utils.flatten(samples))
    out_dir = dd.get_work_dir(samples[0])
    logger.info("summarize metadata")
    out_file = os.path.join(out_dir, "metadata.csv")
    sample_metrics = collections.defaultdict(dict)
    for s in samples:
        m = tz.get_in(['metadata'], s)
        if isinstance(m, six.string_types):
            m = json.loads(m)
        if m:
            for me in list(m.keys()):
                if isinstance(m[me], list) or isinstance(m[me], dict) or isinstance(m[me], tuple):
                    m.pop(me, None)
            sample_metrics[dd.get_sample_name(s)].update(m)
    pd.DataFrame(sample_metrics).transpose().to_csv(out_file)
    return out_file

def _other_pipeline_samples(summary_file, cur_samples):
    """Retrieve samples produced previously by another pipeline in the summary output.
    """
    cur_descriptions = set([s[0]["description"] for s in cur_samples])
    out = []
    if utils.file_exists(summary_file):
        with open(summary_file) as in_handle:
            for s in yaml.safe_load(in_handle).get("samples", []):
                if s["description"] not in cur_descriptions:
                    out.append(s)
    return out

def _save_fields(sample):
    to_save = ["dirs", "genome_resources", "genome_build", "sam_ref", "metadata",
               "description"]
    saved = {k: sample[k] for k in to_save if k in sample}
    if "summary" in sample and "metrics" in sample["summary"]:
        saved["summary"] = {"metrics": sample["summary"]["metrics"]}
    return saved

# ## Generate researcher specific summaries

def _add_researcher_summary(samples, summary_yaml):
    """Generate summary files per researcher if organized via a LIMS.
    """
    by_researcher = collections.defaultdict(list)
    for data in (x[0] for x in samples):
        researcher = utils.get_in(data, ("upload", "researcher"))
        if researcher:
            by_researcher[researcher].append(data["description"])
    out_by_researcher = {}
    for researcher, descrs in by_researcher.items():
        out_by_researcher[researcher] = _summary_csv_by_researcher(summary_yaml, researcher,
                                                                   set(descrs), samples[0][0])
    out = []
    for data in (x[0] for x in samples):
        researcher = utils.get_in(data, ("upload", "researcher"))
        if researcher:
            data["summary"]["researcher"] = out_by_researcher[researcher]
        out.append([data])
    return out

def _summary_csv_by_researcher(summary_yaml, researcher, descrs, data):
    """Generate a CSV file with summary information for a researcher on this project.
    """
    out_file = os.path.join(utils.safe_makedir(os.path.join(data["dirs"]["work"], "researcher")),
                            "%s-summary.tsv" % run_info.clean_name(researcher))
    metrics = ["Total_reads", "Mapped_reads", "Mapped_reads_pct", "Duplicates", "Duplicates_pct"]
    with open(summary_yaml) as in_handle:
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle, dialect="excel-tab")
            writer.writerow(["Name"] + metrics)
            for sample in yaml.safe_load(in_handle)["samples"]:
                if sample["description"] in descrs:
                    row = [sample["description"]] + [utils.get_in(sample, ("summary", "metrics", x), "")
                                                     for x in metrics]
                    writer.writerow(row)
    return out_file

# ## Galaxy functionality

def prep_pdf(qc_dir, config):
    """Create PDF from HTML summary outputs in QC directory.

    Requires wkhtmltopdf installed: http://www.msweet.org/projects.php?Z1
    Thanks to: https://www.biostars.org/p/16991/

    Works around issues with CSS conversion on CentOS by adjusting CSS.
    """
    html_file = os.path.join(qc_dir, "fastqc", "fastqc_report.html")
    html_fixed = "%s-fixed%s" % os.path.splitext(html_file)
    try:
        topdf = config_utils.get_program("wkhtmltopdf", config)
    except config_utils.CmdNotFound:
        topdf = None
    if topdf and utils.file_exists(html_file):
        out_file = "%s.pdf" % os.path.splitext(html_file)[0]
        if not utils.file_exists(out_file):
            cmd = ("sed 's/div.summary/div.summary-no/' %s | sed 's/div.main/div.main-no/' > %s"
                   % (html_file, html_fixed))
            do.run(cmd, "Fix fastqc CSS to be compatible with wkhtmltopdf")
            cmd = [topdf, html_fixed, out_file]
            do.run(cmd, "Convert QC HTML to PDF")
        return out_file
