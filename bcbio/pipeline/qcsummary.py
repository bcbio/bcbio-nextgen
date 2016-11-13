"""Quality control and summary metrics for next-gen alignments and analysis.
"""
import collections
import copy
import csv
import os

import yaml
from datetime import datetime

import toolz as tz
from bcbio.bam import sambamba

from bcbio.variation.bedutils import clean_file
from bcbio import bam, utils
from bcbio.log import logger
from bcbio.pipeline import config_utils, run_info
from bcbio.provenance import do
import bcbio.pipeline.datadict as dd
from bcbio.variation import coverage as cov
from bcbio.rnaseq import gtf
from bcbio.variation.coverage import get_average_coverage
from bcbio.variation import bedutils


# ## High level functions to generate summary

def generate_parallel(samples, run_parallel):
    """Provide parallel preparation of summary information for alignment and variant calling.
    """
    to_analyze, extras = _split_samples_by_qc(samples)
    qced = run_parallel("pipeline_summary", to_analyze)
    samples = _combine_qc_samples(qced) + extras
    qsign_info = run_parallel("qsignature_summary", [samples])
    samples = run_parallel("multiqc_summary", [samples])
    summary_file = write_project_summary(samples, qsign_info)
    out = []
    for data in samples:
        if "summary" not in data[0]:
            data[0]["summary"] = {}
        data[0]["summary"]["project"] = summary_file
        if qsign_info:
            data[0]["summary"]["mixup_check"] = qsign_info[0]["out_dir"]
        out.append(data)
    out = _add_researcher_summary(out, summary_file)
    return out

def pipeline_summary(data):
    """Provide summary information on processing sample.
    """
    data = utils.to_single_data(data)
    work_bam = data.get("align_bam")
    if data["analysis"].lower().startswith("smallrna-seq"):
        work_bam = data["clean_fastq"]
    elif data["analysis"].lower().startswith("chip-seq"):
        work_bam = data["raw_bam"]
    elif not work_bam.endswith(".bam"):
        work_bam = None
    if dd.get_ref_file(data) is not None and work_bam:
        logger.info("QC: %s %s" % (dd.get_sample_name(data), ", ".join(dd.get_algorithm_qc(data))))
        data["summary"] = _run_qc_tools(work_bam, data)
    return [[data]]

def get_qc_tools(data):
    """Retrieve a list of QC tools to use based on configuration and analysis type.

    Uses defaults if previously set.
    """
    if dd.get_algorithm_qc(data):
        return dd.get_algorithm_qc(data)
    analysis = data["analysis"].lower()
    to_run = []
    if "fastqc" not in dd.get_tools_off(data):
        to_run.append("fastqc")
    if any([tool in dd.get_tools_on(data)
            for tool in ["qualimap", "qualimap_full"]]):
        to_run.append("qualimap")
    if analysis.startswith("rna-seq"):
        if gtf.is_qualimap_compatible(dd.get_gtf_file(data)):
            to_run.append("qualimap_rnaseq")
        else:
            logger.debug("GTF not compatible with Qualimap, skipping.")
    if analysis.startswith("smallrna-seq"):
        to_run.append("small-rna")
    if not analysis.startswith("smallrna-seq"):
        to_run.append("samtools")
        to_run.append("gemini")
        if tz.get_in(["config", "algorithm", "kraken"], data):
            to_run.append("kraken")
    if analysis.startswith(("standard", "variant", "variant2")):
        to_run += ["qsignature", "coverage", "variants", "picard"]
    if dd.get_umi_file(data):
        to_run += ["umi"]
    return to_run

def _run_qc_tools(bam_file, data):
    """Run a set of third party quality control tools, returning QC directory and metrics.

        :param bam_file: alignments in bam format
        :param data: dict with all configuration information

        :returns: dict with output of different tools
    """
    from bcbio.qc import fastqc, gemini, kraken, qsignature, qualimap, samtools, picard, srna, umi
    tools = {"fastqc": fastqc.run,
             "small-rna": srna.run,
             "samtools": samtools.run,
             "qualimap": qualimap.run,
             "qualimap_rnaseq": qualimap.run_rnaseq,
             "gemini": gemini.run,
             "qsignature": qsignature.run,
             "coverage": _run_coverage_qc,
             "variants": _run_variants_qc,
             "kraken": kraken.run,
             "picard": picard.run,
             "umi": umi.run}
    qc_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "qc", data["description"]))
    metrics = {}
    qc_out = {}
    for program_name in tz.get_in(["config", "algorithm", "qc"], data):
        qc_fn = tools[program_name]
        cur_qc_dir = os.path.join(qc_dir, program_name)
        out = qc_fn(bam_file, data, cur_qc_dir)
        qc_files = None
        if out and isinstance(out, dict):
            metrics.update(out)
        elif out and isinstance(out, basestring) and os.path.exists(out):
            qc_files = {"base": out, "secondary": []}
        if not qc_files:
            qc_files = _organize_qc_files(program_name, cur_qc_dir)
        if qc_files:
            qc_out[program_name] = qc_files

    bam.remove("%s-downsample%s" % os.path.splitext(bam_file))

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
            if os.path.isfile(fname):
                out_files.append(fname)
            elif os.path.isdir(fname) and not fname.endswith("tx"):
                for root, dirs, files in os.walk(fname):
                    out_files.extend([os.path.join(root, x) for x in files])
        if len(out_files) > 0:
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
        if not dd.get_align_bam(data) or not qcs:
            extras.append([data])
        else:
            for qc in qcs:
                add = copy.deepcopy(data)
                add["config"]["algorithm"]["qc"] = [qc]
                to_process.append([add])
    return to_process, extras

def _combine_qc_samples(samples):
    """Combine split QC analyses into single samples based on BAM files.
    """
    by_bam = collections.defaultdict(list)
    for data in [utils.to_single_data(x) for x in samples]:
        by_bam[dd.get_align_bam(data)].append(data)
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

def _other_pipeline_samples(summary_file, cur_samples):
    """Retrieve samples produced previously by another pipeline in the summary output.
    """
    cur_descriptions = set([s[0]["description"] for s in cur_samples])
    out = []
    if utils.file_exists(summary_file):
        with open(summary_file) as in_handle:
            for s in yaml.load(in_handle).get("samples", []):
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

# ## Coverage

def _run_coverage_qc(bam_file, data, out_dir):
    """Run coverage QC analysis"""
    out = dict()

    total_reads = sambamba.number_of_reads(data, bam_file)
    out['Total_reads'] = total_reads
    mapped = sambamba.number_of_mapped_reads(data, bam_file)
    out['Mapped_reads'] = mapped
    if total_reads:
        out['Mapped_reads_pct'] = 100.0 * mapped / total_reads
    if mapped:
        mapped_unique = sambamba.number_of_mapped_reads(data, bam_file, keep_dups=False)
        out['Mapped_unique_reads'] = mapped
        mapped_dups = mapped - mapped_unique
        out['Duplicates'] = mapped_dups
        out['Duplicates_pct'] = 100.0 * mapped_dups / mapped

        if dd.get_coverage(data):
            cov_bed_file = clean_file(dd.get_coverage(data), data, prefix="cov-", simple=True)
            merged_bed_file = bedutils.merge_overlaps(cov_bed_file, data)
            target_name = "coverage"
        else:
            merged_bed_file = dd.get_variant_regions_merged(data)
            target_name = "variant_regions"

        ontarget = sambamba.number_mapped_reads_on_target(
            data, merged_bed_file, bam_file, keep_dups=False, target_name=target_name)
        if mapped_unique:
            out["Ontarget_unique_reads"] = ontarget
            out["Ontarget_pct"] = 100.0 * ontarget / mapped_unique
            out['Offtarget_pct'] = 100.0 * (mapped_unique - ontarget) / mapped_unique
            padded_bed_file = bedutils.get_padded_bed_file(merged_bed_file, 200, data)
            ontarget_padded = sambamba.number_mapped_reads_on_target(
                data, padded_bed_file, bam_file, keep_dups=False, target_name=target_name + "_padded")
            out["Ontarget_padded_pct"] = 100.0 * ontarget_padded / mapped_unique
        if total_reads:
            out['Usable_pct'] = 100.0 * ontarget / total_reads

        avg_coverage = get_average_coverage(data, bam_file, merged_bed_file, target_name)
        out['Avg_coverage'] = avg_coverage

    priority = cov.priority_coverage(data, out_dir)
    cov.priority_total_coverage(data, out_dir)
    region_coverage_file = cov.coverage_region_detailed_stats(data, out_dir)
    # Re-enable with annotations from internally installed
    # problem region directory
    # if priority:
    #    annotated = cov.decorate_problem_regions(
    #        priority, problem_regions, data)

    return out

def _run_variants_qc(bam_file, data, out_dir):
    """Run variants QC analysis"""
    cov.variants(data, out_dir)
    return None

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
