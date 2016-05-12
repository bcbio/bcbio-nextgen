"""Quality control and summary metrics for next-gen alignments and analysis.  """
import collections
import contextlib
import csv
import os
import glob
import shutil
import subprocess

import pandas as pd
import lxml.html
import yaml
from datetime import datetime
from collections import defaultdict

try:
    from fadapa import Fadapa
except ImportError:
    Fadapa = None
import pybedtools
import pysam
import toolz as tz
import toolz.dicttoolz as dtz

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.log import logger
from bcbio.pipeline import config_utils, run_info
from bcbio.install import _get_data_dir
from bcbio.provenance import do
import bcbio.pipeline.datadict as dd
from bcbio.variation import bedutils
from bcbio.variation import coverage as cov
from bcbio.ngsalign.postalign import dedup_bam
from bcbio.rnaseq import gtf

# ## High level functions to generate summary

def generate_parallel(samples, run_parallel):
    """Provide parallel preparation of summary information for alignment and variant calling.
    """
    samples = run_parallel("pipeline_summary", samples)
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
        data["summary"] = _run_qc_tools(work_bam, data)
    elif data["analysis"].lower().startswith("chip-seq"):
        work_bam = data["raw_bam"]
        data["summary"] = _run_qc_tools(work_bam, data)
    elif dd.get_ref_file(data) is not None and work_bam and work_bam.endswith(".bam"):
        logger.info("Generating summary files: %s" % dd.get_sample_name(data))
        data["summary"] = _run_qc_tools(work_bam, data)
    return [[data]]

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

def _run_qc_tools(bam_file, data):
    """Run a set of third party quality control tools, returning QC directory and metrics.

        :param bam_file: alignments in bam format
        :param data: dict with all configuration information

        :returns: dict with output of different tools
    """
    metrics = {}
    to_run = []
    if "fastqc" not in tz.get_in(("config", "algorithm", "tools_off"), data, []):
        to_run.append(("fastqc", _run_fastqc))
    if data["analysis"].lower().startswith("rna-seq"):
        to_run += [("samtools", _run_samtools_stats)]
        if gtf.is_qualimap_compatible(dd.get_gtf_file(data)):
            to_run.append(("qualimap", _rnaseq_qualimap))
    elif data["analysis"].lower().startswith("chip-seq"):
        bam_file = bam_file.replace(".unique", "")
        to_run += [("samtools", _run_samtools_stats)]
    elif not data["analysis"].lower().startswith("smallrna-seq"):
        to_run += [("samtools", _run_samtools_stats), ("gemini", _run_gemini_stats)]

    if data["analysis"].lower().startswith(("standard", "variant", "variant2")):
        to_run.append(["qsignature", _run_qsignature_generator])
        to_run.append(["coverage", _run_coverage_qc])
        to_run.append(["variants", _run_variants_qc])
        if any([tool in tz.get_in(("config", "algorithm", "tools_on"), data, []) for tool in ["qualimap", "qualimap_full"]]):
            to_run.append(("qualimap", _run_qualimap))
    qc_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "qc", data["description"]))
    metrics = {}
    qc_out = {}
    for program_name, qc_fn in to_run:
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
    if data['config']["algorithm"].get("kraken", None):
        if data["analysis"].lower().startswith("smallrna-seq"):
            logger.info("Kraken is not compatible with srnaseq pipeline yet.")
        else:
            ratio = bam.get_aligned_reads(bam_file, data)
            cur_metrics = _run_kraken(data, ratio)
            metrics.update(cur_metrics)

    bam.remove("%s-downsample%s" % os.path.splitext(bam_file))

    metrics["Name"] = dd.get_sample_name(data)
    metrics["Quality format"] = dd.get_quality_format(data).lower()
    return {"qc": qc_out, "metrics": metrics}

def _organize_qc_files(program, qc_dir):
    """Organize outputs from quality control runs into a base file and secondary outputs.

    Provides compatibility with CWL output. Returns None if no files created during processing.
    """
    base_files = {"fastqc": "fastqc_report.html", "qualimap": "qualimapReport.html"}
    if os.path.exists(qc_dir):
        out_files = []
        for fname in [os.path.join(qc_dir, x) for x in os.listdir(qc_dir)]:
            if os.path.isfile(fname):
                out_files.append(fname)
            elif os.path.isdir(fname) and not fname.endswith("tx"):
                out_files.extend([os.path.join(fname, x) for x in os.listdir(fname)])
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
    if os.path.exists(summary_file):
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
        # check if disambiguation was run
        if "disambiguate" in sample:
            if utils.file_exists(sample["disambiguate"]["summary"]):
                disambigStats = _parse_disambiguate(sample["disambiguate"]["summary"])
                saved["summary"]["metrics"]["Disambiguated %s reads" % str(sample["genome_build"])] = disambigStats[0]
                disambigGenome = (sample["config"]["algorithm"]["disambiguate"][0]
                                  if isinstance(sample["config"]["algorithm"]["disambiguate"], (list, tuple))
                                  else sample["config"]["algorithm"]["disambiguate"])
                saved["summary"]["metrics"]["Disambiguated %s reads" % disambigGenome] = disambigStats[1]
                saved["summary"]["metrics"]["Disambiguated ambiguous reads"] = disambigStats[2]
    return saved

def _parse_disambiguate(disambiguatestatsfilename):
    """Parse disambiguation stats from given file.
    """
    disambig_stats = [0, 0, 0]
    with open(disambiguatestatsfilename, "r") as in_handle:
        for i, line in enumerate(in_handle):
            fields = line.strip().split("\t")
            if i == 0:
                assert fields == ['sample', 'unique species A pairs', 'unique species B pairs', 'ambiguous pairs']
            else:
                disambig_stats = [x + int(y) for x, y in zip(disambig_stats, fields[1:])]
    return disambig_stats

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
    metrics = ["Total reads", "Mapped reads", "Mapped reads pct", "Duplicates", "Duplicates pct"]
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

# ## Run and parse read information from FastQC

class FastQCParser:
    def __init__(self, base_dir, sample=None):
        self._dir = base_dir
        self.sample = sample

    def get_fastqc_summary(self):
        ignore = set(["Total Sequences", "Filtered Sequences",
                      "Filename", "File type", "Encoding"])
        stats = {}
        for stat_line in self._fastqc_data_section("Basic Statistics")[1:]:
            k, v = stat_line.split("\t")[:2]
            if k not in ignore:
                stats[k] = v
        return stats

    def _fastqc_data_section(self, section_name):
        out = []
        in_section = False
        data_file = os.path.join(self._dir, "fastqc_data.txt")
        if os.path.exists(data_file):
            with open(data_file) as in_handle:
                for line in in_handle:
                    if line.startswith(">>%s" % section_name):
                        in_section = True
                    elif in_section:
                        if line.startswith(">>END"):
                            break
                        out.append(line.rstrip("\r\n"))
        return out

    def save_sections_into_file(self):

        data_file = os.path.join(self._dir, "fastqc_data.txt")
        if os.path.exists(data_file) and Fadapa:
            parser = Fadapa(data_file)
            module = [m[1] for m in parser.summary()][2:9]
            for m in module:
                out_file = os.path.join(self._dir, m.replace(" ", "_") + ".tsv")
                dt = self._get_module(parser, m)
                dt.to_csv(out_file, sep="\t", index=False)

    def _get_module(self, parser, module):
        """
        Get module using fadapa package
        """
        dt = []
        lines = parser.clean_data(module)
        header = lines[0]
        for data in lines[1:]:
            if data[0].startswith("#"): #some modules have two headers
                header = data
                continue
            if data[0].find("-") > -1: # expand positions 1-3 to 1, 2, 3
                f, s = map(int, data[0].split("-"))
                for pos in range(f, s):
                    dt.append([str(pos)] + data[1:])
            else:
                dt.append(data)
        dt = pd.DataFrame(dt)
        dt.columns = [h.replace(" ", "_") for h in header]
        dt['sample'] = self.sample
        return dt

def _run_gene_coverage(bam_file, data, out_dir):
    out_file = os.path.join(out_dir, "gene_coverage.pdf")
    ref_file = utils.get_in(data, ("genome_resources", "rnaseq", "transcripts"))
    count_file = data["count_file"]
    if utils.file_exists(out_file):
        return out_file
    with file_transaction(data, out_file) as tx_out_file:
        plot_gene_coverage(bam_file, ref_file, count_file, tx_out_file)
    return {"gene_coverage": out_file}

def _run_kraken(data, ratio):
    """Run kraken, generating report in specified directory and parsing metrics.
       Using only first paired reads.
    """
    # logger.info("Number of aligned reads < than 0.60 in %s: %s" % (dd.get_sample_name(data), ratio))
    logger.info("Running kraken to determine contaminant: %s" % dd.get_sample_name(data))
    qc_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "qc", data["description"]))
    kraken_out = os.path.join(qc_dir, "kraken")
    out = out_stats = None
    db = data['config']["algorithm"]["kraken"]
    kraken_cmd = config_utils.get_program("kraken", data["config"])
    if db == "minikraken":
        db = os.path.join(_get_data_dir(), "genomes", "kraken", "minikraken")

    if not os.path.exists(db):
        logger.info("kraken: no database found %s, skipping" % db)
        return {"kraken_report": "null"}

    if not os.path.exists(os.path.join(kraken_out, "kraken_out")):
        work_dir = os.path.dirname(kraken_out)
        utils.safe_makedir(work_dir)
        num_cores = data["config"]["algorithm"].get("num_cores", 1)
        fn_file = data["files"][0]
        if fn_file.endswith("bam"):
            logger.info("kraken: need fasta files as input")
            return {"kraken_report": "null"}
        with tx_tmpdir(data, work_dir) as tx_tmp_dir:
            with utils.chdir(tx_tmp_dir):
                out = os.path.join(tx_tmp_dir, "kraken_out")
                out_stats = os.path.join(tx_tmp_dir, "kraken_stats")
                cat = "zcat" if fn_file.endswith(".gz") else "cat"
                cl = ("{cat} {fn_file} | {kraken_cmd} --db {db} --quick "
                      "--preload --min-hits 2 "
                      "--threads {num_cores} "
                      "--out {out} --fastq-input /dev/stdin  2> {out_stats}").format(**locals())
                do.run(cl, "kraken: %s" % dd.get_sample_name(data))
                if os.path.exists(kraken_out):
                    shutil.rmtree(kraken_out)
                shutil.move(tx_tmp_dir, kraken_out)
    metrics = _parse_kraken_output(kraken_out, db, data)
    return metrics

def _parse_kraken_output(out_dir, db, data):
    """Parse kraken stat info comming from stderr,
       generating report with kraken-report
    """
    in_file = os.path.join(out_dir, "kraken_out")
    stat_file = os.path.join(out_dir, "kraken_stats")
    out_file = os.path.join(out_dir, "kraken_summary")
    kraken_cmd = config_utils.get_program("kraken-report", data["config"])
    classify = unclassify = None
    with open(stat_file, 'r') as handle:
        for line in handle:
            if line.find(" classified") > -1:
                classify = line[line.find("(") + 1:line.find(")")]
            if line.find(" unclassified") > -1:
                unclassify = line[line.find("(") + 1:line.find(")")]
    if os.path.getsize(in_file) > 0 and not os.path.exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cl = ("{kraken_cmd} --db {db} {in_file} > {tx_out_file}").format(**locals())
            do.run(cl, "kraken report: %s" % dd.get_sample_name(data))
    kraken = {"kraken_clas": classify, "kraken_unclas": unclassify}
    kraken_sum = _summarize_kraken(out_file)
    kraken.update(kraken_sum)
    return kraken

def _summarize_kraken(fn):
    """get the value at species level"""
    kraken = {}
    list_sp, list_value = [], []
    with open(fn) as handle:
        for line in handle:
            cols = line.strip().split("\t")
            sp = cols[5].strip()
            if len(sp.split(" ")) > 1 and not sp.startswith("cellular"):
                list_sp.append(sp)
                list_value.append(cols[0])
    kraken = {"kraken_sp": list_sp, "kraken_value": list_value}
    return kraken

def _run_fastqc(bam_file, data, fastqc_out):
    """Run fastqc, generating report in specified directory and parsing metrics.

    Downsamples to 10 million reads to avoid excessive processing times with large
    files, unless we're running a Standard/smallRNA-seq/QC pipeline.

    Handles fastqc 0.11+, which use a single HTML file and older versions that use
    a directory of files + images. The goal is to eventually move to only 0.11+
    """
    sentry_file = os.path.join(fastqc_out, "fastqc_report.html")
    if not os.path.exists(sentry_file):
        work_dir = os.path.dirname(fastqc_out)
        utils.safe_makedir(work_dir)
        ds_bam = (bam.downsample(bam_file, data, 1e7)
                  if data.get("analysis", "").lower() not in ["standard", "smallrna-seq"]
                  else None)
        bam_file = ds_bam if ds_bam else bam_file
        frmt = "bam" if bam_file.endswith("bam") else "fastq"
        fastqc_name = utils.splitext_plus(os.path.basename(bam_file))[0]
        fastqc_clean_name =  dd.get_sample_name(data)
        num_cores = data["config"]["algorithm"].get("num_cores", 1)
        with tx_tmpdir(data, work_dir) as tx_tmp_dir:
            with utils.chdir(tx_tmp_dir):
                cl = [config_utils.get_program("fastqc", data["config"]),
                      "-d", tx_tmp_dir,
                      "-t", str(num_cores), "--extract", "-o", tx_tmp_dir, "-f", frmt, bam_file]
                do.run(cl, "FastQC: %s" % dd.get_sample_name(data))
                tx_fastqc_out = os.path.join(tx_tmp_dir, "%s_fastqc" % fastqc_name)
                tx_combo_file = os.path.join(tx_tmp_dir, "%s_fastqc.html" % fastqc_name)
                if not os.path.exists(sentry_file) and os.path.exists(tx_combo_file):
                    utils.safe_makedir(fastqc_out)
                    # Use sample name for reports instead of bam file name
                    with open(os.path.join(tx_fastqc_out, "fastqc_data.txt"), 'r') as fastqc_bam_name, \
                            open(os.path.join(tx_fastqc_out, "_fastqc_data.txt"), 'w') as fastqc_sample_name:
                        for line in fastqc_bam_name:
                            fastqc_sample_name.write(line.replace(os.path.basename(bam_file), fastqc_clean_name))
                    shutil.move(os.path.join(tx_fastqc_out, "_fastqc_data.txt"), os.path.join(fastqc_out, 'fastqc_data.txt'))
                    shutil.move(tx_combo_file, sentry_file)
                    if os.path.exists("%s.zip" % tx_fastqc_out):
                        shutil.move("%s.zip" % tx_fastqc_out, os.path.join(fastqc_out, "%s.zip" % fastqc_clean_name))
                elif not os.path.exists(sentry_file):
                    if os.path.exists(fastqc_out):
                        shutil.rmtree(fastqc_out)
                    shutil.move(tx_fastqc_out, fastqc_out)
    parser = FastQCParser(fastqc_out, dd.get_sample_name(data))
    stats = parser.get_fastqc_summary()
    parser.save_sections_into_file()
    return stats

# ## Qualimap

def _parse_num_pct(k, v):
    num, pct = v.split(" / ")
    return {k: num.replace(",", "").strip(), "%s pct" % k: pct.strip()}

def _parse_qualimap_globals(table):
    """Retrieve metrics of interest from globals table.
    """
    out = {}
    want = {"Mapped reads": _parse_num_pct,
            "Duplication rate": lambda k, v: {k: v}}
    for row in table.xpath("table/tr"):
        col, val = [x.text for x in row.xpath("td")]
        if col in want:
            out.update(want[col](col, val))
    return out

def _parse_qualimap_globals_inregion(table):
    """Retrieve metrics from the global targeted region table.
    """
    out = {}
    for row in table.xpath("table/tr"):
        col, val = [x.text for x in row.xpath("td")]
        if col == "Mapped reads":
            out.update(_parse_num_pct("%s (in regions)" % col, val))
    return out

def _parse_qualimap_coverage(table):
    """Parse summary qualimap coverage metrics.
    """
    out = {}
    for row in table.xpath("table/tr"):
        col, val = [x.text for x in row.xpath("td")]
        if col == "Mean":
            out["Coverage (Mean)"] = val
    return out

def _parse_qualimap_insertsize(table):
    """Parse insert size metrics.
    """
    out = {}
    for row in table.xpath("table/tr"):
        col, val = [x.text for x in row.xpath("td")]
        if col == "Median":
            out["Insert size (Median)"] = val
    return out

def _parse_qualimap_metrics(report_file):
    """Extract useful metrics from the qualimap HTML report file.
    """
    if not utils.file_exists(report_file):
        return {}
    out = {}
    parsers = {"Globals": _parse_qualimap_globals,
               "Globals (inside of regions)": _parse_qualimap_globals_inregion,
               "Coverage": _parse_qualimap_coverage,
               "Coverage (inside of regions)": _parse_qualimap_coverage,
               "Insert size": _parse_qualimap_insertsize,
               "Insert size (inside of regions)": _parse_qualimap_insertsize}
    root = lxml.html.parse(report_file).getroot()
    for table in root.xpath("//div[@class='table-summary']"):
        header = table.xpath("h3")[0].text
        if header in parsers:
            out.update(parsers[header](table))
    new_names = []
    for metric in out:
        new_names.append(metric + "_qualimap_1e7reads_est")
    out = dict(zip(new_names, out.values()))
    return out

def _bed_to_bed6(orig_file, out_dir):
    """Convert bed to required bed6 inputs.
    """
    bed6_file = os.path.join(out_dir, "%s-bed6%s" % os.path.splitext(os.path.basename(orig_file)))
    if not utils.file_exists(bed6_file):
        with open(bed6_file, "w") as out_handle:
            for i, region in enumerate(list(x) for x in pybedtools.BedTool(orig_file)):
                region = [x for x in list(region) if x]
                fillers = [str(i), "1.0", "+"]
                full = region + fillers[:6 - len(region)]
                out_handle.write("\t".join(full) + "\n")
    return bed6_file

def _run_qualimap(bam_file, data, out_dir):
    """Run qualimap to assess alignment quality metrics.
    """
    resources = config_utils.get_resources("qualimap", data["config"])
    options = " ".join(resources.get("options", ""))
    report_file = os.path.join(out_dir, "qualimapReport.html")
    pdf_file = "qualimapReport.pdf"
    if not utils.file_exists(report_file) and not utils.file_exists(os.path.join(out_dir, pdf_file)):
        if "qualimap_full" in tz.get_in(("config", "algorithm", "tools_on"), data, []):
            logger.info("Full qualimap analysis for %s may be slow." % bam_file)
            ds_bam = bam_file
        else:
            ds_bam = bam.downsample(bam_file, data, 1e7)
            bam_file = ds_bam if ds_bam else bam_file
        if options.find("PDF") > -1:
            options = "%s -outfile %s" % (options, pdf_file)
        utils.safe_makedir(out_dir)
        num_cores = data["config"]["algorithm"].get("num_cores", 1)
        qualimap = config_utils.get_program("qualimap", data["config"])
        max_mem = config_utils.adjust_memory(resources.get("memory", "1G"),
                                             num_cores)
        cmd = ("unset DISPLAY && {qualimap} bamqc -bam {bam_file} -outdir {out_dir} "
               "-nt {num_cores} --java-mem-size={max_mem} {options}")
        species = tz.get_in(("genome_resources", "aliases", "ensembl"), data, "")
        if species in ["HUMAN", "MOUSE"]:
            cmd += " -gd {species}"
        regions = bedutils.merge_overlaps(dd.get_variant_regions(data), data)
        if regions:
            bed6_regions = _bed_to_bed6(regions, out_dir)
            cmd += " -gff {bed6_regions}"
        do.run(cmd.format(**locals()), "Qualimap: %s" % dd.get_sample_name(data))

    return _parse_qualimap_metrics(report_file)

# ## RNAseq Qualimap

def _parse_metrics(metrics):
    # skipped metrics can sometimes be in unicode, replace unicode with NA if it exists
    metrics = dtz.valmap(lambda x: 'nan' if isinstance(x, unicode) else x, metrics)

    missing = set(["Genes Detected", "Transcripts Detected",
                   "Mean Per Base Cov."])
    correct = set(["rRNA", "rRNA_rate"])
    percentages = set(["Intergenic pct", "Intronic pct", "Exonic pct"])
    to_change = dict({"5'-3' bias": 1, "Intergenic pct": "Intergenic Rate",
                      "Intronic pct": "Intronic Rate",
                      "Exonic pct": "Exonic Rate",
                      "Duplication Rate of Mapped": 1,
                      "Average insert size": 1,
                      })
    total = ["Not aligned", "Aligned to genes", "No feature assigned"]

    out = {}
    total_reads = sum([int(metrics[name]) for name in total])
    out['Mapped'] = sum([int(metrics[name]) for name in total[1:]])
    out['Mapping Rate'] = 1.0 * int(out['Mapped']) / total_reads
    # [out.update({name: 0}) for name in missing]
    out.update({key: val for key, val in metrics.iteritems() if key in correct})
    [metrics.update({name: 1.0 * float(metrics[name]) / 100}) for name in
     percentages]

    for name in to_change:
        if not to_change[name]:
            continue
        try:
            if to_change[name] == 1:
                out.update({name: float(metrics[name])})
            else:
                out.update({to_change[name]: float(metrics[name])})
        # if we can't convert metrics[name] to float (?'s or other non-floats)
        except ValueError:
            continue
    return out

def _detect_duplicates(bam_file, out_dir, data):
    """
    count duplicate percentage
    """
    out_file = os.path.join(out_dir, "dup_metrics.txt")
    if not utils.file_exists(out_file):
        dup_align_bam = dedup_bam(bam_file, data)
        num_cores = dd.get_num_cores(data)
        with file_transaction(out_file) as tx_out_file:
            sambamba = config_utils.get_program("sambamba", data, default="sambamba")
            dup_count = ("{sambamba} view --nthreads {num_cores} --count "
                         "-F 'duplicate and not unmapped' "
                         "{dup_align_bam} >> {tx_out_file}")
            message = "Counting duplicates in {bam_file}.".format(bam_file=bam_file)
            do.run(dup_count.format(**locals()), message)
            tot_count = ("{sambamba} view --nthreads {num_cores} --count "
                         "-F 'not unmapped' "
                         "{dup_align_bam} >> {tx_out_file}")
            message = "Counting reads in {bam_file}.".format(bam_file=bam_file)
            do.run(tot_count.format(**locals()), message)
    with open(out_file) as in_handle:
        dupes = float(in_handle.next().strip())
        total = float(in_handle.next().strip())
    return {"Duplication Rate of Mapped": dupes / total}

def _transform_browser_coor(rRNA_interval, rRNA_coor):
    """
    transform interval format to browser coord: chr:start-end
    """
    with open(rRNA_coor, 'w') as out_handle:
        with open(rRNA_interval, 'r') as in_handle:
            for line in in_handle:
                c, bio, source, s, e = line.split("\t")[:5]
                if bio.startswith("rRNA"):
                    out_handle.write(("{0}:{1}-{2}\n").format(c, s, e))

def _detect_rRNA(data):
    sample = dd.get_sample_name(data)
    gtf_file = dd.get_gtf_file(data)
    tidy_file = dd.get_sailfish_tidy(data)
    rrna_features = gtf.get_rRNA(gtf_file)
    transcripts = set([x[1] for x in rrna_features if x])
    if not transcripts:
        return {'rRNA': "NA", "rRNA_rate": "NA"}
    count_table = pd.read_csv(tidy_file, sep="\t")
    sample_table = count_table[count_table["sample"].isin([sample])]
    rrna_exp = map(float, sample_table[sample_table["id"].isin(transcripts)]["numreads"])
    total_exp = map(float, sample_table["numreads"])
    rrna = sum(rrna_exp)
    rrna_rate = float(rrna) / sum(total_exp)
    return {'rRNA': str(rrna), 'rRNA_rate': str(rrna_rate)}

def _parse_qualimap_rnaseq(table):
    """
    Retrieve metrics of interest from globals table.
    """
    out = {}
    for row in table.xpath("table/tr"):
        col, val = [x.text for x in row.xpath("td")]
        col = col.replace(":", "").strip()
        val = val.replace(",", "")
        m = {col: val}
        if val.find("/") > -1:
            m = _parse_num_pct(col, val.replace("%", ""))
        out.update(m)
    return out

def _parse_rnaseq_qualimap_metrics(report_file):
    """Extract useful metrics from the qualimap HTML report file.
    """
    out = {}
    parsers = ["Reads alignment", "Reads genomic origin", "Transcript coverage profile"]
    root = lxml.html.parse(report_file).getroot()
    for table in root.xpath("//div[@class='table-summary']"):
        header = table.xpath("h3")[0].text
        if header in parsers:
            out.update(_parse_qualimap_rnaseq(table))
    return out

def _rnaseq_qualimap(bam_file, data, out_dir):
    """
    Run qualimap for a rnaseq bam file and parse results
    """
    strandedness = {"firststrand": "strand-specific-reverse",
                    "secondstrand": "strand-specific-forward",
                    "unstranded": "non-strand-specific"}
    report_file = os.path.join(out_dir, "qualimapReport.html")
    config = data["config"]
    gtf_file = dd.get_gtf_file(data)
    ref_file = dd.get_ref_file(data)
    single_end = not bam.is_paired(bam_file)
    library = strandedness[dd.get_strandedness(data)]
    if not utils.file_exists(report_file):
        utils.safe_makedir(out_dir)
        bam.index(bam_file, config)
        cmd = _rnaseq_qualimap_cmd(config, bam_file, out_dir, gtf_file, single_end, library)
        do.run(cmd, "Qualimap for {}".format(dd.get_sample_name(data)))
    metrics = _parse_rnaseq_qualimap_metrics(report_file)
    metrics.update(_detect_duplicates(bam_file, out_dir, data))
    metrics.update(_detect_rRNA(data))
    metrics.update({"Average insert size": bam.estimate_fragment_size(bam_file)})
    metrics = _parse_metrics(metrics)
    return metrics

def _rnaseq_qualimap_cmd(config, bam_file, out_dir, gtf_file=None, single_end=None, library="non-strand-specific"):
    """
    Create command lines for qualimap
    """
    qualimap = config_utils.get_program("qualimap", config)
    resources = config_utils.get_resources("qualimap", config)
    num_cores = resources.get("cores", 1)
    max_mem = config_utils.adjust_memory(resources.get("memory", "4G"),
                                         num_cores)
    cmd = ("unset DISPLAY && {qualimap} rnaseq -outdir {out_dir} -a proportional -bam {bam_file} -p {library} "
           "-gtf {gtf_file} --java-mem-size={max_mem}").format(**locals())
    return cmd

# ## Lightweight QC approaches

def _parse_samtools_stats(stats_file):
    out = {}
    want = {"raw total sequences": "Total reads", "reads mapped": "Mapped reads",
            "reads duplicated": "Duplicates", "insert size average": "Average insert size"}
    with open(stats_file) as in_handle:
        for line in in_handle:
            if not line.startswith("SN"):
                continue
            parts = line.split("\t")
            metric, stat_str = parts[1:3]
            metric = metric.replace(":", "").strip()
            if metric in want:
                stat = float(stat_str.strip())
                out[want[metric]] = stat
                if metric in ["reads mapped", "reads duplicated"]:
                    out["%s pct" % want[metric]] = stat / out["Total reads"]
    return out

def _parse_offtargets(bam_file):
    """
    Add to metrics off-targets reads if it exitst
    """
    off_target = bam_file.replace(".bam", "-offtarget-stats.yaml")
    if os.path.exists(off_target):
        res = yaml.load(open(off_target))
        res['offtarget_pct'] = "%.3f" % (float(res['offtarget']) / float(res['mapped']))
        return res
    return {}

def _run_bamtools_stats(bam_file, data, out_dir):
    """Run bamtools stats with reports on mapped reads, duplicates and insert sizes.
    To be replaced by samtools. Kept in case of back-compatibility issues.
    """
    stats_file = os.path.join(out_dir, "bamtools_stats.txt")
    if not utils.file_exists(stats_file):
        utils.safe_makedir(out_dir)
        bamtools = config_utils.get_program("bamtools", data["config"])
        with file_transaction(data, stats_file) as tx_out_file:
            cmd = "{bamtools} stats -in {bam_file}"
            if bam.is_paired(bam_file):
                cmd += " -insert"
            cmd += " > {tx_out_file}"
            do.run(cmd.format(**locals()), "bamtools stats", data)
    out = _parse_bamtools_stats(stats_file)
    out.update(_parse_offtargets(bam_file))
    return out

def _run_samtools_stats(bam_file, data, out_dir):
    """Run samtools stats with reports on mapped reads, duplicates and insert sizes.
    """
    stats_file = os.path.join(out_dir, "%s.txt" % dd.get_sample_name(data))
    if not utils.file_exists(stats_file):
        utils.safe_makedir(out_dir)
        samtools = config_utils.get_program("samtools", data["config"])
        with file_transaction(data, stats_file) as tx_out_file:
            cmd = "{samtools} stats {bam_file}"
            cmd += " > {tx_out_file}"
            do.run(cmd.format(**locals()), "samtools stats", data)
    out = _parse_samtools_stats(stats_file)
    out.update(_parse_offtargets(bam_file))
    return out

## Variant statistics from gemini

def _run_gemini_stats(bam_file, data, out_dir):
    """Retrieve high level variant statistics from Gemini.
    """
    out = {}
    gemini_dbs = [d for d in
                  [tz.get_in(["population", "db"], x) for x in data.get("variants", [])] if d]
    if len(gemini_dbs) > 0:
        gemini_db = gemini_dbs[0]
        gemini_stat_file = "%s-stats.yaml" % os.path.splitext(gemini_db)[0]
        if not utils.file_uptodate(gemini_stat_file, gemini_db):
            gemini = config_utils.get_program("gemini", data["config"])
            tstv = subprocess.check_output([gemini, "stats", "--tstv", gemini_db])
            gt_counts = subprocess.check_output([gemini, "stats", "--gts-by-sample", gemini_db])
            dbsnp_count = subprocess.check_output([gemini, "query", gemini_db, "-q",
                                                   "SELECT count(*) FROM variants WHERE in_dbsnp==1"])
            out["Transition/Transversion"] = tstv.split("\n")[1].split()[-1]
            for line in gt_counts.split("\n"):
                parts = line.rstrip().split()
                if len(parts) > 0 and parts[0] != "sample":
                    name, hom_ref, het, hom_var, _, total = parts
                    out[name] = {}
                    out[name]["Variations (heterozygous)"] = int(het)
                    out[name]["Variations (homozygous)"] = int(hom_var)
                    # same total variations for all samples, keep that top level as well.
                    out["Variations (total)"] = int(total)
            out["Variations (in dbSNP)"] = int(dbsnp_count.strip())
            if out.get("Variations (total)") > 0:
                out["Variations (in dbSNP) pct"] = "%.1f%%" % (out["Variations (in dbSNP)"] /
                                                               float(out["Variations (total)"]) * 100.0)
            with open(gemini_stat_file, "w") as out_handle:
                yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
        else:
            with open(gemini_stat_file) as in_handle:
                out = yaml.safe_load(in_handle)
    else:
        vcf_file = dd.get_vrn_file(data)
        if isinstance(vcf_file, list):
            vcf_file = vcf_file[0]
        if vcf_file:
            out_file = "%s-bcfstats.tsv" % utils.splitext_plus(vcf_file)[0]
            bcftools = config_utils.get_program("bcftools", data["config"])
            if not utils.file_exists(out_file):
                cmd = ("{bcftools} stats -f PASS {vcf_file} > {out_file}")
                do.run(cmd.format(**locals()), "basic vcf stats %s" % dd.get_sample_name(data))
            with open(out_file) as in_handle:
                for line in in_handle:
                    if line.startswith("SN") and line.find("records") > -1:
                        cols = line.split()
                        out["Variations (total)"] = cols[-1]

    res = {}
    for k, v in out.iteritems():
        if not isinstance(v, dict):
            res.update({k: v})
        if k == dd.get_sample_name(data):
            res.update(v)
    return res

## multiqc

def _check_multiqc_input(path):
    """Check if dir exists, and return empty if it doesn't"""
    if len(glob.glob(path)) > 0:
        return path
    return ""

def multiqc_summary(*samples):
    """Summarize all quality metrics together"""
    samples = utils.unpack_worlds(samples)
    work_dir = dd.get_work_dir(samples[0])
    multiqc = config_utils.get_program("multiqc", samples[0]["config"])
    if not multiqc:
        logger.debug("multiqc not found. Update bcbio_nextgen.py tools to fix this issue.")
    folders = []
    opts = ""
    out_dir = os.path.join(work_dir, "multiqc")
    out_data = os.path.join(work_dir, "multiqc", "multiqc_data")
    out_file = os.path.join(out_dir, "multiqc_report.html")
    samples = report_summary(samples, os.path.join(out_dir, "report"))
    for data in samples:
        for program, pfiles in tz.get_in(["summary", "qc"], data, {}).iteritems():
            if isinstance(pfiles, dict):
                pfiles = pfiles["base"]
            folders.append(os.path.dirname(pfiles))
    # XXX temporary workaround until we can handle larger inputs through MultiQC
    if len(folders) > 250:
        logger.warning("Too many samples for MultiQC, only using first 250 entries.")
        folders = folders[:250]
        opts = "--flat"
    # Back compatible -- to migrate to explicit specifications in input YAML
    folders += ["trimmed", "htseq-count/*summary"]
    if not utils.file_exists(out_file):
        with utils.chdir(work_dir):
            input_dir = " ".join([_check_multiqc_input(d) for d in folders])
            if input_dir.strip():
                cmd = "{multiqc} -f {input_dir} -o {tx_out} {opts}"
                with tx_tmpdir(data, work_dir) as tx_out:
                    do.run(cmd.format(**locals()), "Run multiqc")
                    shutil.move(os.path.join(tx_out, "multiqc_report.html"), out_file)
                    shutil.move(os.path.join(tx_out, "multiqc_data"), out_data)
    out = []
    for i, data in enumerate(samples):
        if i == 0:
            if utils.file_exists(out_file):
                data_files = glob.glob(os.path.join(out_dir, "multiqc_data", "*.txt"))
                data_files += glob.glob(os.path.join(out_dir, "report", "*", "*.bed"))
                data_files += glob.glob(os.path.join(out_dir, "report", "*", "*.txt"))
                data_files += glob.glob(os.path.join(out_dir, "report", "*", "*.tsv"))
                data_files += glob.glob(os.path.join(out_dir, "report", "*.R*"))
                if "summary" not in data:
                    data["summary"] = {}
                data["summary"]["multiqc"] = {"base": out_file, "secondary": data_files}
        out.append(data)
    return [[d] for d in out]

## qsignature

def _run_qsignature_generator(bam_file, data, out_dir):
    """ Run SignatureGenerator to create normalize vcf that later will be input of qsignature_summary

    :param bam_file: (str) path of the bam_file
    :param data: (list) list containing the all the dictionary
                     for this sample
    :param out_dir: (str) path of the output

    :returns: (dict) dict with the normalize vcf file
    """
    qsig = config_utils.get_program("qsignature", data["config"])
    if not qsig:
        logger.info("There is no qsignature tool. Skipping...")
        return None

    position = dd.get_qsig_file(data)
    mixup_check = dd.get_mixup_check(data)
    if mixup_check and mixup_check.startswith("qsignature"):
        utils.safe_makedir(out_dir)
        if not position:
            logger.info("There is no qsignature for this species: %s"
                        % tz.get_in(['genome_build'], data))
            return None
        jvm_opts = "-Xms750m -Xmx2g"
        if mixup_check == "qsignature_full":
            jvm_opts = "-Xms750m -Xmx8g"
            down_bam = bam_file
        else:
            down_bam = _slice_bam_chr21(bam_file, data)
            position = _slice_vcf_chr21(position, out_dir)

        out_name = os.path.basename(down_bam).replace("bam", "qsig.vcf")
        out_file = os.path.join(out_dir, out_name)
        log_file = os.path.join(out_dir, "qsig.log")
        cores = dd.get_cores(data)
        base_cmd = ("{qsig} {jvm_opts} "
                    "org.qcmg.sig.SignatureGenerator "
                    "--noOfThreads {cores} "
                    "-log {log_file} -i {position} "
                    "-i {down_bam} ")
        if not os.path.exists(out_file):
            file_qsign_out = "{0}.qsig.vcf".format(down_bam)
            do.run(base_cmd.format(**locals()), "qsignature vcf generation: %s" % dd.get_sample_name(data))
            if os.path.exists(file_qsign_out):
                with file_transaction(data, out_file) as file_txt_out:
                    shutil.move(file_qsign_out, file_txt_out)
            else:
                raise IOError("File doesn't exist %s" % file_qsign_out)
        return out_file
    return None

def qsignature_summary(*samples):
    """Run SignatureCompareRelatedSimple module from qsignature tool.

    Creates a matrix of pairwise comparison among samples. The
    function will not run if the output exists

    :param samples: list with only one element containing all samples information
    :returns: (dict) with the path of the output to be joined to summary
    """
    warnings, similar = [], []
    qsig = config_utils.get_program("qsignature", samples[0][0]["config"])
    if not qsig:
        return [[]]
    jvm_opts = "-Xms750m -Xmx8g"
    work_dir = samples[0][0]["dirs"]["work"]
    count = 0
    for data in samples:
        data = data[0]
        vcf = tz.get_in(["summary", "qc", "qsignature", "base"], data)
        if vcf:
            count += 1
            vcf_name = dd.get_sample_name(data) + ".qsig.vcf"
            out_dir = utils.safe_makedir(os.path.join(work_dir, "qsignature"))
            if not os.path.lexists(os.path.join(out_dir, vcf_name)):
                os.symlink(vcf, os.path.join(out_dir, vcf_name))
    if count > 0:
        qc_out_dir = utils.safe_makedir(os.path.join(work_dir, "qc", "qsignature"))
        out_file = os.path.join(qc_out_dir, "qsignature.xml")
        out_ma_file = os.path.join(qc_out_dir, "qsignature.ma")
        out_warn_file = os.path.join(qc_out_dir, "qsignature.warnings")
        log = os.path.join(work_dir, "qsignature", "qsig-summary.log")
        if not os.path.exists(out_file):
            with file_transaction(samples[0][0], out_file) as file_txt_out:
                base_cmd = ("{qsig} {jvm_opts} "
                            "org.qcmg.sig.SignatureCompareRelatedSimple "
                            "-log {log} -dir {out_dir} "
                            "-o {file_txt_out} ")
                do.run(base_cmd.format(**locals()), "qsignature score calculation")
        error, warnings, similar = _parse_qsignature_output(out_file, out_ma_file,
                                                            out_warn_file, samples[0][0])
        return [{'total samples': count,
                 'similar samples pairs': len(similar),
                 'warnings samples pairs': len(warnings),
                 'error samples': list(error),
                 'out_dir': qc_out_dir}]
    else:
        return []

def _parse_qsignature_output(in_file, out_file, warning_file, data):
    """ Parse xml file produced by qsignature

    :param in_file: (str) with the path to the xml file
    :param out_file: (str) with the path to output file
    :param warning_file: (str) with the path to warning file

    :returns: (list) with samples that could be duplicated

    """
    name = {}
    error, warnings, similar = set(), set(), set()
    same, replicate, related = 0, 0.1, 0.18
    mixup_check = dd.get_mixup_check(data)
    if mixup_check == "qsignature_full":
        same, replicate, related = 0, 0.01, 0.061
    with open(in_file, 'r') as in_handle:
        with file_transaction(data, out_file) as out_tx_file:
            with file_transaction(data, warning_file) as warn_tx_file:
                with open(out_tx_file, 'w') as out_handle:
                    with open(warn_tx_file, 'w') as warn_handle:
                        et = lxml.etree.parse(in_handle)
                        for i in list(et.iter('file')):
                            name[i.attrib['id']] = os.path.basename(i.attrib['name']).replace(".qsig.vcf", "")
                        for i in list(et.iter('comparison')):
                            msg = None
                            pair = "-".join([name[i.attrib['file1']], name[i.attrib['file2']]])
                            out_handle.write("%s\t%s\t%s\n" %
                                             (name[i.attrib['file1']], name[i.attrib['file2']], i.attrib['score']))
                            if float(i.attrib['score']) == same:
                                msg = 'qsignature ERROR: read same samples:%s\n'
                                error.add(pair)
                            elif float(i.attrib['score']) < replicate:
                                msg = 'qsignature WARNING: read similar/replicate samples:%s\n'
                                warnings.add(pair)
                            elif float(i.attrib['score']) < related:
                                msg = 'qsignature NOTE: read relative samples:%s\n'
                                similar.add(pair)
                            if msg:
                                logger.info(msg % pair)
                                warn_handle.write(msg % pair)
    return error, warnings, similar

def _slice_bam_chr21(in_bam, data):
    """
    return only one BAM file with only chromosome 21
    """
    sambamba = config_utils.get_program("sambamba", data["config"])
    out_file = "%s-chr%s" % os.path.splitext(in_bam)
    if not utils.file_exists(out_file):
        bam.index(in_bam, data['config'])
        with contextlib.closing(pysam.Samfile(in_bam, "rb")) as bamfile:
            bam_contigs = [c["SN"] for c in bamfile.header["SQ"]]
        chromosome = "21"
        if "chr21" in bam_contigs:
            chromosome = "chr21"
        with file_transaction(data, out_file) as tx_out_file:
            cmd = ("{sambamba} slice -o {tx_out_file} {in_bam} {chromosome}").format(**locals())
            out = subprocess.check_output(cmd, shell=True)
    return out_file

def _slice_vcf_chr21(vcf_file, out_dir):
    """
    Slice chr21 of qsignature SNPs to reduce computation time
    """
    tmp_file = os.path.join(out_dir, "chr21_qsignature.vcf")
    if not utils.file_exists(tmp_file):
        cmd = ("grep chr21 {vcf_file} > {tmp_file}").format(**locals())
        out = subprocess.check_output(cmd, shell=True)
    return tmp_file

## report and coverage
def report_summary(samples, out_dir):
    """
    Run coverage report with bcbiocov package
    """
    try:
        import bcbreport.prepare as bcbreport
    except ImportError:
        logger.info("skipping report. No bcbreport installed.")
        return samples
    # samples = utils.unpack_worlds(samples)
    work_dir = dd.get_work_dir(samples[0])
    parent_dir = utils.safe_makedir(out_dir)
    with utils.chdir(parent_dir):
        logger.info("copy qsignature")
        qsignature_fn = os.path.join(work_dir, "qc", "qsignature", "qsignature.ma")
        if qsignature_fn: # this need to be inside summary/qc dict
            if utils.file_exists(qsignature_fn) and not utils.file_exists("qsignature.ma"):
                shutil.copy(qsignature_fn, "bcbio_qsignature.ma")

        out_dir = utils.safe_makedir("fastqc")
        logger.info("summarize fastqc")
        with utils.chdir(out_dir):
            _merge_fastqc(samples)

        logger.info("summarize metrics")
        samples = _merge_metrics(samples)

        out_dir = utils.safe_makedir("coverage")
        logger.info("summarize coverage")
        for data in samples:
            pfiles = tz.get_in(["summary", "qc", "coverage"], data, [])
            if isinstance(pfiles, dict):
                pfiles = [pfiles["base"]] + pfiles["secondary"]
            elif pfiles:
                pfiles = [pfiles]
            for fn in pfiles:
                if os.path.basename(fn).find("coverage_fixed") > -1:
                    utils.copy_plus(fn, os.path.join(out_dir, os.path.basename(fn)))

        out_dir = utils.safe_makedir("variants")
        logger.info("summarize variants")
        for data in samples:
            pfiles = tz.get_in(["summary", "qc", "variants"], data, [])
            if isinstance(pfiles, dict):
                pfiles = [pfiles["base"]] + pfiles["secondary"]
            elif pfiles:
                pfiles = [pfiles]
            for fn in pfiles:
                if os.path.basename(fn).find("gc-depth-parse.tsv") > -1:
                    utils.copy_plus(fn, os.path.join(out_dir, os.path.basename(fn)))
        bcbreport.report(parent_dir)
        out_report = os.path.join(parent_dir, "qc-coverage-report.html")
        if not utils.file_exists(out_report):
            rmd_file = os.path.join(parent_dir, "report-ready.Rmd")
            run_file = "%s-run.R" % (os.path.splitext(out_report)[0])
            with open(run_file, "w") as out_handle:
                out_handle.write("""library(rmarkdown)\nrender("%s")\n""" % rmd_file)
            cmd = "%s %s" % (utils.Rscript_cmd(), run_file)
            # Skip automated generation of coverage report to avoid error
            # messages. We need to generalize coverage reporting and re-include.
            # try:
            #     do.run(cmd, "Prepare coverage summary", log_error=False)
# except subprocess.CalledProcessError, msg:
            #     logger.info("Skipping generation of coverage report: %s" % (str(msg)))
            if utils.file_exists("report-ready.html"):
                shutil.move("report-ready.html", out_report)
    return samples

def _run_coverage_qc(bam_file, data, out_dir):
    """Run coverage QC analysis"""
    priority = cov.priority_coverage(data, out_dir)
    cov.priority_total_coverage(data, out_dir)
    coverage = cov.coverage(data, out_dir)
    problem_regions = dd.get_problem_region_dir(data)
    annotated = None
    if problem_regions and priority:
        annotated = cov.decorate_problem_regions(priority, problem_regions)
    return None

def _run_variants_qc(bam_file, data, out_dir):
    """Run variants QC analysis"""
    return cov.variants(data, out_dir)

def _get_coverage_per_region(name):
    """
    Parse coverage file if it exists to get average value.
    """
    fn = os.path.join("coverage", name + "_coverage_fixed.bed")
    if utils.file_exists(fn):
        try:
            dt = pd.read_csv(fn, sep="\t", index_col=False)
            if len(dt["meanCoverage"]) > 0:
                return "%.3f" % (sum(map(float, dt['meanCoverage'])) / len(dt['meanCoverage']))
        except TypeError:
            logger.debug("%s has no lines in coverage.bed" % name)
    return "NA"

def _merge_metrics(samples):
    """
    parse project.yaml file to get metrics for each bam
    """
    out_file = os.path.join("metrics", "metrics.tsv")
    dt_together = []
    cov = {}
    with file_transaction(out_file) as out_tx:
        for s in samples:
            sample_name = dd.get_sample_name(s)
            if sample_name in cov:
                continue
            m = tz.get_in(['summary', 'metrics'], s)
            sample_file = os.path.abspath(os.path.join("metrics", "%s_bcbio.txt" % sample_name))
            if not tz.get_in(['summary', 'qc'], s):
                s['summary'] = {"qc": {}}
            if m:
                for me in m:
                    if isinstance(m[me], list) or isinstance(m[me], dict) or isinstance(m[me], tuple):
                        continue
                dt = pd.DataFrame(m, index=['1'])
                dt['avg_coverage_per_region'] = _get_coverage_per_region(sample_name)
                cov[sample_name] = dt['avg_coverage_per_region'][0]
                dt.columns = [k.replace(" ", "_").replace("(", "").replace(")", "") for k in dt.columns]
                dt['sample'] = sample_name
                dt.transpose().to_csv(sample_file, sep="\t", header=False)
                dt_together.append(dt)
                s['summary']['qc'].update({'bcbio':{'base': sample_file, 'secondary': []}})
        if len(dt_together) > 0:
            dt_together = utils.rbind(dt_together)
            dt_together.to_csv(out_tx, index=False, sep="\t")

    out = []
    for s in samples:
        if sample_name in cov:
            s['summary']['metrics']['avg_coverage_per_region'] = cov[sample_name]
        out.append(s)
    return samples

def _merge_fastqc(samples):
    """
    merge all fastqc samples into one by module
    """
    fastqc_list = defaultdict(list)
    seen = set()
    for data in samples:
        name = dd.get_sample_name(data)
        if name in seen:
            continue
        seen.add(name)
        fns = glob.glob(os.path.join(dd.get_work_dir(data), "qc", dd.get_sample_name(data), "fastqc") + "/*")
        for fn in fns:
            if fn.endswith("tsv"):
                metric = os.path.basename(fn)
                fastqc_list[metric].append([name, fn])

    for metric in fastqc_list:
        dt_by_sample = []
        for fn in fastqc_list[metric]:
            dt = pd.read_csv(fn[1], sep="\t")
            dt['sample'] = fn[0]
            dt_by_sample.append(dt)
        dt = utils.rbind(dt_by_sample)
        dt.to_csv(metric, sep="\t", index=False, mode ='w')
    return samples
