"""Quality control and summary metrics for next-gen alignments and analysis.
"""
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
# allow graceful during upgrades
try:
    import matplotlib
    matplotlib.use('Agg', force=True)
    import matplotlib.pyplot as plt
    plt.ioff()
except ImportError:
    plt = None
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
import bcbio.rnaseq.qc
import bcbio.pipeline.datadict as dd
from bcbio.variation import bedutils
from bcbio import broad
from bcbio.variation import coverage_experimental as cov
from bcbio.variation.coverage import decorate_problem_regions
# ## High level functions to generate summary


def generate_parallel(samples, run_parallel):
    """Provide parallel preparation of summary information for alignment and variant calling.
    """
    sum_samples = run_parallel("pipeline_summary", samples)
    qsign_info = run_parallel("qsignature_summary", [sum_samples])
    summary_file = write_project_summary(sum_samples, qsign_info)
    samples = []
    for data in sum_samples:
        if "summary" not in data[0]:
            data[0]["summary"] = {}
        data[0]["summary"]["project"] = summary_file
        if qsign_info:
            data[0]["summary"]["mixup_check"] = qsign_info[0]["out_dir"]
        samples.append(data)
    samples = _add_researcher_summary(samples, summary_file)
    return samples

def pipeline_summary(data):
    """Provide summary information on processing sample.
    """
    work_bam = data.get("work_bam")
    if data["sam_ref"] is not None and work_bam and work_bam.endswith(".bam"):
        logger.info("Generating summary files: %s" % str(data["name"]))
        data["summary"] = _run_qc_tools(work_bam, data)
    elif data["analysis"].lower().startswith("smallrna-seq"):
        work_bam = data["clean_fastq"]
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
        to_run.append(("qualimap", _rnaseq_qualimap))
    elif data["analysis"].lower().startswith("chip-seq"):
        to_run.append(["bamtools", _run_bamtools_stats])
    elif not data["analysis"].lower().startswith("smallrna-seq"):
        to_run += [("bamtools", _run_bamtools_stats), ("gemini", _run_gemini_stats)]
    if data["analysis"].lower().startswith(("standard", "variant2")):
        to_run.append(["qsignature", _run_qsignature_generator])
        if "qualimap" in tz.get_in(("config", "algorithm", "tools_on"), data, []):
            to_run.append(("qualimap", _run_qualimap))
    qc_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "qc", data["description"]))
    metrics = {}
    for program_name, qc_fn in to_run:
        cur_qc_dir = os.path.join(qc_dir, program_name)
        cur_metrics = qc_fn(bam_file, data, cur_qc_dir)
        metrics.update(cur_metrics)
    if data['config']["algorithm"].get("kraken", None):
        ratio = bam.get_aligned_reads(bam_file, data)
        cur_metrics = _run_kraken(data, ratio)
        metrics.update(cur_metrics)

    bam.remove("%s-downsample%s" % os.path.splitext(bam_file))

    metrics["Name"] = data["name"][-1]
    metrics["Quality format"] = utils.get_in(data,
                                             ("config", "algorithm",
                                              "quality_format"),
                                             "standard").lower()
    return {"qc": qc_dir, "metrics": metrics}

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
    if "summary" in sample:
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
    # logger.info("Number of aligned reads < than 0.60 in %s: %s" % (str(data["name"]), ratio))
    logger.info("Running kraken to determine contaminant: %s" % str(data["name"]))
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
                do.run(cl, "kraken: %s" % data["name"][-1])
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
            do.run(cl, "kraken report: %s" % data["name"][-1])
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
        num_cores = data["config"]["algorithm"].get("num_cores", 1)
        with tx_tmpdir(data, work_dir) as tx_tmp_dir:
            with utils.chdir(tx_tmp_dir):
                cl = [config_utils.get_program("fastqc", data["config"]),
                      "-t", str(num_cores), "--extract", "-o", tx_tmp_dir, "-f", frmt, bam_file]
                do.run(cl, "FastQC: %s" % data["name"][-1])
                tx_fastqc_out = os.path.join(tx_tmp_dir, "%s_fastqc" % fastqc_name)
                tx_combo_file = os.path.join(tx_tmp_dir, "%s_fastqc.html" % fastqc_name)
                if not os.path.exists(sentry_file) and os.path.exists(tx_combo_file):
                    utils.safe_makedir(fastqc_out)
                    shutil.move(os.path.join(tx_fastqc_out, "fastqc_data.txt"), fastqc_out)
                    shutil.move(tx_combo_file, sentry_file)
                    if os.path.exists("%s.zip" % tx_fastqc_out):
                        shutil.move("%s.zip" % tx_fastqc_out, fastqc_out)
                elif not os.path.exists(sentry_file):
                    if os.path.exists(fastqc_out):
                        shutil.rmtree(fastqc_out)
                    shutil.move(tx_fastqc_out, fastqc_out)
    parser = FastQCParser(fastqc_out, data["name"][-1])
    stats = parser.get_fastqc_summary()
    parser.save_sections_into_file()
    return stats

def _run_complexity(bam_file, data, out_dir):
    try:
        import pandas as pd
        import statsmodels.formula.api as sm
    except ImportError:
        return {"Unique Starts Per Read": "NA"}

    SAMPLE_SIZE = 1000000
    base, _ = os.path.splitext(os.path.basename(bam_file))
    utils.safe_makedir(out_dir)
    out_file = os.path.join(out_dir, base + ".pdf")
    df = bcbio.rnaseq.qc.starts_by_depth(bam_file, data["config"], SAMPLE_SIZE)
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tmp_out_file:
            df.plot(x='reads', y='starts', title=bam_file + " complexity")
            fig = plt.gcf()
            fig.savefig(tmp_out_file)

    print "file saved as", out_file
    print "out_dir is", out_dir
    return bcbio.rnaseq.qc.estimate_library_complexity(df)


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
    report_file = os.path.join(out_dir, "qualimapReport.html")
    if not os.path.exists(report_file):
        ds_bam = bam.downsample(bam_file, data, 1e7)
        bam_file = ds_bam if ds_bam else bam_file
        utils.safe_makedir(out_dir)
        num_cores = data["config"]["algorithm"].get("num_cores", 1)
        qualimap = config_utils.get_program("qualimap", data["config"])
        resources = config_utils.get_resources("qualimap", data["config"])
        max_mem = config_utils.adjust_memory(resources.get("memory", "1G"),
                                             num_cores)
        cmd = ("unset DISPLAY && {qualimap} bamqc -bam {bam_file} -outdir {out_dir} "
               "-nt {num_cores} --java-mem-size={max_mem}")
        species = data["genome_resources"]["aliases"].get("ensembl", "").upper()
        if species in ["HUMAN", "MOUSE"]:
            cmd += " -gd {species}"
        regions = bedutils.merge_overlaps(dd.get_variant_regions(data), data)
        if regions:
            bed6_regions = _bed_to_bed6(regions, out_dir)
            cmd += " -gff {bed6_regions}"
        do.run(cmd.format(**locals()), "Qualimap: %s" % data["name"][-1])
    return _parse_qualimap_metrics(report_file)

# ## RNAseq Qualimap

def _parse_metrics(metrics):
    # skipped metrics can sometimes be in unicode, replace unicode with NA if it exists
    metrics = dtz.valmap(lambda x: 'nan' if isinstance(x, unicode) else x, metrics)

    missing = set(["Genes Detected", "Transcripts Detected",
                   "Mean Per Base Cov."])
    correct = set(["Intergenic pct", "Intronic pct", "Exonic pct"])
    to_change = dict({"5'-3' bias": 1, "Intergenic pct": "Intergenic Rate",
                      "Intronic pct": "Intronic Rate", "Exonic pct": "Exonic Rate",
                      "Not aligned": 0, 'Aligned to genes': 0, 'Non-unique alignment': 0,
                      "No feature assigned": 0, "Duplication Rate of Mapped": 1,
                      "Fragment Length Mean": 1,
                      "rRNA": 1, "Ambiguou alignment": 0})
    total = ["Not aligned", "Aligned to genes", "No feature assigned"]

    out = {}
    total_reads = sum([int(metrics[name]) for name in total])
    out['rRNA rate'] = 1.0 * int(metrics["rRNA"]) / total_reads
    out['Mapped'] = sum([int(metrics[name]) for name in total[1:]])
    out['Mapping Rate'] = 1.0 * int(out['Mapped']) / total_reads
    [out.update({name: 0}) for name in missing]
    [metrics.update({name: 1.0 * float(metrics[name]) / 100}) for name in correct]

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

def _detect_duplicates(bam_file, out_dir, config):
    """
    Detect duplicates metrics with Picard
    """
    out_file = os.path.join(out_dir, "dup_metrics")
    if not utils.file_exists(out_file):
        broad_runner = broad.runner_from_config(config)
        (dup_align_bam, metrics_file) = broad_runner.run_fn("picard_mark_duplicates", bam_file, remove_dups=True)
        shutil.move(metrics_file, out_file)
    metrics = []
    with open(out_file) as in_handle:
        reader = csv.reader(in_handle, dialect="excel-tab")
        for line in reader:
            if line and not line[0].startswith("#"):
                metrics.append(line)
    metrics = dict(zip(metrics[0], metrics[1]))
    return {"Duplication Rate of Mapped": metrics["PERCENT_DUPLICATION"]}

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

def _detect_rRNA(config, bam_file, rRNA_file, ref_file, out_dir, single_end):
    """
    Calculate rRNA with gatk-framework
    """
    if not utils.file_exists(rRNA_file):
        return {'rRNA': 0}
    out_file = os.path.join(out_dir, "rRNA.counts")
    if not utils.file_exists(out_file):
        out_file = _count_rRNA_reads(bam_file, out_file, ref_file, rRNA_file, single_end, config)
    with open(out_file) as in_handle:
        for line in in_handle:
            if line.find("CountReads counted") > -1:
                rRNA_reads = line.split()[6]
                break
    return {'rRNA': rRNA_reads}

def _count_rRNA_reads(in_bam, out_file, ref_file, rRNA_interval, single_end, config):
    """Use GATK counter to count reads in rRNA genes
    """
    bam.index(in_bam, config)
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            rRNA_coor = os.path.join(os.path.dirname(out_file), "rRNA.list")
            _transform_browser_coor(rRNA_interval, rRNA_coor)
            params = ["-T", "CountReads",
                      "-R", ref_file,
                      "-I", in_bam,
                      "-log", tx_out_file,
                      "-L", rRNA_coor,
                      "--filter_reads_with_N_cigar",
                      "-allowPotentiallyMisencodedQuals"]
            jvm_opts = broad.get_gatk_framework_opts(config)
            cmd = [config_utils.get_program("gatk-framework", config)] + jvm_opts + params
            do.run(cmd, "counts rRNA for %s" % in_bam)
        return out_file

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
    report_file = os.path.join(out_dir, "qualimapReport.html")
    config = data["config"]
    gtf_file = dd.get_gtf_file(data)
    ref_file = dd.get_ref_file(data)
    single_end = not bam.is_paired(bam_file)
    if not utils.file_exists(report_file):
        utils.safe_makedir(out_dir)
        bam.index(bam_file, config)
        cmd = _rnaseq_qualimap_cmd(config, bam_file, out_dir, gtf_file, single_end)
        do.run(cmd, "Qualimap for {}".format(data["name"][-1]))
    metrics = _parse_rnaseq_qualimap_metrics(report_file)
    metrics.update(_detect_duplicates(bam_file, out_dir, config))
    metrics.update(_detect_rRNA(config, bam_file, gtf_file, ref_file, out_dir, single_end))
    metrics.update({"Fragment Length Mean": bam.estimate_fragment_size(bam_file)})
    metrics = _parse_metrics(metrics)
    return metrics

def _rnaseq_qualimap_cmd(config, bam_file, out_dir, gtf_file=None, single_end=None):
    """
    Create command lines for qualimap
    """
    qualimap = config_utils.get_program("qualimap", config)
    resources = config_utils.get_resources("qualimap", config)
    num_cores = resources.get("cores", 1)
    max_mem = config_utils.adjust_memory(resources.get("memory", "4G"),
                                         num_cores)
    cmd = ("unset DISPLAY && {qualimap} rnaseq -outdir {out_dir} -a proportional -bam {bam_file} "
           "-gtf {gtf_file} --java-mem-size={max_mem}").format(**locals())
    return cmd

# ## Lightweight QC approaches

def _parse_bamtools_stats(stats_file):
    out = {}
    want = set(["Total reads", "Mapped reads", "Duplicates", "Median insert size"])
    with open(stats_file) as in_handle:
        for line in in_handle:
            parts = line.split(":")
            if len(parts) == 2:
                metric, stat_str = parts
                metric = metric.split("(")[0].strip()
                if metric in want:
                    stat_parts = stat_str.split()
                    if len(stat_parts) == 2:
                        stat, pct = stat_parts
                        pct = pct.replace("(", "").replace(")", "")
                    else:
                        stat = stat_parts[0]
                        pct = None
                    out[metric] = stat
                    if pct:
                        out["%s pct" % metric] = pct
    return out

def _parse_offtargets(bam_file):
    """
    Add to metrics off-targets reads if it exitst
    """
    off_target = bam_file.replace(".bam", "-offtarget-stats.yaml")
    if os.path.exists(off_target):
        res = yaml.load(open(off_target))
        return res
    return {}

def _run_bamtools_stats(bam_file, data, out_dir):
    """Run bamtools stats with reports on mapped reads, duplicates and insert sizes.
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

    res = {}
    for k, v in out.iteritems():
        if not isinstance(v, dict):
            res.update({k: v})
        if k == data['name'][-1]:
            res.update(v)
    return res

## qsignature

def _run_qsignature_generator(bam_file, data, out_dir):
    """ Run SignatureGenerator to create normalize vcf that later will be input of qsignature_summary

    :param bam_file: (str) path of the bam_file
    :param data: (list) list containing the all the dictionary
                     for this sample
    :param out_dir: (str) path of the output

    :returns: (dict) dict with the normalize vcf file
    """
    position = dd.get_qsig_file(data)
    mixup_check = dd.get_mixup_check(data)
    if mixup_check and mixup_check.startswith("qsignature"):
        if not position:
            logger.info("There is no qsignature for this species: %s"
                        % tz.get_in(['genome_build'], data))
            return {}
        jvm_opts = "-Xms750m -Xmx2g"
        limit_reads = 20000000
        if mixup_check == "qsignature_full":
            slice_bam = bam_file
            jvm_opts = "-Xms750m -Xmx8g"
            limit_reads = 100000000
        else:
            slice_bam = _slice_chr22(bam_file, data)
        qsig = config_utils.get_program("qsignature", data["config"])
        if not qsig:
            return {}
        utils.safe_makedir(out_dir)
        out_name = os.path.basename(slice_bam).replace("bam", "qsig.vcf")
        out_file = os.path.join(out_dir, out_name)
        log_file = os.path.join(out_dir, "qsig.log")
        cores = dd.get_cores(data)
        base_cmd = ("{qsig} {jvm_opts} "
                    "org.qcmg.sig.SignatureGenerator "
                    "--noOfThreads {cores} "
                    "-log {log_file} -i {position} "
                    "-i {down_file} ")
        if not os.path.exists(out_file):
            down_file = bam.downsample(slice_bam, data, limit_reads)
            if not down_file:
                down_file = slice_bam
            file_qsign_out = "{0}.qsig.vcf".format(down_file)
            do.run(base_cmd.format(**locals()), "qsignature vcf generation: %s" % data["name"][-1])
            if os.path.exists(file_qsign_out):
                with file_transaction(data, out_file) as file_txt_out:
                    shutil.move(file_qsign_out, file_txt_out)
            else:
                raise IOError("File doesn't exist %s" % file_qsign_out)
        return {'qsig_vcf': out_file}
    return {}

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
        vcf = tz.get_in(["summary", "metrics", "qsig_vcf"], data)
        if vcf:
            count += 1
            vcf_name = data["name"][-1] + ".qsig.vcf"
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

def _slice_chr22(in_bam, data):
    """
    return only one BAM file with only chromosome 22
    """
    sambamba = config_utils.get_program("sambamba", data["config"])
    out_file = "%s-chr%s" % os.path.splitext(in_bam)
    if not utils.file_exists(out_file):
        bam.index(in_bam, data['config'])
        with contextlib.closing(pysam.Samfile(in_bam, "rb")) as bamfile:
            bam_contigs = [c["SN"] for c in bamfile.header["SQ"]]
        chromosome = "22"
        if "chr22" in bam_contigs:
            chromosome = "chr22"
        with file_transaction(data, out_file) as tx_out_file:
            cmd = ("{sambamba} slice -o {tx_out_file} {in_bam} {chromosome}").format(**locals())
            out = subprocess.check_output(cmd, shell=True)
    return out_file

def report_summary(samples, run_parallel):
    """
    Run coverage report with bcbiocov package
    """
    work_dir = dd.get_work_dir(samples[0][0])
    yaml_file = os.path.join(work_dir, "project-summary.yaml")

    parent_dir = utils.safe_makedir(os.path.join(work_dir,"report"))
    qsignature_fn = os.path.join(work_dir, "qc", "qsignature", "qsignature.ma")
    with utils.chdir(parent_dir):

        logger.info("copy qsignature")
        if qsignature_fn:
            if utils.file_exists(qsignature_fn) and not utils.file_exists("qsignature.ma"):
                shutil.copy(qsignature_fn, "qsignature.ma")
        logger.info("summarize metrics")
        _merge_metrics(yaml.load(open(yaml_file)))

        out_dir = utils.safe_makedir("fastqc")
        logger.info("summarize fastqc")
        with utils.chdir(out_dir):
            _merge_fastqc(samples)

        out_dir = utils.safe_makedir("coverage")
        out_dir = utils.safe_makedir("variants")
        samples = run_parallel("coverage_report", samples)

        try:
            import bcbreport.prepare as bcbreport
            bcbreport.report(parent_dir)
        except:
            logger.info("skipping report. No bcbreport installed.")
            pass

    return samples

## report and coverage

def coverage_report(data):
    """
    Run heavy coverage and variants process in parallel
    """
    data = cov.coverage(data)
    data = cov.variants(data)
    data = cov.priority_coverage(data)
    problem_regions = dd.get_problem_region_dir(data)
    name = dd.get_sample_name(data)
    if "coverage" in data:
        coverage = data['coverage']
        annotated = None
        if problem_regions and coverage:
             annotated = decorate_problem_regions(coverage, problem_regions)
        data['coverage'] = {'all': coverage, 'problems': annotated}

    return [[data]]

def _merge_metrics(yaml_data):
    """
    parse project.yaml file to get metrics for each bam
    """
    project = yaml_data
    out_file = os.path.join("metrics", "metrics.tsv")
    dt_together = []
    with file_transaction(out_file) as out_tx:
        for s in project['samples']:
            m = tz.get_in(['summary', 'metrics'], s)
            if m:
                for me in m:
                    if isinstance(m[me], list):
                        m[me] = ":".join(m[me])
                dt = pd.DataFrame(m, index=['1'])
                # dt = pd.DataFrame.from_dict(m)
                dt.columns = [k.replace(" ", "_").replace("(", "").replace(")", "") for k in dt.columns]
                dt['sample'] = s['description']
                dt_together.append(dt)
        dt_together = utils.rbind(dt_together)
        dt_together.to_csv(out_tx, index=False, sep="\t")

def _merge_fastqc(data):
    """
    merge all fastqc samples into one by module
    """
    fastqc_list = defaultdict(list)
    for sample in data:
        name = dd.get_sample_name(sample[0])
        fns = glob.glob(os.path.join(dd.get_work_dir(sample[0]), "qc", dd.get_sample_name(sample[0]), "fastqc") + "/*")
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
        dt.to_csv(metric, sep="\t", index=False, mode = 'w')
    return [data]
