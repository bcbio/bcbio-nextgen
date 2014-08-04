"""Quality control and summary metrics for next-gen alignments and analysis.
"""
import collections
import csv
import os
import shutil
import subprocess

import lxml.html
import yaml
from datetime import datetime

# allow graceful during upgrades
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    plt = None
import pysam
import contextlib

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import config_utils, run_info
from bcbio.install import _get_data_dir
from bcbio.provenance import do
import bcbio.rnaseq.qc
from bcbio.variation.realign import has_aligned_reads
from bcbio.rnaseq.coverage import plot_gene_coverage
import bcbio.pipeline.datadict as dd

# ## High level functions to generate summary

def generate_parallel(samples, run_parallel):
    """Provide parallel preparation of summary information for alignment and variant calling.
    """
    sum_samples = run_parallel("pipeline_summary", samples)
    qsign_info = run_parallel("qsignature_summary",[sum_samples])
    summary_file = write_project_summary(sum_samples,qsign_info)
    samples = []
    for data in sum_samples:
        if "summary" not in data[0]:
            data[0]["summary"] = {}
        data[0]["summary"]["project"] = summary_file
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
    to_run = [("fastqc", _run_fastqc)]
    if data["analysis"].lower().startswith("rna-seq"):
        to_run.append(("rnaseqc", bcbio.rnaseq.qc.sample_summary))
#        to_run.append(("coverage", _run_gene_coverage))
        to_run.append(("complexity", _run_complexity))
    elif data["analysis"].lower().startswith("chip-seq"):
        to_run.append(["bamtools", _run_bamtools_stats])
    else:
        to_run += [("bamtools", _run_bamtools_stats), ("gemini", _run_gemini_stats)]
    if data["analysis"].lower().startswith("standard"):
        to_run.append(["qsignature",_run_qsignature_generator])
    qc_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "qc", data["description"]))
    metrics = {}
    for program_name, qc_fn in to_run:
        cur_qc_dir = os.path.join(qc_dir, program_name)
        cur_metrics = qc_fn(bam_file, data, cur_qc_dir)
        metrics.update(cur_metrics)
    ratio = bam.get_aligned_reads(bam_file, data)
    if ratio < 0.60 and data['config']["algorithm"].get("kraken", False) and data["analysis"].lower() == "rna-seq":
        cur_metrics = _run_kraken(data, ratio)
        metrics.update(cur_metrics)
    metrics["Name"] = data["name"][-1]
    metrics["Quality format"] = utils.get_in(data,
                                             ("config", "algorithm",
                                              "quality_format"),
                                             "standard").lower()
    return {"qc": qc_dir, "metrics": metrics}

# ## Generate project level QC summary for quickly assessing large projects

def write_project_summary(samples,qsign_info = False):
    """Write project summary information on the provided samples.
    write out dirs, genome resources,

    """
    work_dir = samples[0][0]["dirs"]["work"]
    out_file = os.path.join(work_dir, "project-summary.yaml")
    upload_dir = (os.path.join(work_dir, samples[0][0]["upload"]["dir"])
                  if "dir" in samples[0][0]["upload"] else "")
    test_run = samples[0][0].get("test_run", False)
    date = str(datetime.now())
    prev_samples = _other_pipeline_samples(out_file, samples)
    with open(out_file, "w") as out_handle:
        yaml.safe_dump({"date": date}, out_handle,
                       default_flow_style=False, allow_unicode=False)
        if test_run:
            yaml.safe_dump({"test_run": True}, out_handle, default_flow_style=False,
                           allow_unicode=False)
        if qsign_info:
            yaml.safe_dump({"qsignature": qsign_info}, out_handle, default_flow_style=False,
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
    disambig_stats = [-1, -1, -1]
    with open(disambiguatestatsfilename, "r") as in_handle:
        header = in_handle.readline().strip().split("\t")
        if header == ['sample', 'unique species A pairs', 'unique species B pairs', 'ambiguous pairs']:
            disambig_stats_tmp = in_handle.readline().strip().split("\t")[1:]
            if len(disambig_stats_tmp) == 3:
                disambig_stats = [int(x) for x in disambig_stats_tmp]
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
    def __init__(self, base_dir):
        self._dir = base_dir

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

def _run_gene_coverage(bam_file, data, out_dir):
    out_file = os.path.join(out_dir, "gene_coverage.pdf")
    ref_file = utils.get_in(data, ("genome_resources", "rnaseq", "transcripts"))
    count_file = data["count_file"]
    if utils.file_exists(out_file):
        return out_file
    with file_transaction(out_file) as tx_out_file:
        plot_gene_coverage(bam_file, ref_file, count_file, tx_out_file)
    return {"gene_coverage": out_file}



def _run_kraken(data,ratio):
    """Run kraken, generating report in specified directory and parsing metrics.
       Using only first paired reads.
    """
    logger.info("Number of aligned reads < than 0.60 in %s: %s" % (str(data["name"]),ratio))
    logger.info("Running kraken to determine contaminant: %s" % str(data["name"]))
    qc_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "qc", data["description"]))
    kraken_out = os.path.join(qc_dir, "kraken")
    stats = out = out_stats = None
    db = data['config']["algorithm"]["kraken"] 
    if db == "minikraken":
        db = os.path.join(_get_data_dir(),"genome","kraken","minikraken")
    else:
        if not os.path.exists(db):
            logger.info("kraken: no database found %s, skipping" % db)
            return {"kraken_report" : "null"}
    if not os.path.exists(os.path.join(kraken_out,"kraken_out")):
        work_dir = os.path.dirname(kraken_out)
        utils.safe_makedir(work_dir)
        num_cores = data["config"]["algorithm"].get("num_cores", 1)
        files = data["files"]        
        with utils.curdir_tmpdir(data, work_dir) as tx_tmp_dir:
            with utils.chdir(tx_tmp_dir):
                out = os.path.join(tx_tmp_dir,"kraken_out")
                out_stats = os.path.join(tx_tmp_dir,"kraken_stats")
                cl = (" ").join([config_utils.get_program("kraken", data["config"]),
                      "--db",db,"--quick",
                      "--preload","--min-hits","2","--threads",str(num_cores), 
                      "--out", out, files[0]," 2>",out_stats])
                do.run(cl,"kraken: %s" % data["name"][-1])
                if os.path.exists(kraken_out):
                    shutil.rmtree(kraken_out)
                shutil.move(tx_tmp_dir, kraken_out)
    metrics = _parse_kraken_output(kraken_out,db,data)
    return metrics

def _parse_kraken_output(out_dir, db, data):
    """Parse kraken stat info comming from stderr, 
       generating report with kraken-report
    """
    in_file = os.path.join(out_dir,"kraken_out")
    stat_file = os.path.join(out_dir,"kraken_stats")
    out_file = os.path.join(out_dir, "kraken_summary")
    classify = unclassify = None
    with open(stat_file,'r') as handle:
        for line in handle:
            if line.find(" classified") > -1:
                classify = line[line.find("(")+1:line.find(")")]
            if line.find(" unclassified") > -1:
                unclassify = line[line.find("(")+1:line.find(")")]
    if os.path.getsize(in_file)>0:
        with file_transaction(out_file) as tx_out_file:
            cl = (" ").join([config_utils.get_program("kraken-report", data["config"]),
                          "--db",db,in_file,">",tx_out_file])
            do.run(cl, "kraken report: %s" % data["name"][-1])
    return {"kraken_report" : out_file,"kraken_clas": classify,"kraken_unclas": unclassify}

def _run_fastqc(bam_file, data, fastqc_out):
    """Run fastqc, generating report in specified directory and parsing metrics.

    Downsamples to 10 million reads to avoid excessive processing times with large
    files, unless we're running a Standard/QC pipeline.

    Handles fastqc 0.11+, which use a single HTML file and older versions that use
    a directory of files + images. The goal is to eventually move to only 0.11+
    """
    sentry_file = os.path.join(fastqc_out, "fastqc_report.html")
    if not os.path.exists(sentry_file):
        work_dir = os.path.dirname(fastqc_out)
        utils.safe_makedir(work_dir)
        ds_bam = (bam.downsample(bam_file, data, 1e7)
                  if data.get("analysis", "").lower() not in ["standard"]
                  else None)
        bam_file = ds_bam if ds_bam else bam_file
        fastqc_name = os.path.splitext(os.path.basename(bam_file))[0]
        num_cores = data["config"]["algorithm"].get("num_cores", 1)
        with utils.curdir_tmpdir(data, work_dir) as tx_tmp_dir:
            with utils.chdir(tx_tmp_dir):
                cl = [config_utils.get_program("fastqc", data["config"]),
                      "-t", str(num_cores), "--extract", "-o", tx_tmp_dir, "-f", "bam", bam_file]
                do.run(cl, "FastQC: %s" % data["name"][-1])
                tx_fastqc_out = os.path.join(tx_tmp_dir, "%s_fastqc" % fastqc_name)
                tx_combo_file = os.path.join(tx_tmp_dir, "%s_fastqc.html" % fastqc_name)
                if os.path.exists("%s.zip" % tx_fastqc_out):
                    os.remove("%s.zip" % tx_fastqc_out)
                if not os.path.exists(sentry_file) and os.path.exists(tx_combo_file):
                    utils.safe_makedir(fastqc_out)
                    shutil.copy(os.path.join(tx_fastqc_out, "fastqc_data.txt"), fastqc_out)
                    shutil.move(tx_combo_file, sentry_file)
                elif not os.path.exists(sentry_file):
                    if os.path.exists(fastqc_out):
                        shutil.rmtree(fastqc_out)
                    shutil.move(tx_fastqc_out, fastqc_out)
        if ds_bam and os.path.exists(ds_bam):
            os.remove(ds_bam)
    parser = FastQCParser(fastqc_out)
    stats = parser.get_fastqc_summary()
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
        with file_transaction(out_file) as tmp_out_file:
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
    return out

def _bed_to_bed6(orig_file, out_dir):
    """Convert bed to required bed6 inputs.
    """
    import pybedtools
    bed6_file = os.path.join(out_dir, "%s-bed6%s" % os.path.splitext(os.path.basename(orig_file)))
    if not utils.file_exists(bed6_file):
        with open(bed6_file, "w") as out_handle:
            for i, region in enumerate(list(x) for x in pybedtools.BedTool(orig_file)):
                fillers = [str(i), "1.0", "+"]
                full = region + fillers[:6 - len(region)]
                out_handle.write("\t".join(full) + "\n")
    return bed6_file

def _run_qualimap(bam_file, data, out_dir):
    """Run qualimap to assess alignment quality metrics.
    """
    report_file = os.path.join(out_dir, "qualimapReport.html")
    if not os.path.exists(report_file):
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
        regions = data["config"]["algorithm"].get("variant_regions")
        if regions:
            bed6_regions = _bed_to_bed6(regions, out_dir)
            cmd += " -gff {bed6_regions}"
        do.run(cmd.format(**locals()), "Qualimap: %s" % data["name"][-1])
    return _parse_qualimap_metrics(report_file)

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

def _run_bamtools_stats(bam_file, data, out_dir):
    """Run bamtools stats with reports on mapped reads, duplicates and insert sizes.
    """
    stats_file = os.path.join(out_dir, "bamtools_stats.txt")
    if not utils.file_exists(stats_file):
        utils.safe_makedir(out_dir)
        bamtools = config_utils.get_program("bamtools", data["config"])
        with file_transaction(stats_file) as tx_out_file:
            cmd = "{bamtools} stats -in {bam_file}"
            if bam.is_paired(bam_file):
                cmd += " -insert"
            cmd += " > {tx_out_file}"
            do.run(cmd.format(**locals()), "bamtools stats", data)
    return _parse_bamtools_stats(stats_file)

## Variant statistics from gemini

def _run_gemini_stats(bam_file, data, out_dir):
    """Retrieve high level variant statistics from Gemini.
    """
    out = {}
    gemini_db = (data.get("variants", [{}])[0].get("population", {}).get("db") 
                 if data.get("variants") else None)
    if gemini_db:
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
                if len(parts) > 0 and parts[0] == data["name"][-1]:
                    _, hom_ref, het, hom_var, _, total = parts
                    out["Variations (total)"] = int(total)
                    out["Variations (heterozygous)"] = int(het)
                    out["Variations (homozygous)"] = int(hom_var)
                    break
            out["Variations (in dbSNP)"] = int(dbsnp_count.strip())
            if out.get("Variations (total)") > 0:
                out["Variations (in dbSNP) pct"] = "%.1f%%" % (out["Variations (in dbSNP)"] /
                                                               float(out["Variations (total)"]) * 100.0)
            with open(gemini_stat_file, "w") as out_handle:
                yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
        else:
            with open(gemini_stat_file) as in_handle:
                out = yaml.safe_load(in_handle)
    return out


## qsignature

def _run_qsignature_generator(bam_file,data,out_dir):
    """ Run SignatureGenerator to create normalize vcf
    that later will be input of qsignature_summary

    :param bam_file: (str) path of the bam_file
    :param data: (list) list containing the all the dictionary
                     for this sample
    :param out_dir: (str) path of the output
    
    :returns: (dict) dict with the normalize vcf file 

    """
    position = dd.get_qsig_file(data)
    if position:
        slice_bam = _slice_chr22(bam_file, data)
        resources = config_utils.get_resources("qsignature", data["config"])
        qsig = config_utils.get_program("qsignature", data["config"])
        jvm_opts = "-Xms750m -Xmx8g"
        cores = resources.get("cores", 1)
        utils.safe_makedir(out_dir)
        out_name = os.path.basename(slice_bam).replace("bam","qsig.vcf")
        out_file = os.path.join(out_dir, out_name)
        log_file = os.path.join(out_dir, "qsig.log")
        base_cmd = ("{qsig} {jvm_opts} "
                    "org.qcmg.sig.SignatureGenerator "
                    "--noOfThreads {cores} "
                    "-log {log_file} -i {position} "
                    "-i {down_file} ")
        if not os.path.exists(out_file):    
            down_file = bam.downsample(slice_bam, data, 20000000)
            if not down_file:
                down_file = slice_bam
            file_qsign_out = down_file.replace("bam","bam.qsig.vcf")
            do.run(base_cmd.format(**locals()),"qsignature 1: %s" % data["name"][-1])
            if os.path.exists(file_qsign_out):
                with file_transaction(out_file) as file_txt_out:
                    shutil.move(file_qsign_out,file_txt_out)
            else:
    			raise IOError("File doesn't exist %s" % file_qsign_out)
        return {'qsig_vcf':out_file}
    else:
        logger.info("There is no qsignature for this species: %s" 
            % ['config']['algorithm']['genome_build'])

def qsignature_summary(*samples):
    """ Run SignatureCompareRelatedSimple module from 
    qsignature tool to creata a matrix of pairwise 
    comparison among samples. The function
    will not run if the output exisits

    :param samples: list with only one element containing 
        all samples information
    :returns: (dict) with the path of the output to be joined to summary

    """ 
    count = 0
    warnings = []
    qsig = config_utils.get_program("qsignature", samples[0][0]["config"])
    jvm_opts = "-Xms750m -Xmx8g"
    out_dir = utils.safe_makedir(os.path.join(samples[0][0]["dirs"]["work"],"qsignature"))
    log = os.path.join(samples[0][0]["dirs"]["work"], "qsig.log")
    out_file = os.path.join(samples[0][0]["dirs"]["work"], "qc","qsignature.xml")
    out_ma_file = os.path.join(samples[0][0]["dirs"]["work"], "qc","qsignature.ma")
    out_warn_file = os.path.join(samples[0][0]["dirs"]["work"], "qc","qsignature.warnings")
    for data in samples:
        data = data[0]
        if data['summary']['metrics'].get('qsig_vcf',False):
            count += 1
            vcf = data['summary']['metrics']['qsig_vcf']
            vcf_name = os.path.basename(vcf)
            if not os.path.lexists(os.path.join(out_dir,vcf_name)):
                os.symlink(vcf,os.path.join(out_dir,vcf_name))
    if count > 0:
        if not os.path.exists(out_file):
            with file_transaction(out_file) as file_txt_out:
                base_cmd = ("{qsig} {jvm_opts} "
    	                    "org.qcmg.sig.SignatureCompareRelatedSimple "
    	                    "-log {log} -dir {out_dir} "
    	                    "-o {file_txt_out} ")
                do.run(base_cmd.format(**locals()),"qsignature 2")
        warnings = _parse_qsignature_output(out_file,out_ma_file,out_warn_file)
        return [{'qsig_matrix': out_ma_file,
                 'qsig_warnings': out_warn_file,
                 'warnings samples': list(warnings)}]
 
def _parse_qsignature_output(in_file,out_file,warning_file):
    """ Parse xml file produced by qsignature

    :param in_file: (str) with the path to the xml file
    :param out_file: (str) with the path to output file
    :param warning_file: (str) with the path to warning file

    :returns: (list) with samples that could be duplicated

    """ 
    name = {}
    score = {}
    warnings = set()
    with open(in_file,'r') as in_handle:
        with file_transaction(out_file) as out_tx_file:
            with file_transaction(warning_file) as warn_tx_file:
                with open(out_tx_file,'w') as out_handle:
                    with open(warn_tx_file,'w') as warn_handle:
                        ET = lxml.etree.parse(in_handle)
                        for i in list(ET.iter('file')):
                            name[i.attrib['id']] = os.path.basename(i.attrib['name']).replace(".bam.qsig.vcf","")
                        for i in list(ET.iter('comparison')):
                            out_handle.write("%s\t%s\t%s\n" % 
                            (name[i.attrib['file1']],name[i.attrib['file2']],i.attrib['score']))
                            if float(i.attrib['score']) < 0.1:
                                logger.info('qsignature WARNING: risk of duplicated samples:%s' %
                                    (' '.join([name[i.attrib['file1']],name[i.attrib['file2']]])))
                                warn_handle.write('qsignature WARNING: risk of duplicated samples:%s\n' %
                                    (' '.join([name[i.attrib['file1']],name[i.attrib['file2']]])))
                                warnings.add(name[i.attrib['file1']])
                                warnings.add(name[i.attrib['file2']])
        return warnings


def _slice_chr22(in_bam, data):
    """
    return only one BAM file with only chromosome 22
    """
    sambamba = config_utils.get_program("sambamba", data["config"])
    out_file = "%s-chr%s" % os.path.splitext(in_bam)
    if not utils.file_exists(out_file):
        with contextlib.closing(pysam.Samfile(in_bam, "rb")) as bamfile:
            bam_contigs = [c["SN"] for c in bamfile.header["SQ"]]
        chromosome = "22"
        if "chr22" in bam_contigs:
            chromosome = "chr22"
        with file_transaction(out_file) as tx_out_file:
            cmd = ("{sambamba} slice -o {tx_out_file} {in_bam} {chromosome}" ).format(**locals())
            out = subprocess.check_output(cmd, shell=True)
    return out_file
