"""Quality control and summary metrics for next-gen alignments and analysis.
"""
import os
import shutil
import subprocess

import lxml.html
import pybedtools
import yaml
from datetime import datetime

# allow graceful during upgrades
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import config_utils
from bcbio.provenance import do
import bcbio.rnaseq.qc
from bcbio.variation.realign import has_aligned_reads

# ## High level functions to generate summary PDF

def generate_parallel(samples, run_parallel):
    """Provide parallel preparation of summary information for alignment and variant calling.
    """
    sum_samples = run_parallel("pipeline_summary", samples)
    summary_file = write_project_summary(sum_samples)
    samples = []
    for data in sum_samples:
        if "summary" not in data[0]:
            data[0]["summary"] = {}
        data[0]["summary"]["project"] = summary_file
        samples.append(data)
    return samples

def pipeline_summary(data):
    """Provide summary information on processing sample.
    """
    work_bam = (data.get("work_bam")
                if data["config"]["algorithm"].get("merge_bamprep", True)
                else data.get("callable_bam"))
    if data["sam_ref"] is not None and work_bam and has_aligned_reads(work_bam):
        logger.info("Generating summary files: %s" % str(data["name"]))
        data["summary"] = _run_qc_tools(work_bam, data)
    return [[data]]

def _run_qc_tools(bam_file, data):
    """Run a set of third party quality control tools, returning QC directory and metrics.
    """
    to_run = [("fastqc", _run_fastqc)]
    if data["analysis"].lower() == "rna-seq":
        to_run.append(("rnaseqc", bcbio.rnaseq.qc.sample_summary))
        to_run.append(("complexity", _run_complexity))
    elif data["analysis"].lower() == "chip-seq":
        to_run.append(["bamtools", _run_bamtools_stats])
    else:
        to_run += [("bamtools", _run_bamtools_stats), ("gemini", _run_gemini_stats)]
    qc_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "qc", data["name"][-1]))
    metrics = {}
    for program_name, qc_fn in to_run:
        cur_qc_dir = os.path.join(qc_dir, program_name)
        cur_metrics = qc_fn(bam_file, data, cur_qc_dir)
        metrics.update(cur_metrics)
    metrics["Name"] = data["name"][-1]
    metrics["Quality format"] = utils.get_in(data,
                                             ("config", "algorithm",
                                              "quality_format"),
                                             "standard").lower()
    return {"qc": qc_dir, "metrics": metrics}

# ## Generate project level QC summary for quickly assessing large projects

def write_project_summary(samples):
    """Write project summary information on the provided samples.
    write out dirs, genome resources,

    """
    work_dir = samples[0][0]["dirs"]["work"]
    out_file = os.path.join(work_dir, "project-summary.yaml")
    upload_dir = (os.path.join(work_dir, samples[0][0]["upload"]["dir"])
                  if "dir" in samples[0][0]["upload"] else "")
    test_run = samples[0][0].get("test_run", False)
    date = str(datetime.now())
    with open(out_file, "w") as out_handle:
        yaml.dump({"date": date}, out_handle,
                  default_flow_style=False, allow_unicode=False)
        if test_run:
            yaml.dump({"test_run": True}, out_handle, default_flow_style=False,
                      allow_unicode=False)
        yaml.dump({"upload": upload_dir}, out_handle,
                  default_flow_style=False, allow_unicode=False)
        yaml.dump({"bcbio_system": samples[0][0]["config"].get("bcbio_system", "")}, out_handle,
                  default_flow_style=False, allow_unicode=False)
        yaml.dump({"samples": [_save_fields(sample[0]) for sample in samples]}, out_handle,
                  default_flow_style=False, allow_unicode=False)
    return out_file

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
                saved["summary"]["metrics"]["Disambiguated %s reads" % str(sample["algorithm"]["disambiguate"][0])] = disambigStats[1]
                saved["summary"]["metrics"]["Disambiguated ambiguous reads"] = disambigStats[2]
    return saved

def _parse_disambiguate(disambiguatestatsfilename):
    """Parse disambiguation stats from given file.
    """
    disambigStats = [-1 -1 -1]
    with open(disambiguatestatsfilename, "r") as in_handle:
        header = in_handle.readline().strip().split("\t")
        if header == ['sample','unique species A pairs','unique species B pairs','ambiguous pairs']:
            disambigStatsTemp = in_handle.readline().strip().split("\t")[1:]
            if len(disambigStatsTemp) == 3:
                disambigStats = [int(x) for x in disambigStatsTemp]
    return disambigStats


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

def _run_fastqc(bam_file, data, fastqc_out):
    """Run fastqc, generating report in specified directory and parsing metrics.

    Downsamples to 10 million reads to avoid excessive processing times with large
    files.
    """
    if not os.path.exists(os.path.join(fastqc_out, "fastqc_report.html")):
        work_dir = os.path.dirname(fastqc_out)
        utils.safe_makedir(work_dir)
        ds_bam = bam.downsample(bam_file, data, 1e7)
        num_cores = data["config"]["algorithm"].get("num_cores", 1)
        bam_file = ds_bam if ds_bam else bam_file
        with utils.curdir_tmpdir(work_dir) as tx_tmp_dir:
            with utils.chdir(tx_tmp_dir):
                cl = [config_utils.get_program("fastqc", data["config"]),
                      "-t", str(num_cores), "-o", tx_tmp_dir, "-f", "bam", bam_file]
                do.run(cl, "FastQC: %s" % data["name"][-1])
                fastqc_outdir = os.path.join(tx_tmp_dir,
                                             "%s_fastqc" % os.path.splitext(os.path.basename(bam_file))[0])
                if os.path.exists("%s.zip" % fastqc_outdir):
                    os.remove("%s.zip" % fastqc_outdir)
                shutil.move(fastqc_outdir, fastqc_out)
        if ds_bam and os.path.exists(ds_bam):
            os.remove(ds_bam)
    parser = FastQCParser(fastqc_out)
    stats = parser.get_fastqc_summary()
    return stats

def _run_complexity(bam_file, data, out_dir):
    base, _ = os.path.splitext(os.path.basename(bam_file))
    utils.safe_makedir(out_dir)
    out_file = os.path.join(out_dir, base + ".pdf")
    df = bcbio.rnaseq.qc.starts_by_depth(bam_file, data["config"])
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
    parsers = {"Globals" : _parse_qualimap_globals,
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
    bed6_file = os.path.join(out_dir, "%s-bed6%s" % os.path.splitext(os.path.basename(orig_file)))
    if not utils.file_exists(bed6_file):
        with open(bed6_file, "w") as out_handle:
            for i, region in enumerate(list(x) for x in pybedtools.BedTool(orig_file)):
                fillers = [str(i), "1.0", "+"]
                full = region + fillers[:6-len(region)]
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
    gemini_db = data.get("variants", [{}])[0].get("population", {}).get("db")
    if gemini_db:
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
    return out
