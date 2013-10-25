"""Quality control and summary metrics for next-gen alignments and analysis.
"""
import copy
import csv
import os
import shutil

import yaml

from bcbio import utils
from bcbio.log import logger
from bcbio.pipeline import config_utils
from bcbio.provenance import do
import bcbio.rnaseq.qc

# ## High level functions to generate summary PDF

def generate_parallel(samples, run_parallel):
    """Provide parallel preparation of summary information for alignment and variant calling.
    """
    need_summary = False
    for data in samples:
        if data[0]["config"]["algorithm"].get("write_summary", True):
            need_summary = True
    if need_summary:
        sum_samples = run_parallel("pipeline_summary", samples)
        summary_csv = write_project_summary(sum_samples)
        samples = []
        for data in sum_samples:
            if summary_csv:
                if "summary" not in data[0]:
                    data[0]["summary"] = {}
                data[0]["summary"]["project"] = summary_csv
            samples.append(data)
    return samples

def pipeline_summary(data):
    """Provide summary information on processing sample.
    """
    work_bam = (data.get("work_bam")
                if data["config"]["algorithm"].get("merge_bamprep", True)
                else data.get("callable_bam"))
    if data["sam_ref"] is not None and work_bam:
        logger.info("Generating summary files: %s" % str(data["name"]))
        data["summary"] = _run_qc_tools(work_bam, data)
    return [[data]]

def _run_qc_tools(bam_file, data):
    """Run a set of third party quality control tools, returning QC directory and metrics.
    """
    to_run = [("fastqc", _run_fastqc)]
    if data["analysis"].lower() == "rna-seq":
        to_run.append(("rnaseqc", bcbio.rnaseq.qc.sample_summary))
    else:
        to_run.append(("qualimap", _run_qualimap))
    qc_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "qc", data["name"][-1]))
    metrics = {}
    for program_name, qc_fn in to_run:
        cur_qc_dir = os.path.join(qc_dir, program_name)
        cur_metrics = qc_fn(bam_file, data, cur_qc_dir)
        metrics.update(cur_metrics)
    return {"qc": qc_dir, "metrics": metrics}

# ## Generate project level QC summary for quickly assessing large projects

def write_project_summary(samples):
    """Write project summary information on the provided samples.
    """
    def _nocommas(x):
        return x.replace(",", "")
    def _percent(x):
        return x.replace("(", "").replace(")", "").replace("\\", "")
    if len(samples) > 0:
        out_file = os.path.join(samples[0][0]["dirs"]["work"], "project-summary.csv")
        sample_info = _get_sample_summaries(samples)
        header = ["Total", "Aligned", "Pair duplicates", "Insert size",
                  "On target bases", "Mean target coverage", "10x coverage targets",
                  "Zero coverage targets", "Total variations", "In dbSNP",
                  "Transition/Transversion (all)", "Transition/Transversion (dbSNP)",
                  "Transition/Transversion (novel)"]
        select = [(0, _nocommas), (1, _percent), (1, _percent), (0, None),
                  (1, _percent), (0, None), (0, _percent),
                  (0, _percent), (0, None), (0, _percent),
                  (0, None), (0, None), (0, None)]
        rows = [["Sample"] + header]
        for name, info in sample_info:
            cur = [name]
            for col, (i, prep_fn) in zip(header, select):
                val = info.get(col, ["", ""])[i]
                if prep_fn and val:
                    val = prep_fn(val)
                cur.append(val)
            rows.append(cur)
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle)
            for row in rows:
                writer.writerow(row)
        return out_file

def _get_sample_summaries(samples):
    """Retrieve high level summary information for each sample.
    """
    out = []
    for sample in (x[0] for x in samples):
        summary = sample.get("summary", {}).get("metrics")
        if summary:
            sample_info = {}
            for xs in summary:
                n = xs[0]
                if n is not None:
                    sample_info[n] = xs[1:]
            sample_name = ";".join([x for x in sample["name"] if x])
            out.append((sample_name, sample_info))
    return out

# ## Run and parse read information from FastQC

class FastQCParser:
    def __init__(self, base_dir):
        self._dir = base_dir
        self._max_seq_size = 45
        self._max_overrep = 20

    def get_fastqc_summary(self):
        stats = {}
        for stat_line in self._fastqc_data_section("Basic Statistics")[1:]:
            k, v = stat_line.split("\t")[:2]
            stats[k] = v
        over_rep = []
        for line in self._fastqc_data_section("Overrepresented sequences")[1:]:
            parts = line.split("\t")
            over_rep.append(parts)
            over_rep[-1][0] = self._splitseq(over_rep[-1][0])
        return stats, over_rep[:self._max_overrep]

    def _splitseq(self, seq):
        pieces = []
        cur_piece = []
        for s in seq:
            if len(cur_piece) >= self._max_seq_size:
                pieces.append("".join(cur_piece))
                cur_piece = []
            cur_piece.append(s)
        pieces.append("".join(cur_piece))
        return " ".join(pieces)

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
    """
    if not os.path.exists(os.path.join(fastqc_out, "fastqc_report.html")):
        work_dir = os.path.dirname(fastqc_out)
        utils.safe_makedir(work_dir)
        num_cores = data["config"]["algorithm"].get("num_cores", 1)
        cl = [config_utils.get_program("fastqc", data["config"]),
              "-t", str(num_cores), "-o", work_dir, "-f", "bam", bam_file]
        do.run(cl, "FastQC: %s" % data["name"][-1])
        fastqc_outdir = os.path.join(work_dir, "%s_fastqc" % os.path.splitext(os.path.basename(bam_file))[0])
        if os.path.exists("%s.zip" % fastqc_outdir):
            os.remove("%s.zip" % fastqc_outdir)
        shutil.move(fastqc_outdir, fastqc_out)
    parser = FastQCParser(fastqc_out)
    stats, _ = parser.get_fastqc_summary()
    return stats

# ## Qualimap

def _run_qualimap(bam_file, data, out_dir):
    """Run qualimap to assess alignment quality metrics.
    """
    if not os.path.exists(os.path.join(out_dir, "qualimapReport.html")):
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
        # Qualimap requires BED6 inputs. Potentially convert BED3 over.
        #if regions:
        #    cmd += " -gff {regions}"
        do.run(cmd.format(**locals()), "Qualimap: %s" % data["name"][-1])
    # XXX Need to parse out qualimap metrics
    return {}

# ## High level summary in YAML format for loading into Galaxy.

def write_metrics(run_info, dirs):
    """Write an output YAML file containing high level sequencing metrics.
    """
    lane_stats, sample_stats, tab_metrics = summary_metrics(run_info,
            dirs["work"], dirs["fastq"])
    out_file = os.path.join(dirs["work"], "run_summary.yaml")
    with open(out_file, "w") as out_handle:
        metrics = dict(lanes=lane_stats, samples=sample_stats)
        yaml.dump(metrics, out_handle, default_flow_style=False)
    if dirs["flowcell"]:
        tab_out_file = os.path.join(dirs["flowcell"], "run_summary.tsv")
        try:
            with open(tab_out_file, "w") as out_handle:
                writer = csv.writer(out_handle, dialect="excel-tab")
                for info in tab_metrics:
                    writer.writerow(info)
        # If on NFS mounted directory can fail due to filesystem or permissions
        # errors. That's okay, we'll just not write the file.
        except IOError:
            pass
    return out_file

def summary_metrics(items, analysis_dir, fastq_dir):
    """Reformat run and analysis statistics into a YAML-ready format.

    XXX Needs a complete re-working to be less sequencer specific and
    work with more general output structure in work directory.
    """
    tab_out = []
    lane_info = []
    sample_info = []
    for run in items:
        tab_out.append([run["lane"], run.get("researcher", ""),
            run.get("name", ""), run.get("description")])
        base_info = dict(
                researcher = run.get("researcher_id", ""),
                sample = run.get("sample_id", ""),
                lane = run["lane"],
                request = run["upload"].get("run_id"))
        cur_lane_info = copy.deepcopy(base_info)
        cur_lane_info["metrics"] = {}
        lane_info.append(cur_lane_info)
    return lane_info, sample_info, tab_out
