"""Quality control and summary metrics for next-gen alignments and analysis.
"""
import contextlib
import copy
import csv
import glob
import os
import subprocess
import xml.etree.ElementTree as ET

import yaml
from mako.template import Template
import pysam

from bcbio import utils
from bcbio.broad import runner_from_config
from bcbio.broad.metrics import PicardMetrics, PicardMetricsParser
from bcbio.log import logger
from bcbio.pipeline import config_utils

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
    if data["sam_ref"] is not None and data["work_bam"]:
        logger.info("Generating summary files: %s" % str(data["name"]))
        data["summary"] = generate_align_summary(data["work_bam"],
                                                 data["sam_ref"], data["name"],
                                                 data["config"], data["dirs"])
    return [[data]]

def generate_align_summary(bam_file, sam_ref, sample_name,
                           config, dirs):
    """Run alignment summarizing script to produce a pdf with align details.
    """
    with utils.chdir(dirs["work"]):
        with utils.curdir_tmpdir() as tmp_dir:
            graphs, summary, overrep = \
                    _graphs_and_summary(bam_file, sam_ref, tmp_dir, config)
        return {"pdf": _generate_pdf(graphs, summary, overrep, bam_file, sample_name,
                                     dirs, config)}

def _safe_latex(to_fix):
    """Escape characters that make LaTeX unhappy.
    """
    chars = ["%", "_", "&", "#"]
    for char in chars:
        to_fix = to_fix.replace(char, "\\%s" % char)
    return to_fix

def _generate_pdf(graphs, summary, overrep, bam_file, sample_name,
                  dirs, config):
    base = os.path.splitext(os.path.basename(bam_file))[0]
    sample_name = base if sample_name is None else " : ".join(sample_name)
    tmpl = Template(_section_template)
    sample_name = "%s (%s)" % (_safe_latex(sample_name),
                               _safe_latex(base))
    recal_plots = sorted(glob.glob(os.path.join(dirs["work"], "reports", "images",
                                                "%s*-plot.pdf" % base)))
    section = tmpl.render(name=sample_name, summary=None,
                          summary_table=summary,
                          figures=[(f, c, i) for (f, c, i) in graphs if f],
                          overrep=overrep,
                          recal_figures=recal_plots)
    out_file = os.path.join(dirs["work"], "%s-summary.tex" % base)
    out_tmpl = Template(_base_template)
    with open(out_file, "w") as out_handle:
        out_handle.write(out_tmpl.render(parts=[section]))
    pdf_file = "%s.pdf" % os.path.splitext(out_file)[0]
    if (config["algorithm"].get("write_summary", True) and
         not utils.file_exists(pdf_file)):
        cl = [config_utils.get_program("pdflatex", config), out_file]
        subprocess.check_call(cl)
    return pdf_file

def _graphs_and_summary(bam_file, sam_ref, tmp_dir, config):
    """Prepare picard/FastQC graphs and summary details.
    """
    bait = config["algorithm"].get("hybrid_bait", None)
    target = config["algorithm"].get("hybrid_target", None)
    if target is not None:
        assert os.path.exists(target), (target, "does not exist!")
    broad_runner = runner_from_config(config)
    metrics = PicardMetrics(broad_runner, tmp_dir)
    summary_table, metrics_graphs = \
                   metrics.report(bam_file, sam_ref, is_paired(bam_file), bait, target)
    metrics_graphs = [(p, c, 0.75) for p, c in metrics_graphs]
    fastqc_graphs, fastqc_stats, fastqc_overrep = \
                   fastqc_report(bam_file, config)
    all_graphs = fastqc_graphs + metrics_graphs
    summary_table = _update_summary_table(summary_table, sam_ref, fastqc_stats)
    return all_graphs, summary_table, fastqc_overrep

def _update_summary_table(summary_table, ref_file, fastqc_stats):
    stats_want = []
    summary_table[0] = (summary_table[0][0], summary_table[0][1],
            "%sbp %s" % (fastqc_stats.get("Sequence length", "0"), summary_table[0][-1]))
    for stat in stats_want:
        summary_table.insert(0, (stat, fastqc_stats.get(stat, ""), ""))
    ref_org = os.path.splitext(os.path.split(ref_file)[-1])[0]
    summary_table.insert(0, ("Reference organism",
        ref_org.replace("_", " "), ""))
    return summary_table

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
    with utils.curdir_tmpdir() as tmp_dir:
        for sample in (x[0] for x in samples):
            if sample["work_bam"] is not None:
                _, summary, _ = _graphs_and_summary(sample["work_bam"], sample["sam_ref"],
                                                    tmp_dir, sample["config"])
                sample_info = {}
                for xs in summary:
                    n = xs[0]
                    if n is not None:
                        sample_info[n] = xs[1:]
                sample_name = ";".join([x for x in sample["name"] if x])
                out.append((sample_name, sample_info))
    return out

def is_paired(bam_file):
    """Determine if a BAM file has paired reads.
    """
    with contextlib.closing(pysam.Samfile(bam_file, "rb")) as in_pysam:
        for read in in_pysam:
            return read.is_paired

# ## Run and parse read information from FastQC

def fastqc_report(bam_file, config):
    """Calculate statistics about a read using FastQC.
    """
    out_dir = _run_fastqc(bam_file, config)
    parser = FastQCParser(out_dir)
    graphs = parser.get_fastqc_graphs()
    stats, overrep = parser.get_fastqc_summary()
    return graphs, stats, overrep

class FastQCParser:
    def __init__(self, base_dir):
        self._dir = base_dir
        self._max_seq_size = 45
        self._max_overrep = 20

    def get_fastqc_graphs(self):
        graphs = (("per_base_quality.png", "", 1.0),
                  ("per_base_sequence_content.png", "", 0.85),
                  ("per_sequence_gc_content.png", "", 0.85),
                  ("kmer_profiles.png", "", 0.85),)
        final_graphs = []
        for f, caption, size in graphs:
            full_f = os.path.join(self._dir, "Images", f)
            if os.path.exists(full_f):
                final_graphs.append((full_f, caption, size))
        return final_graphs

    def get_fastqc_summary(self):
        stats = {}
        for stat_line in self._fastqc_data_section("Basic Statistics")[1:]:
            k, v = [_safe_latex(x) for x in stat_line.split("\t")[:2]]
            stats[k] = v
        over_rep = []
        for line in self._fastqc_data_section("Overrepresented sequences")[1:]:
            parts = [_safe_latex(x) for x in line.split("\t")]
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


def _run_fastqc(bam_file, config):
    out_base = "fastqc"
    utils.safe_makedir(out_base)
    fastqc_out = os.path.join(out_base, "%s_fastqc" %
                              os.path.splitext(os.path.basename(bam_file))[0])
    if not os.path.exists(fastqc_out):
        cl = [config_utils.get_program("fastqc", config),
              "-o", out_base, "-f", "bam", bam_file]
        subprocess.check_call(cl)
    if os.path.exists("%s.zip" % fastqc_out):
        os.remove("%s.zip" % fastqc_out)
    return fastqc_out


# ## High level summary in YAML format for loading into Galaxy.

def write_metrics(run_info, fc_name, fc_date, dirs):
    """Write an output YAML file containing high level sequencing metrics.
    """
    lane_stats, sample_stats, tab_metrics = summary_metrics(run_info,
            dirs["work"], fc_name, fc_date, dirs["fastq"])
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

def summary_metrics(run_info, analysis_dir, fc_name, fc_date, fastq_dir):
    """Reformat run and analysis statistics into a YAML-ready format.
    """
    tab_out = []
    lane_info = []
    sample_info = []
    for lane_xs in run_info["details"]:
        run = lane_xs[0]
        tab_out.append([run["lane"], run.get("researcher", ""),
            run.get("name", ""), run.get("description")])
        base_info = dict(
                researcher = run.get("researcher_id", ""),
                sample = run.get("sample_id", ""),
                lane = run["lane"],
                request = run_info["run_id"])
        cur_lane_info = copy.deepcopy(base_info)
        cur_lane_info["metrics"] = _bustard_stats(run["lane"], fastq_dir,
                                                  fc_date, analysis_dir)
        lane_info.append(cur_lane_info)
        for lane_x in lane_xs:
            cur_name = "%s_%s_%s" % (run["lane"], fc_date, fc_name)
            if lane_x["barcode_id"]:
                cur_name = "%s_%s-" % (cur_name, lane_x["barcode_id"])
            stats = _metrics_from_stats(_lane_stats(cur_name, analysis_dir))
            if stats:
                cur_run_info = copy.deepcopy(base_info)
                cur_run_info["metrics"] = stats
                cur_run_info["barcode_id"] = str(lane_x["barcode_id"])
                cur_run_info["barcode_type"] = str(lane_x.get("barcode_type", ""))
                sample_info.append(cur_run_info)
    return lane_info, sample_info, tab_out

def _metrics_from_stats(stats):
    """Remap Broad metrics names to our local names.
    """
    if stats:
        s_to_m = dict(
                AL_MEAN_READ_LENGTH = 'Read length',
                AL_TOTAL_READS = 'Reads',
                AL_PF_READS_ALIGNED = 'Aligned',
                DUP_READ_PAIR_DUPLICATES = 'Pair duplicates'
                )
        metrics = dict()
        for stat_name, metric_name in s_to_m.iteritems():
            metrics[metric_name] = stats.get(stat_name, 0)
        return metrics

def _bustard_stats(lane_num, fastq_dir, fc_date, analysis_dir):
    """Extract statistics about the flow cell from Bustard outputs.
    """
    stats = dict()
    if fastq_dir:
        sum_file = os.path.join(fastq_dir, os.pardir, "BustardSummary.xml")
        #sum_file = os.path.join(fc_dir, "Data", "Intensities", "BaseCalls",
        #        "BustardSummary.xml")
        if os.path.exists(sum_file):
            with open(sum_file) as in_handle:
                results = ET.parse(in_handle).getroot().find("TileResultsByLane")
                for lane in results:
                    if lane.find("laneNumber").text == str(lane_num):
                        stats = _collect_cluster_stats(lane)
    read_stats = _calc_fastq_stats(analysis_dir, lane_num, fc_date)
    stats.update(read_stats)
    return stats

def _calc_fastq_stats(analysis_dir, lane_num, fc_date):
    """Grab read length from fastq; could provide distribution if non-equal.
    """
    stats = dict()
    fastqc_dirs = glob.glob(os.path.join(analysis_dir, "fastqc",
                                         "%s_%s*" % (lane_num, fc_date)))
    if len(fastqc_dirs) > 0:
        parser = FastQCParser(sorted(fastqc_dirs)[-1])
        fastqc_stats, _ = parser.get_fastqc_summary()
        stats["Read length"] = fastqc_stats["Sequence length"]
    return stats

def _collect_cluster_stats(lane):
    """Retrieve total counts on cluster statistics.
    """
    stats = {"Clusters" : 0, "Clusters passed": 0}
    for tile in lane.find("Read").findall("Tile"):
        stats["Clusters"] += int(tile.find("clusterCountRaw").text)
        stats["Clusters passed"] += int(tile.find("clusterCountPF").text)
    return stats

def _lane_stats(cur_name, work_dir):
    """Parse metrics information from files in the working directory.
    """
    parser = PicardMetricsParser()
    metrics_files = glob.glob(os.path.join(work_dir, "%s*metrics" % cur_name))
    metrics = parser.extract_metrics(metrics_files)
    return metrics

# ## LaTeX templates for output PDF

_section_template = r"""
\subsection*{${name}}

% if summary_table:
    \begin{table}[h]
    \centering
    \begin{tabular}{|l|rr|}
    \hline
    % for label, val, extra in summary_table:
        %if label is not None:
            ${label} & ${val} & ${extra} \\%
        %else:
            \hline
        %endif
    %endfor
    \hline
    \end{tabular}
    \caption{Summary of lane results}
    \end{table}
% endif

% if summary:
    \begin{verbatim}
    ${summary}
    \end{verbatim}
% endif

% for i, (figure, caption, size) in enumerate(figures):
    \begin{figure}[htbp]
      \centering
      \includegraphics[width=${size}\linewidth] {${figure}}
      \caption{${caption}}
    \end{figure}
% endfor

% if len(overrep) > 0:
    \begin{table}[htbp]
    \centering
    \begin{tabular}{|p{8cm}rrp{4cm}|}
    \hline
    Sequence & Count & Percent & Match \\%
    \hline
    % for seq, count, percent, match in overrep:
        \texttt{${seq}} & ${count} & ${"%.2f" % float(percent)} & ${match} \\%
    % endfor
    \hline
    \end{tabular}
    \caption{Overrepresented read sequences}
    \end{table}
% endif

\FloatBarrier
% if len(recal_figures) > 0:
    \subsubsection*{Quality score recalibration}
    % for figure in recal_figures:
        \begin{figure}[htbp]
          \centering
          \includegraphics[width=0.48\linewidth]{${figure}}
        \end{figure}
    % endfor
% endif
\FloatBarrier
"""

_base_template = r"""
\documentclass{article}
\usepackage{fullpage}
\usepackage{graphicx}
\usepackage{placeins}

\begin{document}
% for part in parts:
    ${part}
% endfor
\end{document}
"""
