#!/usr/bin/env python
"""Generate alignment plots and tables for a sequencing run and alignment.

Requires:
- Picard reporting tools
- fastx toolkit
- SolexaQA
- LaTeX -- pdflatex
- ps2pdf

Generate a LaTeX report on a sequencing run and it's alignment to the 
reference genome.

Usage:
    align_summary_report.py <picard dir> <input bam> <ref file> <fastq one> <fastq two>
    --bait=<bait file>
    --target=<target file>
    --config=<YAML configuration file specifying executable information>
    --sort -- Sort the file first
"""
import os
import sys
import glob
import subprocess
import shutil
import tempfile
from optparse import OptionParser

import yaml
from mako.template import Template

from bcbio.picard import PicardRunner
from bcbio.picard.metrics import PicardMetrics

PARAM_DEFAULT = dict(
        fastx_stats = "fastx_quality_stats",
        fastx_graph = "fastq_quality_boxplot_graph.sh",
        pdflatex = "pdflatex",
        ps2pdf = "ps2pdf",
        solexaqa = "SolexaQA.pl",
        )

def main(picard_dir, align_bam, ref_file, fastq_one, fastq_pair=None,
        bait_file=None, target_file=None, do_sort=False, sample_name="",
        config=None):
    tmp_dir = _make_tmpdir()
    work_dir = os.getcwd()
    if config:
        with open(config) as in_handle:
            params = yaml.load(in_handle)["program"]
    else:
        params = PARAM_DEFAULTS
    picard = PicardRunner(picard_dir)
    if do_sort:
        align_bam = picard_sort(picard, align_bam, tmp_dir)

    metrics = PicardMetrics(picard, tmp_dir)
    summary_table, metrics_graphs = metrics.report(
            align_bam, ref_file, fastq_pair is not None,
            bait_file, target_file)
    base, ext = os.path.splitext(align_bam)
    base = base.replace(".", "-")
    total_count, read_size, fastq_graphs = plot_fastq_stats(
            [fastq_one, fastq_pair], base, params)
    qa_graphs = solexaqa_plots([fastq_one, fastq_pair], params, work_dir)

    # add read_size to the total summary table
    summary_table[0] = (summary_table[0][0], summary_table[0][1],
            "%sbp %s" % (read_size, summary_table[0][-1]))
    ref_org = os.path.splitext(os.path.split(ref_file)[-1])[0]
    summary_table.insert(0, ("Reference organism",
        ref_org.replace("_", " "), ""))
    tmpl = Template(section_template)
    sample_name = "%s (%s)" % (sample_name.replace("_", "\_"),
            base.replace("_", " "))
    section = tmpl.render(name=sample_name, summary=None,
            summary_table=summary_table,
            figures=[(f, c) for (f, c) in metrics_graphs + fastq_graphs +
                     qa_graphs if f],
            recal_figures=_get_recal_plots(work_dir, align_bam))
    out_file = "%s-summary.tex" % base
    out_tmpl = Template(base_template)
    with open(out_file, "w") as out_handle:
        out_handle.write(out_tmpl.render(parts=[section]))
    run_pdflatex(out_file, params)
    shutil.rmtree(tmp_dir)

def plot_fastq_stats(fastq_files, out_base, params):
    """Use fastx toolkit to prepare a plot of quality statistics.
    """
    print "Drawing plot of fastq quality scores"
    fastq_files = [f for f in fastq_files if f]
    graphs = []
    for i, fastq_file in enumerate(fastq_files):
        if len(fastq_files) > 1:
            cur_out_base = "%s_%s" % (out_base, i+1)
        else:
            cur_out_base = out_base
        stat_file = "%s_fastq_stats.txt" % cur_out_base
        if not os.path.exists(stat_file):
            stats_cl = [params["fastx_stats"], "-i", fastq_file,
                        "-o", stat_file]
            child = subprocess.Popen(stats_cl)
            child.wait()
        with open(stat_file) as in_handle:
            # look at the last line in the stats file
            for line in in_handle:
                pass
            parts = line.split()
            count = int(parts[1])
            read_size = int(parts[0])
        graph_file = "%s_fastq_qual.pdf" % cur_out_base
        if not os.path.exists(graph_file):
            ps_file = "%s.ps" % os.path.splitext(graph_file)[0]
            boxplot_cl = [params["fastx_graph"], "-i", stat_file,
                         "-p",
                         "-t", fastq_file,
                         "-o", ps_file]
            child = subprocess.Popen(boxplot_cl)
            child.wait()
            topdf_cl = [params["ps2pdf"], ps_file, graph_file]
            child = subprocess.Popen(topdf_cl)
            child.wait()
        graphs.append((graph_file,
            "Quality score distribution per base for read %s" % (i + 1)))
    return count, read_size, graphs

def solexaqa_plots(fastq_files, params, work_dir):
    print "SolexaQA plots of fastq error distributions"
    graphs = []
    for i, fastq_file in enumerate(f for f in fastq_files if f):
        orig_tile_graph, tile_graph = _sqa_file(fastq_file, "png", work_dir)
        orig_qual_graph, qual_graph  = _sqa_file(fastq_file, "quality.pdf",
                                                   work_dir)
        if not os.path.exists(tile_graph) or not os.path.exists(qual_graph):
            cl = [params["solexaqa"], fastq_file]
            subprocess.check_call(cl)
            os.rename(orig_tile_graph, tile_graph)
            os.rename(orig_qual_graph, qual_graph)
        #graphs.append((tile_graph,
        #    "Error distribution per position and tile for read %s. "\
        #    "Darker squares correspond to poor quality scores." % (i + 1)))
        graphs.append((qual_graph,
            "Mean error probability per read position and tile for read %s. "\
            "Ideal flowcells will have a tight range of values "\
            "for all tiles." % (i + 1)))
    return graphs

def _sqa_file(fname, ext, work_dir):
    """Naming for SolexaQA.pl output. Handles local files and fixing extensions
    """
    orig_file = os.path.join(work_dir, "%s.%s" % (os.path.basename(fname), ext))
    # replace all '.' extensions except for last so pdflatex can handle the file
    safe_file = orig_file.replace(".", "_", orig_file.count(".") - 1)
    return orig_file, safe_file

def picard_sort(picard, align_bam, tmp_dir):
    base, ext = os.path.splitext(align_bam)
    out_file = "%s-sort%s" % (base, ext)
    if not os.path.exists(out_file):
        opts = [("INPUT", align_bam),
                ("OUTPUT", out_file),
                ("TMP_DIR", tmp_dir),
                ("SORT_ORDER", "coordinate")]
        picard.run("SortSam", opts)
    return out_file

def run_pdflatex(tex_file, params):
    cl = [params["pdflatex"], tex_file]
    with open(os.devnull, "w") as ignore:
        child = subprocess.Popen(cl, stdout=ignore, stderr=ignore)
    child.wait()

def _get_recal_plots(work_dir, align_bam):
    """Retrieve any recalibration report plots for display.
    """
    (base, _) = os.path.splitext(align_bam)
    reports = glob.glob(os.path.join(work_dir, "reports", "images",
        "%s*-plot.pdf" % base))
    reports.sort()
    return reports

def _make_tmpdir():
    tmp_dir_base = os.path.join(os.getcwd(), "tmp")
    if not os.path.exists(tmp_dir_base):
        os.makedirs(tmp_dir_base)
    tmp_dir = tempfile.mkdtemp(dir=tmp_dir_base)
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    return tmp_dir

section_template = r"""
\subsection*{${name}}

% if summary_table:
    \begin{table}[h]
    \centering
    \begin{tabular}{|l|rr|}
    \hline
    % for label, val, extra in summary_table:
        %if label is not None:
            ${label} & ${val} & ${extra} \\ 
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

% for i, (figure, caption) in enumerate(figures):
    \begin{figure}[htbp]
      \centering
      \includegraphics[width=
      % if i < 1:
        0.52
      % else:
        0.75
      %endif
      \linewidth]{${figure}}
      \caption{${caption}}
    \end{figure}
% endfor

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

base_template = r"""
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

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-b", "--bait", dest="bait_file")
    parser.add_option("-t", "--target", dest="target_file")
    parser.add_option("-n", "--name", dest="sample_name")
    parser.add_option("-s", "--sort", dest="do_sort", action="store_true",
            default=False)
    parser.add_option("-c", "--config", dest="config")
    (options, args) = parser.parse_args()
    if len(args) < 4:
        print __doc__
        sys.exit()
    kwargs = dict(
        bait_file=options.bait_file,
        target_file=options.target_file,
        sample_name=options.sample_name,
        config=options.config,
        do_sort=options.do_sort)
    main(*args, **kwargs)
