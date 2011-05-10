#!/usr/bin/env python
"""Generate PDF report with plots and tables for a sequencing run and alignment.

Requires:
- Picard (http://picard.sourceforge.net/)
- FastQC (http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/)
- LaTeX -- pdflatex

Usage:
    align_summary_report.py <picard dir> <input bam> <ref file>
    --paired
    --bait=<bait file>
    --target=<target file>
    --config=<YAML configuration file specifying executable information>
    --sort
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

from bcbio.broad import BroadRunner
from bcbio.broad.metrics import PicardMetrics
from bcbio import utils

def main(picard_dir, align_bam, ref_file, is_paired, bait_file=None,
         target_file=None, do_sort=False, sample_name="", config=None):
    with utils.curdir_tmpdir() as tmp_dir:
        work_dir = os.getcwd()
        params = {}
        java_memory = ""
        if config:
            with open(config) as in_handle:
                info = yaml.load(in_handle)
                params = info["program"]
                java_memory = info["algorithm"].get("java_memory", "")
        picard = BroadRunner(picard_dir, max_memory=java_memory)
        if do_sort:
            align_bam = picard_sort(picard, align_bam, tmp_dir)

        metrics = PicardMetrics(picard, tmp_dir)
        summary_table, metrics_graphs = metrics.report(
                align_bam, ref_file, is_paired, bait_file, target_file)
        metrics_graphs = [(p, c, 0.75) for p, c in metrics_graphs]
        base, ext = os.path.splitext(align_bam)
        base = base.replace(".", "-")
        fastqc_graphs, fastqc_stats, fastqc_overrep = \
                       fastqc_report(align_bam, params)

        all_graphs = fastqc_graphs + metrics_graphs
        summary_table = _update_summary_table(summary_table, ref_file, fastqc_stats)
        tmpl = Template(section_template)
        if sample_name is None:
            sample_name = fastqc_stats["Filename"]
        sample_name = "%s (%s)" % (sample_name.replace("_", "\_"),
                base.replace("_", " "))
        section = tmpl.render(name=sample_name, summary=None,
                              summary_table=summary_table,
                              figures=[(f, c, i) for (f, c, i) in all_graphs if f],
                              overrep=fastqc_overrep,
                              recal_figures=_get_recal_plots(work_dir, align_bam))
        out_file = "%s-summary.tex" % base
        out_tmpl = Template(base_template)
        with open(out_file, "w") as out_handle:
            out_handle.write(out_tmpl.render(parts=[section]))
        run_pdflatex(out_file, params)

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
            k, v = [self._safe_latex(x) for x in stat_line.split("\t")[:2]]
            stats[k] = v
        over_rep = []
        for line in self._fastqc_data_section("Overrepresented sequences")[1:]:
            parts = [self._safe_latex(x) for x in line.split("\t")]
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

    def _safe_latex(self, to_fix):
        """Escape characters that make LaTeX unhappy.
        """
        chars = ["%", "_", "&"]
        for char in chars:
            to_fix = to_fix.replace(char, "\\%s" % char)
        return to_fix

def _run_fastqc(bam_file, config):
    out_base = "fastqc"
    utils.safe_makedir(out_base)
    fastqc_out = os.path.join(out_base, "%s_fastqc" %
                              os.path.splitext(os.path.basename(bam_file))[0])
    if not os.path.exists(fastqc_out):
        cl = [config.get("fastqc", "fastqc"), "-o", out_base, "-f", "bam", bam_file]
        subprocess.check_call(cl)
    if os.path.exists("%s.zip" % fastqc_out):
        os.remove("%s.zip" % fastqc_out)
    return fastqc_out

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
    cl = [params.get("pdflatex", "pdflatex"), tex_file]
    subprocess.check_call(cl)

def _get_recal_plots(work_dir, align_bam):
    """Retrieve any recalibration report plots for display.
    """
    (base, _) = os.path.splitext(align_bam)
    reports = glob.glob(os.path.join(work_dir, "reports", "images",
        "%s*-plot.pdf" % base))
    reports.sort()
    return reports

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
    Sequence & Count & Percent & Match \\ 
    \hline
    % for seq, count, percent, match in overrep:
        \texttt{${seq}} & ${count} & ${"%.2f" % float(percent)} & ${match} \\ 
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
    parser.add_option("-p", "--paired", dest="is_paired", action="store_true",
            default=False)
    parser.add_option("-c", "--config", dest="config")
    (options, args) = parser.parse_args()
    if len(args) < 3:
        print __doc__
        sys.exit()
    kwargs = dict(
        is_paired=options.is_paired,
        bait_file=options.bait_file,
        target_file=options.target_file,
        sample_name=options.sample_name,
        config=options.config,
        do_sort=options.do_sort)
    main(*args, **kwargs)
