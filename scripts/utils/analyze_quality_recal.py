#!/usr/bin/env python
"""Provide plots summarizing recalibration of quality scores.

Usage:
    analyze_quality_recal.py <recal_bam> <input_fastq1> <input_fastq2>
    --chunk_size=25 --input_format=fastq-illumina
    --dbdir=/tmp/chapmanb

    <recal_bam> is a BAM alignment file containing recalibrarted quality scores
    <input_fastq> are the initial fastq files with quality scores
    --chunk_size -- How many positions to read in at once. Higher scores are
    faster but require more memory.
    --input_format -- Quality score encoding of input fastq files.
    --dbdir -- Where to store database files. This is needed for cluster
    jobs on NFS since sqlite can behave strangely on NFS with lock errors.
    --workdir -- The working directory to write output files to. Defaults to the current
      directory

Requirements:
    sqlite
    biopython
    pysam
    R with ggplot2, plyr and sqldf
    rpy2
    mako
    latex (texlive)
"""
import sys
import os
import csv
import glob
import collections
import subprocess
from optparse import OptionParser
try:
    from pysqlite2 import dbapi2 as sqlite3
except ImportError:
    import sqlite3

from Bio import SeqIO
from Bio import Seq
import pysam
from mako.template import Template
try:
    import rpy2.robjects as robjects
except (ImportError, LookupError):
    robjects = None

def main(recal_bam, fastq1, fastq2=None, chunk_size=None, input_format=None,
        db_dir=None, work_dir=None):
    if not _are_libraries_installed():
        print "R libraries or rpy2 not installed. Not running recalibration plot."
        return
    if work_dir is None:
        work_dir = os.getcwd()
    report_dir = os.path.join(work_dir, "reports")
    image_dir = os.path.join(report_dir, "images")
    if db_dir is None:
        db_dir = work_dir
    if not os.path.exists(image_dir):
        # avoid error with creating directories simultaneously on two threads
        try:
            os.makedirs(image_dir)
        except OSError:
            assert os.path.isdir(image_dir)
    base = os.path.splitext(os.path.basename(recal_bam))[0]
    orig_files = {1: fastq1, 2: fastq2}

    section_info = []
    pairs = ([1] if fastq2 is None else [1, 2])
    for pair in pairs:
        plots = []
        db_file = os.path.join(db_dir, "%s_%s-qualities.sqlite" % (base, pair))
        if not os.path.exists(db_file):
            print "Converting BAM alignment to fastq files"
            recal_fastq1, recal_fastq2 = bam_to_fastq(recal_bam, len(pairs) > 1)
            recal_files = {1: recal_fastq1, 2: recal_fastq2}
            print "Normalizing and sorting fastq files"
            orig = sort_csv(fastq_to_csv(orig_files[pair], input_format,
                work_dir))
            recal = sort_csv(fastq_to_csv(recal_files[pair], "fastq", work_dir))
            print "Summarizing remapped qualities for pair", pair
            summarize_qualities(db_file, orig, recal, chunk_size)
        print "Plotting for pair", pair
        for position_select, pname in _positions_to_examine(db_file):
            title = "Pair %s; Position: %s" % (pair, position_select)
            plot_file = os.path.join(image_dir,
                    "%s_%s_%s-plot.pdf" % (base, pair, pname))
            draw_quality_plot(db_file, plot_file, position_select, title)
            plots.append(plot_file)
        section_info.append(("Pair %s" % pair, plots))

    run_latex_report(base, report_dir, section_info)
    _clean_intermediates(recal_bam, fastq1, fastq2, report_dir)

def _are_libraries_installed():
    if robjects is None:
        print "rpy2 not installed: http://rpy.sourceforge.net/rpy2.html"
        return False
    import rpy2.rinterface
    try:
        robjects.r('''
          library(sqldf)
          library(plyr)
          library(ggplot2)
        ''')
    except rpy2.rinterface.RRuntimeError:
        print "Some R libraries not installed"
        return False
    return True

def draw_quality_plot(db_file, plot_file, position_select, title):
    """Draw a plot of remapped qualities using ggplot2.

    Remapping information is pulled from the sqlite3 database using sqldf
    according to the position select attribute, which is a selection phrase like
    '> 50' or '=28'.

    plyr is used to summarize data by the original and remapped score for all
    selected positions.

    ggplot2 plots a heatmap of remapped counts at each (original, remap)
    coordinate, with a x=y line added for reference.
    """
    robjects.r.assign('db.file', db_file)
    robjects.r.assign('plot.file', plot_file)
    robjects.r.assign('position.select', position_select)
    robjects.r.assign('title', title)
    robjects.r('''
      library(sqldf)
      library(plyr)
      library(ggplot2)
      sql <- paste("select * from data WHERE position", position.select, sep=" ")
      exp.data <- sqldf(sql, dbname=db.file)
      remap.data <- ddply(exp.data, c("orig", "remap"), transform, count=sum(count))
      p <- ggplot(remap.data, aes(orig, remap)) +
           geom_tile(aes(fill = count)) +
           scale_fill_gradient(low = "white", high = "steelblue", trans="log") +
           opts(panel.background = theme_rect(fill = "white"),
                title=title) +
           geom_abline(intercept=0, slope=1)
      ggsave(plot.file, p, width=6, height=6)
    ''')

def _positions_to_examine(db_file):
    """Determine how to sub-divide recalibration analysis based on read length.
    """
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("""SELECT MAX(position) FROM data""")
    position = cursor.fetchone()[0]
    if position is not None:
        position = int(position)
    cursor.close()
    split_at = 50
    if position is None:
        return []
    elif position < split_at:
        return [("<= %s" % position, "lt%s" % position)]
    else:
        return [("< %s" % split_at, "lt%s" % split_at),
                (">= %s" % split_at, "gt%s" % split_at)]

def summarize_qualities(db_file, orig_file, cmp_file, chunk_size):
    out_conn = sqlite3.connect(db_file)
    out_cursor = out_conn.cursor()
    out_cursor.execute("""create table data
                (position integer, orig integer,
                 remap integer, count integer)""")
    try:
        cur_pos = 1
        for pos, orig_val, final_val, count in _organize_by_position(orig_file,
                cmp_file, chunk_size):
            out_cursor.execute("INSERT INTO data VALUES (?,?,?,?)",
                    [pos, orig_val, final_val, count])
            if pos != cur_pos:
                cur_pos = pos
                out_conn.commit()
    finally:
        out_conn.commit()
        out_cursor.close()

def _organize_by_position(orig_file, cmp_file, chunk_size):
    """Read two CSV files of qualities, organizing values by position.
    """
    with open(orig_file) as in_handle:
        reader1 = csv.reader(in_handle)
        positions = len(next(reader1)) - 1
    for positions in _chunks(range(positions), chunk_size):
        with open(orig_file) as orig_handle:
            with open(cmp_file) as cmp_handle:
                orig_reader = csv.reader(orig_handle)
                cmp_reader = csv.reader(cmp_handle)
                for item in _counts_at_position(positions,
                        orig_reader, cmp_reader):
                    yield item

def _chunks(l, n):
    """ Yield successive n-sized chunks from l.

    http://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks-in-python
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def _counts_at_position(positions, orig_reader, cmp_reader):
    """Combine orignal and new qualities at each position, generating counts.
    """
    pos_counts = collections.defaultdict(lambda:
                 collections.defaultdict(lambda:
                 collections.defaultdict(int)))
    for orig_parts in orig_reader:
        cmp_parts = next(cmp_reader)
        for pos in positions:
            try:
                pos_counts[pos][int(orig_parts[pos+1])][int(cmp_parts[pos+1])] += 1
            except IndexError:
                pass
    for pos, count_dict in pos_counts.iteritems():
        for orig_val, cmp_dict in count_dict.iteritems():
            for cmp_val, count in cmp_dict.iteritems():
                yield pos+1, orig_val, cmp_val, count

def sort_csv(in_file):
    """Sort a CSV file by read name, allowing direct comparison.
    """
    out_file = "%s.sort" % in_file
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        cl = ["sort", "-k", "1,1", in_file]
        with open(out_file, "w") as out_handle:
            child = subprocess.Popen(cl, stdout=out_handle)
            child.wait()
    return out_file

def fastq_to_csv(in_file, fastq_format, work_dir):
    """Convert a fastq file into a CSV of phred quality scores.
    """
    out_file = "%s.csv" % (os.path.splitext(os.path.basename(in_file))[0])
    out_file = os.path.join(work_dir, out_file)
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        with open(in_file) as in_handle:
            with open(out_file, "w") as out_handle:
                writer = csv.writer(out_handle)
                for rec in SeqIO.parse(in_handle, fastq_format):
                    writer.writerow([rec.id] + rec.letter_annotations["phred_quality"])
    return out_file

def bam_to_fastq(bam_file, is_paired):
    """Convert a BAM file to fastq files.
    """
    out_files, out_handles = _get_fastq_handles(bam_file,
            is_paired)
    if len(out_handles) > 0:
        in_bam = pysam.Samfile(bam_file, mode='rb')
        for read in in_bam:
            num = 1 if (not read.is_paired or read.is_read1) else 2
            # reverse the sequence and quality if mapped to opposite strand
            if read.is_reverse:
                seq = str(Seq.reverse_complement(Seq.Seq(read.seq)))
                qual = "".join(reversed(read.qual))
            else:
                seq = read.seq
                qual = read.qual
            out_handles[num].write("@%s\n%s\n+\n%s\n" % (read.qname,
                seq, qual))
    [h.close() for h in out_handles.values()]
    return out_files

def _get_fastq_handles(bam_file, is_paired):
    (base, _) = os.path.splitext(bam_file)
    out_files = []
    out_handles = dict()
    if is_paired:
        for index in [1, 2]:
            cur_file = "%s_%s_fastq.txt" % (base, index)
            out_files.append(cur_file)
            if not (os.path.exists(cur_file) and os.path.getsize(cur_file) > 0):
                out_handles[index] = open(cur_file, "w")
    else:
        cur_file = "%s_fastq.txt" % base
        out_files.append(cur_file)
        out_files.append(None)
        if not(os.path.exists(cur_file) and os.path.getsize(cur_file) > 0):
            out_handles[1] = open(cur_file, "w")
    return out_files, out_handles

def _clean_intermediates(bam_file, fastq1, fastq2, report_dir):
    base = os.path.splitext(bam_file)[0]
    for bam_rem in glob.glob("%s_*fastq*" % base):
        os.remove(bam_rem)
    for fastq in (fastq1, fastq2):
        if fastq:
            for fastq_rem in glob.glob("%s.csv*" %
                    os.path.splitext(os.path.basename(fastq))[0]):
                os.remove(fastq_rem)
    for latex_ext in ["aux", "log"]:
        for latex_rem in glob.glob(os.path.join(report_dir, "%s*.%s" %
                            (os.path.basename(base), latex_ext))):
            os.remove(latex_rem)

def run_latex_report(base, report_dir, section_info):
    """Generate a pdf report with plots using latex.
    """
    out_name = "%s_recal_plots.tex" % base
    out = os.path.join(report_dir, out_name)
    with open(out, "w") as out_handle:
        out_tmpl = Template(out_template)
        out_handle.write(out_tmpl.render(sections=section_info))
    start_dir = os.getcwd()
    try:
        os.chdir(report_dir)
        cl = ["pdflatex", out_name]
        child = subprocess.Popen(cl)
        child.wait()
    finally:
        os.chdir(start_dir)

out_template = r"""
\documentclass{article}
\usepackage{fullpage}
\usepackage[top=0.5in,right=1in,left=1in,bottom=0.5in]{geometry}
\usepackage{graphicx}
\usepackage{placeins}

\begin{document}
% for section, figures in sections:
    \subsubsection*{${section}}
    % for figure in figures:
        \begin{figure}[htbp]
          \centering
          \includegraphics[width=0.57\linewidth]{${figure}}
        \end{figure}
    % endfor
   \FloatBarrier
   \newpage
% endfor
\end{document}
"""

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-c", "--chunk_size", dest="chunk_size", default=25)
    parser.add_option("-i", "--input_format", dest="input_format",
            default="fastq-illumina")
    parser.add_option("-d", "--dbdir", dest="dbdir",
            default=None)
    parser.add_option("-w", "--workdir", dest="workdir",
            default=None)
    (options, args) = parser.parse_args()
    kwargs = dict(chunk_size = int(options.chunk_size),
                  input_format = options.input_format,
                  db_dir = options.dbdir, work_dir = options.workdir)
    main(*args, **kwargs)
