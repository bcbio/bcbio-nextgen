import os
import string
import argparse
import os.path as op
import sys
import shutil
from collections import namedtuple

import pysam

try:
    from seqcluster import prepare_data as prepare
    from seqcluster import make_clusters as main_cluster
    from seqcluster.libs import parse
    from seqcluster import templates as template_seqcluster
except ImportError:
    pass

from bcbio.utils import file_exists, safe_makedir
from bcbio.provenance import do
from bcbio.distributed.transaction import tx_tmpdir, file_transaction
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.pipeline.sample import process_alignment


def run_prepare(*data):
    """
    Run seqcluster prepare to merge all samples in one file
    """
    out_dir = os.path.join(dd.get_work_dir(data[0][0]), "seqcluster", "prepare")
    out_dir = os.path.abspath(safe_makedir(out_dir))
    prepare_dir = os.path.join(out_dir, "prepare")
    fn = []
    for sample in data:
        name = sample[0]["rgnames"]['sample']
        fn.append("%s\t%s" % (sample[0]['collapse'], name))
    args = namedtuple('args', 'debug print_debug minc minl maxl out')
    args = args(False, False, 2, 17, 40, out_dir)
    ma_out = op.join(out_dir, "seqs.ma")
    seq_out = op.join(out_dir, "seqs.fastq")
    min_shared = max(int(len(fn) / 10.0), 1)
    if not file_exists(ma_out):
        seq_l, sample_l = prepare._read_fastq_files(fn, args)
        with file_transaction(ma_out) as ma_tx:
            with open(ma_tx, 'w') as ma_handle:
                with open(seq_out, 'w') as seq_handle:
                    prepare._create_matrix_uniq_seq(sample_l, seq_l, ma_handle, seq_handle, min_shared)

    return data

def run_align(*data):
    """
    Prepare data to run alignment step, only once for each project
    """
    work_dir = dd.get_work_dir(data[0][0])
    out_dir = os.path.join(work_dir, "seqcluster", "prepare")
    seq_out = op.join(out_dir, "seqs.fastq")
    bam_dir = os.path.join(work_dir, "align")
    new_bam_file = op.join(bam_dir, "seqs.bam")
    if not file_exists(new_bam_file):
        sample = process_alignment(data[0][0], [seq_out, None])
        # data = data[0][0]
        bam_file = dd.get_work_bam(sample[0][0])
        shutil.move(bam_file, new_bam_file)
        shutil.move(bam_file + ".bai", new_bam_file + ".bai")
        shutil.rmtree(op.join(bam_dir, sample[0][0]["rgnames"]['sample']))
    return data

def run_cluster(*data):
    """
    Run seqcluster cluster to detect smallRNA clusters
    """
    work_dir = dd.get_work_dir(data[0][0])
    out_dir = os.path.join(work_dir, "seqcluster", "cluster")
    out_dir = os.path.abspath(safe_makedir(out_dir))
    out_file = os.path.join(out_dir, "seqcluster.json")
    prepare_dir = op.join(work_dir, "seqcluster", "prepare")
    bam_file = op.join(work_dir, "align", "seqs.bam")
    cluster_dir = _cluster(bam_file, prepare_dir, out_dir, dd.get_ref_file(data[0][0]), dd.get_srna_gtf_file(data[0][0]))
    report_file = _report(data[0][0], dd.get_ref_file(data[0][0]))
    for sample in data:
        sample[0]["seqcluster"] = out_dir
    return data

def _cluster(bam_file, prepare_dir, out_dir, reference, annotation_file=None):
    """
    Connect to seqcluster to run cluster with python directly
    """
    seqcluster = os.path.join(os.path.dirname(sys.executable), "seqcluster")
    ma_file = op.join(prepare_dir, "seqs.ma")
    # cl = ["cluster", "-o", out_dir, "-m", ma_file, "-a", bam_file, "-r", reference]
    if annotation_file:
        annotation_file = "-g " + annotation_file
    else:
        annotation_file = ""

    if not file_exists(op.join(out_dir, "counts.tsv")):
        cmd = ("{seqcluster} cluster -o {out_dir} -m {ma_file} -a {bam_file} -r {reference} {annotation_file}")
        do.run(cmd.format(**locals()), "Running seqcluster.")
    return out_dir

def _report(data, reference):
    """
    Run report of seqcluster to get browser options for results
    """
    seqcluster = os.path.join(os.path.dirname(sys.executable), "seqcluster")
    work_dir = dd.get_work_dir(data)
    out_dir = safe_makedir(os.path.join(work_dir, "seqcluster", "report"))
    out_file = op.join(out_dir, "seqcluster.db")
    json = op.join(work_dir, "seqcluster", "cluster", "seqcluster.json")
    cmd = ("{seqcluster} report -o {out_dir} -r {reference} -j {json}")
    if not file_exists(out_file):
        do.run(cmd.format(**locals()), "Run report on clusters")
    return out_file

def report(data):
    """Create a Rmd report for small RNAseq analysis"""
    work_dir = dd.get_work_dir(data[0][0])
    out_dir = os.path.join(work_dir, "report")
    safe_makedir(out_dir)
    summary_file = os.path.join(out_dir, "summary.csv")
    with file_transaction(summary_file) as out_tx:
        with open(out_tx, 'w') as out_handle:
            print >>out_handle, "sample_id,size_stats,miraligner,group"
            for sample in data:
                info = sample[0]
                group = _guess_group(info)
                print >>out_handle, ",".join([dd.get_sample_name(info),
                                              info['size_stats'],
                                              info['seqbuster'], group])
    _create_rmd(summary_file)
    return summary_file

def _guess_group(info):
    """Add the first group to get report with some factor"""
    value = "fake"
    if "metadata" in info:
        if info["metadata"]:
            key, value = info['metadata'].popitem()
    return value

def _create_rmd(summary_fn):
    """Create relatie path files for Rmd report"""
    root_path, fn = os.path.split(os.path.abspath(summary_fn))
    out_file = os.path.join(root_path, fn.replace(".csv", "_re.csv"))
    with open(summary_fn) as in_handle:
        with open(out_file, 'w') as out_handle:
            for line in in_handle:
                cols = line.strip().split(",")
                fix_line = ",".join([os.path.relpath(c, root_path) if os.path.exists(c) else c for c in cols])
                print >>out_handle, fix_line
    report_file = _modify_report(root_path, out_file)

    return out_file, report_file

def _modify_report(summary_path, summary_fn):
    """Read Rmd template and dump with project path."""
    summary_path = os.path.abspath(summary_path)
    template = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(template_seqcluster.__file__)), "report.rmd"))
    content = open(template).read()
    out_content = string.Template(content).safe_substitute({'path_abs': summary_path,
                                                            'path_summary': os.path.join(summary_path, summary_fn)})
    out_file = os.path.join(os.path.dirname(summary_fn), "ready_report.rmd")
    with open(out_file, 'w') as out_handle:
        print >>out_handle, out_content

    return out_file
