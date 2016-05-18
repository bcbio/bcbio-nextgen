import os
import string
import os.path as op
import sys
import shutil
from collections import namedtuple

try:
    from seqcluster import prepare_data as prepare
    from seqcluster import templates as template_seqcluster
    from seqcluster.seqbuster import _create_counts, _read_miraligner, _tab_output
except ImportError:
    pass

from bcbio.utils import file_exists, safe_makedir, move_safe, append_stem
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.pipeline.sample import process_alignment
from bcbio.srna import mirdeep

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
    out_dir = op.join(work_dir, "seqcluster", "prepare")
    seq_out = op.join(out_dir, "seqs.fastq")
    bam_dir = op.join(work_dir, "align")
    new_bam_file = op.join(bam_dir, "seqs.bam")
    if not file_exists(new_bam_file):
        sample = process_alignment(data[0][0], [seq_out, None])
        bam_file = dd.get_work_bam(sample[0][0])
        shutil.move(bam_file, new_bam_file)
        shutil.move(bam_file + ".bai", new_bam_file + ".bai")
        shutil.rmtree(op.join(bam_dir, sample[0][0]["rgnames"]['sample']))
    return data

def run_cluster(*data):
    """
    Run seqcluster cluster to detect smallRNA clusters
    """
    sample = data[0][0]
    tools_off = dd.get_tools_off(data[0][0])
    work_dir = dd.get_work_dir(sample)
    out_dir = op.join(work_dir, "seqcluster", "cluster")
    out_dir = op.abspath(safe_makedir(out_dir))
    prepare_dir = op.join(work_dir, "seqcluster", "prepare")
    bam_file = op.join(work_dir, "align", "seqs.bam")
    if "seqcluster" not in tools_off:
        cluster_dir = _cluster(bam_file, prepare_dir, out_dir, dd.get_ref_file(sample), dd.get_srna_gtf_file(sample))
        sample["report"] = _report(sample, dd.get_ref_file(sample))
        sample["seqcluster"] = out_dir

    out_mirna = _make_isomir_counts(data, out_dir=op.join(work_dir, "mirbase"))
    if out_mirna:
        sample = dd.set_mirna_counts(sample, out_mirna[0])
        sample = dd.set_isomir_counts(sample, out_mirna[1])

    out_novel = _make_isomir_counts(data, "seqbuster_novel", op.join(work_dir, "mirdeep2"), "_novel")
    novel_db = mirdeep.run(data)
    if out_novel:
        sample = dd.set_novel_mirna_counts(sample, out_novel[0])
        sample = dd.set_novel_isomir_counts(sample, out_novel[1])
    data[0][0] = sample
    return data

def _cluster(bam_file, prepare_dir, out_dir, reference, annotation_file=None):
    """
    Connect to seqcluster to run cluster with python directly
    """
    seqcluster = op.join(os.path.dirname(sys.executable), "seqcluster")
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
    seqcluster = op.join(os.path.dirname(sys.executable), "seqcluster")
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
    out_dir = op.join(work_dir, "report")
    safe_makedir(out_dir)
    summary_file = op.join(out_dir, "summary.csv")
    with file_transaction(summary_file) as out_tx:
        with open(out_tx, 'w') as out_handle:
            print >>out_handle, "sample_id,%s" % _guess_header(data[0][0])
            for sample in data:
                info = sample[0]
                group = _guess_group(info)
                files = info["seqbuster"] if "seqbuster" in info else "None"
                print >>out_handle, ",".join([dd.get_sample_name(info),
                                              group])
    _modify_report(work_dir, out_dir)
    return summary_file

def _guess_header(info):
    """Add the first group to get report with some factor"""
    value = "group"
    if "metadata" in info:
        if info["metadata"]:
            return ",".join(info["metadata"].keys())
    return value

def _guess_group(info):
    """Add the first group to get report with some factor"""
    value = "fake"
    if "metadata" in info:
        if info["metadata"]:
            return ",".join(info["metadata"].values())
    return value

def _modify_report(summary_path, out_dir):
    """Read Rmd template and dump with project path."""
    summary_path = op.abspath(summary_path)
    template = op.normpath(op.join(op.dirname(op.realpath(template_seqcluster.__file__)), "report.rmd"))
    content = open(template).read()
    out_content = string.Template(content).safe_substitute({'path_abs': summary_path})
    out_file = op.join(out_dir, "srna_report.rmd")
    with open(out_file, 'w') as out_handle:
        print >>out_handle, out_content
    return out_file

def _make_isomir_counts(data, srna_type="seqbuster", out_dir=None, stem=""):
    """
    Parse miraligner files to create count matrix.
    """
    work_dir = dd.get_work_dir(data[0][0])
    if not out_dir:
        out_dir = op.join(work_dir, "mirbase")
    out_novel_isomir = append_stem(op.join(out_dir, "counts.tsv"), stem)
    out_novel_mirna = append_stem(op.join(out_dir, "counts_mirna.tsv"), stem)
    logger.debug("Create %s count data at %s." % (srna_type, out_dir))
    if file_exists(out_novel_mirna):
        return [out_novel_mirna, out_novel_isomir]
    out_dts = []
    for sample in data:
        if sample[0].get(srna_type):
            miraligner_fn = sample[0][srna_type]
            reads = _read_miraligner(miraligner_fn)
            if reads:
                out_file, dt, dt_pre = _tab_output(reads, miraligner_fn + ".back", dd.get_sample_name(sample[0]))
                out_dts.append(dt)
            else:
                logger.debug("WARNING::%s has NOT miRNA annotated for %s. Check if fasta files is small or species value." % (dd.get_sample_name(sample[0]), srna_type))
    if out_dts:
        out_files = _create_counts(out_dts, out_dir)
        out_files = [move_safe(out_files[0], out_novel_isomir), move_safe(out_files[1], out_novel_mirna)]
        return out_files
    else:
        logger.debug("WARNING::any samples have miRNA annotated for %s. Check if fasta files is small or species value." % srna_type)
