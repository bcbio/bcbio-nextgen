"""
perform exon-level counting using DEXSeq
"""
import sys
import os
from bcbio.utils import R_package_path, file_exists, safe_makedir
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio import bam
from bcbio.log import logger
import bcbio.pipeline.datadict as dd

def bcbio_run(data):
    out_dir = os.path.join(dd.get_work_dir(data), "dexseq")
    safe_makedir(out_dir)
    sample_name = dd.get_sample_name(data)
    out_file = os.path.join(out_dir, sample_name + ".dexseq")
    bam_file = dd.get_work_bam(data)
    dexseq_gff = dd.get_dexseq_gff(data)
    stranded = dd.get_strandedness(data)
    return run_count(bam_file, dexseq_gff, stranded, out_file, data)

def run_count(bam_file, dexseq_gff, stranded, out_file, data):
    """
    run dexseq_count on a BAM file
    """
    assert file_exists(bam_file), "%s does not exist." % bam_file
    sort_order = bam._get_sort_order(bam_file, {})
    assert sort_order, "Cannot determine sort order of %s." % bam_file
    strand_flag = _strand_flag(stranded)
    assert strand_flag, "%s is not a valid strandedness value." % stranded
    if not file_exists(dexseq_gff):
        logger.info("%s was not found, so exon-level counting is being "
                    "skipped." % dexseq_gff)
        return None

    dexseq_count = _dexseq_count_path()
    if not dexseq_count:
        logger.info("DEXseq is not installed, skipping exon-level counting.")
        return None

    sort_flag = "name" if sort_order == "queryname" else "pos"
    is_paired = bam.is_paired(bam_file)
    paired_flag = "yes" if is_paired else "no"
    bcbio_python = sys.executable

    if file_exists(out_file):
        return out_file
    cmd = ("{bcbio_python} {dexseq_count} -f bam -r {sort_flag} -p {paired_flag} "
           "-s {strand_flag} {dexseq_gff} {bam_file} {tx_out_file}")
    message = "Counting exon-level counts with %s and %s." % (bam_file, dexseq_gff)
    with file_transaction(data, out_file) as tx_out_file:
        do.run(cmd.format(**locals()), message)
    return out_file

def _strand_flag(stranded):
    strand_flag = {"unstranded": "no",
                   "firststrand": "reverse",
                   "secondstrand": "yes"}
    return strand_flag.get(stranded, None)

def _dexseq_count_path():
    package_path = R_package_path("DEXSeq")
    if not package_path:
        return None
    return os.path.join(package_path, "python_scripts", "dexseq_count.py")

def _dexseq_gtf_path(genome_dir):
    return os.path.join(genome_dir, "rnaseq", "ref-transcripts.dexseq.gff")
