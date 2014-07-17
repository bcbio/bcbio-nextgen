"""
perform exon-level counting using DEXSeq
"""
import os
import sys
from bcbio.utils import R_package_path, file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio import bam

def run_count(bam_file, dexseq_gff, stranded, out_file):
    """
    run dexseq_count on a BAM file
    """
    assert file_exists(bam_file), "%s does not exist." % bam_file
    assert file_exists(dexseq_gff), "%s does not exist." % dexseq_gff
    sort_order = bam._get_sort_order(bam_file, {})
    assert sort_order, "Cannot determine sort order of %s." % bam_file
    strand_flag = _strand_flag(stranded)
    assert strand_flag, "%s is not a valid strandedness value." % stranded

    dexseq_count = _dexseq_count_path()
    if not dexseq_count:
        logger.error("DEXseq is not installed, aborting.")
        sys.exit(1)

    sort_flag = "name" if sort_order == "queryname" else "pos"
    is_paired = bam.is_paired(bam_file)
    paired_flag = "yes" if is_paired else "no"

    if file_exists(out_file):
        return out_file
    cmd = ("python {dexseq_count} -f bam -r {sort_flag} -p {paired_flag} "
           "-s {strand_flag} {dexseq_gff} {bam_file} {tx_out_file}")
    message = "Counting exon-level counts with %s and %s." % (bam_file, dexseq_gff)
    with file_transaction(out_file) as tx_out_file:
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
