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
    counts = run_count(bam_file, dexseq_gff, stranded, out_file, data)
    data = dd.set_dexseq_counts(data, counts)
    return data

def run_count(bam_file, dexseq_gff, stranded, out_file, data):
    """
    run dexseq_count on a BAM file
    """
    assert file_exists(bam_file), "%s does not exist." % bam_file
    sort_order = bam._get_sort_order(bam_file, {})
    assert sort_order, "Cannot determine sort order of %s." % bam_file
    strand_flag = _strand_flag(stranded)
    assert strand_flag, "%s is not a valid strandedness value." % stranded
    if not dexseq_gff:
        logger.info("No DEXSeq GFF file was found, skipping exon-level counting.")
        return None
    elif not file_exists(dexseq_gff):
        logger.info("%s was not found, so exon-level counting is being "
                    "skipped." % dexseq_gff)
        return None

    dexseq_count = _dexseq_count_path()
    if not dexseq_count:
        logger.info("DEXseq is not installed, skipping exon-level counting.")
        return None

    if dd.get_aligner(data) == "bwa":
        logger.info("Can't use DEXSeq with bwa alignments, skipping exon-level counting.")
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

def create_dexseq_annotation(gff, count_file):
    """
    Create an easy data frame to allow easy annotation
    during differential expression analysis i.e
    gene:exon_id chr start end strand
    """
    out_file = count_file + ".ann"
    if file_exists(out_file):
        return out_file
    with file_transaction(out_file) as tx_out:
        with open(tx_out, 'w') as out_handle:
            with open(gff) as in_handle:
                for line in in_handle:
                    cols = line.strip().split("\t")
                    if cols[2] == "exonic_part":
                        exon = [f for f in cols[8].split(";") if f.strip().startswith("exonic_part_number")]
                        gene = [f for f in cols[8].split(";") if f.strip().startswith("gene_id")]
                        exon = exon[0].replace("\"", "").split()[1]
                        gene = gene[0].replace("\"", "").split()[1]
                        length = int(cols[4]) - int(cols[3]) + 1
                        line = "%s:%s\t%s\t%s\t%s\t%s\t%s\t%s" % (gene, exon, gene,
                                                                  cols[0], cols[3],
                                                                  cols[4],
                                                                  length, cols[6])
                        print >>out_handle, line
