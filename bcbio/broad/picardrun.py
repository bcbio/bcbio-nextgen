"""Convenience functions for running common Picard utilities.
"""
import os
import collections
from contextlib import closing

import pysam

from bcbio.utils import curdir_tmpdir, file_exists
from bcbio.distributed.transaction import file_transaction


def picard_rnaseq_metrics(picard, align_bam, ref, ribo="null", out_file=None):
    """ Collect RNASeq metrics for a bam file """
    base, ext = os.path.splitext(align_bam)
    if out_file is None:
        out_file = "%s.metrics" % (base)
    if not file_exists(out_file):
        with curdir_tmpdir() as tmp_dir:
            with file_transaction(out_file) as tx_out_file:
                opts = [("INPUT", align_bam),
                        ("OUTPUT", tx_out_file),
                        ("TMP_DIR", tmp_dir),
                        ("REF_FLAT", ref),
                        ("STRAND_SPECIFICITY", "NONE"),
                        ("ASSUME_SORTED", "True"),
                        ("RIBOSOMAL_INTERVALS", ribo)]

                picard.run("CollectRnaSeqMetrics", opts)
    return out_file


def picard_sort(picard, align_bam, sort_order="coordinate",
                out_file=None, compression_level=None, pipe=False):
    """Sort a BAM file by coordinates.
    """
    base, ext = os.path.splitext(align_bam)
    if out_file is None:
        out_file = "%s-sort%s" % (base, ext)
    if not file_exists(out_file):
        with curdir_tmpdir() as tmp_dir:
            with file_transaction(out_file) as tx_out_file:
                opts = [("INPUT", align_bam),
                        ("OUTPUT", out_file if pipe else tx_out_file),
                        ("TMP_DIR", tmp_dir),
                        ("SORT_ORDER", sort_order)]
                if compression_level:
                    opts.append(("COMPRESSION_LEVEL", compression_level))
                picard.run("SortSam", opts, pipe=pipe)
    return out_file

def picard_merge(picard, in_files, out_file=None,
                 merge_seq_dicts=False):
    """Merge multiple BAM files together with Picard.
    """
    if out_file is None:
        out_file = "%smerge.bam" % os.path.commonprefix(in_files)
    if not file_exists(out_file):
        with curdir_tmpdir() as tmp_dir:
            with file_transaction(out_file) as tx_out_file:
                opts = [("OUTPUT", tx_out_file),
                        ("SORT_ORDER", "coordinate"),
                        ("MERGE_SEQUENCE_DICTIONARIES",
                         "true" if merge_seq_dicts else "false"),
                        ("USE_THREADING", "true"),
                        ("TMP_DIR", tmp_dir)]
                for in_file in in_files:
                    opts.append(("INPUT", in_file))
                picard.run("MergeSamFiles", opts)
    return out_file

def picard_index(picard, in_bam):
    index_file = "%s.bai" % in_bam
    if not file_exists(index_file):
        with file_transaction(index_file) as tx_index_file:
            opts = [("INPUT", in_bam),
                    ("OUTPUT", tx_index_file)]
            picard.run("BuildBamIndex", opts)
    return index_file

def picard_index_ref(picard, ref_file):
    """Provide a Picard style dict index file for a reference genome.
    """
    dict_file = "%s.dict" % os.path.splitext(ref_file)[0]
    if not file_exists(dict_file):
        with file_transaction(dict_file) as tx_dict_file:
            opts = [("REFERENCE", ref_file),
                    ("OUTPUT", tx_dict_file)]
            picard.run("CreateSequenceDictionary", opts)
    return dict_file

def picard_fastq_to_bam(picard, fastq_one, fastq_two, out_dir,
                        platform, sample_name="", rg_name="", pu_name="",
                        qual_format=None):
    """Convert fastq file(s) to BAM, adding sample, run group and platform information.
    """
    qual_formats = {"illumina": "Illumina"}
    if qual_format is None:
        try:
            qual_format = qual_formats[platform.lower()]
        except KeyError:
            raise ValueError("Need to specify quality format for %s" % platform)
    out_bam = os.path.join(out_dir, "%s-fastq.bam" %
                           os.path.splitext(os.path.basename(fastq_one))[0])
    if not file_exists(out_bam):
        with curdir_tmpdir() as tmp_dir:
            with file_transaction(out_bam) as tx_out_bam:
                opts = [("FASTQ", fastq_one),
                        ("QUALITY_FORMAT", qual_format),
                        ("READ_GROUP_NAME", rg_name),
                        ("SAMPLE_NAME", sample_name),
                        ("PLATFORM_UNIT", pu_name),
                        ("PLATFORM", platform),
                        ("TMP_DIR", tmp_dir),
                        ("OUTPUT", tx_out_bam)]
                if fastq_two:
                    opts.append(("FASTQ2", fastq_two))
                picard.run("FastqToSam", opts)
    return out_bam

def picard_bam_to_fastq(picard, in_bam, fastq_one, fastq_two=None):
    """Convert BAM file to fastq.
    """
    if not file_exists(fastq_one):
        with curdir_tmpdir() as tmp_dir:
            with file_transaction(fastq_one) as tx_out1:
                opts = [("INPUT", in_bam),
                        ("FASTQ", tx_out1),
                        ("TMP_DIR", tmp_dir)]
                if fastq_two is not None:
                    opts += [("SECOND_END_FASTQ", fastq_two)]
                picard.run("SamToFastq", opts)
    return (fastq_one, fastq_two)

def picard_sam_to_bam(picard, align_sam, fastq_bam, ref_file,
                      is_paired=False):
    """Convert SAM to BAM, including unmapped reads from fastq BAM file.
    """
    if align_sam.endswith(".sam"):
        out_bam = "%s.bam" % os.path.splitext(align_sam)[0]
    elif align_sam.endswith("-align.bam"):
        out_bam = "%s.bam" % align_sam.replace("-align.bam", "")
    else:
        raise NotImplementedError("Input format not recognized")
    if not file_exists(out_bam):
        with curdir_tmpdir() as tmp_dir:
            with file_transaction(out_bam) as tx_out_bam:
                opts = [("UNMAPPED", fastq_bam),
                        ("ALIGNED", align_sam),
                        ("OUTPUT", tx_out_bam),
                        ("REFERENCE_SEQUENCE", ref_file),
                        ("TMP_DIR", tmp_dir),
                        ("PAIRED_RUN", ("true" if is_paired else "false")),
                        ]
                picard.run("MergeBamAlignment", opts)
    return out_bam

def picard_formatconverter(picard, align_sam):
    """Convert aligned SAM file to BAM format.
    """
    out_bam = "%s.bam" % os.path.splitext(align_sam)[0]
    if not file_exists(out_bam):
        with curdir_tmpdir() as tmp_dir:
            with file_transaction(out_bam) as tx_out_bam:
                opts = [("INPUT", align_sam),
                        ("OUTPUT", tx_out_bam)]
                picard.run("SamFormatConverter", opts)
    return out_bam

def picard_mark_duplicates(picard, align_bam, remove_dups=False):
    base, ext = os.path.splitext(align_bam)
    base = base.replace(".", "-")
    dup_bam = "%s-dup%s" % (base, ext)
    dup_metrics = "%s-dup.dup_metrics" % base
    if not file_exists(dup_bam):
        with curdir_tmpdir() as tmp_dir:
            with file_transaction(dup_bam, dup_metrics) as (tx_dup_bam, tx_dup_metrics):
                opts = [("INPUT", align_bam),
                        ("OUTPUT", tx_dup_bam),
                        ("TMP_DIR", tmp_dir),
                        ("REMOVE_DUPLICATES", "true" if remove_dups else "false"),
                        ("METRICS_FILE", tx_dup_metrics)]
                if picard.get_picard_version("MarkDuplicates") >= 1.82:
                    opts += [("PROGRAM_RECORD_ID", "null")]
                picard.run("MarkDuplicates", opts)
    return dup_bam, dup_metrics

def picard_fixmate(picard, align_bam):
    """Run Picard's FixMateInformation generating an aligned output file.
    """
    base, ext = os.path.splitext(align_bam)
    out_file = "%s-sort%s" % (base, ext)
    if not file_exists(out_file):
        with curdir_tmpdir() as tmp_dir:
            with file_transaction(out_file) as tx_out_file:
                opts = [("INPUT", align_bam),
                        ("OUTPUT", tx_out_file),
                        ("TMP_DIR", tmp_dir),
                        ("SORT_ORDER", "coordinate")]
                picard.run("FixMateInformation", opts)
    return out_file

def picard_idxstats(picard, align_bam):
    """Retrieve alignment stats from picard using BamIndexStats.
    """
    opts = [("INPUT", align_bam)]
    stdout = picard.run("BamIndexStats", opts, get_stdout=True)
    out = []
    AlignInfo = collections.namedtuple("AlignInfo", ["contig", "length", "aligned", "unaligned"])
    for line in stdout.split("\n"):
        if line:
            parts = line.split()
            if len(parts) == 2:
                _, unaligned = parts
                out.append(AlignInfo("nocontig", 0, 0, int(unaligned)))
            else:
                contig, _, length, _, aligned, _, unaligned = parts
                out.append(AlignInfo(contig, int(length), int(aligned), int(unaligned)))
    return out

def bed2interval(align_file, bed, out_file=None):
    """Converts a bed file to an interval file for use with some of the
    Picard tools by grabbing the header from the alignment file, reording
    the bed file columns and gluing them together.

    align_file can be in BAM or SAM format.
    bed needs to be in bed12 format:
    http://genome.ucsc.edu/FAQ/FAQformat.html#format1.5

    """

    base, ext = os.path.splitext(align_file)
    if out_file is None:
        out_file = base + ".interval"

    with closing(pysam.Samfile(align_file, "r" if ext.endswith(".sam") else "rb")) as in_bam:
        header = in_bam.text

    def reorder_line(line):
        splitline = line.strip().split("\t")
        reordered = "\t".join([splitline[0], splitline[1], splitline[2],
                               splitline[5], splitline[3]])
        return reordered + "\n"

    with file_transaction(out_file) as tx_out_file:
        with open(bed) as bed_handle:
            with open(tx_out_file, "w") as out_handle:
                out_handle.write(header)
                for line in bed_handle:
                    out_handle.write(reorder_line(line))
    return out_file
