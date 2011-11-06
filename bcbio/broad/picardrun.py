"""Convenience functions for running common Picard utilities.
"""
import os

from bcbio.utils import curdir_tmpdir, file_exists
from bcbio.distributed.transaction import file_transaction

def picard_sort(picard, align_bam):
    """Sort a BAM file by coordinates.
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
                picard.run("SortSam", opts)
    return out_file

def picard_merge(picard, in_files, out_file=None):
    """Merge multiple BAM files together with Picard.
    """
    if out_file is None:
        out_file = "%smerge.bam" % os.path.commonprefix(in_files)
    if not file_exists(out_file):
        with curdir_tmpdir() as tmp_dir:
            with file_transaction(out_file) as tx_out_file:
                opts = [("OUTPUT", tx_out_file),
                        ("SORT_ORDER", "coordinate"),
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
    out_bam = os.path.join(out_dir, "%s.bam" %
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

def picard_sam_to_bam(picard, align_sam, fastq_bam, ref_file,
                      is_paired=False):
    """Convert SAM to BAM, including unmapped reads from fastq BAM file.
    """
    out_bam = "%s.bam" % os.path.splitext(align_sam)[0]
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

def picard_mark_duplicates(picard, align_bam):
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
                        ("METRICS_FILE", tx_dup_metrics)]
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
