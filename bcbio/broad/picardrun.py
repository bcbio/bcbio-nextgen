"""Convenience functions for running common Picard utilities.
"""
import os

from bcbio.utils import curdir_tmpdir, file_transaction

def picard_sort(picard, align_bam):
    """Sort a BAM file by coordinates.
    """
    base, ext = os.path.splitext(align_bam)
    out_file = "%s-sort%s" % (base, ext)
    if not os.path.exists(out_file):
        with curdir_tmpdir() as tmp_dir:
            opts = [("INPUT", align_bam),
                    ("OUTPUT", out_file),
                    ("TMP_DIR", tmp_dir),
                    ("SORT_ORDER", "coordinate")]
            with file_transaction(out_file):
                picard.run("SortSam", opts)
    return out_file

def picard_merge(picard, in_files, out_file=None):
    """Merge multiple BAM files together with Picard.
    """
    if out_file is None:
        out_file = "%smerge.bam" % os.path.commonprefix(in_files)
    if not os.path.exists(out_file):
        with curdir_tmpdir() as tmp_dir:
            opts = [("OUTPUT", out_file),
                    ("SORT_ORDER", "coordinate"),
                    ("TMP_DIR", tmp_dir)]
            for in_file in in_files:
                opts.append(("INPUT", in_file))
            with file_transaction(out_file):
                picard.run("MergeSamFiles", opts)
    return out_file

def picard_index(picard, in_bam):
    index_file = "%s.bai" % in_bam
    if not os.path.exists(index_file):
        opts = [("INPUT", in_bam),
                ("OUTPUT", index_file)]
        with file_transaction(index_file):
            picard.run("BuildBamIndex", opts)
    return index_file

def picard_index_ref(picard, ref_file):
    """Provide a Picard style dict index file for a reference genome.
    """
    dict_file = "%s.dict" % os.path.splitext(ref_file)[0]
    if not os.path.exists(dict_file):
        opts = [("REFERENCE", ref_file),
                ("OUTPUT", dict_file)]
        with file_transaction(dict_file):
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
    if not (os.path.exists(out_bam) and os.path.getsize(out_bam) > 0):
        with curdir_tmpdir() as tmp_dir:
            opts = [("FASTQ", fastq_one),
                    ("QUALITY_FORMAT", qual_format),
                    ("READ_GROUP_NAME", rg_name),
                    ("SAMPLE_NAME", sample_name),
                    ("PLATFORM_UNIT", pu_name),
                    ("PLATFORM", platform),
                    ("TMP_DIR", tmp_dir),
                    ("OUTPUT", out_bam)]
            if fastq_two:
                opts.append(("FASTQ2", fastq_two))
            with file_transaction(out_bam):
                picard.run("FastqToSam", opts)
    return out_bam

def picard_sam_to_bam(picard, align_sam, fastq_bam, ref_file,
                      is_paired=False):
    """Convert SAM to BAM, including unmapped reads from fastq BAM file.
    """
    out_bam = "%s.bam" % os.path.splitext(align_sam)[0]
    if not os.path.exists(out_bam):
        with curdir_tmpdir() as tmp_dir:
            opts = [("UNMAPPED", fastq_bam),
                    ("ALIGNED", align_sam),
                    ("OUTPUT", out_bam),
                    ("REFERENCE_SEQUENCE", ref_file),
                    ("TMP_DIR", tmp_dir),
                    ("PAIRED_RUN", ("true" if is_paired else "false")),
                    ]
            with file_transaction(out_bam):
                picard.run("MergeBamAlignment", opts)
    return out_bam

def picard_mark_duplicates(picard, align_bam):
    base, ext = os.path.splitext(align_bam)
    base = base.replace(".", "-")
    dup_bam = "%s-dup%s" % (base, ext)
    dup_metrics = "%s-dup.dup_metrics" % base
    if not os.path.exists(dup_bam):
        with curdir_tmpdir() as tmp_dir:
            opts = [("INPUT", align_bam),
                    ("OUTPUT", dup_bam),
                    ("TMP_DIR", tmp_dir),
                    ("METRICS_FILE", dup_metrics)]
        with file_transaction(dup_bam, dup_metrics):
            picard.run("MarkDuplicates", opts)
    return dup_bam, dup_metrics

def picard_fixmate(picard, align_bam):
    """Run Picard's FixMateInformation generating an aligned output file.
    """
    base, ext = os.path.splitext(align_bam)
    out_file = "%s-sort%s" % (base, ext)
    if not os.path.exists(out_file):
        with curdir_tmpdir() as tmp_dir:
            opts = [("INPUT", align_bam),
                    ("OUTPUT", out_file),
                    ("TMP_DIR", tmp_dir),
                    ("SORT_ORDER", "coordinate")]
            with file_transaction(out_file):
                picard.run("FixMateInformation", opts)
    return out_file
