"""Convenience functions for running common Picard and GATK utilities.
"""
import os

from bcbio.picard.utils import curdir_tmpdir

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
            picard.run("SortSam", opts)
    return out_file

def picard_merge(picard, in_files):
    """Merge multiple BAM files together with Picard.
    """
    out_file = "%smerge.bam" % os.path.commonprefix(in_files)
    if not os.path.exists(out_file):
        with curdir_tmpdir() as tmp_dir:
            opts = [("OUTPUT", out_file),
                    ("SORT_ORDER", "coordinate"),
                    ("TMP_DIR", tmp_dir)]
            for in_file in in_files:
                opts.append(("INPUT", in_file))
            picard.run("MergeSamFiles", opts)
    return out_file

def picard_index(picard, in_bam):
    index_file = "%s.bai" % in_bam
    if not os.path.exists(index_file):
        opts = [("INPUT", in_bam),
                ("OUTPUT", index_file)]
        picard.run("BuildBamIndex", opts)
    return index_file
