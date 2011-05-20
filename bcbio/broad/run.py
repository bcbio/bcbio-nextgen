"""Convenience functions for running common Picard and GATK utilities.
"""
import os

from bcbio.utils import curdir_tmpdir

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

def picard_index_ref(picard, ref_file):
    """Provide a Picard style dict index file for a reference genome.
    """
    dict_file = "%s.dict" % os.path.splitext(ref_file)[0]
    if not os.path.exists(dict_file):
        opts = [("REFERENCE", ref_file),
                ("OUTPUT", dict_file)]
        picard.run("CreateSequenceDictionary", opts)
    return dict_file

def gatk_realigner_targets(runner, align_bam, ref_file, dbsnp=None):
    """Generate a list of interval regions for realignment around indels.
    """
    out_file = "%s-realign.intervals" % os.path.splitext(align_bam)[0]
    params = ["-T", "RealignerTargetCreator",
              "-I", align_bam,
              "-R", ref_file,
              "-o", out_file,
              "-l", "INFO",
              ]
    if dbsnp:
        params += ["-B:dbsnp,VCF", dbsnp]
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        runner.run_gatk(params)
    return out_file

def gatk_indel_realignment(runner, align_bam, ref_file, intervals):
    """Perform realignment of BAM file in specified regions
    """
    out_file = "%s-realign.bam" % os.path.splitext(align_bam)[0]
    params = ["-T", "IndelRealigner",
              "-I", align_bam,
              "-R", ref_file,
              "-targetIntervals", intervals,
              "-o", out_file,
              "-l", "INFO",
              ]
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        with curdir_tmpdir() as tmp_dir:
            runner.run_gatk(params, tmp_dir)
    return out_file

def gatk_realigner(runner, ref_file, align_bam, dbsnp=None):
    """Realign a BAM file around indels using GATK.
    """
    picard_index_ref(runner, ref_file)
    realign_target_file = gatk_realigner_targets(runner, align_bam,
                                                 ref_file, dbsnp)
    realign_bam = gatk_indel_realignment(runner, align_bam, ref_file,
                                         realign_target_file)
    return realign_bam

