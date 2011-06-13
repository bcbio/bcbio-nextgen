"""Convenience functions for running common GATK utilities.
"""
import os

from bcbio.utils import curdir_tmpdir

def gatk_realigner_targets(runner, align_bam, ref_file, dbsnp=None,
                           deep_coverage=False):
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
    if deep_coverage:
        params += ["--mismatchFraction", "0.30",
                   "--maxIntervalSize", "650"]
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        runner.run_gatk(params)
    return out_file

def gatk_indel_realignment(runner, align_bam, ref_file, intervals,
                           deep_coverage=False):
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
    if deep_coverage:
        params += ["--maxReadsInMemory", "300000",
                   "--maxReadsForRealignment", str(int(5e5)),
                   "--maxReadsForConsensuses", "500",
                   "--maxConsensuses", "100"]
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        with curdir_tmpdir() as tmp_dir:
            runner.run_gatk(params, tmp_dir)
    return out_file

def gatk_realigner(runner, align_bam, ref_file, dbsnp=None,
                   deep_coverage=False):
    """Realign a BAM file around indels using GATK.
    """
    picard_index_ref(runner, ref_file)
    realign_target_file = gatk_realigner_targets(runner, align_bam,
                                                 ref_file, dbsnp, deep_coverage)
    realign_bam = gatk_indel_realignment(runner, align_bam, ref_file,
                                         realign_target_file, deep_coverage)
    return realign_bam
