"""Perform realignment of BAM files around indels using the GATK toolkit.
"""
import os

import pysam

from bcbio import broad
from bcbio.utils import curdir_tmpdir, file_transaction

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
        with file_transaction(out_file):
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
            with file_transaction(out_file):
                runner.run_gatk(params, tmp_dir)
    return out_file

def gatk_realigner(align_bam, ref_file, config, dbsnp=None,
                   deep_coverage=False):
    """Realign a BAM file around indels using GATK, returning sorted BAM.
    """
    runner = broad.runner_from_config(config)
    runner.run_fn("picard_index", align_bam)
    runner.run_fn("picard_index_ref", ref_file)
    if not os.path.exists("%s.fai" % ref_file):
        pysam.faidx(ref_file)
    realign_target_file = gatk_realigner_targets(runner, align_bam,
                                                 ref_file, dbsnp, deep_coverage)
    realign_bam = gatk_indel_realignment(runner, align_bam, ref_file,
                                         realign_target_file, deep_coverage)
    # No longer required in recent GATK (> Feb 2011) -- now done on the fly
    # realign_sort_bam = runner.run_fn("picard_fixmate", realign_bam)
    return realign_bam
