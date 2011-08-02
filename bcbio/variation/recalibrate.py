"""Perform quality score recalibration with the GATK toolkit.

Read quality scores can be corrected post-alignment to provide better estimates of
actual error rates based on alignments to the reference genome.

http://www.broadinstitute.org/gsa/wiki/index.php/Base_quality_score_recalibration
"""
import os
import shutil

from bcbio import broad
from bcbio.utils import curdir_tmpdir, file_transaction

def gatk_recalibrate(align_bam, ref_file, config, snp_file=None):
    """Perform a GATK recalibration of the sorted aligned BAM, producing recalibrated BAM.
    """
    picard = broad.runner_from_config(config)
    platform = config["algorithm"]["platform"]
    picard.run_fn("picard_index_ref", ref_file)
    (dup_align_bam, _) = picard.run_fn("picard_mark_duplicates", align_bam)
    recal_file = _gatk_count_covariates(picard, dup_align_bam, ref_file, platform,
            snp_file)
    return _gatk_table_recalibrate(picard, dup_align_bam, ref_file, recal_file,
                                   platform)

def _gatk_table_recalibrate(picard, dup_align_bam, ref_file, recal_file, platform):
    """Step 2 of GATK recalibration -- use covariates to re-write output file.
    """
    out_file = "%s-gatkrecal.bam" % os.path.splitext(dup_align_bam)[0]
    params = ["-T", "TableRecalibration",
              "-recalFile", recal_file,
              "-R", ref_file,
              "-I", dup_align_bam,
              "--out", out_file,
              "-baq",  "RECALCULATE",
              "-l", "INFO",
              "-U",
              "-OQ",
              "--default_platform", platform,
              ]
    if not os.path.exists(out_file):
        if _recal_available(recal_file):
            with curdir_tmpdir() as tmp_dir:
                with file_transaction(out_file):
                    picard.run_gatk(params, tmp_dir)
        else:
            shutil.copy(dup_align_bam, out_file)
    return out_file

def _recal_available(recal_file):
    """Determine if it's possible to do a recalibration; do we have data?
    """
    if os.path.exists(recal_file):
        with open(recal_file) as in_handle:
            while 1:
                line = in_handle.next()
                if not line.startswith("#"):
                    break
            test_line = in_handle.next()
            if test_line and not test_line.startswith("EOF"):
                return True
    return False

def _gatk_count_covariates(picard, dup_align_bam, ref_file, platform,
        snp_file):
    """Step 1 of GATK recalibration process -- counting covariates.
    """
    out_file = "%s.recal" % os.path.splitext(dup_align_bam)[0]
    params = ["-T", "CountCovariates",
              "-cov", "ReadGroupCovariate",
              "-cov", "QualityScoreCovariate",
              "-cov", "CycleCovariate",
              "-cov", "DinucCovariate",
              "-recalFile", out_file,
              "-I", dup_align_bam,
              "-R", ref_file,
              "-l", "INFO",
              "-U",
              "-OQ",
              "--default_platform", platform,
              ]
    if snp_file:
        params += ["-B:dbsnp,VCF", snp_file]
    if not os.path.exists(out_file):
        with curdir_tmpdir() as tmp_dir:
            with file_transaction(out_file):
                picard.run_gatk(params, tmp_dir)
    return out_file

