"""Perform quality score recalibration with the GATK toolkit.

Corrects read quality scores post-alignment to provide improved estimates of
error rates based on alignments to the reference genome.

http://www.broadinstitute.org/gsa/wiki/index.php/Base_quality_score_recalibration
"""
import os
import shutil

from bcbio import broad
from bcbio.utils import curdir_tmpdir, file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.variation.realign import has_aligned_reads

def gatk_recalibrate(align_bam, ref_file, config, snp_file=None):
    """Perform a GATK recalibration of the sorted aligned BAM, producing recalibrated BAM.
    """
    broad_runner = broad.runner_from_config(config)
    platform = config["algorithm"]["platform"]
    broad_runner.run_fn("picard_index_ref", ref_file)
    (dup_align_bam, _) = broad_runner.run_fn("picard_mark_duplicates", align_bam, remove_dups=True)
    broad_runner.run_fn("picard_index", dup_align_bam)
    intervals = config["algorithm"].get("variant_regions", None)
    recal_file = _gatk_base_recalibrator(broad_runner, dup_align_bam, ref_file, platform,
                                         snp_file, intervals)
    recal_bam = _gatk_recalibrate(broad_runner, dup_align_bam, ref_file, recal_file,
                                  platform, intervals)
    broad_runner.run_fn("picard_index", recal_bam)
    return recal_bam

def _gatk_recalibrate(broad_runner, dup_align_bam, ref_file, recal_file,
                      platform, intervals):
    """Step 2 of GATK recalibration -- use covariates to re-write output file.
    """
    out_file = "%s-gatkrecal.bam" % os.path.splitext(dup_align_bam)[0]
    if not file_exists(out_file):
        if _recal_available(recal_file):
            with curdir_tmpdir() as tmp_dir:
                with file_transaction(out_file) as tx_out_file:
                    params = ["-T", "PrintReads",
                              "-BQSR", recal_file,
                              "-R", ref_file,
                              "-I", dup_align_bam,
                              "--out", tx_out_file,
                              ]
                    if intervals:
                        params += ["-L", intervals, "--interval_set_rule", "INTERSECTION"]
                    broad_runner.run_gatk(params, tmp_dir)
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

def _gatk_base_recalibrator(broad_runner, dup_align_bam, ref_file, platform,
        snp_file, intervals):
    """Step 1 of GATK recalibration process, producing table of covariates.
    """
    out_file = "%s.grp" % os.path.splitext(dup_align_bam)[0]
    if not file_exists(out_file):
        if has_aligned_reads(dup_align_bam):
            with curdir_tmpdir() as tmp_dir:
                with file_transaction(out_file) as tx_out_file:
                    params = ["-T", "BaseRecalibrator",
                              "-o", tx_out_file,
                              "-I", dup_align_bam,
                              "-R", ref_file,
                              ]
                    # GATK-lite does not have support for
                    # insertion/deletion quality modeling
                    if not broad_runner.has_gatk_full():
                        params += ["--disable_indel_quals"]
                    if snp_file:
                        params += ["--knownSites", snp_file]
                    if intervals:
                        params += ["-L", intervals, "--interval_set_rule", "INTERSECTION"]
                    broad_runner.run_gatk(params, tmp_dir)
        else:
            with open(out_file, "w") as out_handle:
                out_handle.write("# No aligned reads")
    return out_file

