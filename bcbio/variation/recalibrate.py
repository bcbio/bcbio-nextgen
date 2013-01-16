"""Perform quality score recalibration with the GATK toolkit.

Corrects read quality scores post-alignment to provide improved estimates of
error rates based on alignments to the reference genome.

http://www.broadinstitute.org/gsa/wiki/index.php/Base_quality_score_recalibration
"""
import os
import shutil
from contextlib import closing

import pysam

from bcbio import broad
from bcbio.log import logger
from bcbio.utils import curdir_tmpdir, file_exists
from bcbio.distributed.split import parallel_split_combine
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline.shared import (configured_ref_file, process_bam_by_chromosome,
                                   subset_bam_by_region, write_nochr_reads)
from bcbio.variation.realign import has_aligned_reads

def prep_recal(data):
    """Perform a GATK recalibration of the sorted aligned BAM, producing recalibrated BAM.
    """
    if data["config"]["algorithm"].get("recalibrate", True):
        logger.info("Recalibrating %s with GATK" % str(data["name"]))
        ref_file = data["sam_ref"]
        config = data["config"]
        dbsnp_file = configured_ref_file("dbsnp", config, ref_file)
        broad_runner = broad.runner_from_config(config)
        platform = config["algorithm"]["platform"]
        broad_runner.run_fn("picard_index_ref", ref_file)
        if config["algorithm"].get("mark_duplicates", True):
            (dup_align_bam, _) = broad_runner.run_fn("picard_mark_duplicates", data["work_bam"],
                                                     remove_dups=True)
        else:
            dup_align_bam = data["work_bam"]
        broad_runner.run_fn("picard_index", dup_align_bam)
        intervals = config["algorithm"].get("variant_regions", None)
        data["work_bam"] = dup_align_bam
        data["prep_recal"] = _gatk_base_recalibrator(broad_runner, dup_align_bam, ref_file,
                                                     platform, dbsnp_file, intervals)
    return [[data]]

# ## Identify recalibration information

def _get_downsample_pct(runner, in_bam):
    """Calculate a downsampling percent to use for large BAM files.

    Large whole genome BAM files take an excessively long time to recalibrate and
    the extra inputs don't help much beyond a certain point. See the 'Downsampling analysis'
    plots in the GATK documentation:

    http://gatkforums.broadinstitute.org/discussion/44/base-quality-score-recalibrator#latest

    This identifies large files and calculates the fraction to downsample to.
    """
    target_counts = 1e8 # 100 million reads per read group, 20x the plotted max
    total = sum(x.aligned for x in runner.run_fn("picard_idxstats", in_bam))
    with closing(pysam.Samfile(in_bam, "rb")) as work_bam:
        n_rgs = max(1, len(work_bam.header["RG"]))
    rg_target = n_rgs * target_counts
    if total > rg_target:
        return float(rg_target) / float(total)

def _gatk_base_recalibrator(broad_runner, dup_align_bam, ref_file, platform,
        dbsnp_file, intervals):
    """Step 1 of GATK recalibration process, producing table of covariates.
    """
    out_file = "%s.grp" % os.path.splitext(dup_align_bam)[0]
    plot_file = "%s-plots.pdf" % os.path.splitext(dup_align_bam)[0]
    if not file_exists(out_file):
        if has_aligned_reads(dup_align_bam):
            with curdir_tmpdir() as tmp_dir:
                with file_transaction(out_file) as tx_out_file:
                    params = ["-T", "BaseRecalibrator",
                              "-o", tx_out_file,
                              "--plot_pdf_file", plot_file,
                              "-I", dup_align_bam,
                              "-R", ref_file,
                              ]
                    downsample_pct = _get_downsample_pct(broad_runner, dup_align_bam)
                    if downsample_pct:
                        params += ["--downsample_to_fraction", str(downsample_pct),
                                   "--downsampling_type", "ALL_READS"]
                    # GATK-lite does not have support for
                    # insertion/deletion quality modeling
                    if not broad_runner.has_gatk_full():
                        params += ["--disable_indel_quals"]
                    if dbsnp_file:
                        params += ["--knownSites", dbsnp_file]
                    if intervals:
                        params += ["-L", intervals, "--interval_set_rule", "INTERSECTION"]
                    broad_runner.run_gatk(params, tmp_dir)
        else:
            with open(out_file, "w") as out_handle:
                out_handle.write("# No aligned reads")
    return out_file

# ## Create recalibrated BAM

def parallel_write_recal_bam(xs, parallel_fn):
    """Rewrite a recalibrated BAM file in parallel, working off each chromosome.
    """
    to_process = []
    finished = []
    for x in xs:
        if x[0]["config"]["algorithm"].get("recalibrate", True):
            to_process.append(x)
        else:
            finished.append(x)
    if len(to_process) > 0:
        file_key = "work_bam"
        split_fn = process_bam_by_chromosome("-gatkrecal.bam", file_key,
                                           default_targets=["nochr"])
        processed = parallel_split_combine(to_process, split_fn, parallel_fn,
                                           "write_recal_bam", "combine_bam",
                                           file_key, ["config"])
        finished.extend(processed)
        # Save diskspace from original to recalibrated
        #save_diskspace(data["work_bam"], "Recalibrated to %s" % recal_bam,
        #               data["config"])
    return finished

def write_recal_bam(data, region=None, out_file=None):
    """Step 2 of GATK recalibration -- use covariates to re-write output file.
    """
    config = data["config"]
    if out_file is None:
        out_file = "%s-gatkrecal.bam" % os.path.splitext(data["work_bam"])[0]
    logger.info("Writing recalibrated BAM for %s to %s" % (data["name"], out_file))
    if region == "nochr":
        out_bam = write_nochr_reads(data["work_bam"], out_file)
    else:
        out_bam = _run_recal_bam(data["work_bam"], data["prep_recal"],
                                 region, data["sam_ref"], out_file, config)
    data["work_bam"] = out_bam
    return [data]

def _run_recal_bam(dup_align_bam, recal_file, region, ref_file, out_file, config):
    """Run BAM recalibration with the given input
    """
    if not file_exists(out_file):
        if _recal_available(recal_file):
            broad_runner = broad.runner_from_config(config)
            intervals = config["algorithm"].get("variant_regions", None)
            with curdir_tmpdir() as tmp_dir:
                with file_transaction(out_file) as tx_out_file:
                    params = ["-T", "PrintReads",
                              "-BQSR", recal_file,
                              "-R", ref_file,
                              "-I", dup_align_bam,
                              "--out", tx_out_file,
                              ]
                    if region:
                        params += ["-L", region]
                    if intervals:
                        params += ["-L", intervals]
                    if params and intervals:
                        params += ["--interval_set_rule", "INTERSECTION"]
                    broad_runner.run_gatk(params, tmp_dir)
        elif region:
            subset_bam_by_region(dup_align_bam, region, out_file)
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
