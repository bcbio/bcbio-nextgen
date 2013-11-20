"""Perform quality score recalibration with the GATK toolkit.

Corrects read quality scores post-alignment to provide improved estimates of
error rates based on alignments to the reference genome.

http://www.broadinstitute.org/gsa/wiki/index.php/Base_quality_score_recalibration
"""
import os
import shutil

from bcbio import bam, broad, utils
from bcbio.bam import cram
from bcbio.log import logger
from bcbio.utils import curdir_tmpdir, file_exists
from bcbio.distributed.split import parallel_split_combine
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import (process_bam_by_chromosome,
                                   subset_bam_by_region, write_nochr_reads,
                                   subset_variant_regions)
from bcbio.variation.realign import has_aligned_reads

# ## BAMutil recalibration

def bamutil_dedup_recal_cl(in_file, out_file, data, do_recal):
    """Prepare commandline for running deduplication and recalibration with bamutil.
    http://genome.sph.umich.edu/wiki/BamUtil:_dedup
    """
    raise NotImplementedError("Not functional for piped BAM analysis")
    config = data["config"]
    bam_cmd = config_utils.get_program("bam", config)
    ref_file = data["sam_ref"]
    dbsnp_file = data["genome_resources"]["variation"]["dbsnp"]

    cmd = "{bam_cmd} dedup --in {in_file} --out {out_file} --oneChrom"
    if do_recal:
        cmd += " --recab --refFile {ref_file} --dbsnp {dbsnp_file}"
    return cmd.format(**locals())

# ## GATK recalibration

def prep_recal(data):
    """Perform a GATK recalibration of the sorted aligned BAM, producing recalibrated BAM.
    """
    if data["config"]["algorithm"].get("recalibrate", True) in [True, "gatk"]:
        logger.info("Recalibrating %s with GATK" % str(data["name"]))
        ref_file = data["sam_ref"]
        config = data["config"]
        dbsnp_file = data["genome_resources"]["variation"]["dbsnp"]
        broad_runner = broad.runner_from_config(config)
        platform = config["algorithm"]["platform"]
        broad_runner.run_fn("picard_index_ref", ref_file)
        if config["algorithm"].get("mark_duplicates", True):
            (dup_align_bam, _) = broad_runner.run_fn("picard_mark_duplicates", data["work_bam"])
        else:
            dup_align_bam = data["work_bam"]
        bam.index(dup_align_bam, config)
        intervals = config["algorithm"].get("variant_regions", None)
        data["work_bam"] = dup_align_bam
        data["prep_recal"] = _gatk_base_recalibrator(broad_runner, dup_align_bam, ref_file,
                                                     platform, dbsnp_file, intervals)
    return [[data]]

# ## Identify recalibration information

def _gatk_base_recalibrator(broad_runner, dup_align_bam, ref_file, platform,
        dbsnp_file, intervals):
    """Step 1 of GATK recalibration process, producing table of covariates.

    Large whole genome BAM files take an excessively long time to recalibrate and
    the extra inputs don't help much beyond a certain point. See the 'Downsampling analysis'
    plots in the GATK documentation:

    http://gatkforums.broadinstitute.org/discussion/44/base-quality-score-recalibrator#latest

    This identifies large files and calculates the fraction to downsample to.

    TODO: Use new GATK 2.6+ AnalyzeCovariates tool to plot recalibration results.
    """
    target_counts = 1e8 # 100 million reads per read group, 20x the plotted max
    out_file = "%s.grp" % os.path.splitext(dup_align_bam)[0]
    if not file_exists(out_file):
        if has_aligned_reads(dup_align_bam, intervals):
            with curdir_tmpdir() as tmp_dir:
                with file_transaction(out_file) as tx_out_file:
                    params = ["-T", "BaseRecalibrator",
                              "-o", tx_out_file,
                              "-I", dup_align_bam,
                              "-R", ref_file,
                              ]
                    downsample_pct = bam.get_downsample_pct(broad_runner, dup_align_bam, target_counts)
                    if downsample_pct:
                        params += ["--downsample_to_fraction", str(downsample_pct),
                                   "--downsampling_type", "ALL_READS"]
                    if platform.lower() == "solid":
                        params += ["--solid_nocall_strategy", "PURGE_READ",
                                   "--solid_recal_mode", "SET_Q_ZERO_BASE_N"]
                    # GATK-lite does not have support for
                    # insertion/deletion quality modeling
                    if broad_runner.gatk_type() == "lite":
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
        out_bam = write_nochr_reads(data["work_bam"], out_file, data["config"])
    else:
        out_bam = _run_recal_bam(data["work_bam"], data["prep_recal"],
                                 region, data["sam_ref"], out_file, config)
    qual_bin = config["algorithm"].get("quality_bin", None)
    if ((qual_bin is True or qual_bin == "postrecal" or
         isinstance(qual_bin, list) and "postrecal" in qual_bin)
         and has_aligned_reads(out_bam)):
        binned_bam = cram.illumina_qual_bin(out_bam, data["sam_ref"],
                                         os.path.dirname(out_bam), config)
        shutil.move(out_bam, out_bam + ".binned")
        shutil.move(binned_bam, out_bam)
        utils.save_diskspace(out_bam + ".binned",
                             "Quality binned to %s" % out_bam, config)
    data["work_bam"] = out_bam
    return [data]

def _run_recal_bam(dup_align_bam, recal_file, region, ref_file, out_file, config):
    """Run BAM recalibration with the given input
    """
    if not file_exists(out_file):
        if _recal_available(recal_file):
            broad_runner = broad.runner_from_config(config)
            with curdir_tmpdir() as tmp_dir:
                with file_transaction(out_file) as tx_out_file:
                    params = ["-T", "PrintReads",
                              "-BQSR", recal_file,
                              "-R", ref_file,
                              "-I", dup_align_bam,
                              "--out", tx_out_file,
                              ]
                    base_bed = config["algorithm"].get("variant_regions", None)
                    region_bed = subset_variant_regions(base_bed, region, tx_out_file)
                    if region_bed:
                        params += ["-L", region_bed, "--interval_set_rule", "INTERSECTION"]
                    elif region:
                        params += ["-L", region, "--interval_set_rule", "INTERSECTION"]
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
