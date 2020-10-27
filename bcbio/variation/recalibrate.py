"""Perform quality score recalibration with the GATK toolkit.

Corrects read quality scores post-alignment to provide improved estimates of
error rates based on alignments to the reference genome.

http://www.broadinstitute.org/gsa/wiki/index.php/Base_quality_score_recalibration
"""
import os
from distutils.version import LooseVersion

import toolz as tz

from bcbio import bam, broad, utils
from bcbio.log import logger
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.variation.realign import has_aligned_reads
from bcbio.variation import sentieon

def prep_recal(data):
    """Do pre-BQSR recalibration, calculation of recalibration tables.
    """
    if dd.get_recalibrate(data) in [True, "gatk"]:
        logger.info("Prepare BQSR tables with GATK: %s " % str(dd.get_sample_name(data)))
        dbsnp_file = tz.get_in(("genome_resources", "variation", "dbsnp"), data)
        if not dbsnp_file:
            logger.info("Skipping GATK BaseRecalibrator because no VCF file of known variants was found.")
            return data
        broad_runner = broad.runner_from_config(data["config"])
        data["prep_recal"] = _gatk_base_recalibrator(broad_runner, dd.get_align_bam(data),
                                                     dd.get_ref_file(data), dd.get_platform(data),
                                                     dbsnp_file,
                                                     dd.get_variant_regions(data) or dd.get_sample_callable(data),
                                                     data)
    elif dd.get_recalibrate(data) == "sentieon":
        logger.info("Prepare BQSR tables with sentieon: %s " % str(dd.get_sample_name(data)))
        data["prep_recal"] = sentieon.bqsr_table(data)
    elif dd.get_recalibrate(data):
        raise NotImplementedError("Unsupported recalibration type: %s" % (dd.get_recalibrate(data)))
    return data

def apply_recal(data):
    """Apply recalibration tables to the sorted aligned BAM, producing recalibrated BAM.
    """
    orig_bam = dd.get_align_bam(data) or dd.get_work_bam(data)
    had_work_bam = "work_bam" in data
    if dd.get_recalibrate(data) in [True, "gatk"]:
        if data.get("prep_recal"):
            logger.info("Applying BQSR recalibration with GATK: %s " % str(dd.get_sample_name(data)))
            data["work_bam"] = _gatk_apply_bqsr(data)
    elif dd.get_recalibrate(data) == "sentieon":
        if data.get("prep_recal"):
            logger.info("Applying BQSR recalibration with sentieon: %s " % str(dd.get_sample_name(data)))
            data["work_bam"] = sentieon.apply_bqsr(data)
    elif dd.get_recalibrate(data):
        raise NotImplementedError("Unsupported recalibration type: %s" % (dd.get_recalibrate(data)))
    # CWL does not have work/alignment BAM separation
    if not had_work_bam and dd.get_work_bam(data):
        data["align_bam"] = dd.get_work_bam(data)
    if orig_bam != dd.get_work_bam(data) and orig_bam != dd.get_align_bam(data):
        utils.save_diskspace(orig_bam, "BAM recalibrated to %s" % dd.get_work_bam(data), data["config"])
    return data

# ## GATK recalibration

def _gatk_base_recalibrator(broad_runner, dup_align_bam, ref_file, platform,
                            dbsnp_file, intervals, data):
    """Step 1 of GATK recalibration process, producing table of covariates.

    For GATK 4 we use local multicore spark runs:
    https://github.com/broadinstitute/gatk/issues/2345

    For GATK3, Large whole genome BAM files take an excessively long time to recalibrate and
    the extra inputs don't help much beyond a certain point. See the 'Downsampling analysis'
    plots in the GATK documentation:

    http://gatkforums.broadinstitute.org/discussion/44/base-quality-score-recalibrator#latest

    This identifies large files and calculates the fraction to downsample to.

    spark host and timeout settings help deal with runs on restricted systems
    where we encounter network and timeout errors
    """
    target_counts = 1e8  # 100 million reads per read group, 20x the plotted max
    out_file = os.path.join(dd.get_work_dir(data), "align", dd.get_sample_name(data),
                            "%s-recal.grp" % utils.splitext_plus(os.path.basename(dup_align_bam))[0])
    if not utils.file_exists(out_file):
        if has_aligned_reads(dup_align_bam, intervals):
            with file_transaction(data, out_file) as tx_out_file:
                gatk_type = broad_runner.gatk_type()
                assert gatk_type in ["restricted", "gatk4"], \
                    "Require full version of GATK 2.4+ or GATK4 for BQSR"
                params = ["-I", dup_align_bam]
                cores = dd.get_num_cores(data)
                if gatk_type == "gatk4":
                    resources = config_utils.get_resources("gatk-spark", data["config"])
                    spark_opts = [str(x) for x in resources.get("options", [])]
                    params += ["-T", "BaseRecalibratorSpark",
                               "--output", tx_out_file, "--reference", dd.get_ref_file(data)]
                    if spark_opts:
                        params += spark_opts
                    else:
                        params += ["--spark-master", "local[%s]" % cores,
                                   "--conf", "spark.driver.host=localhost", "--conf", "spark.network.timeout=800",
                                   "--conf", "spark.executor.heartbeatInterval=100",
                                   "--conf", "spark.local.dir=%s" % os.path.dirname(tx_out_file)]
                    if dbsnp_file:
                        params += ["--known-sites", dbsnp_file]
                    if intervals:
                        params += ["-L", intervals, "--interval-set-rule", "INTERSECTION"]
                else:
                    params += ["-T", "BaseRecalibrator",
                                "-o", tx_out_file, "-R", ref_file]
                    downsample_pct = bam.get_downsample_pct(dup_align_bam, target_counts, data)
                    if downsample_pct:
                        params += ["--downsample_to_fraction", str(downsample_pct),
                                   "--downsampling_type", "ALL_READS"]
                    if platform.lower() == "solid":
                        params += ["--solid_nocall_strategy", "PURGE_READ",
                                   "--solid_recal_mode", "SET_Q_ZERO_BASE_N"]
                    if dbsnp_file:
                        params += ["--knownSites", dbsnp_file]
                    if intervals:
                        params += ["-L", intervals, "--interval_set_rule", "INTERSECTION"]
                memscale = {"magnitude": 0.9 * cores, "direction": "increase"} if cores > 1 else None
                broad_runner.run_gatk(params, os.path.dirname(tx_out_file), memscale=memscale,
                                      parallel_gc=True)
        else:
            with open(out_file, "w") as out_handle:
                out_handle.write("# No aligned reads")
    return out_file

def _gatk_apply_bqsr(data):
    """Parallel BQSR support for GATK4.

    Normalized qualities to 3 bin outputs at 10, 20 and 30 based on pipeline standard
    recommendations, which will help with output file sizes:
    https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md#base-quality-score-binning-scheme
    https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/5585cdf7877104f2c61b2720ddfe7235f2fad577/PairedEndSingleSampleWf.gatk4.0.wdl#L1081

    spark host and timeout settings help deal with runs on restricted systems
    where we encounter network and timeout errors
    """
    in_file = dd.get_align_bam(data) or dd.get_work_bam(data)
    out_file = os.path.join(dd.get_work_dir(data), "align", dd.get_sample_name(data),
                            "%s-recal.bam" % utils.splitext_plus(os.path.basename(in_file))[0])
    if not utils.file_uptodate(out_file, in_file):
        with file_transaction(data, out_file) as tx_out_file:
            broad_runner = broad.runner_from_config(data["config"])
            gatk_type = broad_runner.gatk_type()
            cores = dd.get_num_cores(data)
            if gatk_type == "gatk4":
                resources = config_utils.get_resources("gatk-spark", data["config"])
                spark_opts = [str(x) for x in resources.get("options", [])]
                params = ["-T", "ApplyBQSRSpark",
                          "--input", in_file, "--output", tx_out_file, "--bqsr-recal-file", data["prep_recal"],
                          "--static-quantized-quals", "10", "--static-quantized-quals", "20",
                          "--static-quantized-quals", "30"]
                if spark_opts:
                    params += spark_opts
                else:
                    params += ["--spark-master", "local[%s]" % cores,
                               "--conf", "spark.local.dir=%s" % os.path.dirname(tx_out_file),
                               "--conf", "spark.driver.host=localhost", "--conf", "spark.network.timeout=800"]
                # Avoid problems with StreamClosedErrors on GATK 4.1+
                # https://github.com/bcbio/bcbio-nextgen/issues/2806#issuecomment-492504497
                params += ["--create-output-bam-index", "false"]
            else:
                params = ["-T", "PrintReads", "-R", dd.get_ref_file(data), "-I", in_file,
                          "-BQSR", data["prep_recal"], "-o", tx_out_file]
            # Avoid problems with intel deflater for GATK 3.8 and GATK4
            # https://github.com/bcbio/bcbio-nextgen/issues/2145#issuecomment-343095357
            if gatk_type == "gatk4":
                params += ["--jdk-deflater", "--jdk-inflater"]
            elif LooseVersion(broad_runner.gatk_major_version()) > LooseVersion("3.7"):
                params += ["-jdk_deflater", "-jdk_inflater"]
            memscale = {"magnitude": 0.9 * cores, "direction": "increase"} if cores > 1 else None
            broad_runner.run_gatk(params, os.path.dirname(tx_out_file), memscale=memscale,
                                  parallel_gc=True)
    bam.index(out_file, data["config"])
    return out_file
