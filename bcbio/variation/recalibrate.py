"""Perform quality score recalibration with the GATK toolkit.

Corrects read quality scores post-alignment to provide improved estimates of
error rates based on alignments to the reference genome.

http://www.broadinstitute.org/gsa/wiki/index.php/Base_quality_score_recalibration
"""
import os

import toolz as tz

from bcbio import bam, broad
from bcbio.log import logger
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.pipeline import datadict as dd
from bcbio.variation.realign import has_aligned_reads

# ## GATK recalibration

def prep_recal(data):
    """Perform a GATK recalibration of the sorted aligned BAM, producing recalibrated BAM.
    """
    if dd.get_recalibrate(data) in [True, "gatk"]:
        logger.info("Recalibrating %s with GATK" % str(dd.get_sample_name(data)))
        dbsnp_file = tz.get_in(("genome_resources", "variation", "dbsnp"), data)
        if not dbsnp_file:
            logger.info("Skipping GATK BaseRecalibrator because no VCF file of known variants was found.")
            return [[data]]
        broad_runner = broad.runner_from_config(data["config"])
        data["prep_recal"] = _gatk_base_recalibrator(broad_runner, dd.get_align_bam(data),
                                                     dd.get_ref_file(data), dd.get_platform(data),
                                                     dbsnp_file, dd.get_variant_regions(data), data)
    return data

# ## Identify recalibration information

def _gatk_base_recalibrator(broad_runner, dup_align_bam, ref_file, platform,
                            dbsnp_file, intervals, data):
    """Step 1 of GATK recalibration process, producing table of covariates.

    Large whole genome BAM files take an excessively long time to recalibrate and
    the extra inputs don't help much beyond a certain point. See the 'Downsampling analysis'
    plots in the GATK documentation:

    http://gatkforums.broadinstitute.org/discussion/44/base-quality-score-recalibrator#latest

    This identifies large files and calculates the fraction to downsample to.
    """
    target_counts = 1e8  # 100 million reads per read group, 20x the plotted max
    out_file = "%s.grp" % os.path.splitext(dup_align_bam)[0]
    if not file_exists(out_file):
        if has_aligned_reads(dup_align_bam, intervals):
            with tx_tmpdir(data) as tmp_dir:
                with file_transaction(data, out_file) as tx_out_file:
                    gatk_type = broad_runner.gatk_type()
                    assert gatk_type in ["restricted", "gatk4"], \
                        "Require full version of GATK 2.4+, or GATK4 for BQSR"
                    params = ["-T", "BaseRecalibrator",
                              "-I", dup_align_bam,
                              "-R", ref_file,
                              ]
                    downsample_pct = bam.get_downsample_pct(dup_align_bam, target_counts, data)
                    if downsample_pct:
                        params += ["--downsample_to_fraction", str(downsample_pct),
                                   "--downsampling_type", "ALL_READS"]
                    if gatk_type == "gatk4":
                        params += ["--output", tx_out_file]
                    else:
                        params += ["-o", tx_out_file]
                        if platform.lower() == "solid":
                            params += ["--solid_nocall_strategy", "PURGE_READ",
                                    "--solid_recal_mode", "SET_Q_ZERO_BASE_N"]
                    if dbsnp_file:
                        params += ["--knownSites", dbsnp_file]
                    if intervals:
                        params += ["-L", intervals, "--interval_set_rule", "INTERSECTION"]
                    broad_runner.run_gatk(params, tmp_dir)
        else:
            with open(out_file, "w") as out_handle:
                out_handle.write("# No aligned reads")
    return out_file
