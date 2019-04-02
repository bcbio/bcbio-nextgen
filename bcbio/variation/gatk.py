"""GATK variant calling -- HaplotypeCaller and UnifiedGenotyper.
"""
import os
from distutils.version import LooseVersion
import shutil
import subprocess

import toolz as tz

from bcbio import bam, broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.pipeline import datadict as dd
from bcbio.variation import annotation, bamprep, bedutils, ploidy, vcfutils

def standard_cl_params(items):
    """Shared command line parameters for GATK programs.

    Handles no removal of duplicate reads for amplicon or
    non mark duplicate experiments. If we have pre-aligned inputs we
    ignore the value or mark duplicates (since they may already be
    marked in the input BAM).
    """
    out = []
    def _skip_duplicates(data):
        return (dd.get_coverage_interval(data) == "amplicon" or
                (dd.get_aligner(data) and not dd.get_mark_duplicates(data)))
    if any(_skip_duplicates(d) for d in items):
        broad_runner = broad.runner_from_config(items[0]["config"])
        gatk_type = broad_runner.gatk_type()
        if gatk_type == "gatk4":
            out += ["--disable-read-filter", "NotDuplicateReadFilter"]
        elif LooseVersion(broad_runner.gatk_major_version()) >= LooseVersion("3.5"):
            out += ["-drf", "DuplicateRead"]
    return out

def _shared_gatk_call_prep(align_bams, items, ref_file, region, out_file, num_cores=1):
    """Shared preparation work for GATK variant calling.
    """
    data = items[0]
    config = data["config"]
    broad_runner = broad.runner_from_config(config)
    gatk_type = broad_runner.gatk_type()
    for x in align_bams:
        bam.index(x, config)
    picard_runner = broad.runner_from_path("picard", config)
    picard_runner.run_fn("picard_index_ref", ref_file)
    params = ["-R", ref_file]
    coverage_depth_min = tz.get_in(["algorithm", "coverage_depth_min"], config)
    if coverage_depth_min and coverage_depth_min < 4:
        confidence = "4.0"
        params += ["--standard_min_confidence_threshold_for_calling", confidence]
    for a in annotation.get_gatk_annotations(config):
        params += ["--annotation", a]
    for x in align_bams:
        params += ["-I", x]
    variant_regions = bedutils.population_variant_regions(items)
    region = subset_variant_regions(variant_regions, region, out_file, items)
    if region:
        if gatk_type == "gatk4":
            params += ["-L", bamprep.region_to_gatk(region), "--interval-set-rule", "INTERSECTION"]
        else:
            params += ["-L", bamprep.region_to_gatk(region), "--interval_set_rule", "INTERSECTION"]
    params += standard_cl_params(items)
    return broad_runner, params

def unified_genotyper(align_bams, items, ref_file, assoc_files,
                       region=None, out_file=None):
    """Perform SNP genotyping on the given alignment file.
    """
    if out_file is None:
        out_file = "%s-variants.vcf.gz" % utils.splitext_plus(align_bams[0])[0]
    if not utils.file_exists(out_file):
        broad_runner, params = \
            _shared_gatk_call_prep(align_bams, items, ref_file, region, out_file)
        with file_transaction(items[0], out_file) as tx_out_file:
            params += ["-T", "UnifiedGenotyper",
                       "-o", tx_out_file,
                       "-ploidy", (str(ploidy.get_ploidy(items, region))
                                   if broad_runner.gatk_type() == "restricted" else "2"),
                       "--genotype_likelihoods_model", "BOTH"]
            resources = config_utils.get_resources("gatk", items[0]["config"])
            if "options" in resources:
                params += [str(x) for x in resources.get("options", [])]
            broad_runner.run_gatk(params)
    return vcfutils.bgzip_and_index(out_file, items[0]["config"])

def _joint_calling(items):
    """Determine if this call feeds downstream into joint calls.
    """
    jointcaller = tz.get_in(("config", "algorithm", "jointcaller"), items[0])
    if jointcaller:
        assert len(items) == 1, "Can only do joint calling preparation with GATK with single samples"
        assert tz.get_in(("metadata", "batch"), items[0]) is not None, \
            "Joint calling requires batched samples, %s has no metadata batch." % dd.get_sample_name(items[0])
    return jointcaller

def _use_spark(num_cores, gatk_type, items, opts):
    return ((len(items) == 1 and num_cores > 1 and gatk_type == "gatk4") or
            "--spark-master" in opts)

def haplotype_caller(align_bams, items, ref_file, assoc_files,
                       region=None, out_file=None):
    """Call variation with GATK's HaplotypeCaller.

    This requires the full non open-source version of GATK.
    """
    if out_file is None:
        out_file = "%s-variants.vcf.gz" % utils.splitext_plus(align_bams[0])[0]
    if not utils.file_exists(out_file):
        num_cores = dd.get_num_cores(items[0])
        broad_runner, params = \
            _shared_gatk_call_prep(align_bams, items, ref_file, region, out_file, num_cores)
        gatk_type = broad_runner.gatk_type()
        assert gatk_type in ["restricted", "gatk4"], \
            "Require full version of GATK 2.4+, or GATK4 for haplotype calling"
        with file_transaction(items[0], out_file) as tx_out_file:
            resources = config_utils.get_resources("gatk-spark", items[0]["config"])
            spark_opts = [str(x) for x in resources.get("options", [])]
            if _use_spark(num_cores, gatk_type, items, spark_opts):
                params += ["-T", "HaplotypeCallerSpark"]
                if spark_opts:
                    params += spark_opts
                else:
                    params += ["--spark-master", "local[%s]" % num_cores,
                               "--conf", "spark.local.dir=%s" % os.path.dirname(tx_out_file),
                               "--conf", "spark.driver.host=localhost", "--conf", "spark.network.timeout=800",
                               "--conf", "spark.executor.heartbeatInterval=100"]
            else:
                params += ["-T", "HaplotypeCaller"]
            params += ["--annotation", "ClippingRankSumTest",
                       "--annotation", "DepthPerSampleHC"]
            # Enable hardware based optimizations in GATK 3.1+
            if LooseVersion(broad_runner.gatk_major_version()) >= LooseVersion("3.1"):
                if _supports_avx():
                    # Scale down HMM thread default to avoid overuse of cores
                    # https://github.com/bcbio/bcbio-nextgen/issues/2442
                    if gatk_type == "gatk4":
                        params += ["--native-pair-hmm-threads", "1"]
                    # GATK4 selects the right HMM optimization automatically with FASTEST_AVAILABLE
                    # GATK3 needs to be explicitly set
                    else:
                        params += ["--pair_hmm_implementation", "VECTOR_LOGLESS_CACHING"]
            resources = config_utils.get_resources("gatk-haplotype", items[0]["config"])
            if "options" in resources:
                params += [str(x) for x in resources.get("options", [])]
            # Prepare gVCFs if doing joint calling
            is_joint = False
            if _joint_calling(items) or any("gvcf" in dd.get_tools_on(d) for d in items):
                is_joint = True
                # If joint calling parameters not set in user options
                if not any([x in ["--emit-ref-confidence", "-ERC", "--emitRefConfidence"] for x in params]):
                    if gatk_type == "gatk4":
                        params += ["--emit-ref-confidence", "GVCF"]
                    else:
                        params += ["--emitRefConfidence", "GVCF"]
                        params += ["--variant_index_type", "LINEAR", "--variant_index_parameter", "128000"]
                # Set GQ banding to not be single GQ resolution
                # No recommended default but try to balance resolution and size
                # http://gatkforums.broadinstitute.org/gatk/discussion/7051/recommendation-best-practices-gvcf-gq-bands

                if not any([x in ["-GQB"] for x in params]):
                    for boundary in [10, 20, 30, 40, 60, 80]:
                        params += ["-GQB", str(boundary)]
            # Enable non-diploid calling in GATK 3.3+
            if LooseVersion(broad_runner.gatk_major_version()) >= LooseVersion("3.3"):
                params += ["-ploidy", str(ploidy.get_ploidy(items, region))]
            if gatk_type == "gatk4":
                # GATK4 Spark calling does not support bgzipped output, use plain VCFs
                if is_joint and _use_spark(num_cores, gatk_type, items, spark_opts):
                    tx_out_file = tx_out_file.replace(".vcf.gz", ".vcf")
                params += ["--output", tx_out_file]
            else:
                params += ["-o", tx_out_file]
            broad_runner.new_resources("gatk-haplotype")
            memscale = {"magnitude": 0.9 * num_cores, "direction": "increase"} if num_cores > 1 else None
            try:
                broad_runner.run_gatk(params, os.path.dirname(tx_out_file), memscale=memscale,
                                      parallel_gc=_use_spark(num_cores, gatk_type, items, spark_opts))
            except subprocess.CalledProcessError as msg:
                # Spark failing on regions without any reads, write an empty VCF instead
                # https://github.com/broadinstitute/gatk/issues/4234
                if (_use_spark(num_cores, gatk_type, items, spark_opts) and
                      str(msg).find("java.lang.UnsupportedOperationException: empty collection") >= 0 and
                      str(msg).find("at org.apache.spark.rdd.RDD") >= 0):
                    vcfutils.write_empty_vcf(tx_out_file, samples=[dd.get_sample_name(d) for d in items])
                else:
                    raise
            if tx_out_file.endswith(".vcf"):
                vcfutils.bgzip_and_index(tx_out_file, items[0]["config"])


    # avoid bug in GATK where files can get output as non-compressed
    if out_file.endswith(".gz") and not os.path.exists(out_file + ".tbi"):
        with open(out_file, "r") as in_handle:
            is_plain_text = in_handle.readline().startswith("##fileformat")
        if is_plain_text:
            text_out_file = out_file
            out_file = out_file.replace(".vcf.gz", ".vcf")
            shutil.move(text_out_file, out_file)
    return vcfutils.bgzip_and_index(out_file, items[0]["config"])

def _supports_avx():
    """Check for support for Intel AVX acceleration."""
    if os.path.exists("/proc/cpuinfo"):
        with open("/proc/cpuinfo") as in_handle:
            for line in in_handle:
                if line.startswith("flags") and line.find("avx") > 0:
                    return True
