"""GATK variant calling -- HaplotypeCaller and UnifiedGenotyper.
"""
from distutils.version import LooseVersion

import toolz as tz

from bcbio import bam, broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.pipeline import datadict as dd
from bcbio.variation import annotation, bamprep, ploidy

def _shared_gatk_call_prep(align_bams, items, ref_file, dbsnp, region, out_file):
    """Shared preparation work for GATK variant calling.
    """
    data = items[0]
    config = data["config"]
    broad_runner = broad.runner_from_config(config)
    broad_runner.run_fn("picard_index_ref", ref_file)
    for x in align_bams:
        bam.index(x, config)
    params = ["-R", ref_file]
    if dd.is_set_coverage_depth_max(data):
        coverage_depth_max = dd.get_coverage_depth_max(data)
        # GATK can only downsample to a minimum of 200
        coverage_depth_max = max([200, coverage_depth_max])
        params += ["--downsample_to_coverage", str(coverage_depth_max),
                   "--downsampling_type", "BY_SAMPLE"]
    coverage_depth_min = tz.get_in(["algorithm", "coverage_depth_min"], config)
    if coverage_depth_min and coverage_depth_min < 4:
        confidence = "4.0"
        params += ["--standard_min_confidence_threshold_for_calling", confidence,
                   "--standard_min_confidence_threshold_for_emitting", confidence]
    for a in annotation.get_gatk_annotations(config):
        params += ["--annotation", a]
    for x in align_bams:
        params += ["-I", x]
    if dbsnp:
        params += ["--dbsnp", dbsnp]
    variant_regions = tz.get_in(["algorithm", "variant_regions"], config)
    region = subset_variant_regions(variant_regions, region, out_file, items)
    if region:
        params += ["-L", bamprep.region_to_gatk(region), "--interval_set_rule", "INTERSECTION"]
    return broad_runner, params

def unified_genotyper(align_bams, items, ref_file, assoc_files,
                       region=None, out_file=None):
    """Perform SNP genotyping on the given alignment file.
    """
    if out_file is None:
        out_file = "%s-variants.vcf.gz" % utils.splitext_plus(align_bams[0])[0]
    if not utils.file_exists(out_file):
        broad_runner, params = \
            _shared_gatk_call_prep(align_bams, items,
                                   ref_file, assoc_files.get("dbsnp"),
                                   region, out_file)
        with file_transaction(items[0], out_file) as tx_out_file:
            params += ["-T", "UnifiedGenotyper",
                       "-o", tx_out_file,
                       "-ploidy", (str(ploidy.get_ploidy(items, region))
                                   if broad_runner.gatk_type() == "restricted" else "2"),
                       "--genotype_likelihoods_model", "BOTH"]
            broad_runner.run_gatk(params)
    return out_file

def _joint_calling(items):
    """Determine if this call feeds downstream into joint calls.
    """
    jointcaller = tz.get_in(("config", "algorithm", "jointcaller"), items[0])
    if jointcaller:
        assert len(items) == 1, "Can only do joint calling preparation with GATK with single samples"
        assert tz.get_in(("metadata", "batch"), items[0]) is not None, \
            "Joint calling requires batched samples, %s has no metadata batch." % dd.get_sample_name(items[0])
    return jointcaller

def haplotype_caller(align_bams, items, ref_file, assoc_files,
                       region=None, out_file=None):
    """Call variation with GATK's HaplotypeCaller.

    This requires the full non open-source version of GATK.
    """
    if out_file is None:
        out_file = "%s-variants.vcf.gz" % utils.splitext_plus(align_bams[0])[0]
    if not utils.file_exists(out_file):
        broad_runner, params = \
            _shared_gatk_call_prep(align_bams, items,
                                   ref_file, assoc_files.get("dbsnp"),
                                   region, out_file)
        assert broad_runner.gatk_type() == "restricted", \
            "Require full version of GATK 2.4+ for haplotype calling"
        with file_transaction(items[0], out_file) as tx_out_file:
            params += ["-T", "HaplotypeCaller",
                       "-o", tx_out_file,
                       "--annotation", "ClippingRankSumTest",
                       "--annotation", "DepthPerSampleHC"]
            # Enable hardware based optimizations in GATK 3.1+
            if LooseVersion(broad_runner.gatk_major_version()) >= LooseVersion("3.1"):
                params += ["--pair_hmm_implementation", "VECTOR_LOGLESS_CACHING"]
            # Enable non-diploid calling in GATK 3.3+
            if LooseVersion(broad_runner.gatk_major_version()) >= LooseVersion("3.3"):
                params += ["-ploidy", str(ploidy.get_ploidy(items, region))]
            if _joint_calling(items):  # Prepare gVCFs if doing joint calling
                params += ["--emitRefConfidence", "GVCF", "--variant_index_type", "LINEAR",
                           "--variant_index_parameter", "128000"]
            resources = config_utils.get_resources("gatk-haplotype", items[0]["config"])
            if "options" in resources:
                params += [str(x) for x in resources.get("options", [])]
            broad_runner.new_resources("gatk-haplotype")
            broad_runner.run_gatk(params)
    return out_file
