"""Provide SNP and indel calling using GATK genotyping tools.
"""
import os
import collections
import copy
from distutils.version import LooseVersion

from bcbio import bam, broad, utils
from bcbio.utils import file_exists, safe_makedir
from bcbio.distributed.transaction import file_transaction
from bcbio.distributed.split import grouped_parallel_split_combine
from bcbio.log import logger
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import (process_bam_by_chromosome, subset_variant_regions)
from bcbio.variation.realign import has_aligned_reads
from bcbio.variation import annotation, bamprep, multi, phasing, ploidy, vcfutils, vfilter

# ## GATK Genotype calling

def _shared_gatk_call_prep(align_bams, ref_file, config, dbsnp, region, out_file):
    """Shared preparation work for GATK variant calling.
    """
    broad_runner = broad.runner_from_config(config)
    broad_runner.run_fn("picard_index_ref", ref_file)
    for x in align_bams:
        bam.index(x, config)
    # GATK can only downsample to a minimum of 200
    coverage_depth_max = max(200, utils.get_in(config, ("algorithm", "coverage_depth_max"), 10000))
    coverage_depth_min = utils.get_in(config, ("algorithm", "coverage_depth_min"), 4)
    variant_regions = config["algorithm"].get("variant_regions", None)
    confidence = "4.0" if coverage_depth_min < 4 else "30.0"
    region = subset_variant_regions(variant_regions, region, out_file)

    params = ["-R", ref_file,
              "--standard_min_confidence_threshold_for_calling", confidence,
              "--standard_min_confidence_threshold_for_emitting", confidence,
              "--downsample_to_coverage", str(coverage_depth_max),
              "--downsampling_type", "BY_SAMPLE",
              ]
    for a in annotation.get_gatk_annotations(config):
        params += ["--annotation", a]
    for x in align_bams:
        params += ["-I", x]
    if dbsnp:
        params += ["--dbsnp", dbsnp]
    if region:
        params += ["-L", bamprep.region_to_gatk(region), "--interval_set_rule", "INTERSECTION"]
    return broad_runner, params

def unified_genotyper(align_bams, items, ref_file, assoc_files,
                       region=None, out_file=None):
    """Perform SNP genotyping on the given alignment file.
    """
    if out_file is None:
        out_file = "%s-variants.vcf.gz" % os.path.splitext(align_bams[0])[0]
    if not file_exists(out_file):
        config = items[0]["config"]
        broad_runner, params = \
            _shared_gatk_call_prep(align_bams, ref_file, items[0]["config"], assoc_files["dbsnp"],
                                   region, out_file)
        if (not isinstance(region, (list, tuple)) and
                not all(has_aligned_reads(x, region) for x in align_bams)):
            vcfutils.write_empty_vcf(out_file, config)
        else:
            with file_transaction(out_file) as tx_out_file:
                params += ["-T", "UnifiedGenotyper",
                           "-o", tx_out_file,
                           "-ploidy", (str(ploidy.get_ploidy(items, region))
                                       if broad_runner.gatk_type() == "restricted" else "2"),
                           "--genotype_likelihoods_model", "BOTH"]
                broad_runner.run_gatk(params)
    return out_file

def haplotype_caller(align_bams, items, ref_file, assoc_files,
                       region=None, out_file=None):
    """Call variation with GATK's HaplotypeCaller.

    This requires the full non open-source version of GATK.
    """
    if out_file is None:
        out_file = "%s-variants.vcf.gz" % os.path.splitext(align_bams[0])[0]
    if not file_exists(out_file):
        config = items[0]["config"]
        broad_runner, params = \
            _shared_gatk_call_prep(align_bams, ref_file, items[0]["config"], assoc_files["dbsnp"],
                                   region, out_file)
        assert broad_runner.gatk_type() == "restricted", \
            "Require full version of GATK 2.4+ for haplotype calling"
        if not all(has_aligned_reads(x, region) for x in align_bams):
            vcfutils.write_empty_vcf(out_file, config)
        else:
            with file_transaction(out_file) as tx_out_file:
                params += ["-T", "HaplotypeCaller",
                           "-o", tx_out_file,
                           "--annotation", "ClippingRankSumTest",
                           "--annotation", "DepthPerSampleHC"]
                # Enable hardware based optimizations in GATK 3.1+
                if LooseVersion(broad_runner.gatk_major_version()) >= LooseVersion("3.1"):
                    params += ["--pair_hmm_implementation", "VECTOR_LOGLESS_CACHING"]
                broad_runner.new_resources("gatk-haplotype")
                broad_runner.run_gatk(params)
    return out_file

# ## Variant filtration -- shared functionality

def variant_filtration(call_file, ref_file, vrn_files, data):
    """Filter variant calls using Variant Quality Score Recalibration.

    Newer GATK with Haplotype calling has combined SNP/indel filtering.
    """
    caller = data["config"]["algorithm"].get("variantcaller")
    call_file = ploidy.filter_vcf_by_sex(call_file, data)
    if caller in ["freebayes"]:
        return vfilter.freebayes(call_file, ref_file, vrn_files, data)
    # no additional filtration for callers that filter as part of call process
    elif caller in ["samtools", "varscan", "mutect"]:
        return call_file
    else:
        config = data["config"]
        snp_file, indel_file = vcfutils.split_snps_indels(call_file, ref_file, config)
        snp_filter_file = _variant_filtration_snp(snp_file, ref_file, vrn_files, data)
        indel_filter_file = _variant_filtration_indel(indel_file, ref_file, vrn_files, data)
        orig_files = [snp_filter_file, indel_filter_file]
        out_file = "{base}combined.vcf.gz".format(base=os.path.commonprefix(orig_files))
        return vcfutils.combine_variant_files(orig_files, out_file, ref_file, config)

def _apply_variant_recal(broad_runner, snp_file, ref_file, recal_file,
                         tranch_file, filter_type):
    """Apply recalibration details, returning filtered VCF file.
    """
    base, ext = utils.splitext_plus(snp_file)
    out_file = "{base}-{filter}filter{ext}".format(base=base, ext=ext,
                                                   filter=filter_type)
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            params = ["-T", "ApplyRecalibration",
                      "-R", ref_file,
                      "--input", snp_file,
                      "--out", tx_out_file,
                      "--tranches_file", tranch_file,
                      "--recal_file", recal_file,
                      "--mode", filter_type]
            broad_runner.run_gatk(params)
    return out_file

def _shared_variant_filtration(filter_type, snp_file, ref_file, vrn_files, variantcaller):
    """Share functionality for filtering variants.
    """
    recal_file = "{base}.recal".format(base=utils.splitext_plus(snp_file)[0])
    tranches_file = "{base}.tranches".format(base=utils.splitext_plus(snp_file)[0])
    params = ["-T", "VariantRecalibrator",
              "-R", ref_file,
              "--input", snp_file,
              "--mode", filter_type,
              "-an", "DP",
              "-an", "FS",
              "-an", "ReadPosRankSum",
              "-an", "MQRankSum"]
    if filter_type in ["SNP", "BOTH"]:
        # Haplotype Score no longer calculated for indels as of GATK 2.4
        # and only used for GATK Unified Genotyper calls
        if variantcaller == "gatk":
            params.extend(["-an", "HaplotypeScore"])
        for name, train_info in [("train_hapmap", "known=false,training=true,truth=true,prior=15.0"),
                                 ("train_1000g_omni", "known=false,training=true,truth=false,prior=12.0"),
                                 ("dbsnp", "known=true,training=false,truth=false,prior=8.0")]:
            if name in vrn_files:
                params.extend(["-resource:%s,VCF,%s" % (name.replace("train_", ""), train_info),
                               vrn_files[name]])
    if filter_type in ["INDEL", "BOTH"]:
        params.extend(
            ["-resource:mills,VCF,known=true,training=true,truth=true,prior=12.0",
             vrn_files["train_indels"]])
    return params, recal_file, tranches_file

# ## SNP specific variant filtration

def _variant_filtration_snp(snp_file, ref_file, vrn_files, data):
    """Filter SNP variant calls using GATK best practice recommendations.
    """
    config = data["config"]
    broad_runner = broad.runner_from_config(config)
    filter_type = "SNP"
    variantcaller = config["algorithm"].get("variantcaller", "gatk")
    filters = ["QD < 2.0", "MQ < 40.0", "FS > 60.0",
               "MQRankSum < -12.5", "ReadPosRankSum < -8.0"]
    # GATK Haplotype caller (v2.2) appears to have much larger HaplotypeScores
    # resulting in excessive filtering, so avoid this metric
    if variantcaller not in ["gatk-haplotype"]:
        filters.append("HaplotypeScore > 13.0")
    if not config_utils.use_vqsr([config["algorithm"]]):
        return vfilter.hard_w_expression(snp_file, " || ".join(filters), data, filter_type)
    else:
        # also check if we've failed recal and needed to do strict filtering
        filter_file = "{base}-filter{ext}.vcf.gz".format(base=utils.splitext_plus(snp_file)[0], ext=filter_type)
        if file_exists(filter_file):
            config["algorithm"]["coverage_interval"] = "regional"
            return _variant_filtration_snp(snp_file, ref_file, vrn_files, data)
        assert "train_hapmap" in vrn_files and "train_1000g_omni" in vrn_files, \
            "Need HapMap and 1000 genomes training files"
        params, recal_file, tranches_file = _shared_variant_filtration(
            filter_type, snp_file, ref_file, vrn_files, variantcaller)
        if not file_exists(recal_file):
            with file_transaction(recal_file, tranches_file) as (tx_recal, tx_tranches):
                params.extend(["--recal_file", tx_recal,
                               "--tranches_file", tx_tranches])
                try:
                    broad_runner.new_resources("gatk-vqsr")
                    broad_runner.run_gatk(params, log_error=False)
                # Can fail to run if not enough values are present to train. Rerun with regional
                # filtration approach instead
                except:
                    logger.info("VQSR failed due to lack of training data. Using hard filtering.")
                    config["algorithm"]["coverage_interval"] = "regional"
                    return _variant_filtration_snp(snp_file, ref_file, vrn_files, data)
        return _apply_variant_recal(broad_runner, snp_file, ref_file, recal_file,
                                    tranches_file, filter_type)

# ## Indel specific variant filtration

def _variant_filtration_indel(snp_file, ref_file, vrn_files, data):
    """Filter indel variant calls using GATK best practice recommendations.
    """
    config = data["config"]
    broad_runner = broad.runner_from_config(config)
    filter_type = "INDEL"
    variantcaller = config["algorithm"].get("variantcaller", "gatk")
    if not config_utils.use_vqsr([config["algorithm"]]):
        filterexp = " || ".join(["QD < 2.0", "ReadPosRankSum < -20.0", "FS > 200.0"])
        return vfilter.hard_w_expression(snp_file, filterexp, data, filter_type)
    else:
        # also check if we've failed recal and needed to do strict filtering
        filter_file = "{base}-filter{ext}.vcf.gz".format(base=utils.splitext_plus(snp_file)[0], ext=filter_type)
        if file_exists(filter_file):
            config["algorithm"]["coverage_interval"] = "regional"
            return _variant_filtration_indel(snp_file, ref_file, vrn_files, data)
        assert "train_indels" in vrn_files, "Need indel training file specified"
        params, recal_file, tranches_file = _shared_variant_filtration(
            filter_type, snp_file, ref_file, vrn_files, variantcaller)
        if not file_exists(recal_file):
            with file_transaction(recal_file, tranches_file) as (tx_recal, tx_tranches):
                params.extend(["--recal_file", tx_recal,
                               "--tranches_file", tx_tranches])
                if LooseVersion(broad_runner.gatk_major_version()) >= LooseVersion("2.7"):
                    params.extend(["--numBadVariants", "3000"])
                try:
                    broad_runner.new_resources("gatk-vqsr")
                    broad_runner.run_gatk(params, log_error=False)
                except:
                    logger.info("VQSR failed due to lack of training data. Using hard filtering.")
                    config["algorithm"]["coverage_interval"] = "regional"
                    return _variant_filtration_indel(snp_file, ref_file, vrn_files, data)
        return _apply_variant_recal(broad_runner, snp_file, ref_file, recal_file,
                                    tranches_file, filter_type)

# ## High level functionality to run genotyping in parallel

def get_variantcaller(data):
    return data["config"]["algorithm"].get("variantcaller", "gatk")

def combine_multiple_callers(data):
    """Collapse together variant calls from multiple approaches into variants
    """
    by_bam = collections.defaultdict(list)
    for x in data:
        by_bam[x[0]["work_bam"]].append(x[0])
    out = []
    for grouped_calls in by_bam.itervalues():
        ready_calls = [{"variantcaller": get_variantcaller(x),
                        "vrn_file": x.get("vrn_file"),
                        "validate": x.get("validate")}
                       for x in grouped_calls]
        final = grouped_calls[0]
        def orig_variantcaller_order(x):
            return final["config"]["algorithm"]["orig_variantcaller"].index(x["variantcaller"])
        if len(ready_calls) > 1 and "orig_variantcaller" in final["config"]["algorithm"]:
            final["variants"] = sorted(ready_calls, key=orig_variantcaller_order)
        else:
            final["variants"] = ready_calls
        out.append([final])
    return out

def handle_multiple_variantcallers(data):
    """Split samples that potentially require multiple variant calling approaches.
    """
    assert len(data) == 1
    callers = get_variantcaller(data[0])
    if isinstance(callers, basestring):
        return [data]
    elif not callers:
        return []
    else:
        out = []
        for caller in callers:
            base = copy.deepcopy(data[0])
            base["config"]["algorithm"]["orig_variantcaller"] = \
              base["config"]["algorithm"]["variantcaller"]
            base["config"]["algorithm"]["variantcaller"] = caller
            out.append([base])
        return out

def parallel_variantcall(sample_info, parallel_fn):
    """Provide sample genotyping, running in parallel over individual chromosomes.
    """
    to_process = []
    finished = []
    for x in sample_info:
        if get_variantcaller(x[0]):
            to_process.extend(handle_multiple_variantcallers(x))
        else:
            finished.append(x)
    if len(to_process) > 0:
        split_fn = process_bam_by_chromosome("-variants.vcf.gz", "work_bam",
                                             dir_ext_fn=get_variantcaller)
        processed = grouped_parallel_split_combine(
            to_process, split_fn, multi.group_batches, parallel_fn,
            "variantcall_sample", "split_variants_by_sample", "combine_variant_files",
            "vrn_file", ["sam_ref", "config"])
        finished.extend(processed)
    return finished

def get_variantcallers():
    from bcbio.variation import freebayes, cortex, samtools, varscan, mutect
    return {"gatk": unified_genotyper,
            "gatk-haplotype": haplotype_caller,
            "freebayes": freebayes.run_freebayes,
            "cortex": cortex.run_cortex,
            "samtools": samtools.run_samtools,
            "varscan": varscan.run_varscan,
            "mutect": mutect.mutect_caller}

def variantcall_sample(data, region=None, out_file=None):
    """Parallel entry point for doing genotyping of a region of a sample.
    """
    if out_file is None or not os.path.exists(out_file) or not os.path.lexists(out_file):
        safe_makedir(os.path.dirname(out_file))
        sam_ref = data["sam_ref"]
        config = data["config"]
        caller_fns = get_variantcallers()
        caller_fn = caller_fns[config["algorithm"].get("variantcaller", "gatk")]
        if isinstance(data["work_bam"], basestring):
            align_bams = [data["work_bam"]]
            items = [data]
        else:
            align_bams = data["work_bam"]
            items = data["work_items"]
        call_file = "%s-raw%s" % utils.splitext_plus(out_file)
        call_file = caller_fn(align_bams, items, sam_ref,
                              data["genome_resources"]["variation"],
                              region, call_file)
        if data["config"]["algorithm"].get("phasing", False) == "gatk":
            call_file = phasing.read_backed_phasing(call_file, align_bams, sam_ref, region, config)
        utils.symlink_plus(call_file, out_file)
        if "work_items" in data:
            del data["work_items"]
    data["vrn_file"] = out_file
    return [data]
