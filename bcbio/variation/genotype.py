"""Provide SNP, indel calling and variation analysis using GATK genotyping tools.

Genotyping:

http://www.broadinstitute.org/gsa/wiki/index.php/Best_Practice_Variant_Detection_with_the_GATK_v3
http://www.broadinstitute.org/gsa/wiki/index.php/Unified_genotyper

Variant Evaluation:

http://www.broadinstitute.org/gsa/wiki/index.php/VariantEval
"""
import os
import collections
import copy
from distutils.version import LooseVersion
import itertools

from bcbio import bam, broad
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
    coverage_depth = config["algorithm"].get("coverage_depth", "high").lower()
    variant_regions = config["algorithm"].get("variant_regions", None)
    confidence = "4.0" if coverage_depth in ["low"] else "30.0"
    region = subset_variant_regions(variant_regions, region, out_file)

    params = ["-R", ref_file,
              "--standard_min_confidence_threshold_for_calling", confidence,
              "--standard_min_confidence_threshold_for_emitting", confidence,
              "--downsample_to_coverage", "250",
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
        out_file = "%s-variants.vcf" % os.path.splitext(align_bams[0])[0]
    if not file_exists(out_file):
        broad_runner, params = \
            _shared_gatk_call_prep(align_bams, ref_file, items[0]["config"], assoc_files["dbsnp"],
                                   region, out_file)
        if (not isinstance(region, (list, tuple)) and
                not all(has_aligned_reads(x, region) for x in align_bams)):
            vcfutils.write_empty_vcf(out_file)
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
        out_file = "%s-variants.vcf" % os.path.splitext(align_bams[0])[0]
    if not file_exists(out_file):
        broad_runner, params = \
            _shared_gatk_call_prep(align_bams, ref_file, items[0]["config"], assoc_files["dbsnp"],
                                   region, out_file)
        assert broad_runner.gatk_type() == "restricted", \
            "Require full version of GATK 2.4+ for haplotype calling"
        if not all(has_aligned_reads(x, region) for x in align_bams):
            vcfutils.write_empty_vcf(out_file)
        else:
            with file_transaction(out_file) as tx_out_file:
                params += ["-T", "HaplotypeCaller",
                           "-o", tx_out_file]
                #params = _gatk_location_hack(params)
                broad_runner.new_resources("gatk-haplotype")
                broad_runner.run_gatk(params)
    return out_file

def _gatk_location_hack(args):
    """Temporary work around for issues in GATK 2.4-9 and 2.5-2 HaplotypeCaller.

    Fixes:
    - softclipped reads at end of chromosomes.
      Pads these regions to avoid working exclusively with the end. Needs to be
      fixed properly in upstream GATK.
    - Problematic assembly around repeat regions. Excludes these regions.
    """
    region_idxs = [i + 1 for i, x in enumerate(args) if x == "-L"]
    # padding
    problem_chrs = ["GL000195.1"]
    pad_start = 250
    # exclusion
    exclude_args = {"20": ["-XL", "20:33972777-33973070"]}
    extra_args = []
    for ridx in region_idxs:
        if os.path.isfile(args[ridx]):
            with open(args[ridx]) as in_handle:
                chrom = in_handle.readline().split()[0]
        else:
            chrom = args[ridx].split(":")[0]
        if chrom in problem_chrs and args[ridx].find(":") > 0:
            chrom, rest = args[ridx].split(":")
            start, end = rest.split("-")
            new_start = max([pad_start, int(start)])
            args[ridx] = "%s:%s-%s" % (chrom, new_start, end)
        elif chrom in exclude_args:
            extra_args.extend(exclude_args[chrom])
    return args + extra_args

# ## Variant filtration -- shared functionality

def variant_filtration(call_file, ref_file, vrn_files, data):
    """Filter variant calls using Variant Quality Score Recalibration.

    Newer GATK with Haplotype calling has combined SNP/indel filtering.
    """
    caller = data["config"]["algorithm"].get("variantcaller")
    call_file = ploidy.filter_vcf_by_sex(call_file, data)
    if caller in ["freebayes"]:
        return vfilter.freebayes(call_file, ref_file, vrn_files, data["config"])
    # no additional filtration for callers that filter as part of call process
    elif caller in ["samtools", "varscan"]:
        return call_file
    else:
        config = data["config"]
        snp_file, indel_file = vcfutils.split_snps_indels(call_file, ref_file, config)
        snp_filter_file = _variant_filtration_snp(snp_file, ref_file, vrn_files, config)
        indel_filter_file = _variant_filtration_indel(indel_file, ref_file, vrn_files, config)
        orig_files = [snp_filter_file, indel_filter_file]
        out_file = "{base}combined.vcf".format(base=os.path.commonprefix(orig_files))
        return vcfutils.combine_variant_files(orig_files, out_file, ref_file, config)

def _apply_variant_recal(broad_runner, snp_file, ref_file, recal_file,
                         tranch_file, filter_type):
    """Apply recalibration details, returning filtered VCF file.
    """
    base, ext = os.path.splitext(snp_file)
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
    recal_file = "{base}.recal".format(base=os.path.splitext(snp_file)[0])
    tranches_file = "{base}.tranches".format(base=os.path.splitext(snp_file)[0])
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

def _variant_filtration_snp(snp_file, ref_file, vrn_files, config):
    """Filter SNP variant calls using GATK best practice recommendations.
    """
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
        return vfilter.jexl_hard(broad_runner, snp_file, ref_file, filter_type, filters)
    else:
        # also check if we've failed recal and needed to do strict filtering
        filter_file = "{base}-filter{ext}.vcf".format(base=os.path.splitext(snp_file)[0], ext=filter_type)
        if file_exists(filter_file):
            config["algorithm"]["coverage_interval"] = "regional"
            return _variant_filtration_snp(snp_file, ref_file, vrn_files, config)
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
                    return _variant_filtration_snp(snp_file, ref_file, vrn_files, config)
        return _apply_variant_recal(broad_runner, snp_file, ref_file, recal_file,
                                    tranches_file, filter_type)

# ## Indel specific variant filtration

def _variant_filtration_indel(snp_file, ref_file, vrn_files, config):
    """Filter indel variant calls using GATK best practice recommendations.
    """
    broad_runner = broad.runner_from_config(config)
    filter_type = "INDEL"
    variantcaller = config["algorithm"].get("variantcaller", "gatk")
    if not config_utils.use_vqsr([config["algorithm"]]):
        return vfilter.jexl_hard(broad_runner, snp_file, ref_file, filter_type,
                                 ["QD < 2.0", "ReadPosRankSum < -20.0", "FS > 200.0"])
    else:
        # also check if we've failed recal and needed to do strict filtering
        filter_file = "{base}-filter{ext}.vcf".format(base=os.path.splitext(snp_file)[0], ext=filter_type)
        if file_exists(filter_file):
            config["algorithm"]["coverage_interval"] = "regional"
            return _variant_filtration_indel(snp_file, ref_file, vrn_files, config)
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
                    return _variant_filtration_indel(snp_file, ref_file, vrn_files, config)
        return _apply_variant_recal(broad_runner, snp_file, ref_file, recal_file,
                                    tranches_file, filter_type)

# ## Variant evaluation

def gatk_evaluate_variants(vcf_file, ref_file, config, dbsnp=None):
    """Evaluate variants, return SNP counts and Transition/Transversion ratios.
    """
    runner = broad.runner_from_config(config)
    eval_file = variant_eval(vcf_file, ref_file, dbsnp, runner)
    stats = _extract_eval_stats(eval_file)
    return _format_stats(stats['called'])

def _format_stats(stats):
    """Convert statistics into high level summary of major variables.
    """
    total = sum(itertools.chain.from_iterable(s.itervalues() for s in stats.itervalues()))
    if total > 0:
        dbsnp = sum(stats['known'].itervalues()) / float(total) * 100.0
    else:
        dbsnp = -1.0
    tv_dbsnp = stats['known']['tv']
    ti_dbsnp = stats['known']['ti']
    tv_novel = stats['novel']['tv']
    ti_novel = stats['novel']['ti']
    if tv_novel > 0 and tv_dbsnp > 0:
        titv_all = float(ti_novel + ti_dbsnp) / float(tv_novel + tv_dbsnp)
        titv_dbsnp = float(ti_dbsnp) / float(tv_dbsnp)
        titv_novel = float(ti_novel) / float(tv_novel)
    else:
        titv_all, titv_dbsnp, titv_novel = (-1.0, -1.0, -1.0)
    return dict(total=total, dbsnp_pct=dbsnp, titv_all=titv_all,
                titv_dbsnp=titv_dbsnp, titv_novel=titv_novel)

def _extract_eval_stats(eval_file):
    """Parse statistics of interest from GATK output file.
    """
    stats = dict()
    for snp_type in ['called', 'filtered']:
        stats[snp_type] = dict()
        for dbsnp_type in ['known', 'novel']:
            stats[snp_type][dbsnp_type] = dict(ti=0, tv=0)
    for line in _eval_analysis_type(eval_file, "Ti/Tv Variant Evaluator"):
        if line[1:3] == ['dbsnp', 'eval']:
            snp_type = line[3]
            dbsnp_type = line[5]
            try:
                cur = stats[snp_type][dbsnp_type]
            except KeyError:
                cur = None
            if cur:
                stats[snp_type][dbsnp_type]["ti"] = int(line[6])
                stats[snp_type][dbsnp_type]["tv"] = int(line[7])
    return stats

def _eval_analysis_type(in_file, analysis_name):
    """Retrieve data lines associated with a particular analysis.
    """
    supported_versions = ["v0.2", "v1.0", "v1.1"]
    with open(in_file) as in_handle:
        # read until we reach the analysis
        for line in in_handle:
            if line.startswith(("##:GATKReport", "#:GATKReport")):
                version = line.split()[0].split(".", 1)[-1].split(":")[0]
                assert version in supported_versions, \
                       "Unexpected GATKReport version: {0}".format(version)
                if line.find(analysis_name) > 0:
                    break
        # read off header lines
        for _ in range(1):
            in_handle.next()
        # read the table until a blank line
        for line in in_handle:
            if not line.strip():
                break
            parts = line.rstrip("\n\r").split()
            yield parts

def variant_eval(vcf_in, ref_file, dbsnp, picard):
    """Evaluate variants in comparison with dbSNP reference.
    """
    out_file = "%s.eval" % os.path.splitext(vcf_in)[0]
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            params = ["-T", "VariantEval",
                      "-R", ref_file,
                      "--eval", vcf_in,
                      "-L", vcf_in,
                      "--dbsnp", dbsnp,
                      "-ST", "Filter",
                      "-o", tx_out_file,
                      "-l", "INFO",
                      "--doNotUseAllStandardModules",
                      "--evalModule", "CompOverlap",
                      "--evalModule", "CountVariants",
                      "--evalModule", "TiTvVariantEvaluator",
                      "--evalModule", "ValidationReport",
                      "--stratificationModule", "Filter"]
            picard.run_gatk(params)
    return out_file

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
        if len(ready_calls) > 1:
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
        split_fn = process_bam_by_chromosome("-variants.vcf", "work_bam",
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
    call_file = "%s-raw%s" % os.path.splitext(out_file)
    call_file = caller_fn(align_bams, items, sam_ref,
                          data["genome_resources"]["variation"],
                          region, call_file)
    if data["config"]["algorithm"].get("phasing", False) == "gatk":
        call_file = phasing.read_backed_phasing(call_file, align_bams, sam_ref, region, config)
    for ext in ["", ".idx"]:
        if not os.path.exists(out_file + ext):
            if os.path.exists(call_file + ext):
                try:
                    os.symlink(call_file + ext, out_file + ext)
                except OSError, msg:
                    if  str(msg).find("File exists") == -1:
                        raise
    if "work_items" in data:
        del data["work_items"]
    data["vrn_file"] = out_file
    return [data]
