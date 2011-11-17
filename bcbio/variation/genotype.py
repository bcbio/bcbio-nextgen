"""Provide SNP, indel calling and variation analysis using GATK genotyping tools.

Genotyping:

http://www.broadinstitute.org/gsa/wiki/index.php/Best_Practice_Variant_Detection_with_the_GATK_v3
http://www.broadinstitute.org/gsa/wiki/index.php/Unified_genotyper

Variant Evaluation:

http://www.broadinstitute.org/gsa/wiki/index.php/VariantEval
"""
import os
import itertools

from bcbio import broad
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.distributed.split import parallel_split_combine
from bcbio.pipeline.shared import (split_bam_by_chromosome, configured_ref_file)
from bcbio.variation.realign import has_aligned_reads

# ## GATK Genotype calling

def unified_genotyper(align_bam, ref_file, config, dbsnp=None,
                       region=None, out_file=None):
    """Perform SNP genotyping on the given alignment file.
    """
    broad_runner = broad.runner_from_config(config)
    broad_runner.run_fn("picard_index_ref", ref_file)
    broad_runner.run_fn("picard_index", align_bam)
    coverage_depth = config["algorithm"].get("coverage_depth", "high").lower()
    if coverage_depth in ["low"]:
        confidence = "4.0"
    else:
        confidence = "30.0"
    if out_file is None:
        out_file = "%s-variants.vcf" % os.path.splitext(align_bam)[0]
    if not file_exists(out_file):
        if has_aligned_reads(align_bam, region):
            with file_transaction(out_file) as tx_out_file:
                params = ["-T", "UnifiedGenotyper",
                          "-I", align_bam,
                          "-R", ref_file,
                          "-o", tx_out_file,
                          "--annotation", "QualByDepth",
                          "--annotation", "HaplotypeScore",
                          "--annotation", "MappingQualityRankSumTest",
                          "--annotation", "ReadPosRankSumTest",
                          "--annotation", "FisherStrand",
                          "--annotation", "RMSMappingQuality",
                          "--annotation", "DepthOfCoverage",
                          "--genotype_likelihoods_model", "BOTH",
                          "--standard_min_confidence_threshold_for_calling", confidence,
                          "--standard_min_confidence_threshold_for_emitting", confidence,
                          "-l", "INFO",
                          ]
                if dbsnp:
                    params += ["--dbsnp", dbsnp]
                if region:
                    params += ["-L", region]
                broad_runner.run_gatk(params)
        else:
            with open(out_file, "w") as out_handle:
                out_handle.write("##fileformat=VCFv4.1\n"
                                 "## No variants; no reads aligned in region\n"
                                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    return out_file

# ## Utility functions for dealing with VCF files

def split_snps_indels(broad_runner, orig_file, ref_file):
    """Split a variant call file into SNPs and INDELs for processing.
    """
    base, ext = os.path.splitext(orig_file)
    snp_file = "{base}-snp{ext}".format(base=base, ext=ext)
    indel_file = "{base}-indel{ext}".format(base=base, ext=ext)
    params = ["-T", "SelectVariants",
              "-R", ref_file,
              "--variant", orig_file]
    for out_file, select_type in [(snp_file, "SNP"), (indel_file, "INDEL")]:
        if not file_exists(out_file):
            with file_transaction(out_file) as tx_out_file:
                cur_params = params + ["--out", tx_out_file,
                                       "--selectTypeToInclude", select_type]
                broad_runner.run_gatk(cur_params)
    return snp_file, indel_file

def combine_variant_files(orig_files, out_file, ref_file, config):
    """Combine multiple VCF files into a single output file.
    """
    broad_runner = broad.runner_from_config(config)
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            params = ["-T", "CombineVariants",
                      "-R", ref_file,
                      "--out", tx_out_file]
            priority_order = []
            for orig_file in orig_files:
                name = os.path.splitext(os.path.basename(orig_file))[0]
                params.extend(["--variant:{name}".format(name=name), orig_file])
                priority_order.append(name)
            params.extend(["--rod_priority_list", ",".join(priority_order)])
            broad_runner.run_gatk(params)
    return out_file

# ## Variant filtration -- shared functionality

def variant_filtration(call_file, ref_file, vrn_files, config):
    """Filter variant calls using Variant Quality Score Recalibration.
    """
    broad_runner = broad.runner_from_config(config)
    snp_file, indel_file = split_snps_indels(broad_runner, call_file, ref_file)
    snp_filter_file = _variant_filtration_snp(broad_runner, snp_file, ref_file,
                                              vrn_files, config)
    indel_filter_file = _variant_filtration_indel(broad_runner, indel_file,
                                                  ref_file, vrn_files, config)
    orig_files = [snp_filter_file, indel_filter_file]
    out_file = "{base}combined.vcf".format(base=os.path.commonprefix(orig_files))
    return combine_variant_files(orig_files, out_file, ref_file, config)

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

def _shared_variant_filtration(filter_type, snp_file, ref_file):
    """Share functionality for filtering variants.
    """
    recal_file = "{base}.recal".format(base = os.path.splitext(snp_file)[0])
    tranches_file = "{base}.tranches".format(base = os.path.splitext(snp_file)[0])
    params = ["-T", "VariantRecalibrator",
              "-R", ref_file,
              "--input", snp_file,
              "--mode", filter_type]
    return params, recal_file, tranches_file

def _variant_filtration_no_recal(broad_runner, snp_file, ref_file, filter_type,
                                 expressions):
    """Perform hard filtering if coverage is in limited regions.

    Variant quality score recalibration will not work on some regions; it
    requires enough positions to train the model.
    """
    base, ext = os.path.splitext(snp_file)
    out_file = "{base}-filter{ftype}{ext}".format(base=base, ext=ext,
                                                  ftype=filter_type)
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            params = ["-T", "VariantFiltration",
                      "-R", ref_file,
                      "--out", tx_out_file,
                      "--variant", snp_file]
            for exp in expressions:
                params.extend(["--filterName", "GATKStandard{e}".format(e=exp.split()[0]),
                               "--filterExpression", exp])
            broad_runner.run_gatk(params)
    return out_file

# ## SNP specific variant filtration

def _variant_filtration_snp(broad_runner, snp_file, ref_file, vrn_files,
                            config):
    """Filter SNP variant calls using GATK best practice recommendations.
    """
    filter_type = "SNP"
    cov_interval = config["algorithm"].get("coverage_interval", "exome").lower()
    params, recal_file, tranches_file = _shared_variant_filtration(
        filter_type, snp_file, ref_file)
    assert vrn_files.train_hapmap and vrn_files.train_1000g_omni, \
           "Need HapMap and 1000 genomes training files"
    if cov_interval == "regional":
        return _variant_filtration_no_recal(broad_runner, snp_file, ref_file, filter_type,
                                            ["QD < 5.0", "HRun > 5", "FS > 200.0"])
    else:
        params.extend(
            ["-resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0",
             vrn_files.train_hapmap,
             "-resource:omni,VCF,known=false,training=true,truth=false,prior=12.0",
             vrn_files.train_1000g_omni,
             "-resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0",
             vrn_files.dbsnp,
              "-an", "QD",
              "-an", "HaplotypeScore",
              "-an", "MQRankSum",
              "-an", "ReadPosRankSum",
              "-an", "FS",
              "-an", "MQ"])
        if cov_interval == "exome":
            params.extend(["--maxGaussians", "4", "--percentBadVariants", "0.05"])
        else:
            params.extend(["-an", "DP"])
        if not file_exists(recal_file):
            with file_transaction(recal_file, tranches_file) as (tx_recal, tx_tranches):
                params.extend(["--recal_file", tx_recal,
                               "--tranches_file", tx_tranches])
                try:
                    broad_runner.run_gatk(params)
                # Can fail to run if not enough values are present to train. Rerun with regional
                # filtration approach instead
                except:
                    config["algorithm"]["coverage_interval"] = "regional"
                    return _variant_filtration_snp(broad_runner, snp_file, ref_file, vrn_files,
                                                   config)
        return _apply_variant_recal(broad_runner, snp_file, ref_file, recal_file,
                                    tranches_file, filter_type)

# ## Indel specific variant filtration

def _variant_filtration_indel(broad_runner, snp_file, ref_file, vrn_files,
                              config):
    """Filter indel variant calls using GATK best practice recommendations.
    """
    filter_type = "INDEL"
    cov_interval = config["algorithm"].get("coverage_interval", "exome").lower()
    params, recal_file, tranches_file = _shared_variant_filtration(
        filter_type, snp_file, ref_file)
    if cov_interval in ["exome", "regional"]:
        return _variant_filtration_no_recal(broad_runner, snp_file, ref_file, filter_type,
                                            ["QD < 2.0", "ReadPosRankSum < -20.0", "FS > 200.0"])
    else:
        assert vrn_files.train_indels, \
               "Need indel training file specified"
        params.extend(
            ["-resource:mills,VCF,known=true,training=true,truth=true,prior=12.0",
             vrn_files.train_indels,
             "-an", "QD",
             "-an", "FS",
             "-an", "HaplotypeScore",
             "-an", "ReadPosRankSum"])
        if not file_exists(recal_file):
            with file_transaction(recal_file, tranches_file) as (tx_recal, tx_tranches):
                params.extend(["--recal_file", tx_recal,
                               "--tranches_file", tx_tranches])
                broad_runner.run_gatk(params)
        return _apply_variant_recal(broad_runner, snp_file, ref_file, recal_file,
                                    tranches_file, filter_type)

# ## Variant evaluation

def gatk_evaluate_variants(vcf_file, ref_file, config, dbsnp=None, intervals=None):
    """Evaluate variants, return SNP counts and Transition/Transversion ratios.
    """
    runner = broad.runner_from_config(config)
    eval_file = variant_eval(vcf_file, ref_file, dbsnp, intervals, runner)
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
    return dict(total=total, dbsnp_pct = dbsnp, titv_all=titv_all,
                titv_dbsnp=titv_dbsnp, titv_novel=titv_novel)

def _extract_eval_stats(eval_file):
    """Parse statistics of interest from GATK output file.
    """
    stats = dict()
    for snp_type in ['called', 'filtered']:
        stats[snp_type]  = dict()
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
    with open(in_file) as in_handle:
        # read until we reach the analysis
        for line in in_handle:
            if (line.startswith("##:GATKReport") and
                line.find(analysis_name) > 0):
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

def variant_eval(vcf_in, ref_file, dbsnp, target_intervals, picard):
    """Evaluate variants in comparison with dbSNP reference.
    """
    out_file = "%s.eval" % os.path.splitext(vcf_in)[0]
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            params = ["-T", "VariantEval",
                      "-R", ref_file,
                      "--eval", vcf_in,
                      "--dbsnp", dbsnp,
                      "-ST", "Filter",
                      "-o", tx_out_file,
                      "-l", "INFO"
                      ]
            if target_intervals:
                # BED file target intervals are explicit with GATK 1.3
                # http://getsatisfaction.com/gsa/topics/
                # gatk_v1_3_and_bed_interval_file_must_be_parsed_through_tribble
                if _is_bed_file(target_intervals):
                    flag = "-L:bed"
                else:
                    flag = "-L"
                params.extend([flag, target_intervals])
            picard.run_gatk(params)
    return out_file

def _is_bed_file(fname):
    """Simple check if a file is in BED format.
    """
    if fname.lower().endswith(".bed"):
        return True
    with open(fname) as in_handle:
        for line in in_handle:
            if not line.startswith("#"):
                parts = line.split("\t")
                if len(parts) > 3:
                    try:
                        int(parts[1])
                        int(parts[2])
                        return True
                    except ValueError:
                        pass
                break
    return False

# ## High level functionality to run genotyping in parallel

def parallel_unified_genotyper(sample_info, parallel_fn):
    """Realign samples, running in parallel over individual chromosomes.
    """
    to_process = []
    finished = []
    for x in sample_info:
        if x[0]["config"]["algorithm"]["snpcall"]:
            to_process.append(x)
        else:
            finished.append(x)
    if len(to_process) > 0:
        split_fn = split_bam_by_chromosome("-variants.vcf", "work_bam")
        processed = parallel_split_combine(to_process, split_fn, parallel_fn,
                                           "unified_genotyper_sample",
                                           "combine_variant_files",
                                           "vrn_file", ["sam_ref", "config"])
        finished.extend(processed)
    return finished

def unified_genotyper_sample(data, region=None, out_file=None):
    """Parallel entry point for doing genotyping of a region of a sample.
    """
    if data["config"]["algorithm"]["snpcall"]:
        sam_ref = data["sam_ref"]
        config = data["config"]
        data["vrn_file"] = unified_genotyper(data["work_bam"], sam_ref, config,
                                             configured_ref_file("dbsnp", config, sam_ref),
                                             region, out_file)
    return [data]
