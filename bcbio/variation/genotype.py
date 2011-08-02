"""Provide SNP, indel calling and variation analysis using GATK genotyping tools.

Genotyping:

http://www.broadinstitute.org/gsa/wiki/index.php/Unified_genotyper
http://www.broadinstitute.org/gsa/wiki/index.php/Local_realignment_around_indels
http://www.broadinstitute.org/gsa/wiki/index.php/IndelGenotyper
http://www.broadinstitute.org/gsa/wiki/index.php/VariantFiltrationWalker

Variant Evaluation:

http://www.broadinstitute.org/gsa/wiki/index.php/VariantEval
"""
import os
import itertools

from bcbio import broad
from bcbio.utils import file_transaction

# ## SNP Genotyping

def gatk_genotyper(align_bam, ref_file, config, dbsnp=None):
    """Perform genotyping and filtration on a sorted aligned BAM file.
    """
    picard = broad.runner_from_config(config)
    picard.run_fn("picard_index_ref", ref_file)
    picard.run_fn("picard_index", align_bam)
    snp_file = _unified_genotyper(picard, align_bam, ref_file, dbsnp)
    filter_snp = _variant_filtration(picard, snp_file, ref_file)
    return filter_snp

def _unified_genotyper(picard, align_bam, ref_file, dbsnp=None):
    """Perform SNP genotyping on the given alignment file.
    """
    out_file = "%s-snp.vcf" % os.path.splitext(align_bam)[0]
    params = ["-T", "UnifiedGenotyper",
              "-I", align_bam,
              "-R", ref_file,
              "-o", out_file,
              "-A", "DepthOfCoverage",
              "-A", "AlleleBalance",
              "-A", "HomopolymerRun",
              "-A", "QualByDepth",
              "--genotype_likelihoods_model", "SNP",
              "-baq", "CALCULATE_AS_NECESSARY",
              "--standard_min_confidence_threshold_for_calling", "10.0",
              "--standard_min_confidence_threshold_for_emitting", "10.0",
              #"--trigger_min_confidence_threshold_for_calling", "10.0",
              #"--trigger_min_confidence_threshold_for_emitting", "10.0",
              "--downsample_to_coverage", 10000,
              "--min_base_quality_score", 20,
              "-l", "INFO",
              ]
    if dbsnp:
        params += ["-B:dbsnp,VCF", dbsnp]
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        with file_transaction(out_file):
            picard.run_gatk(params)
    return out_file

def _variant_filtration(picard, snp_file, ref_file):
    """Filter out problematic SNP calls.

    Recommended Broad hard filtering for deep coverage exomes:
        QUAL < 30.0 || AB > 0.75 && DP > 40 || QD < 5.0 || HRun > 5 || SB > -0.10
    """
    out_file = "%s-filter%s" % os.path.splitext(snp_file)
    params = ["-T", "VariantFiltration",
              "-R", ref_file,
              "-o", out_file,
              "-B:variant,VCF", snp_file,
              "--filterName", "QUALFilter",
              "--filterExpression", "QUAL <= 50.0",
              "--filterName", "QDFilter",
              "--filterExpression", "QD < 5.0",
              "--filterName", "ABFilter",
              "--filterExpression", "AB > 0.75 && DP > 40",
              "--filterName", "HRunFilter",
              "--filterExpression", "HRun > 3.0",
              "--filterName", "SBFilter",
              "--filterExpression", "SB > -0.10",
              "-l", "INFO",
              ]
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        with file_transaction(out_file):
            picard.run_gatk(params)
    return out_file

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
    params = ["-T", "VariantEval",
              "-R", ref_file,
              "-B:eval,VCF", vcf_in,
              "-B:dbsnp,VCF", dbsnp,
              "-ST", "Filter",
              "-o", out_file,
              "-l", "INFO"
              ]
    if target_intervals:
        params.extend(["-L", target_intervals])
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        with file_transaction(out_file):
            picard.run_gatk(params)
    return out_file
