"""Provide SNP and indel calling using GATK genotyping tools.

http://www.broadinstitute.org/gsa/wiki/index.php/Unified_genotyper
http://www.broadinstitute.org/gsa/wiki/index.php/Local_realignment_around_indels
http://www.broadinstitute.org/gsa/wiki/index.php/IndelGenotyper
http://www.broadinstitute.org/gsa/wiki/index.php/VariantFiltrationWalker
"""
import os

from bcbio.broad import BroadRunner

def gatk_genotyper(align_bam, ref_file, config, dbsnp=None):
    """Perform genotyping and filtration on a sorted aligned BAM file.
    """
    picard = BroadRunner(config["program"]["picard"],
                         config["program"].get("gatk", ""),
                         max_memory=config["algorithm"].get("java_memory", ""))
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
        picard.run_gatk(params)
    return out_file
