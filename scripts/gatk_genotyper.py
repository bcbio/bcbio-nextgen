#!/usr/bin/env python
"""Provide SNP and indel calling using GATK genotyping tools.

http://www.broadinstitute.org/gsa/wiki/index.php/Unified_genotyper
http://www.broadinstitute.org/gsa/wiki/index.php/Local_realignment_around_indels
http://www.broadinstitute.org/gsa/wiki/index.php/IndelGenotyper
http://www.broadinstitute.org/gsa/wiki/index.php/VariantFiltrationWalker


Usage:
    gatk_genotyper.py <config file> <reference file> <align BAM file>
                      <dbsnp file>

Requires:
  - Picard
  - GATK
"""
import os
import sys
import subprocess

import yaml

from bcbio.broad import BroadRunner
from bcbio.utils import curdir_tmpdir

def main(config_file, ref_file, align_bam, dbsnp=None):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    picard = BroadRunner(config["program"]["picard"],
                         config["program"].get("gatk", ""),
                         max_memory=config["algorithm"].get("java_memory", ""))
    ref_dict = picard.run_fn("picard_index_ref", ref_file)
    picard.run_fn("picard_index", align_bam)
    realign_bam = picard.run_fn("gatk_realigner", align_bam, ref_file, dbsnp)
    realign_sort_bam = picard_fixmate(picard, realign_bam)
    picard.run_fn("picard_index", realign_sort_bam)
    snp_file = unified_genotyper(picard, realign_sort_bam, ref_file,
                                 dbsnp)
    filter_snp = variant_filtration(picard, snp_file, ref_file)

def unified_genotyper(picard, align_bam, ref_file, dbsnp=None):
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

def variant_filtration(picard, snp_file, ref_file):
    """Filter out problematic SNP calls.

    XXX missing:
        interval list

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

def picard_fixmate(picard, align_bam):
    """Run Picard's FixMateInformation generating an aligned output file.
    """
    base, ext = os.path.splitext(align_bam)
    out_file = "%s-sort%s" % (base, ext)
    if not os.path.exists(out_file):
        with curdir_tmpdir() as tmp_dir:
            opts = [("INPUT", align_bam),
                    ("OUTPUT", out_file),
                    ("TMP_DIR", tmp_dir),
                    ("SORT_ORDER", "coordinate")]
            picard.run("FixMateInformation", opts)
    return out_file

if __name__ == "__main__":
    main(*sys.argv[1:])
