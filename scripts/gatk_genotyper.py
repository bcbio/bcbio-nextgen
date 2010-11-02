#!/usr/bin/env python
"""Provide SNP and indel calling using GATK genotyping tools.

http://www.broadinstitute.org/gsa/wiki/index.php/Unified_genotyper
http://www.broadinstitute.org/gsa/wiki/index.php/Local_realignment_around_indels
http://www.broadinstitute.org/gsa/wiki/index.php/IndelGenotyper
http://www.broadinstitute.org/gsa/wiki/index.php/VariantFiltrationWalker
http://www.broadinstitute.org/gsa/wiki/index.php/VariantEval_v2.0


Usage:
    gatk_genotyper.py <config file> <reference file> <align BAM file>
                      <dbsnp file>

Requires:
  - Picard
  - GATK
  - samtools
"""
import os
import sys
import subprocess

import yaml

from bcbio.picard import PicardRunner
from bcbio.picard.utils import curdir_tmpdir

def main(config_file, ref_file, align_bam, dbsnp=None):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    picard = PicardRunner(config["program"]["picard"])
    ref_dict = index_ref_file(picard, ref_file)
    index_bam(align_bam, config["program"]["samtools"])
    realign_target_file = realigner_targets(picard, align_bam,
            ref_file, dbsnp)
    realign_bam = indel_realignment(picard, align_bam, ref_file,
            realign_target_file)
    realign_sort_bam = picard_fixmate(picard, realign_bam)
    index_bam(realign_sort_bam, config["program"]["samtools"])
    snp_file = unified_genotyper(picard, realign_sort_bam, ref_file,
            config["algorithm"]["platform"], dbsnp)
    filter_snp = variant_filtration(picard, snp_file, ref_file)
    #eval_snp = variant_eval(picard, filter_snp, ref_file, dbsnp)

def unified_genotyper(picard, align_bam, ref_file, platform, dbsnp=None):
    """Perform SNP genotyping on the given alignment file.
    """
    out_file = "%s-snp.vcf" % os.path.splitext(align_bam)[0]
    params = ["-T", "UnifiedGenotyper",
              "-I", align_bam,
              "-R", ref_file,
              "-o", out_file,
              "--genotype_model", "JOINT_ESTIMATE",
              "--base_model", "EMPIRICAL",
              "--standard_min_confidence_threshold_for_calling", "10.0",
              "--standard_min_confidence_threshold_for_emitting", "10.0",
              "--trigger_min_confidence_threshold_for_calling", "10.0",
              "--trigger_min_confidence_threshold_for_emitting", "10.0",
              "--downsample_to_coverage", 10000,
              "--min_base_quality_score", 20,
              "--platform", platform,
              "-l", "INFO",
              ]
    if dbsnp:
        params += ["-B", "dbsnp,VCF,%s" % dbsnp]
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        picard.run_gatk(params)
    return out_file

def variant_filtration(picard, snp_file, ref_file):
    """Filter out problematic SNP calls.

    XXX missing:
        interval list
    """
    out_file = "%s-filter%s" % os.path.splitext(snp_file)
    params = ["-T", "VariantFiltration",
              "-R", ref_file,
              "-o", out_file,
              "-B", "variant,VCF,%s" % snp_file,
              "--filterName", "QUALFilter",
              "--filterExpression", "QUAL <= 50.0",
              "--filterName", "QDFilter",
              "--filterExpression", "QD < 5.0",
              "--filterName", "ABFilter",
              "--filterExpression", "AB > 0.75",
              "--filterName", "HRunFilter",
              "--filterExpression", "HRun > 3.0",
              "-l", "INFO",
              ]
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        picard.run_gatk(params)
    return out_file

def NOTUSED_variant_eval(picard, filter_snp, ref_file, dbsnp):
    """Provide summary evaluating called variants.

    XXX missing:
        interval list
    """
    out_file = "%s-eval.txt" % os.path.splitext(filter_snp)[0]
    params = ["-T", "VariantEval",
              "-B", "eval,VCF,%s" % filter_snp,
              "-R", ref_file,
              "-o", out_file,
              "-l", "INFO",
              ]
    if dbsnp:
        params += ["-B", "comp,VCF,%s" % dbsnp]
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        picard.run_gatk(params)
    return out_file

def realigner_targets(picard, align_bam, ref_file, dbsnp=None):
    """Generate a list of interval regions for realignment around indels.
    """
    out_file = "%s-realign.intervals" % os.path.splitext(align_bam)[0]
    params = ["-T", "RealignerTargetCreator",
              "-I", align_bam,
              "-R", ref_file,
              "-o", out_file,
              "-l", "INFO",
              ]
    if dbsnp:
        params += ["-B", "dbsnp,VCF,%s" % dbsnp]
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        picard.run_gatk(params)
    return out_file

def indel_realignment(picard, align_bam, ref_file, intervals):
    """Perform realignment of BAM file in specified regions
    """
    out_file = "%s-realign.bam" % os.path.splitext(align_bam)[0]
    params = ["-T", "IndelRealigner",
              "-I", align_bam,
              "-R", ref_file,
              "-targetIntervals", intervals,
              "-o", out_file,
              "-l", "INFO",
              ]
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        with curdir_tmpdir() as tmp_dir:
            picard.run_gatk(params, tmp_dir)
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

def index_ref_file(picard, ref_file):
    """Provide a Picard style dict index file for a reference genome.
    """
    dict_file = "%s.dict" % os.path.splitext(ref_file)[0]
    if not os.path.exists(dict_file):
        opts = [("REFERENCE", ref_file),
                ("OUTPUT", dict_file)]
        picard.run("CreateSequenceDictionary", opts)
    return dict_file

def index_bam(bam_file, samtools_cmd):
    index_file = "%s.bai" % bam_file
    if not os.path.exists(index_file):
        cl = [samtools_cmd, "index", bam_file]
        subprocess.check_call(cl)
    return index_file

if __name__ == "__main__":
    main(*sys.argv[1:])
