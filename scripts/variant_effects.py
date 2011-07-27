#!/usr/bin/env python
"""Compute the effects of SNP modifications using snpEff.

http://sourceforge.net/projects/snpeff/

Usage:
    variant_effects.py <snpEff jar> <VCF input> <genome name> [<interval file>]
"""
import os
import sys
import glob

from bcbio.variation.effects import snpeff_effects

def main(snpeff_jar, vcf_ref, genome, interval_file=None):
    if os.path.isdir(vcf_ref):
        vcf_files = sorted(glob.glob(os.path.join(vcf_ref, "*-snp-filter.vcf")))
    else:
        vcf_files = [vcf_ref]
    for vcf_file in vcf_files:
        snpeff_effects(snpeff_jar, vcf_file, genome, interval_file)

if __name__ == "__main__":
    main(*sys.argv[1:])
