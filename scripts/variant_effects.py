#!/usr/bin/env python
"""Compute the effects of SNP modifications using snpEff.

http://sourceforge.net/projects/snpeff/

Usage:
    variant_effects.py <snpEff jar> <VCF input> <genome name> [<interval file>]
"""
import os
import sys
import csv
import glob
import subprocess

# remap Galaxy genome names to the ones used by snpEff
GENOME_REMAP = {
        "GRCh37": "hg37.60",
        "hg19" : "hg37.60",
        "mm9" : "mm37.60",
        }

def main(snpeff_jar, vcf_ref, genome, interval_file=None):
    snpeff_config = "%s.config" % os.path.splitext(snpeff_jar)[0]
    genome = GENOME_REMAP[genome]
    if os.path.isdir(vcf_ref):
        vcf_files = sorted(glob.glob(os.path.join(vcf_ref, "*-snp-filter.vcf")))
    else:
        vcf_files = [vcf_ref]
    for vcf_in in (v for v in vcf_files if _vcf_has_items(v)):
        se_interval = (convert_to_snpeff_interval(interval_file, vcf_ref)
                       if interval_file else None)
        try:
            out_file = run_snpeff(vcf_in, genome, snpeff_jar, snpeff_config,
                    se_interval)
        finally:
            for fname in [se_interval]:
                if fname and os.path.exists(fname):
                    os.remove(fname)

def _vcf_has_items(in_file):
    with open(in_file) as in_handle:
        while 1:
            line = in_handle.next()
            if not line.startswith("#"):
                break
        line = in_handle.next()
        if line:
            return True
    return False

def run_snpeff(snp_in, genome, snpeff_jar, snpeff_config, se_interval):
    out_file = "%s-effects.tsv" % (os.path.splitext(snp_in)[0])
    if not os.path.exists(out_file):
        cl = ["java", "-jar", snpeff_jar, "-1", "-vcf4", "-pass", "-c", snpeff_config,
              genome, snp_in]
        if se_interval:
            cl.extend(["-filterInterval", se_interval])
        print " ".join(cl)
        with open(out_file, "w") as out_handle:
            subprocess.check_call(cl, stdout=out_handle)
    return out_file

def convert_to_snpeff_interval(in_file, base_file):
    out_file = "%s-snpeff-intervals.bed" % os.path.splitext(base_file)[0]
    if not os.path.exists(out_file):
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle, dialect="excel-tab")
            with open(in_file) as in_handle:
                for line in (l for l in in_handle if not l.startswith("@")):
                    parts = line.split()
                    writer.writerow(parts[:3])
    return out_file

if __name__ == "__main__":
    main(*sys.argv[1:])
