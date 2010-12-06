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
import collections

# remap Galaxy genome names to the ones used by snpEff
GENOME_REMAP = {
        "GRCh37": "hg37.60",
        }

def main(snpeff_jar, vcf_ref, genome, interval_file=None):
    snpeff_config = "%s.config" % os.path.splitext(snpeff_jar)[0]
    genome = GENOME_REMAP[genome]
    #allowed_chroms = snpeff_allowed_chroms(genome, snpeff_config)
    #allowed_bases = get_allowed_bases(interval_file) if interval_file else None
    if os.path.isdir(vcf_ref):
        vcf_files = sorted(glob.glob(os.path.join(vcf_ref, "*-snp-filter.vcf")))
    else:
        vcf_files = [vcf_ref]
    for vcf_in in vcf_files:
        se_interval = (convert_to_snpeff_interval(interval_file, vcf_ref)
                       if interval_file else None)
        #snp_in = vcf_to_snpeff(vcf_in, allowed_chroms, allowed_bases)
        try:
            out_file = run_snpeff(vcf_in, genome, snpeff_jar, snpeff_config,
                    se_interval)
        finally:
            for fname in [se_interval]:
                if fname and os.path.exists(fname):
                    os.remove(fname)

def run_snpeff(snp_in, genome, snpeff_jar, snpeff_config, se_interval):
    out_file = "%s-effects.tsv" % (os.path.splitext(snp_in)[0])
    if not os.path.exists(out_file):
        cl = ["java", "-jar", snpeff_jar, "-1", "-vcf4", "-pass", "-c", snpeff_config,
              genome, snp_in]
        if se_interval:
            cl.extend(["-filterInterval", se_interval])
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

def snpeff_allowed_chroms(genome, snpeff_config):
    """Retrieve allowed chromosomes from snpEff configuration file.
    """
    with open(snpeff_config) as in_handle:
        for line in (l.lstrip() for l in in_handle):
            if line.startswith("%s.chromosomes" % genome):
                _, chroms = line.split(" : ")
                return [c.strip() for c in chroms.split(",")]
    raise ValueError("Did not find chromosomes for %s in %s" % (genome,
        snpeff_config))

def get_allowed_bases(interval_file):
    allowed = collections.defaultdict(list)
    with open(interval_file) as in_handle:
        for line in (l for l in in_handle if not l.startswith("@")):
            chrom, start, end, _, _ = line.rstrip("\r\n").split("\t")
            allowed[chrom].extend(range(int(start)-1, int(end)+1))
    final = dict()
    for chrom, bases in allowed.iteritems():
        final[chrom] = set(bases)
    return final

def vcf_to_snpeff(vcf_in, allowed_chroms, allowed_bases):
    """Convert a VCF input file into snpEff's input format.
    """
    out_file = "%s.insnpeff" % os.path.splitext(vcf_in)[0]
    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle, dialect="excel-tab")
        for chrom, pos, ref, alt in _vcf_items(vcf_in):
            if chrom in allowed_chroms:
                if allowed_bases is None or pos in allowed_bases[chrom]:
                    writer.writerow([chrom, pos, ref, alt, "+"])
    return out_file

def _vcf_items(vcf_in):
    with open(vcf_in) as in_handle:
        for line in (l for l in in_handle if not l.startswith("#")):
            parts = line.split("\t")
            if parts[6].upper() == "PASS":
                yield (parts[0], int(parts[1]), parts[3], parts[4])

if __name__ == "__main__":
    main(*sys.argv[1:])
