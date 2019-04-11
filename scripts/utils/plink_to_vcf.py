#!/usr/bin/env python
"""Convert Plink ped/map files into VCF format using plink and Plink/SEQ.

Latest version available as part of bcbio-nextgen:
https://github.com/bcbio/bcbio-nextgen/blob/master/scripts/plink_to_vcf.py

Requires:

plink: http://pngu.mgh.harvard.edu/~purcell/plink/
PLINK/SEQ: https://atgu.mgh.harvard.edu/plinkseq/
bx-python: https://github.com/bxlab/bx-python

You also need the genome reference file in 2bit format:
http://genome.ucsc.edu/FAQ/FAQformat.html#format7
using faToTwoBit:
http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/

Usage:
  plink_to_vcf.py <ped file> <map file> <UCSC reference file in 2bit format)

"""
import os
import sys
import subprocess

from bx.seq import twobit

def main(ped_file, map_file, ref_file):
    base_dir = os.getcwd()
    pbed_prefix = convert_to_plink_bed(ped_file, map_file, base_dir)
    vcf_file = convert_bed_to_vcf(pbed_prefix, ped_file, base_dir)
    fix_nonref_positions(vcf_file, ref_file)

def convert_to_plink_bed(ped_file, map_file, base_dir):
    # from ubuntu package, 'plink' otherwise
    for plink_cl in ["p-link", "plink", "plink2"]:
        try:
            subprocess.check_call([plink_cl, "--help"])
            break
        except:
            pass
    plink_prefix = os.path.splitext(os.path.basename(ped_file))[0].replace(".", "_")
    work_dir = os.path.join(base_dir, "vcfconvert")
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    out_base = os.path.join(work_dir, plink_prefix)
    if not os.path.exists("{0}.bed".format(out_base)):
        subprocess.check_call([plink_cl, "--noweb",
                               "--ped", ped_file, "--map", map_file,
                               "--make-bed", "--out", out_base])
    return out_base

def convert_bed_to_vcf(pbed_prefix, ped_file, base_dir):
    out_file = os.path.join(base_dir,
                            "{0}-raw.vcf".format(os.path.splitext(os.path.basename(ped_file))[0]))
    if not os.path.exists(out_file):
        subprocess.check_call(["pseq", pbed_prefix, "new-project"])
        subprocess.check_call(["pseq", pbed_prefix, "load-plink",
                               "--file", pbed_prefix, "--id", "vcfconvert"])
        with open(out_file, "w") as out_handle:
            subprocess.check_call(["pseq", pbed_prefix, "write-vcf"],
                                  stdout=out_handle)
    return out_file

def fix_line_problems(parts):
    """Fix problem alleles and reference/variant bases in VCF line.
    """
    varinfo = parts[:9]
    genotypes = []
    # replace haploid calls
    for x in parts[9:]:
        if len(x) == 1:
            x = "./."
        genotypes.append(x)
    if varinfo[3] == "0": varinfo[3] = "N"
    if varinfo[4] == "0": varinfo[4] = "N"
    return varinfo, genotypes

def fix_vcf_line(parts, ref_base):
    """Orient VCF allele calls with respect to reference base.

    Handles cases with ref and variant swaps. strand complements.
    """
    swap = {"1/1": "0/0", "0/1": "0/1", "0/0": "1/1", "./.": "./."}
    complements = {"G": "C", "A": "T", "C": "G", "T": "A", "N": "N"}
    varinfo, genotypes = fix_line_problems(parts)
    ref, var = varinfo[3:5]
    # non-reference regions or non-informative, can't do anything
    if ref_base in [None, "N"] or set(genotypes) == set(["./."]):
        varinfo = None
    # matching reference, all good
    elif ref_base == ref:
        assert ref_base == ref, (ref_base, parts)
    # swapped reference and alternate regions
    elif ref_base == var or ref in ["N", "0"]:
        varinfo[3] = var
        varinfo[4] = ref
        genotypes = [swap[x] for x in genotypes]
    # reference is on alternate strand
    elif ref_base != ref and complements.get(ref) == ref_base:
        varinfo[3] = complements[ref]
        varinfo[4] = ",".join([complements[v] for v in var.split(",")])
    # unspecified alternative base
    elif ref_base != ref and var in ["N", "0"]:
        varinfo[3] = ref_base
        varinfo[4] = ref
        genotypes = [swap[x] for x in genotypes]
    # swapped and on alternate strand
    elif ref_base != ref and complements.get(var) == ref_base:
        varinfo[3] = complements[var]
        varinfo[4] = ",".join([complements[v] for v in ref.split(",")])
        genotypes = [swap[x] for x in genotypes]
    else:
        print "Did not associate ref {0} with line: {1}".format(
            ref_base, varinfo)
    if varinfo is not None:
        return varinfo + genotypes

def fix_nonref_positions(in_file, ref_file):
    """Fix Genotyping VCF positions where the bases are all variants.

    The plink/pseq output does not handle these correctly, and
    has all reference/variant bases reversed.
    """
    ignore_chrs = ["."]
    ref2bit = twobit.TwoBitFile(open(ref_file))
    out_file = in_file.replace("-raw.vcf", ".vcf")

    with open(in_file) as in_handle:
        with open(out_file, "w") as out_handle:
            for line in in_handle:
                if line.startswith("#"):
                    out_handle.write(line)
                else:
                    parts = line.rstrip("\r\n").split("\t")
                    pos = int(parts[1])
                    # handle chr/non-chr naming
                    if parts[0] not in ref2bit.keys() and parts[0].replace("chr", "") in ref2bit.keys():
                        parts[0] = parts[0].replace("chr", "")
                    # handle X chromosome
                    elif parts[0] not in ref2bit.keys() and parts[0] == "23":
                        for test in ["X", "chrX"]:
                            if test in ref2bit.keys():
                                parts[0] == test
                    ref_base = None
                    if parts[0] not in ignore_chrs:
                        try:
                            ref_base = ref2bit[parts[0]].get(pos-1, pos).upper()
                        except Exception as msg:
                            print "Skipping line. Failed to retrieve reference base for %s\n%s" % (str(parts), msg)
                    parts = fix_vcf_line(parts, ref_base)
                    if parts is not None:
                        out_handle.write("\t".join(parts) + "\n")
        return out_file

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "Incorrect arguments"
        print __doc__
        sys.exit(1)
    main(*sys.argv[1:])
