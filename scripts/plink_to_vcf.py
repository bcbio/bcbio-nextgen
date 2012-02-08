#!/usr/bin/env python
"""Convert Plink ped/map files into VCF format using plink and Plink/SEQ.

Requires:

plink: http://pngu.mgh.harvard.edu/~purcell/plink/
PLINK/SEQ: http://atgu.mgh.harvard.edu/plinkseq/
bx-python: https://bitbucket.org/james_taylor/bx-python/wiki/Home

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
    pbed_prefix = convert_to_plink_bed(ped_file, map_file)
    vcf_file = convert_bed_to_vcf(pbed_prefix, ped_file)
    fix_nonref_positions(vcf_file, ref_file)

def convert_to_plink_bed(ped_file, map_file):
    plink_cl = "p-link" # from ubuntu package, 'plink' otherwise
    plink_prefix = os.path.splitext(os.path.basename(ped_file))[0]
    work_dir = os.path.join(os.path.dirname(ped_file), "vcfconvert")
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    out_base = os.path.join(work_dir, plink_prefix)
    if not os.path.exists("{}.bed".format(out_base)):
        subprocess.check_call([plink_cl, "--ped", ped_file, "--map", map_file,
                               "--make-bed", "--out", out_base])
    return out_base

def convert_bed_to_vcf(pbed_prefix, ped_file):
    out_file = "{}.vcf".format(os.path.splitext(ped_file)[0])
    if not os.path.exists(out_file):
        subprocess.check_call(["pseq", pbed_prefix, "new-project"])
        subprocess.check_call(["pseq", pbed_prefix, "load-plink",
                               "--file", pbed_prefix, "--id", "vcfconvert"])
        with open(out_file, "w") as out_handle:
            subprocess.check_call(["pseq", pbed_prefix, "write-vcf"],
                                  stdout=out_handle)
    return out_file

def fix_nonref_positions(in_file, ref_file):
    """Fix Genotyping VCF positions where the bases are all variants.

    The plink/pseq output does not handle these correctly, and
    has all reference/variant bases reversed.
    """
    swap = {"1/1": "0/0", "0/1": "0/1", "0/0": "1/1", "./.": "./."}
    complements = {"G": "C", "A": "T", "C": "G", "T": "A"}
    ref2bit = twobit.TwoBitFile(open(ref_file))
    out_file = apply("{}-fix{}".format, os.path.splitext(in_file))

    with open(in_file) as in_handle:
        with open(out_file, "w") as out_handle:
            for line in in_handle:
                if line.startswith("#"):
                    out_handle.write(line)
                else:
                    parts = line.rstrip("\t\n").split("\t")
                    pos = int(parts[1])
                    ref_base = ref2bit[parts[0]].get(pos-1, pos).upper()
                    ref, var = parts[3:5]
                    if ref_base == var or ref == "N":
                        base = parts[:9]
                        base[3] = var
                        base[4] = ref
                        swap_genotypes = [swap[x] for x in parts[9:]]
                        parts = base + swap_genotypes
                        #print "swapped", pos
                    elif ref_base == ref:
                        #print "same", pos
                        assert ref_base == ref, (ref_base, parts)
                    elif ref_base != ref and var in ["N", "0"]:
                        parts[3] = ref_base
                    elif ref_base != ref and complements[ref] == ref_base:
                        parts[3] = complements[ref]
                        parts[4] = complements[var]
                    else:
                        raise ValueError("Cannot associate ref {} with line: {}".format(
                            ref_base, parts))
                    out_handle.write("\t".join(parts) + "\n")
        return out_file

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "Incorrect arguments"
        print __doc__
        sys.exit(1)
    main(*sys.argv[1:])
