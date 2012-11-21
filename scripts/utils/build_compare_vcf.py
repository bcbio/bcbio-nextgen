#!/usr/bin/env python
"""Build a test comparison dataset from an existing VCF file.

Produces a slightly modified comparison set to use for testing comparison
software.

Usage:
  build_compare_vcf.py <in_vcf_file>
"""
import os
import sys
import random

import vcf

def main(in_file):
    out_file = apply("{0}-cmp{1}".format, os.path.splitext(in_file))
    with open(in_file) as in_handle:
        with open(out_file, "w") as out_handle:
            rdr = vcf.Reader(in_handle)
            wtr = vcf.Writer(out_handle, rdr)
            for rec in rdr:
                out_rec = adjust_variant(rec)
                if out_rec:
                    wtr.write_record(out_rec)

def adjust_variant(rec):
    do_change = random.random()
    if do_change < 0.2:
        return None
    elif do_change < 0.5:
        return rec
    else:
        rec.samples = [adjust_genotype(g) for g in rec.samples]
        return rec

def adjust_genotype(g):
    alts = ["0", "1"]
    do_change = random.random()
    if do_change < 0.7:
        new_gt = None
    elif do_change < 0.9:
        new_gt = g.gt_phase_char().join(["."] * (len(g.gt_alleles)))
    else:
        new_gt = g.gt_phase_char().join([random.choice(alts) for x in g.gt_alleles])
    if new_gt:
        g.data = g.data._replace(GT=new_gt)
    return g

if __name__ == "__main__":
    main(sys.argv[1])
