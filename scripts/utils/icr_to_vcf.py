#!/usr/bin/env python
"""Convert supplemental table 1 from the ICR142 validation series into VCF and BED.

Creates a VCF truth set and BED region file for each sample for running
validations with standard tools.

http://f1000research.com/articles/5-386/v1

"""
import collections
import csv
import math
import sys

import pandas as pd

def main(in_file, hom_file):
    call_buffer = 1000
    val_buffer = 3
    homs = _read_hom_file(hom_file)
    df = pd.read_csv(in_file, sep="\t")
    grouped = df.groupby(["Sample"])
    grouped.apply(_write_sample_vcf_bed(call_buffer, val_buffer, homs))

def _read_hom_file(hom_file):
    """Read file of expected homozygotes based on initial screen. Rest are hets.
    """
    homs = collections.defaultdict(set)
    with open(hom_file) as in_handle:
        reader = csv.reader(in_handle)
        for sample, chrom, pos, ref, alt in reader:
            homs[sample].add((chrom, pos, ref, alt))
    return homs

def _write_vcf_header(handle, sample):
    handle.write("##fileformat=VCFv4.1\n")
    handle.write("""##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n""")
    handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample)

def chrom_sort(row):
    try:
        return int(row[0])
    except ValueError:
        return row[0]

def _write_sample_vcf_bed(call_buffer, val_buffer, homs):
    def run(df):
        df = df.sort_values(["Chr", "EvaluatedPosition"])
        sample = list(df["Sample"])[0]
        vcf_file = "ICR142-%s.vcf" % list(df["Sample"])[0]
        call_file = "ICR142-%s-call.bed" % list(df["Sample"])[0]
        val_file = "ICR142-%s-validate.bed" % list(df["Sample"])[0]
        with open(vcf_file, "w") as vcf_handle:
            vcf_rows = []
            with open(call_file, "w") as call_handle:
                with open(val_file, "w") as val_handle:
                    for _, row in df.iterrows():
                        pos = int(row.EvaluatedPosition if math.isnan(row.POS) else row.POS) - 1
                        call_handle.write("%s\t%s\t%s\t%s\n" % (row.Chr, pos - call_buffer, pos + call_buffer, row.Gene))
                        val_handle.write("%s\t%s\t%s\t%s\n" % (row.Chr, pos - val_buffer, pos + val_buffer, row.Gene))
                        if row.SangerCall != "No":
                            if math.isnan(row.POS):
                                print list(df["Sample"])[0], len(df), row.Gene
                                pos, ref, alt = None, None, None
                            else:
                                pos, ref, alt = (int(row.POS), row.REF, row.ALT)
                            if pos is not None:
                                key = (row.Chr, str(pos), ref, alt)
                                gt = ("1/1" if key in homs[sample] else "0/1")
                                vcf_rows.append([row.Chr, str(pos), ".", ref, alt, ".", "PASS", ".", "GT", gt])
            _write_vcf_header(vcf_handle, list(df["Sample"])[0])
            for vcf_row in sorted(vcf_rows, key=chrom_sort):
                vcf_handle.write("\t".join(vcf_row) + "\n")
    return run

if __name__ == "__main__":
    main(*sys.argv[1:])
