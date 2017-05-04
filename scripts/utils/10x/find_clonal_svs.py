#!/usr/bin/env python
"""Find 10x structural variants present uniquely in parent or clones.

Includes fixes for 10x SVs to allow usage with snpEff and downstream VCF
parsing libraries.
"""
import math
import os
import sys

import cyvcf2

def main(in_file):
    in_cyvcf = cyvcf2.VCF(in_file)
    parental_file = "%s-parental.vcf" % os.path.splitext(in_file)[0]
    parental_handle = open(parental_file, "w")
    clones_file = "%s-clones.vcf" % os.path.splitext(in_file)[0]
    clones_handle = open(clones_file, "w")
    parental_handle.write(_fix_header(in_cyvcf.raw_header))
    clones_handle.write(_fix_header(in_cyvcf.raw_header))
    for rec in in_cyvcf:
        calls = [parse_name(s) for s, gt in zip(in_cyvcf.samples, rec.gt_bases) if _has_call(gt)]
        if len(calls) == 1 and calls[0] == "parental":
            parental_handle.write(_fix_start_end(rec))
            _summarize("parental", rec)
        elif len(calls) > 0 and "parental" not in calls:
            clones_handle.write(_fix_start_end(rec))
            _summarize("clones", rec)

def _fix_header(header):
    fix = "##FORMAT=<ID=DR,"
    out = []
    for line in header.split("\n"):
        if line.startswith(fix):
            line = line.replace("Number=1", "Number=.")
        out.append(line)
    return "\n".join(out)

def _fix_start_end(rec):
    line = str(rec)
    parts = line.split("\t")
    start = int(parts[1])
    info = dict(p.split("=") for p in parts[7].split(";"))
    if "SVLEN" in info and abs(int(math.round(info["SVLEN"]))) > 1e7:
        return ""
    end = int(info["END"])
    if end < start:
        info["END"] = str(start)
        parts[7] = ";".join(["%s=%s" % (k, v) for (k, v) in info.items()])
        parts[1] = str(end)
        return "\t".join(parts)
    else:
        return line

def _has_call(g):
    gts = g.split("|")
    if len(gts) == 1:
        gts = g.split("/")
    return any([gt not in [".", "N"] for gt in gts])

def _summarize(classification, rec):
    print(classification, rec.CHROM, rec.POS, rec.INFO["SVTYPE"], rec.INFO["AVGLEN"])

def parse_name(x):
    return os.path.basename(x).split("-")[0]


if __name__ == "__main__":
    main(*sys.argv[1:])
