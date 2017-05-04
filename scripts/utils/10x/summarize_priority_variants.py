#!/usr/bin/env python
"""Summarize priority calls in annotated structural variants.
"""
import csv
import os
import sys

import cyvcf2

def main(in_file):
    in_cyvcf = cyvcf2.VCF(in_file)
    writer = csv.writer(sys.stdout)
    writer.writerow(["chrom", "start", "end", "svtype", "samples", "size", "gene", "annotation", "detail"])
    for rec in in_cyvcf:
        calls = [parse_name(s) for s, gt in zip(in_cyvcf.samples, rec.gt_bases) if _has_call(gt)]
        anns = [x.split("|") for x in rec.INFO.get("SIMPLE_ANN", "").split(",")]
        svtypes = set([])
        all_genes = set([])
        annotations = set([])
        details = set([])
        for svtype, annotation, genes, _, detail, _ in (x for x in anns if x and len(x) > 1):
            if detail != "NOT_PRIORITISED":
                for c in "'[]' ":
                    svtype = svtype.replace(c, "")
                svtypes.add(svtype)
                all_genes.add(genes)
                annotations.add(annotation)
                details.add(detail)
        if svtypes:
            start = int(rec.POS)
            end = rec.INFO.get("END")
            size = abs(rec.INFO.get("SVLEN", end - start))
            writer.writerow([rec.CHROM, start, end, _combine(svtypes), size, ";".join(calls), _combine(all_genes),
                             _combine(annotations), _combine(details)])

def _combine(xs):
    return ";".join(sorted(list(xs)))

def _has_call(g):
    gts = g.split("|")
    if len(gts) == 1:
        gts = g.split("/")
    return any([gt not in [".", "N"] for gt in gts])

def parse_name(x):
    return os.path.basename(x).split("-")[0]

if __name__ == "__main__":
    main(*sys.argv[1:])
