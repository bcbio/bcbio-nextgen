#!/usr/bin/env python
"""Convert Complete Genomics SvEvents file of structural variants to VCF.

Handles:

  inversion/probable-inversion -> INV
  deletion -> DEL
  tandem-duplication -> DUP
  distal-duplication -> Breakends (BND)

Does not convert: complex

Requirements:

bx-python: https://bitbucket.org/james_taylor/bx-python/wiki/Home

Usage:
  cg_svevents_to_vcf.py <SV Events TSV file> <Genome in UCSC 2bit format>
"""
import sys
import csv
from collections import namedtuple

from bx.seq import twobit

def main(svevents_file, genome_file):
    genome_2bit = twobit.TwoBitFile(open(genome_file))
    for event in svevent_reader(svevents_file):
        for vcf_line in _svevent_to_vcf(event):
            print vcf_line

# ## Convert different types of svEvents into VCF info

VcfLine = namedtuple('VcfLine', ["chrom", "pos", "id", "ref", "alt", "info"])

def _svevent_to_vcf(event):
    if event["Type"] in ["inversion", "probable-inversion"]:
        out = _convert_event_inv(event)
    elif event["Type"] in ["deletion"]:
        out = _convert_event_del(event)
    elif event["Type"] in ["tandem-duplication"]:
        out = _convert_event_dup(event)
    elif event["Type"] in ["distal-duplication"]:
        out = _convert_event_bnd(event)
    elif event["Type"] in ["complex"]:
        out = [] # ignore complex events
    else:
        raise ValueError("Unexpected event type %s" % event["Type"])
    return out

def _convert_event_inv(event):
    print event
    return []

def _convert_event_del(event):
    return []

def _convert_event_dup(event):
    return []

def _convert_event_bnd(event):
    return []

def svevent_reader(in_file):
    """Lazy generator of SV events, returned as dictionary of parts.
    """
    with open(in_file) as in_handle:
        while 1:
            line = next(in_handle)
            if line.startswith(">"):
                break
        header = line[1:].rstrip().split("\t")
        reader = csv.reader(in_handle, dialect="excel-tab")
        for parts in reader:
            out = {}
            for h, p in zip(header, parts):
                out[h] = p
            yield out

if __name__ == "__main__":
    main(*sys.argv[1:])
