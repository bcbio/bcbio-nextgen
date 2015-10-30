#!/usr/bin/env python
"""Collapse HLAs present in hg38 1000 genomes distribution to p-groups.

P-groups have identical proteins in antigen binding domains:

http://hla.alleles.org/alleles/p_groups.html
"""
import collections
import pprint
import sys

def main(fasta_fai, pgroup_file):
    hlas = read_hlas(fasta_fai)
    pgroups = read_pgroups(pgroup_file)

    groups = collections.defaultdict(list)
    for orig_hla in hlas:
        #print list(hla_choices(orig_hla))
        found = False
        for cur_hla in hla_choices(orig_hla[:]):
            if cur_hla in pgroups:
                found = True
                groups[pgroups[cur_hla]].append(orig_hla)
                break
        if not found:
            groups[""].append(orig_hla)
            print("Did not find group for %s" % orig_hla)
    final = {}
    for group in sorted(groups.keys()):
        for hla in groups[group]:
            print group, hla
            final[hla] = group
    pprint.pprint(final)

def hla_choices(orig_hla, min_parts=2):
    """Provide a range of options for HLA type, with decreassing resolution.
    """
    yield orig_hla
    try:
        int(orig_hla[-1])
    except ValueError:
        yield orig_hla[:-1]
    hla_parts = orig_hla.split(":")
    for sub_i in range(len(hla_parts) - min_parts + 1):
        yield ":".join(hla_parts[:len(hla_parts) - sub_i])

def read_pgroups(in_file):
    """Read HLAs and the pgroups they fall in.
    """
    out = {}
    with open(in_file) as in_handle:
        for line in (l for l in in_handle if not l.startswith("#")):
            locus, alleles, group = line.strip().split(";")
            for allele in alleles.split("/"):
                out["HLA-%s%s" % (locus, allele)] = group
    return out

def read_hlas(fasta_fai):
    """Get HLA alleles from the hg38 fasta fai file.
    """
    out = []
    with open(fasta_fai) as in_handle:
        for line in in_handle:
            if line.startswith("HLA"):
                out.append(line.split()[0])
    return out

if __name__ == "__main__":
    main(*sys.argv[1:])