#!/usr/bin/env python
"""Format DREAM challenge truth sets to contain BED files of covered regions and SVs.

Prepares inputs for evaluation with bcbio-nextgen.
Usage:
  format_dream_truthset.py <truth_file> <genome FAI file>

"""
import glob
import gzip
import subprocess
import sys

import pybedtools

from bcbio import utils

def main(vcf_file, fai_file):
    out, out_vcf, out_bed = _prep_out_handles(vcf_file)
    nocall_regions = []
    with gzip.open(vcf_file) as in_handle:
        for line in in_handle:
            if line.startswith("#"):
                out["vcf"].write(line)
            else:
                nocall_region = _write_event(line, out)
                if nocall_region:
                    nocall_regions.append(nocall_region)
    for h in out.values():
        h.close()
    _write_outbed(nocall_regions, fai_file, out_bed)
    subprocess.check_call(["bgzip", "-f", out_vcf])
    subprocess.check_call(["tabix", "-f", "-p", "vcf", out_vcf + ".gz"])
    base = utils.splitext_plus(out_vcf)[0]
    subprocess.check_call(["rm", "-f", "%s.tar.gz" % base])
    subprocess.check_call(["tar", "-czvpf", "%s.tar.gz" % base] +
                          glob.glob("%s*.vcf.gz*" % base) + glob.glob("%s*.bed" % base))

def _write_outbed(nocall_regions, fai_file, out_file):
    ref_lines = []
    with open(fai_file) as in_handle:
        for line in in_handle:
            parts = line.split()
            ref_lines.append("%s\t0\t%s" % (parts[0], parts[1]))
    ref_bedtool = pybedtools.BedTool("\n".join(ref_lines), from_string=True)
    nocall_bedtool = pybedtools.BedTool("\n".join("\t".join(r) for r in nocall_regions), from_string=True)
    ref_bedtool.subtract(nocall_bedtool).saveas(out_file)

def _write_event(line, out):
    parts = line.split("\t")
    if parts[4].startswith("<"):
        event = parts[4][1:-1]
        if event in out:
            bed_parts = _parts_to_bed(parts)
            out[event].write("%s\n" % "\t".join(bed_parts))
        else:
            assert event in ["MSK", "IGN"], event
            return _parts_to_bed(parts)
    else:
        out["vcf"].write(line)

def _parts_to_bed(parts):
    chrom = parts[0]
    start = str(int(parts[1]) - 1)
    ends = [x for x in parts[7].split(";") if x.startswith("END=")]
    assert len(ends) == 1, parts[7]
    end = str(int(ends[0].replace("END=", "")))
    return [chrom, start, end]

def _prep_out_handles(base_file):
    base = utils.splitext_plus(base_file)[0].replace(".", "_")
    out_bed = "%s_regions.bed" % base
    out_vcf = "%s.vcf" % base
    fs = {"vcf": out_vcf}
    for event in ["DEL", "DUP", "INS", "INV"]:
        fs[event] = "%s_sv_%s.bed" % (base, event)
    out = {}
    for key, f in fs.items():
        out[key] = open(f, "w")
    return out, out_vcf, out_bed

if __name__ == "__main__":
    main(*sys.argv[1:])
