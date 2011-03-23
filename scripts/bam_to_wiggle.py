#!/usr/bin/env python
"""Convert BAM files to BigWig file format in a specified region.

Usage:
    bam_to_wiggle.py <BAM file> [<YAML config>]
    [--outfile=<output file name>
     --chrom=<chrom>
     --start=<start>
     --end=<end>]

chrom start and end are optional, in which case they default to everything.

The config file is in YAML format and specifies the location of the wigToBigWig
program from UCSC:

program:
  ucsc_bigwig: wigToBigWig

If not specified, these will be assumed to be present in the system path.

The script requires:
    pysam (http://code.google.com/p/pysam/)
    wigToBigWig from UCSC (http://hgdownload.cse.ucsc.edu/admin/exe/)
If a configuration file is used, then PyYAML is also required (http://pyyaml.org/)
"""
import os
import sys
import subprocess
from optparse import OptionParser
from contextlib import contextmanager

import pysam

def main(bam_file, config_file=None, chrom='all', start=0, end=None,
         outfile=None):
    if config_file:
        import yaml
        with open(config_file) as in_handle:
            config = yaml.load(in_handle)
    else:
        config = {"program": {"ucsc_bigwig" : "wigToBigWig"}}
    if outfile is None:
        outfile = "%s.bigwig" % os.path.splitext(bam_file)[0]
    if start > 0:
        start = int(start) - 1
    if end is not None:
        end = int(end)
    regions = [(chrom, start, end)]
    if not os.path.exists(outfile):
        wig_file = "%s.wig" % os.path.splitext(bam_file)[0]
        with open(wig_file, "w") as out_handle:
            chr_sizes, wig_valid = write_bam_track(bam_file, regions, config, out_handle)
        try:
            if wig_valid:
                convert_to_bigwig(wig_file, chr_sizes, config, outfile)
        finally:
            os.remove(wig_file)

@contextmanager
def indexed_bam(bam_file, config):
    if not os.path.exists(bam_file + ".bai"):
        pysam.index(bam_file)
    sam_reader = pysam.Samfile(bam_file, "rb")
    yield sam_reader
    sam_reader.close()

def write_bam_track(bam_file, regions, config, out_handle):
    out_handle.write("track %s\n" % " ".join(["type=wiggle_0",
        "name=%s" % os.path.splitext(os.path.split(bam_file)[-1])[0],
        "visibility=full",
        ]))
    sizes = []
    is_valid = False
    with indexed_bam(bam_file, config) as work_bam:
        for ref_info in work_bam.header.get("SQ", []):
            sizes.append((ref_info["SN"], ref_info["LN"]))
        if len(regions) == 1 and regions[0][0] == "all":
            regions = []
            for ref_info in work_bam.header.get("SQ", []):
                regions.append((ref_info["SN"], 0, None))
        for chrom, start, end in regions:
            if end is None:
                for ref_info in work_bam.header.get("SQ", []):
                    if ref_info["SN"] == chrom:
                        end = int(ref_info["LN"])
                        break
            assert end is not None, "Could not find %s in header" % chrom
            out_handle.write("variableStep chrom=%s\n" % chrom)
            for col in work_bam.pileup(chrom, start, end):
                out_handle.write("%s %s\n" % (col.pos+1, col.n))
                is_valid = True
    return sizes, is_valid

def convert_to_bigwig(wig_file, chr_sizes, config, bw_file=None):
    if not bw_file:
        bw_file = "%s.bigwig" % (os.path.splitext(wig_file)[0])
    size_file = "%s-sizes.txt" % (os.path.splitext(wig_file)[0])
    with open(size_file, "w") as out_handle:
        for chrom, size in chr_sizes:
            out_handle.write("%s\t%s\n" % (chrom, size))
    try:
        cl = [config["program"]["ucsc_bigwig"], wig_file, size_file, bw_file]
        subprocess.check_call(cl)
    finally:
        os.remove(size_file)
    return bw_file

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-o", "--outfile", dest="outfile")
    parser.add_option("-c", "--chrom", dest="chrom")
    parser.add_option("-s", "--start", dest="start")
    parser.add_option("-e", "--end", dest="end")
    (options, args) = parser.parse_args()
    if len(args) not in [1, 2]:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    kwargs = dict(
        outfile=options.outfile,
        chrom=options.chrom or 'all',
        start=options.start or 0,
        end=options.end)
    main(*args, **kwargs)
