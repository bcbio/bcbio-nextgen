#!/usr/bin/env python
"""Convert BAM files to BigWig file format in a specified region.

Usage:
    bam_to_wiggle.py <BAM file> [<YAML config>]
    [--outfile=<output file name>
     --chrom=<chrom>
     --start=<start>
     --end=<end>
     --normalize]

chrom start and end are optional, in which case they default to everything.
The normalize flag adjusts counts to reads per million.

The config file is in YAML format and specifies the location of the wigToBigWig
program from UCSC:

program:
  ucsc_bigwig: wigToBigWig

If not specified, these will be assumed to be present in the system path.

The script requires:
    pysam (http://code.google.com/p/pysam/)
    wigToBigWig from UCSC (http://hgdownload.cse.ucsc.edu/admin/exe/)
If a configuration file is used, then PyYAML is also required (http://pyyaml.org/)
along with bcbio
"""
import os
import sys
import subprocess
import tempfile
from optparse import OptionParser
from contextlib import contextmanager, closing

import pysam

def main(bam_file, config_file=None, chrom='all', start=0, end=None,
         outfile=None, normalize=False, use_tempfile=False):
    if config_file:
        from bcbio.pipeline.config_utils import load_config, get_program
        config = load_config(config_file)
    else:
        config = {"program": {"ucsc_bigwig": "wigToBigWig"}}
        get_program = None
    if outfile is None:
        outfile = "%s.bigwig" % os.path.splitext(bam_file)[0]
    if start > 0:
        start = int(start) - 1
    if end is not None:
        end = int(end)
    regions = [(chrom, start, end)]
    if os.path.abspath(bam_file) == os.path.abspath(outfile):
        sys.stderr.write("Bad arguments, input and output files are the same.\n")
        sys.exit(1)
    if not (os.path.exists(outfile) and os.path.getsize(outfile) > 0):
        if use_tempfile:
            #Use a temp file to avoid any possiblity of not having write permission
            out_handle = tempfile.NamedTemporaryFile(delete=False)
            wig_file = out_handle.name
        else:
            wig_file = "%s.wig" % os.path.splitext(outfile)[0]
            out_handle = open(wig_file, "w")
        with closing(out_handle):
            chr_sizes, wig_valid = write_bam_track(bam_file, regions, config, out_handle,
                                                   normalize)
        try:
            if wig_valid:
                convert_to_bigwig(wig_file, chr_sizes, config, outfile, get_program=get_program)
        finally:
            os.remove(wig_file)

@contextmanager
def indexed_bam(bam_file, config):
    if not os.path.exists(bam_file + ".bai"):
        pysam.index(bam_file)
    sam_reader = pysam.Samfile(bam_file, "rb")
    yield sam_reader
    sam_reader.close()

def write_bam_track(bam_file, regions, config, out_handle, normalize):
    out_handle.write("track %s\n" % " ".join(["type=wiggle_0",
        "name=%s" % os.path.splitext(os.path.split(bam_file)[-1])[0],
        "visibility=full",
        ]))
    normal_scale = 1e6
    is_valid = False
    with indexed_bam(bam_file, config) as work_bam:
        total = sum(1 for r in work_bam.fetch() if not r.is_unmapped) if normalize else None
        sizes = zip(work_bam.references, work_bam.lengths)
        if len(regions) == 1 and regions[0][0] == "all":
            regions = [(name, 0, length) for name, length in sizes]
        for chrom, start, end in regions:
            if end is None and chrom in work_bam.references:
                end = work_bam.lengths[work_bam.references.index(chrom)]
            assert end is not None, "Could not find %s in header" % chrom
            out_handle.write("variableStep chrom=%s\n" % chrom)
            for col in work_bam.pileup(chrom, start, end):
                if normalize:
                    n = float(col.n) / total * normal_scale
                else:
                    n = col.n
                out_handle.write("%s %.1f\n" % (col.pos+1, n))
                is_valid = True
    return sizes, is_valid

def convert_to_bigwig(wig_file, chr_sizes, config, bw_file=None, get_program=None):
    if not bw_file:
        bw_file = "%s.bigwig" % (os.path.splitext(wig_file)[0])
    size_file = "%s-sizes.txt" % (os.path.splitext(wig_file)[0])
    with open(size_file, "w") as out_handle:
        for chrom, size in chr_sizes:
            out_handle.write("%s\t%s\n" % (chrom, size))
    cmd = get_program("ucsc_bigwig", config, default="wigToBigWig") if get_program else "wigToBigWig"
    try:
        cl = [cmd, wig_file, size_file, bw_file]
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
    parser.add_option("-n", "--normalize", dest="normalize",
                      action="store_true", default=False)
    parser.add_option("-t", "--tempfile", dest="use_tempfile",
                      action="store_true", default=False)
    (options, args) = parser.parse_args()
    if len(args) not in [1, 2]:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    kwargs = dict(
        outfile=options.outfile,
        chrom=options.chrom or 'all',
        start=options.start or 0,
        end=options.end,
        normalize=options.normalize,
        use_tempfile=options.use_tempfile)
    main(*args, **kwargs)
