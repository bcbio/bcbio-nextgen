#!/usr/bin/env python
"""Prepare paired end fastq files from a chromosome region in an aligned input BAM file.

Useful for preparing test or other example files with subsets of aligned data.

Usage:
  bam_to_fastq_region.py <YAML config> <BAM input> <chromosome> <start> <end>
"""
import os
import sys
import contextlib

import yaml
import pysam
from Bio import Seq

from bcbio import broad

def main(config_file, in_file, space, start, end):
    with open(config_file) as in_handle:
        config = yaml.safe_load(in_handle)
    runner = broad.runner_from_config(config)
    target_region = (space, int(start), int(end))
    for pair in [1, 2]:
        out_file = "%s_%s-%s.fastq" % (os.path.splitext(os.path.basename(in_file))[0],
                                          pair, target_region[0])
        with open(out_file, "w") as out_handle:
            for name, seq, qual in bam_to_fastq_pair(in_file, target_region, pair):
                out_handle.write("@%s/%s\n%s\n+\n%s\n" % (name, pair, seq, qual))
        sort_fastq(out_file, runner)

def bam_to_fastq_pair(in_file, target_region, pair):
    """Generator to convert BAM files into name, seq, qual in a region.
    """
    space, start, end = target_region
    bam_file = pysam.Samfile(in_file, "rb")
    for read in bam_file:
        if (not read.is_unmapped and not read.mate_is_unmapped
                and bam_file.getrname(read.tid) == space
                and bam_file.getrname(read.mrnm) == space
                and read.pos >= start and read.pos <= end
                and read.mpos >= start and read.mpos <= end
                and not read.is_secondary
                and read.is_paired and getattr(read, "is_read%s" % pair)):
            seq = Seq.Seq(read.seq)
            qual = list(read.qual)
            if read.is_reverse:
                seq = seq.reverse_complement()
                qual.reverse()
            yield read.qname, str(seq), "".join(qual)

@contextlib.contextmanager
def fastq_to_bam(in_file, runner):
    bam_file = "%s.bam" % os.path.splitext(in_file)[0]
    try:
        opts = [("FASTQ", in_file),
                ("OUTPUT", bam_file),
                ("QUALITY_FORMAT", "Standard"),
                ("SAMPLE_NAME", "t")]
        runner.run("FastqToSam", opts)
        yield bam_file
    finally:
        if os.path.exists(bam_file):
            os.remove(bam_file)

@contextlib.contextmanager
def sort_bam(in_file, runner):
    base, ext = os.path.splitext(in_file)
    out_file = "%s-sort%s" % (base, ext)
    try:
        opts = [("INPUT", in_file),
                ("OUTPUT", out_file),
                ("SORT_ORDER", "queryname")]
        runner.run("SortSam", opts)
        yield out_file
    finally:
        if os.path.exists(out_file):
            os.remove(out_file)

def bam_to_fastq(in_file, runner):
    out_file = "%s.fastq" % os.path.splitext(in_file)[0]
    opts = [("INPUT", in_file),
            ("FASTQ", out_file)]
    runner.run("SamToFastq", opts)

def sort_fastq(in_file, runner):
    with fastq_to_bam(in_file, runner) as bam_file:
        with sort_bam(bam_file, runner) as sort_bam_file:
            bam_to_fastq(sort_bam_file, runner)

if __name__ == "__main__":
    main(*sys.argv[1:])

