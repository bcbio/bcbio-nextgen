#!/usr/bin/env python
"""Convert a SAM alignment file into a sorted BAM file with all reads.

This creates an output BAM file with both aligned and unaligned reads
in a format ready for downstream analyses.

Usage:
    picard_sam_to_bam.py <config_file> <alignment sam file>
                         <reference file>
                         <fastq_reads> [<fastq pair>]
"""
import sys
import os
import tempfile
from optparse import OptionParser

import yaml

from bcbio.broad import BroadRunner
from bcbio.utils import curdir_tmpdir, save_diskspace

def main(config_file, align_sam, ref_file, fastq_one, fastq_pair=None,
        sample_name="", rg_name="", pu_name=""):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    picard = BroadRunner(config["program"]["picard"],
                         max_memory=config["algorithm"].get("java_memory", ""))
    platform = config["algorithm"]["platform"]
    picard.run_fn("picard_index_ref", ref_file)
    base_dir = os.path.dirname(align_sam)
    out_fastq_bam = picard.run_fn("picard_fastq_to_bam", fastq_one, fastq_pair,
                                  base_dir, platform, sample_name, rg_name, pu_name)
    out_bam = picard.run_fn("picard_sam_to_bam", align_sam, out_fastq_bam, ref_file,
                            fastq_pair is not None)
    sort_bam = picard.run_fn("picard_sort", out_bam)
    save_diskspace(out_fastq_bam, "Combined into output BAM %s" % out_bam, config)
    save_diskspace(out_bam, "Sorted to %s" % sort_bam, config)

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-n", "--name", dest="sample_name", default="")
    parser.add_option("-r", "--rg", dest="rg_name", default="A")
    parser.add_option("-p", "--pu", dest="pu_name", default="")
    (options, args) = parser.parse_args()
    kwargs = dict(
            sample_name=options.sample_name,
            rg_name=options.rg_name,
            pu_name=options.pu_name)
    main(*args, **kwargs)
