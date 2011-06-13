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
    if platform.lower() == "illumina":
        qual_format = "Illumina"
    else:
        raise ValueError("Need to specify quality format for %s" % platform)
    picard.run_fn("picard_index_ref", ref_file)
    base_dir = os.path.split(align_sam)[0]
    with curdir_tmpdir() as tmp_dir:
        out_fastq_bam = picard_fastq_to_bam(picard, fastq_one, fastq_pair,
                base_dir, platform, qual_format, sample_name, rg_name, pu_name,
                tmp_dir)
        out_bam = picard_merge_bam(picard, align_sam, out_fastq_bam,
                ref_file, tmp_dir, fastq_pair is not None)
        sort_bam = picard.run_fn("picard_sort", out_bam)
    save_diskspace(out_fastq_bam, "Combined into output BAM %s" % out_bam, config)
    save_diskspace(out_bam, "Sorted to %s" % sort_bam, config)

def picard_merge_bam(picard, align_sam, fastq_bam, ref_file,
        tmp_dir, is_paired):
    out_bam = "%s.bam" % os.path.splitext(align_sam)[0]
    if not os.path.exists(out_bam):
        opts = [("UNMAPPED", fastq_bam),
                ("ALIGNED", align_sam),
                ("OUTPUT", out_bam),
                ("REFERENCE_SEQUENCE", ref_file),
                ("TMP_DIR", tmp_dir),
                ("PAIRED_RUN", ("true" if is_paired else "false")),
                ]
        picard.run("MergeBamAlignment", opts)
    return out_bam

def picard_fastq_to_bam(picard, fastq_one, fastq_two, base_dir,
        platform, qual_format, sample_name, rg_name, pu_name, tmp_dir):
    out_bam = os.path.join(base_dir, "%s.bam" %
            os.path.splitext(os.path.split(fastq_one)[1])[0])
    if not (os.path.exists(out_bam) and os.path.getsize(out_bam) > 0):
        opts = [("FASTQ", fastq_one),
                ("QUALITY_FORMAT", qual_format),
                ("READ_GROUP_NAME", rg_name),
                ("SAMPLE_NAME", sample_name),
                ("PLATFORM_UNIT", pu_name),
                ("PLATFORM", platform),
                ("TMP_DIR", tmp_dir),
                ("OUTPUT", out_bam),
                ]
        if fastq_two:
            opts.append(("FASTQ2", fastq_two))
        picard.run("FastqToSam", opts)
    return out_bam

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
