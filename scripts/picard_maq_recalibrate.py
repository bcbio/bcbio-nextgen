#!/usr/bin/env python
"""Use Picard tools to run Maq with recalibrated Fastq quality scores.

This starts with input fastq files, either single or paired end, and generates
an output alignment file to feed into downstream analyses.

Usage:
    picard_maq_recalibrate.py <config file> <sample name> <reference file>
                              <fastq 1> [<fastq 2>]

Requires:
    Picard
    maq
"""
import os
import sys
import glob
import subprocess
from optparse import OptionParser

import yaml

from bcbio.picard import PicardRunner
from bcbio.utils import curdir_tmpdir

def main(config_file, out_base, ref_file, read1, read2=None, sample_name=""):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    barcode, lane = _get_sample_details(read1)
    maq_cmd = config["program"]["maq"]
    stringency = config["algorithm"]["stringency"]

    picard = PicardRunner(config["program"]["picard"])
    bam_reads = fastq_to_bam(picard, sample_name,
            config["algorithm"]["quality_format"], read1, read2)
    base_align = picard_run_maq(picard, maq_cmd, bam_reads, ref_file, barcode,
            lane, out_base, stringency, read2 is not None, limit=1e6)
    cal_bam_reads = calibrate_scores(picard, bam_reads, base_align, ref_file)
    final_align = picard_run_maq(picard, maq_cmd, cal_bam_reads, ref_file,
            barcode, lane, out_base, stringency, read2 is not None, ext="-cal")

def _get_sample_details(in_file):
    read_parts = os.path.split(in_file)[-1].split("_")
    lane = str(int(read_parts[0]))
    barcode = [p for p in read_parts if p.endswith("XX")][0]
    print "Running alignment for %s, lane %s" % (barcode, lane)
    return barcode, lane

def calibrate_scores(picard, input_bam, base_align, ref_file):
    base, ext = os.path.splitext(input_bam)
    cal_bam = "%s-cal.bam" % base
    table_file = "%s-cal.table" % base
    if not os.path.exists(cal_bam):
        opts = [("ALIGNED_SAM", base_align),
                ("INPUT", input_bam),
                ("OUTPUT", cal_bam),
                ("REFERENCE_SEQUENCE", ref_file),
                ("TABLE", table_file)
                ]
        picard.run("CalibrateQualityScores", opts)
    return cal_bam

def picard_run_maq(picard, maq_cmd, input_bam, ref_file, barcode, lane,
        out_base, stringency, is_paired=False, limit=None, ext=""):
    out_dir = "%s-maq%s" % (out_base, ext)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    bam_out_file = "%s.bam" % (out_dir)
    with curdir_tmpdir() as tmp_dir:
        std_opts = [("INPUT", input_bam),
                    ("ANALYSIS_DIR", out_dir),
                    ("FLOWCELL_BARCODE", barcode),
                    ("LANE", lane),
                    ("REFERENCE_SEQUENCE", ref_file),
                    ("TMP_DIR", tmp_dir),
                    ("PAIRED_RUN", ("true" if is_paired else "false"))]
        # Convert fastq to Maq ready files
        if len(glob.glob(
               os.path.join(out_dir, "%s.%s*bfq" % (barcode, lane)))) == 0:
            opts = std_opts + [
                    ("PREPARE", "true"),
                    ]
            if limit:
                opts.append(("READS_TO_ALIGN", int(limit)))
            picard.run("RunMaq", opts)
        # actually run Maq. Use python as Picard is failing with same parameters
        if len(glob.glob(
               os.path.join(out_dir, "%s.%s*out*.map" % (barcode, lane)))) == 0:
            #opts = std_opts + [
            #        ("STRINGENCY", stringency),
            #        ("ALIGN", "true"),
            #        ]
            #picard.run("RunMaq", opts)
            run_maq(maq_cmd, stringency, out_dir, ref_file, barcode, lane)
        # Convert the output file to BAM aligned
        if not os.path.exists(bam_out_file):
            opts = std_opts + [
                    ("OUTPUT", bam_out_file),
                    ("BAM_OUTPUT", "true")
                    ]
            index_file = index_ref_file(picard, ref_file)
            picard.run("RunMaq", opts)
            #convert_map_to_bam(picard, out_dir, bam_out_file, ref_file,
            #        barcode, lane)
    return bam_out_file

def convert_map_to_bam(picard, out_dir, bam_out_file, ref_file, barcode, lane):
    """Covert Maq map files into an output BAM file.
    """
    opts = [("OUTPUT", bam_out_file),
            ("SEQUENCE_DICTIONARY", index_file)]
    for map_file in map_files(out_dir, barcode, lane):
        opts.append(("MAP_FILE", map_file))
    picard.run("MaqToSam", opts)

def convert_map_to_bam_old(picard, out_dir, bam_out_file, ref_file, barcode, lane):
    """Convert an output file from Maq to BAM format.
    """
    sam_file = None
    for i, map_file in enumerate(map_files(out_dir, barcode, lane)):
        cl = ["maq2sam-long", map_file]
        print " ".join(cl)
        if sam_file is None:
            sam_file = "%s.sam" % os.path.splitext(map_file)[0]
        with open(sam_file, ("w" if i == 0 else "a")) as out_handle:
            child = subprocess.Popen(cl, stdout=out_handle)
            child.wait()
    assert sam_file is not None
    # XXX Picard does not write the reference header
    #opts = [("INPUT", sam_file), ("OUTPUT", bam_out_file)]
    #picard.run("SamFormatConverter", opts)
    ref_dir, ref_base = os.path.split(ref_file)
    sam_ref_file = os.path.join(ref_dir, os.pardir, "seq", "%s.fa.fai" %
            os.path.splitext(ref_base)[0])
    cl = ["samtools", "view", "-bt", sam_ref_file, sam_file]
    with open(bam_out_file, "w") as out_handle:
        print " ".join(cl)
        child = subprocess.Popen(cl, stdout=out_handle)
        child.wait()

def map_files(out_dir, barcode, lane):
    map_files = glob.glob(os.path.join(out_dir, 
                                       "%s.%s*.map" % (barcode, lane)))
    map_files.sort()
    for m in map_files:
        yield m

def run_maq(maq_cmd, stringency, out_dir, ref_file, barcode, lane):
    """Run Maq directly as would be done by Picard.
    """
    # XXX Don't exactly correlate to Picard params
    s_params = dict(
            high = ["-a", "500", "-e", "100"],
            low = ["-a", "1500", "-e", "2100"])
    for input_num, input_files in bfq_files(out_dir, barcode, lane):
        out_file = os.path.join(out_dir, "%s.%s.%s.out.aln.map" % (barcode,
            lane, input_num))
        if not os.path.exists(out_file):
            cl = ["maq", "map", "-s", "0"] + s_params[stringency] + [
                    out_file,
                    "%s.bfa" % os.path.splitext(ref_file)[0]
                    ] + input_files
            err_file = os.path.join(out_dir, "%s.%s.%s.out.map.log" % (barcode,
                lane, input_num))
            print " ".join(cl)
            with open(err_file, "w") as out_handle:
                child = subprocess.Popen(cl, stderr=out_handle)
                child.wait()

def bfq_files(out_dir, barcode, lane):
    """Generator supplying bfq files in a directory, split by number.
    """
    input_files = glob.glob(os.path.join(out_dir,
                                         "%s.%s*.bfq" % (barcode, lane)))
    input_files.sort()
    cur_num = 0
    cur_files = []
    for in_file in input_files:
        (_, _, num, _, _) = os.path.split(in_file)[-1].split(".")
        if num == cur_num:
            cur_files.append(in_file)
        else:
            if len(cur_files) > 0:
                yield cur_num, cur_files
            cur_num = num
            cur_files = [in_file]
    if len(cur_files) > 0:
        yield cur_num, cur_files

def fastq_to_bam(picard, sample_name, quality_format, read1, read2):
    base, ext = os.path.splitext(os.path.basename(read1))
    out_file = "%s.bam" % base
    if not os.path.exists(out_file):
        with curdir_tmpdir() as tmp_dir:
            opts = [("FASTQ", read1),
                    ("TMP_DIR", tmp_dir),
                    ("QUALITY_FORMAT", quality_format),
                    ("SAMPLE_NAME", sample_name),
                    ("OUTPUT", out_file)]
            if read2:
                opts.append(("FASTQ2", read2))
            picard.run("FastqToSam", opts)
    return out_file

def index_ref_file(picard, ref_file):
    """Provide a Picard style dict index file for a reference genome.
    """
    dict_file = "%s.dict" % os.path.splitext(ref_file)[0]
    if not os.path.exists(dict_file):
        opts = [("REFERENCE", ref_file),
                ("OUTPUT", dict_file)]
        picard.run("CreateSequenceDictionary", opts)
    return dict_file

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-n", "--name", dest="sample_name")
    (options, args) = parser.parse_args()
    kwargs = dict(
            sample_name=options.sample_name)
    main(*args, **kwargs)
