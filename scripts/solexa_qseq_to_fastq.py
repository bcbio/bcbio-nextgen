#!/usr/bin/env python
"""Convert output solexa qseq files into fastq format.

Works with qseq output from Illumina's on-machine base caller in:

Data/Intensities/BaseCalls/

or from the offline base caller in:

Data/*_Firecrest*/Bustard*

Usage:
    solexa_qseq_to_fastq.py <run name> <list of lane numbers>

The lane numbers should be separated by commas, so to build fastq files for all
lanes, you should pass:

    1,2,3,4,5,6,7,8

Output files will be in the fastq directory as <lane>_<run_name>_fastq.txt

Optional arguments:
    --failed (-f): Also write out reads failing the Illumina quality checks.
"""
from __future__ import with_statement
import os
import sys
import glob
from optparse import OptionParser

def main(run_name, lane_nums, do_fail=False):
    startdir = os.getcwd()
    outdir = os.path.join(startdir, "fastq")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if do_fail:
        fail_dir = os.path.join(outdir, "failed")
        if not os.path.exists(fail_dir):
            os.makedirs(fail_dir)
    else:
        fail_dir = None
    for lane_num in lane_nums:
        lane_prefix = "s_%s" % lane_num
        out_prefix = "%s_%s" % (lane_num, run_name)
        write_lane(lane_prefix, out_prefix, outdir, fail_dir)

def write_lane(lane_prefix, out_prefix, outdir, fail_dir):
    qseq_files = glob.glob("%s_*qseq.txt" % lane_prefix)
    one_files, two_files = _split_paired(qseq_files)
    is_paired = len(two_files) > 0
    out_files = _get_outfiles(out_prefix, outdir, is_paired)
    fail_files = (_get_outfiles(out_prefix, fail_dir, is_paired)
                  if fail_dir else None)
    for (num, files) in [("1", one_files), ("2", two_files)]:
        for fname in files:
            convert_qseq_to_fastq(fname, num, out_files, fail_files)

def convert_qseq_to_fastq(fname, num, out_files, fail_files=None):
    """Convert a qseq file into the appropriate fastq output.
    """
    for basename, seq, qual, passed in _qseq_iterator(fname):
        name = "%s/%s" % (basename, num)
        out = "@%s\n%s\n+%s\n%s\n" % (name, seq, name, qual)
        if passed:
            out_files[num].write(out)
        elif fail_files:
            fail_files[num].write(out)

def _qseq_iterator(fname):
    """Return the name, sequence, quality, and pass info of qseq reads.

    Names look like:

    HWI-EAS264:4:1:1111:3114#0/1
    """
    with open(fname) as qseq_handle:
        for line in qseq_handle:
            parts = line.strip().split("\t")
            passed = int(parts[-1]) == 1
            name = ":".join([parts[0]] +  parts[2:6]) + "#" + parts[6]
            seq = parts[8].replace(".", "N")
            qual = parts[9]
            assert len(seq) == len(qual)
            yield name, seq, qual, passed

def _get_outfiles(out_prefix, outdir, has_paired_files):
    out_files = {}
    if has_paired_files:
        for num in ("1", "2"):
            out_files[num] = os.path.join(outdir, "%s_%s_fastq.txt" % (
                out_prefix, num))
    else:
        out_files["1"] = os.path.join(outdir, "%s_fastq.txt" % out_prefix)
    for index, fname in out_files.items():
        out_files[index] = open(fname, "w")
    return out_files

def _split_paired(files):
    one = []
    two = []
    for f in files:
        parts = f.split("_")
        if parts[2] == "1":
            one.append(f)
        else:
            two.append(f)
    one.sort()
    two.sort()
    return one, two

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-f", "--failed", dest="do_fail", action="store_true",
            default=False)
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print __doc__
        sys.exit()
    main(args[0], args[1].split(","), options.do_fail)
