#!/usr/bin/env python
"""Small wrapper around Illumina pipeline scripts to build fastq files.

Usage:
    solexa_qseq_to_fastq.py <run name> <list of lane numbers>

The lane numbers should be separated by commas, so to build fastq files for all
lanes, you should pass:

    1,2,3,4,5,6,7,8
"""
from __future__ import with_statement
import os
import sys
import glob
import re
import subprocess
from optparse import OptionParser

from Mgh.Solexa import Config

def main(run_name, lane_nums):
    startdir = os.getcwd()
    outdir = os.path.join(startdir, "fastq")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for lane_num in lane_nums:
        lane_prefix = "s_%s" % lane_num
        out_prefix = "%s_%s" % (lane_num, run_name)
        write_lane(lane_prefix, out_prefix, outdir)

def write_lane(lane_prefix, out_prefix, outdir):
    qseq_files = glob.glob("%s_*qseq.txt" % lane_prefix)
    one_files, two_files = _split_paired(qseq_files)
    out_files = _get_outfiles(out_prefix, outdir, len(two_files) > 0)
    for (num, files) in [("1", one_files), ("2", two_files)]:
        for fname in files:
            convert_qseq_to_fastq(fname, num, out_files)

def convert_qseq_to_fastq(fname, num, out_files):
    """Convert a qseq file into the appropriate fastq output.
    """
    for basename, seq, qual in _qseq_iterator(fname):
        assert len(seq) == len(qual)
        split_seqs = [(num, seq, qual)]
        for snum, sseq, squal in split_seqs:
            name = "%s/%s" % (basename, snum)
            out_files[snum].write("@%s\n%s\n+%s\n%s\n" %
                    (name, sseq, name, squal))

def _qseq_iterator(fname):
    """Return the name, sequence, and quality of passing qseq reads.

    Names look like:

    HWI-EAS264:4:1:1111:3114#0/1
    """
    with open(fname) as qseq_handle:
        for line in qseq_handle:
            parts = line.strip().split("\t")
            if int(parts[-1]) == 1:
                name = ":".join([parts[0]] +  parts[2:6]) + "#" + parts[6]
                seq = parts[8].replace(".", "N")
                qual = parts[9]
                yield name, seq, qual
   
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

# --- Old code using Solexa tools, here for double checking verification
# --- purposes, if needed.

def old_main(run_lane, lane_nums):
    """The old main function, which uses the solexa pipeline.
    """
    config = Config.from_hostname()
    mask_script = os.path.join(config.pipeline_dir, "bin", "maskQseq.pl")
    build_script = os.path.join(config.pipeline_dir, "bin", "buildSeq.pl")
    
    startdir = os.getcwd()
    outdir = os.path.join(startdir, "fastq")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for lane_num in lane_nums:
        lane_prefix = "s_%s" % lane_num
        out_prefix = "%s_%s" % (lane_num, run_name)
        write_lane_old(lane_prefix, build_script, mask_script, out_prefix, outdir)
    # clean up intermediate files
    for to_remove in ["_ub_qseq.txt", "_tiles.txt"]:
        for rem_file in glob.glob(os.path.join(outdir, "*%s" % to_remove)):
            os.remove(rem_file)

def write_lane_old(lane_prefix, build_script, mask_script, out_prefix, outdir):
    tiles_file, tiles = make_tiles_file(lane_prefix, outdir)
    print "Writing mask files for", lane_prefix
    has_paired_end = make_mask_files(lane_prefix, tiles, outdir, mask_script)
    if has_paired_end:
        ends = [1, 2]
    else:
        ends = [1]
    start_dir = os.getcwd()
    for end in ends:
        try:
            os.chdir(outdir)
            cl = "%s --qualityFilter --fastq %s %s %s _ub_qseq.txt" % (build_script,
                    tiles_file, lane_prefix, end)
            print cl
            if len(ends) > 1:
                out_file = "%s_%s_fastq.txt" % (out_prefix, end)
            else:
                out_file = "%s_fastq.txt" % (out_prefix)
            with open(out_file, "w") as out_handle:
                subprocess.call(cl.split(), stdout=out_handle)
            print "Written to", out_file
        finally:
            os.chdir(start_dir)

def make_mask_files(lane_prefix, tiles, outdir, mask_script):
    """Generate masked files for all qseq files with the given prefix.
    """
    has_paired_end = False
    for tile in tiles:
        for paired_end in [1, 2]:
            base_name = "%s_%s_%04d" % (lane_prefix, paired_end, tile)
            in_file = "%s_qseq.txt" % (base_name)
            out_file = os.path.join(outdir, "%s_ub_qseq.txt" % base_name)
            if os.path.exists(in_file):
                if not has_paired_end and paired_end == 2:
                    has_paired_end = True
                if not os.path.exists(out_file):
                    with open(in_file) as in_handle:
                        line = in_handle.readline().split()
                        assert len(line[8]) == len(line[9])
                        read_length = len(line[8])
                        base_str = "Y" * read_length
                    cl = "%s --use_bases %s --orig_read_length %s %s --output %s" % (
                            mask_script, base_str, read_length, in_file, out_file)
                    print cl
                    subprocess.call(cl.split())
    return has_paired_end

def make_tiles_file(lane_prefix, outdir):
    match_pat = re.compile("%s_\d_(?P<tile>\d+)_.*" % lane_prefix)
    files = glob.glob("%s_*qseq*.txt" % lane_prefix)
    tiles = []
    for fname in files:
        match = match_pat.search(fname)
        tiles.append(int(match.group('tile')))
    tiles = list(set(tiles)) # remove duplicates for paired end reads
    tiles.sort()
    final_tiles = ["%s_%04d" % (lane_prefix, t) for t in tiles]
    out_file = os.path.join(outdir, "%s_tiles.txt" % lane_prefix)
    with open(out_file, "w") as out_handle:
        out_handle.write(" ".join(final_tiles) + "\n")
    return out_file, tiles

# ----

if __name__ == "__main__":
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print __doc__
        sys.exit()
    main(args[0], args[1].split(","))
