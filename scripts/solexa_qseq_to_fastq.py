#!/usr/bin/env python
"""Convert output solexa qseq files into fastq format, handling multiplexing.

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

Illumina barcoded samples contain barcodes in a separate qseq lane, which are
identified by being much shorter than the primary read. Barcodes are added to
the 3' end of the first sequence to remain consistent with other homebrew
barcoding methods.

Optional arguments:
    --failed (-f): Write out reads failing the Illumina quality checks instead.
    --outdir (-o): Write out fastq files to different output directory; defaults
                   to a directory named fastq in the current directory.
"""
from __future__ import with_statement
import os
import sys
import glob
from optparse import OptionParser

def main(run_name, lane_nums, do_fail=False, outdir=None):
    if outdir is None:
        outdir = os.path.join(os.getcwd(), "fastq")
    if not os.path.exists(outdir):
        try:
            os.makedirs(outdir)
        except OSError:
            assert os.path.isdir(outdir)
    if do_fail:
        fail_dir = os.path.join(outdir, "failed")
        if not os.path.exists(fail_dir):
            try:
                os.makedirs(fail_dir)
            except OSError:
                assert os.path.isdir(fail_dir)
    else:
        fail_dir = None
    for lane_num in lane_nums:
        lane_prefix = "s_%s" % lane_num
        out_prefix = "%s_%s" % (lane_num, run_name)
        write_lane(lane_prefix, out_prefix, outdir, fail_dir)

def write_lane(lane_prefix, out_prefix, outdir, fail_dir):
    qseq_files = glob.glob("%s_*qseq.txt" % lane_prefix)
    one_files, two_files, bc_files = _split_paired(qseq_files)
    is_paired = len(two_files) > 0
    out_files = (_get_outfiles(out_prefix, outdir, is_paired)
                 if not fail_dir else None)
    fail_files = (_get_outfiles(out_prefix, fail_dir, is_paired)
                  if fail_dir else None)
    for (num, files) in [("1", one_files), ("2", two_files)]:
        for i, fname in enumerate(files):
            bc_file = _get_associated_barcode(num, i, fname, bc_files)
            convert_qseq_to_fastq(fname, num, bc_file, out_files, fail_files)

def _get_associated_barcode(read_num, file_num, fname, bc_files):
    """Get barcodes for the first read if present.
    """
    if read_num == "1" and len(bc_files) > 0:
        bc_file = bc_files[file_num]
        bc_parts = bc_file.split("_")
        read_parts = fname.split("_")
        assert (bc_parts[1] == read_parts[1] and
                bc_parts[3] == read_parts[3]), (bc_parts, read_parts)
        return bc_file
    return None

def convert_qseq_to_fastq(fname, num, bc_file, out_files, fail_files=None):
    """Convert a qseq file into the appropriate fastq output.
    """
    bc_iterator = _qseq_iterator(bc_file, fail_files is None) if bc_file else None
    for basename, seq, qual, passed in _qseq_iterator(fname, fail_files is None):
        # if we have barcodes, add them to the 3' end of the sequence
        if bc_iterator:
            (_, bc_seq, bc_qual, _) = bc_iterator.next()
            seq = "%s%s" % (seq, bc_seq)
            qual = "%s%s" % (qual, bc_qual)
        name = "%s/%s" % (basename, num)
        out = "@%s\n%s\n+\n%s\n" % (name, seq, qual)
        if passed:
            out_files[num].write(out)
        elif fail_files:
            fail_files[num].write(out)

def _qseq_iterator(fname, pass_wanted):
    """Return the name, sequence, quality, and pass info of qseq reads.

    Names look like:

    HWI-EAS264:4:1:1111:3114#0/1
    """
    with open(fname) as qseq_handle:
        for line in qseq_handle:
            parts = line.strip().split("\t")
            passed = int(parts[-1]) == 1
            if passed is pass_wanted:
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
    """Identify first read, second read and barcode sequences in qseqs.

    Barcoded sequences are identified by being much shorter than reads
    in the first lane.
    """
    files.sort()
    one = []
    two = []
    bcs = []
    ref_size = None
    for f in files:
        parts = f.split("_")
        if parts[2] == "1":
            one.append(f)
            if ref_size is None:
                ref_size = _get_qseq_seq_size(f) // 2
        elif parts[2] == "2":
            cur_size = _get_qseq_seq_size(f)
            assert ref_size is not None
            if cur_size < ref_size:
                bcs.append(f)
            else:
                two.append(f)
        elif parts[2] == "3":
            two.append(f)
        else:
            raise ValueError("Unexpected part: %s" % f)
    one.sort()
    two.sort()
    bcs.sort()
    if len(two) > 0: assert len(two) == len(one)
    if len(bcs) > 0: assert len(bcs) == len(one)
    return one, two, bcs

def _get_qseq_seq_size(fname):
    with open(fname) as in_handle:
        parts = in_handle.readline().split("\t")
        return len(parts[8])

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-f", "--failed", dest="do_fail", action="store_true",
                      default=False)
    parser.add_option("-o", "--outdir", dest="outdir", action="store",
                      default=None)
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print __doc__
        sys.exit()
    main(args[0], args[1].split(","), options.do_fail, options.outdir)
