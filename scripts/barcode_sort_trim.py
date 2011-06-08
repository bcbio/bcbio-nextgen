#!/usr/bin/env python
"""Identify fastq reads with barcodes, trimming and sorting for downstream work.

Given a fastq file or pair of fastq files containing barcodes, this identifies
the barcode in each read and writes it into a unique file. Mismatches are
allowed and barcode position within the reads can be specified.

Usage:
    barcode_sort_trim.py <barcode file> <out format> <in file> [<pair file>]
        --mismatch=n (number of allowed mismatches, default 1)
        --second (barcode is on the second read, defaults to first)
        --five (barcode is on the 5' end of the sequence, default to 3')
        --noindel (disallow insertion/deletions on barcode matches)
        --quiet (do not print out summary information on tags)

<barcode file> is a text file of:
    <name> <sequence>
for all barcodes present in the fastq multiplex.

<out format> specifies how the output files should be written:
    1_100721_FC626DUAAX_--b--_--r--_fastq.txt
  It should contain two values for substitution:
    --b-- Location of the barcode identifier
    --r-- Location of the read number (1 or 2)
  This can be used to specify any output location:
    /your/output/dir/out_--b--_--r--.txt

Requires:
    Python -- versions 2.6 or 2.7
    Biopython -- http://biopython.org
"""
from __future__ import with_statement
import sys
import os
import itertools
import unittest
import collections
import csv
from optparse import OptionParser

from Bio import pairwise2
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def main(barcode_file, out_format, in1, in2, mismatch,
         first_read, three_end, allow_indels,
         metrics_file, verbose):
    barcodes = read_barcodes(barcode_file)
    stats = collections.defaultdict(int)
    out_writer = output_to_fastq(out_format)
    for (name1, seq1, qual1), (name2, seq2, qual2) in itertools.izip(
            read_fastq(in1), read_fastq(in2)):
        end_gen = end_generator(seq1, seq2, first_read, three_end)
        bc_name, bc_seq, match_seq = best_match(end_gen, barcodes, mismatch,
                                                allow_indels)
        seq1, qual1, seq2, qual2 = remove_barcode(seq1, qual1, seq2, qual2,
                match_seq, first_read, three_end)
        out_writer(bc_name, name1, seq1, qual1, name2, seq2, qual2)
        stats[bc_name] += 1

    sort_bcs = []
    for bc in stats.keys():
        try:
            sort_bc = float(bc)
        except ValueError:
            sort_bc = str(bc)
        sort_bcs.append((sort_bc, bc))
    sort_bcs.sort()
    sort_bcs = [s[1] for s in sort_bcs]
    if verbose:
        print "% -10s %s" % ("barcode", "count")
        for bc in sort_bcs:
            print "% -10s %s" % (bc, stats[bc])
        print "% -10s %s" % ("total", sum(stats.values()))
    if metrics_file:
        with open(metrics_file, "w") as out_handle:
            writer = csv.writer(out_handle, dialect="excel-tab")
            for bc in sort_bcs:
                writer.writerow([bc, stats[bc]])

def best_match(end_gen, barcodes, mismatch, allow_indels=True):
    """Identify barcode best matching to the test sequence, with mismatch.

    Returns the barcode id, barcode sequence and match sequence.
    unmatched is returned for items which can't be matched to a barcode within
    the provided parameters.
    """
    if len(barcodes) == 1 and barcodes.values() == ["trim"]:
        size = len(barcodes.keys()[0])
        test_seq = end_gen(size)
        return barcodes.values()[0], test_seq, test_seq
    # easiest, fastest case -- exact match
    sizes = list(set(len(b) for b in barcodes.keys()))
    for s in sizes:
        test_seq = end_gen(s)
        try:
            bc_id = barcodes[test_seq]
            return bc_id, test_seq, test_seq
        except KeyError:
            pass
    # check for best approximate match within mismatch values
    match_info = []
    if mismatch > 0:
        for bc_seq, bc_id in barcodes.iteritems():
            test_seq = end_gen(len(bc_seq))
            aligns = pairwise2.align.globalms(bc_seq, test_seq,
                    5.0, -4.0, -9.0, -0.5, one_alignment_only=True)
            (abc_seq, atest_seq) = aligns[0][:2] if len(aligns) == 1 else ("", "")
            matches = sum(1 for i, base in enumerate(abc_seq)
                          if base == atest_seq[i])
            gaps = abc_seq.count("-")
            cur_mismatch = len(test_seq) - matches + gaps
            if cur_mismatch <= mismatch and (allow_indels or gaps == 0):
                match_info.append((cur_mismatch, bc_id, abc_seq, atest_seq))
    if len(match_info) > 0:
        match_info.sort()
        name, bc_seq, test_seq = match_info[0][1:]
        return name, bc_seq.replace("-", ""), test_seq.replace("-", "")
    else:
        return "unmatched", "", ""

def end_generator(seq1, seq2=None, first_read=True, three_end=True):
    """Function which pulls a barcode of a provided size from paired seqs.

    This respects the provided details about location of the barcode, returning
    items of the specified size to check against the read.
    """
    seq = seq1 if first_read else seq2
    assert seq is not None
    def _get_end(size):
        assert size > 0
        if three_end:
            return seq[-size:]
        else:
            return seq[:size]
    return _get_end

def _remove_from_end(seq, qual, match_seq, three_end):
    if match_seq:
        if three_end:
            assert seq[-len(match_seq):] == match_seq
            seq = seq[:-len(match_seq)]
            qual = qual[:-len(match_seq)]
        else:
            assert seq[:len(match_seq)] == match_seq
            seq = seq[len(match_seq):]
            qual = qual[len(match_seq):]
    return seq, qual

def remove_barcode(seq1, qual1, seq2, qual2, match_seq, first_read, three_end):
    """Trim found barcode from the appropriate sequence end.
    """
    if first_read:
        seq1, qual1 = _remove_from_end(seq1, qual1, match_seq, three_end)
    else:
        assert seq2 and qual2
        seq2, qual2 = _remove_from_end(seq2, qual2, match_seq, three_end)
    return seq1, qual1, seq2, qual2

def _write_to_handles(name, seq, qual, fname, out_handles):
    try:
        out_handle = out_handles[fname]
    except KeyError:
        out_handle = open(fname, "w")
        out_handles[fname] = out_handle
    out_handle.write("@%s\n%s\n+\n%s\n" % (name, seq, qual))

def output_to_fastq(output_base):
    """Write a set of paired end reads as fastq, managing output handles.
    """
    work_dir = os.path.dirname(output_base)
    if not os.path.exists(work_dir) and work_dir:
        try:
            os.makedirs(work_dir)
        except OSError:
            assert os.path.isdir(work_dir)
    out_handles = dict()
    def write_reads(barcode, name1, seq1, qual1, name2, seq2, qual2):
        read1name = output_base.replace("--r--", "1").replace("--b--", barcode)
        _write_to_handles(name1, seq1, qual1, read1name, out_handles)
        if name2:
            read2name = output_base.replace("--r--", "2").replace("--b--", barcode)
            _write_to_handles(name2, seq2, qual2, read2name, out_handles)
    return write_reads

def read_barcodes(fname):
    barcodes = {}
    with open(fname) as in_handle:
        for line in (l for l in in_handle if not l.startswith("#")):
            name, seq = line.rstrip("\r\n").split()
            barcodes[seq] = name
    return barcodes

def read_fastq(fname):
    """Provide read info from fastq file, potentially not existing.
    """
    if fname:
        with open(fname) as in_handle:
            for info in FastqGeneralIterator(in_handle):
                yield info
    else:
        for info in itertools.repeat((None, None, None)):
            yield info

# --- Testing code: run with 'nosetests -v -s barcode_sort_trim.py'

class BarcodeTest(unittest.TestCase):
    """Test identification and removal of barcodes with local alignments.
    """
    def setUp(self):
        self.barcodes = {"CGATGT": "2", "CAGATC": "7"}

    def test_1_end_generator(self):
        """Ensure the proper end is returned for sequences.
        """
        seq1, seq2 = ("AAATTT", "GGGCCC")
        end_gen = end_generator(seq1, seq2, True, True)
        assert end_gen(3) == "TTT"
        end_gen = end_generator(seq1, seq2, True, False)
        assert end_gen(3) == "AAA"
        assert end_gen(4) == "AAAT"
        end_gen = end_generator(seq1, seq2, False, True)
        assert end_gen(3) == "CCC"
        end_gen = end_generator(seq1, seq2, False, False)
        assert end_gen(3) == "GGG"

    def test_2_identical_match(self):
        """Ensure we can identify identical barcode matches.
        """
        bc_id, seq, _ = best_match(end_generator("CGATGT"), self.barcodes, 0)
        assert bc_id == "2"
        assert seq == "CGATGT"

    def test_3_allowed_mismatch(self):
        """Identify barcodes with the allowed number of mismatches.
        """
        # 1 and 2 mismatches
        (bc_id, _, _) = best_match(end_generator("CGTTGT"), self.barcodes, 1)
        assert bc_id == "2"
        (bc_id, _, _) = best_match(end_generator("GCATGT"), self.barcodes, 2)
        assert bc_id == "2"
        # single gap insertion
        (bc_id, _, _) = best_match(end_generator("GATTGT"), self.barcodes, 1)
        # single gap deletion
        (bc_id, _, _) = best_match(end_generator("GCGAGT"), self.barcodes, 1)
        assert bc_id == "unmatched"
        (bc_id, _, _) = best_match(end_generator("GCGAGT"), self.barcodes, 2)
        assert bc_id == "2"
        (bc_id, _, _) = best_match(end_generator("GCGAGT"), self.barcodes, 2, False)
        assert bc_id == "unmatched"
        # too many errors
        (bc_id, _, _) = best_match(end_generator("GCATGT"), self.barcodes, 1)
        assert bc_id == "unmatched"
        (bc_id, _, _) = best_match(end_generator("GCTTGT"), self.barcodes, 2)
        assert bc_id == "unmatched"

    def test_4_custom_barcodes(self):
        """ Detect longer non-standard custom barcodes, trimming
        """
        # Use the custom long barcode
        custom_barcode = dict((bc_seq, bc_id) for bc_id, bc_seq in self.barcodes.iteritems())
        # Simulate an arbitrary read, attach barcode and remove it from the 3' end
        seq = "GATTACA"*5+custom_barcode["7"]
        (bc_id, bc_seq, match_seq) = best_match(end_generator(seq), self.barcodes, 1)
        (removed, _, _, _) = remove_barcode(seq, "B"*9, seq, "g"*9, match_seq, True, True)
        # Was the barcode properly identified and removed with 1 mismatch allowed ?
        assert bc_id == "7"
        assert bc_seq == match_seq
        assert removed == "GATTACA"*5

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-s", "--second", dest="first_read",
                      action="store_false", default=True)
    parser.add_option("-f", "--five", dest="three_end",
                      action="store_false", default=True)
    parser.add_option("-i", "--noindel", dest="indels",
                      action="store_false", default=True)
    parser.add_option("-q", "--quiet", dest="verbose",
                      action="store_false", default=True)
    parser.add_option("-m", "--mismatch", dest="mismatch", default=1)
    parser.add_option("-o", "--metrics", dest="metrics_file", default=None)
    options, args = parser.parse_args()
    if len(args) == 3:
        barcode_file, out_format, in1 = args
        in2 = None
    elif len(args) == 4:
        barcode_file, out_format, in1, in2 = args
    else:
        print __doc__
        sys.exit()
    main(barcode_file, out_format, in1, in2, int(options.mismatch),
         options.first_read, options.three_end, options.indels,
         options.metrics_file, options.verbose)
