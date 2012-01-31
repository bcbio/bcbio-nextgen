#!/usr/bin/env python
"""Identify fastq reads with barcodes, trimming and sorting for downstream work.

Given a fastq file or pair of fastq files containing barcodes, this identifies
the barcode in each read and writes it into a unique file. Mismatches are
allowed and barcode position within the reads can be specified.

Usage:
    barcode_sort_trim.py <barcode file> <out format> <in file> [<pair file>]
        --mismatch=n (number of allowed mismatches, default 1)
        --bc_offset=n (an offset into the read where the barcode starts (5' barcode)
                       or ends (3' barcode))
        --read=n Integer read number containing the barcode (default to 1)
        --five (barcode is on the 5' end of the sequence, default to 3')
        --noindel (disallow insertion/deletions on barcode matches)
        --quiet (do not print out summary information on tags)
        --tag_title (append matched barcode to sequence header)

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

def main(barcode_file, out_format, in1, in2, in3, mismatch, bc_offset,
         bc_read_i, three_end, allow_indels,
         metrics_file, verbose, tag_title):
    barcodes = read_barcodes(barcode_file)
    stats = collections.defaultdict(int)
    out_writer = output_to_fastq(out_format)
    for (name1, seq1, qual1), (name2, seq2, qual2), (name3, seq3, qual3) in itertools.izip(
            read_fastq(in1), read_fastq(in2), read_fastq(in3)):
        end_gen = end_generator(seq1, seq2, seq3, bc_read_i, three_end, bc_offset)
        bc_name, bc_seq, match_seq = best_match(end_gen, barcodes, mismatch,
                                                allow_indels)
        seq1, qual1, seq2, qual2, seq3, qual3 = remove_barcode(
                seq1, qual1, seq2, qual2, seq3, qual3,
                match_seq, bc_read_i, three_end, bc_offset)
        if tag_title:
            name1 += " %s" % match_seq
            name2 += " %s" % match_seq
            name3 += " %s" % match_seq
        out_writer(bc_name, name1, seq1, qual1, name2, seq2, qual2,
                   name3, seq3, qual3)
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
    if mismatch > 0 or _barcode_has_ambiguous(barcodes):
        for bc_seq, bc_id in barcodes.iteritems():
            test_seq = end_gen(len(bc_seq))
            if _barcode_very_ambiguous(barcodes):
                gapopen_penalty = -18.0
            else:
                gapopen_penalty = -9.0
            aligns = pairwise2.align.globalms(bc_seq, test_seq,
                    5.0, -4.0, gapopen_penalty, -0.5, one_alignment_only=True)
            (abc_seq, atest_seq) = aligns[0][:2] if len(aligns) == 1 else ("", "")
            matches = sum(1 for i, base in enumerate(abc_seq)
                          if (base == atest_seq[i] or base == "N"))
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

def _barcode_very_ambiguous(barcodes):
    max_size = max(len(x) for x in barcodes.keys())
    max_ns = max(x.count("N") for x in barcodes.keys())
    return float(max_ns) / float(max_size) > 0.5

def _barcode_has_ambiguous(barcodes):
    for seq in barcodes.keys():
        if "N" in seq:
            return True
    return False

def end_generator(seq1, seq2=None, seq3=None, bc_read_i=1, three_end=True, bc_offset=0):
    """Function which pulls a barcode of a provided size from paired seqs.

    This respects the provided details about location of the barcode, returning
    items of the specified size to check against the read.
    """
    seq_choice = {1: seq1, 2: seq2, 3: seq3}
    seq = seq_choice[bc_read_i]
    assert seq is not None

    def _get_end(size):
        assert size > 0
        if three_end:
            return seq[-size-bc_offset:len(seq)-bc_offset]
        else:
            return seq[bc_offset:size+bc_offset]
    return _get_end

def _remove_from_end(seq, qual, match_seq, three_end, bc_offset):
    if match_seq:
        if three_end:
            assert seq[-len(match_seq)-bc_offset:len(seq)-bc_offset] == match_seq
            seq = seq[:-len(match_seq)-bc_offset]
            qual = qual[:-len(match_seq)-bc_offset]
        else:
            assert seq[bc_offset:len(match_seq)+bc_offset] == match_seq
            seq = seq[len(match_seq)+bc_offset:]
            qual = qual[len(match_seq)+bc_offset:]
    return seq, qual

def remove_barcode(seq1, qual1, seq2, qual2, seq3, qual3,
                   match_seq, bc_read_i, three_end, bc_offset=0):
    """Trim found barcode from the appropriate sequence end.
    """
    if bc_read_i == 1:
        seq1, qual1 = _remove_from_end(seq1, qual1, match_seq, three_end, bc_offset)
    elif bc_read_i == 2:
        assert seq2 and qual2
        seq2, qual2 = _remove_from_end(seq2, qual2, match_seq, three_end, bc_offset)
    else:
        assert bc_read_i == 3
        assert seq3 and qual3
        seq3, qual3 = _remove_from_end(seq3, qual3, match_seq, three_end, bc_offset)
    return seq1, qual1, seq2, qual2, seq3, qual3

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

    def write_reads(barcode, name1, seq1, qual1, name2, seq2, qual2,
                    name3, seq3, qual3):
        read1name = output_base.replace("--r--", "1").replace("--b--", barcode)
        _write_to_handles(name1, seq1, qual1, read1name, out_handles)
        if seq2:
            read2name = output_base.replace("--r--", "2").replace("--b--", barcode)
            _write_to_handles(name2, seq2, qual2, read2name, out_handles)
        if seq3:
            read3name = output_base.replace("--r--", "3").replace("--b--", barcode)
            _write_to_handles(name3, seq3, qual3, read3name, out_handles)
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
        for info in itertools.repeat(("", None, None)):
            yield info


# --- Testing code: run with 'nosetests -v -s barcode_sort_trim.py'

class BarcodeTest(unittest.TestCase):
    """Test identification and removal of barcodes with local alignments.
    """
    def setUp(self):
        self.barcodes = {"CGATGT": "2", "CAGATC": "7", "TTAGGCATC": "8"}

    def test_1_end_generator(self):
        """Ensure the proper end is returned for sequences.
        """
        seq1, seq2 = ("AAATTT", "GGGCCC")
        end_gen = end_generator(seq1, seq2, None, 1, True)
        assert end_gen(3) == "TTT"
        end_gen = end_generator(seq1, seq2, None, 1, False)
        assert end_gen(3) == "AAA"
        assert end_gen(4) == "AAAT"
        end_gen = end_generator(seq1, seq2, None, 2, True)
        assert end_gen(3) == "CCC"
        end_gen = end_generator(seq1, seq2, None, 2, False)
        assert end_gen(3) == "GGG"
        # Test end generation with an offset
        end_gen = end_generator(seq1, seq2, None, 1, True,1)
        assert end_gen(3) == "ATT"
        end_gen = end_generator(seq1, seq2, None, 1, False,1)
        assert end_gen(3) == "AAT"
        assert end_gen(4) == "AATT"
        end_gen = end_generator(seq1, seq2, None, 2, True,1)
        assert end_gen(3) == "GCC"
        end_gen = end_generator(seq1, seq2, None, 2, False,1)

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
        # with indels permitted, accepts 2 mismatches, even if "1" is specified
        (bc_id, _, _) = best_match(end_generator("CGAAGT"), self.barcodes, 1)
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
        seq = "GATTACA" * 5 + custom_barcode["8"]
        (bc_id, bc_seq, match_seq) = best_match(end_generator(seq), self.barcodes, 1)
        (removed, _, _, _, _, _) = remove_barcode(seq, "B" * 9, seq, "g" * 9, None, None,
                                                  match_seq, True, True)
        # Was the barcode properly identified and removed with 1 mismatch allowed ?
        assert bc_id == "8"
        assert bc_seq == match_seq
        assert removed == "GATTACA" * 5

    def test_5_illumina_barcodes(self):
        """ Test that Illumina reads with a trailing A are demultiplexed correctly
        """
        # Use the first barcode
        for bc_seq, bc_id in self.barcodes.items():
            if bc_id == "2":
                break
            
        # Simulate an arbitrary read, attach barcode and add a trailing A
        seq = "GATTACA" * 5 + bc_seq + "A"
        (bc_id, bc_seq, match_seq) = best_match(end_generator(seq,None,None,1,True,1), self.barcodes, 1)
        (removed, _, _, _, _, _) = remove_barcode(seq, "B" * 9, seq, "g" * 9, None, None,
                                                  match_seq, True, True, 1)
        # Was the barcode properly identified and removed with 1 mismatch allowed ?
        assert bc_id == "2"
        assert bc_seq == match_seq
        assert removed == "GATTACA" * 5

    def test_6_ambiguous_barcodes(self):
        """Allow mismatch N characters in specified barcodes.
        """
        bcs = {"CGATGN": "2", "CAGATC": "7"}
        (bc_id, _, _) = best_match(end_generator("CGATGT"), bcs, 0, False)
        assert bc_id == "2", bc_id
        (bc_id, _, _) = best_match(end_generator("CGATGN"), bcs, 0, False)
        assert bc_id == "2", bc_id
        (bc_id, _, _) = best_match(end_generator("CGATNT"), bcs, 0, False)
        assert bc_id == "unmatched", bc_id
        (bc_id, _, _) = best_match(end_generator("CGATNT"), bcs, 1, False)
        assert bc_id == "2", bc_id

    def test_7_very_ambiguous_barcodes(self):
        """Matching with highly ambiguous barcodes used for sorting."""
        bcs = {"ANNNNNN": "A", "CNNNNNN": "C", "GNNNNNN": "G", "TNNNNNN": "T"}
        (bc_id, _, _) = best_match(end_generator("CGGGAGA", bc_offset=0), bcs, 2, True)
        assert bc_id == "C", bc_id
        (bc_id, _, _) = best_match(end_generator("GCGGGAG", bc_offset=0), bcs, 0, False)
        assert bc_id == "G", bc_id

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-s", "--second", dest="deprecated_first_read",
                      action="store_false", default=True)
    parser.add_option("-r", "--read", dest="bc_read_i", action="store",
                      default=1)
    parser.add_option("-f", "--five", dest="three_end",
                      action="store_false", default=True)
    parser.add_option("-i", "--noindel", dest="indels",
                      action="store_false", default=True)
    parser.add_option("-q", "--quiet", dest="verbose",
                      action="store_false", default=True)
    parser.add_option("-m", "--mismatch", dest="mismatch", default=1)
    parser.add_option("-b", "--bc_offset", dest="bc_offset", default=0)
    parser.add_option("-o", "--metrics", dest="metrics_file", default=None)
    parser.add_option("-t", "--tag_title", dest="tag_title",
                      action="store_true", default=False)
    options, args = parser.parse_args()
    in2, in3 = (None, None)
    if len(args) == 3:
        barcode_file, out_format, in1 = args
    elif len(args) == 4:
        barcode_file, out_format, in1, in2 = args
    elif len(args) == 5:
        barcode_file, out_format, in1, in2, in3 = args
    else:
        print __doc__
        sys.exit()
    # handle deprecated less general options
    if options.deprecated_first_read is False:
        options.bc_read_i = 2
    main(barcode_file, out_format, in1, in2, in3, int(options.mismatch), int(options.bc_offset),
         int(options.bc_read_i), options.three_end, options.indels,
         options.metrics_file, options.verbose, options.tag_title)
