"""Utilities for working with fastq files.
"""

import difflib
from itertools import izip

from Bio import SeqIO

from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio import utils


@utils.memoize_outfile(stem=".groom")
def groom(in_file, in_qual="fastq-sanger", out_dir=None, out_file=None):
    """
    Grooms a FASTQ file into sanger format, if it is not already in that
    format. Use fastq-illumina for Illumina 1.3-1.7 qualities and
    fastq-solexa for the original solexa qualities. When in doubt, your
    sequences are probably fastq-sanger.

    """
    if in_qual == "fastq-sanger":
        logger.info("%s is already in Sanger format." % (in_file))
        return out_file
    with file_transaction(out_file) as tmp_out_file:
        count = SeqIO.convert(in_file, in_qual, tmp_out_file, "fastq-sanger")
    logger.info("Converted %d reads in %s to %s." % (count, in_file, out_file))
    return out_file

@utils.memoize_outfile(stem=".fixed")
def filter_single_reads_by_length(in_file, quality_format, min_length=20,
                                  out_file=None):
    """
    removes reads from a fastq file which are shorter than a minimum
    length

    """
    logger.info("Removing reads in %s thare are less than %d bases."
                % (in_file, min_length))
    in_iterator = SeqIO.parse(in_file, quality_format)
    out_iterator = (record for record in in_iterator if
                    len(record.seq) > min_length)
    with file_transaction(out_file) as tmp_out_file:
        with open(tmp_out_file, "w") as out_handle:
            SeqIO.write(out_iterator, out_handle, quality_format)
    return out_file

def filter_reads_by_length(fq1, fq2, quality_format, min_length=20):
    """
    removes reads from a pair of fastq files that are shorter than
    a minimum length. removes both ends of a read if one end falls
    below the threshold while maintaining the order of the reads

    """

    logger.info("Removing reads in %s and %s that "
                "are less than %d bases." % (fq1, fq2, min_length))
    fq1_out = utils.append_stem(fq1, ".fixed")
    fq2_out = utils.append_stem(fq2, ".fixed")
    fq1_single = utils.append_stem(fq1, ".singles")
    fq2_single = utils.append_stem(fq2, ".singles")
    if all(map(utils.file_exists, [fq1_out, fq2_out, fq2_single, fq2_single])):
        return [fq1_out, fq2_out]

    fq1_in = SeqIO.parse(fq1, quality_format)
    fq2_in = SeqIO.parse(fq2, quality_format)

    with open(fq1_out, 'w') as fq1_out_handle, open(fq2_out, 'w') as fq2_out_handle, open(fq1_single, 'w') as fq1_single_handle, open(fq2_single, 'w') as fq2_single_handle:
        for fq1_record, fq2_record in izip(fq1_in, fq2_in):
            if len(fq1_record.seq) >= min_length and len(fq2_record.seq) >= min_length:
                fq1_out_handle.write(fq1_record.format(quality_format))
                fq2_out_handle.write(fq2_record.format(quality_format))
            else:
                if len(fq1_record.seq) > min_length:
                    fq1_single_handle.write(fq1_record.format(quality_format))
                if len(fq2_record.seq) > min_length:
                    fq2_single_handle.write(fq2_record.format(quality_format))

    return [fq1_out, fq2_out]

def combine_pairs(input_files):
    """ calls files pairs if they are completely the same except
    for one has _1 and the other has _2 returns a list of tuples
    of pairs or singles.
    From bipy.utils (https://github.com/roryk/bipy/blob/master/bipy/utils.py)"""
    PAIR_FILE_IDENTIFIERS = ["1", "2"]

    pairs = []
    used = []
    for in_file in input_files:
        if in_file in used:
            continue
        for comp_file in input_files:
            if comp_file in used:
                continue
            s = difflib.SequenceMatcher(a=in_file, b=comp_file)
            blocks = s.get_matching_blocks()
            # length 3 means on match in the middle of the string
            if len(s.get_matching_blocks()) is not 3:
                continue
            if comp_file[blocks[0][2]] in PAIR_FILE_IDENTIFIERS:
                # e.g. _R1, _R2 or _1, _2
                if comp_file[blocks[0][2] - 1] in ("R", "_"):
                    used.append(in_file)
                    used.append(comp_file)
                    pairs.append([in_file, comp_file])
                    break
        if in_file not in used:
            pairs.append([in_file])
            used.append(in_file)

    return pairs
