"""Utilities for working with fastq files.
"""

from itertools import izip, product
import os
import random
import gzip

from Bio import SeqIO
from bcbio.distributed import objectstore
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio import utils
from bcbio.utils import open_possible_gzip
from bcbio.pipeline import config_utils
from bcbio.provenance import do

@utils.memoize_outfile(stem=".groom")
def groom(in_file, data, in_qual="illumina", out_dir=None, out_file=None):
    """
    Grooms a FASTQ file from Illumina 1.3/1.5 quality scores into
    sanger format, if it is not already in that format.
    """
    seqtk = config_utils.get_program("seqtk", data["config"])
    if in_qual == "fastq-sanger":
        logger.info("%s is already in Sanger format." % in_file)
        return out_file
    with file_transaction(out_file) as tmp_out_file:
        cmd = "{seqtk} seq -Q64 {in_file} | gzip > {tmp_out_file}".format(**locals())
        do.run(cmd, "Converting %s to Sanger format." % in_file)
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

    out_files = [fq1_out, fq2_out, fq1_single, fq2_single]

    with file_transaction(out_files) as tmp_out_files:
        fq1_out_handle = open(tmp_out_files[0], "w")
        fq2_out_handle = open(tmp_out_files[1], "w")
        fq1_single_handle = open(tmp_out_files[2], "w")
        fq2_single_handle = open(tmp_out_files[3], "w")

        for fq1_record, fq2_record in izip(fq1_in, fq2_in):
            if len(fq1_record.seq) >= min_length and len(fq2_record.seq) >= min_length:
                fq1_out_handle.write(fq1_record.format(quality_format))
                fq2_out_handle.write(fq2_record.format(quality_format))
            else:
                if len(fq1_record.seq) > min_length:
                    fq1_single_handle.write(fq1_record.format(quality_format))
                if len(fq2_record.seq) > min_length:
                    fq2_single_handle.write(fq2_record.format(quality_format))
        fq1_out_handle.close()
        fq2_out_handle.close()
        fq1_single_handle.close()
        fq2_single_handle.close()

    return [fq1_out, fq2_out]

def rstrip_extra(fname):
    """Strip extraneous, non-discriminative filename info from the end of a file.
    """
    to_strip = ("_R", "_", "fastq", ".", "-")
    while fname.endswith(to_strip):
        for x in to_strip:
            if fname.endswith(x):
                fname = fname[:len(fname) - len(x)]
                break
    return fname

def combine_pairs(input_files):
    """ calls files pairs if they are completely the same except
    for one has _1 and the other has _2 returns a list of tuples
    of pairs or singles.
    From bipy.utils (https://github.com/roryk/bipy/blob/master/bipy/utils.py)
    Adjusted to allow different input paths or extensions for matching files.
    """
    PAIR_FILE_IDENTIFIERS = set(["1", "2"])

    pairs = []
    used = set([])
    for in_file in input_files:
        if in_file in used:
            continue
        for comp_file in input_files:
            if comp_file in used or comp_file == in_file:
                continue
            a = rstrip_extra(utils.splitext_plus(os.path.basename(in_file))[0])
            b = rstrip_extra(utils.splitext_plus(os.path.basename(comp_file))[0])
            if len(a) != len(b):
                continue
            s = dif(a,b)
            if len(s) > 1:
                continue #there is only 1 difference
            if (a[s[0]] in PAIR_FILE_IDENTIFIERS and
                  b[s[0]] in PAIR_FILE_IDENTIFIERS):
                # if the 1/2 isn't the last digit before a separator, skip
                # this skips stuff like 2P 2A, often denoting replicates, not
                # read pairings
                if len(b) > (s[0] + 1):
                    if (b[s[0]+1] not in ("_", "-", ".")):
                        continue
                # if the 1/2 is not a separator or prefaced with R, skip
                if b[s[0]- 1] in ("R", "_", "-", "."):
                    used.add(in_file)
                    used.add(comp_file)
                    if b[s[0]] == "2":
                        pairs.append([in_file, comp_file])
                    else:
                        pairs.append([comp_file, in_file])
                    break
        if in_file not in used:
            pairs.append([in_file])
            used.add(in_file)
    return pairs

def dif(a, b):
    """ copy from http://stackoverflow.com/a/8545526 """
    return [i for i in range(len(a)) if a[i] != b[i]]

def is_fastq(in_file, bzip=True):
    fastq_ends = [".txt", ".fq", ".fastq"]
    zip_ends = [".gzip", ".gz"]
    if bzip:
        zip_ends += [".bz2", ".bzip2"]
    base, first_ext = os.path.splitext(in_file)
    second_ext = os.path.splitext(base)[1]
    if first_ext in fastq_ends:
        return True
    elif (second_ext, first_ext) in product(fastq_ends, zip_ends):
        return True
    else:
        return False

def downsample(f1, f2, data, N, quick=False):
    """ get N random headers from a fastq file without reading the
    whole thing into memory
    modified from: http://www.biostars.org/p/6544/
    quick=True will just grab the first N reads rather than do a true
    downsampling
    """
    if quick:
        rand_records = range(N)
    else:
        records = sum(1 for _ in open(f1)) / 4
        N = records if N > records else N
        rand_records = random.sample(xrange(records), N)

    fh1 = open_possible_gzip(f1)
    fh2 = open_possible_gzip(f2) if f2 else None
    outf1 = os.path.splitext(f1)[0] + ".subset" + os.path.splitext(f1)[1]
    outf2 = os.path.splitext(f2)[0] + ".subset" + os.path.splitext(f2)[1] if f2 else None

    if utils.file_exists(outf1):
        if not outf2:
            return outf1, outf2
        elif utils.file_exists(outf2):
            return outf1, outf2

    out_files = (outf1, outf2) if outf2 else (outf1)

    with file_transaction(out_files) as tx_out_files:
        if isinstance(tx_out_files, basestring):
            tx_out_f1 = tx_out_files
        else:
            tx_out_f1, tx_out_f2 = tx_out_files
        sub1 = open_possible_gzip(tx_out_f1, "w")
        sub2 = open_possible_gzip(tx_out_f2, "w") if outf2 else None
        rec_no = - 1
        for rr in rand_records:
            while rec_no < rr:
                rec_no += 1
                for i in range(4): fh1.readline()
                if fh2:
                    for i in range(4): fh2.readline()
            for i in range(4):
                sub1.write(fh1.readline())
                if sub2:
                    sub2.write(fh2.readline())
            rec_no += 1
        fh1.close()
        sub1.close()
        if f2:
            fh2.close()
            sub2.close()

    return outf1, outf2

def estimate_read_length(fastq_file, quality_format="fastq-sanger", nreads=1000):
    """
    estimate average read length of a fastq file
    """

    in_handle = SeqIO.parse(fastq_file, quality_format)
    read = in_handle.next()
    average = len(read.seq)
    for _ in range(nreads):
        try:
            average = (average + len(in_handle.next().seq)) / 2
        except StopIteration:
            break
    in_handle.close()
    return average

def open_fastq(in_file):
    """ open a fastq file, using gzip if it is gzipped
    """
    if objectstore.is_remote(in_file):
        return objectstore.open(in_file)
    _, ext = os.path.splitext(in_file)
    if ext == ".gz":
        return gzip.open(in_file, 'rb')
    if ext in [".fastq", ".fq"]:
        return open(in_file, 'r')
    # default to just opening it
    return open(in_file, "r")
