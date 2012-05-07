"""Split input FASTQ files into pieces to allow parallel cluster processing.

This is useful for speeding up alignments on a cluster at the price of
temporary increased disk usage.
"""
import os
import glob
from bcbio.bam.trim import _save_diskspace

from Bio.SeqIO.QualityIO import FastqGeneralIterator

def _find_current_split(in_fastq, out_dir):
    """Check for existing split files to avoid re-splitting.
    """
    base = os.path.join(out_dir,
                        os.path.splitext(os.path.basename(in_fastq))[0])
    return sorted(glob.glob("{0}*".format(base)))

def _split_by_size(in_fastq, split_size, out_dir):
    """Split FASTQ files by a specified number of records.
    """
    existing = _find_current_split(in_fastq, out_dir)
    if len(existing) > 0:
        return existing
    def new_handle(num):
        base, ext = os.path.splitext(os.path.basename(in_fastq))
        fname = os.path.join(out_dir, "{base}_{num}{ext}".format(
            base=base, num=num, ext=ext))
        return fname, open(fname, "w")
    cur_index = 0
    cur_count = 0
    out_fname, out_handle = new_handle(cur_index)
    out_files = [out_fname]
    with open(in_fastq) as in_handle:
        for name, seq, qual in FastqGeneralIterator(in_handle):
            if cur_count < split_size:
                cur_count += 1
            else:
                cur_count = 0
                cur_index += 1
                out_handle.close()
                out_fname, out_handle = new_handle(cur_index)
                out_files.append(out_fname)
            out_handle.write("@%s\n%s\n+\n%s\n" % (name, seq, qual))
    out_handle.close()
    return out_files

def split_fastq_files(fastq1, fastq2, split_size, out_dir, config):
    """Split paired end FASTQ files into pieces for parallel analysis.
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    split_fastq1 = _split_by_size(fastq1, split_size, out_dir)
    _save_diskspace(fastq1, split_fastq1[0], config)
    if fastq2:
        split_fastq2 = _split_by_size(fastq2, split_size, out_dir)
        _save_diskspace(fastq2, split_fastq2[0], config)
    else:
        split_fastq2 = [None] * len(split_fastq1)
    return zip(split_fastq1, split_fastq2, [None] + [x+1 for x in range(len(split_fastq1) - 1)])
