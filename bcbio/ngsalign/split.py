"""Split input FASTQ files into pieces to allow parallel cluster processing.

This is useful for speeding up alignments on a cluster at the price of
temporary increased disk usage.
"""
import os
import glob
import itertools
import operator
import time

import pysam
from Bio import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from bcbio.bam.trim import _save_diskspace
from bcbio.bam import cram
from bcbio import utils, broad
from bcbio.pipeline import config_utils, alignment

def _find_current_split(in_fastq, out_dir):
    """Check for existing split files to avoid re-splitting.
    """
    base = os.path.join(out_dir,
                        os.path.splitext(os.path.basename(in_fastq))[0])
    def get_splitnum(fname):
        """Number from filename like: NA12878-E2-XPR855_2_69.fastq
        """
        base = os.path.splitext(os.path.basename(fname))[0]
        _, num = base.rsplit("_", 1)
        return int(num)
    return sorted(glob.glob("{0}*".format(base)), key=get_splitnum)

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

def _get_seq_qual(read):
    if read.is_reverse:
        seq = str(Seq.Seq(read.seq).reverse_complement())
        tmp = list(read.qual)
        tmp.reverse()
        qual = "".join(tmp)
    else:
        seq = read.seq
        qual = read.qual
    return seq, qual

def _find_current_bam_split(bam_file, out_dir):
    """Check for existing split files from BAM inputs, to avoid re-splitting.
    """
    base = os.path.join(out_dir,
                        os.path.splitext(os.path.basename(bam_file))[0])
    def get_pair_and_splitnum(fname):
        base = os.path.splitext(os.path.basename(fname))[0]
        _, pair, num = base.rsplit("_", 2)
        return int(num), int(pair)
    xs = []
    for fname in glob.glob("{0}_*".format(base)):
        num, pair = get_pair_and_splitnum(fname)
        xs.append((num, pair, fname))
    out = []
    for num, g in itertools.groupby(sorted(xs), operator.itemgetter(0)):
        f1, f2 = [x[-1] for x in sorted(g)]
        split = num if num > 0 else None
        out.append((f1, f2, split))
    return out

def split_bam_file(bam_file, split_size, out_dir, config):
    """Split a BAM file into paired end fastq splits based on split size.

    XXX Need to generalize for non-paired end inputs.
    """
    existing = _find_current_bam_split(bam_file, out_dir)
    if len(existing) > 0:
        return existing
    pipe = True

    utils.safe_makedir(out_dir)
    broad_runner = broad.runner_from_config(config)
    out_files = []
    def new_handle(num):
        out = []
        for pair in [1, 2]:
            fname = os.path.join(out_dir, "{base}_{pair}_{num}.fastq".format(
                base=os.path.splitext(os.path.basename(bam_file))[0], pair=pair, num=num))
            out += [fname, open(fname, "w")]
        return out
    with utils.curdir_tmpdir(base_dir=config_utils.get_resources("tmp", config).get("dir")) as tmp_dir:
        if pipe:
            sort_file = os.path.join(tmp_dir, "%s-sort.bam" %
                                     os.path.splitext(os.path.basename(bam_file))[0])
            os.mkfifo(sort_file)
            broad_runner.run_fn("picard_sort", bam_file, "queryname", sort_file,
                                compression_level=0, pipe=True)
        else:
            sort_file = os.path.join(out_dir, "%s-sort.bam" %
                                     os.path.splitext(os.path.basename(bam_file))[0])
            broad_runner.run_fn("picard_sort", bam_file, "queryname", sort_file)

        samfile = pysam.Samfile(sort_file, "rb")
        i = 0
        num = 0
        f1, out_handle1, f2, out_handle2 = new_handle(num)
        out_files.append([f1, f2, None])
        for x1, x2 in utils.partition_all(2, samfile):
            x1_seq, x1_qual = _get_seq_qual(x1)
            out_handle1.write("@%s/1\n%s\n+\n%s\n" % (i, x1_seq, x1_qual))
            x2_seq, x2_qual = _get_seq_qual(x2)
            out_handle2.write("@%s/2\n%s\n+\n%s\n" % (i, x2_seq, x2_qual))
            i += 1
            if i % split_size == 0:
                num += 1
                out_handle1.close()
                out_handle2.close()
                f1, out_handle1, f2, out_handle2 = new_handle(num)
                out_files.append([f1, f2, num])
        out_handle1.close()
        out_handle2.close()
        samfile.close()
        if pipe:
            os.unlink(sort_file)
        else:
            utils.save_diskspace(sort_file, "Split to {}".format(out_files[0][0]), config)
    return out_files

def split_read_files(fastq1, fastq2, item, split_size, out_dir, dirs, config):
    """Split input reads for parallel processing, dispatching on input type.
    """
    if fastq1.endswith(".bam") and fastq2 is None:
        qual_bin_method = config["algorithm"].get("quality_bin")
        if (qual_bin_method == "prealignment" or
             (isinstance(qual_bin_method, list) and "prealignment" in qual_bin_method)):
            _, sam_ref = alignment.get_genome_ref(item["genome_build"], None, dirs["galaxy"])
            out_bindir = utils.safe_makedir(os.path.join(out_dir, "qualbin"))
            fastq1 = cram.illumina_qual_bin(fastq1, sam_ref, out_bindir, config)
        return split_bam_file(fastq1, split_size, out_dir, config)
    else:
        return split_fastq_files(fastq1, fastq2, split_size, out_dir, config)
