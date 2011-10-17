"""Next-gen sequencing alignment with Novoalign: http://www.novocraft.com
"""
import os
import subprocess

from bcbio.utils import (memoize_outfile, file_exists)
from bcbio.distributed.transaction import file_transaction

@memoize_outfile(".ndx")
def refindex(ref_file, kmer_size=None, step_size=None, out_file=None):
    cl = ["novoindex"]
    if kmer_size:
        cl += ["-k", str(kmer_size)]
    if step_size:
        cl += ["-s", str(step_size)]
    cl += [out_file, ref_file]
    subprocess.check_call(cl)

def _get_base_filename(fname):
    fname = os.path.splitext(os.path.basename(fname))[0]
    to_replace = ["_fastq", "-unique"]
    for rep in to_replace:
        fname = fname.replace(rep, "")
    test_fname, ext = fname.rsplit("_", 1)
    try:
        int(ext)
        fname = test_fname
    except ValueError:
        pass
    return fname

def align(out_dir, ref_index, fastq1, fastq2=None, qual_format=None):
    out_file = os.path.join(out_dir, "%s.sam" % _get_base_filename(fastq1))
    if not file_exists(out_file):
        cl = ["novoalign", "-o", "SAM", "-r", "None", "-d", ref_index, "-f", fastq1]
        if fastq2:
            cl.append(fastq2)
        if qual_format:
            cl += ["-F", qual_format]
        print " ".join(cl)
        with file_transaction(out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                subprocess.check_call(cl, stdout=out_handle)
    return out_file
