"""Utilities to manipulate VCF files.
"""
import os

from bcbio.utils import file_exists, replace_suffix, append_stem
from bcbio.distributed.transaction import file_transaction

import sh

def split_vcf(in_file, config, out_dir=None):
    """
    split a VCF file into separate files by chromosome
    requires tabix to be installed

    """
    if out_dir is None:
        out_dir = os.path.join(os.path.dirname(in_file), "split")

    fasta_file = config["ref"]["fasta"]
    fasta_index = fasta_file + ".fai"
    samtools_path = config["program"].get("samtools", "samtools")
    tabix_path = config["program"].get("tabix", "tabix")

    if not file_exists(fasta_index):
        samtools = sh.Command(samtools_path)
        samtools.faidx(fasta_file)

    # if in_file is not compressed, compress it
    (_, ext) = os.path.splitext(in_file)
    if ext is not ".gz":
        gzip_file = in_file + ".gz"
        if not file_exists(gzip_file):
            sh.bgzip("-c", in_file, _out=gzip_file)
        in_file = gzip_file

    # create tabix index
    tabix_index(in_file)

    # find the chromosome names from the fasta index file
    chroms = str(sh.cut("-f1", fasta_index)).split()

    # make outfile from chromosome name
    def chr_out(chrom):
        out_file = replace_suffix(append_stem(in_file, chrom), ".vcf")
        return os.path.join(out_dir, os.path.basename(out_file))

    # run tabix to break up the vcf file
    def run_tabix(chrom):
        tabix = sh.Command(tabix_path)
        out_file = chr_out(chrom)
        if file_exists(out_file):
            return out_file
        with file_transaction(out_file) as tmp_out_file:
            tabix("-h", in_file, chrom, _out=tmp_out_file)
        return out_file

    out_files = map(run_tabix, chroms)
    return out_files


def tabix_index(in_file, preset="vcf", config=None):
    """
    index a file using tabix

    """
    if config:
        tabix_path = config["program"].get("tabix", "tabix")
    else:
        tabix_path = sh.which("tabix")
    tabix = sh.Command(tabix_path)
    out_file = in_file + ".tbi"
    if file_exists(out_file):
        return out_file
    tabix("-p", preset, in_file)

    return out_file
