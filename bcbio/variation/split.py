"""Utilities for manipulating VCF files.
"""
import os

from bcbio.utils import file_exists, replace_suffix, append_stem
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils, tools
from bcbio.provenance import do
from bcbio.variation import bamprep

def split_vcf(in_file, ref_file, config, out_dir=None):
    """Split a VCF file into separate files by chromosome.
    """
    if out_dir is None:
        out_dir = os.path.join(os.path.dirname(in_file), "split")
    out_files = []
    with open(fasta_idx(ref_file, config)) as in_handle:
        for line in in_handle:
            chrom, size = line.split()[:2]
            out_file = os.path.join(out_dir,
                                    os.path.basename(replace_suffix(append_stem(in_file, "-%s" % chrom), ".vcf")))
            subset_vcf(in_file, (chrom, 0, size), out_file, config)
            out_files.append(out_file)
    return out_files

def subset_vcf(in_file, region, out_file, config):
    """Subset VCF in the given region, handling bgzip and indexing of input.
    """
    work_file = bgzip_and_index(in_file, config)
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            bcftools = config_utils.get_program("bcftools", config)
            region_str = bamprep.region_to_gatk(region)
            cmd = "{bcftools} subset -r {region_str} {work_file} > {tx_out_file}"
            do.run(cmd.format(**locals()), "subset %s: %s" % (os.path.basename(work_file), region_str))
    return out_file

def fasta_idx(in_file, config):
    """Retrieve fasta index
    """
    fasta_index = in_file + ".fai"
    if not file_exists(fasta_index):
        samtools = config_utils.get_program("samtools", config)
        cmd = "{samtools} faidx {in_file}"
        do.run(cmd.format(**locals()), "samtools faidx")
    return fasta_index

def bgzip_and_index(in_file, config):
    """bgzip and tabix index an input VCF file.
    """
    out_file = in_file if in_file.endswith(".gz") else in_file + ".gz"
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            bgzip = tools.get_bgzip_cmd(config)
            cmd = "{bgzip} -c {in_file} > {tx_out_file}"
            do.run(cmd.format(**locals()), "bgzip %s" % os.path.basename(in_file))
        os.remove(in_file)
    tabix_index(out_file, config)
    return out_file

def tabix_index(in_file, config, preset="vcf"):
    """Index a file using tabix.
    """
    in_file = os.path.abspath(in_file)
    out_file = in_file + ".tbi"
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            tabix = tools.get_tabix_cmd(config)
            tx_in_file = os.path.splitext(tx_out_file)[0]
            os.symlink(in_file, tx_in_file)
            cmd = "{tabix} -p {preset} {tx_in_file}"
            do.run(cmd.format(**locals()), "tabix index %s" % os.path.basename(in_file))
    return out_file
