"""Utilities for manipulating VCF files.
"""
import os

from bcbio.bam import ref
from bcbio.utils import file_exists, replace_suffix, append_stem
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.variation import bamprep, vcfutils

def split_vcf(in_file, ref_file, config, out_dir=None):
    """Split a VCF file into separate files by chromosome.
    """
    if out_dir is None:
        out_dir = os.path.join(os.path.dirname(in_file), "split")
    out_files = []
    with open(ref.fasta_idx(ref_file, config)) as in_handle:
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
    work_file = vcfutils.bgzip_and_index(in_file, config)
    if not file_exists(out_file):
        with file_transaction(config, out_file) as tx_out_file:
            bcftools = config_utils.get_program("bcftools", config)
            region_str = bamprep.region_to_gatk(region)
            cmd = "{bcftools} view -r {region_str} {work_file} > {tx_out_file}"
            do.run(cmd.format(**locals()), "subset %s: %s" % (os.path.basename(work_file), region_str))
    return out_file
