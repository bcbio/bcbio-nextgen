"""Hacks to handle subsetting of chromosomes for heterogeneity analysis.

This puts ugly chromosome naming assumptions that restrict heterogeneity estimations
to autosomal chromosomes in a single place.
"""
import os

from bcbio import utils
from bcbio.distributed.transaction import file_transaction

def is_autosomal(chrom):
    """Keep chromosomes that are a digit 1-22, or chr prefixed digit chr1-chr22
    """
    try:
        int(chrom)
        return True
    except ValueError:
        try:
            int(str(chrom.lower().replace("chr", "").replace("_", "").replace("-", "")))
            return True
        except ValueError:
            return False

def is_sex(chrom):
    return chrom in ["X", "chrX", "Y", "chrY"]

def is_mitochondrial(chrom):
    return chrom in ["MT", "chrM", "chrMT"]

def is_autosomal_or_x(chrom):
    return is_autosomal(chrom) or chrom in ["X", "chrX"]

def is_autosomal_or_sex(chrom):
    return is_autosomal(chrom) or is_sex(chrom)

def is_nonalt(chrom):
    """Check that a chromosome is on 1-22, X, Y, MT.
    """
    return is_autosomal_or_sex(chrom) or is_mitochondrial(chrom)

def bed_to_standardonly(in_file, data, headers=None, include_sex_chroms=False, out_dir=None):
    out_file = "%s-stdchrs%s" % utils.splitext_plus(in_file)
    if out_dir:
        out_file = os.path.join(out_dir, os.path.basename(out_file))
    checkfn = is_autosomal_or_sex if include_sex_chroms else is_autosomal
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(in_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    for line in in_handle:
                        if checkfn(line.split()[0]) or (headers and line.startswith(headers)):
                            out_handle.write(line)
    return out_file
