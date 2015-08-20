"""Hacks to handle subsetting of chromosomes for heterogeneity analysis.

This puts ugly chromosome naming assumptions that restrict heterogeneity estimations
to autosomal chromosomes in a single place.
"""
from bcbio import utils
from bcbio.bam import ref
from bcbio.distributed.transaction import file_transaction

def is_autosomal(chrom):
    """Keep chromosomes that are a digit 1-22, or chr prefixed digit chr1-chr22
    """
    try:
        int(chrom)
        return True
    except ValueError:
        try:
            int(str(chrom.replace("chr", "")))
            return True
        except ValueError:
            return False

def is_autosomal_or_x(chrom):
    return is_autosomal(chrom) or chrom in ["X", "chrX"]

def autosomal_or_x_coords(ref_file):
    out = []
    for contig in ref.file_contigs(ref_file):
        if is_autosomal_or_x(contig.name):
            out.append((contig.name, 0, contig.size))
    return out

def bed_to_standardonly(in_file, data, headers=None):
    out_file = "%s-stdchrs%s" % utils.splitext_plus(in_file)
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(in_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    for line in in_handle:
                        if is_autosomal(line.split()[0]) or (headers and line.startswith(headers)):
                            out_handle.write(line)
    return out_file
