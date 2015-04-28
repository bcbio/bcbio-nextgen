"""Annotate structural variant calls with associated genes.
"""
import os

from bcbio import utils
from bcbio.bam import ref
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do

import pybedtools

def add_genes(in_file, data, max_distance=10000):
    """Add gene annotations to a BED file from pre-prepared RNA-seq data.

    max_distance -- only keep annotations within this distance of event
    """
    gene_file = dd.get_gene_bed(data)
    if gene_file and utils.file_exists(in_file):
        out_file = "%s-annotated.bed" % utils.splitext_plus(in_file)[0]
        if not utils.file_uptodate(out_file, in_file):
            input_rec = iter(pybedtools.BedTool(in_file)).next()
            # keep everything after standard chrom/start/end, 1-based
            extra_fields = range(4, len(input_rec.fields) + 1)
            # keep the new gene annotation
            gene_index = len(input_rec.fields) + 4
            extra_fields.append(gene_index)
            columns = ",".join([str(x) for x in extra_fields])
            max_column = max(extra_fields) + 1
            ops = ",".join(["distinct"] * len(extra_fields))
            with file_transaction(data, out_file) as tx_out_file:
                # swap over gene name to '.' if beyond maximum distance
                # cut removes the last distance column which can cause issues
                # with bedtools merge: 'ERROR: illegal character '.' found in integer conversion of string'
                distance_filter = (r"""awk -F$'\t' -v OFS='\t' '{if ($NF > %s) $%s = "."} {print}'""" %
                                   (max_distance, gene_index))
                cmd = ("sort -k1,1 -k2,2n {in_file} | "
                       "bedtools closest -d -t all -a - -b {gene_file} | "
                       "{distance_filter} | cut -f 1-{max_column} | "
                       "bedtools merge -i - -c {columns} -o {ops} -delim ',' > {tx_out_file}")
                do.run(cmd.format(**locals()), "Annotate BED file with gene info")
        return out_file
    else:
        return in_file

def subset_by_genes(in_file, data, out_dir, pad):
    """Subset BED file of regions to only those within pad of the final output.
    """
    gene_file = dd.get_gene_bed(data)
    fai_file = ref.fasta_idx(dd.get_ref_file(data))
    if not gene_file or not utils.file_exists(in_file):
        return in_file
    else:
        out_file = os.path.join(out_dir, "%s-geneonly.bed" % utils.splitext_plus(os.path.basename(in_file))[0])
        if not utils.file_uptodate(out_file, in_file):
            with file_transaction(data, out_file) as tx_out_file:
                want_region_file = "%s-targetregions%s" % utils.splitext_plus(out_file)
                pybedtools.BedTool(gene_file).slop(g=fai_file, b=pad).merge().saveas(want_region_file)
                pybedtools.BedTool(in_file).intersect(b=want_region_file).sort().saveas(tx_out_file)
        return out_file