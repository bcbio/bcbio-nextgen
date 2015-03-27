"""Annotate structural variant calls with associated genes.
"""
from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do

def add_genes(in_file, data):
    """Add gene annotations to a BED file from pre-prepared RNA-seq data.
    """
    import pybedtools
    gene_file = dd.get_gene_bed(data)
    if gene_file:
        out_file = "%s-annotated.bed" % utils.splitext_plus(in_file)[0]
        if not utils.file_uptodate(out_file, in_file):
            input_rec = iter(pybedtools.BedTool(in_file)).next()
            # keep everything after standard chrom/start/end, 1-based
            extra_fields = range(4, len(input_rec.fields) + 1)
            # keep the new gene annotation
            extra_fields.append(len(input_rec.fields) + 4)
            columns = ",".join([str(x) for x in extra_fields])
            ops = ",".join(["distinct"] * len(extra_fields))
            with file_transaction(data, out_file) as tx_out_file:
                cmd = ("sort -k1,1 -k2,2n {in_file} | "
                       "bedtools closest -t all -a - -b {gene_file} | "
                       "bedtools merge -i - -c {columns} -o {ops} > {tx_out_file}")
                do.run(cmd.format(**locals()), "Annotate BED file with gene info")
        return out_file
    else:
        return in_file
