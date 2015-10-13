"""
count number of reads mapping to features of transcripts

"""
import os
import sys
import itertools
import pandas as pd
import gffutils

from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio import bam
import bcbio.pipeline.datadict as dd

def combine_count_files(files, out_file=None, ext=".fpkm"):
    """
    combine a set of count files into a single combined file
    """
    assert all([file_exists(x) for x in files]), \
        "Some count files in %s do not exist." % files
    for f in files:
        assert file_exists(f), "%s does not exist or is empty." % f
    col_names = [os.path.basename(x.replace(ext, "")) for x in files]
    if not out_file:
        out_dir = os.path.join(os.path.dirname(files[0]))
        out_file = os.path.join(out_dir, "combined.counts")

    if file_exists(out_file):
        return out_file

    for i, f in enumerate(files):
        if i == 0:
            df = pd.io.parsers.read_table(f, sep="\t", index_col=0, header=None,
                                          names=[col_names[0]])
        else:
            df = df.join(pd.io.parsers.read_table(f, sep="\t", index_col=0,
                                                  header=None,
                                                  names=[col_names[i]]))

    df.to_csv(out_file, sep="\t", index_label="id")
    return out_file

def annotate_combined_count_file(count_file, gtf_file, out_file=None):
    dbfn = gtf_file + ".db"
    if not file_exists(dbfn):
        return None

    if not gffutils:
        return None

    db = gffutils.FeatureDB(dbfn, keep_order=True)

    if not out_file:
        out_dir = os.path.dirname(count_file)
        out_file = os.path.join(out_dir, "annotated_combined.counts")

    # if the genes don't have a gene_id or gene_name set, bail out
    try:
        symbol_lookup = {f['gene_id'][0]: f['gene_name'][0] for f in
                         db.features_of_type('exon')}
    except KeyError:
        return None

    df = pd.io.parsers.read_table(count_file, sep="\t", index_col=0, header=0)

    df['symbol'] = df.apply(lambda x: symbol_lookup.get(x.name, ""), axis=1)
    df.to_csv(out_file, sep="\t", index_label="id")
    return out_file
