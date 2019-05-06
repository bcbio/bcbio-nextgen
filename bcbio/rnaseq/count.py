"""
count number of reads mapping to features of transcripts

"""
import os
import pandas as pd
from collections import defaultdict
import gffutils
from bcbio.log import logger

from bcbio.utils import file_exists

def combine_count_files(files, out_file=None, ext=".fpkm"):
    """
    combine a set of count files into a single combined file
    """
    files = list(files)
    if not files:
        return None
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
    logger.info("Combining count files into %s." % out_file)
    row_names = []
    col_vals = defaultdict(list)
    for i, f in enumerate(files):
        vals = []
        if i == 0:
            with open(f) as in_handle:
                for line in in_handle:
                    if not line.strip().startswith("#"):
                        rname, val = line.strip().split("\t")
                        row_names.append(rname)
                        vals.append(val)
        else:
            with open(f) as in_handle:
                for line in in_handle:
                    if not line.strip().startswith("#"):
                        try:
                            _, val = line.strip().split("\t")
                        except ValueError:
                            print(f, line)
                            raise
                        vals.append(val)
        col_vals[col_names[i]] = vals

    df = pd.DataFrame(col_vals, index=row_names)
    df.to_csv(out_file, sep="\t", index_label="id")
    return out_file

def annotate_combined_count_file(count_file, gtf_file, out_file=None):
    if not count_file:
        return None
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

    df = pd.io.parsers.read_csv(count_file, sep="\t", index_col=0, header=0)

    df['symbol'] = df.apply(lambda x: symbol_lookup.get(x.name, ""), axis=1)
    df.to_csv(out_file, sep="\t", index_label="id")
    return out_file
