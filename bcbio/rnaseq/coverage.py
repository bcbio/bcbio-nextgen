"""
Functions to handle plotting coverage across genes
"""
try:
    import gffutils
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg', force=True)
    import matplotlib.pyplot as plt
    plt.ioff()
except ImportError:
    gffutils, pd, plt = None, None, None

import random
import numpy as np
from collections import defaultdict, Counter

from bcbio.utils import file_exists

def _select_random_nonzero_genes(count_file):
    """
    given a count file with rows of gene_ids and columns of counts
    return a random set of genes with non-zero counts
    """
    MIN_READS = 100
    DEFAULT_NUM_SAMPLES = 100
    df = pd.io.parsers.read_csv(count_file, delimiter="\t", header=0, index_col=0)
    means = pd.DataFrame({"mean": df.mean(1)}, index=df.index)
    means = means[means['mean'] > MIN_READS]
    NUM_SAMPLES = min(DEFAULT_NUM_SAMPLES, len(means))
    rows = random.sample(means.index, NUM_SAMPLES)
    return list(means.ix[rows].index)

def _plot_coverage(df, out_file):
    fig = plt.gcf()
    df.plot(x='distance', y='depths', subplots=True)
    fig.savefig(out_file)
    plt.close(fig)
    return out_file

def _normalize_coverage(read_depths):
    """
    given a list of read depths for a gene, scales read depth to
    a 100 bp faux-gene so multiple genes can be averaged
    together to get an overall view of coverage for a set of genes
    """
    gene_length = len(read_depths)
    norm_dist = [100 * float(x) / gene_length for x in range(gene_length)]
    df = pd.DataFrame({"distance": norm_dist, "depths": read_depths})
    return df


def _gene_depth(dbfn, bamfn, gene):
    """
    takes a gffutils db (dbfn), a BAM file (bamfn) and a gene_id (gene)
    and returns a 5' -> 3' list of per-base read depths for the exons of
    the gene
    """
    from chanjo import bam
    db = gffutils.FeatureDB(dbfn, keep_order=True)
    read_depths = []
    bam_handle = bam.CoverageAdapter(bamfn)
    for exon in db.children(gene, featuretype="exon", order_by='start'):
        strand = exon.strand
        coord = [exon.start, exon.end]
        read_depths += bam_handle.read(exon.seqid, min(coord), max(coord)).tolist()
    # return a list of depths going in the 5' -> 3' direction
    if strand == "-":
        read_depths = read_depths[::-1]
    return read_depths

def plot_gene_coverage(bam_file, ref_file, count_file, out_file):
    if file_exists(out_file):
        return out_file
    coverage = pd.DataFrame()
    ref_db = ref_file + ".db"
    for gene in _select_random_nonzero_genes(count_file):
        depth = _gene_depth(ref_db, bam_file, gene)
        coverage = coverage.append(_normalize_coverage(depth))

    # group into 100 bins for 0->100% along the transcript
    if coverage.empty:
        return None

    groups = coverage.groupby(np.digitize(coverage.distance, range(100)))
    out_file = _plot_coverage(groups.mean(), out_file)
    return out_file

def estimate_library_content(bam_file, ref_file):
    from chanjo import bam
    ref_db = ref_file + ".db"
    library_content = defaultdict(Counter)
    db = gffutils.FeatureDB(ref_db, keep_order=True)
    with bam.open_samfile(bam_file) as bam_handle:
        for read in bam_handle:
            name = read.getrname(read.tid)
            start = read.pos
            end = read.aend
            overlapped = db.region((name, start, end))
