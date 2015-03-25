"""
calculate coverage across a list of regions
"""
import tempfile
import subprocess
import sys
import pandas as pd
if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from bcbio.utils import rbind, file_exists
from bcbio.distributed.transaction import file_transaction

def _calc_regional_coverage(in_bam, chrom, start, end, samplename):
    """
    given a BAM and a region, calculate the coverage for each base in that
    region. returns a pandas dataframe of the format:

    chrom position coverage name

    where the samplename column is the coverage at chrom:position
    """
    region_file = tempfile.NamedTemporaryFile(delete=False).name
    with open(region_file, "w") as out_handle:
        out_handle.write("%s\t%s\t%s\n" % (chrom, start, end))

    cmd = ("bedtools coverage -abam {in_bam} -b {region_file} -d")
    out = subprocess.check_output(cmd.format(**locals()), shell=True)
    names = ["chom", "start", "end", "offset", "coverage"]
    df = pd.io.parsers.read_table(StringIO(out), sep="\t", header=None,
                                  names=names)
    df["sample"] = samplename
    df["chrom"] = chrom
    df["position"] = df["start"] + df["offset"] - 1
    return df[["chrom", "position", "coverage", "sample"]]

def _combine_regional_coverage(in_bams, samplenames, chrom, start, end):
    """
    given a list of bam files, sample names and a region, calculate the
    coverage in the region for each of the samples and return a tidy pandas
    dataframe of the format:

    chrom position coverage name
    """
    dfs = [_calc_regional_coverage(bam, chrom, start, end, sample) for bam, sample
           in zip(in_bams, samplenames)]
    return rbind(dfs)

def _bed_to_regions(bed_file):
    regions = []
    with open(bed_file) as in_handle:
        for line in in_handle:
            tokens = line.split()
            regions.append((tokens[0], tokens[1], tokens[2]))
    return regions

def plot_multiple_regions_coverage(samples, out_file, region_bed=None, regions=None):
    """
    given a list of bcbio samples and a bed file of regions or a list of tuples
    of regions the form (chrom, start, end) make a plot of the coverage in the
    regions for the set of samples
    """
    assert (region_bed or regions), ("A region bed file or a list of tuples of "
                                     "regions must be provided")
    if file_exists(out_file):
        return out_file
    in_bams = [dd.get_work_bam(x) for x in samples]
    samplenames = [dd.get_sample_name(x) for x in samples]
    if not regions:
        regions = _bed_to_regions(region_bed)
    with file_transaction(out_file) as tx_out_file:
        with PdfPages(tx_out_file) as pdf_out:
            sns.despine()
            for chrom, start, end in regions:
                df = _combine_regional_coverage(in_bams, samplenames, chrom,
                                                start, end)
                plot = sns.tsplot(df, time="position", unit="chrom",
                                  value="coverage", condition="sample")
                plt.title("{chrom}:{start}-{end}".format(**locals()))
                pdf_out.savefig(plot.get_figure())
                plt.close()
    return out_file
