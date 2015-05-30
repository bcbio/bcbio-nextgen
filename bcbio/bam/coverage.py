"""
calculate coverage across a list of regions
"""
import os

import six
import matplotlib as mpl
mpl.use('Agg', force=True)
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
import seaborn as sns
import pandas as pd
import pybedtools

from bcbio.utils import rbind, file_exists
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
import bcbio.pipeline.datadict as dd
from pylab import stem, setp
from collections import defaultdict
from itertools import repeat

def _calc_regional_coverage(in_bam, chrom, start, end, samplename, work_dir):
    """
    given a BAM and a region, calculate the coverage for each base in that
    region. returns a pandas dataframe of the format:

    chrom position coverage name

    where the samplename column is the coverage at chrom:position
    """
    region_bt = pybedtools.BedTool("%s\t%s\t%s\n" % (chrom, start, end), from_string=True).saveas()
    region_file = region_bt.fn
    coords = "%s:%s-%s" % (chrom, start, end)
    tx_tmp_file = os.path.join(work_dir, "coverage-%s-%s.txt" % (samplename, coords.replace(":", "_")))
    cmd = ("samtools view -b {in_bam} {coords} | "
           "bedtools coverage -a {region_file} -b - -d > {tx_tmp_file}")
    do.run(cmd.format(**locals()), "Plotting coverage for %s %s" % (samplename, coords))
    names = ["chom", "start", "end", "offset", "coverage"]
    df = pd.io.parsers.read_table(tx_tmp_file, sep="\t", header=None,
                                  names=names).dropna()
    os.remove(tx_tmp_file)
    df["sample"] = samplename
    df["chrom"] = chrom
    df["position"] = df["start"] + df["offset"] - 1
    return df[["chrom", "position", "coverage", "sample"]]

def _combine_regional_coverage(in_bams, samplenames, chrom, start, end, work_dir):
    """
    given a list of bam files, sample names and a region, calculate the
    coverage in the region for each of the samples and return a tidy pandas
    dataframe of the format:

    chrom position coverage name
    """
    dfs = [_calc_regional_coverage(bam, chrom, start, end, sample, work_dir) for bam, sample
           in zip(in_bams, samplenames)]
    return rbind(dfs)

def _get_caller_colormap(callers):
    colors = matplotlib.colors.ColorConverter.colors.keys()
    return {caller: colors[index] for index, caller in enumerate(callers)}

def _get_caller_heights(callers, plot):
    max_y = plot.get_ylim()[1] * 0.2
    spacing = max_y / len(callers)
    return {caller: spacing + spacing * index for index, caller in enumerate(callers)}

def _get_stems_by_callers(intervals):
    stems = defaultdict(list)
    for interval in intervals:
        pos = interval.start
        caller = interval.fields[3]
        stems[caller].append(pos)
    return stems

def _add_stems_to_plot(interval, stem_bed, samples, plot):
    stems = _get_stems_by_callers(stem_bed.tabix_intervals(interval))
    callers = sorted(stems.keys())
    caller_colormap = _get_caller_colormap(callers)
    caller_heights = _get_caller_heights(callers, plot)
    for caller in callers:
        stem_color = caller_colormap[caller]
        caller_stems = stems[caller]
        stem_heights = list(repeat(caller_heights[caller], len(caller_stems)))
        markerline, _, baseline = stem(caller_stems, stem_heights, '-.',
                                       label=caller)
        setp(markerline, 'markerfacecolor', stem_color)
        setp(baseline, 'color', 'r', 'linewidth', 0)
        plt.legend()

def _split_regions(chrom, start, end):
    """Split regions longer than 100kb into smaller sections.
    """
    window_size = 1e5
    if end - start < window_size * 5:
        return [(chrom, start, end)]
    else:
        out = []
        for r in pybedtools.BedTool().window_maker(w=window_size,
                                                   b=pybedtools.BedTool("%s\t%s\t%s" % (chrom, start, end),
                                                                        from_string=True)):
            out.append((r.chrom, r.start, r.end))
        return out

def plot_multiple_regions_coverage(samples, out_file, region_bed=None, stem_bed=None):
    """
    given a list of bcbio samples and a bed file or BedTool of regions,
    makes a plot of the coverage in the regions for the set of samples

    if given a bed file or BedTool of locations in stem_bed with a label,
    plots lollipops at those locations
    """
    PAD = 100
    if file_exists(out_file):
        return out_file
    in_bams = [dd.get_align_bam(x) for x in samples]
    samplenames = [dd.get_sample_name(x) for x in samples]
    if isinstance(region_bed, six.string_types):
        region_bed = pybedtools.BedTool(region_bed)
    if isinstance(stem_bed, six.string_types):
        stem_bed = pybedtools.BedTool(stem_bed)
    if stem_bed is not None:  # tabix indexed bedtools eval to false
        stem_bed = stem_bed.tabix()
    plt.clf()
    plt.cla()
    with file_transaction(out_file) as tx_out_file:
        with PdfPages(tx_out_file) as pdf_out:
            sns.despine()
            for line in region_bed:
                for chrom, start, end in _split_regions(line.chrom, max(line.start - PAD, 0),
                                                        line.end + PAD):
                    df = _combine_regional_coverage(in_bams, samplenames, chrom,
                                                    start, end, os.path.dirname(tx_out_file))
                    plot = sns.tsplot(df, time="position", unit="chrom",
                                      value="coverage", condition="sample")
                    if stem_bed is not None:  # tabix indexed bedtools eval to false
                        interval = pybedtools.Interval(chrom, start, end)
                        _add_stems_to_plot(interval, stem_bed, samples, plot)
                    plt.title("{chrom}:{start}-{end}".format(**locals()))
                    pdf_out.savefig(plot.get_figure())
                    plt.close()
    return out_file
