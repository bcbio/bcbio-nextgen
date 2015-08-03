import argparse
import os
from bcbio.rnaseq import qc
from collections import Counter
import bcbio.bam as bam
import bcbio.utils as utils
from itertools import ifilter
import bcbio.pipeline.datadict as dd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def count_duplicate_starts(bam_file, sample_size=10000000):
    """
    Return a set of x, y points where x is the number of reads sequenced and
    y is the number of unique start sites identified
    If sample size < total reads in a file the file will be downsampled.
    """
    count = Counter()
    with bam.open_samfile(bam_file) as samfile:
        # unmapped reads should not be counted
        filtered = ifilter(lambda x: not x.is_unmapped, samfile)
        def read_parser(read):
            return ":".join([str(read.tid), str(read.pos)])
        samples = utils.reservoir_sample(filtered, sample_size, read_parser)

    count.update(samples)
    return count


if __name__ == "__main__":
    description = ("Create reads sequenced vs unique start sites graph for "
                   "examining the quality of a library. The idea for this "
                   " metric was borrowed from: "
                   "https://github.com/mbusby/ComplexityByStartPosition."
                   "This can also be used to generate counts for use with "
                   "the preseq tool: "
                   "http://smithlab.usc.edu/plone/software/librarycomplexity.")

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("alignment_file", help="Alignment file to process,"
                        "can be SAM or BAM format.")
    parser.add_argument("--complexity", default=False, action='store_true',
                        help="Rough estimate of library complexity")
    parser.add_argument("--histogram", default=False, action='store_true',
                        help="Output a histogram of the unique reads vs. counts.")
    parser.add_argument("--counts", default=False, action='store_true',
                        help="Output a counts of each start site")
    parser.add_argument("--figure", default=None, help="Generate a figure for the complexity")
    parser.add_argument("--sample-size", default=50000000, type=int,
                        help="Number of reads to sample.")
    args = parser.parse_args()
    data = dd.set_config({}, {})
    df = qc.starts_by_depth(args.alignment_file, data, args.sample_size)
    if args.figure:
        df.plot(x='reads', y='starts')
        fig = plt.gcf()
        fig.savefig(args.figure)

    if args.histogram:
        base, _ = os.path.splitext(args.alignment_file)
        df.to_csv(base + ".histogram", sep="\t", header=True, index=False)
    if args.complexity:
        print qc.estimate_library_complexity(df)
    if args.counts:
        c = count_duplicate_starts(args.alignment_file)
        for item in c.items():
            print "{0}\t{1}".format(item[0], item[1])
