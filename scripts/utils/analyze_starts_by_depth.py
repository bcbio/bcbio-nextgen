import argparse

from bcbio.rnaseq import qc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


if __name__ == "__main__":
    description = ("Create reads sequenced vs unique start sites graph for "
                   "examining the quality of a library. The idea for this "
                   " metric was borrowed from: "
                   "https://github.com/mbusby/ComplexityByStartPosition")
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("alignment_file", help="Alignment file to process,"
                        "can be SAM or BAM format.")
    parser.add_argument("--complexity", default=False, action='store_true',
                        help="Rough estimate of library complexity")
    parser.add_argument("--figure", default=None, help="Generate a figure for the complexity")
    args = parser.parse_args()
    if args.figure:
        df = qc.starts_by_depth(args.alignment_file)
        df.plot(x='reads', y='starts')
        fig = plt.gcf()
        fig.savefig(args.figure)

    if args.complexity:
        print qc.estimate_library_complexity(args.alignment_file)
