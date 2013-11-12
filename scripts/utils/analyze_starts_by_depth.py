import argparse
import pandas as pd

from bcbio import bam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def starts_by_depth(bam_file):
    """
    Return a set of x, y points where x is the number of reads sequenced and
    y is the number of unique start sites identified
    """
    BINSIZE_IN_READS = 100
    seen_starts = set()
    counted = 0
    num_reads = []
    starts = []
    buffer = []
    df = pd.DataFrame(columns=('reads', 'starts'))
    with bam.open_samfile(bam_file) as samfile:
        for read in samfile:
            counted += 1
            buffer.append(":".join([str(read.tid), str(read.pos)]))
            if counted % BINSIZE_IN_READS == 0:
                seen_starts.update(buffer)
                buffer = []
                num_reads.append(counted)
                starts.append(len(seen_starts))
        seen_starts.update(buffer)
        num_reads.append(counted)
        starts.append(len(seen_starts))
    return pd.DataFrame({"reads": num_reads, "starts": starts})

if __name__ == "__main__":
    description = ("Create reads sequenced vs unique start sites graph for "
                   "examining the quality of a library. The idea for this "
                   " metric was borrowed from: "
                   "https://github.com/mbusby/ComplexityByStartPosition")
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("alignment_file", help="Alignment file to process,"
                        "can be SAM or BAM format.")
    parser.add_argument("--out_file", help="Name of output figure.")
    args = parser.parse_args()
    df = starts_by_depth(args.alignment_file)
    df.plot(x='reads', y='starts')
    fig = plt.gcf()
    fig.savefig("test.pdf")
