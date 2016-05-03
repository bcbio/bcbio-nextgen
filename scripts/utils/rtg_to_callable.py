#!/bin/env python
"""Convert RTG coverage statistics into a BED file of callable regions.

Provides a supplement to VCF files indicating callable regions where
a no-call indicates reference.

Defaults to requiring more than 4 reads of coverage to be callable.
"""
import sys
import gzip

def main(cov_file):
    min_cov = 4
    out_file = cov_file.replace(".bed.gz", "-callable.bed")
    with gzip.open(cov_file) as in_handle:
        with open(out_file, "w") as out_handle:
            cur_chrom, cur_start, cur_end = None, None, None
            for line in in_handle:
                if not line.startswith("#"):
                    chrom, start, end, _, coverage = line.strip().split()
                    if int(coverage) > min_cov:
                        if chrom == cur_chrom:
                            cur_end = end
                        else:
                            if cur_chrom:
                                out_handle.write("%s\t%s\t%s\n" % (cur_chrom, cur_start, cur_end))
                            cur_chrom, cur_start, cur_end = chrom, start, end
                    elif cur_chrom:
                        out_handle.write("%s\t%s\t%s\n" % (cur_chrom, cur_start, cur_end))
                        cur_chrom, cur_start, cur_end = (None, None, None)

if __name__ == "__main__":
    main(*sys.argv[1:])
