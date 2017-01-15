"""Utilities to examine BAM counts in defined regions.

These are useful for plotting comparisons between BAM files to look at
differences in defined or random regions.
"""
from __future__ import print_function
import random
import collections

import pysam

class NormalizedBam:
    """Prepare and query an alignment BAM file for normalized read counts.
    """
    def __init__(self, name, fname, picard, quick=False):
        self.name = name
        self._bam = pysam.Samfile(fname, "rb")
        picard.run_fn("picard_index", fname)
        if quick:
            self._total = 1e6
        else:
            self._total = sum(1 for r in self._bam.fetch() if not r.is_unmapped)
            print(name, self._total)

    def all_regions(self):
        """Get a tuple of all chromosome, start and end regions.
        """
        regions = []
        for sq in self._bam.header["SQ"]:
            regions.append((sq["SN"], 1, int(sq["LN"])))
        return regions

    def read_count(self, space, start, end):
        """Retrieve the normalized read count in the provided region.
        """
        read_counts = 0
        for read in self._bam.fetch(space, start, end):
            read_counts += 1
        return self._normalize(read_counts, self._total)

    def coverage_pileup(self, space, start, end):
        """Retrieve pileup coverage across a specified region.
        """
        return ((col.pos, self._normalize(col.n, self._total))
                for col in self._bam.pileup(space, start, end))

    def _normalize(self, count, total):
        """Normalize to reads per million.
        """
        return float(count) / float(total) * 1e6

def random_regions(base, n, size):
    """Generate n random regions of 'size' in the provided base spread.
    """
    spread = size // 2
    base_info = collections.defaultdict(list)
    for space, start, end in base:
        base_info[space].append(start + spread)
        base_info[space].append(end - spread)
    regions = []
    for _ in range(n):
        space = random.choice(base_info.keys())
        pos = random.randint(min(base_info[space]), max(base_info[space]))
        regions.append([space, pos-spread, pos+spread])
    return regions

