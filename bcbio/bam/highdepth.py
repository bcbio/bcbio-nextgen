"""Identify windows with very high depth for potential filtering.

In non-targeted experiments, high depth regions are often due to collapsed repeats
or other structure which can create long run times and incorrect results in
small and structural variant calling.
"""
import os
import sys

import yaml

from bcbio import utils
from bcbio.pipeline import datadict as dd
from bcbio.variation import bedutils

def _get_files(data):
    work_bam = dd.get_align_bam(data) or dd.get_work_bam(data)
    out_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "align", dd.get_sample_name(data)))
    out_file = "%s-highdepth.bed" % os.path.join(out_dir, utils.splitext_plus(os.path.basename(work_bam))[0])
    stats_file = "%s-stats.yaml" % utils.splitext_plus(out_file)[0]
    return work_bam, out_file, stats_file

def combine_file_rename(x):
    return x.replace("-callable.bed", "-highdepth.bed")

def combine_callable_bed(in_files, out_file, config):
    bedutils.combine(filter(lambda x: utils.file_exists(x), [combine_file_rename(x) for x in in_files]),
                     combine_file_rename(out_file), config)
    callable_bed = bedutils.combine(in_files, out_file, config)
    return callable_bed

def bin_depths(min_cov, max_cov, window_size, callable_out, highdepth_out):
    """Provide bins of covered regions, including a separate file of high depth regions.
    """
    last = (None, None)
    last_window = -1
    depth_cache = []
    cache = []

    def mean(vals):
        return sum(vals) / float(window_size)

    with open(callable_out, "w") as callable_handle:
        with open(highdepth_out, "w") as highdepth_handle:
            for chrom, pos, depth in (l.rstrip("\r\n").split("\t", 3) for l in sys.stdin):
                depth, pos = int(depth), int(pos)
                if pos / window_size != last_window:
                    if mean(depth_cache) > max_cov:
                        highdepth_handle.write("%s\t%d\t%d\t%.2f\n" %
                                               (chrom, last_window * window_size,
                                                last_window * window_size + window_size, mean(depth_cache)))
                    depth_cache = depth_cache[:0]
                    last_window = pos / window_size
                depth_cache.append(depth)

                key = (chrom, "NO_COVERAGE" if depth == 0 else "LOW_COVERAGE" if depth < min_cov else "CALLABLE")
                if key != last or pos != cache[-1][1] + 1:
                    if last[0] is not None:
                        callable_handle.write("\t".join((last[0], str(int(cache[0][1]) - 1),
                                                         str(cache[-1][1]), last[1])) + "\n")
                    last = key
                    cache = cache[:0]
                cache.append((chrom, pos, depth))
            if cache:
                callable_handle.write("\t".join((chrom, str(int(cache[0][1]) - 1),
                                                 str(cache[-1][1]), last[1])) + "\n")

def get_stats_file(data):
    return _get_files(data)[-1]

def get_median_coverage(data):
    stats_file = get_stats_file(data)
    if not utils.file_exists(stats_file):
        return 0
    else:
        with open(stats_file) as in_handle:
            stats = yaml.safe_load(in_handle)
        return stats["median_cov"]
