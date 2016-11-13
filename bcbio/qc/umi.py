"""Calculate quality control metrics for UMI tags and consensus generation.
"""
import collections
import math
import os

import numpy as np
import pysam
import yaml

from bcbio import bam, utils
from bcbio.pipeline import datadict as dd

def run(_, data, out_dir):
    stats_file = os.path.join(utils.safe_makedir(out_dir), "%s_umi_stats.yaml" % dd.get_sample_name(data))
    if not utils.file_uptodate(stats_file, dd.get_align_bam(data)):
        out = {}
        counts = collections.defaultdict(lambda: collections.defaultdict(int))
        total = 0
        mapped = 0
        duplicates = 0
        with pysam.AlignmentFile(data["umi_bam"], "rb", check_sq=False) as bam_iter:
            for rec in bam_iter:
                total += 1
                umi = rec.get_tag("RX")
                if umi and not rec.is_unmapped:
                    mapped += 1
                    if rec.is_duplicate:
                        duplicates += 1
                    chrom = bam_iter.getrname(rec.reference_id)
                    pos = rec.reference_start
                    key = (chrom, pos)
                counts[key][umi] += 1
        umi_reductions = []
        umi_counts = collections.defaultdict(int)
        for key in sorted(counts.keys()):
            for c in counts[key].values():
                umi_counts[c] += 1
            total_seqs = sum(counts[key].values())
            umi_count = len(counts[key])
            umi_reductions.append(float(total_seqs) / umi_count)
        consensus_count = sum([x.aligned for x in bam.idxstats(dd.get_align_bam(data), data)])
        out["umi_baseline_all"] = total
        out["umi_baseline_mapped"] = mapped
        out["umi_baseline_duplicate_pct"] = float(duplicates) / float(mapped) * 100.0
        out["umi_consensus_mapped"] = consensus_count
        out["umi_consensus_pct"] = (100.0 - float(consensus_count) / float(mapped) * 100.0)
        out["umi_reduction_median"] = int(math.ceil(np.median(umi_reductions)))
        out["umi_reduction_max"] = int(max(umi_reductions))
        out["umi_counts"] = dict(umi_counts)
        with open(stats_file, "w") as out_handle:
            yaml.safe_dump({dd.get_sample_name(data): out}, out_handle,
                           default_flow_style=False, allow_unicode=False)
    return stats_file
