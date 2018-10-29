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
        total = 0
        mapped = 0
        duplicates = 0
        umi_reductions = []
        umi_counts = collections.defaultdict(int)
        with pysam.AlignmentFile(data["umi_bam"], "rb", check_sq=False) as bam_iter:
            cur_counts = collections.defaultdict(int)
            cur_key = None
            for rec in bam_iter:
                total += 1
                umi = _get_umi_tag(rec)
                if umi and not rec.is_unmapped:
                    mapped += 1
                    if rec.is_duplicate:
                        duplicates += 1
                    chrom = bam_iter.getrname(rec.reference_id)
                    pos = rec.reference_start
                    key = (chrom, pos)
                    if key != cur_key:
                        # update counts
                        if cur_counts:
                            for c in cur_counts.values():
                                umi_counts[c] += 1
                            total_seqs = sum(cur_counts.values())
                            umi_count = len(cur_counts)
                            umi_reductions.append(float(total_seqs) / umi_count)
                        # update current keys
                        cur_key = key
                        cur_counts = collections.defaultdict(int)
                    cur_counts[umi] += 1
            if cur_counts:
                for c in cur_counts.values():
                    umi_counts[c] += 1
                total_seqs = sum(cur_counts.values())
                umi_count = len(cur_counts)
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
        out["umi_raw_avg_cov"] = data["config"]["algorithm"].get("rawumi_avg_cov", 0)
        with open(stats_file, "w") as out_handle:
            yaml.safe_dump({dd.get_sample_name(data): out}, out_handle,
                           default_flow_style=False, allow_unicode=False)
    return stats_file

def _get_umi_tag(rec):
    """Handle UMI and duplex tag retrieval.
    """
    for tag in ["RX", "XC"]:
        try:
            return rec.get_tag(tag)
        except KeyError:
            pass
