#!/usr/bin/env python
"""Collect alignment summary metrics from multiple lanes and summarize as CSV.

Usage:
    collect_metrics_to_csv.py <comma separated list of lanes>

"""
import sys
import os
import csv
import glob
import collections

import pysam

from Mgh.Picard.metrics import PicardMetricsParser

WANT_METRICS = [
"AL_TOTAL_READS",
"AL_PF_READS_ALIGNED",
"AL_PF_HQ_ALIGNED_Q20_BASES",
"AL_PCT_READS_ALIGNED_IN_PAIRS",
"DUP_READ_PAIR_DUPLICATES",
"DUP_PERCENT_DUPLICATION",
"DUP_ESTIMATED_LIBRARY_SIZE",
"HS_BAIT_SET",
"HS_GENOME_SIZE",
"HS_LIBRARY_SIZE",
"HS_BAIT_TERRITORY",
"HS_TARGET_TERRITORY",
"HS_PF_UQ_BASES_ALIGNED",
"HS_ON_BAIT_BASES",
"HS_ON_TARGET_BASES",
"HS_PCT_SELECTED_BASES",
"HS_MEAN_TARGET_COVERAGE",
"HS_FOLD_ENRICHMENT",
"HS_ZERO_CVG_TARGETS_PCT",
"HS_FOLD_80_BASE_PENALTY",
"HS_PCT_TARGET_BASES_2X",
"HS_PCT_TARGET_BASES_10X",
"HS_PCT_TARGET_BASES_20X",
"HS_PENALTY_20X",
# ToDo
"SNP_TOTAL_SNPS",
"SNP_PCT_DBSNP",
"Lane IC PCT Mean RD1 Err Rate",
"Lane IC PCT Mean RD2 Err Rate",
]

def main(run_name):
    work_dir = os.getcwd()
    parser = PicardMetricsParser()

    base_bams = _get_base_bams(work_dir, run_name)
    samples = [_get_sample_name(b) for (_, _, b) in base_bams]
    metrics = [lane_stats(l, b, run_name, work_dir, parser)
               for (l, b, _) in base_bams]
    header_counts = _get_header_counts(samples, metrics)
    header = [m for m in WANT_METRICS if header_counts[m] > 0]
    out_file = "%s-summary.csv" % (run_name)
    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["sample"] + header)
        for i, sample in enumerate(samples):
            info = [metrics[i].get(m, "") for m in header]
            writer.writerow([sample] + info)

def _get_header_counts(samples, metrics):
    header_counts = collections.defaultdict(int)
    for i, _ in enumerate(samples):
        for metric in WANT_METRICS:
            try:
                metrics[i][metric]
                header_counts[metric] += 1
            except KeyError:
                pass
    return header_counts

def _get_sample_name(in_file):
    bam_file = pysam.Samfile(in_file, "rb")
    name = bam_file.header["RG"][0]["SM"]
    bam_file.close()
    return name

def _get_base_bams(work_dir, run_name):
    bam_files = glob.glob(os.path.join(work_dir, "*_%s*-sort.bam" % run_name))
    lane_info = dict()
    for cur_file in bam_files:
        lane_name = os.path.basename(cur_file).split("-")[0]
        lane_parts = lane_name.split("_")
        assert "_".join(lane_parts[1:3]) == run_name
        bc_id = lane_parts[3] if len(lane_parts) == 4 else None
        try:
            bc_id = int(bc_id)
        except ValueError:
            pass
        lane_id = lane_parts[0]
        try:
            lane_id = int(lane_id)
        except ValueError:
            pass
        lane_info[(lane_id, bc_id)] = cur_file
    final = []
    for key in sorted(lane_info.keys()):
        lane, bc = key
        final.append((lane, bc, lane_info[key]))
    return final

def lane_stats(lane, bc_id, run_name, work_dir, parser):
    base_name = "%s_%s" % (lane, run_name)
    if bc_id:
        base_name += "_%s-" % bc_id
    metrics_files = glob.glob(os.path.join(work_dir, "%s*metrics" % base_name))
    metrics = parser.extract_metrics(metrics_files)
    return metrics

if __name__ == "__main__":
    main(*sys.argv[1:])

