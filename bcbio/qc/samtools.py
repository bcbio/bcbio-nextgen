"""Quality control metrics from samtools.
"""
import os

import toolz as tz

from bcbio.distributed.transaction import file_transaction
from bcbio import utils
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do

def run(_, data, out_dir=None):
    """Run samtools stats with reports on mapped reads, duplicates and insert sizes.
    """
    stats_file, idxstats_file = _get_stats_files(data, out_dir)
    samtools = config_utils.get_program("samtools", data["config"])
    bam_file = dd.get_align_bam(data) or dd.get_work_bam(data)
    if not utils.file_exists(stats_file):
        utils.safe_makedir(out_dir)
        with file_transaction(data, stats_file) as tx_out_file:
            cores = dd.get_num_cores(data)
            cmd = "{samtools} stats -@ {cores} {bam_file}"
            cmd += " > {tx_out_file}"
            do.run(cmd.format(**locals()), "samtools stats", data)
    if not utils.file_exists(idxstats_file):
        utils.safe_makedir(out_dir)
        with file_transaction(data, idxstats_file) as tx_out_file:
            cmd = "{samtools} idxstats {bam_file}"
            cmd += " > {tx_out_file}"
            do.run(cmd.format(**locals()), "samtools index stats", data)
    out = {"base": idxstats_file, "secondary": [stats_file]}
    out["metrics"] = _parse_samtools_stats(stats_file)
    return out

def run_and_save(data):
    """Run QC, saving file outputs in data dictionary.
    """
    run(None, data)
    stats_file, idxstats_file = _get_stats_files(data)
    data = tz.update_in(data, ["depth", "samtools", "stats"], lambda x: stats_file)
    data = tz.update_in(data, ["depth", "samtools", "idxstats"], lambda x: idxstats_file)
    return data

def _get_stats_files(data, out_dir=None):
    """Retrieve stats files from pre-existing dictionary or filesystem.
    """
    if not out_dir:
        out_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data),
                                                  "qc", dd.get_sample_name(data), "samtools"))
    stats_file = tz.get_in(["depth", "samtools", "stats"], data)
    idxstats_file = tz.get_in(["depth", "samtools", "idxstats"], data)
    if not stats_file:
        stats_file = os.path.join(out_dir, "%s.txt" % dd.get_sample_name(data))
    if not idxstats_file:
        idxstats_file = os.path.join(out_dir, "%s-idxstats.txt" % dd.get_sample_name(data))
    return stats_file, idxstats_file

def _parse_samtools_stats(stats_file):
    out = {}
    want = {"raw total sequences": "Total_reads",
            "reads mapped": "Mapped_reads",
            "reads mapped and paired": "Mapped_paired_reads",
            "reads duplicated": "Duplicates",
            "insert size average": "Average_insert_size",
            "average length": "Average_read_length",
            }
    with open(stats_file) as in_handle:
        for line in in_handle:
            if not line.startswith("SN"):
                continue
            parts = line.split("\t")
            metric, stat_str = parts[1:3]
            metric = metric.replace(":", "").strip()
            if metric in want:
                stat = float(stat_str.strip())
                out[want[metric]] = stat
    # Ensure we have zero values for any metrics not present in stats output
    for metric in want.values():
        if metric not in out:
            out[metric] = 0
    return out
