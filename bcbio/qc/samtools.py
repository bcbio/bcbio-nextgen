"""Quality control metrics from samtools.
"""
import os

import yaml

from bcbio.distributed.transaction import file_transaction
from bcbio import utils
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do

def run(bam_file, data, out_dir):
    """Run samtools stats with reports on mapped reads, duplicates and insert sizes.
    """
    stats_file = os.path.join(out_dir, "%s.txt" % dd.get_sample_name(data))
    if not utils.file_exists(stats_file):
        utils.safe_makedir(out_dir)
        samtools = config_utils.get_program("samtools", data["config"])
        with file_transaction(data, stats_file) as tx_out_file:
            cmd = "{samtools} stats {bam_file}"
            cmd += " > {tx_out_file}"
            do.run(cmd.format(**locals()), "samtools stats", data)
    out = _parse_samtools_stats(stats_file)
    out.update(_parse_offtargets(bam_file))
    return out

def _parse_samtools_stats(stats_file):
    out = {}
    want = {"raw total sequences": "Total reads", "reads mapped": "Mapped reads",
            "reads duplicated": "Duplicates", "insert size average": "Average insert size"}
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
                if metric in ["reads mapped", "reads duplicated"]:
                    out["%s pct" % want[metric]] = stat / out["Total reads"]
    return out

def _parse_offtargets(bam_file):
    """
    Add to metrics off-targets reads if it exitst
    """
    off_target = bam_file.replace(".bam", "-offtarget-stats.yaml")
    if os.path.exists(off_target):
        res = yaml.load(open(off_target))
        res['offtarget_pct'] = "%.3f" % (float(res['offtarget']) / float(res['mapped']))
        return res
    return {}
