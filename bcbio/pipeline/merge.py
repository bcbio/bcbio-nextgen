"""Handle multiple samples present on a single flowcell

Merges samples located in multiple lanes on a flowcell. Unique sample names identify
items to combine within a group.
"""
import os
import shutil
import collections

from bcbio import utils, broad
from bcbio.pipeline.fastq import get_fastq_files

def combine_fastq_files(in_files, work_dir):
    if len(in_files) == 1:
        return in_files[0]
    else:
        cur1, cur2 = in_files[0]
        out1 = os.path.join(work_dir, os.path.basename(cur1))
        out2 = os.path.join(work_dir, os.path.basename(cur2)) if cur2 else None
        if not os.path.exists(out1):
            with open(out1, "a") as out_handle:
                for (cur1, _) in in_files:
                    with open(cur1) as in_handle:
                        shutil.copyfileobj(in_handle, out_handle)
        if out2 and not os.path.exists(out2):
            with open(out2, "a") as out_handle:
                for (_, cur2) in in_files:
                    with open(cur2) as in_handle:
                        shutil.copyfileobj(in_handle, out_handle)
        return out1, out2

def organize_samples(items, dirs, config_file):
    """Organize BAM output files by sample name.
    """
    def _sort_by_lane_barcode(x):
        """Index a sample by lane and barcode.
        """
        return (x["info"]["lane"], x["info"]["barcode_id"])
    items_by_name = collections.defaultdict(list)
    for item in items:
        name = (item["info"].get("name", ""), item["info"]["description"])
        items_by_name[name].append(item)
    out = []
    for name, item_group in items_by_name.iteritems():
        fastq_files = [x["fastq"] for x in item_group]
        bam_files = [x["out_bam"] for x in item_group]
        item_group.sort(key=_sort_by_lane_barcode)

        out.append({"name": name, "info": item_group[0]["info"],
                    "fastq_files": fastq_files, "bam_files": bam_files,
                    "dirs": dirs, "config": item_group[0]["config"],
                    "config_file": config_file})
    out.sort(key=_sort_by_lane_barcode)
    out = [[x] for x in out]
    return out

def merge_bam_files(bam_files, work_dir, config):
    """Merge multiple BAM files from a sample into a single BAM for processing.
    """
    out_file = os.path.join(work_dir, os.path.basename(bam_files[0]))
    picard = broad.runner_from_config(config)
    picard.run_fn("picard_merge", bam_files, out_file)
    for b in bam_files:
        utils.save_diskspace(b, "BAM merged to %s" % out_file, config)
    return out_file

