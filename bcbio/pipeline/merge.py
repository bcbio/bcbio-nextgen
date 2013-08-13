"""Handle multiple samples present on a single flowcell

Merges samples located in multiple lanes on a flowcell. Unique sample names identify
items to combine within a group.
"""
import copy
import collections
import os
import shutil

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.provenance import do

def combine_fastq_files(in_files, work_dir, config):
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
        for f1, f2 in in_files:
            utils.save_diskspace(f1, "fastq merged to %s" % out1, config)
            if f2:
                utils.save_diskspace(f2, "fastq merged to %s" % out2, config)
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
        bam_files = [x["work_bam"] for x in item_group]
        item_group.sort(key=_sort_by_lane_barcode)

        out.append({"name": name, "info": item_group[0]["info"],
                    "fastq_files": fastq_files, "bam_files": bam_files,
                    "dirs": copy.deepcopy(dirs), "config": item_group[0]["config"],
                    "config_file": config_file})
    out.sort(key=_sort_by_lane_barcode)
    out = [[x] for x in out]
    return out

def merge_bam_files(bam_files, work_dir, config, out_file=None):
    """Merge multiple BAM files from a sample into a single BAM for processing.

    Uses bamtools for merging, which handles large numbers of input BAMs.
    """
    if len(bam_files) == 1:
        return bam_files[0]
    else:
        if out_file is None:
            out_file = os.path.join(work_dir, os.path.basename(sorted(bam_files)[0]))
        if not utils.file_exists(out_file) or not utils.file_exists(out_file + ".bai"):
            bamtools = config_utils.get_program("bamtools", config)
            samtools = config_utils.get_program("samtools", config)
            resources = config_utils.get_resources("samtools", config)
            # Could use multiple cores if we build multicore steps in here
            #num_cores = config["algorithm"].get("num_cores", 1)
            num_cores = 1
            max_mem = resources.get("memory", "2048M")
            with file_transaction(out_file) as tx_out_file:
                tx_out_prefix = os.path.splitext(tx_out_file)[0]
                with utils.tmpfile(dir=work_dir, prefix="bammergelist") as bam_file_list:
                    with open(bam_file_list, "w") as out_handle:
                        for f in sorted(bam_files):
                            out_handle.write("%s\n" % f)
                    cmd = ("{bamtools} merge -list {bam_file_list} | "
                           "{samtools} sort -@ {num_cores} -m {max_mem} - {tx_out_prefix}")
                    do.run(cmd.format(**locals()), "Merge bam files", None)
            for b in bam_files:
                utils.save_diskspace(b, "BAM merged to %s" % out_file, config)
        picard = broad.runner_from_config(config)
        picard.run_fn("picard_index", out_file)
        return out_file
