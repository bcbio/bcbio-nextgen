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

def _organize_bam_by_lane(align_items):
    out = {}
    for a in align_items:
        out[a["lane"]] = a["out_bam"]
    return out

def organize_samples(dirs, fc_name, fc_date, run_items, align_items):
    """Organize BAM output files by sample name, handling multiplexing.
    """
    bams_by_lane = _organize_bam_by_lane(align_items)
    bams_by_sample = collections.defaultdict(list)
    sample_info = dict()
    fastq_by_sample = collections.defaultdict(list)
    for lane_info in run_items:
        multiplex = lane_info.get("multiplex", None)
        if multiplex:
            mfastq_dir = os.path.join(dirs["work"], "%s_%s_%s_barcode" %
                    (lane_info["lane"], fc_date, fc_name))
            for multi in multiplex:
                name = (lane_info.get("name", ""), lane_info["description"],
                        multi["name"])
                base = "%s_%s_%s_%s" % (lane_info["lane"], fc_date, fc_name, multi["barcode_id"])
                fname = os.path.join(dirs["align"], "%s-sort.bam" % base)
                has_bam = False
                if os.path.exists(fname):
                    has_bam = True
                    bams_by_sample[name].append(fname)
                elif bams_by_lane.has_key(base):
                    has_bam = True
                    bams_by_sample[name].append(bams_by_lane[base])
                else:
                    pass # Not all barcodes may exist; would like a way to check here
                if has_bam:
                    sample_info[name] = lane_info
                    fastq_by_sample[name].append(get_fastq_files(mfastq_dir, lane_info,
                                                                 fc_name, multi["barcode_id"]))
        else:
            name = (lane_info.get("name", ""), lane_info["description"])
            base = "%s_%s_%s" % (lane_info["lane"], fc_date, fc_name)
            fname = os.path.join(dirs["align"], "%s-sort.bam" % base)
            if os.path.exists(fname):
                bams_by_sample[name].append(fname)
            elif bams_by_lane.has_key(base):
                bams_by_sample[name].append(bams_by_lane[base])
            else:
                raise ValueError("Did not find BAM files for %s" % lane_info)
            sample_info[name] = lane_info
            fastq_by_sample[name].append(get_fastq_files(dirs["fastq"],
                                                         lane_info, fc_name))
    return sorted(bams_by_sample.items()), dict(fastq_by_sample), sample_info

def merge_bam_files(bam_files, work_dir, config):
    """Merge multiple BAM files from a sample into a single BAM for processing.
    """
    out_file = os.path.join(work_dir, os.path.basename(bam_files[0]))
    picard = broad.runner_from_config(config)
    picard.run_fn("picard_merge", bam_files, out_file)
    for b in bam_files:
        utils.save_diskspace(b, "BAM merged to %s" % out_file, config)
    return out_file

