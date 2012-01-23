"""Top level driver functionality for processing a sequencing lane.
"""
import os
import copy

from bcbio.log import logger
from bcbio.pipeline.fastq import get_fastq_files
from bcbio.pipeline.demultiplex import split_by_barcode
from bcbio.pipeline.alignment import align_to_sort_bam
from bcbio.bam.trim import brun_trim_fastq

def process_lane(lane_items, fc_name, fc_date, dirs, config):
    """Prepare lanes, potentially splitting based on barcodes.
    """
    lane_name = "%s_%s_%s" % (lane_items[0]['lane'], fc_date, fc_name)
    logger.info("Demulitplexing %s" % lane_name)
    full_fastq1, full_fastq2 = get_fastq_files(dirs["fastq"], dirs["work"],
                                               lane_items[0], fc_name, config=config)
    bc_files = split_by_barcode(full_fastq1, full_fastq2, lane_items,
                                lane_name, dirs, config)
    out = []
    for item in lane_items:
        config = _update_config_w_custom(config, item)
        # Can specify all barcodes but might not have actual sequences
        # Would be nice to have a good way to check this is okay here.
        if bc_files.has_key(item["barcode_id"]):
            fastq1, fastq2 = bc_files[item["barcode_id"]]
            cur_lane_name = lane_name
            cur_lane_desc = item["description"]
            if item.get("name", "") and config["algorithm"].get("include_short_name", True):
                cur_lane_desc = "%s : %s" % (item["name"], cur_lane_desc)
            if item["barcode_id"] is not None:
                cur_lane_name += "_%s" % (item["barcode_id"])
            if config["algorithm"].get("trim_reads", False):
                trim_info = brun_trim_fastq([x for x in [fastq1, fastq2] if x is not None],
                                            dirs, config)
                fastq1 = trim_info[0]
                if fastq2 is not None:
                    fastq2 = trim_info[1]
            out.append((fastq1, fastq2, item, cur_lane_name, cur_lane_desc,
                        dirs, config))
    return out

def process_alignment(fastq1, fastq2, info, lane_name, lane_desc,
                      dirs, config):
    """Do an alignment of fastq files, preparing a sorted BAM output file.
    """
    aligner = config["algorithm"].get("aligner", None)
    out_bam = ""
    if os.path.exists(fastq1) and aligner:
        logger.info("Aligning lane %s with %s aligner" % (lane_name, aligner))
        out_bam = align_to_sort_bam(fastq1, fastq2, info["genome_build"], aligner,
                                    lane_name, lane_desc, dirs, config)
    return [{"fastq": [fastq1, fastq2], "out_bam": out_bam, "info": info,
             "config": config}]

def _update_config_w_custom(config, lane_info):
    """Update the configuration for this lane if a custom analysis is specified.
    """
    config = copy.deepcopy(config)
    analysis_type = lane_info.get("analysis", "")
    custom = config["custom_algorithms"].get(analysis_type, None)
    if custom:
        for key, val in custom.iteritems():
            config["algorithm"][key] = val
    # apply any algorithm details specified with the lane
    for key, val in lane_info.get("algorithm", {}).iteritems():
        config["algorithm"][key] = val
    return config

