"""Top level driver functionality for processing a sequencing lane.
"""
import copy
import os
import re

from bcbio import utils, broad
from bcbio.log import logger
from bcbio.bam import callable
from bcbio.bam.trim import brun_trim_fastq, trim_read_through
from bcbio.pipeline.fastq import get_fastq_files, needs_fastq_conversion
from bcbio.pipeline.demultiplex import split_by_barcode
from bcbio.pipeline.alignment import align_to_sort_bam, get_genome_ref
from bcbio.pipeline import cleanbam, merge, shared, sample
from bcbio.ngsalign.split import split_read_files
from bcbio.variation import recalibrate

def _item_needs_compute(lanes):
    """Determine if any item needs computing resources to spin up a cluster.
    """
    for lane_items, _, _, _, config in lanes:
        # check if multiplexed
        if len(lane_items) > 1 or lane_items[0]["barcode_id"] is not None:
            return True
        # check if we need to process the input by splitting or conversion
        item = lane_items[0]
        config = shared.update_config_w_custom(config, item)
        split_size = config.get("distributed", {}).get("align_split_size",
                                                       config["algorithm"].get("align_split_size", None))
        if split_size is not None:
            return True
        if needs_fastq_conversion(item, config):
            return True
    return False

def _prep_fastq_files(item, bc_files, dirs, config):
    """Potentially prepare input FASTQ files for processing.
    """
    fastq1, fastq2 = bc_files[item["barcode_id"]]
    split_size = config.get("distributed", {}).get("align_split_size",
                                                   config["algorithm"].get("align_split_size", None))
    if split_size:
        split_dir = utils.safe_makedir(os.path.join(dirs["work"], "align_splitprep", item["description"]))
        return split_read_files(fastq1, fastq2, item, split_size, split_dir, dirs, config)
    else:
        return [[fastq1, fastq2, None]]

def process_all_lanes(lanes, run_parallel):
    """Process all input lanes, avoiding starting a cluster if not needed.
    """
    lanes = list(lanes)
    if _item_needs_compute(lanes):
        return run_parallel("process_lane", lanes)
    else:
        return [apply(process_lane, xs)[0] for xs in lanes]

def process_lane(lane_items, fc_name, fc_date, dirs, config):
    """Prepare lanes, potentially splitting based on barcodes.
    """
    lane_name = "%s_%s_%s" % (lane_items[0]['lane'], fc_date, fc_name)
    logger.info("Preparing %s" % lane_name)
    full_fastq1, full_fastq2 = get_fastq_files(dirs["fastq"],
                                               dirs["work"], lane_items[0], fc_name, dirs=dirs,
                                               config=shared.update_config_w_custom(config, lane_items[0]))
    bc_files = split_by_barcode(full_fastq1, full_fastq2, lane_items,
                                lane_name, dirs, config)
    out = []
    for item in lane_items:
        config = shared.update_config_w_custom(config, item)
        # Can specify all barcodes but might not have actual sequences
        # Would be nice to have a good way to check this is okay here.
        if item["barcode_id"] in bc_files:
            for fastq1, fastq2, lane_ext in _prep_fastq_files(item, bc_files, dirs, config):
                cur_lane_name = lane_name
                cur_lane_desc = item["description"]
                if item.get("name", "") and config["algorithm"].get("include_short_name", True):
                    cur_lane_desc = "%s : %s" % (item["name"], cur_lane_desc)
                if item["barcode_id"] is not None:
                    cur_lane_name += "_%s" % (item["barcode_id"])
                if lane_ext is not None:
                    cur_lane_name += "_s{0}".format(lane_ext)
                out.append((fastq1, fastq2, item, cur_lane_name, cur_lane_desc,
                            dirs, config))
    return out


def trim_lane(fastq1, fastq2, info, lane_name, lane_desc, dirs, config):
    """
    if trim_reads is set with no trimmer specified, default to B-run trimming
    only. if trimmer is set to a supported type, perform that trimming
    instead.

    """
    to_trim = [x for x in [fastq1, fastq2] if x is not None]
    # this block is to maintain legacy configuration files
    trim_reads = config["algorithm"].get("trim_reads", False)
    if not trim_reads:
        logger.info("Skipping trimming of %s." % (", ".join(to_trim)))
        return [(fastq1, fastq2, info, lane_name, lane_desc, dirs, config)]

    # swap the default to None if trim_reads gets deprecated



    if trim_reads == "low_quality" or trim_reads == "true":
        logger.info("Trimming low quality ends from %s."
                    % (", ".join(to_trim)))
        out_files = brun_trim_fastq(to_trim, dirs, config)

    if trim_reads == "read_through":
        logger.info("Trimming low quality ends and read through adapter "
                    "sequence from %s." % (", ".join(to_trim)))
        out_files = trim_read_through(to_trim, dirs, config)

    else:
        logger.info("Trimming low quality ends from %s."
                    % (", ".join(to_trim)))
        out_files = brun_trim_fastq(to_trim, dirs, config)

    fastq1 = out_files[0]
    if fastq2 is not None:
        fastq2 = out_files[1]

    return [(fastq1, fastq2, info, lane_name, lane_desc, dirs, config)]

# ## Alignment

def rg_names(lane_name, lane_desc, config):
    return {"rg": lane_name.split("_")[0],
            "sample": lane_desc,
            "lane": lane_name,
            "pl": config["algorithm"]["platform"].lower(),
            "pu": (lane_name.rsplit("_", 1)[0]
                   if re.search(r"_s\d+$", lane_name) is not None
                   else lane_name)}

def link_bam_file(orig_file, new_dir):
    """Provide symlinks of BAM file and existing indexes.
    """
    new_dir = utils.safe_makedir(new_dir)
    update_files = []
    for fname in (orig_file, "%s.bai" % orig_file,
                  "%s.bai" % os.path.splitext(orig_file)[0]):
        if utils.file_exists(fname):
            sym_file = os.path.join(new_dir, os.path.basename(fname))
            if not os.path.exists(sym_file):
                os.symlink(fname, sym_file)
            update_files.append(sym_file)
    return update_files[0]

def process_alignment(fastq1, fastq2, info, lane_name, lane_desc,
                      dirs, config):
    """Do an alignment of fastq files, preparing a sorted BAM output file.
    """
    aligner = config["algorithm"].get("aligner", None)
    out_bam = ""
    names = rg_names(lane_name, lane_desc, config)
    _, ref_file = get_genome_ref(info["genome_build"], None, dirs["galaxy"])
    if os.path.exists(fastq1) and aligner:
        logger.info("Aligning lane %s with %s aligner" % (lane_name, aligner))
        out_bam, ref_file = align_to_sort_bam(fastq1, fastq2, names,
                                              info["genome_build"], aligner,
                                              dirs, config)
    elif os.path.exists(fastq1) and fastq1.endswith(".bam"):
        sort_method = config["algorithm"].get("bam_sort")
        bamclean = config["algorithm"].get("bam_clean")
        if sort_method:
            runner = broad.runner_from_config(config)
            out_file = os.path.join(dirs["work"], "{}-sort.bam".format(
                os.path.splitext(os.path.basename(fastq1))[0]))
            out_bam = runner.run_fn("picard_sort", fastq1, sort_method, out_file)
        elif bamclean is True or bamclean == "picard":
            out_bam = cleanbam.picard_prep(fastq1, names, ref_file, dirs, config)
        else:
            out_bam = link_bam_file(fastq1, os.path.join(dirs["work"], "prealign",
                                                         names["sample"]))
    if not out_bam and not os.path.exists(fastq1):
        raise ValueError("Could not find input file: %s" % fastq1)
    return [{"fastq": [fastq1, fastq2], "work_bam": out_bam, "info": info,
             "sam_ref": ref_file, "config": config}]

def align_prep_full(fastq1, fastq2, info, lane_name, lane_desc,
                    dirs, config, config_file):
    """Perform alignment and post-processing required on full BAM files.
    Prepare list of callable genome regions allowing subsequent parallelization.
    """
    if fastq1 is None and "vrn_file" in info:
        _, ref_file = get_genome_ref(info["genome_build"], None, dirs["galaxy"])
        config["algorithm"]["variantcaller"] = ""
        data = {"info": info, "sam_ref": ref_file,
                "work_bam": None,
                "genome_build": info["genome_build"],
                "name": ("", lane_desc),
                "vrn_file": info["vrn_file"],
                "dirs": copy.deepcopy(dirs), "config": config}
    else:
        align_out = process_alignment(fastq1, fastq2, info, lane_name, lane_desc,
                                      dirs, config)[0]
        data = _organize_merge_samples(align_out, dirs, config_file)
        callable_region_bed, analysis_regions = callable.block_regions(data["work_bam"],
                                                                       data["sam_ref"], config)
        data["regions"] = analysis_regions
        if (os.path.exists(callable_region_bed) and
                not data["config"]["algorithm"].get("variant_regions")):
            data["config"]["algorithm"]["variant_regions"] = callable_region_bed
        data["callable_bam"] = data["work_bam"]
        data = _recal_no_markduplicates(data)
    return [data]

def _recal_no_markduplicates(data):
    orig_config = copy.deepcopy(data["config"])
    data["config"]["algorithm"]["mark_duplicates"] = False
    data = recalibrate.prep_recal(data)[0][0]
    data["config"] = orig_config
    return data

def _organize_merge_samples(align_out, dirs, config_file):
    """Back compatibility handling organizing and merging samples.
    """
    samples = merge.organize_samples([align_out], dirs, config_file)
    assert len(samples) == 1 and len(samples[0]) == 1
    sample_data = samples[0][0]
    sample_data["dirs"]["work"] = os.path.dirname(align_out["work_bam"])
    samples2 = sample.merge_sample(sample_data)
    assert len(samples2) == 1 and len(samples2[0]) == 1
    data = samples2[0][0]
    data["dirs"] = copy.deepcopy(dirs)
    return data
