"""Top level driver functionality for processing a sequencing lane.
"""
import contextlib
import copy
import itertools
import os

import pysam

from bcbio import utils, broad
from bcbio.log import logger
from bcbio.bam import callable
from bcbio.bam.trim import brun_trim_fastq, trim_read_through
from bcbio.pipeline.fastq import get_fastq_files, needs_fastq_conversion
from bcbio.pipeline.alignment import align_to_sort_bam, get_genome_ref
from bcbio.pipeline import cleanbam
from bcbio.ngsalign.split import split_read_files
from bcbio.variation import recalibrate, vcfutils

def _item_needs_compute(lane_items):
    """Determine if any item needs computing resources to spin up a cluster.
    """
    for item in lane_items:
        config = item["config"]
        split_size = config.get("distributed", {}).get("align_split_size",
                                                       config["algorithm"].get("align_split_size", None))
        if split_size is not None:
            return True
        if needs_fastq_conversion(item, item["config"]):
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
        return run_parallel("process_lane", [[x] for x in lanes])
    else:
        return [process_lane(x)[0] for x in lanes]

def process_lane(item):
    """Prepare lanes, potentially splitting based on barcodes.
    """
    logger.debug("Preparing %s" % item["rgnames"]["lane"])
    item["files"] = get_fastq_files(item)
    return [item]

def trim_lane(item):
    """
    if trim_reads is set with no trimmer specified, default to B-run trimming
    only. if trimmer is set to a supported type, perform that trimming
    instead.

    """
    to_trim = [x for x in item["files"] if x is not None]
    dirs = item["dirs"]
    config = item["config"]
    # this block is to maintain legacy configuration files
    trim_reads = config["algorithm"].get("trim_reads", False)
    if not trim_reads:
        logger.info("Skipping trimming of %s." % (", ".join(to_trim)))
        return item

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
    item["files"] = out_files
    return [[item]]

# ## Alignment

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

def _check_prealigned_bam(in_bam, ref_file, config):
    """Ensure a pre-aligned BAM file matches the expected reference genome.
    """
    ref_contigs = [c["SN"] for c in vcfutils.ref_file_contigs(ref_file, config)]
    with contextlib.closing(pysam.Samfile(in_bam, "rb")) as bamfile:
        bam_contigs = [c["SN"] for c in bamfile.header["SQ"]]
    problems = []
    warnings = []
    for bc, rc in itertools.izip_longest(bam_contigs, ref_contigs):
        if bc != rc:
            if bc and rc:
                problems.append("Reference mismatch. BAM: %s Reference: %s" (bc, rc))
            elif bc:
                problems.append("Extra BAM chromosomes: %s" % bc)
            elif rc:
                warnings.append("Extra reference chromosomes: %s" % rc)
    if problems:
        raise ValueError("Unexpected order, name or contig mismatches between input BAM and reference file:\n%s\n"
                         % "\n".join(problems))
    if warnings:
        print("*** Potential problems in input BAM compared to reference:\n%s\n" %
              "\n".join(warnings))

def process_alignment(data):
    """Do an alignment of fastq files, preparing a sorted BAM output file.
    """
    fastq1, fastq2 = data["files"]
    config = data["config"]
    aligner = config["algorithm"].get("aligner", None)
    out_bam = ""
    if os.path.exists(fastq1) and aligner:
        logger.info("Aligning lane %s with %s aligner" % (data["rgnames"]["lane"], aligner))
        out_bam, ref_file = align_to_sort_bam(fastq1, fastq2, data["rgnames"],
                                              data["genome_build"], aligner,
                                              data["dirs"], data["config"])
    elif os.path.exists(fastq1) and fastq1.endswith(".bam"):
        sort_method = config["algorithm"].get("bam_sort")
        bamclean = config["algorithm"].get("bam_clean")
        if sort_method:
            runner = broad.runner_from_config(config)
            out_file = os.path.join(data["dirs"]["work"], "{}-sort.bam".format(
                os.path.splitext(os.path.basename(fastq1))[0]))
            out_bam = runner.run_fn("picard_sort", fastq1, sort_method, out_file)
        elif bamclean is True or bamclean == "picard":
            out_bam = cleanbam.picard_prep(fastq1, data["rgnames"], data["sam_ref"], data["dirs"], config)
        else:
            out_bam = link_bam_file(fastq1, os.path.join(data["dirs"]["work"], "prealign",
                                                         data["rgnames"]["sample"]))
        _check_prealigned_bam(fastq1, data["sam_ref"], config)
    if not out_bam and not os.path.exists(fastq1):
        raise ValueError("Could not find input file: %s" % fastq1)
    data["sam_ref"] = get_genome_ref(data["genome_build"], None, data["dirs"]["galaxy"])[-1]
    data["work_bam"] = out_bam
    return [[data]]

def align_prep_full(data, config_file):
    """Perform alignment and post-processing required on full BAM files.
    Prepare list of callable genome regions allowing subsequent parallelization.
    """
    _, ref_file = get_genome_ref(data["genome_build"], None, data["dirs"]["galaxy"])
    data["sam_ref"] = ref_file
    if data["files"][0] is None and "vrn_file" in data:
        data["config"]["algorithm"]["variantcaller"] = ""
        data["work_bam"] = None
    else:
        data = process_alignment(data)[0][0]
        callable_region_bed, nblock_bed = callable.block_regions(data["work_bam"],
                                                                 data["sam_ref"], data["config"])
        data["regions"] = {"nblock": nblock_bed}
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
