"""Top level driver functionality for processing a sequencing lane.
"""
import copy
import os

from bcbio import bam, broad, utils
from bcbio.log import logger
from bcbio.bam import callable, fastq
from bcbio.bam.trim import trim_adapters
from bcbio.pipeline.fastq import get_fastq_files
from bcbio.pipeline.alignment import align_to_sort_bam
from bcbio.pipeline import cleanbam
from bcbio.variation import bedutils, recalibrate

def process_lane(item):
    """Prepare lanes, potentially splitting based on barcodes and reducing the
    number of reads for a test run
    """
    NUM_DOWNSAMPLE = 10000
    logger.debug("Preparing %s" % item["rgnames"]["lane"])
    file1, file2 = get_fastq_files(item)
    if item.get("test_run", False):
        if bam.is_bam(file1):
            file1 = bam.downsample(file1, item, NUM_DOWNSAMPLE)
            file2 = None
        else:
            file1, file2 = fastq.downsample(file1, file2, item,
                                            NUM_DOWNSAMPLE, quick=True)
    item["files"] = [file1, file2]
    return [[item]]

def trim_lane(item):
    """Trim reads with the provided trimming method.
    Support methods: read_through.
    """
    to_trim = [x for x in item["files"] if x is not None]
    dirs = item["dirs"]
    config = item["config"]
    # this block is to maintain legacy configuration files
    trim_reads = config["algorithm"].get("trim_reads", False)
    if not trim_reads:
        logger.info("Skipping trimming of %s." % (", ".join(to_trim)))
        return [[item]]

    if trim_reads == "read_through":
        logger.info("Trimming low quality ends and read through adapter "
                    "sequence from %s." % (", ".join(to_trim)))
        out_files = trim_adapters(to_trim, dirs, config)
    item["files"] = out_files
    return [[item]]

# ## Alignment

def link_bam_file(orig_file, new_dir):
    """Provide symlinks of BAM file and existing indexes.
    """
    new_dir = utils.safe_makedir(new_dir)
    sym_file = os.path.join(new_dir, os.path.basename(orig_file))
    utils.symlink_plus(orig_file, sym_file)
    return sym_file

def _add_supplemental_bams(data):
    """Add supplemental files produced by alignment, useful for structural variant calling.
    """
    file_key = "work_bam"
    if data.get(file_key):
        for supext in ["disc", "sr"]:
            base, ext = os.path.splitext(data[file_key])
            test_file = "%s-%s%s" % (base, supext, ext)
            if os.path.exists(test_file):
                sup_key = file_key + "-plus"
                if not sup_key in data:
                    data[sup_key] = {}
                data[sup_key][supext] = test_file
    return data

def process_alignment(data):
    """Do an alignment of fastq files, preparing a sorted BAM output file.
    """
    if "files" not in data:
        fastq1, fastq2 = None, None
    elif len(data["files"]) == 2:
        fastq1, fastq2 = data["files"]
    else:
        assert len(data["files"]) == 1, data["files"]
        fastq1, fastq2 = data["files"][0], None
    config = data["config"]
    aligner = config["algorithm"].get("aligner", None)
    if fastq1 and os.path.exists(fastq1) and aligner:
        logger.info("Aligning lane %s with %s aligner" % (data["rgnames"]["lane"], aligner))
        data = align_to_sort_bam(fastq1, fastq2, aligner, data)
        data = _add_supplemental_bams(data)
    elif fastq1 and os.path.exists(fastq1) and fastq1.endswith(".bam"):
        sort_method = config["algorithm"].get("bam_sort")
        bamclean = config["algorithm"].get("bam_clean")
        if bamclean is True or bamclean == "picard":
            if sort_method and sort_method != "coordinate":
                raise ValueError("Cannot specify `bam_clean: picard` with `bam_sort` other than coordinate: %s"
                                 % sort_method)
            out_bam = cleanbam.picard_prep(fastq1, data["rgnames"], data["sam_ref"], data["dirs"],
                                           config)
        elif sort_method:
            runner = broad.runner_from_config(config)
            out_file = os.path.join(data["dirs"]["work"], "{}-sort.bam".format(
                os.path.splitext(os.path.basename(fastq1))[0]))
            out_bam = runner.run_fn("picard_sort", fastq1, sort_method, out_file)
        else:
            out_bam = link_bam_file(fastq1, os.path.join(data["dirs"]["work"], "prealign",
                                                         data["rgnames"]["sample"]))
        bam.check_header(out_bam, data["rgnames"], data["sam_ref"], data["config"])
        data["work_bam"] = out_bam
    elif fastq1 is None and "vrn_file" in data:
        data["config"]["algorithm"]["variantcaller"] = ""
        data["work_bam"] = None
    else:
        raise ValueError("Could not process input file: %s" % fastq1)
    return [[data]]

def postprocess_alignment(data):
    """Perform post-processing steps required on full BAM files.
    Prepares list of callable genome regions allowing subsequent parallelization.
    Cleans input BED files to avoid issues with overlapping input segments.
    """
    data = bedutils.clean_inputs(data)
    if data["work_bam"]:
        callable_region_bed, nblock_bed, callable_bed = \
            callable.block_regions(data["work_bam"], data["sam_ref"], data["config"])
        data["regions"] = {"nblock": nblock_bed, "callable": callable_bed}
        if (os.path.exists(callable_region_bed) and
                not data["config"]["algorithm"].get("variant_regions")):
            data["config"]["algorithm"]["variant_regions"] = callable_region_bed
            data = bedutils.clean_inputs(data)
        data = _recal_no_markduplicates(data)
    return [data]

def _recal_no_markduplicates(data):
    orig_config = copy.deepcopy(data["config"])
    data["config"]["algorithm"]["mark_duplicates"] = False
    data = recalibrate.prep_recal(data)[0][0]
    data["config"] = orig_config
    return data
