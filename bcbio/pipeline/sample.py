"""High level entry point for processing a sample.

Samples may include multiple lanes, or barcoded subsections of lanes,
processed together.
"""
import collections
import copy
import glob
import os
import re

import toolz as tz

from bcbio import utils, bam, broad
from bcbio.log import logger
from bcbio.distributed import objectstore
from bcbio.pipeline.merge import merge_bam_files
from bcbio.bam import callable, highdepth, trim
from bcbio.hla import optitype
from bcbio.ngsalign import postalign
from bcbio.pipeline.fastq import get_fastq_files
from bcbio.pipeline.alignment import align_to_sort_bam
from bcbio.pipeline import cleanbam
from bcbio.variation import bedutils, coverage, recalibrate
from bcbio.variation import multi as vmulti
import bcbio.pipeline.datadict as dd
from bcbio.pipeline.fastq import merge as fq_merge
from bcbio.bam import merge as bam_merge
from bcbio.bam import skewer

def prepare_sample(data):
    """Prepare a sample to be run, potentially converting from BAM to
    FASTQ and/or downsampling the number of reads for a test run
    """
    logger.debug("Preparing %s" % data["rgnames"]["sample"])
    file1, file2 = get_fastq_files(data)
    data["files"] = [file1, file2]
    return [[data]]

def trim_sample(data):
    """Trim from a sample with the provided trimming method.
    Support methods: read_through.
    """
    config = data["config"]
    # this block is to maintain legacy configuration files
    trim_reads = config["algorithm"].get("trim_reads", False)
    sample_name = dd.get_sample_name(data)
    if not trim_reads:
        logger.info("Skipping trimming of %s." % sample_name)
        return [[data]]

    if trim_reads == "read_through":
        if "skewer" in dd.get_tools_on(data):
            trim_adapters = skewer.trim_adapters
        else:
            trim_adapters = trim.trim_adapters
        out_files = trim_adapters(data)
        data["files"] = out_files
    return [[data]]

# ## Alignment

def link_bam_file(orig_file, new_dir):
    """Provide symlinks of BAM file and existing indexes.
    """
    new_dir = utils.safe_makedir(new_dir)
    sym_file = os.path.join(new_dir, os.path.basename(orig_file))
    utils.symlink_plus(orig_file, sym_file)
    return sym_file

def _add_supplemental_bams(data):
    """Add supplemental files produced by alignment, useful for structural
    variant calling.
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

def _add_hla_files(data):
    """Add extracted fastq files of HLA alleles for typing.
    """
    if "hla" not in data:
        data["hla"] = {}
    align_file = dd.get_align_bam(data)
    hla_dir = os.path.join(os.path.dirname(align_file), "hla")
    if not os.path.exists(hla_dir):
        hla_dir = None
    data["hla"]["fastq"] = hla_dir
    return data

def process_alignment(data, alt_input=None):
    """Do an alignment of fastq files, preparing a sorted BAM output file.
    """
    data = utils.to_single_data(data)
    fastq1, fastq2 = dd.get_input_sequence_files(data)
    if alt_input:
        fastq1, fastq2 = alt_input
    config = data["config"]
    aligner = config["algorithm"].get("aligner", None)
    if fastq1 and objectstore.file_exists_or_remote(fastq1) and aligner:
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
                                           data)
        elif sort_method:
            runner = broad.runner_from_path("picard", config)
            out_file = os.path.join(data["dirs"]["work"], "{}-sort.bam".format(
                os.path.splitext(os.path.basename(fastq1))[0]))
            out_bam = runner.run_fn("picard_sort", fastq1, sort_method, out_file)
        else:
            out_bam = link_bam_file(fastq1, os.path.join(data["dirs"]["work"], "prealign",
                                                         data["rgnames"]["sample"]))
        bam.index(out_bam, data["config"])
        bam.check_header(out_bam, data["rgnames"], data["sam_ref"], data["config"])
        dedup_bam = postalign.dedup_bam(out_bam, data)
        bam.index(dedup_bam, data["config"])
        data["work_bam"] = dedup_bam
    elif fastq1 and objectstore.file_exists_or_remote(fastq1) and fastq1.endswith(".cram"):
        data["work_bam"] = fastq1
    elif fastq1 is None and "vrn_file" in data:
        data["config"]["algorithm"]["variantcaller"] = False
        data["work_bam"] = None
    elif not fastq1:
        raise ValueError("No 'files' specified for input sample: %s" % dd.get_sample_name(data))
    else:
        raise ValueError("Could not process input file from sample configuration. \n" +
                         fastq1 +
                         "\nIs the path to the file correct or is empty?\n" +
                         "If it is a fastq file (not pre-aligned BAM or CRAM), "
                         "is an aligner specified in the input configuration?")
    if data.get("work_bam"):
        # Add stable 'align_bam' target to use for retrieving raw alignment
        data["align_bam"] = data["work_bam"]
        data = _add_hla_files(data)
    return [[data]]

def prep_samples(*items):
    """Handle any global preparatory steps for samples with potentially shared data.

    Avoids race conditions in postprocess alignment when performing prep tasks
    on shared files between multiple similar samples.

    Cleans input BED files to avoid issues with overlapping input segments.
    """
    out = []
    for data in (utils.to_single_data(x) for x in items):
        data = bedutils.clean_inputs(data)
        out.append([data])
    return out

def postprocess_alignment(data):
    """Perform post-processing steps required on full BAM files.
    Prepares list of callable genome regions allowing subsequent parallelization.
    """
    data = utils.to_single_data(data)
    bam_file = data.get("align_bam") or data.get("work_bam")
    if vmulti.bam_needs_processing(data) and bam_file and bam_file.endswith(".bam"):
        ref_file = dd.get_ref_file(data)
        out_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "align",
                                                  dd.get_sample_name(data)))
        bam_file_ready = os.path.join(out_dir, os.path.basename(bam_file))
        if not utils.file_exists(bam_file_ready):
            utils.symlink_plus(bam_file, bam_file_ready)
        bam.index(bam_file_ready, data["config"])
        callable_region_bed, nblock_bed, callable_bed = \
            callable.block_regions(bam_file_ready, ref_file, data)
        highdepth_bed = highdepth.identify(data)
        sample_callable = callable.sample_callable_bed(bam_file_ready, ref_file, data)
        offtarget_stats = callable.calculate_offtarget(bam_file_ready, ref_file, data)
        data["regions"] = {"nblock": nblock_bed, "callable": callable_bed, "highdepth": highdepth_bed,
                           "sample_callable": sample_callable,
                           "offtarget_stats": offtarget_stats}
        data = coverage.assign_interval(data)
        if (os.path.exists(callable_region_bed) and
                not data["config"]["algorithm"].get("variant_regions")):
            data["config"]["algorithm"]["variant_regions"] = callable_region_bed
            data = bedutils.clean_inputs(data)
        data = _recal_no_markduplicates(data)
    return [[data]]

def _recal_no_markduplicates(data):
    orig_config = copy.deepcopy(data["config"])
    data["config"]["algorithm"]["mark_duplicates"] = False
    data = recalibrate.prep_recal(data)[0][0]
    data["config"] = orig_config
    return data

def _merge_out_from_infiles(in_files):
    """Generate output merged file name from set of input files.

    Handles non-shared filesystems where we don't know output path when setting
    up split parts.
    """
    fname = os.path.commonprefix([os.path.basename(f) for f in in_files])
    while fname.endswith(("-", "_", ".")):
        fname = fname[:-1]
    ext = os.path.splitext(in_files[0])[-1]
    dirname = os.path.dirname(in_files[0])
    while dirname.endswith(("split", "merge")):
        dirname = os.path.dirname(dirname)
    return os.path.join(dirname, "%s%s" % (fname, ext))

def delayed_bam_merge(data):
    """Perform a merge on previously prepped files, delayed in processing.

    Handles merging of associated split read and discordant files if present.
    """
    if data.get("combine"):
        assert len(data["combine"].keys()) == 1
        file_key = data["combine"].keys()[0]
        extras = []
        for x in data["combine"][file_key].get("extras", []):
            if isinstance(x, (list, tuple)):
                extras.extend(x)
            else:
                extras.append(x)
        if file_key in data:
            extras.append(data[file_key])
        in_files = sorted(list(set(extras)))
        out_file = tz.get_in(["combine", file_key, "out"], data, _merge_out_from_infiles(in_files))
        sup_exts = data.get(file_key + "-plus", {}).keys()
        for ext in sup_exts + [""]:
            merged_file = None
            if os.path.exists(utils.append_stem(out_file, "-" + ext)):
                cur_out_file, cur_in_files = out_file, []
            if ext:
                cur_in_files = list(filter(os.path.exists, (utils.append_stem(f, "-" + ext) for f in in_files)))
                cur_out_file = utils.append_stem(out_file, "-" + ext) if len(cur_in_files) > 0 else None
            else:
                cur_in_files, cur_out_file = in_files, out_file
            if cur_out_file:
                config = copy.deepcopy(data["config"])
                config["algorithm"]["save_diskspace"] = False
                if len(cur_in_files) > 0:
                    merged_file = merge_bam_files(cur_in_files, os.path.dirname(cur_out_file), config,
                                                  out_file=cur_out_file)
                else:
                    assert os.path.exists(cur_out_file)
                    merged_file = cur_out_file
            if merged_file:
                if ext:
                    data[file_key + "-plus"][ext] = merged_file
                else:
                    data[file_key] = merged_file
        data.pop("region", None)
        data.pop("combine", None)
    return [[data]]

def merge_split_alignments(data):
    """Merge split BAM inputs generated by common workflow language runs.
    """
    assert isinstance(data, (list, tuple)) and len(data) == 1
    data = data[0]
    data = _merge_align_bams(data)
    data = _merge_hla_fastq_inputs(data)
    return [[data]]

def _merge_align_bams(data):
    """Merge multiple alignment BAMs, including split and discordant reads.
    """
    for key in (["work_bam"], ["work_bam-plus", "disc"], ["work_bam-plus", "sr"]):
        in_files = tz.get_in(key, data)
        if in_files:
            if not isinstance(in_files, (list, tuple)):
                in_files = [in_files]
            ext = "-%s" % key[-1] if len(key) > 1 else ""
            out_file = os.path.join(dd.get_work_dir(data), "align", dd.get_sample_name(data),
                                    "%s-sort%s.bam" % (dd.get_sample_name(data), ext))
            merged_file = merge_bam_files(in_files, utils.safe_makedir(os.path.dirname(out_file)),
                                          data["config"], out_file=out_file)
            data = tz.update_in(data, key, lambda x: merged_file)
    if "align_bam" in data and "work_bam" in data:
        data["align_bam"] = data["work_bam"]
    return data

def _merge_hla_fastq_inputs(data):
    """Merge HLA inputs from a split initial alignment.
    """
    hla_key = ["hla", "fastq"]
    hla_dirs = tz.get_in(hla_key, data)
    hla_outdir = None
    if hla_dirs:
        if not isinstance(hla_dirs, (list, tuple)):
            hla_dirs = [hla_dirs]
        hla_files = collections.defaultdict(list)
        for hla_dir in hla_dirs:
            if hla_dir and hla_dir != "None":
                for hla_file in sorted(list(glob.glob(os.path.join(hla_dir, ".*.fq")))):
                    hlatype = re.search(".hla.(?P<hlatype>[\w-]+).fq", hla_file)
                    hla_files[hlatype].append(hla_file)
        if len(hla_files) > 0:
            hla_outdir = os.path.join(dd.get_work_dir(data), "align", dd.get_sample_name(data), "hla")
            for hlatype, files in hla_files:
                out_file = os.path.join(hla_outdir, "%s-%s.fq" % (dd.get_sample_name(data), hlatype))
                optitype.combine_hla_fqs(files, out_file, data)
    data = tz.update_in(data, hla_key, lambda x: hla_outdir)
    return data

def prepare_bcbio_samples(sample):
    """
    Function that will use specific function to merge input files
    """
    if sample['fn'] == "fq_merge":
        out_file = fq_merge(sample['files'], sample['out_file'], sample['config'])
    elif sample['fn'] == "bam_merge":
        out_file = bam_merge(sample['files'], sample['out_file'], sample['config'])
    sample['out_file'] = out_file
    return [sample]
