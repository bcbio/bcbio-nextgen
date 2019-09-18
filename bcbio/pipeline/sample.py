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
from bcbio.cwl import cwlutils
from bcbio.log import logger
from bcbio.distributed import objectstore
from bcbio.pipeline.merge import merge_bam_files
from bcbio.bam import callable, readstats, trim
from bcbio.hla import optitype
from bcbio.ngsalign import postalign
from bcbio.pipeline.fastq import get_fastq_files
from bcbio.pipeline.alignment import align_to_sort_bam
from bcbio.pipeline import cleanbam
from bcbio.variation import coverage, recalibrate
from bcbio.variation import multi as vmulti
import bcbio.pipeline.datadict as dd
from bcbio.pipeline.fastq import merge as fq_merge
from bcbio.bam import merge as bam_merge
from bcbio.pipeline.sra import query_gsm, query_srr
from bcbio.bam import skewer
from bcbio.structural.seq2c import prep_seq2c_bed
from bcbio.variation.bedutils import clean_file, merge_overlaps
from bcbio.structural import get_svcallers, regions
from bcbio.qc import samtools

def prepare_sample(data):
    """Prepare a sample to be run, potentially converting from BAM to
    FASTQ and/or downsampling the number of reads for a test run
    """
    data = utils.to_single_data(data)
    logger.debug("Preparing %s" % data["rgnames"]["sample"])
    data["files"] = get_fastq_files(data)
    # get_fastq_files swaps over quality scores to standard, unless trimming
    if not(dd.get_trim_reads(data)):
        data = dd.set_quality_format(data, "standard")
    return [[data]]

def trim_sample(data):
    """Trim from a sample with the provided trimming method.
    Support methods: read_through.
    """
    data = utils.to_single_data(data)
    trim_reads = dd.get_trim_reads(data)
    # this block is to maintain legacy configuration files
    if not trim_reads:
        logger.info("Skipping trimming of %s." % dd.get_sample_name(data))
    else:
        if "skewer" in dd.get_tools_on(data) or trim_reads == "skewer":
            trim_adapters = skewer.trim_adapters
        else:
            trim_adapters = trim.trim_adapters
        out_files = trim_adapters(data)
        data["files"] = out_files
    return [[data]]

# ## Alignment

def _link_bam_file(in_file, new_dir, data):
    """Provide symlinks of BAM file and existing indexes if needed.
    """
    new_dir = utils.safe_makedir(new_dir)
    out_file = os.path.join(new_dir, os.path.basename(in_file))
    if not utils.file_exists(out_file):
        out_file = os.path.join(new_dir, "%s-prealign.bam" % dd.get_sample_name(data))
    if data.get("cwl_keys"):
        # Has indexes, we're okay to go with the original file
        if utils.file_exists(in_file + ".bai"):
            out_file = in_file
        else:
            utils.copy_plus(in_file, out_file)
    else:
        utils.symlink_plus(in_file, out_file)
    return out_file

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
                sup_key = file_key + "_plus"
                if sup_key not in data:
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
        hla_files = None
    else:
        hla_files = sorted(list(glob.glob(os.path.join(hla_dir, "%s.*.fq" % os.path.basename(align_file)))))
    data["hla"]["fastq"] = hla_files
    return data

def process_alignment(data, alt_input=None):
    """Do an alignment of fastq files, preparing a sorted BAM output file.
    """
    data = cwlutils.normalize_missing(utils.to_single_data(data))
    data = cwlutils.unpack_tarballs(data, data)
    fastq1, fastq2 = dd.get_input_sequence_files(data)
    if alt_input:
        fastq1, fastq2 = alt_input
    config = data["config"]
    aligner = config["algorithm"].get("aligner", None)
    if fastq1 and objectstore.file_exists_or_remote(fastq1) and aligner:
        logger.info("Aligning lane %s with %s aligner" % (data["rgnames"]["lane"], aligner))
        data = align_to_sort_bam(fastq1, fastq2, aligner, data)
        if dd.get_correct_umis(data):
            data["work_bam"] = postalign.correct_umis(data)
        if dd.get_umi_consensus(data):
            data["umi_bam"] = dd.get_work_bam(data)
            if fastq2:
                f1, f2, avg_cov = postalign.umi_consensus(data)
                data["config"]["algorithm"]["rawumi_avg_cov"] = avg_cov
                del data["config"]["algorithm"]["umi_type"]
                data["config"]["algorithm"]["mark_duplicates"] = False
                data = align_to_sort_bam(f1, f2, aligner, data)
            else:
                raise ValueError("Single fastq input for UMI processing; fgbio needs paired reads: %s" %
                                 dd.get_sample_name(data))
        data = _add_supplemental_bams(data)
    elif fastq1 and objectstore.file_exists_or_remote(fastq1) and fastq1.endswith(".bam"):
        sort_method = config["algorithm"].get("bam_sort")
        bamclean = config["algorithm"].get("bam_clean")
        if bamclean is True or bamclean == "picard":
            if sort_method and sort_method != "coordinate":
                raise ValueError("Cannot specify `bam_clean: picard` with `bam_sort` other than coordinate: %s"
                                 % sort_method)
            out_bam = cleanbam.picard_prep(fastq1, data["rgnames"], dd.get_ref_file(data), data["dirs"],
                                           data)
        elif bamclean == "fixrg":
            out_bam = cleanbam.fixrg(fastq1, data["rgnames"], dd.get_ref_file(data), data["dirs"], data)
        elif bamclean == "remove_extracontigs":
            out_bam = cleanbam.remove_extracontigs(fastq1, data)
        elif sort_method:
            runner = broad.runner_from_path("picard", config)
            out_file = os.path.join(data["dirs"]["work"], "{}-sort.bam".format(
                os.path.splitext(os.path.basename(fastq1))[0]))
            if not utils.file_exists(out_file):
                work_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "bamclean",
                                                           dd.get_sample_name(data)))
                out_file = os.path.join(work_dir, "{}-sort.bam".format(dd.get_sample_name(data)))
            out_bam = runner.run_fn("picard_sort", fastq1, sort_method, out_file)
        else:
            out_bam = _link_bam_file(fastq1, os.path.join(dd.get_work_dir(data), "prealign",
                                                          dd.get_sample_name(data)), data)
        bam.index(out_bam, data["config"])
        bam.check_header(out_bam, data["rgnames"], dd.get_ref_file(data), data["config"])
        dedup_bam = postalign.dedup_bam(out_bam, data)
        bam.index(dedup_bam, data["config"])
        data["work_bam"] = dedup_bam
    elif fastq1 and objectstore.file_exists_or_remote(fastq1) and fastq1.endswith(".cram"):
        data["work_bam"] = fastq1
    elif fastq1 is None and not dd.get_aligner(data):
        data["config"]["algorithm"]["variantcaller"] = False
        data["work_bam"] = None
    elif not fastq1:
        raise ValueError("No 'files' specified for input sample: %s" % dd.get_sample_name(data))
    elif "kraken" in config["algorithm"]:  # kraken doesn's need bam
        pass
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
        data = cwlutils.normalize_missing(data)
        data = cwlutils.unpack_tarballs(data, data)
        data = clean_inputs(data)
        out.append([data])
    return out

def clean_inputs(data):
    """Clean BED input files to avoid overlapping segments that cause downstream issues.

    Per-merges inputs to avoid needing to call multiple times during later parallel steps.
    """
    if not utils.get_in(data, ("config", "algorithm", "variant_regions_orig")):
        data["config"]["algorithm"]["variant_regions_orig"] = dd.get_variant_regions(data)
    clean_vr = clean_file(dd.get_variant_regions(data), data, prefix="cleaned-")
    merged_vr = merge_overlaps(clean_vr, data)
    data["config"]["algorithm"]["variant_regions"] = clean_vr
    data["config"]["algorithm"]["variant_regions_merged"] = merged_vr

    if dd.get_coverage(data):
        if not utils.get_in(data, ("config", "algorithm", "coverage_orig")):
            data["config"]["algorithm"]["coverage_orig"] = dd.get_coverage(data)
        clean_cov_bed = clean_file(dd.get_coverage(data), data, prefix="cov-", simple=True)
        merged_cov_bed = merge_overlaps(clean_cov_bed, data)
        data["config"]["algorithm"]["coverage"] = clean_cov_bed
        data["config"]["algorithm"]["coverage_merged"] = merged_cov_bed

    if 'seq2c' in get_svcallers(data):
        seq2c_ready_bed = prep_seq2c_bed(data)
        if not seq2c_ready_bed:
            logger.warning("Can't run Seq2C without a svregions or variant_regions BED file")
        else:
            data["config"]["algorithm"]["seq2c_bed_ready"] = seq2c_ready_bed
    elif regions.get_sv_bed(data):
        dd.set_sv_regions(data, clean_file(regions.get_sv_bed(data), data, prefix="svregions-"))
    return data

def postprocess_alignment(data):
    """Perform post-processing steps required on full BAM files.
    Prepares list of callable genome regions allowing subsequent parallelization.
    """
    data = cwlutils.normalize_missing(utils.to_single_data(data))
    data = cwlutils.unpack_tarballs(data, data)
    bam_file = data.get("align_bam") or data.get("work_bam")
    ref_file = dd.get_ref_file(data)
    if vmulti.bam_needs_processing(data) and bam_file and bam_file.endswith(".bam"):
        out_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "align",
                                                  dd.get_sample_name(data)))
        bam_file_ready = os.path.join(out_dir, os.path.basename(bam_file))
        if not utils.file_exists(bam_file_ready):
            utils.symlink_plus(bam_file, bam_file_ready)
        bam.index(bam_file_ready, data["config"])
        covinfo = callable.sample_callable_bed(bam_file_ready, ref_file, data)
        callable_region_bed, nblock_bed = \
            callable.block_regions(covinfo.raw_callable, bam_file_ready, ref_file, data)
        data["regions"] = {"nblock": nblock_bed,
                           "callable": covinfo.raw_callable,
                           "sample_callable": covinfo.callable,
                           "mapped_stats": readstats.get_cache_file(data)}
        data["depth"] = covinfo.depth_files
        data = coverage.assign_interval(data)
        data = samtools.run_and_save(data)
        data = recalibrate.prep_recal(data)
        data = recalibrate.apply_recal(data)
    elif dd.get_variant_regions(data):
        callable_region_bed, nblock_bed = \
            callable.block_regions(dd.get_variant_regions(data), bam_file, ref_file, data)
        data["regions"] = {"nblock": nblock_bed, "callable": dd.get_variant_regions(data),
                           "sample_callable": dd.get_variant_regions(data)}
    return [[data]]

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
        file_key = list(data["combine"].keys())[0]
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
        sup_exts = data.get(file_key + "_plus", {}).keys()
        for ext in list(sup_exts) + [""]:
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
                if len(cur_in_files) > 0:
                    merged_file = merge_bam_files(cur_in_files, os.path.dirname(cur_out_file), data,
                                                  out_file=cur_out_file)
                else:
                    assert os.path.exists(cur_out_file)
                    merged_file = cur_out_file
            if merged_file:
                if ext:
                    data[file_key + "_plus"][ext] = merged_file
                else:
                    data[file_key] = merged_file
        data.pop("region", None)
        data.pop("combine", None)
    return [[data]]

def merge_split_alignments(data):
    """Merge split BAM inputs generated by common workflow language runs.
    """
    data = utils.to_single_data(data)
    data = _merge_align_bams(data)
    data = _merge_hla_fastq_inputs(data)
    return [[data]]

def _merge_align_bams(data):
    """Merge multiple alignment BAMs, including split and discordant reads.
    """
    for key in (["work_bam"], ["work_bam_plus", "disc"], ["work_bam_plus", "sr"], ["umi_bam"]):
        in_files = tz.get_in(key, data, [])
        if not isinstance(in_files, (list, tuple)):
            in_files = [in_files]
        in_files = [x for x in in_files if x and x != "None"]
        if in_files:
            ext = "-%s" % key[-1] if len(key) > 1 else ""
            out_file = os.path.join(dd.get_work_dir(data), "align", dd.get_sample_name(data),
                                    "%s-sort%s.bam" % (dd.get_sample_name(data), ext))
            merged_file = merge_bam_files(in_files, utils.safe_makedir(os.path.dirname(out_file)),
                                          data, out_file=out_file)
            data = tz.update_in(data, key, lambda x: merged_file)
        else:
            data = tz.update_in(data, key, lambda x: None)
    if "align_bam" in data and "work_bam" in data:
        data["align_bam"] = data["work_bam"]
    return data

def _merge_hla_fastq_inputs(data):
    """Merge HLA inputs from a split initial alignment.
    """
    hla_key = ["hla", "fastq"]
    hla_sample_files = [x for x in (tz.get_in(hla_key, data) or []) if x and x != "None"]
    merged_hlas = None
    if hla_sample_files:
        out_files = collections.defaultdict(list)
        for hla_file in utils.flatten(hla_sample_files):
            rehla = re.search(r".hla.(?P<hlatype>[\w-]+).fq", hla_file)
            if rehla:
                hlatype = rehla.group("hlatype")
                out_files[hlatype].append(hla_file)
        if len(out_files) > 0:
            hla_outdir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "align",
                                                         dd.get_sample_name(data), "hla"))
            merged_hlas = []
            for hlatype, files in out_files.items():
                out_file = os.path.join(hla_outdir, "%s-%s.fq" % (dd.get_sample_name(data), hlatype))
                optitype.combine_hla_fqs([(hlatype, f) for f in files], out_file, data)
                merged_hlas.append(out_file)
    data = tz.update_in(data, hla_key, lambda x: merged_hlas)
    return data

def prepare_bcbio_samples(sample):
    """
    Function that will use specific function to merge input files
    """
    logger.info("Preparing %s files %s to merge into %s." % (sample['name'], sample['files'], sample['out_file']))
    if sample['fn'] == "fq_merge":
        out_file = fq_merge(sample['files'], sample['out_file'], sample['config'])
    elif sample['fn'] == "bam_merge":
        out_file = bam_merge(sample['files'], sample['out_file'], sample['config'])
    elif sample['fn'] == "query_gsm":
        out_file = query_gsm(sample['files'], sample['out_file'], sample['config'])
    elif sample['fn'] == "query_srr":
        out_file = query_srr(sample['files'], sample['out_file'], sample['config'])
    sample['out_file'] = out_file
    return [sample]
