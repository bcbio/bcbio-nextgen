"""Support for Copy Number Variations (CNVs) with GATK4

https://software.broadinstitute.org/gatk/documentation/article?id=11682
https://gatkforums.broadinstitute.org/dsde/discussion/11683/
"""
import glob
import os
import shutil

import numpy as np
import toolz as tz

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.variation import bedutils, vcfutils

def run(items, background=None):
    """Detect copy number variations from batched set of samples using GATK4 CNV calling.

    TODO: implement germline calling with DetermineGermlineContigPloidy and GermlineCNVCaller
    """
    if not background: background = []
    paired = vcfutils.get_paired(items + background)
    if paired:
        out = _run_paired(paired)
    else:
        out = items
        logger.warn("GATK4 CNV calling currently only available for somatic samples: %s" %
                    ", ".join([dd.get_sample_name(d) for d in items + background]))
    return out

def _run_paired(paired):
    """Run somatic variant calling pipeline.
    """
    from bcbio.structural import titancna
    work_dir = _sv_workdir(paired.tumor_data)
    seg_files = model_segments(tz.get_in(["depth", "bins", "normalized"], paired.tumor_data),
                               work_dir, paired)
    call_file = call_copy_numbers(seg_files["seg"], work_dir, paired.tumor_data)
    out = []
    if paired.normal_data:
        out.append(paired.normal_data)
    if "sv" not in paired.tumor_data:
        paired.tumor_data["sv"] = []
    paired.tumor_data["sv"].append({"variantcaller": "gatk-cnv",
                                    "call_file": call_file,
                                    "vrn_file": titancna.to_vcf(call_file, "GATK4-CNV", _get_seg_header,
                                                                _seg_to_vcf, paired.tumor_data),
                                    "seg": seg_files["seg"],
                                    "plot": plot_model_segments(seg_files, work_dir, paired.tumor_data)})
    out.append(paired.tumor_data)
    return out

def call_copy_numbers(seg_file, work_dir, data):
    """Call copy numbers from a normalized and segmented input file.
    """
    out_file = os.path.join(work_dir, "%s-call.seg" % dd.get_sample_name(data))
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            params = ["-T", "CallCopyRatioSegments",
                      "-I", seg_file, "-O", tx_out_file]
            _run_with_memory_scaling(params, tx_out_file, data)
    return out_file

def plot_model_segments(seg_files, work_dir, data):
    """Diagnostic plots of segmentation and inputs.
    """
    from bcbio.heterogeneity import chromhacks
    out_file = os.path.join(work_dir, "%s.modeled.png" % dd.get_sample_name(data))
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            dict_file = utils.splitext_plus(dd.get_ref_file(data))[0] + ".dict"
            plot_dict = os.path.join(os.path.dirname(tx_out_file), os.path.basename(dict_file))
            with open(dict_file) as in_handle:
                with open(plot_dict, "w") as out_handle:
                    for line in in_handle:
                        if line.startswith("@SQ"):
                            cur_chrom = [x.split(":", 1)[1].strip()
                                         for x in line.split("\t") if x.startswith("SN:")][0]
                            if chromhacks.is_autosomal_or_sex(cur_chrom):
                                out_handle.write(line)
                        else:
                            out_handle.write(line)
            params = ["-T", "PlotModeledSegments",
                      "--denoised-copy-ratios", tz.get_in(["depth", "bins", "normalized"], data),
                      "--segments", seg_files["final_seg"],
                      "--allelic-counts", seg_files["tumor_hets"],
                      "--sequence-dictionary", plot_dict,
                      "--minimum-contig-length", "10",
                      "--output-prefix", dd.get_sample_name(data),
                      "-O", os.path.dirname(tx_out_file)]
            _run_with_memory_scaling(params, tx_out_file, data)
    return {"seg": out_file}

def model_segments(copy_file, work_dir, paired):
    """Perform segmentation on input copy number log2 ratio file.
    """
    out_file = os.path.join(work_dir, "%s.cr.seg" % dd.get_sample_name(paired.tumor_data))
    tumor_counts, normal_counts = heterogzygote_counts(paired)
    if not utils.file_exists(out_file):
        with file_transaction(paired.tumor_data, out_file) as tx_out_file:
            params = ["-T", "ModelSegments",
                      "--denoised-copy-ratios", copy_file,
                      "--allelic-counts", tumor_counts,
                      "--output-prefix", dd.get_sample_name(paired.tumor_data),
                      "-O", os.path.dirname(tx_out_file)]
            if normal_counts:
                params += ["--normal-allelic-counts", normal_counts]
            _run_with_memory_scaling(params, tx_out_file, paired.tumor_data)
            for tx_fname in glob.glob(os.path.join(os.path.dirname(tx_out_file),
                                                   "%s*" % dd.get_sample_name(paired.tumor_data))):
                shutil.copy(tx_fname, os.path.join(work_dir, os.path.basename(tx_fname)))
    return {"seg": out_file, "tumor_hets": out_file.replace(".cr.seg", ".hets.tsv"),
            "final_seg": out_file.replace(".cr.seg", ".modelFinal.seg")}

def denoise(data, pon, work_dir):
    """Normalize read counts using panel of normal background or GC/mappability
    """
    std_file = os.path.join(work_dir, "%s-crstandardized.tsv" % dd.get_sample_name(data))
    denoise_file = os.path.join(work_dir, "%s-crdenoised.tsv" % dd.get_sample_name(data))
    if not utils.file_exists(std_file):
        with file_transaction(data, std_file, denoise_file) as (tx_std_file, tx_denoise_file):
            params = ["-T", "DenoiseReadCounts",
                      "-I", tz.get_in(["depth", "bins", "target"], data),
                      "--standardized-copy-ratios", tx_std_file,
                      "--denoised-copy-ratios", tx_denoise_file]
            if pon:
                params += ["--count-panel-of-normals", pon]
            else:
                params += ["--annotated-intervals", tz.get_in(["regions", "bins", "gcannotated"], data)]
            _run_with_memory_scaling(params, tx_std_file, data)
    return denoise_file if pon else std_file

def create_panel_of_normals(items, group_id, work_dir):
    """Create a panel of normals from one or more background read counts.
    """
    out_file = os.path.join(work_dir, "%s-%s-pon.hdf5" % (dd.get_sample_name(items[0]), group_id))
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            params = ["-T", "CreateReadCountPanelOfNormals",
                      "-O", tx_out_file,
                      "--annotated-intervals", tz.get_in(["regions", "bins", "gcannotated"], items[0])]
            for data in items:
                params += ["-I", tz.get_in(["depth", "bins", "target"], data)]
            _run_with_memory_scaling(params, tx_out_file, items[0], ld_preload=True)
    return out_file

def pon_to_bed(pon_file, out_dir, data):
    """Extract BED intervals from a GATK4 hdf5 panel of normal file.
    """
    out_file = os.path.join(out_dir, "%s-intervals.bed" % (utils.splitext_plus(os.path.basename(pon_file))[0]))
    if not utils.file_uptodate(out_file, pon_file):
        import h5py
        with file_transaction(data, out_file) as tx_out_file:
            with h5py.File(pon_file, "r") as f:
                with open(tx_out_file, "w") as out_handle:
                    intervals = f["original_data"]["intervals"]
                    for i in range(len(intervals["transposed_index_start_end"][0])):
                        chrom = intervals["indexed_contig_names"][intervals["transposed_index_start_end"][0][i]]
                        start = int(intervals["transposed_index_start_end"][1][i]) - 1
                        end = int(intervals["transposed_index_start_end"][2][i])
                        out_handle.write("%s\t%s\t%s\n" % (chrom, start, end))
    return out_file

def prepare_intervals(data, region_file, work_dir):
    """Prepare interval regions for targeted and gene based regions.
    """
    target_file = os.path.join(work_dir, "%s-target.interval_list" % dd.get_sample_name(data))
    if not utils.file_uptodate(target_file, region_file):
        with file_transaction(data, target_file) as tx_out_file:
            params = ["-T", "PreprocessIntervals", "-R", dd.get_ref_file(data),
                      "--interval-merging-rule", "OVERLAPPING_ONLY",
                      "-O", tx_out_file]
            if dd.get_coverage_interval(data) == "genome":
                params += ["--bin-length", "1000", "--padding", "0"]
            else:
                params += ["-L", region_file, "--bin-length", "0", "--padding", "250"]
            _run_with_memory_scaling(params, tx_out_file, data)
    return target_file

def annotate_intervals(target_file, data):
    """Provide GC annotated intervals for error correction during panels and denoising.

    TODO: include mappability and segmentation duplication inputs
    """
    out_file = "%s-gcannotated.tsv" % utils.splitext_plus(target_file)[0]
    if not utils.file_uptodate(out_file, target_file):
        with file_transaction(data, out_file) as tx_out_file:
            params = ["-T", "AnnotateIntervals", "-R", dd.get_ref_file(data),
                      "-L", target_file,
                      "--interval-merging-rule", "OVERLAPPING_ONLY",
                      "-O", tx_out_file]
            _run_with_memory_scaling(params, tx_out_file, data)
    return out_file

def collect_read_counts(data, work_dir):
    """Count reads in defined bins using CollectReadCounts.
    """
    out_file = os.path.join(work_dir, "%s-target-coverage.hdf5" % dd.get_sample_name(data))
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            params = ["-T", "CollectReadCounts", "-I", dd.get_align_bam(data),
                      "-L", tz.get_in(["regions", "bins", "target"], data),
                      "--interval-merging-rule", "OVERLAPPING_ONLY",
                      "-O", tx_out_file, "--format", "HDF5"]
            _run_with_memory_scaling(params, tx_out_file, data)
    return out_file

def heterogzygote_counts(paired):
    """Provide tumor/normal counts at population heterozyogte sites with CollectAllelicCounts.
    """
    work_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(paired.tumor_data), "structural", "counts"))
    key = "germline_het_pon"
    het_bed = tz.get_in(["genome_resources", "variation", key], paired.tumor_data)
    vr = bedutils.population_variant_regions([x for x in [paired.tumor_data, paired.normal_data] if x])
    cur_het_bed = bedutils.intersect_two(het_bed, vr, work_dir, paired.tumor_data)
    tumor_counts = _run_collect_allelic_counts(cur_het_bed, key, work_dir, paired.tumor_data)
    normal_counts = (_run_collect_allelic_counts(cur_het_bed, key, work_dir, paired.normal_data)
                     if paired.normal_data else None)
    if normal_counts:
        tumor_counts, normal_counts = _filter_by_normal(tumor_counts, normal_counts, paired.tumor_data)
    return tumor_counts, normal_counts

def _filter_by_normal(tumor_counts, normal_counts, data):
    """Filter count files based on normal frequency and median depth, avoiding high depth regions.

    For frequency, restricts normal positions to those between 0.4 and 0.65

    For depth, matches approach used in AMBER to try and avoid problematic genomic regions
    with high count in the normal:
    https://github.com/hartwigmedical/hmftools/tree/master/amber#usage
    """
    from bcbio.heterogeneity import bubbletree
    fparams = bubbletree.NORMAL_FILTER_PARAMS
    tumor_out = "%s-normfilter%s" % utils.splitext_plus(tumor_counts)
    normal_out = "%s-normfilter%s" % utils.splitext_plus(normal_counts)
    if not utils.file_uptodate(tumor_out, tumor_counts):
        with file_transaction(data, tumor_out, normal_out) as (tx_tumor_out, tx_normal_out):
            median_depth = _get_normal_median_depth(normal_counts)
            min_normal_depth = median_depth * fparams["min_depth_percent"]
            max_normal_depth = median_depth * fparams["max_depth_percent"]
            with open(tumor_counts) as tumor_handle:
                with open(normal_counts) as normal_handle:
                    with open(tx_tumor_out, "w") as tumor_out_handle:
                        with open(tx_normal_out, "w") as normal_out_handle:
                            header = None
                            for t, n in zip(tumor_handle, normal_handle):
                                if header is None:
                                    if not n.startswith("@"):
                                        header = n.strip().split()
                                    tumor_out_handle.write(t)
                                    normal_out_handle.write(n)
                                elif (_normal_passes_depth(header, n, min_normal_depth, max_normal_depth) and
                                      _normal_passes_freq(header, n, fparams)):
                                    tumor_out_handle.write(t)
                                    normal_out_handle.write(n)
    return tumor_out, normal_out

def _normal_passes_freq(header, line, fparams):
    vals = dict(zip(header, line.strip().split()))
    cur_depth = float(vals["REF_COUNT"]) + int(vals["ALT_COUNT"])
    if cur_depth > 0:
        cur_freq = float(vals["ALT_COUNT"]) / cur_depth
    else:
        cur_freq = 0.0
    return cur_freq >= fparams["min_freq_narrow"] and cur_freq <= fparams["max_freq_narrow"]

def _normal_passes_depth(header, line, min_normal_depth, max_normal_depth):
    vals = dict(zip(header, line.strip().split()))
    cur_depth = int(vals["REF_COUNT"]) + int(vals["ALT_COUNT"])
    return cur_depth >= min_normal_depth and cur_depth <= max_normal_depth

def _get_normal_median_depth(normal_counts):
    depths = []
    with open(normal_counts) as in_handle:
        header = None
        for line in in_handle:
            if header is None and not line.startswith("@"):
                header = line.strip().split()
            elif header:
                n_vals = dict(zip(header, line.strip().split()))
                depths.append(int(n_vals["REF_COUNT"]) + int(n_vals["ALT_COUNT"]))
    return np.median(depths)

def _run_collect_allelic_counts(pos_file, pos_name, work_dir, data):
    """Counts by alleles for a specific sample and set of positions.
    """
    out_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "structural", "counts"))
    out_file = os.path.join(out_dir, "%s-%s-counts.tsv" % (dd.get_sample_name(data), pos_name))
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            params = ["-T", "CollectAllelicCounts", "-L", pos_file, "-I", dd.get_align_bam(data),
                      "-R", dd.get_ref_file(data), "-O", tx_out_file]
            _run_with_memory_scaling(params, tx_out_file, data)
    return out_file

def _run_with_memory_scaling(params, tx_out_file, data, ld_preload=False):
    num_cores = dd.get_num_cores(data)
    memscale = {"magnitude": 0.9 * num_cores, "direction": "increase"} if num_cores > 1 else None
    # Ignore tools_off: [gatk4], since it doesn't apply to GATK CNV calling
    config = utils.deepish_copy(data["config"])
    if "gatk4" in dd.get_tools_off({"config": config}):
        config["algorithm"]["tools_off"].remove("gatk4")
    broad_runner = broad.runner_from_config(config)
    broad_runner.run_gatk(params, os.path.dirname(tx_out_file), memscale=memscale, ld_preload=ld_preload)

# ## VCF output

def _get_seg_header(in_handle):
    for line in in_handle:
        if not line.startswith("@"):
            break
    return line.strip().split("\t"), in_handle

def _seg_to_vcf(vals):
    """Convert GATK CNV calls seg output to a VCF line.
    """
    call_to_cn = {"+": 3, "-": 1}
    call_to_type = {"+": "DUP", "-": "DEL"}
    if vals["CALL"] not in ["0"]:
        info = ["FOLD_CHANGE_LOG=%s" % vals["MEAN_LOG2_COPY_RATIO"],
                "PROBES=%s" % vals["NUM_POINTS_COPY_RATIO"],
                "SVTYPE=%s" % call_to_type[vals["CALL"]],
                "SVLEN=%s" % (int(vals["END"]) - int(vals["START"])),
                "END=%s" % vals["END"],
                "CN=%s" % call_to_cn[vals["CALL"]]]
        return [vals["CONTIG"], vals["START"], ".", "N", "<%s>" % call_to_type[vals["CALL"]], ".",
                ".", ";".join(info), "GT", "0/1"]

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(dd.get_work_dir(data), "structural",
                                           dd.get_sample_name(data), "gatk-cnv"))
