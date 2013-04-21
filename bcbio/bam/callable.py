"""Examine callable regions following genome mapping of short reads.

Identifies callable analysis regions surrounded by larger regions lacking
aligned bases. This allows parallelization of smaller chromosome chunks
through post-processing and variant calling, with each sub-section
mapping handled separately.
"""
import contextlib
import os
import shutil

import pybedtools
import pysam
from py_descriptive_statistics import Enum as Stats

from bcbio import utils, broad
from bcbio.log import logger
from bcbio.distributed.messaging import parallel_runner
from bcbio.distributed.split import parallel_split_combine
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import shared

def parallel_callable_loci(in_bam, ref_file, config):
    num_cores = config["algorithm"].get("num_cores", 1)
    data = {"work_bam": in_bam, "sam_ref": ref_file, "config": config}
    parallel = {"type": "local", "cores": num_cores, "module": "bcbio.distributed"}
    runner = parallel_runner(parallel, {}, config)
    split_fn = shared.process_bam_by_chromosome("-callable.bed", "work_bam")
    out = parallel_split_combine([[data]], split_fn, runner,
                                 "calc_callable_loci", "combine_bed",
                                 "callable_bed", ["config"])[0]
    return out[0]["callable_bed"]

def combine_bed(in_files, out_file, config):
    """Combine multiple BED files into a single output.
    """
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                for in_file in in_files:
                    with open(in_file) as in_handle:
                        shutil.copyfileobj(in_handle, out_handle)
    return out_file

def calc_callable_loci(data, region=None, out_file=None):
    """Determine callable bases for input BAM using Broad's CallableLoci walker.

    http://www.broadinstitute.org/gatk/gatkdocs/
    org_broadinstitute_sting_gatk_walkers_coverage_CallableLoci.html
    """
    broad_runner = broad.runner_from_config(data["config"])
    if out_file is None:
        out_file = "%s-callable.bed" % os.path.splitext(data["work_bam"])[0]
    out_summary = "%s-callable-summary.txt" % os.path.splitext(data["work_bam"])[0]
    variant_regions = data["config"]["algorithm"].get("variant_regions", None)
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            broad_runner.run_fn("picard_index", data["work_bam"])
            params = ["-T", "CallableLoci",
                      "-R", data["sam_ref"],
                      "-I", data["work_bam"],
                      "--out", tx_out_file,
                      "--summary", out_summary]
            ready_region = shared.subset_variant_regions(variant_regions, region, tx_out_file)
            if ready_region:
                params += ["-L", ready_region]
            if ((variant_regions and ready_region and os.path.isfile(ready_region))
                 or not variant_regions or not region):
                broad_runner.run_gatk(params)
            else:
                with open(out_file, "w") as out_handle:
                    for tregion in get_ref_bedtool(data["sam_ref"], data["config"]):
                        if tregion.chrom == region:
                            out_handle.write("%s\t%s\t%s\tNO_COVERAGE\n" %
                                             (tregion.chrom, tregion.start, tregion.stop))
    return [{"callable_bed": out_file, "config": data["config"], "work_bam": data["work_bam"]}]

def get_ref_bedtool(ref_file, config):
    """Retrieve a pybedtool BedTool object with reference sizes from input reference.
    """
    broad_runner = broad.runner_from_config(config)
    ref_dict = broad_runner.run_fn("picard_index_ref", ref_file)
    ref_lines = []
    with contextlib.closing(pysam.Samfile(ref_dict, "r")) as ref_sam:
        for sq in ref_sam.header["SQ"]:
            ref_lines.append("%s\t%s\t%s" % (sq["SN"], 0, sq["LN"]))
    return pybedtools.BedTool("\n".join(ref_lines), from_string=True)

def _get_nblock_regions(in_file, min_n_size):
    """Retrieve coordinates of regions in reference genome with no mapping.
    These are potential breakpoints for parallelizing analysis.
    """
    out_lines = []
    with open(in_file) as in_handle:
        for line in in_handle:
            contig, start, end, ctype = line.rstrip().split()
            if (ctype in ["REF_N", "NO_COVERAGE"] and
                  int(end) - int(start) > min_n_size):
                out_lines.append("%s\t%s\t%s\n" % (contig, start, end))
    return pybedtools.BedTool("\n".join(out_lines), from_string=True)

def _combine_regions(all_regions, ref_regions):
    """Combine multiple BEDtools regions of regions into sorted final BEDtool.
    """
    chrom_order = {}
    for i, x in enumerate(ref_regions):
        chrom_order[x.chrom] = i
    def wchrom_key(x):
        chrom, start, end = x
        return (chrom_order[chrom], start, end)
    all_intervals = []
    for region_group in all_regions:
        for region in region_group:
            all_intervals.append((region.chrom, int(region.start), int(region.stop)))
    all_intervals.sort(key=wchrom_key)
    bed_lines = ["%s\t%s\t%s" % (c, s, e) for (c, s, e) in all_intervals]
    return pybedtools.BedTool("\n".join(bed_lines), from_string=True)

def _add_config_regions(nblock_regions, ref_regions, config):
    """Add additional nblock regions based on configured regions to call.
    Identifies user defined regions which we should not be analyzing.
    """
    input_regions_bed = config["algorithm"].get("variant_regions", None)
    if input_regions_bed:
        input_regions = pybedtools.BedTool(input_regions_bed)
        # work around problem with single region not subtracted correctly.
        if len(input_regions) == 1:
            str_regions = str(input_regions[0]).strip()
            input_regions = pybedtools.BedTool("%s\n%s" % (str_regions, str_regions),
                                               from_string=True)
        input_nblock = ref_regions.subtract(input_regions)
        all_intervals = _combine_regions([input_nblock, nblock_regions], ref_regions)
        return all_intervals.merge()
    else:
        return nblock_regions

def _avoid_small_regions(regions, min_size, ref_regions):
    """Expand regions less than min_size to merge with nearby regions.
    This avoids large numbers of very small regions, which are
    problematic for parallelizing.
    """
    chromsizes = {}
    for r in ref_regions:
        chromsizes[r.chrom] = (r.start, r.stop)
    small_regions = regions.filter(lambda b: (b.stop - b.start) < min_size)
    expand_small_regions = small_regions.slop(g=chromsizes, b=min_size * 4)
    if len(expand_small_regions) > 0:
        return regions.cat(expand_small_regions, postmerge=True)
    else:
        return regions

def block_regions(in_bam, ref_file, config):
    """Find blocks of regions for analysis from mapped input BAM file.

    Identifies islands of callable regions, surrounding by regions
    with no read support, that can be analyzed independently.
    """
    min_n_size = int(config["algorithm"].get("nomap_split_size", 5000))
    callable_bed = parallel_callable_loci(in_bam, ref_file, config)
    block_bed = "%s-analysisblocks%s" % os.path.splitext(callable_bed)
    if utils.file_uptodate(block_bed, callable_bed):
        ready_regions = pybedtools.BedTool(block_bed)
    else:
        ref_regions = get_ref_bedtool(ref_file, config)
        nblock_regions = _get_nblock_regions(callable_bed, min_n_size)
        nblock_regions = _add_config_regions(nblock_regions, ref_regions, config)
        ready_regions = ref_regions.subtract(nblock_regions)
        ready_regions = ready_regions.merge(d=min_n_size)
        ready_regions = _avoid_small_regions(ready_regions, min_n_size, ref_regions)
        ready_regions.saveas(block_bed)
    return [(r.chrom, int(r.start), int(r.stop)) for r in ready_regions]

def _write_bed_regions(sample, final_regions):
    work_dir = sample["dirs"]["work"]
    ref_regions = get_ref_bedtool(sample["sam_ref"], sample["config"])
    noanalysis_regions = ref_regions.subtract(final_regions)
    out_file = os.path.join(work_dir, "analysis_blocks.bed")
    out_file_ref = os.path.join(work_dir, "noanalysis_blocks.bed")
    final_regions.saveas(out_file)
    noanalysis_regions.saveas(out_file_ref)
    return out_file, out_file_ref

def _analysis_block_stats(regions):
    """Provide statistics on sizes and number of analysis blocks.
    """
    prev = None
    between_sizes = []
    region_sizes = []
    for region in regions:
        if prev and prev.chrom == region.chrom:
            between_sizes.append(region.start - prev.end)
        region_sizes.append(region.end - region.start)
        prev = region
    def descriptive_stats(xs):
        if len(xs) < 2:
            return xs
        calc = Stats(xs)
        parts = ["min: %s" % min(xs),
                 "5%%: %s" % calc.percentile(5),
                 "25%%: %s" % calc.percentile(25),
                 "median: %s" % calc.percentile(50),
                 "75%%: %s" % calc.percentile(75),
                 "95%%: %s" % calc.percentile(95),
                 "99%%: %s" % calc.percentile(99),
                 "max: %s" % max(xs)]
        return "\n".join(["  " + x for x in parts])
    logger.info("Identified %s parallel analysis blocks\n" % len(region_sizes) +
                "Block sizes:\n%s\n" % descriptive_stats(region_sizes) +
                "Between block sizes:\n%s\n" % descriptive_stats(between_sizes))
    if len(region_sizes) == 0:
        raise ValueError("No callable analysis regions found in all samples")

def combine_sample_regions(samples):
    """Combine islands of callable regions from multiple samples.
    Creates a global set of callable samples usable across a
    project with multi-sample calling.
    """
    min_n_size = int(samples[0]["config"]["algorithm"].get("nomap_split_size", 5000))
    final_regions = None
    all_regions = []
    for regions in (x["regions"] for x in samples if "regions" in x):
        bed_lines = ["%s\t%s\t%s" % (c, s, e) for (c, s, e) in regions]
        all_regions.append(pybedtools.BedTool("\n".join(bed_lines), from_string=True))
    if len(all_regions) == 0:
        final_regions = []
    elif len(all_regions) == 1:
        final_regions = all_regions[0]
    else:
        ref_bedtool = get_ref_bedtool(samples[0]["sam_ref"], samples[0]["config"])
        combo_regions = _combine_regions(all_regions, ref_bedtool)
        final_regions = combo_regions.merge(d=min_n_size)
    if len(all_regions) > 0:
        _analysis_block_stats(final_regions)
        analysis_file, no_analysis_file = _write_bed_regions(samples[0], final_regions)
    else:
        analysis_file, no_analysis_file = None, None
    regions = {"analysis": [(r.chrom, int(r.start), int(r.stop)) for r in final_regions],
               "noanalysis": no_analysis_file,
               "analysis_bed": analysis_file}
    return regions
