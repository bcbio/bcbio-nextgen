"""Examine callable regions following genome mapping of short reads.

Identifies callable analysis regions surrounded by larger regions lacking
aligned bases. This allows parallelization of smaller chromosome chunks
through post-processing and variant calling, with each sub-section
mapping handled separately.

Regions are split to try to maintain relative uniformity across the
genome and avoid extremes of large blocks or large numbers of
small blocks.
"""
import contextlib
import copy
import operator
import os
import shutil

import pybedtools
import pysam
from py_descriptive_statistics import Enum as Stats

from bcbio import utils, broad
from bcbio.log import logger, setup_local_logging
from bcbio.distributed.messaging import parallel_runner
from bcbio.distributed.split import parallel_split_combine
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import shared

def parallel_callable_loci(in_bam, ref_file, config):
    num_cores = config["algorithm"].get("num_cores", 1)
    config = copy.deepcopy(config)
    config["algorithm"]["memory_adjust"] = {"direction": "decrease", "magnitude": 2}
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
    if data["config"].get("parallel", {}).get("log_queue"):
        handler = setup_local_logging(data["config"], data["config"]["parallel"])
    else:
        handler = None
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
    if handler and hasattr(handler, "close"):
        handler.close()
    return [{"callable_bed": out_file, "config": data["config"], "work_bam": data["work_bam"]}]

def sample_callable_bed(bam_file, ref_file, config):
    """Retrieve callable regions for a sample subset by defined analysis regions.
    """
    out_file = "%s-callable_sample.bed" % os.path.splitext(bam_file)[0]
    callable_bed = parallel_callable_loci(bam_file, ref_file, config)
    input_regions_bed = config["algorithm"].get("variant_regions", None)
    if not utils.file_uptodate(out_file, callable_bed):
        with file_transaction(out_file) as tx_out_file:
            callable_regions = pybedtools.BedTool(callable_bed)
            filter_regions = callable_regions.filter(lambda x: x.name == "CALLABLE")
            if input_regions_bed:
                if not utils.file_uptodate(out_file, input_regions_bed):
                    input_regions = pybedtools.BedTool(input_regions_bed)
                    filter_regions.intersect(input_regions).saveas(tx_out_file)
            else:
                filter_regions.saveas(tx_out_file)
    return out_file

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
        if input_nblock == ref_regions:
            raise ValueError("Input variant_region file (%s) "
                             "excludes all genomic regions. Do the chromosome names "
                             "in the BED file match your genome (chr1 vs 1)?" % input_regions_bed)
        all_intervals = _combine_regions([input_nblock, nblock_regions], ref_regions)
        return all_intervals.merge()
    else:
        return nblock_regions

class NBlockRegionPicker:
    """Choose nblock regions reasonably spaced across chromosomes.

    This avoids excessively large blocks and also large numbers of tiny blocks
    by splitting to a defined number of blocks.

    Assumes to be iterating over an ordered input file and needs re-initiation
    with each new file processed as it keeps track of previous blocks to
    maintain the splitting.
    """
    def __init__(self, ref_regions, config):
        self._chr_last_blocks = {}
        target_blocks = int(config["algorithm"].get("nomap_split_targets", 2000))
        self._target_size = self._get_target_size(target_blocks, ref_regions)

    def _get_target_size(self, target_blocks, ref_regions):
        size = 0
        for x in ref_regions:
            size += (x.end - x.start)
        return size // target_blocks

    def include_block(self, x):
        """Check for inclusion of block based on distance from previous.
        """
        last_pos = self._chr_last_blocks.get(x.chrom, 0)
        if (x.start - last_pos) > self._target_size:
            self._chr_last_blocks[x.chrom] = x.stop
            return True
        else:
            return False

def block_regions(in_bam, ref_file, config):
    """Find blocks of regions for analysis from mapped input BAM file.

    Identifies islands of callable regions, surrounding by regions
    with no read support, that can be analyzed independently.
    """
    min_n_size = int(config["algorithm"].get("nomap_split_size", 100))
    callable_bed = parallel_callable_loci(in_bam, ref_file, config)
    nblock_bed = "%s-nblocks%s" % os.path.splitext(callable_bed)
    callblock_bed = "%s-callableblocks%s" % os.path.splitext(callable_bed)
    if not utils.file_uptodate(nblock_bed, callable_bed):
        ref_regions = get_ref_bedtool(ref_file, config)
        nblock_regions = _get_nblock_regions(callable_bed, min_n_size)
        nblock_regions = _add_config_regions(nblock_regions, ref_regions, config)
        nblock_regions.saveas(nblock_bed)
        ref_regions.subtract(nblock_bed).merge(d=min_n_size).saveas(callblock_bed)
    return callblock_bed, nblock_bed

def _write_bed_regions(sample, final_regions, out_file, out_file_ref):
    ref_regions = get_ref_bedtool(sample["sam_ref"], sample["config"])
    noanalysis_regions = ref_regions.subtract(final_regions)
    final_regions.saveas(out_file)
    noanalysis_regions.saveas(out_file_ref)

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

def _needs_region_update(out_file, samples):
    """Check if we need to update BED file of regions, supporting back compatibility.
    """
    nblock_files = [x["regions"]["nblock"] for x in samples if "regions" in x]
    # For older approaches and do not create a new set of analysis
    # regions, since the new algorithm will re-do all BAM and variant
    # steps with new regions
    for nblock_file in nblock_files:
        test_old = nblock_file.replace("-nblocks", "-analysisblocks")
        if os.path.exists(test_old):
            return False
    # Check if any of the local files have changed so we need to refresh
    for noblock_file in nblock_files:
        if not utils.file_uptodate(out_file, noblock_file):
            return True
    return False

def combine_sample_regions(samples):
    """Create global set of callable regions for multi-sample calling.

    Intersects all non-callable (nblock) regions from all samples,
    producing a global set of callable regions.
    """
    config = samples[0]["config"]
    work_dir = samples[0]["dirs"]["work"]
    analysis_file = os.path.join(work_dir, "analysis_blocks.bed")
    no_analysis_file = os.path.join(work_dir, "noanalysis_blocks.bed")
    min_n_size = int(config["algorithm"].get("nomap_split_size", 100))

    if not utils.file_exists(analysis_file) or _needs_region_update(analysis_file, samples):
        # Combine all nblocks into a final set of intersecting regions
        # without callable bases. HT @brentp for intersection approach
        # https://groups.google.com/forum/?fromgroups#!topic/bedtools-discuss/qA9wK4zN8do
        nblock_regions = reduce(operator.add,
                                (pybedtools.BedTool(x["regions"]["nblock"])
                                 for x in samples if "regions" in x))
        ref_regions = get_ref_bedtool(samples[0]["sam_ref"], config)
        block_filter = NBlockRegionPicker(ref_regions, config)
        nblock_size_filtered = nblock_regions.filter(block_filter.include_block).saveas()
        if len(nblock_size_filtered) > len(ref_regions):
            final_nblock_regions = nblock_size_filtered
        else:
            final_nblock_regions = nblock_regions
        final_regions = ref_regions.subtract(final_nblock_regions).merge(d=min_n_size)
        _write_bed_regions(samples[0], final_regions, analysis_file, no_analysis_file)
    else:
        final_regions = pybedtools.BedTool(analysis_file)
    _analysis_block_stats(final_regions)
    regions = {"analysis": [(r.chrom, int(r.start), int(r.stop)) for r in final_regions],
               "noanalysis": no_analysis_file,
               "analysis_bed": analysis_file}
    return regions
