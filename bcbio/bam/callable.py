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
import sys

import numpy
import pysam
try:
    import pybedtools
except ImportError:
    pybedtools = None

from bcbio import bam, broad, utils
from bcbio.log import logger
from bcbio.distributed import multi, prun
from bcbio.distributed.split import parallel_split_combine
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import shared

def parallel_callable_loci(in_bam, ref_file, config):
    num_cores = config["algorithm"].get("num_cores", 1)
    config = copy.deepcopy(config)
    config["algorithm"]["memory_adjust"] = {"direction": "decrease", "magnitude": 2}
    data = {"work_bam": in_bam, "sam_ref": ref_file, "config": config}
    parallel = {"type": "local", "cores": num_cores, "module": "bcbio.distributed"}
    items = [[data]]
    with prun.start(parallel, items, config) as runner:
        split_fn = shared.process_bam_by_chromosome("-callable.bed", "work_bam")
        out = parallel_split_combine(items, split_fn, runner,
                                     "calc_callable_loci", "combine_bed",
                                     "callable_bed", ["config"])[0]
    return out[0]["callable_bed"]

@multi.zeromq_aware_logging
def calc_callable_loci(data, region=None, out_file=None):
    """Determine callable bases for an input BAM in the given region.

    We also identify super high depth regions (7x more than the set maximum depth) to
    avoid calling in since these are repetitive centromere and telomere regions that spike
    memory usage.
    """
    if out_file is None:
        out_file = "%s-callable.bed" % os.path.splitext(data["work_bam"])[0]
    max_depth = utils.get_in(data, ("config", "algorithm", "coverage_depth_max"), 10000)
    depth = {"max": max_depth * 7 if max_depth > 0 else sys.maxint - 1,
             "min": utils.get_in(data, ("config", "algorithm", "coverage_depth_min"), 4)}
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            bam.index(data["work_bam"], data["config"])
            with contextlib.closing(pysam.Samfile(data["work_bam"], "rb")) as in_bam:
                with open(tx_out_file, "w") as out_handle:
                    for r in _regions_for_coverage(data, region, tx_out_file):
                        for chrom, start, end, info in _get_coverage_in_region(in_bam, r, depth):
                            out_handle.write("%s\t%s\t%s\t%s\n" % (chrom, start, end, info))
    return [{"callable_bed": out_file, "config": data["config"], "work_bam": data["work_bam"]}]

def _get_coverage_in_region(in_bam, region, depth):
    """Retrieve summary of coverage in a region.
    Uses chajo's approach of pre-allocating a numpy array for coverage before collapsing
    to cleanly handle non-covered positions. This uses 2Gb memory for human chr1. If memory
    requirements become an issue, could look at splitting long chromosomes at smart places.
    XXX Can replace `positions` with chanjo functionality when it accepts max_depth keyword.
    """
    # special case, do not calculate if we are in a chromosome not covered by BED file
    if region.attrs.get("no_coverage"):
        yield region.chrom, region.start, region.end, "NO_COVERAGE"
    else:
        positions = numpy.zeros(region.end + 1 - region.start)
        for col in in_bam.pileup(str(region.chrom), region.start, region.end + 1, stepper="all",
                                 max_depth=depth["max"] + 10, truncate=True):
            positions[col.pos - region.start] = col.n
        cur_ctype = None
        cur_start = None
        for i, count in enumerate(positions):
            ctype = _get_ctype(count, depth)
            if cur_ctype is None:
                cur_ctype = ctype
                cur_start = i
            elif cur_ctype != ctype:
                yield region.chrom, region.start + cur_start, region.start + i, cur_ctype
                cur_ctype = ctype
                cur_start = i
        if i > cur_start:
            yield region.chrom, region.start + cur_start, region.start + i, cur_ctype

def _get_ctype(count, depth):
    if count == 0:
        return "NO_COVERAGE"
    elif count < depth["min"]:
        return "LOW_COVERAGE"
    elif count > depth["max"]:
        return "EXCESSIVE_COVERAGE"
    else:
        return "CALLABLE"

def _regions_for_coverage(data, region, out_file):
    """Retrieve BedTool iterator over regions we need to calculate coverage in.
    """
    variant_regions = utils.get_in(data, ("config", "algorithm", "variant_regions"))
    ready_region = shared.subset_variant_regions(variant_regions, region, out_file)
    if not ready_region:
        return get_ref_bedtool(data["sam_ref"], data["config"])
    elif os.path.isfile(ready_region):
        return pybedtools.BedTool(ready_region).intervals
    elif isinstance(ready_region, (list, tuple)):
        c, s, e = ready_region
        return [pybedtools.Interval(c, s, e)]
    else:
        assert isinstance(ready_region, basestring)
        out = []
        for r in [x for x in get_ref_bedtool(data["sam_ref"], data["config"])
                  if x.chrom == ready_region]:
            # If we have variant regions but none in this region, don't calculate coverage
            r.attrs["no_coverage"] = variant_regions is not None
            out.append(r)
        return out

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
            if (ctype in ["REF_N", "NO_COVERAGE", "EXCESSIVE_COVERAGE", "LOW_COVERAGE"] and
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
        self._ref_sizes = {x.chrom: x.stop for x in ref_regions}

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
        # fills an entire chromosome, handles smaller decoy and haplotype chromosomes
        elif last_pos == 0 and x.stop >= self._ref_sizes.get(x.chrom, 0) - 100:
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
    return callblock_bed, nblock_bed, callable_bed

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
        parts = ["min: %s" % min(xs),
                 "5%%: %s" % numpy.percentile(xs, 5),
                 "25%%: %s" % numpy.percentile(xs, 25),
                 "median: %s" % numpy.percentile(xs, 50),
                 "75%%: %s" % numpy.percentile(xs, 75),
                 "95%%: %s" % numpy.percentile(xs, 95),
                 "99%%: %s" % numpy.percentile(xs, 99),
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

def _combine_excessive_coverage(samples, ref_regions, min_n_size):
    """Provide a global set of regions with excessive coverage to avoid.
    """
    flag = "EXCESSIVE_COVERAGE"
    ecs = (pybedtools.BedTool(x["regions"]["callable"]).filter(lambda x: x.name == flag)
           for x in samples if "regions" in x)
    merge_ecs = _combine_regions(ecs, ref_regions).saveas()
    if len(merge_ecs) > 0:
        return merge_ecs.merge(d=min_n_size).filter(lambda x: x.stop - x.start > min_n_size).saveas()
    else:
        return merge_ecs

def combine_sample_regions(*samples):
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
        ec_regions = _combine_excessive_coverage(samples, ref_regions, min_n_size)
        block_filter = NBlockRegionPicker(ref_regions, config)
        nblock_size_filtered = nblock_regions.filter(block_filter.include_block).saveas()
        if len(nblock_size_filtered) >= len(ref_regions):
            final_nblock_regions = nblock_size_filtered
        else:
            final_nblock_regions = nblock_regions
        final_regions = ref_regions.subtract(final_nblock_regions)
        if len(ec_regions) > 0:
            final_regions = final_regions.subtract(ec_regions)
        final_regions.merge(d=min_n_size)
        _write_bed_regions(samples[0], final_regions, analysis_file, no_analysis_file)
    else:
        final_regions = pybedtools.BedTool(analysis_file)
    _analysis_block_stats(final_regions)
    regions = {"analysis": [(r.chrom, int(r.start), int(r.stop)) for r in final_regions],
               "noanalysis": no_analysis_file,
               "analysis_bed": analysis_file}
    return [regions]
