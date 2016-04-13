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
import os
import shutil
import subprocess

import numpy
import pybedtools
import pysam
import toolz as tz
import yaml

from bcbio import broad, utils
from bcbio.bam import ref
from bcbio.log import logger
from bcbio.distributed import multi, prun
from bcbio.distributed.split import parallel_split_combine
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils, shared
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import bedutils, realign
from bcbio.variation import multi as vmulti

def parallel_callable_loci(in_bam, ref_file, data):
    config = copy.deepcopy(data["config"])
    num_cores = config["algorithm"].get("num_cores", 1)
    data = {"work_bam": in_bam, "config": config,
            "reference": data["reference"]}
    parallel = {"type": "local", "cores": num_cores, "module": "bcbio.distributed"}
    items = [[data]]
    with prun.start(parallel, items, config, multiplier=int(num_cores)) as runner:
        split_fn = shared.process_bam_by_chromosome("-callable.bed", "work_bam", remove_alts=True)
        out = parallel_split_combine(items, split_fn, runner,
                                     "calc_callable_loci", "combine_bed",
                                     "callable_bed", ["config"])[0]
    return out[0]["callable_bed"]

@multi.zeromq_aware_logging
def calc_callable_loci(data, region=None, out_file=None):
    """Determine callable bases for an input BAM in the given region.
    """
    if out_file is None:
        out_file = "%s-callable.bed" % os.path.splitext(data["work_bam"])[0]
    depth = {"min": dd.get_coverage_depth_min(data)}
    if not utils.file_exists(out_file):
        ref_file = tz.get_in(["reference", "fasta", "base"], data)
        region_file, calc_callable = _regions_for_coverage(data, region, ref_file, out_file)
        if calc_callable:
            _group_by_ctype(_get_coverage_file(data["work_bam"], ref_file, region, region_file, depth,
                                               out_file, data),
                            depth, region, region_file, out_file, data)
        # special case, do not calculate if we are in a chromosome not covered by BED file
        else:
            with file_transaction(data, out_file) as tx_out_file:
                shutil.move(region_file, tx_out_file)
    return [{"callable_bed": out_file, "config": data["config"], "work_bam": data["work_bam"]}]

def _group_by_ctype(bed_file, depth, region, region_file, out_file, data):
    """Group adjacent callable/uncallble regions into defined intervals.

    Uses tips from bedtools discussion:
    https://groups.google.com/d/msg/bedtools-discuss/qYDE6XF-GRA/2icQtUeOX_UJ
    https://gist.github.com/arq5x/b67196a46db5b63bee06
    """
    with file_transaction(data, out_file) as tx_out_file:
        min_cov = depth["min"]
        sort_cmd = bedutils.get_sort_cmd()
        cmd = (r"""cat {bed_file} | awk '{{if ($4 == 0) {{print $0"\tNO_COVERAGE"}} """
               r"""else if ($4 < {min_cov}) {{print $0"\tLOW_COVERAGE"}} """
               r"""else {{print $0"\tCALLABLE"}} }}' | """
               "bedtools groupby -prec 21 -g 1,5 -c 1,2,3,5 -o first,first,max,first | "
               "cut -f 3-6 | "
               "bedtools intersect -nonamecheck -a - -b {region_file} | "
               "{sort_cmd} -k1,1 -k2,2n  > {tx_out_file}")
        do.run(cmd.format(**locals()), "bedtools groupby coverage: %s" % (str(region)), data)

def _get_coverage_file(in_bam, ref_file, region, region_file, depth, base_file, data):
    """Retrieve summary of coverage in a region.
    Requires positive non-zero mapping quality at a position, matching GATK's
    CallableLoci defaults.
    """
    out_file = "%s-genomecov.bed" % utils.splitext_plus(base_file)[0]
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            fai_file = ref.fasta_idx(ref_file, data["config"])
            sambamba = config_utils.get_program("sambamba", data["config"])
            bedtools = config_utils.get_program("bedtools", data["config"])
            cmd = ("{sambamba} view -F 'mapping_quality > 0' -L {region_file} -f bam -l 1 {in_bam} | "
                   "{bedtools} genomecov -split -ibam stdin -bga -g {fai_file} "
                   "> {tx_out_file}")
            do.run(cmd.format(**locals()), "bedtools genomecov: %s" % (str(region)), data)
    # Empty output file, no coverage for the whole contig
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                for feat in get_ref_bedtool(ref_file, data["config"], region):
                    out_handle.write("%s\t%s\t%s\t%s\n" % (feat.chrom, feat.start, feat.end, 0))
    return out_file

def _regions_for_coverage(data, region, ref_file, out_file):
    """Retrieve BED file of regions we need to calculate coverage in.

    Checks for variant region specifications that do not overlap contigs
    (in which case we do not calculate coverage) and regions smaller than
    callable_min_size (in which case we assign everything as callable).
    callable_min_size avoids calculations for small chromosomes we won't
    split on later, saving computation and disk IO.
    """
    variant_regions = dd.get_variant_regions_merged(data)
    ready_region = shared.subset_variant_regions(variant_regions, region, out_file)
    custom_file = "%s-coverageregions.bed" % utils.splitext_plus(out_file)[0]
    region_size = _get_region_size(ref_file, data, region)
    if variant_regions is None and region_size is not None and region_size < dd.get_callable_min_size(data):
        coverage_str = "CALLABLE" if realign.has_aligned_reads(dd.get_work_bam(data), region) else "NO_COVERAGE"
        custom_file = _write_all_chrom_file(coverage_str, custom_file, ref_file, region, data)
        return custom_file, False
    elif not ready_region:
        get_ref_bedtool(ref_file, data["config"]).saveas(custom_file)
        return custom_file, True
    elif os.path.isfile(ready_region):
        return ready_region, True
    elif isinstance(ready_region, (list, tuple)):
        c, s, e = ready_region
        pybedtools.BedTool("%s\t%s\t%s\n" % (c, s, e), from_string=True).saveas(custom_file)
        return custom_file, True
    else:
        custom_file = _write_all_chrom_file("NO_COVERAGE", custom_file, ref_file, region, data)
        return custom_file, variant_regions is None

def _write_all_chrom_file(coverage_str, out_file, ref_file, region, data):
    """Write a coverage calculation for an entire chromosome, skipping work of running bedtools.
    """
    with file_transaction(data, out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            for feat in get_ref_bedtool(ref_file, data["config"], region):
                out_handle.write("%s\t%s\t%s\t%s\n" % (feat.chrom, feat.start, feat.end, coverage_str))
    return out_file

def _get_region_size(ref_file, data, region=None):
    """Retrieve size of a region, potentially returning None if not set.
    """
    if region:
        for contig in ref.file_contigs(ref_file, data["config"]):
            if contig.name == region:
                return contig.size

def sample_callable_bed(bam_file, ref_file, data):
    """Retrieve callable regions for a sample subset by defined analysis regions.
    """
    config = data["config"]
    out_file = "%s-callable_sample.bed" % os.path.splitext(bam_file)[0]
    with shared.bedtools_tmpdir({"config": config}):
        callable_bed = parallel_callable_loci(bam_file, ref_file, data)
        input_regions_bed = config["algorithm"].get("variant_regions", None)
        if not utils.file_uptodate(out_file, callable_bed):
            with file_transaction(config, out_file) as tx_out_file:
                callable_regions = pybedtools.BedTool(callable_bed)
                filter_regions = callable_regions.filter(lambda x: x.name == "CALLABLE")
                if input_regions_bed:
                    if not utils.file_uptodate(out_file, input_regions_bed):
                        input_regions = pybedtools.BedTool(input_regions_bed)
                        filter_regions.intersect(input_regions, nonamecheck=True).saveas(tx_out_file)
                else:
                    filter_regions.saveas(tx_out_file)
    return out_file

def calculate_offtarget(bam_file, ref_file, data):
    """Generate file of offtarget read counts for inputs with variant regions.
    """
    vrs_file = dd.get_variant_regions_merged(data)
    if vrs_file:
        out_file = "%s-offtarget-stats.yaml" % os.path.splitext(bam_file)[0]
        if not utils.file_exists(out_file):
            with file_transaction(data, out_file) as tx_out_file:
                offtarget_regions = "%s-regions.bed" % utils.splitext_plus(out_file)[0]
                ref_bed = get_ref_bedtool(ref_file, data["config"])
                ref_bed.subtract(pybedtools.BedTool(vrs_file), nonamecheck=True).saveas(offtarget_regions)
                samtools = config_utils.get_program("samtools", data["config"])
                bedtools = config_utils.get_program("bedtools", data["config"])
                cmd = ("{samtools} view -u {bam_file} -L {offtarget_regions} | "
                       "{bedtools} intersect -abam - -b {offtarget_regions} -f 1.0 -bed | wc -l")
                offtarget_count = int(subprocess.check_output(cmd.format(**locals()), shell=True))
                cmd = "{samtools} idxstats {bam_file} | awk '{{s+=$3}} END {{print s}}'"
                mapped_count = int(subprocess.check_output(cmd.format(**locals()), shell=True))
                with open(tx_out_file, "w") as out_handle:
                    yaml.safe_dump({"mapped": mapped_count, "offtarget": offtarget_count}, out_handle,
                                   allow_unicode=False, default_flow_style=False)
        return out_file

def get_ref_bedtool(ref_file, config, chrom=None):
    """Retrieve a pybedtool BedTool object with reference sizes from input reference.
    """
    broad_runner = broad.runner_from_path("picard", config)
    ref_dict = broad_runner.run_fn("picard_index_ref", ref_file)
    ref_lines = []
    with contextlib.closing(pysam.Samfile(ref_dict, "r")) as ref_sam:
        for sq in ref_sam.header["SQ"]:
            if not chrom or sq["SN"] == chrom:
                ref_lines.append("%s\t%s\t%s" % (sq["SN"], 0, sq["LN"]))
    return pybedtools.BedTool("\n".join(ref_lines), from_string=True)

def _get_nblock_regions(in_file, min_n_size, ref_regions):
    """Retrieve coordinates of regions in reference genome with no mapping.
    These are potential breakpoints for parallelizing analysis.
    """
    out_lines = []
    called_contigs = set([])
    with open(in_file) as in_handle:
        for line in in_handle:
            contig, start, end, ctype = line.rstrip().split()
            called_contigs.add(contig)
            if (ctype in ["REF_N", "NO_COVERAGE", "EXCESSIVE_COVERAGE", "LOW_COVERAGE"] and
                  int(end) - int(start) > min_n_size):
                out_lines.append("%s\t%s\t%s\n" % (contig, start, end))
    for refr in ref_regions:
        if refr.chrom not in called_contigs:
            out_lines.append("%s\t%s\t%s\n" % (refr.chrom, 0, refr.stop))
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
        input_nblock = ref_regions.subtract(input_regions, nonamecheck=True)
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
    def __init__(self, ref_regions, config, min_n_size):
        self._end_buffer = 250 if min_n_size > 50 else 0
        self._chr_last_blocks = {}
        target_blocks = int(config["algorithm"].get("nomap_split_targets", 200))
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
        # Region excludes an entire chromosome, typically decoy/haplotypes
        if last_pos <= self._end_buffer and x.stop >= self._ref_sizes.get(x.chrom, 0) - self._end_buffer:
            return True
        # Do not split on smaller decoy and haplotype chromosomes
        elif self._ref_sizes.get(x.chrom, 0) <= self._target_size:
            return False
        elif (x.start - last_pos) > self._target_size:
            self._chr_last_blocks[x.chrom] = x.stop
            return True
        else:
            return False

    def expand_block(self, feat):
        """Expand any blocks which are near the start or end of a contig.
        """
        chrom_end = self._ref_sizes.get(feat.chrom)
        if chrom_end:
            if feat.start < self._end_buffer:
                feat.start = 0
            if feat.stop >= chrom_end - self._end_buffer:
                feat.stop = chrom_end
        return feat

def block_regions(in_bam, ref_file, data):
    """Find blocks of regions for analysis from mapped input BAM file.

    Identifies islands of callable regions, surrounding by regions
    with no read support, that can be analyzed independently.
    """
    config = data["config"]
    min_n_size = int(config["algorithm"].get("nomap_split_size", 250))
    with shared.bedtools_tmpdir({"config": config}):
        callable_bed = parallel_callable_loci(in_bam, ref_file, data)
        nblock_bed = "%s-nblocks%s" % os.path.splitext(callable_bed)
        callblock_bed = "%s-callableblocks%s" % os.path.splitext(callable_bed)
        if not utils.file_uptodate(nblock_bed, callable_bed):
            ref_regions = get_ref_bedtool(ref_file, config)
            nblock_regions = _get_nblock_regions(callable_bed, min_n_size, ref_regions)
            nblock_regions = _add_config_regions(nblock_regions, ref_regions, config)
            nblock_regions.filter(lambda r: len(r) > min_n_size).saveas(nblock_bed)
            if len(ref_regions.subtract(nblock_regions, nonamecheck=True)) > 0:
                ref_regions.subtract(nblock_bed, nonamecheck=True).merge(d=min_n_size).saveas(callblock_bed)
            else:
                raise ValueError("No callable regions found from BAM file. Alignment regions might "
                                 "not overlap with regions found in your `variant_regions` BED: %s" % in_bam)
    return callblock_bed, nblock_bed, callable_bed

def _write_bed_regions(data, final_regions, out_file, out_file_ref):
    ref_file = tz.get_in(["reference", "fasta", "base"], data)
    ref_regions = get_ref_bedtool(ref_file, data["config"])
    noanalysis_regions = ref_regions.subtract(final_regions, nonamecheck=True)
    final_regions.saveas(out_file)
    noanalysis_regions.saveas(out_file_ref)

def _analysis_block_stats(regions, samples):
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
        raise ValueError("No callable regions found in: %s" %
                         (", ".join([dd.get_sample_name(x) for x in samples])))

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

def combine_sample_regions(*samples):
    """Create batch-level sets of callable regions for multi-sample calling.

    Intersects all non-callable (nblock) regions from all samples in a batch,
    producing a global set of callable regions.
    """
    samples = utils.unpack_worlds(samples)
    # back compatibility -- global file for entire sample set
    global_analysis_file = os.path.join(samples[0]["dirs"]["work"], "analysis_blocks.bed")
    if utils.file_exists(global_analysis_file) and not _needs_region_update(global_analysis_file, samples):
        global_no_analysis_file = os.path.join(os.path.dirname(global_analysis_file), "noanalysis_blocks.bed")
    else:
        global_analysis_file = None
    out = []
    analysis_files = []
    batches = []
    with shared.bedtools_tmpdir(samples[0]):
        for batch, items in vmulti.group_by_batch(samples, require_bam=False).items():
            batches.append(items)
            if global_analysis_file:
                analysis_file, no_analysis_file = global_analysis_file, global_no_analysis_file
            else:
                analysis_file, no_analysis_file = _combine_sample_regions_batch(batch, items)
            for data in items:
                vr_file = dd.get_variant_regions(data)
                if analysis_file:
                    analysis_files.append(analysis_file)
                    data["config"]["algorithm"]["callable_regions"] = analysis_file
                    data["config"]["algorithm"]["non_callable_regions"] = no_analysis_file
                    data["config"]["algorithm"]["callable_count"] = pybedtools.BedTool(analysis_file).count()
                elif vr_file:
                    data["config"]["algorithm"]["callable_count"] = pybedtools.BedTool(vr_file).count()
                highdepth_bed = tz.get_in(["regions", "highdepth"], data)
                if highdepth_bed:
                    data["config"]["algorithm"]["highdepth_regions"] = highdepth_bed
                # attach a representative sample for calculating callable region
                if not data.get("work_bam"):
                    for x in items:
                        if x.get("work_bam"):
                            data["work_bam_callable"] = x["work_bam"]
                out.append([data])
        assert len(out) == len(samples)
        if len(analysis_files) > 0:
            final_regions = pybedtools.BedTool(analysis_files[0])
            _analysis_block_stats(final_regions, batches[0])
    return out

def _combine_sample_regions_batch(batch, items):
    """Combine sample regions within a group of batched samples.
    """
    config = items[0]["config"]
    work_dir = utils.safe_makedir(os.path.join(items[0]["dirs"]["work"], "regions"))
    analysis_file = os.path.join(work_dir, "%s-analysis_blocks.bed" % batch)
    no_analysis_file = os.path.join(work_dir, "%s-noanalysis_blocks.bed" % batch)
    if not utils.file_exists(analysis_file) or _needs_region_update(analysis_file, items):
        # Combine all nblocks into a final set of intersecting regions
        # without callable bases. HT @brentp for intersection approach
        # https://groups.google.com/forum/?fromgroups#!topic/bedtools-discuss/qA9wK4zN8do
        bed_regions = [pybedtools.BedTool(x["regions"]["nblock"])
                       for x in items if "regions" in x]
        if len(bed_regions) == 0:
            analysis_file, no_analysis_file = None, None
        else:
            with file_transaction(items[0], analysis_file, no_analysis_file) as (tx_afile, tx_noafile):
                def intersect_two(a, b):
                    return a.intersect(b, nonamecheck=True)
                nblock_regions = reduce(intersect_two, bed_regions).saveas(
                    "%s-nblock%s" % utils.splitext_plus(tx_afile))
                ref_file = tz.get_in(["reference", "fasta", "base"], items[0])
                ref_regions = get_ref_bedtool(ref_file, config)
                min_n_size = int(config["algorithm"].get("nomap_split_size", 250))
                block_filter = NBlockRegionPicker(ref_regions, config, min_n_size)
                final_nblock_regions = nblock_regions.filter(
                    block_filter.include_block).saveas().each(block_filter.expand_block).saveas(
                        "%s-nblockfinal%s" % utils.splitext_plus(tx_afile))
                final_regions = ref_regions.subtract(final_nblock_regions, nonamecheck=True).merge(d=min_n_size)
                _write_bed_regions(items[0], final_regions, tx_afile, tx_noafile)
    if analysis_file and utils.file_exists(analysis_file):
        return analysis_file, no_analysis_file
    else:
        return None, None
