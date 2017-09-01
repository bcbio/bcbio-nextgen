"""Examine and query coverage in sequencing experiments.

Provides estimates of coverage intervals based on callable regions
"""
import csv
import itertools
import os
import shutil
import yaml

import pybedtools
import numpy as np
import pysam

from bcbio.variation.bedutils import clean_file
from bcbio.utils import (file_exists, append_stem, copy_plus)
from bcbio import bam, utils
from bcbio.bam import ref, sambamba
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.pipeline import shared

GENOME_COV_THRESH = 0.40  # percent of genome covered for whole genome analysis
OFFTARGET_THRESH = 0.01  # percent of offtarget reads required to be capture (not amplification) based

def assign_interval(data):
    """Identify coverage based on percent of genome covered and relation to targets.

    Classifies coverage into 3 categories:
      - genome: Full genome coverage
      - regional: Regional coverage, like exome capture, with off-target reads
      - amplicon: Amplication based regional coverage without off-target reads
    """
    if not dd.get_coverage_interval(data):
        vrs = dd.get_variant_regions_merged(data)
        callable_file = dd.get_sample_callable(data)
        if vrs:
            callable_size = pybedtools.BedTool(vrs).total_coverage()
        else:
            callable_size = pybedtools.BedTool(callable_file).total_coverage()
        total_size = sum([c.size for c in ref.file_contigs(dd.get_ref_file(data), data["config"])])
        genome_cov_pct = callable_size / float(total_size)
        if genome_cov_pct > GENOME_COV_THRESH:
            cov_interval = "genome"
            offtarget_pct = 0.0
        elif not vrs:
            cov_interval = "regional"
            offtarget_pct = 0.0
        else:
            offtarget_pct = _count_offtarget(data, dd.get_align_bam(data) or dd.get_work_bam(data),
                                             vrs or callable_file, "variant_regions")
            if offtarget_pct > OFFTARGET_THRESH:
                cov_interval = "regional"
            else:
                cov_interval = "amplicon"
        logger.info("%s: Assigned coverage as '%s' with %.1f%% genome coverage and %.1f%% offtarget coverage"
                    % (dd.get_sample_name(data), cov_interval, genome_cov_pct * 100.0, offtarget_pct * 100.0))
        data["config"]["algorithm"]["coverage_interval"] = cov_interval
    return data

def _count_offtarget(data, bam_file, bed_file, target_name):
    mapped_unique = sambamba.number_of_mapped_reads(data, bam_file, keep_dups=False)
    ontarget = sambamba.number_of_mapped_reads(
        data, bam_file, keep_dups=False, bed_file=bed_file, target_name=target_name)
    if mapped_unique:
        return float(mapped_unique - ontarget) / mapped_unique
    else:
        return 0.0

def calculate(bam_file, data):
    """Calculate coverage in parallel using samtools depth through goleft.

    samtools depth removes duplicates and secondary reads from the counts:
    if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
    """
    params = {"window_size": 5000, "parallel_window_size": 1e5, "min": dd.get_coverage_depth_min(data)}
    prefix = os.path.join(
        utils.safe_makedir(os.path.join(dd.get_work_dir(data), "align", dd.get_sample_name(data))),
        "%s-coverage" % (dd.get_sample_name(data)))
    depth_file = prefix + ".depth.bed"
    callable_file = prefix + ".callable.bed"
    variant_regions = dd.get_variant_regions_merged(data)
    if not utils.file_uptodate(callable_file, bam_file):
        cmd = ["goleft", "depth", "--q", "1", "--mincov", str(params["min"]),
               "--processes", str(dd.get_num_cores(data)), "--ordered"]
        with file_transaction(data, depth_file) as tx_depth_file:
            with utils.chdir(os.path.dirname(tx_depth_file)):
                tx_callable_file = tx_depth_file.replace(".depth.bed", ".callable.bed")
                prefix = tx_depth_file.replace(".depth.bed", "")
                bam_ref_file = "%s-bamref.fa" % utils.splitext_plus(bam_file)[0]
                bam.fai_from_bam(dd.get_ref_file(data), bam_file, bam_ref_file + ".fai", data)
                cmd += ["--reference", bam_ref_file]
                cmd += ["--prefix", prefix, bam_file]
                bcbio_env = utils.get_bcbio_env()
                msg = "Calculate coverage: %s" % dd.get_sample_name(data)
                do.run(cmd, msg, env=bcbio_env)
                shutil.move(tx_callable_file, callable_file)
    final_callable = _subset_to_variant_regions(callable_file, variant_regions, data)
    return depth_file, final_callable

def _create_genome_regions(callable_file, data):
    """Create whole genome contigs we want to process, only non-alts.

    Skips problem contigs like HLAs for downstream analysis.
    """
    variant_regions = "%s-genome.bed" % utils.splitext_plus(callable_file)[0]
    with file_transaction(data, variant_regions) as tx_variant_regions:
        with open(tx_variant_regions, "w") as out_handle:
            for c in shared.get_noalt_contigs(data):
                out_handle.write("%s\t%s\t%s\n" % (c.name, 0, c.size))
    return variant_regions

def _subset_to_variant_regions(callable_file, variant_regions, data):
    """Subset output callable file to only variant regions of interest.
    """
    out_file = "%s-vrsubset.bed" % utils.splitext_plus(callable_file)[0]
    if not utils.file_uptodate(out_file, callable_file):
        if not variant_regions:
            variant_regions = _create_genome_regions(callable_file, data)
        with file_transaction(data, out_file) as tx_out_file:
            pybedtools.BedTool(callable_file).intersect(variant_regions).saveas(tx_out_file)
    return out_file

def _get_cache_file(data, target_name):
    prefix = os.path.join(
        utils.safe_makedir(os.path.join(dd.get_work_dir(data), "align", dd.get_sample_name(data))),
        "%s-coverage" % (dd.get_sample_name(data)))
    cache_file = prefix + "-" + target_name + "-stats.yaml"
    return cache_file

def _read_cache(cache_file, reuse_cmp_files):
    reuse_cmp_file = [fn for fn in reuse_cmp_files if fn]
    if all(utils.file_uptodate(cache_file, fn) for fn in reuse_cmp_file):
        with open(cache_file) as in_handle:
            return yaml.safe_load(in_handle)
    return dict()

def _write_cache(cache, cache_file):
    with open(cache_file, "w") as out_handle:
        yaml.safe_dump(cache, out_handle, default_flow_style=False, allow_unicode=False)

def get_average_coverage(target_name, bed_file, data, bam_file=None):
    if not bam_file:
        bam_file = dd.get_align_bam(data) or dd.get_work_bam(data)
    cache_file = _get_cache_file(data, target_name)
    cache = _read_cache(cache_file, [bam_file, bed_file])
    if "avg_coverage" in cache:
        return int(cache["avg_coverage"])

    if bed_file:
        avg_cov = _average_bed_coverage(bed_file, data)
    else:
        avg_cov = _average_genome_coverage(data, bam_file)

    cache["avg_coverage"] = int(avg_cov)
    _write_cache(cache, cache_file)
    return int(avg_cov)

def _average_genome_coverage(data, bam_file):
    """Quickly calculate average coverage for whole genome files using indices.

    Includes all reads, with duplicates.
    """
    total = sum([c.size for c in ref.file_contigs(dd.get_ref_file(data), data["config"])])
    read_counts = sum(x.aligned for x in bam.idxstats(bam_file, data))
    with pysam.Samfile(bam_file, "rb") as pysam_bam:
        read_size = np.median(list(itertools.islice((a.query_length for a in pysam_bam.fetch()), 1e5)))
    avg_cov = float(read_counts * read_size) / total
    return avg_cov

def _average_bed_coverage(bed_file, data):
    sambamba_depth_file = regions_coverage(bed_file, data)
    avg_covs = []
    total_len = 0
    with open(sambamba_depth_file) as fh:
        for line_tokens in (l.rstrip().split() for l in fh if not l.startswith("#")):
            line_tokens = [x for x in line_tokens if x.strip()]
            start, end = map(int, line_tokens[1:3])
            size = end - start
            avg_covs.append(float(line_tokens[-1]) * size)
            total_len += size
    avg_cov = sum(avg_covs) / total_len if total_len > 0 else 0
    return avg_cov

def _calculate_percentiles(cov_file, dist_file, cutoffs, out_dir, data):
    """Calculate percentage over over specified cutoff range.

    XXX Does not calculate the per-bin coverage estimations which we had
    earlier with sambamba depth. Instead has a global metric of percent coverage
    which provides a more defined look at coverage changes by depth.
    """
    if not utils.file_exists(dist_file) or not utils.file_exists(cov_file):
        return []
    sample = dd.get_sample_name(data)
    out_total_file = append_stem(dist_file, "_total_summary")
    if not utils.file_exists(out_total_file):
        with file_transaction(data, out_total_file) as tx_file:
            with open(tx_file, 'w') as out_handle:
                writer = csv.writer(out_handle, dialect="excel-tab")
                writer.writerow(["cutoff_reads", "bases_pct", "sample"])
                with open(dist_file) as in_handle:
                    for line in in_handle:
                        count, pct = line.strip().split()
                        count = int(count)
                        pct = "%.1f" % (float(pct) * 100.0)
                        if count >= min(cutoffs) and count <= max(cutoffs):
                            writer.writerow(["percentage%s" % count, pct, sample])
                    if min(cutoffs) < count:
                        writer.writerow(["percentage%s" % min(cutoffs), pct, sample])
    # To move metrics to multiqc, will remove older files
    # when bcbreport accepts these one, to avoid errors
    # while porting everything to multiqc
    # These files will be copied to final
    out_total_fixed = os.path.join(os.path.dirname(out_total_file), "%s_bcbio_coverage_avg.txt" % sample)
    copy_plus(out_total_file, out_total_fixed)
    return [out_total_fixed]

def regions_coverage(bed_file, data):
    """Generate coverage over regions of interest using mosdepth.
    """
    cov_file, _ = _run_mosdepth(bed_file, data)
    return cov_file

def _run_mosdepth(bed_file, data):
    """Run mosdepth for a specific BED file generating coverage and distribution.
    """
    bam_file = dd.get_align_bam(data) or dd.get_work_bam(data)
    work_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "coverage", dd.get_sample_name(data)))
    target_name = utils.splitext_plus(os.path.basename(bed_file))[0]
    out_file = os.path.join(work_dir, "%s-coverage.bed" % target_name)
    dist_file = os.path.join(work_dir, "%s-distribution.txt" % target_name)
    if not utils.file_uptodate(out_file, bam_file) or not utils.file_uptodate(out_file, bed_file):
        with file_transaction(data, out_file) as tx_out_file:
            with file_transaction(data, dist_file) as tx_dist_file:
                num_cores = dd.get_cores(data)
                cmd = ("mosdepth -t {num_cores} -F 1804 -b {bed_file} -d {tx_dist_file} {bam_file} > {tx_out_file}")
                message = "Calculating regions coverage of {target_name} in {bam_file}"
                do.run(cmd.format(**locals()), message.format(**locals()))
    return out_file, dist_file

def coverage_region_detailed_stats(data, out_dir, extra_cutoffs=None):
    """
    Calculate coverage at different completeness cutoff
    for region in coverage option.
    """
    bed_file = dd.get_coverage(data)
    if not bed_file or not utils.file_exists(bed_file):
        return []
    else:
        bed_file = clean_file(bed_file, data, prefix="cov-", simple=True)
        cov_file, dist_file = _run_mosdepth(bed_file, data)
        cutoffs = {1, 5, 10, 20, 50, 100, 250, 500, 1000, 5000, 10000, 50000}
        if extra_cutoffs:
            cutoffs = sorted(list(cutoffs | extra_cutoffs))
        out_files = _calculate_percentiles(cov_file, dist_file, cutoffs, out_dir, data)
        return [os.path.abspath(x) for x in out_files]
