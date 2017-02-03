"""Examine and query coverage in sequencing experiments.

Provides estimates of coverage intervals based on callable regions
"""
import itertools
import os
import shutil
import yaml

import pybedtools
import pandas as pd
import numpy as np
import pysam

from bcbio.variation.bedutils import clean_file
from bcbio.utils import (file_exists, chdir, safe_makedir,
                         append_stem, copy_plus)
from bcbio import utils
from bcbio.bam import ref, sambamba
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.pipeline import shared

def assign_interval(data):
    """Identify coverage based on percent of genome covered and relation to targets.

    Classifies coverage into 3 categories:
      - genome: Full genome coverage
      - regional: Regional coverage, like exome capture, with off-target reads
      - amplicon: Amplication based regional coverage without off-target reads
    """
    genome_cov_thresh = 0.40  # percent of genome covered for whole genome analysis
    offtarget_thresh = 0.01  # percent of offtarget reads required to be capture (not amplification) based
    if not dd.get_coverage_interval(data):
        vrs = dd.get_variant_regions_merged(data)
        callable_file = dd.get_sample_callable(data)
        if vrs:
            callable_size = pybedtools.BedTool(vrs).total_coverage()
        else:
            callable_size = pybedtools.BedTool(callable_file).total_coverage()
        total_size = sum([c.size for c in ref.file_contigs(dd.get_ref_file(data), data["config"])])
        genome_cov_pct = callable_size / float(total_size)
        if genome_cov_pct > genome_cov_thresh:
            cov_interval = "genome"
            offtarget_pct = 0.0
        elif not vrs:
            cov_interval = "regional"
            offtarget_pct = 0.0
        else:
            offtarget_pct = _count_offtarget(data, dd.get_align_bam(data) or dd.get_work_bam(data),
                                             vrs or callable_file, "variant_regions")
            if offtarget_pct > offtarget_thresh:
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
    return 0.0

def calculate(bam_file, data):
    """Calculate coverage in parallel using samtools depth through goleft.

    samtools depth removes duplicates and secondary reads from the counts:
    if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
    """
    params = {"window_size": 5000, "parallel_window_size": 1e5, "min": dd.get_coverage_depth_min(data),
              "high_multiplier": 20}
    prefix = os.path.join(
        utils.safe_makedir(os.path.join(dd.get_work_dir(data), "align", dd.get_sample_name(data))),
        "%s-coverage" % (dd.get_sample_name(data)))
    depth_file = prefix + ".depth.bed"
    callable_file = prefix + ".callable.bed"
    variant_regions = dd.get_variant_regions_merged(data)
    variant_regions_avg_cov = get_average_coverage(data, bam_file, variant_regions, "variant_regions")
    if not utils.file_uptodate(callable_file, bam_file):
        ref_file = dd.get_ref_file(data)
        cmd = ["goleft", "depth", "--q", "1",
               "--mincov", str(params["min"]), "--reference", ref_file,
               "--processes", str(dd.get_num_cores(data)), "--ordered"]
        max_depth = _get_max_depth(variant_regions_avg_cov, params, data)
        if max_depth:
            cmd += ["--maxmeandepth", str(int(max_depth))]
        with file_transaction(data, depth_file) as tx_depth_file:
            with utils.chdir(os.path.dirname(tx_depth_file)):
                tx_callable_file = tx_depth_file.replace(".depth.bed", ".callable.bed")
                prefix = tx_depth_file.replace(".depth.bed", "")
                cmd += ["--prefix", prefix, bam_file]
                bcbio_env = utils.get_bcbio_env()
                msg = "Calculate coverage: %s" % dd.get_sample_name(data)
                do.run(cmd, msg, env=bcbio_env)
                shutil.move(tx_callable_file, callable_file)
    final_callable = _subset_to_variant_regions(callable_file, variant_regions, data)
    return depth_file, final_callable, _extract_highdepth(final_callable, data), variant_regions_avg_cov

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

def _extract_highdepth(callable_file, data):
    out_file = "%s-highdepth.bed" % utils.splitext_plus(callable_file)[0]
    if not utils.file_uptodate(out_file, callable_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(callable_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    for line in in_handle:
                        parts = line.strip().split("\t")
                        if "EXCESSIVE_COVERAGE" in parts:
                            out_handle.write("\t".join(parts[:3] + ["highdepth"]) + "\n")
    return out_file

def _get_max_depth(average_coverage, params, data):
    """Calculate maximum depth based on a rough multiplier of average coverage.
    """
    if dd.get_coverage_interval(data) == "genome":
        avg_cov = min(30.0, average_coverage)
        return avg_cov * params["high_multiplier"]

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

def get_average_coverage(data, bam_file, bed_file=None, target_name="genome"):
    cache_file = _get_cache_file(data, target_name)
    cache = _read_cache(cache_file, [bam_file, bed_file])
    if "avg_coverage" in cache:
        return cache["avg_coverage"]

    if bed_file:
        avg_cov = _average_bed_coverage(data, bed_file, bam_file, target_name=target_name)
    else:
        avg_cov = _average_genome_coverage(data, bam_file)

    cache["avg_coverage"] = avg_cov
    _write_cache(cache, cache_file)
    return avg_cov

def _average_genome_coverage(data, bam_file):
    total = sum([c.size for c in ref.file_contigs(dd.get_ref_file(data), data["config"])])
    read_counts = sambamba.number_of_mapped_reads(data, bam_file, keep_dups=False)
    with pysam.Samfile(bam_file, "rb") as pysam_bam:
        read_size = np.median(list(itertools.islice((a.query_length for a in pysam_bam.fetch()), 1e5)))
    avg_cov = float(read_counts * read_size) / total
    return avg_cov

def _average_bed_coverage(data, bed_file, bam_file, target_name):
    sambamba_depth_file = regions_coverage(data, bed_file, bam_file, target_name)
    avg_covs = []
    mean_cov_col = None
    total_len = 0
    with open(sambamba_depth_file) as fh:
        for line in fh:
            if line.startswith('#'):
                mean_cov_col = line.split('\t').index('meanCoverage')
                continue
            line_tokens = line.replace('\n', '').split()
            start, end = map(int, line_tokens[1:3])
            size = end - start
            avg_covs.append(float(line_tokens[mean_cov_col]) * size)
            total_len += size
    avg_cov = sum(avg_covs) / total_len if total_len > 0 else 0
    return avg_cov

def checkpoint(stem):
    def check_file(f):
        def wrapper(*args, **kwargs):
            out_file = append_stem(args[0], stem)
            if file_exists(out_file):
                logger.debug("Skipping %s" % out_file)
                return out_file
            return f(*args, **kwargs)
        return wrapper
    return check_file

def _calculate_percentiles(in_file, sample, data=None, cutoffs=None):
    """
    Parse pct bases per region to summarize it in
    7 different pct of regions points with pct bases covered
    higher than a completeness cutoff (5, 10, 20, 50 ...)
    """
    has_data = False
    with open(in_file) as in_handle:
        for i, line in enumerate(in_handle):
            if i > 0:
                has_data = True
                break
    if not has_data:
        return []
    out_file = append_stem(in_file, "_summary")
    out_total_file = append_stem(in_file, "_total_summary")
    if not utils.file_exists(out_file) or not utils.file_exists(out_total_file):
        dt = pd.read_csv(in_file, sep="\t", index_col=False)
        pct = dict()
        pct_bases = dict()
        size = np.array(dt["chromEnd"]) - np.array(dt["chromStart"])
        for cutoff in [h for h in list(dt) if h.startswith("percentage")]:
            if cutoffs and int(cutoff.split("percentage")[1]) in cutoffs:
                a = np.array(dt[cutoff])
                for p_point in [0.01, 10, 25, 50, 75, 90, 99.9]:
                    q = np.percentile(a, p_point)
                    pct[(cutoff, p_point)] = q
                pct_bases[cutoff] = sum(size * a) / float(sum(size))

        with file_transaction(data, out_total_file) as tx_file:
            with open(tx_file, 'w') as out_handle:
                print >>out_handle, "cutoff_reads\tbases_pct\tsample"
                for k in pct_bases:
                    print >>out_handle, "\t".join(map(str, [k, pct_bases[k], sample]))
        with file_transaction(data, out_file) as tx_file:
            with open(tx_file, 'w') as out_handle:
                print >>out_handle, "cutoff_reads\tregion_pct\tbases_pct\tsample"
                for k in pct:
                    print >>out_handle, "\t".join(map(str, [k[0], k[1], pct[k], sample]))
    # To move metrics to multiqc, will remove older files
    # when bcbreport accepts these one, to avoid errors
    # while porting everything to multiqc
    # These files will be copied to final
    out_file_fixed = os.path.join(os.path.dirname(out_file), "%s_bcbio_coverage.txt" % sample)
    out_total_fixed = os.path.join(os.path.dirname(out_file), "%s_bcbio_coverage_avg.txt" % sample)
    copy_plus(out_file, out_file_fixed)
    copy_plus(out_total_file, out_total_fixed)
    return [out_file_fixed, out_total_fixed]

def _read_regions(fn):
    """
    Save in a dict the position of regions with
    the information of the coverage stats.
    """
    regions = {}
    with open(fn) as in_handle:
        for line in in_handle:
            if line.startswith("chrom"):
                regions["header"] = line.strip()
                continue
            idx = "".join(line.split("\t")[:2])
            regions[idx] = line.strip()
    return regions

@checkpoint("_fixed")
def _add_high_covered_regions(in_file, bed_file, sample, data=None):
    """
    Add regions with higher coverage than the limit
    as fully covered.
    """
    out_file = append_stem(in_file, "_fixed")
    regions = _read_regions(in_file)
    with file_transaction(data, out_file) as out_tx:
        with open(bed_file) as in_handle:
            with open(out_tx, 'w') as out_handle:
                if "header" in regions:
                    print >>out_handle, regions["header"]
                for line in in_handle:
                    idx = "".join(line.split("\t")[:2])
                    if idx not in regions:
                        print >>out_handle, "%s\t1000\t1000\t100\t100\t100\t100\t100\t100\t100\t100\t100\t100\t%s" % (line.strip(), sample)
                    else:
                        print >>out_handle, regions[idx]
    return out_file

def _summary_variants(in_file, out_file, data=None):
    """Parse GC and depth variant file
       to be ready for multiqc.
    """
    dt = pd.read_csv(in_file, sep="\t", index_col=False,
                     dtype={"CG": np.float64, "depth": np.float64}, na_values=["."]).dropna()
    row = list()
    with file_transaction(data, out_file) as out_tx:
        cg = dt["CG"]
        d = dt["depth"]
        for p_point in [0.01, 10, 25, 50, 75, 90, 99.9, 100]:
            if len(cg) > 0:
                q_cg = np.percentile(cg, p_point)
            else:
                q_cg = 0
            if len(d) > 0:
                q_d = np.percentile(d, p_point)
            else:
                q_d = 0
            row.append([p_point, q_d, q_cg])
        pd.DataFrame(row).to_csv(out_tx, header=["pct_variants", "depth", "cg"], index=False, sep="\t")

def regions_coverage(data, bed_file, bam_file, target_name, depth_thresholds=None):
    work_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "coverage", dd.get_sample_name(data)))
    out_file = os.path.join(work_dir, target_name + "_regions_depth.bed")
    if utils.file_uptodate(out_file, bam_file) and utils.file_uptodate(out_file, bed_file):
        return out_file
    with file_transaction(data, out_file) as tx_out_file:
        cmdl = sambamba.make_command(data, "depth region", bam_file, bed_file, depth_thresholds=depth_thresholds)
        cmdl += " -o " + tx_out_file
        message = "Calculating regions coverage of {target_name} in {bam_file}"
        do.run(cmdl, message.format(**locals()))
    return out_file

def coverage_region_detailed_stats(data, out_dir, extra_cutoffs=None):
    """
    Calculate coverage at different completeness cutoff
    for region in coverage option.
    """
    bed_file = dd.get_coverage(data)
    if not bed_file or not utils.file_exists(bed_file):
        return []
    work_dir = safe_makedir(out_dir)
    cleaned_bed = clean_file(bed_file, data, prefix="cov-", simple=True)

    cutoffs = {1, 5, 10, 20, 50, 100, 250, 500, 1000, 5000, 10000, 50000}

    with chdir(work_dir):
        in_bam = dd.get_align_bam(data) or dd.get_work_bam(data)
        sample = dd.get_sample_name(data)
        logger.debug("doing coverage for %s" % sample)
        parse_file = os.path.join(sample + "_coverage.bed")
        if utils.file_uptodate(parse_file, cleaned_bed) and utils.file_uptodate(parse_file, in_bam):
            pass
        else:
            with file_transaction(data, parse_file) as out_tx:
                depth_thresholds = sorted(list(cutoffs | extra_cutoffs))
                cmdl = sambamba.make_command(data, "depth region", in_bam, cleaned_bed, depth_thresholds=depth_thresholds)
                cmdl += " | sed 's/# chrom/chrom/' > " + out_tx
                do.run(cmdl, "Run coverage regional analysis for {}".format(sample))
        out_files = _calculate_percentiles(os.path.abspath(parse_file), sample, data=data, cutoffs=cutoffs)
    return [os.path.abspath(x) for x in out_files]
