"""Support for Copy Number Variations (CNVs) with GATK4

https://software.broadinstitute.org/gatk/documentation/article?id=11682
https://gatkforums.broadinstitute.org/dsde/discussion/11683/
"""
import os

import toolz as tz

from bcbio import broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.variation import bedutils

def run(items, background=None):
    """Detect copy number variations from batched set of samples using GATK4 CNV calling.
    """
    if not background: background = []
    return items

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
                      "-O", tx_out_file]
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
    normal_counts = _run_collect_allelic_counts(cur_het_bed, key, work_dir, paired.normal_data)
    return tumor_counts, normal_counts

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

def _run_with_memory_scaling(params, tx_out_file, data):
    num_cores = dd.get_num_cores(data)
    memscale = {"magnitude": 0.9 * num_cores, "direction": "increase"} if num_cores > 1 else None
    broad_runner = broad.runner_from_config(data["config"])
    broad_runner.run_gatk(params, os.path.dirname(tx_out_file), memscale=memscale)
