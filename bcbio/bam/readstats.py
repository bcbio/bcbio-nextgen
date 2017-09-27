"""Calculation of mapped reads by BAM counting, currently implemented with samtools.
"""
import json
import os
import toolz as tz

from bcbio import bam, utils
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.distributed.transaction import tx_tmpdir

pybedtools = utils.LazyImport("pybedtools")

def _backcompatible_cache_file(query_flags, bed_file, target_name, data):
    """Back-compatible: retrieve cache file from previous location.
    """
    cmd_id = "num_" + " and ".join(query_flags).replace(" ", "_")
    if bed_file is not None:
        target_name = target_name or os.path.basename(bed_file)
        cmd_id += "_on_" + target_name
    work_dir = os.path.join(dd.get_work_dir(data), "coverage", dd.get_sample_name(data), "sambamba")
    output_file = os.path.join(work_dir, cmd_id)
    if utils.file_exists(output_file):
        return output_file

def get_cache_file(data):
    cache_file = tz.get_in(["regions", "mapped_stats"], data)
    if cache_file:
        return cache_file
    else:
        return os.path.join(utils.safe_makedir(os.path.join(dd.get_work_dir(data), "coverage")),
                            "mapped_stats.txt")

def number_of_mapped_reads(data, bam_file, keep_dups=True, bed_file=None, target_name=None):
    """Count mapped reads, allow adjustment for duplicates and BED regions.

    Since samtools view does not use indexes for BED files
    (https://github.com/samtools/samtools/issues/88)
    we loop over regions in a BED file and add the counts together.

    Uses a global cache file to store counts, making it possible to pass this single
    file for CWL runs. For parallel processes it can have concurrent append writes,
    with the hopes this will avoid race conditions. Need to revisit with some kind
    of global file locking if it becomes an issue.
    """
    # Flag explainer https://broadinstitute.github.io/picard/explain-flags.html
    callable_flags = ["not unmapped", "not mate_is_unmapped", "not secondary_alignment",
                      "not failed_quality_control"]
    if keep_dups:
        query_flags = callable_flags
        flag = 780  # not (read unmapped or mate unmapped or fails QC or secondary alignment)
    else:
        query_flags = callable_flags + ["not duplicate"]
        flag = 1804  # as above plus not duplicate

    # Back compatible cache
    oldcache_file = _backcompatible_cache_file(query_flags, bed_file, target_name, data)
    if oldcache_file:
        with open(oldcache_file) as f:
            return int(f.read().strip())

    # New cache
    key = json.dumps({"flags": sorted(query_flags),
                      "region": os.path.basename(bed_file) if bed_file else "",
                      "sample": dd.get_sample_name(data)},
                     separators=(",", ":"), sort_keys=True)
    cache_file = get_cache_file(data)
    if utils.file_exists(cache_file):
        with open(cache_file) as in_handle:
            for cur_key, cur_val in (l.strip().split("\t") for l in in_handle):
                if cur_key == key:
                    return int(cur_val)

    # Calculate stats
    bam.index(bam_file, data["config"], check_timestamp=False)
    num_cores = dd.get_num_cores(data)
    count = 0
    with tx_tmpdir(data) as cur_tmpdir:
        # Covert to samtools regions (they are 1-based, BED is 0-based)
        regions = (["%s:%s-%s" % (r.chrom, r.start + 1, r.end) for r in pybedtools.BedTool(bed_file)]
                    if bed_file else [None])
        logger.debug("Count mapped reads with samtools view: %s" % (dd.get_sample_name(data)))
        for i, region_group in enumerate(tz.partition_all(10000, regions)):
            if len(region_group) == 1 and not region_group[0]:
                region_str = ""
            else:
                region_in = os.path.join(cur_tmpdir, "%s-regions-%s.bed" % (dd.get_sample_name(data), i))
                with open(region_in, "w") as out_handle:
                    out_handle.write(" ".join(region_group))
                region_str = " `cat %s`" % region_in
            count_out = os.path.join(cur_tmpdir, "%s-count-%s.txt" % (dd.get_sample_name(data), i))
            cmd = "samtools view -c -F {flag} -@ {num_cores} {bam_file}{region_str} > {count_out}"
            do.run(cmd.format(**locals()))
            with open(count_out) as in_handle:
                count += int(in_handle.read().strip())

    # Update cache
    with open(cache_file, "a") as out_handle:
        out_handle.write("%s\t%s\n" % (key, count))
    return count
