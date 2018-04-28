"""Calculation of mapped reads by BAM counting, currently implemented with samtools.
"""
import contextlib
import json
import os
import time
import toolz as tz

from bcbio import bam, utils
from bcbio.bam import ref
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction

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
    so we have a simple file locking mechanism to avoid this.
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
    count_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "coverage",
                                                dd.get_sample_name(data), "counts"))
    if not bed_file:
        bed_file = os.path.join(count_dir, "fullgenome.bed")
        if not utils.file_exists(bed_file):
            with file_transaction(data, bed_file) as tx_out_file:
                with open(tx_out_file, "w") as out_handle:
                    for c in ref.file_contigs(dd.get_ref_file(data), data["config"]):
                        out_handle.write("%s\t%s\t%s\n" % (c.name, 0, c.size))
    count_file = os.path.join(count_dir,
                              "%s-%s-counts.txt" % (os.path.splitext(os.path.basename(bed_file))[0], flag))
    if not utils.file_exists(count_file):
        bam.index(bam_file, data["config"], check_timestamp=False)
        num_cores = dd.get_num_cores(data)
        with file_transaction(data, count_file) as tx_out_file:
            cmd = ("hts_nim_tools count-reads -t {num_cores} -F {flag} {bed_file} {bam_file} > {tx_out_file}")
            do.run(cmd.format(**locals()), "Count mapped reads: %s" % (dd.get_sample_name(data)))
    count = 0
    with open(count_file) as in_handle:
        for line in in_handle:
            count += int(line.rstrip().split()[-1])

    with _simple_lock(cache_file):
        with open(cache_file, "a") as out_handle:
            out_handle.write("%s\t%s\n" % (key, count))
    return count

@contextlib.contextmanager
def _simple_lock(f):
    """Simple file lock, times out after 20 second assuming lock is stale
    """
    lock_file = f + ".lock"
    timeout = 20
    curtime = 0
    interval = 2
    while os.path.exists(lock_file):
        time.sleep(interval)
        curtime += interval
        if curtime > timeout:
            os.remove(lock_file)
    with open(lock_file, "w") as out_handle:
        out_handle.write("locked")
    yield
    if os.path.exists(lock_file):
        os.remove(lock_file)
