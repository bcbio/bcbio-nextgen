"""Examine and query coverage in sequencing experiments.

Provides estimates of coverage intervals based on callable regions and
stores coverage per regions in a database using Chanjo
(https://github.com/robinandeer/chanjo)
"""
import collections
import os
import sys

import toolz as tz
import yaml

from bcbio import utils
from bcbio.bam import ref
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import bedutils

def assign_interval(data):
    """Identify coverage based on percent of genome covered and relation to targets.

    Classifies coverage into 3 categories:
      - genome: Full genome coverage
      - regional: Regional coverage, like exome capture, with off-target reads
      - amplicon: Amplication based regional coverage without off-target reads
    """
    genome_cov_thresh = 0.40  # percent of genome covered for whole genome analysis
    offtarget_thresh = 0.10  # percent of offtarget reads required to be capture (not amplification) based
    import pybedtools
    if not dd.get_coverage_interval(data):
        vrs = dd.get_variant_regions(data)
        callable_file = dd.get_sample_callable(data)
        if vrs:
            seq_size = pybedtools.BedTool(vrs).total_coverage()
        else:
            seq_size = pybedtools.BedTool(callable_file).total_coverage()
        total_size = sum([c.size for c in ref.file_contigs(dd.get_ref_file(data), data["config"])])
        genome_cov_pct = seq_size / float(total_size)
        if genome_cov_pct > genome_cov_thresh:
            cov_interval = "genome"
            offtarget_pct = 0.0
        else:
            offtarget_stat_file = dd.get_offtarget_stats(data)
            if not offtarget_stat_file:
                offtarget_pct = 0.0
            else:
                with open(offtarget_stat_file) as in_handle:
                    stats = yaml.safe_load(in_handle)
                offtarget_pct = stats["offtarget"] / float(stats["mapped"])
            if offtarget_pct > offtarget_thresh:
                cov_interval = "regional"
            else:
                cov_interval = "amplicon"
        logger.info("Assigned coverage as '%s' with %.1f%% genome coverage and %.1f%% offtarget coverage"
                    % (cov_interval, genome_cov_pct * 100.0, offtarget_pct * 100.0))
        data["config"]["algorithm"]["coverage_interval"] = cov_interval
    return data

def summary(items):
    cutoff = 4  # coverage for completeness

    out_dir = utils.safe_makedir(os.path.join(items[0]["dirs"]["work"], "coverage"))
    clean_bed = bedutils.clean_file(tz.get_in(["config", "algorithm", "coverage"], items[0]),
                                    items[0])
    bed_file = _uniquify_bed_names(clean_bed, out_dir, items[0])
    batch = _get_group_batch(items)
    assert batch, "Did not find batch for samples: %s" % ",".join([dd.get_sample_name(x) for x in items])

    out_file = os.path.join(out_dir, "%s-coverage.db" % batch)
    if not utils.file_exists(out_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            chanjo = os.path.join(os.path.dirname(sys.executable), "chanjo")
            cmd = ("{chanjo} --db {tx_out_file} build {bed_file}")
            do.run(cmd.format(**locals()), "Prep chanjo database")
            for data in items:
                sample = dd.get_sample_name(data)
                bam_file = data["work_bam"]
                cmd = ("{chanjo} annotate -s {sample} -g {batch} -c {cutoff} {bam_file} {bed_file} | "
                       "{chanjo} --db {tx_out_file} import")
                do.run(cmd.format(**locals()), "Chanjo coverage", data)
    out = []
    for data in items:
        data["coverage"] = {"summary": out_file}
        out.append([data])
    return out

def _uniquify_bed_names(bed_file, out_dir, data):
    """Chanjo required unique names in the BED file to map to intervals.
    """
    out_file = os.path.join(out_dir, "%s-unames%s" % utils.splitext_plus(os.path.basename(bed_file)))
    if not utils.file_exists(out_file) or not utils.file_uptodate(out_file, bed_file):
        with file_transaction(data, out_file) as tx_out_file:
            with open(bed_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    namecounts = collections.defaultdict(int)
                    for i, line in enumerate(in_handle):
                        parts = line.rstrip("\r\n").split("\t")
                        if len(parts) >= 4:
                            name = parts[3]
                        else:
                            name = str(i)
                        namecount = namecounts.get(name, 0)
                        namecounts[name] += 1
                        if namecount > 0:
                            name = "%s-%s" % (name, namecount)
                        if len(parts) >= 4:
                            parts[3] = name
                        else:
                            assert len(parts) == 3
                            parts.append(name)
                        out_handle.write("\t".join(parts) + "\n")
    return out_file

def _get_group_batch(items):
    out = None
    all_batches = []
    for data in items:
        batches = tz.get_in(("metadata", "batch"), data, [dd.get_sample_name(data)])
        if not isinstance(batches, (list, tuple)):
            batches = [batches]
        batches = [b for b in batches if b]
        all_batches.extend(batches)
        if not out:
            out = set(batches)
        else:
            out = out.intersection(set(batches))
    if len(out) > 0:
        return out.pop()
    else:
        return all_batches[0]

def _handle_multi_batches(prepped, multi_batches):
    """Avoid carrying items present in multiple batches along in analysis.
    """
    out = []
    handled = set([])
    for data in (x[0] for x in prepped):
        name = dd.get_sample_name(data)
        if name in multi_batches:
            if name not in handled:
                out.append([data])
                handled.add(name)
            multi_batches.remove(name)
        elif name not in handled:
            out.append([data])
    assert len(multi_batches) == 0, "Did not find all multi_batch items: %s" % (list(multi_batches))
    return out

def summarize_samples(samples, run_parallel):
    """Back compatibility for existing pipelines. Should be replaced with summary when ready.
    """
    extras = []
    to_run = collections.defaultdict(list)
    multi_batches = set([])
    for data in [x[0] for x in samples]:
        if tz.get_in(["config", "algorithm", "coverage"], data):
            batches = tz.get_in(("metadata", "batch"), data, [dd.get_sample_name(data)])
            if not isinstance(batches, (tuple, list)):
                batches = [batches]
            else:
                multi_batches.add(dd.get_sample_name(data))
            for batch in batches:
                to_run[batch].append(utils.deepish_copy(data))
        else:
            extras.append([data])
    out = run_parallel("coverage_summary", [[xs] for xs in to_run.values()]) if len(to_run) > 0 else []
    out = _handle_multi_batches(out, multi_batches)
    assert len(out + extras) == len(samples), (len(out + extras), len(samples))
    return out + extras
