"""Examine and query coverage in sequencing experiments.

Provides estimates of coverage intervals based on callable regions
"""
import collections
import os
import sys
import shutil
import glob

import toolz as tz
import yaml
import sqlite3
from pybedtools import BedTool
import pybedtools

from bcbio import utils, bed
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

def decorate_problem_regions(query_bed, problem_bed_dir):
    """
    decorate query_bed with percentage covered by BED files of regions specified
    in the problem_bed_dir
    """
    if utils.is_gzipped(query_bed):
        stem, _ = os.path.splitext(query_bed)
        stem, ext = os.path.splitext(stem)
    else:
        stem, ext = os.path.splitext(query_bed)
    out_file = stem + ".problem_annotated" + ext + ".gz"
    if utils.file_exists(out_file):
        return out_file
    bed_files = _find_bed_files(problem_bed_dir)
    bed_file_string = " ".join(bed_files)
    names = [os.path.splitext(os.path.basename(x))[0] for x in bed_files]
    names_string = " ".join(names)
    with utils.open_gzipsafe(query_bed) as in_handle:
        header = map(str, in_handle.next().strip().split())
    header = "\t".join(header + names)
    cmd = ("bedtools annotate -i {query_bed} -files {bed_file_string} "
           "-names {names_string} | sed -s 's/^#.*$/{header}/' | bgzip -c > {tx_out_file}")
    with file_transaction(out_file) as tx_out_file:
        message = "Annotate %s with problem regions." % query_bed
        do.run(cmd.format(**locals()), message)
    return out_file

def _find_bed_files(path):
    """
    recursively walk directories to find all of the BED files in the
    problem regions directory
    """
    bed_files = []
    for dirpath, subdirs, files in os.walk(path):
        for x in files:
            if x.endswith(".bed") or x.endswith(".bed.gz"):
                bed_files.append(os.path.join(dirpath, x))
    return bed_files
