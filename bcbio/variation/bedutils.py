"""Utilities for manipulating BED files.
"""
import os
import shutil
import sys
import subprocess

import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.variation import vcfutils

def get_sort_cmd():
    """Retrieve GNU coreutils sort command, using version-sort if available.

    Recent versions of sort have alpha-numeric sorting, which provides
    more natural sorting of chromosomes (chr1, chr2) instead of (chr1, chr10).
    This also fixes versions of sort, like 8.22 in CentOS 7.1, that have broken
    sorting without version sorting specified.

    https://github.com/chapmanb/bcbio-nextgen/issues/624
    https://github.com/chapmanb/bcbio-nextgen/issues/1017
    """
    has_versionsort = subprocess.check_output("sort --help | grep version-sort; exit 0", shell=True).strip()
    if has_versionsort:
        return "sort -V"
    else:
        return "sort"

def clean_file(in_file, data, prefix="", bedprep_dir=None):
    """Prepare a clean sorted input BED file without headers
    """
    if in_file:
        if not bedprep_dir:
            bedprep_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "bedprep"))
        out_file = os.path.join(bedprep_dir, "%s%s" % (prefix, os.path.basename(in_file))).replace(".gz", "")
        if not utils.file_uptodate(out_file, in_file):
            with file_transaction(data, out_file) as tx_out_file:
                py_cl = os.path.join(os.path.dirname(sys.executable), "py")
                cat_cmd = "zcat" if in_file.endswith(".gz") else "cat"
                sort_cmd = get_sort_cmd()
                cmd = ("{cat_cmd} {in_file} | grep -v ^track | grep -v ^browser | "
                       "grep -v ^# | "
                       "{py_cl} -x 'bcbio.variation.bedutils.remove_bad(x)' | "
                       "{sort_cmd} -k1,1 -k2,2n > {tx_out_file}")
                do.run(cmd.format(**locals()), "Prepare cleaned BED file", data)
        vcfutils.bgzip_and_index(out_file, data.get("config", {}), remove_orig=False)
        return out_file

def sort_merge(in_file, data):
    """Sort and merge a BED file, collapsing gene names.
    """
    out_file = "%s-sort.bed" % os.path.splitext(in_file)[0]
    if not utils.file_uptodate(out_file, in_file):
        with file_transaction(data, out_file) as tx_out_file:
            cat_cmd = "zcat" if in_file.endswith(".gz") else "cat"
            sort_cmd = get_sort_cmd()
            cmd = ("{cat_cmd} {in_file} | {sort_cmd} -k1,1 -k2,2n | "
                   "bedtools merge -i - -c 4 -o distinct > {tx_out_file}")
            do.run(cmd.format(**locals()), "Sort BED file", data)
    return out_file

def remove_bad(line):
    """Remove non-increasing BED lines which will cause variant callers to choke.
    """
    parts = line.strip().split("\t")
    if int(parts[2]) > int(parts[1]):
        return line
    else:
        return None

def merge_overlaps(in_file, data, distance=None, out_dir=None):
    """Merge bed file intervals to avoid overlapping regions.

    Overlapping regions (1:1-100, 1:90-100) cause issues with callers like FreeBayes
    that don't collapse BEDs prior to using them.
    """
    config = data["config"]
    if in_file:
        bedtools = config_utils.get_program("bedtools", config,
                                            default="bedtools")
        work_dir = tz.get_in(["dirs", "work"], data)
        if out_dir:
            bedprep_dir = out_dir
        elif work_dir:
            bedprep_dir = utils.safe_makedir(os.path.join(work_dir, "bedprep"))
        else:
            bedprep_dir = os.path.dirname(in_file)
        out_file = os.path.join(bedprep_dir, "%s-merged.bed" % (utils.splitext_plus(os.path.basename(in_file))[0]))
        if not utils.file_exists(out_file):
            with file_transaction(data, out_file) as tx_out_file:
                distance = "-d %s" % distance if distance else ""
                cmd = "{bedtools} merge {distance} -i {in_file} > {tx_out_file}"
                do.run(cmd.format(**locals()), "Prepare merged BED file", data)
        vcfutils.bgzip_and_index(out_file, data["config"], remove_orig=False)
        return out_file

def clean_inputs(data):
    """Clean BED input files to avoid overlapping segments that cause downstream issues.

    Per-merges inputs to avoid needing to call multiple times during later parallel steps.
    """
    clean_vr = clean_file(utils.get_in(data, ("config", "algorithm", "variant_regions")), data)
    merge_overlaps(clean_vr, data)
    data["config"]["algorithm"]["variant_regions"] = clean_vr
    return data

def combine(in_files, out_file, config):
    """Combine multiple BED files into a single output.
    """
    if not utils.file_exists(out_file):
        with file_transaction(config, out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                for in_file in in_files:
                    with open(in_file) as in_handle:
                        shutil.copyfileobj(in_handle, out_handle)
    return out_file
