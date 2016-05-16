"""Utilities for manipulating BED files.
"""
import os
import shutil
import sys
import subprocess

import toolz as tz

from bcbio import utils
from bcbio.bam import ref
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
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

def check_bed_contigs(in_file, data):
    """Ensure BED file contigs match the reference genome.
    """
    contigs = set([])
    with utils.open_gzipsafe(in_file) as in_handle:
        for line in in_handle:
            if not line.startswith(("#", "track", "browser")) and line.strip():
                contigs.add(line.split()[0])
    ref_contigs = set([x.name for x in ref.file_contigs(dd.get_ref_file(data))])
    if len(contigs - ref_contigs) / float(len(contigs)) > 0.25:
        raise ValueError("Contigs in BED file %s not in reference genome:\n %s\n"
                         % (in_file, list(contigs - ref_contigs)) +
                         "This is typically due to chr1 versus 1 differences in BED file and reference.")

def clean_file(in_file, data, prefix="", bedprep_dir=None):
    """Prepare a clean sorted input BED file without headers
    """
    if in_file:
        if not bedprep_dir:
            bedprep_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "bedprep"))
        out_file = os.path.join(bedprep_dir, "%s%s" % (prefix, os.path.basename(in_file))).replace(".gz", "")
        if not utils.file_uptodate(out_file, in_file):
            check_bed_contigs(in_file, data)
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
    if line.strip() and len(parts) > 2 and int(parts[2]) > int(parts[1]):
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
        if not utils.file_uptodate(out_file, in_file):
            with file_transaction(data, out_file) as tx_out_file:
                distance = "-d %s" % distance if distance else ""
                cmd = "{bedtools} merge {distance} -i {in_file} > {tx_out_file}"
                do.run(cmd.format(**locals()), "Prepare merged BED file", data)
        vcfutils.bgzip_and_index(out_file, data["config"], remove_orig=False)
        return out_file

def population_variant_regions(items):
    """Retrieve the variant region BED file from a population of items.

    If tumor/normal, return the tumor BED file. If a population, return
    the BED file covering the most bases.
    """
    import pybedtools
    if len(items) == 1:
        return dd.get_variant_regions(items[0])
    else:
        paired = vcfutils.get_paired(items)
        if paired:
            return dd.get_variant_regions(paired.tumor_data)
        else:
            vrs = []
            for data in items:
                vr_bed = dd.get_variant_regions(data)
                if vr_bed:
                    vrs.append((pybedtools.BedTool(vr_bed).total_coverage(), vr_bed))
            vrs.sort(reverse=True)
            if vrs:
                return vrs[0][1]

def clean_inputs(data):
    """Clean BED input files to avoid overlapping segments that cause downstream issues.

    Per-merges inputs to avoid needing to call multiple times during later parallel steps.
    """
    clean_vr = clean_file(utils.get_in(data, ("config", "algorithm", "variant_regions")), data)
    merged_vr = merge_overlaps(clean_vr, data)
    data["config"]["algorithm"]["variant_regions"] = clean_vr
    data["config"]["algorithm"]["variant_regions_merged"] = merged_vr
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

def intersect_two(f1, f2, work_dir, data):
    """Intersect two regions, handling cases where either file is not present.
    """
    f1_exists = f1 and utils.file_exists(f1)
    f2_exists = f2 and utils.file_exists(f2)
    if not f1_exists and not f2_exists:
        return None
    elif f1_exists and not f2_exists:
        return f1
    elif f2_exists and not f1_exists:
        return f2
    else:
        out_file = os.path.join(work_dir, "%s-merged.bed" % (utils.splitext_plus(os.path.basename(f1))[0]))
        if not utils.file_exists(out_file):
            with file_transaction(data, out_file) as tx_out_file:
                cmd = "bedtools intersect -a {f1} -b {f2} > {tx_out_file}"
                do.run(cmd.format(**locals()), "Intersect BED files", data)
        return out_file
