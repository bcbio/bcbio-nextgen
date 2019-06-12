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

def get_sort_cmd(tmp_dir=None):
    """Retrieve GNU coreutils sort command, using version-sort if available.

    Recent versions of sort have alpha-numeric sorting, which provides
    more natural sorting of chromosomes (chr1, chr2) instead of (chr1, chr10).
    This also fixes versions of sort, like 8.22 in CentOS 7.1, that have broken
    sorting without version sorting specified.

    https://github.com/bcbio/bcbio-nextgen/issues/624
    https://github.com/bcbio/bcbio-nextgen/issues/1017
    """
    has_versionsort = subprocess.check_output("sort --help | grep version-sort; exit 0", shell=True).strip()
    if has_versionsort:
        cmd = "sort -V"
    else:
        cmd = "sort"
    if tmp_dir and os.path.exists(tmp_dir) and os.path.isdir(tmp_dir):
        cmd += " -T %s" % tmp_dir
    return cmd

def check_bed_contigs(in_file, data):
    """Ensure BED file contigs match the reference genome.
    """
    if not dd.get_ref_file(data):
        return
    contigs = set([])
    with utils.open_gzipsafe(in_file) as in_handle:
        for line in in_handle:
            if not line.startswith(("#", "track", "browser", "@")) and line.strip():
                contigs.add(line.split()[0])
    ref_contigs = set([x.name for x in ref.file_contigs(dd.get_ref_file(data))])
    if contigs and len(contigs - ref_contigs) / float(len(contigs)) > 0.25:
        raise ValueError("Contigs in BED file %s not in reference genome:\n %s\n"
                         % (in_file, list(contigs - ref_contigs)) +
                         "This is typically due to chr1 versus 1 differences in BED file and reference.")

def has_regions(in_file):
    with utils.open_gzipsafe(in_file) as in_handle:
        for line in in_handle:
            if not line.startswith(("#", "track", "browser", "@")) and line.strip():
                return True
    return False

def check_bed_coords(in_file, data):
    """Ensure BED file coordinates match reference genome.

    Catches errors like using a hg38 BED file for an hg19 genome run.
    """
    if dd.get_ref_file(data):
        contig_sizes = {}
        for contig in ref.file_contigs(dd.get_ref_file(data)):
            contig_sizes[contig.name] = contig.size
        with utils.open_gzipsafe(in_file) as in_handle:
            for line in in_handle:
                if not line.startswith(("#", "track", "browser", "@")) and line.strip():
                    parts = line.split()
                    if len(parts) > 3:
                        try:
                            end = int(parts[2])
                        except ValueError:
                            continue
                        contig = parts[0]
                        check_size = contig_sizes.get(contig)
                        if check_size and end > check_size:
                            raise ValueError("Found BED coordinate off the end of the chromosome:\n%s%s\n"
                                             "Is the input BED from the right genome build?" %
                                             (line, in_file))

def clean_file(in_file, data, prefix="", bedprep_dir=None, simple=None):
    """Prepare a clean sorted input BED file without headers
    """
    # Remove non-ascii characters. Used in coverage analysis, to support JSON code in one column
    #   and be happy with sambamba:
    simple = "iconv -c -f utf-8 -t ascii | sed 's/ //g' |" if simple else ""
    if in_file:
        if not bedprep_dir:
            bedprep_dir = utils.safe_makedir(os.path.join(data["dirs"]["work"], "bedprep"))
        # Avoid running multiple times with same prefix
        if prefix and os.path.basename(in_file).startswith(prefix):
            return in_file
        out_file = os.path.join(bedprep_dir, "%s%s" % (prefix, os.path.basename(in_file)))
        out_file = out_file.replace(".interval_list", ".bed")
        if out_file.endswith(".gz"):
            out_file = out_file[:-3]
        if not utils.file_uptodate(out_file, in_file):
            check_bed_contigs(in_file, data)
            check_bed_coords(in_file, data)
            with file_transaction(data, out_file) as tx_out_file:
                bcbio_py = sys.executable
                cat_cmd = "zcat" if in_file.endswith(".gz") else "cat"
                sort_cmd = get_sort_cmd(os.path.dirname(tx_out_file))
                cmd = ("{cat_cmd} {in_file} | grep -v ^track | grep -v ^browser | grep -v ^@ | "
                       "grep -v ^# | {simple} "
                       "{bcbio_py} -c 'from bcbio.variation import bedutils; bedutils.remove_bad()' | "
                       "{sort_cmd} -k1,1 -k2,2n > {tx_out_file}")
                do.run(cmd.format(**locals()), "Prepare cleaned BED file", data)
        vcfutils.bgzip_and_index(out_file, data.get("config", {}), remove_orig=False)
        return out_file

def sort_merge(in_file, data, out_dir=None):
    """Sort and merge a BED file, collapsing gene names.
       Output is a 3 or 4 column file (the 4th column values go comma-separated).
    """
    out_file = "%s-sortmerge.bed" % os.path.splitext(in_file)[0]
    bedtools = config_utils.get_program("bedtools", data, default="bedtools")
    if out_dir:
        out_file = os.path.join(out_dir, os.path.basename(out_file))
    if not utils.file_uptodate(out_file, in_file):
        column_opt = ""
        with utils.open_gzipsafe(in_file) as in_handle:
            for line in in_handle:
                if not line.startswith(("#", "track", "browser", "@")):
                    parts = line.split()
                    if len(parts) >= 4:
                        column_opt = "-c 4 -o distinct"
        with file_transaction(data, out_file) as tx_out_file:
            cat_cmd = "zcat" if in_file.endswith(".gz") else "cat"
            sort_cmd = get_sort_cmd(os.path.dirname(tx_out_file))
            cmd = ("{cat_cmd} {in_file} | {sort_cmd} -k1,1 -k2,2n | "
                   "{bedtools} merge -i - {column_opt} > {tx_out_file}")
            do.run(cmd.format(**locals()), "Sort and merge BED file", data)
    return out_file

def remove_bad():
    """Remove non-increasing BED lines which will cause variant callers to choke.

    Also fixes space separated BED inputs.
    """
    for line in sys.stdin:
        parts = line.strip().split("\t")
        if len(parts) == 1 and len(line.strip().split()) > 1:
            parts = line.strip().split()
        if line.strip() and len(parts) > 2 and int(parts[2]) > int(parts[1]):
            sys.stdout.write("\t".join(parts) + "\n")

def merge_overlaps(in_file, data, distance=None, out_dir=None):
    """Merge bed file intervals to avoid overlapping regions. Output is always a 3 column file.

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

def population_variant_regions(items, merged=False):
    """Retrieve the variant region BED file from a population of items.

    If tumor/normal, return the tumor BED file. If a population, return
    the BED file covering the most bases.
    """
    def _get_variant_regions(data):
        out = dd.get_variant_regions(data) or dd.get_sample_callable(data)
        # Only need to merge for variant region inputs, not callable BED regions which don't overlap
        if merged and dd.get_variant_regions(data):
            merged_out = dd.get_variant_regions_merged(data)
            if merged_out:
                out = merged_out
            else:
                out = merge_overlaps(out, data)
        return out
    import pybedtools
    if len(items) == 1:
        return _get_variant_regions(items[0])
    else:
        paired = vcfutils.get_paired(items)
        if paired:
            return _get_variant_regions(paired.tumor_data)
        else:
            vrs = []
            for data in items:
                vr_bed = _get_variant_regions(data)
                if vr_bed:
                    vrs.append((pybedtools.BedTool(vr_bed).total_coverage(), vr_bed))
            vrs.sort(reverse=True)
            if vrs:
                return vrs[0][1]

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
    bedtools = config_utils.get_program("bedtools", data, default="bedtools")
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
                cmd = "{bedtools} intersect -a {f1} -b {f2} > {tx_out_file}"
                do.run(cmd.format(**locals()), "Intersect BED files", data)
        return out_file

def get_padded_bed_file(out_dir, bed_file, padding, data):
    bedtools = config_utils.get_program("bedtools", data, default="bedtools")
    out_file = os.path.join(out_dir, "%s-padded.bed" % (utils.splitext_plus(os.path.basename(bed_file))[0]))
    if utils.file_uptodate(out_file, bed_file):
        return out_file
    fai_file = ref.fasta_idx(dd.get_ref_file(data))
    with file_transaction(data, out_file) as tx_out_file:
        cmd = "{bedtools} slop -i {bed_file} -g {fai_file} -b {padding} | bedtools merge -i - > {tx_out_file}"
        do.run(cmd.format(**locals()), "Pad BED file", data)
    return out_file

def subset_to_genome(in_file, out_file, data):
    """Subset a BED file to only contain contigs present in the reference genome.
    """
    if not utils.file_uptodate(out_file, in_file):
        contigs = set([x.name for x in ref.file_contigs(dd.get_ref_file(data))])
        with utils.open_gzipsafe(in_file) as in_handle:
            with file_transaction(data, out_file) as tx_out_file:
                with open(tx_out_file, "w") as out_handle:
                    for line in in_handle:
                        parts = line.split()
                        if parts and parts[0] in contigs:
                            out_handle.write(line)
    return out_file
