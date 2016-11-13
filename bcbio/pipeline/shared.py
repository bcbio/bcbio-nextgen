"""Pipeline functionality shared amongst multiple analysis types.
"""
import os
from contextlib import contextmanager
import fileinput
import functools
import tempfile

import pybedtools
import pysam
import toolz as tz

from bcbio import bam, broad, utils
from bcbio.bam import ref
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.utils import file_exists, save_diskspace
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.provenance import do

# ## Split/Combine helpers

def combine_bam(in_files, out_file, config):
    """Parallel target to combine multiple BAM files.
    """
    runner = broad.runner_from_path("picard", config)
    runner.run_fn("picard_merge", in_files, out_file)
    for in_file in in_files:
        save_diskspace(in_file, "Merged into {0}".format(out_file), config)
    bam.index(out_file, config)
    return out_file

def get_noalt_contigs(data):
    """Retrieve contigs without alternatives as defined in bwa *.alts files.

    If no alt files present (when we're not aligning with bwa), work around
    with standard set of alts based on hg38 -- anything with HLA, _alt or
    _decoy in the name.
    """
    alts = set([])
    alt_files = [f for f in tz.get_in(["reference", "bwa", "indexes"], data, []) if f.endswith("alt")]
    if alt_files:
        for alt_file in alt_files:
            with open(alt_file) as in_handle:
                for line in in_handle:
                    if not line.startswith("@"):
                        alts.add(line.split()[0].strip())
    else:
        for contig in ref.file_contigs(dd.get_ref_file(data)):
            if ("_alt" in contig.name or "_decoy" in contig.name or
                  contig.name.startswith("HLA-") or ":" in contig.name):
                alts.add(contig.name)
    return [c for c in ref.file_contigs(dd.get_ref_file(data)) if c.name not in alts]

def write_nochr_reads(in_file, out_file, config):
    """Write a BAM file of reads that are not mapped on a reference chromosome.

    This is useful for maintaining non-mapped reads in parallel processes
    that split processing by chromosome.
    """
    if not file_exists(out_file):
        with file_transaction(config, out_file) as tx_out_file:
            samtools = config_utils.get_program("samtools", config)
            cmd = "{samtools} view -b -f 4 {in_file} > {tx_out_file}"
            do.run(cmd.format(**locals()), "Select unmapped reads")
    return out_file

def write_noanalysis_reads(in_file, region_file, out_file, config):
    """Write a BAM file of reads in the specified region file that are not analyzed.

    We want to get only reads not in analysis regions but also make use of
    the BAM index to perform well on large files. The tricky part is avoiding
    command line limits. There is a nice discussion on SeqAnswers:
    http://seqanswers.com/forums/showthread.php?t=29538
    sambamba supports intersection via an input BED file so avoids command line
    length issues.
    """
    if not file_exists(out_file):
        with file_transaction(config, out_file) as tx_out_file:
            bedtools = config_utils.get_program("bedtools", config)
            sambamba = config_utils.get_program("sambamba", config)
            cl = ("{sambamba} view -f bam -L {region_file} {in_file} | "
                  "{bedtools} intersect -abam - -b {region_file} -f 1.0 -nonamecheck"
                  "> {tx_out_file}")
            do.run(cl.format(**locals()), "Select unanalyzed reads")
    return out_file

def subset_bam_by_region(in_file, region, config, out_file_base=None):
    """Subset BAM files based on specified chromosome region.
    """
    if out_file_base is not None:
        base, ext = os.path.splitext(out_file_base)
    else:
        base, ext = os.path.splitext(in_file)
    out_file = "%s-subset%s%s" % (base, region, ext)
    if not file_exists(out_file):
        with pysam.Samfile(in_file, "rb") as in_bam:
            target_tid = in_bam.gettid(region)
            assert region is not None, \
                   "Did not find reference region %s in %s" % \
                   (region, in_file)
            with file_transaction(config, out_file) as tx_out_file:
                with pysam.Samfile(tx_out_file, "wb", template=in_bam) as out_bam:
                    for read in in_bam:
                        if read.tid == target_tid:
                            out_bam.write(read)
    return out_file

def subset_bed_by_chrom(in_file, chrom, data, out_dir=None):
    """Subset a BED file to only have items from the specified chromosome.
    """
    if out_dir is None:
        out_dir = os.path.dirname(in_file)
    base, ext = os.path.splitext(os.path.basename(in_file))
    out_file = os.path.join(out_dir, "%s-%s%s" % (base, chrom, ext))
    if not utils.file_uptodate(out_file, in_file):
        with file_transaction(data, out_file) as tx_out_file:
            _rewrite_bed_with_chrom(in_file, tx_out_file, chrom)
    return out_file

def _rewrite_bed_with_chrom(in_file, out_file, chrom):
    with open(in_file) as in_handle:
        with open(out_file, "w") as out_handle:
            for line in in_handle:
                if line.startswith("%s\t" % chrom):
                    out_handle.write(line)


def _subset_bed_by_region(in_file, out_file, region, do_merge=True):
    orig_bed = pybedtools.BedTool(in_file)
    region_bed = pybedtools.BedTool("\t".join(str(x) for x in region) + "\n", from_string=True)
    if do_merge:
        orig_bed.intersect(region_bed, nonamecheck=True).filter(lambda x: len(x) > 1).merge().saveas(out_file)
    else:
        orig_bed.intersect(region_bed, nonamecheck=True).filter(lambda x: len(x) > 1).saveas(out_file)

def get_lcr_bed(items):
    lcr_bed = utils.get_in(items[0], ("genome_resources", "variation", "lcr"))
    do_lcr = any([utils.get_in(data, ("config", "algorithm", "remove_lcr"), False)
                  for data in items])
    if do_lcr and lcr_bed and os.path.exists(lcr_bed):
        return lcr_bed

def remove_lcr_regions(orig_bed, items):
    """If configured and available, update a BED file to remove low complexity regions.
    """
    lcr_bed = get_lcr_bed(items)
    if lcr_bed:
        nolcr_bed = os.path.join("%s-nolcr.bed" % (utils.splitext_plus(orig_bed)[0]))
        with file_transaction(items[0], nolcr_bed) as tx_nolcr_bed:
            with bedtools_tmpdir(items[0]):
                pybedtools.BedTool(orig_bed).subtract(pybedtools.BedTool(lcr_bed), nonamecheck=True).\
                    saveas(tx_nolcr_bed)
        # If we have a non-empty file, convert to the LCR subtracted for downstream analysis
        if utils.file_exists(nolcr_bed):
            orig_bed = nolcr_bed
    return orig_bed

def remove_highdepth_regions(in_file, items):
    """Remove high depth regions from a BED file for analyzing a set of calls.

    Tries to avoid spurious errors and slow run times in collapsed repeat regions.

    Also adds ENCODE blacklist regions which capture additional collapsed repeats
    around centromeres.
    """
    from bcbio.variation import bedutils
    highdepth_beds = filter(lambda x: x is not None and utils.file_exists(x),
                            list(set([tz.get_in(["config", "algorithm", "highdepth_regions"], x) for x in items])))
    encode_bed = tz.get_in(["genome_resources", "variation", "encode_blacklist"], items[0])
    if encode_bed and os.path.exists(encode_bed):
        highdepth_beds.append(encode_bed)
    out_file = "%s-glimit%s" % utils.splitext_plus(in_file)
    if not utils.file_uptodate(out_file, in_file):
        with file_transaction(items[0], out_file) as tx_out_file:
            with bedtools_tmpdir(items[0]):
                all_file = "%s-all.bed" % utils.splitext_plus(tx_out_file)[0]
                if len(highdepth_beds) > 0:
                    with open(all_file, "w") as out_handle:
                        for highdepth_bed in highdepth_beds:
                            with utils.open_gzipsafe(highdepth_bed) as in_handle:
                                for line in in_handle:
                                    parts = line.split("\t")
                                    out_handle.write("\t".join(parts[:4]).rstrip() + "\n")
                if utils.file_exists(all_file):
                    to_remove = bedutils.sort_merge(all_file, items[0])
                    cmd = "bedtools subtract -nonamecheck -a {in_file} -b {to_remove} > {tx_out_file}"
                    do.run(cmd.format(**locals()), "Remove high depth regions")
                else:
                    utils.symlink_plus(in_file, out_file)
    return out_file

@contextmanager
def bedtools_tmpdir(data):
    with tx_tmpdir(data) as tmpdir:
        orig_tmpdir = tempfile.gettempdir()
        pybedtools.set_tempdir(tmpdir)
        yield
        if orig_tmpdir and os.path.exists(orig_tmpdir):
            pybedtools.set_tempdir(orig_tmpdir)
        else:
            tempfile.tempdir = None

def subtract_low_complexity(f):
    """Remove low complexity regions from callable regions if available.
    """
    @functools.wraps(f)
    def wrapper(variant_regions, region, out_file, items=None, do_merge=True, data=None):
        region_bed = f(variant_regions, region, out_file, items, do_merge, data)
        if region_bed and isinstance(region_bed, basestring) and os.path.exists(region_bed) and items:
            region_bed = remove_lcr_regions(region_bed, items)
        return region_bed
    return wrapper

@subtract_low_complexity
def subset_variant_regions(variant_regions, region, out_file, items=None, do_merge=True, data=None):
    """Return BED file subset by a specified chromosome region.

    variant_regions is a BED file, region is a chromosome name or tuple
    of (name, start, end) for a genomic region.
    """
    if region is None:
        return variant_regions
    elif variant_regions is None:
        return region
    elif not isinstance(region, (list, tuple)) and region.find(":") > 0:
        raise ValueError("Partial chromosome regions not supported")
    else:
        merge_text = "-unmerged" if not do_merge else ""
        subset_file = "{0}".format(utils.splitext_plus(out_file)[0])
        subset_file += "%s-regions.bed" % (merge_text)
        if not os.path.exists(subset_file):
            config = items[0] if items else data
            with file_transaction(config, subset_file) as tx_subset_file:
                if isinstance(region, (list, tuple)):
                    _subset_bed_by_region(variant_regions, tx_subset_file, region, do_merge = do_merge)
                else:
                    _rewrite_bed_with_chrom(variant_regions, tx_subset_file, region)
        if os.path.getsize(subset_file) == 0:
            return region
        else:
            return subset_file
