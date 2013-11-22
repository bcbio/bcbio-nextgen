"""Pipeline functionality shared amongst multiple analysis types.
"""
import os
from contextlib import closing

import pybedtools
import pysam

from bcbio import bam, broad
from bcbio.pipeline import config_utils
from bcbio.utils import file_exists, safe_makedir, save_diskspace
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do

# ## Split/Combine helpers

def combine_bam(in_files, out_file, config):
    """Parallel target to combine multiple BAM files.
    """
    runner = broad.runner_from_config(config)
    runner.run_fn("picard_merge", in_files, out_file)
    for in_file in in_files:
        save_diskspace(in_file, "Merged into {0}".format(out_file), config)
    bam.index(out_file, config)
    return out_file

def process_bam_by_chromosome(output_ext, file_key, default_targets=None, dir_ext_fn=None):
    """Provide targets to process a BAM file by individual chromosome regions.

    output_ext: extension to supply to output files
    file_key: the key of the BAM file in the input data map
    default_targets: a list of extra chromosome targets to process, beyond those specified
                     in the BAM file. Useful for retrieval of non-mapped reads.
    dir_ext_fn: A function to retrieve a directory naming extension from input data map.
    """
    if default_targets is None:
        default_targets = []
    def _do_work(data):
        bam_file = data[file_key]
        out_dir = os.path.dirname(bam_file)
        if dir_ext_fn:
            out_dir = os.path.join(out_dir, dir_ext_fn(data))

        out_file = os.path.join(out_dir, "{base}{ext}".format(
                base=os.path.splitext(os.path.basename(bam_file))[0],
                ext=output_ext))
        part_info = []
        if not file_exists(out_file):
            work_dir = safe_makedir(
                "{base}-split".format(base=os.path.splitext(out_file)[0]))
            with closing(pysam.Samfile(bam_file, "rb")) as work_bam:
                for chr_ref in list(work_bam.references) + default_targets:
                    chr_out = os.path.join(work_dir,
                                           "{base}-{ref}{ext}".format(
                                               base=os.path.splitext(os.path.basename(bam_file))[0],
                                               ref=chr_ref, ext=output_ext))
                    part_info.append((chr_ref, chr_out))
        return out_file, part_info
    return _do_work

def write_nochr_reads(in_file, out_file, config):
    """Write a BAM file of reads that are not mapped on a reference chromosome.

    This is useful for maintaining non-mapped reads in parallel processes
    that split processing by chromosome.
    """
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
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
    """
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            bedtools = config_utils.get_program("bedtools", config)
            samtools = config_utils.get_program("samtools", config)
            region_str = " ".join("%s:%s-%s" % tuple(r) for r in pybedtools.BedTool(region_file))
            cl = ("{samtools} view -b {in_file} {region_str} | "
                  "{bedtools} intersect -abam - -b {region_file} -f 1.0 "
                  "> {tx_out_file}")
            do.run(cl.format(**locals()), "Select unanalyzed reads")
    return out_file

def subset_bam_by_region(in_file, region, out_file_base=None):
    """Subset BAM files based on specified chromosome region.
    """
    if out_file_base is not None:
        base, ext = os.path.splitext(out_file_base)
    else:
        base, ext = os.path.splitext(in_file)
    out_file = "%s-subset%s%s" % (base, region, ext)
    if not file_exists(out_file):
        with closing(pysam.Samfile(in_file, "rb")) as in_bam:
            target_tid = in_bam.gettid(region)
            assert region is not None, \
                   "Did not find reference region %s in %s" % \
                   (region, in_file)
            with file_transaction(out_file) as tx_out_file:
                with closing(pysam.Samfile(tx_out_file, "wb", template=in_bam)) as out_bam:
                    for read in in_bam:
                        if read.tid == target_tid:
                            out_bam.write(read)
    return out_file

def _rewrite_bed_with_chrom(in_file, out_file, chrom):
    with open(in_file) as in_handle:
        with open(out_file, "w") as out_handle:
            for line in in_handle:
                if line.startswith("%s\t" % chrom):
                    out_handle.write(line)

def _subset_bed_by_region(in_file, out_file, region):
    orig_bed = pybedtools.BedTool(in_file)
    region_bed = pybedtools.BedTool("\t".join(str(x) for x in region) + "\n", from_string=True)
    orig_bed.intersect(region_bed).saveas(out_file)

def subset_variant_regions(variant_regions, region, out_file):
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
        subset_file = "{0}-regions.bed".format(os.path.splitext(out_file)[0])
        if not os.path.exists(subset_file):
            with file_transaction(subset_file) as tx_subset_file:
                if isinstance(region, (list, tuple)):
                    c, s, e = region
                    safe_region = [c, s, e - 2]
                    _subset_bed_by_region(variant_regions, tx_subset_file, safe_region)
                else:
                    _rewrite_bed_with_chrom(variant_regions, tx_subset_file, region)
        if os.path.getsize(subset_file) == 0:
            return region
        else:
            return subset_file
