"""Pipeline functionality shared amongst multiple analysis types.
"""
import os
import collections
from contextlib import closing

import pysam

from bcbio import broad
from bcbio.pipeline.alignment import get_genome_ref
from bcbio.utils import file_exists, safe_makedir, save_diskspace
from bcbio.distributed.transaction import file_transaction

# ## Split/Combine helpers

def combine_bam(in_files, out_file, config):
    """Parallel target to combine multiple BAM files.
    """
    runner = broad.runner_from_config(config)
    runner.run_fn("picard_merge", in_files, out_file)
    for in_file in in_files:
        save_diskspace(in_file, "Merged into {0}".format(out_file), config)
    runner.run_fn("picard_index", out_file)
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

def write_nochr_reads(in_file, out_file):
    """Write a BAM file of reads that are not on a reference chromosome.

    This is useful for maintaining non-mapped reads in parallel processes
    that split processing by chromosome.
    """
    if not file_exists(out_file):
        with closing(pysam.Samfile(in_file, "rb")) as in_bam:
            with file_transaction(out_file) as tx_out_file:
                with closing(pysam.Samfile(tx_out_file, "wb", template=in_bam)) as out_bam:
                    for read in in_bam:
                        if read.tid < 0:
                            out_bam.write(read)
    return out_file

def subset_bam_by_region(in_file, region, out_file_base = None):
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

def subset_variant_regions(variant_regions, region, out_file):
    """Return BED file subset by a specified chromosome region.

    variant_regions is a BED file, region is a chromosome name.
    """
    if region is None:
        return variant_regions
    elif variant_regions is None:
        return region
    elif region.find(":") > 0:
        raise ValueError("Partial chromosome regions not supported")
    else:
        # create an ordered subset file for processing
        subset_file = "{0}-regions.bed".format(os.path.splitext(out_file)[0])
        items = []
        with open(variant_regions) as in_handle:
            for line in in_handle:
                if line.startswith(region) and line.split("\t")[0] == region:
                    start = int(line.split("\t")[1])
                    items.append((start, line))
        if len(items) > 0:
            if not os.path.exists(subset_file):
                with open(subset_file, "w") as out_handle:
                    items.sort()
                    for _, line in items:
                        out_handle.write(line)
            return subset_file
        else:
            return region

# ## Retrieving file information from configuration variables

def configured_ref_file(name, config, sam_ref):
    """Full path to a reference file specified in the configuration.

    Resolves non-absolute paths relative to the base genome reference directory.
    """
    ref_file = config["algorithm"].get(name, None)
    if ref_file:
        if not os.path.isabs(ref_file):
            base_dir = os.path.dirname(os.path.dirname(sam_ref))
            ref_file = os.path.join(base_dir, ref_file)
    return ref_file

def configured_vrn_files(config, sam_ref):
    """Full path to all configured files for variation assessment.
    """
    names = ["dbsnp", "train_hapmap", "train_1000g_omni", "train_indels"]
    VrnFiles = collections.namedtuple("VrnFiles", names)
    return apply(VrnFiles, [configured_ref_file(n, config, sam_ref) for n in names])

def ref_genome_info(info, config, dirs):
    """Retrieve reference genome information from configuration variables.
    """
    genome_build = info.get("genome_build", None)
    (_, sam_ref) = get_genome_ref(genome_build, config["algorithm"]["aligner"],
                                  dirs["galaxy"])
    return genome_build, sam_ref
