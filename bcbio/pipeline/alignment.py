"""Pipeline code to run alignments and prepare BAM files.

This works as part of the lane/flowcell process step of the pipeline.
"""
import os
from collections import namedtuple

from bcbio import utils, broad
from bcbio.ngsalign import bowtie, bwa, tophat

# Define a next-generation sequencing tool to plugin:
# align_fn -- runs an aligner and generates SAM output
# galaxy_loc_file -- name of a Galaxy location file to retrieve
#  the genome index location
# remap_index_fn -- Function that will take the location provided
#  from galaxy_loc_file and find the actual location of the index file.
#  This is useful for indexes that don't have an associated location file
#  but are stored in the same directory structure.
NgsTool = namedtuple("NgsTool", ["align_fn", "galaxy_loc_file",
                                 "remap_index_fn"])
_tools = {
    "bowtie": NgsTool(bowtie.align, bowtie.galaxy_location_file, None),
    "bwa": NgsTool(bwa.align, bwa.galaxy_location_file, None),
    "tophat": NgsTool(tophat.align, tophat.galaxy_location_file, None),
    "samtools": NgsTool(None, "sam_fa_indices.loc", None),
    }

def align_to_sort_bam(fastq1, fastq2, genome_build, aligner,
                      lane_name, sample_name, dirs, config):
    """Align to the named genome build, returning a sorted BAM file.
    """
    utils.safe_makedir(dirs["align"])
    align_ref, sam_ref = get_genome_ref(genome_build, aligner, dirs["galaxy"])
    align_fn = _tools[aligner].align_fn
    sam_file = align_fn(fastq1, fastq2, align_ref, lane_name, dirs["align"], config)
    return sam_to_sort_bam(sam_file, sam_ref, fastq1, fastq2, sample_name,
                           lane_name, config)

def sam_to_sort_bam(sam_file, ref_file, fastq1, fastq2, sample_name,
                    lane_name, config):
    """Convert SAM file to merged and sorted BAM file.
    """
    rg_name = lane_name.split("_")[0]
    picard = broad.runner_from_config(config)
    platform = config["algorithm"]["platform"]
    qual_format = config["algorithm"].get("quality_format", None)
    base_dir = os.path.dirname(sam_file)

    picard.run_fn("picard_index_ref", ref_file)
    out_fastq_bam = picard.run_fn("picard_fastq_to_bam", fastq1, fastq2,
                                  base_dir, platform, sample_name, rg_name, lane_name,
                                  qual_format)
    out_bam = picard.run_fn("picard_sam_to_bam", sam_file, out_fastq_bam, ref_file,
                            fastq2 is not None)
    sort_bam = picard.run_fn("picard_sort", out_bam)

    utils.save_diskspace(sam_file, "SAM converted to BAM", config)
    utils.save_diskspace(out_fastq_bam, "Combined into output BAM %s" % out_bam, config)
    utils.save_diskspace(out_bam, "Sorted to %s" % sort_bam, config)

    return sort_bam

def get_genome_ref(genome_build, aligner, galaxy_base):
    """Retrieve the reference genome file location from galaxy configuration.
    """
    if not aligner or not genome_build:
        return (None, None)
    ref_dir = os.path.join(galaxy_base, "tool-data")
    out_info = []
    for ref_get in [aligner, "samtools"]:
        ref_file = os.path.join(ref_dir, _tools[ref_get].galaxy_loc_file)
        cur_ref = None
        with open(ref_file) as in_handle:
            for line in in_handle:
                if line.strip() and not line.startswith("#"):
                    parts = line.strip().split()
                    if parts[0] == "index":
                        parts = parts[1:]
                    if parts[0] == genome_build:
                        cur_ref = parts[-1]
                        break
        if cur_ref is None:
            raise IndexError("Genome %s not found in %s" % (genome_build,
                ref_file))
        remap_fn = _tools[ref_get].remap_index_fn
        if remap_fn:
            cur_ref = remap_fn(cur_ref)
        out_info.append(utils.add_full_path(cur_ref, ref_dir))

    if len(out_info) != 2:
        raise ValueError("Did not find genome reference for %s %s" %
                (genome_build, aligner))
    else:
        return tuple(out_info)

