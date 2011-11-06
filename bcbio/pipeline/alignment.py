"""Pipeline code to run alignments and prepare BAM files.

This works as part of the lane/flowcell process step of the pipeline.
"""
import os
from collections import namedtuple

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from bcbio import utils, broad
from bcbio.ngsalign import bowtie, bwa, tophat
from bcbio.distributed.transaction import file_transaction

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
    if fastq2 is None and aligner in ["bwa"]:
        fastq1 = _remove_read_number(fastq1, sam_file)
    return sam_to_sort_bam(sam_file, sam_ref, fastq1, fastq2, sample_name,
                           lane_name, config)

def _remove_read_number(in_file, sam_file):
    """Work around problem with MergeBamAlignment with BWA and single end reads.

    Need to remove read number ends from Fastq to match BWA stripping of numbers.

    http://sourceforge.net/mailarchive/forum.php?thread_name=87bosvbbqz.fsf%
    40fastmail.fm&forum_name=samtools-help
    http://sourceforge.net/mailarchive/forum.php?thread_name=4EB03C42.2060405%
    40broadinstitute.org&forum_name=samtools-help
    """
    out_file = os.path.join(os.path.dirname(sam_file),
                            "%s-safe%s" % os.path.splitext(os.path.basename(in_file)))
    if not os.path.exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            with open(in_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    for i, (name, seq, qual) in enumerate(FastqGeneralIterator(in_handle)):
                        if i == 0 and not name.endswith("/1"):
                            out_file = in_file
                            break
                        else:
                            name = name.rsplit("/", 1)[0]
                            out_handle.write("@%s\n%s\n+\n%s\n" % (name, seq, qual))
    return out_file

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
    # merge FASTQ files, only if barcoded samples in the work directory
    if (os.path.commonprefix([fastq1, sort_bam]) == os.path.dirname(sort_bam) and
          not config["algorithm"].get("upload_fastq", True)):
        utils.save_diskspace(fastq1, "Merged into output BAM %s" % out_bam, config)
        if fastq2:
            utils.save_diskspace(fastq2, "Merged into output BAM %s" % out_bam, config)
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

