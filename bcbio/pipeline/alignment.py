"""Pipeline code to run alignments and prepare BAM files.

This works as part of the lane/flowcell process step of the pipeline.
"""
from collections import namedtuple
import os
import sys
import glob

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from bcbio import bam, broad, utils
from bcbio.bam import cram
from bcbio.distributed.transaction import file_transaction
from bcbio.ngsalign import (bowtie, bwa, tophat, bowtie2, mosaik,
                            novoalign, star)
from bcbio.log import logger

# Define a next-generation sequencing tool to plugin:
# align_fn -- runs an aligner and generates SAM output
# galaxy_loc_file -- name of a Galaxy location file to retrieve
#  the genome index location
# bam_align_fn -- runs an aligner on a BAM file
# remap_index_fn -- Function that will take the location provided
#  from galaxy_loc_file and find the actual location of the index file.
#  This is useful for indexes that don't have an associated location file
#  but are stored in the same directory structure.
NgsTool = namedtuple("NgsTool", ["align_fn", "bam_align_fn",
                                 "galaxy_loc_file", "remap_index_fn"])


BASE_LOCATION_FILE = "sam_fa_indices.loc"

TOOLS = {
    "bowtie": NgsTool(bowtie.align, None, bowtie.galaxy_location_file, None),
    "bowtie2": NgsTool(bowtie2.align, None,
                       bowtie2.galaxy_location_file, bowtie2.remap_index_fn),
    "bwa": NgsTool(bwa.align_pipe, bwa.align_bam, bwa.galaxy_location_file, None),
    "mosaik": NgsTool(mosaik.align, None, mosaik.galaxy_location_file, None),
    "novoalign": NgsTool(novoalign.align_pipe, novoalign.align_bam,
                         novoalign.galaxy_location_file, novoalign.remap_index_fn),
    "tophat": NgsTool(tophat.align, None,
                      bowtie2.galaxy_location_file, bowtie2.remap_index_fn),
    "samtools": NgsTool(None, None, BASE_LOCATION_FILE, None),
    "star": NgsTool(star.align, None, None, star.remap_index_fn),
    "tophat2": NgsTool(tophat.align, None,
                       bowtie2.galaxy_location_file, bowtie2.remap_index_fn)}

metadata = {"support_bam": [k for k, v in TOOLS.iteritems() if v.bam_align_fn is not None]}

def align_to_sort_bam(fastq1, fastq2, aligner, data):
    """Align to the named genome build, returning a sorted BAM file.
    """
    names = data["rgnames"]
    align_dir_parts = [data["dirs"]["work"], "align", names["sample"]]
    if data.get("disambiguate"):
        align_dir_parts.append(data["disambiguate"]["genome_build"])
    align_dir = utils.safe_makedir(apply(os.path.join, align_dir_parts))
    if fastq1.endswith(".bam"):
        out_bam = _align_from_bam(fastq1, aligner, data["align_ref"], data["sam_ref"],
                                  names, align_dir, data)
        data["work_bam"] = out_bam
    else:
        data = _align_from_fastq(fastq1, fastq2, aligner, data["align_ref"], data["sam_ref"],
                                 names, align_dir, data)
    if data["work_bam"] and utils.file_exists(data["work_bam"]):
        bam.index(data["work_bam"], data["config"])
    return data

def _align_from_fastq_pipe(fastq1, fastq2, aligner, align_ref, sam_ref, names, align_dir, data):
    """Align longer reads using new piped strategies that avoid disk IO.
    """
    align_fn = TOOLS[aligner].pipe_align_fn
    if align_fn is None:
        raise NotImplementedError("Do not yet support piped alignment with %s" % aligner)
    return align_fn(fastq1, fastq2, align_ref, names, align_dir, data)

def _align_from_bam(fastq1, aligner, align_ref, sam_ref, names, align_dir, data):
    assert not data.get("align_split"), "Do not handle split alignments with BAM yet"
    config = data["config"]
    qual_bin_method = config["algorithm"].get("quality_bin")
    if (qual_bin_method == "prealignment" or
         (isinstance(qual_bin_method, list) and "prealignment" in qual_bin_method)):
        out_dir = utils.safe_makedir(os.path.join(align_dir, "qualbin"))
        fastq1 = cram.illumina_qual_bin(fastq1, sam_ref, out_dir, config)
    align_fn = TOOLS[aligner].bam_align_fn
    if align_fn is None:
        raise NotImplementedError("Do not yet support BAM alignment with %s" % aligner)
    return align_fn(fastq1, align_ref, names, align_dir, config)

def _align_from_fastq(fastq1, fastq2, aligner, align_ref, sam_ref, names,
                      align_dir, data):
    """Align from fastq inputs, producing sorted BAM output.
    """
    config = data["config"]
    align_fn = TOOLS[aligner].align_fn
    out = align_fn(fastq1, fastq2, align_ref, names, align_dir, data)
    if isinstance(out, basestring) and out.endswith(".sam"):
        if fastq2 is None and aligner in ["bwa", "bowtie2", "tophat2"]:
            fastq1 = _remove_read_number(fastq1, out)
        data["work_bam"] = sam_to_sort_bam(out, sam_ref, fastq1, fastq2, names, config)
        return data
    else:
        return out

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
    # file already exists and is zero means we already skipped the removal and
    # are just using the original file
    if os.path.exists(out_file) and os.path.getsize(out_file) == 0:
        return in_file
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

def sam_to_sort_bam(sam_file, ref_file, fastq1, fastq2, names, config):
    """Convert SAM file to merged and sorted BAM file.
    """
    picard = broad.runner_from_config(config)
    base_dir = os.path.dirname(sam_file)

    picard.run_fn("picard_index_ref", ref_file)
    out_fastq_bam = picard.run_fn("picard_fastq_to_bam", fastq1, fastq2, base_dir, names)
    out_bam = picard.run_fn("picard_sam_to_bam", sam_file, out_fastq_bam, ref_file,
                            fastq2 is not None)
    sort_bam = picard.run_fn("picard_sort", out_bam)

    utils.save_diskspace(sam_file, "SAM converted to BAM", config)
    utils.save_diskspace(out_fastq_bam, "Combined into output BAM %s" % out_bam, config)
    utils.save_diskspace(out_bam, "Sorted to %s" % sort_bam, config)
    # merge FASTQ files, only if barcoded samples in the work directory
    if (os.path.commonprefix([fastq1, sort_bam]) ==
             os.path.split(os.path.dirname(sort_bam))[0]
          and not config["algorithm"].get("upload_fastq", True)):
        utils.save_diskspace(fastq1, "Merged into output BAM %s" % out_bam, config)
        if fastq2:
            utils.save_diskspace(fastq2, "Merged into output BAM %s" % out_bam, config)
    return sort_bam
