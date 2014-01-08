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
NgsTool = namedtuple("NgsTool", ["align_fn", "pipe_align_fn", "bam_align_fn",
                                 "galaxy_loc_file", "remap_index_fn", "can_pipe",
                                 "index_fn"])

BASE_LOCATION_FILE = "sam_fa_indices.loc"

TOOLS = {
    "bowtie": NgsTool(bowtie.align, None, None,
                      bowtie.galaxy_location_file, None, None,
                      None),
    "bowtie2": NgsTool(bowtie2.align, None, None,
                       bowtie2.galaxy_location_file, bowtie2.remap_index_fn, None,
                       None),
    "bwa": NgsTool(bwa.align, bwa.align_pipe, bwa.align_bam,
                   bwa.galaxy_location_file, None, bwa.can_pipe,
                   None),
    "mosaik": NgsTool(mosaik.align, None, None,
                      mosaik.galaxy_location_file, None, None,
                      None),
    "novoalign": NgsTool(novoalign.align, novoalign.align_pipe, novoalign.align_bam,
                         novoalign.galaxy_location_file, novoalign.remap_index_fn, novoalign.can_pipe,
                         None),
    "tophat": NgsTool(tophat.align, None, None,
                      bowtie2.galaxy_location_file, bowtie2.remap_index_fn, None,
                      None),
    "samtools": NgsTool(None, None, None, BASE_LOCATION_FILE,
                        None, None, None),
    "star": NgsTool(star.align, None, None,
                    None, star.remap_index_fn, None,
                    star.index),
    "tophat2": NgsTool(tophat.align, None, None,
                       bowtie2.galaxy_location_file, bowtie2.remap_index_fn, None,
                       None)}

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
    elif _can_pipe(aligner, fastq1, data):
        data = _align_from_fastq_pipe(fastq1, fastq2, aligner, data["align_ref"], data["sam_ref"],
                                      names, align_dir, data)
    else:
        out_bam = _align_from_fastq(fastq1, fastq2, aligner, data["align_ref"], data["sam_ref"],
                                    names, align_dir, data)
        data["work_bam"] = out_bam
    if data["work_bam"] and utils.file_exists(data["work_bam"]):
        bam.index(data["work_bam"], data["config"])
    return data

def _can_pipe(aligner, fastq_file, data):
    """Check if current aligner support piping for a particular input fastq file.
    """
    if TOOLS[aligner].can_pipe and TOOLS[aligner].pipe_align_fn:
        return TOOLS[aligner].can_pipe(fastq_file, data)
    return False

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
    assert not data.get("align_split"), "Do not handle split alignments with non-piped fastq yet"
    config = data["config"]
    align_fn = TOOLS[aligner].align_fn
    sam_file = align_fn(fastq1, fastq2, align_ref, names["lane"], align_dir, data,
                        names=names)
    if fastq2 is None and aligner in ["bwa", "bowtie2", "tophat2"]:
        fastq1 = _remove_read_number(fastq1, sam_file)
    return sam_to_sort_bam(sam_file, sam_ref, fastq1, fastq2, names, config)

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


def make_missing_index(item):
    aligner = item["algorithm"].get("aligner", None)
    index_loc = item.get("align_ref", None)
    indexer_fn = TOOLS[aligner].index_fn if aligner in TOOLS else None
    if not index_loc or os.path.exists(index_loc):
        return item
    if not indexer_fn:
        logger.error("Index for %s is missing and the code to generate "
                     "it is not in place. Please open an issue here: "
                     "https://github.com/chapmanb/bcbio-nextgen/issues?state=open"
                     % aligner)
        sys.exit(1)
    else:
        logger.info("Index for %s is missing so it is being generated. This may "
                    "take a couple of hours depending on the index but it will "
                    "only happen the first time bcbio-nextgen is run." % aligner)
        indexer_fn(item)
    return item

def make_missing_indices(lane_items, run_parallel):
    items = [x[0] for x in lane_items]
    for item in items:
        index_loc = item.get("align_ref", None)
        if not _index_exists(index_loc):
            run_parallel("make_missing_index", [[item]])

def _index_exists(index_loc):
    files_match = glob.glob(index_loc + ".*")
    if os.path.exists(index_loc) or files_match:
        return True
    return False
