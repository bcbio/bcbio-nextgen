"""Pipeline code to run alignments and prepare BAM files.

This works as part of the lane/flowcell process step of the pipeline.
"""
from collections import namedtuple
import os

import toolz as tz

from bcbio import bam, utils
from bcbio.bam import cram
from bcbio.ngsalign import (bowtie, bwa, tophat, bowtie2,
                            novoalign, snap, star)

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
    "novoalign": NgsTool(novoalign.align_pipe, novoalign.align_bam,
                         novoalign.galaxy_location_file, novoalign.remap_index_fn),
    "tophat": NgsTool(tophat.align, None,
                      bowtie2.galaxy_location_file, bowtie2.remap_index_fn),
    "samtools": NgsTool(None, None, BASE_LOCATION_FILE, None),
    "snap": NgsTool(snap.align, snap.align_bam, snap.galaxy_location_file, snap.remap_index_fn),
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
    aligner_indexes = os.path.commonprefix(tz.get_in(("reference", aligner, "indexes"), data))
    if aligner_indexes.endswith("."):
        aligner_indexes = aligner_indexes[:-1]
    ref_file = tz.get_in(("reference", "fasta", "base"), data)
    if fastq1.endswith(".bam"):
        data = _align_from_bam(fastq1, aligner, aligner_indexes, ref_file,
                               names, align_dir, data)
    else:
        data = _align_from_fastq(fastq1, fastq2, aligner, aligner_indexes, ref_file,
                                 names, align_dir, data)
    if data["work_bam"] and utils.file_exists(data["work_bam"]):
        if not data.get("align_split"):
            bam.index(data["work_bam"], data["config"])
        for extra in ["-sr", "-disc"]:
            extra_bam = utils.append_stem(data['work_bam'], extra)
            if utils.file_exists(extra_bam):
                bam.index(extra_bam, data["config"])
    return data

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
    out = align_fn(fastq1, align_ref, names, align_dir, data)
    if isinstance(out, dict):
        assert "work_bam" in out
        return out
    else:
        data["work_bam"] = out
        return data

def _align_from_fastq(fastq1, fastq2, aligner, align_ref, sam_ref, names,
                      align_dir, data):
    """Align from fastq inputs, producing sorted BAM output.
    """
    config = data["config"]
    align_fn = TOOLS[aligner].align_fn
    out = align_fn(fastq1, fastq2, align_ref, names, align_dir, data)
    # handle align functions that update the main data dictionary in place
    if isinstance(out, dict):
        assert "work_bam" in out
        return out
    # handle output of raw SAM files that need to be converted to BAM
    else:
        work_bam = bam.sam_to_bam(out, config)
        data["work_bam"] = bam.sort(work_bam, config)
        return data
