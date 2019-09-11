"""Pipeline code to run alignments and prepare BAM files.

This works as part of the lane/flowcell process step of the pipeline.
"""
from collections import namedtuple
import glob
import os

import toolz as tz

from bcbio import bam, utils
from bcbio.bam import cram
from bcbio.ngsalign import (bbmap, bowtie, bwa, tophat, bowtie2, minimap2,
                            novoalign, snap, star, hisat2, bismark, bsmap)
from bcbio.pipeline import datadict as dd

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
    "bbmap": NgsTool(bbmap.align, None, None, bbmap.remap_index_fn),
    "bowtie": NgsTool(bowtie.align, None, bowtie.galaxy_location_file, None),
    "bowtie2": NgsTool(bowtie2.align, None,
                       bowtie2.galaxy_location_file, bowtie2.remap_index_fn),
    "bwa": NgsTool(bwa.align_pipe, bwa.align_bam, bwa.galaxy_location_file, None),
    "bismark": NgsTool(bismark.align, None, None, bismark.remap_index_fn),
    "bsmap": NgsTool(bsmap.align, None, None, None),
    "sentieon-bwa": NgsTool(bwa.align_pipe, bwa.align_bam, bwa.galaxy_location_file, None),
    "minimap2": NgsTool(minimap2.align, None, None, minimap2.remap_index_fn),
    "novoalign": NgsTool(novoalign.align_pipe, novoalign.align_bam,
                         novoalign.galaxy_location_file, novoalign.remap_index_fn),
    "tophat": NgsTool(tophat.align, None,
                      bowtie2.galaxy_location_file, bowtie2.remap_index_fn),
    "samtools": NgsTool(None, None, BASE_LOCATION_FILE, None),
    "snap": NgsTool(snap.align, None, snap.galaxy_location_file, snap.remap_index_fn),
    "star": NgsTool(star.align, None, None, star.remap_index_fn),
    "tophat2": NgsTool(tophat.align, None,
                       bowtie2.galaxy_location_file, bowtie2.remap_index_fn),
    "hisat2": NgsTool(hisat2.align, None, None, hisat2.remap_index_fn)}

metadata = {"support_bam": [k for k, v in TOOLS.items() if v.bam_align_fn is not None]}

def organize_noalign(data):
    """CWL target to skip alignment and organize input data.
    """
    data = utils.to_single_data(data[0])
    work_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "align", dd.get_sample_name(data)))
    work_bam = os.path.join(work_dir, "%s-input.bam" % dd.get_sample_name(data))
    if data.get("files"):
        if data["files"][0].endswith(".cram"):
            work_bam = cram.to_bam(data["files"][0], work_bam, data)
        else:
            assert data["files"][0].endswith(".bam"), data["files"][0]
            utils.copy_plus(data["files"][0], work_bam)
        bam.index(work_bam, data["config"])
    else:
        work_bam = None
    data["align_bam"] = work_bam
    return data

def align_to_sort_bam(fastq1, fastq2, aligner, data):
    """Align to the named genome build, returning a sorted BAM file.
    """
    names = data["rgnames"]
    align_dir_parts = [data["dirs"]["work"], "align", names["sample"]]
    if data.get("disambiguate"):
        align_dir_parts.append(data["disambiguate"]["genome_build"])
    aligner_index = get_aligner_index(aligner, data)
    align_dir = utils.safe_makedir(os.path.join(*align_dir_parts))
    ref_file = tz.get_in(("reference", "fasta", "base"), data)
    if fastq1.endswith(".bam"):
        data = _align_from_bam(fastq1, aligner, aligner_index, ref_file,
                               names, align_dir, data)
    else:
        data = _align_from_fastq(fastq1, fastq2, aligner, aligner_index, ref_file,
                                 names, align_dir, data)
    if data["work_bam"] and utils.file_exists(data["work_bam"]):
        if data.get("align_split") and dd.get_mark_duplicates(data):
            # If merging later with with bamsormadup need query sorted inputs
            # but CWL requires a bai file. Create a fake one to make it happy.
            bam.fake_index(data["work_bam"], data)
        else:
            bam.index(data["work_bam"], data["config"])
        for extra in ["-sr", "-disc"]:
            extra_bam = utils.append_stem(data['work_bam'], extra)
            if utils.file_exists(extra_bam):
                bam.index(extra_bam, data["config"])
    return data

def get_aligner_with_aliases(aligner, data):
    """Retrieve aligner index retriever, including aliases for shared.

    Handles tricky cases like gridss where we need bwa indices even with
    no aligner specified since they're used internally within GRIDSS.
    """
    aligner_aliases = {"sentieon-bwa": "bwa"}
    from bcbio import structural
    if not aligner and "gridss" in structural.get_svcallers(data):
        aligner = "bwa"
    return aligner_aliases.get(aligner) or aligner

def allow_noindices():
    return set(["minimap2"])

def get_aligner_index(aligner, data):
    """Handle multiple specifications of aligner indexes, returning value to pass to aligner.

    Original bcbio case -- a list of indices.
    CWL case: a single file with secondaryFiles staged in the same directory.
    """
    aligner_indexes = tz.get_in(("reference", get_aligner_with_aliases(aligner, data), "indexes"), data)
    # standard bcbio case
    if aligner_indexes and isinstance(aligner_indexes, (list, tuple)):
        aligner_index = os.path.commonprefix(aligner_indexes)
        if aligner_index.endswith("."):
            aligner_index = aligner_index[:-1]
        return aligner_index
    # single file -- check for standard naming or directory
    elif aligner_indexes and os.path.exists(aligner_indexes):
        aligner_dir = os.path.dirname(aligner_indexes)
        aligner_prefix = os.path.splitext(aligner_indexes)[0]
        if len(glob.glob("%s.*" % aligner_prefix)) > 0:
            return aligner_prefix
        else:
            return aligner_dir

    if aligner not in allow_noindices():
        raise ValueError("Did not find reference indices for aligner %s in genome: %s" %
                         (aligner, data["reference"]))

def _align_from_bam(fastq1, aligner, align_ref, sam_ref, names, align_dir, data):
    assert not data.get("align_split"), "Do not handle split alignments with BAM yet"
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
        assert out.get("work_bam"), (dd.get_sample_name(data), out.get("work_bam"))
        return out
    # handle output of raw SAM files that need to be converted to BAM
    else:
        work_bam = bam.sam_to_bam(out, config)
        data["work_bam"] = bam.sort(work_bam, config)
        return data
