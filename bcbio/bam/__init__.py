"""Functionality to query and extract information from aligned BAM files.
"""
import contextlib
import os

import pysam

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.provenance import do

def is_paired(bam_file):
    """Determine if a BAM file has paired reads.
    """
    with contextlib.closing(pysam.Samfile(bam_file, "rb")) as in_pysam:
        for read in in_pysam:
            return read.is_paired

def index(in_bam, config):
    """Index a BAM file, skipping if index present.

    Centralizes BAM indexing providing ability to switch indexing approaches.
    """
    index_file = "%s.bai" % in_bam
    alt_index_file = "%s.bai" % os.path.splitext(in_bam)[0]
    if not utils.file_exists(index_file) and not utils.file_exists(alt_index_file):
        try:
            sambamba = config_utils.get_program("sambamba", config)
        except config_utils.CmdNotFound:
            sambamba = None
        samtools = config_utils.get_program("samtools", config)
        num_cores = config["algorithm"].get("num_cores", 1)
        with file_transaction(index_file) as tx_index_file:
            samtools_cmd = "{samtools} index {in_bam} {tx_index_file}"
            if sambamba:
                cmd = "{sambamba} index -t {num_cores} {in_bam} {tx_index_file}"
            else:
                cmd = samtools_cmd
            # sambamba has intermittent multicore failures. Allow
            # retries with single core
            try:
                do.run(cmd.format(**locals()), "Index BAM file: %s" % os.path.basename(in_bam),
                       log_error=False)
            except:
                do.run(samtools_cmd.format(**locals()),
                       "Index BAM file (single core): %s" % os.path.basename(in_bam))
    return index_file if utils.file_exists(index_file) else alt_index_file

def open_samfile(in_file):
    if is_bam(in_file):
        return pysam.Samfile(in_file, "rb")
    elif is_sam(in_file):
        return pysam.Samfile(in_file, "r")
    else:
        raise IOError("in_file must be either a BAM file or SAM file. Is the "
                      "extension .sam or .bam?")

def is_bam(in_file):
    _, ext = os.path.splitext(in_file)
    if ext == ".bam":
        return True
    else:
        return False


def is_sam(in_file):
    _, ext = os.path.splitext(in_file)
    if ext == ".sam":
        return True
    else:
        return False
