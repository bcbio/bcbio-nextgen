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
            if sambamba:
                cmd = "{sambamba} index -t {num_cores} {in_bam} {tx_index_file}"
            else:
                cmd = "{samtools} index {in_bam} {tx_index_file}"
            do.run(cmd.format(**locals()), "Index BAM file: %s" % os.path.basename(in_bam))
    return index_file if utils.file_exists(index_file) else alt_index_file
