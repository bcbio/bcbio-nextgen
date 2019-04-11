"""Perform archiving of output files for compressed storage.

Handles conversion to CRAM format.
"""

from bcbio import utils
from bcbio.bam import cram
from bcbio.cwl import cwlutils
from bcbio.pipeline import datadict as dd

def to_cram(data):
    """Convert BAM archive files into indexed CRAM.
    """
    data = utils.to_single_data(data)
    cram_file = cram.compress(dd.get_work_bam(data) or dd.get_align_bam(data), data)
    out_key = "archive_bam" if cwlutils.is_cwl_run(data) else "work_bam"
    data[out_key] = cram_file
    return [[data]]

def compress(samples, run_parallel):
    """Perform compression of output files for long term storage.
    """
    to_cram = []
    finished = []
    for data in [x[0] for x in samples]:
        if "cram" in dd.get_archive(data) or "cram-lossless" in dd.get_archive(data):
            to_cram.append([data])
        else:
            finished.append([data])
    crammed = run_parallel("archive_to_cram", to_cram)
    return finished + crammed
