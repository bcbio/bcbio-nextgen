"""Perform archiving of output files for compressed storage.

Handles conversion to CRAM format.
"""

from bcbio import utils
from bcbio.bam import cram

def to_cram(data):
    """Convert BAM archive files into indexed CRAM.
    """
    cram_file = cram.compress(data["work_bam"], data)
    data["work_bam"] = cram_file
    return [[data]]

def compress(samples, run_parallel):
    """Perform compression of output files for long term storage.
    """
    to_cram = []
    finished = []
    for data in [x[0] for x in samples]:
        to_archive = set(utils.get_in(data, ("config", "algorithm", "archive"), []))
        if "cram" in to_archive:
            to_cram.append([data])
        else:
            finished.append([data])
    crammed = run_parallel("archive_to_cram", to_cram)
    return finished + crammed
