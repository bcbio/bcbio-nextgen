"""Functionality to query and extract information from aligned BAM files.
"""
import contextlib
import os

import pysam

from bcbio import broad, utils
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

def get_downsample_pct(runner, in_bam, target_counts):
    """Retrieve percentage of file to downsample to get to target counts.
    """
    total = sum(x.aligned for x in runner.run_fn("picard_idxstats", in_bam))
    with contextlib.closing(pysam.Samfile(in_bam, "rb")) as work_bam:
        n_rgs = max(1, len(work_bam.header["RG"]))
    rg_target = n_rgs * target_counts
    if total > rg_target:
        return float(rg_target) / float(total)

def downsample(in_bam, data, target_counts):
    """Downsample a BAM file to the specified number of target counts.
    """
    broad_runner = broad.runner_from_config(data["config"])
    ds_pct = get_downsample_pct(broad_runner, in_bam, target_counts)
    if ds_pct:
        out_file = "%s-downsample%s" % os.path.splitext(in_bam)
        if not utils.file_exists(out_file):
            with file_transaction(out_file) as tx_out_file:
                args = ["-T", "PrintReads",
                        "-R", data["sam_ref"],
                        "-I", in_bam,
                        "--downsample_to_fraction", "%.3f" % ds_pct,
                        "--out", tx_out_file]
                broad_runner.run_gatk(args)
        return out_file

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
