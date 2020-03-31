"""Handle conversions to/from CRAM reference based compression.
"""
import os
import subprocess

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd, config_utils
from bcbio.provenance import do


def compress(in_bam, data):
    """Compress a BAM file to CRAM, providing indexed CRAM file.

    Does 8 bin compression of quality score and read name removal
    using Staden io_lib if `cram` specified:

    https://github.com/jkbonfield/io_lib

    Otherwise does `cram-lossless` which only converts to CRAM.
    """
    out_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "archive"))
    out_file = os.path.join(out_dir, "%s.cram" % os.path.splitext(os.path.basename(in_bam))[0])
    cores = dd.get_num_cores(data)
    ref_file = dd.get_ref_file(data)
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            compress_type = dd.get_archive(data)
            scramble = config_utils.get_program("scramble", data["config"])
            cmd = [scramble, "-I", "bam", "-O", "cram", "-9", "-X", "archive", "-r", ref_file, "-V", "3.0", "-t", cores,
                   in_bam, tx_out_file]
            compressed = False
            if "cram" in compress_type:
                try:
                    cmd.extend(["-B", "-n"])
                    do.run(cmd, "Compress BAM to CRAM: quality score binning")
                    compressed = True
                except subprocess.CalledProcessError:
                    pass
            if not compressed:
                do.run(cmd, "Compress BAM to CRAM: lossless")
    index(out_file, data["config"])
    return out_file


def index(in_cram, config):
    """Ensure CRAM file has a .crai index file.
    """
    out_file = in_cram + ".crai"
    if not utils.file_uptodate(out_file, in_cram):
        with file_transaction(config, in_cram + ".crai") as tx_out_file:
            tx_in_file = os.path.splitext(tx_out_file)[0]
            utils.symlink_plus(in_cram, tx_in_file)
            cmd = ["cram_index", tx_in_file]
            do.run(cmd, "Index CRAM file")
    return out_file


def to_bam(in_file, out_file, data):
    """Convert CRAM file into BAM.
    """
    if not utils.file_uptodate(out_file, in_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = ["samtools", "view", "-O", "BAM", "-o", tx_out_file, in_file]
            do.run(cmd, "Convert CRAM to BAM")
    bam.index(out_file, data["config"])
    return out_file
