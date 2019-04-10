"""Handle conversions to/from CRAM reference based compression.

http://www.ebi.ac.uk/ena/about/cram_toolkit
"""
import os
import subprocess

from bcbio import bam, utils
from bcbio.provenance import do
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.distributed.transaction import file_transaction

def compress(in_bam, data):
    """Compress a BAM file to CRAM, providing indexed CRAM file.

    Does 8 bin compression of quality score and read name removal
    using bamUtils squeeze if `cram` specified:

    http://genome.sph.umich.edu/wiki/BamUtil:_squeeze

    Otherwise does `cram-lossless` which only converts to CRAM.
    """
    out_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "archive"))
    out_file = os.path.join(out_dir, "%s.cram" % os.path.splitext(os.path.basename(in_bam))[0])
    cores = dd.get_num_cores(data)
    ref_file = dd.get_ref_file(data)
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            compress_type = dd.get_archive(data)
            samtools = config_utils.get_program("samtools", data["config"])
            try:
                bam_cmd = config_utils.get_program("bam", data["config"])
            except config_utils.CmdNotFound:
                bam_cmd = None
            to_cram = ("{samtools} view -T {ref_file} -@ {cores} "
                       "-C -x BD -x BI -o {tx_out_file}")
            compressed = False
            if "cram" in compress_type and bam_cmd:
                try:
                    cmd = ("{bam_cmd} squeeze --in {in_bam} --out -.ubam --keepDups "
                           "--binQualS=2,10,20,25,30,35,70 --binMid | " + to_cram)
                    do.run(cmd.format(**locals()), "Compress BAM to CRAM: quality score binning")
                    compressed = True
                # Retry failures avoiding using bam squeeze which can cause issues
                except subprocess.CalledProcessError:
                    pass
            if not compressed:
                cmd = (to_cram + " {in_bam}")
                do.run(cmd.format(**locals()), "Compress BAM to CRAM: lossless")
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
            cmd = "samtools index {tx_in_file}"
            do.run(cmd.format(**locals()), "Index CRAM file")
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
