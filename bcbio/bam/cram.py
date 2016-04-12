"""Handle conversions to/from CRAM reference based compression.

http://www.ebi.ac.uk/ena/about/cram_toolkit
"""
import os

from bcbio import utils
from bcbio.provenance import do
from bcbio.pipeline import datadict as dd
from bcbio.distributed.transaction import file_transaction

def compress(in_bam, data):
    """Compress a BAM file to CRAM, providing indexed CRAM file.

    XXX Does not 8-bin quality scores on output. We switched to samtools
    CRAM support due to problems with cramtools 3.0, and samtools does not yet
    support quality binning.
    """
    out_file = "%s.cram" % os.path.splitext(in_bam)[0]
    cores = dd.get_num_cores(data)
    ref_file = dd.get_ref_file(data)
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            cmd = ("samtools view -T {ref_file} -@ {cores} "
                   "-C -x BD -x BI "
                   "-o {tx_out_file} {in_bam}")
            do.run(cmd.format(**locals()), "Compress BAM to CRAM")
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
