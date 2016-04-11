"""Handle conversions to/from CRAM reference based compression.

http://www.ebi.ac.uk/ena/about/cram_toolkit
"""
import os
import subprocess

from bcbio import utils
from bcbio.pipeline import config_utils
from bcbio.distributed.transaction import file_transaction

def compress(in_bam, ref_file, config):
    """Compress a BAM file to CRAM, binning quality scores. Indexes CRAM file.
    """
    out_file = "%s.cram" % os.path.splitext(in_bam)[0]
    resources = config_utils.get_resources("cram", config)
    jvm_opts = " ".join(resources.get("jvm_opts", ["-Xms1500m", "-Xmx3g"]))
    if not utils.file_exists(out_file):
        with file_transaction(config, out_file) as tx_out_file:
            cmd = ("cramtools {jvm_opts} cram "
                   "--input-bam-file {in_bam} "
                   "--capture-all-tags "
                   "--ignore-tags 'BD:BI' "
                   "--reference-fasta-file {ref_file} "
                   "--lossy-quality-score-spec '*8' "
                   "--output-cram-file {tx_out_file}")
            subprocess.check_call(cmd.format(**locals()), shell=True)
    index(out_file, config)
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
            subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file
