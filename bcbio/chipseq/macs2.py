import os

import subprocess

from bcbio import utils
from bcbio.log import logger
from bcbio.provenance import do
import bcbio.pipeline.datadict as dd
from bcbio.pipeline import config_utils

HS = {"hg19": "2.7e9",
      "mm10": "1.87e9"}

def run(name, chip_bam, input_bam, genome_build, out_dir, config):
    """
    Run macs2 for chip and input samples avoiding
    errors due to samples.
    """
    out_file = os.path.join(out_dir, name + "_peaks.xls")
    if utils.file_exists(out_file):
        return out_file
    macs2 = config_utils.get_program("macs2", config)
    genome_size = HS[genome_build]
    with utils.chdir(out_dir):
        cmd = _macs2_cmd()
        try:
            do.run(cmd.format(**locals()), "macs2 for %s" % name)
        except subprocess.CalledProcessError:
            logger.debug("macs2 terminated with an error."
                         "please, check the message and report "
                         "error if related to bcbio.")
            # do.run("touch error" , "")
            return "error"
            pass
    return out_file

def _macs2_cmd():
    """Main command for macs2 tool."""
    cmd = ("{macs2} callpeak -t {chip_bam} -c {input_bam}"
            " -g {genome_size} -n {name} -B")
    return cmd
