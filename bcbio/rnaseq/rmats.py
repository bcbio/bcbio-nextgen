import os

import subprocess

from bcbio import utils
from bcbio.log import logger
from bcbio.provenance import do
import bcbio.pipeline.datadict as dd
from bcbio.pipeline import config_utils

def run(name, chip_bam, input_bam, gtf_file, out_dir, read_len, config):
    """
    Run rmats for muatant and control samples avoiding
    errors due to samples.
    """
    # output file name need to have the caller name
    MATS_output = os.path.join(out_dir, "MATS_output")
    rmats_file = os.path.join(out_dir, "summary.txt")
    out_file = os.path.join(out_dir, name + "_summary.txt")
    if utils.file_exists(out_file):
        return out_file
    rmats = config_utils.get_program("rmats", config)
    options = " ".join(config_utils.get_resources("rmats", config).get("options", ""))
    with utils.chdir(out_dir):
        cmd = _rmats_cmd()
        try:
            do.run(cmd.format(**locals()), "rmats for %s" % name)
            utils.move_safe(rmats_file, out_file)
        except subprocess.CalledProcessError:
            raise RuntimeWarning("rMATS terminated with an error.\n"
                                 "Please, check the message and report "
                                 "error if it is related to bcbio.\n"
                                 "You can add specific options for the sample "
                                 "setting resources as explained in docs: "
                                 "https://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#sample-specific-resources")
    return (MATS_output, out_file)

def _rmats_cmd():
    """Main command for rmats tool."""
    cmd = ("{rmats} -b1 {chip_bam} -b2 {input_bam}"
            " -gtf {gtf_file} -o {out_dir} -len {read_len} -analysis U")
    return cmd
