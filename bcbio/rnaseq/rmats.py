import os

import subprocess
import commands
from bcbio import utils
from bcbio.log import logger
from bcbio.provenance import do
import bcbio.pipeline.datadict as dd
from bcbio.pipeline import config_utils

def _get_stranded_flag(config):
    strand_flag = {"unstranded": "fr-unstranded",
                   "firststrand": "fr-firststrand",
                   "secondstrand": "fr-secondstrand"}
    stranded = utils.get_in(config, ("algorithm", "strandedness"), "unstranded").lower()
    assert stranded in strand_flag, ("%s is not a valid strandedness value. "
                                     "Valid values are 'firststrand', "
                                     "'secondstrand' and 'unstranded" % (stranded))
    flag = strand_flag[stranded]
    return flag


def run(name, chip_bam, rep_bam, input_bam, gtf_file, out_dir, rlength, rpair, config):
    """
    Run rmats for muatant and control samples avoiding
    errors due to samples.
    """
    # output file name need to have the caller name
    MATS_output = os.path.join(out_dir, name + "_MATS_output")
    MATS_dir = os.path.join(out_dir, "MATS_output")
    rmats_file = os.path.join(out_dir, "summary.txt")
    out_file = os.path.join(out_dir, name + "_summary.txt")
    libType = _get_stranded_flag(config)
    if rep_bam != "":
            chip_bam = chip_bam + "," + rep_bam
    if utils.file_exists(out_file):
        return out_file
    rmats = config_utils.get_program("rmats", config)
    options = " ".join(config_utils.get_resources("rmats", config).get("options", ""))
    with utils.chdir(out_dir):
        cmd = _rmats_cmd()
        try:
            do.run(cmd.format(**locals()), "rmats for %s" % name)
            utils.move_safe(rmats_file, out_file)
            utils.move_safe(MATS_dir, MATS_output)
            repdir_dir = os.path.join(out_dir,"SAMPLE_1")
            utils.remove_safe(repdir_dir)
            repdir_dir = os.path.join(out_dir,"SAMPLE_2")
            utils.remove_safe(repdir_dir)
            print repdir_dir
        except subprocess.CalledProcessError:
            raise RuntimeWarning("rMATS terminated with an error.\n"
                                 "Please, check the message and report "
                                 "error if it is related to bcbio.\n"
                                 "You can add specific options for the sample "
                                 "setting resources as explained in docs: "
                                 "https://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#sample-specific-resources")
    return (out_file)

def _rmats_cmd():
    """Main command for rmats tool."""
    cmd = ("{rmats} -b1 {chip_bam} -b2 {input_bam}"
            " -gtf {gtf_file} -o {out_dir} -len {rlength} -analysis U -libType {libType} -t {rpair}")
    return cmd
