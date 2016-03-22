import os

import subprocess

from bcbio import utils
from bcbio.log import logger
from bcbio.provenance import do
import bcbio.pipeline.datadict as dd
from bcbio.pipeline import config_utils

HS = {"hg19": "2.7e9",
      "GRCh37": "2.7e9",
      "mm10": "1.87e9"}

def run(name, chip_bam, input_bam, genome_build, out_dir, config):
    """
    Run macs2 for chip and input samples avoiding
    errors due to samples.
    """
    # output file name need to have the caller name
    out_file = os.path.join(out_dir, name + "_peaks_macs2.xls")
    macs2_file = os.path.join(out_dir, name + "_peaks.xls")
    if utils.file_exists(out_file):
        return out_file
    macs2 = config_utils.get_program("macs2", config)
    options = " ".join(config_utils.get_resources("macs2", config).get("options", ""))
    if genome_build not in HS and options.find("-g") == -1:
        raise ValueError("This %s genome doesn't have a pre-set value."
                          "You can add specific values using resources "
                          "option for macs2 in the YAML file (-g genome_size)."
                          "Check Chip-seq configuration in "
                          "bcbio-nextgen documentation.")

    genome_size = "" if options.find("-g") > -1 else "-g %s" % HS[genome_build]
    with utils.chdir(out_dir):
        cmd = _macs2_cmd()
        try:
            do.run(cmd.format(**locals()), "macs2 for %s" % name)
            utils.move_safe(macs2_file, out_file)
        except subprocess.CalledProcessError:
            raise RuntimeWarning("macs2 terminated with an error.\n"
                                 "Please, check the message and report "
                                 "error if it is related to bcbio.\n"
                                 "You can add specific options for the sample "
                                 "setting resources as explained in docs: "
                                 "https://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#sample-specific-resources")
    return out_file

def _macs2_cmd():
    """Main command for macs2 tool."""
    cmd = ("{macs2} callpeak -t {chip_bam} -c {input_bam}"
            " {genome_size} -n {name} -B {options}")
    return cmd
