"""Transfer files from sequencer to remote analysis machine.
"""
import glob
import operator
import os
import subprocess

from bcbio import utils
from bcbio.log import logger

def copy_to_remote(dname, fastq_dir, sample_cfile, config):
    """Copy required files for processing to remote server using rsync.
    """
    with utils.chdir(dname):
        reports = reduce(operator.add,
                         [glob.glob("*.xml"),
                          glob.glob("Data/Intensities/BaseCalls/*.xml"),
                          glob.glob("Data/Intensities/BaseCalls/*.xsl"),
                          glob.glob("Data/Intensities/BaseCalls/*.htm"),
                          ["Data/Intensities/BaseCalls/Plots", "Data/reports",
                           "Data/Status.htm", "Data/Status_Files", "InterOp"]])
        run_info = reduce(operator.add,
                          [glob.glob("run_info.yaml"),
                           glob.glob("*.csv")])
        fastq = glob.glob(os.path.join(fastq_dir.replace(dname + "/", "", 1),
                                       "*.gz"))
        configs = [sample_cfile.replace(dname + "/", "", 1)]
    include_file = os.path.join(dname, "transfer_files.txt")
    with open(include_file, "w") as out_handle:
        out_handle.write("+ */\n")
        for fname in configs + fastq + run_info + reports:
            out_handle.write("+ %s\n" % fname)
        out_handle.write("- *\n")
    cmd = ["rsync", "-akmrtv", "--include-from=%s" % include_file, dname,
           "%s@%s:%s" % (utils.get_in(config, ("process", "username")),
                         utils.get_in(config, ("process", "host")),
                         utils.get_in(config, ("process", "dir")))]
    logger.info("Copying files to analysis machine")
    logger.info(" ".join(cmd))
    subprocess.check_call(cmd)
