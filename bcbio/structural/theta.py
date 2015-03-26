"""Provide estimates of sample purity and subclonal copy number using THetA.

Identifying cellularity and subclonal populations within somatic calling using
tumor normal pairs.

https://github.com/raphael-group/THetA
"""
import os
import sys
import subprocess

import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do

def run(cnv_info, somatic_info):
    """Run THetA analysis given output from CNV caller on a tumor/normal pair.
    """
    cmds = _get_cmds()
    if not cmds:
        logger.info("THetA scripts not found in current PATH. Skipping.")
        return cnv_info
    else:
        work_dir = _sv_workdir(somatic_info.tumor_data)
        exome_input = _create_exome_input(cmds, cnv_info, somatic_info, work_dir)
        _run_theta(cmds, exome_input, somatic_info.tumor_data, work_dir)
        return cnv_info

def _run_theta(cmds, exome_input, data, work_dir):
    """Run theta, calculating subpopulations and normal contamination.

    TODO: test on larger datasets and explore n=3 or higher values to evaluate
    multiple subclones.
    """
    out_dir = os.path.join(work_dir, "raw")
    result_file = os.path.join(out_dir, "%s.n2.results" % utils.splitext_plus(os.path.basename(exome_input))[0])
    if not utils.file_exists(result_file) and not utils.file_exists(result_file + ".skipped"):
        with file_transaction(data, out_dir) as tx_out_dir:
            utils.safe_makedir(tx_out_dir)
            cmd = cmds["run_theta"] + ["-n", "2", "-k", "4", "-m", ".90",
                                       exome_input, "--NUM_PROCESSES", dd.get_cores(data),
                                       "--FORCE", "-d", tx_out_dir]
            try:
                do.run(cmd, "Run THetA to calculate purity", log_error=False)
            except subprocess.CalledProcessError, msg:
                if ("Number of intervals must be greater than 1" in str(msg) or
                      "This sample isn't a good candidate for THetA analysis" in str(msg)):
                    with open(os.path.join(tx_out_dir,
                                           os.path.basename(result_file) + ".skipped"), "w") as out_handle:
                        out_handle.write("Expected TheTA failure, skipping")
                else:
                    raise
    return result_file

def _create_exome_input(cmds, cnv_info, somatic_info, work_dir):
    """Create exome inputs for THetA from existing CNV segmentation inputs.
    """
    out_prefix = os.path.join(work_dir, "%s_exome" % dd.get_sample_name(somatic_info.tumor_data))
    out_file = out_prefix + ".input"
    if not utils.file_exists(out_file):
        with file_transaction(somatic_info.tumor_data, out_file) as tx_out_file:
            tx_out_prefix = tx_out_file.replace(".input", "")
            target_bed = tz.get_in(["config", "algorithm", "variant_regions"], somatic_info.tumor_data)
            cmd = cmds["create_input"] + ["-t", somatic_info.tumor_bam, "-n", somatic_info.normal_bam,
                                          "-s", cnv_info["cns"], "--EXON_FILE", target_bed,
                                          "--FA", dd.get_ref_file(somatic_info.tumor_data),
                                          "--OUTPUT_PREFIX", tx_out_prefix]
            do.run(cmd, "Create exome inputs for THetA")
    return out_file

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "theta"))

def _get_cmds():
    """Retrieve required commands for running THetA with our local bcbio python.
    """
    cmds = {}
    for (name, cmd) in [("run_theta", "RunTHetA.py"), ("create_input", "createTHetAExomeInput.py")]:
        try:
            local_cmd = subprocess.check_output(["which", cmd]).strip()
        except subprocess.CalledProcessError:
            return None
        cmds[name] = [sys.executable, local_cmd]
    return cmds
