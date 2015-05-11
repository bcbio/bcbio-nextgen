"""Provide estimates of sample purity and subclonal copy number using THetA.

Identifying cellularity and subclonal populations within somatic calling using
tumor normal pairs.

https://github.com/raphael-group/THetA
"""
import os
import sys
import subprocess

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.structural import cnvkit

def run(cnv_info, somatic_info):
    """Run THetA analysis given output from CNV caller on a tumor/normal pair.
    """
    cmd = _get_cmd("RunTHeTA.py")
    if not cmd:
        logger.info("THetA scripts not found in current PATH. Skipping.")
        return cnv_info
    else:
        work_dir = _sv_workdir(somatic_info.tumor_data)
        cnv_info = cnvkit.export_theta(cnv_info, somatic_info.tumor_data)
        cnv_info = _run_theta(cnv_info, somatic_info.tumor_data, work_dir)
        return cnv_info

def _run_theta(cnv_info, data, work_dir):
    """Run theta, calculating subpopulations and normal contamination.
    """
    max_cnv = "4"
    n2_result = _safe_run_theta(cnv_info["theta_input"], os.path.join(work_dir, "n2"), ".n2.results",
                                ["-n", "2", "-k", max_cnv], data)
    if n2_result:
        n2_bounds = "%s.withBounds" % os.path.splitext(n2_result)[0]
        n3_result = _safe_run_theta(n2_bounds, os.path.join(work_dir, "n3"), ".n3.results",
                                    ["-n", "3", "-k", max_cnv, "--RESULTS", n2_result], data)
        if n3_result:
            best_result = _select_model(n2_bounds, n2_result, n3_result,
                                        os.path.join(work_dir, "best"), data)
            cnv_info["theta"] = best_result
    return cnv_info

def _select_model(n2_bounds, n2_result, n3_result, out_dir, data):
    """Run final model selection from n=2 and n=3 options.
    """
    out_file = os.path.join(out_dir, _split_theta_ext(n2_bounds) + ".BEST.results")
    if not utils.file_exists(out_file):
        with file_transaction(data, out_dir) as tx_out_dir:
            utils.safe_makedir(tx_out_dir)
            with utils.chdir(tx_out_dir):
                cmd = _get_cmd("ModelSelection.py") + [n2_bounds, n2_result, n3_result]
                do.run(cmd, "Select best THetA model")
    return out_file

def _safe_run_theta(input_file, out_dir, output_ext, args, data):
    """Run THetA, catching and continuing on any errors.
    """
    out_file = _split_theta_ext(input_file) + output_ext
    skip_file = out_file = ".skipped"
    if utils.file_exists(skip_file):
        return None
    if not utils.file_exists(out_file):
        with file_transaction(data, out_dir) as tx_out_dir:
            utils.safe_makedir(tx_out_dir)
            cmd = _get_cmd("RunTHetA.py") + args + \
                  [input_file, "--NUM_PROCESSES", dd.get_cores(data),
                   "--FORCE", "-d", tx_out_dir]
            try:
                do.run(cmd, "Run THetA to calculate purity", log_error=False)
            except subprocess.CalledProcessError, msg:
                if ("Number of intervals must be greater than 1" in str(msg) or
                      "This sample isn't a good candidate for THetA analysis" in str(msg)):
                    with open(os.path.join(tx_out_dir, os.path.basename(skip_file)), "w") as out_handle:
                        out_handle.write("Expected TheTA failure, skipping")
                    return None
                else:
                    raise
    return out_file

def _split_theta_ext(fname):
    base = os.path.splitext(fname)[0]
    if base.endswith(".n2", ".n3"):
        base = os.path.splitext(base)[0]
    return base

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "structural",
                                           dd.get_sample_name(data), "theta"))

def _get_cmd(cmd):
    """Retrieve required commands for running THetA with our local bcbio python.
    """
    check_cmd = "RunTHetA.py"
    try:
        local_cmd = subprocess.check_output(["which", check_cmd]).strip()
    except subprocess.CalledProcessError:
        return None
    return [sys.executable, "%s/%s" % (os.path.dirname(local_cmd), cmd)]
