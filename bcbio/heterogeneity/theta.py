"""Provide estimates of sample purity and subclonal copy number using THetA.

Identifying cellularity and subclonal populations within somatic calling using
tumor normal pairs.

https://github.com/raphael-group/THetA
"""
import os
import sys
import subprocess

import numpy as np

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do

def run(vrn_info, cnvs_by_name, somatic_info):
    """Run THetA analysis given output from CNV caller on a tumor/normal pair.
    """
    cmd = _get_cmd("RunTHeTA.py")
    if not cmd:
        logger.info("THetA scripts not found in current PATH. Skipping.")
    else:
        from bcbio.structural import cnvkit
        work_dir = _sv_workdir(somatic_info.tumor_data)
        assert "cnvkit" in cnvs_by_name, "THetA requires CNVkit calls"
        cnv_info = cnvkit.export_theta(cnvs_by_name["cnvkit"], somatic_info.tumor_data)
        return _run_theta(cnv_info, somatic_info.tumor_data, work_dir)

def _run_theta(cnv_info, data, work_dir):
    """Run theta, calculating subpopulations and normal contamination.
    """
    out = {"caller": "theta"}
    max_normal = "0.9"
    opts = ["-m", max_normal]
    n2_result = _safe_run_theta(cnv_info["theta_input"], os.path.join(work_dir, "n2"), ".n2.results",
                                ["-n", "2"] + opts, data)
    if n2_result:
        out["estimate"] = n2_result
        n2_bounds = "%s.withBounds" % os.path.splitext(n2_result)[0]
        n3_result = _safe_run_theta(n2_bounds, os.path.join(work_dir, "n3"), ".n3.results",
                                    ["-n", "3", "--RESULTS", n2_result] + opts,
                                    data)
        if n3_result:
            best_result = _select_model(n2_bounds, n2_result, n3_result,
                                        os.path.join(work_dir, "n3"), data)
            out["estimate"] = best_result
            out["cnvs"] = _merge_theta_calls(n2_bounds, best_result, cnv_info["vrn_file"], data)
    return out

def _update_with_calls(result_file, cnv_file):
    """Update bounds with calls from CNVkit, inferred copy numbers and p-values from THetA.
    """
    results = {}
    with open(result_file) as in_handle:
        in_handle.readline()  # header
        _, _, cs, ps = in_handle.readline().strip().split()
        for i, (c, p) in enumerate(zip(cs.split(":"), ps.split(","))):
            results[i] = (c, p)
    cnvs = {}
    with open(cnv_file) as in_handle:
        for line in in_handle:
            chrom, start, end, _, count = line.rstrip().split()[:5]
            cnvs[(chrom, start, end)] = count
    def update(i, line):
        parts = line.rstrip().split("\t")
        chrom, start, end = parts[1:4]
        parts += cnvs.get((chrom, start, end), ".")
        parts += list(results[i])
        return "\t".join(parts) + "\n"
    return update

def _merge_theta_calls(bounds_file, result_file, cnv_file, data):
    """Create a final output file with merged CNVkit and THetA copy and population estimates.
    """
    out_file = "%s-merged.txt" % (result_file.replace(".BEST.results", ""))
    if not utils.file_uptodate(out_file, result_file):
        with file_transaction(data, out_file) as tx_out_file:
            updater = _update_with_calls(result_file, cnv_file)
            with open(bounds_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    i = 0
                    for line in in_handle:
                        if line.startswith("#"):
                            parts = line.rstrip().split("\t")
                            parts += ["cnv", "pop_cnvs", "pop_pvals"]
                            out_handle.write("\t".join(parts) + "\n")
                        else:
                            out_handle.write(updater(i, line))
                            i += 1
    return out_file

def _select_model(n2_bounds, n2_result, n3_result, out_dir, data):
    """Run final model selection from n=2 and n=3 options.
    """
    n2_out_file = n2_result.replace(".n2.results", ".BEST.results")
    n3_out_file = n3_result.replace(".n3.results", ".BEST.results")
    if not utils.file_exists(n2_out_file) and not utils.file_exists(n3_out_file):
        cmd = _get_cmd("ModelSelection.py") + [n2_bounds, n2_result, n3_result]
        do.run(cmd, "Select best THetA model")
    if utils.file_exists(n2_out_file):
        return n2_out_file
    else:
        assert utils.file_exists(n3_out_file)
        return n3_out_file

def _safe_run_theta(input_file, out_dir, output_ext, args, data):
    """Run THetA, catching and continuing on any errors.
    """
    out_file = os.path.join(out_dir, _split_theta_ext(input_file) + output_ext)
    skip_file = out_file + ".skipped"
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
    base = os.path.splitext(os.path.basename(fname))[0]
    if base.endswith((".n2", ".n3")):
        base = os.path.splitext(base)[0]
    return base

def _sv_workdir(data):
    return utils.safe_makedir(os.path.join(data["dirs"]["work"], "heterogeneity",
                                           dd.get_sample_name(data), "theta"))

def _get_cmd(cmd):
    """Retrieve required commands for running THetA with our local bcbio python.
    """
    check_cmd = "RunTHetA.py"
    try:
        local_cmd = subprocess.check_output(["which", check_cmd]).strip()
    except subprocess.CalledProcessError:
        return None
    return [sys.executable, "%s/%s" % (os.path.dirname(os.path.realpath(local_cmd)), cmd)]
