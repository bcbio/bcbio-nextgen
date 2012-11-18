#!/usr/bin/env python
"""Perform an automated analysis on a sequencing run using Galaxy information.

Given a directory of solexa output, this retrieves details about the sequencing
run from the Galaxy description, and uses this to perform an initial alignment
and analysis.

Usage:
    automated_initial_analysis.py <YAML config file>
                                  [<flow cell dir>]
                                  [<YAML run information>]

The optional <YAML run information> file specifies details about the
flowcell lanes, instead of retrieving it from Galaxy. An example
configuration file is located in 'config/run_info.yaml'

Workflow:
    - Retrieve details on a run.
    - Align fastq files to reference genome.
    - Perform secondary analyses like SNP calling.
    - Generate summary report.
"""
import os
import sys

from bcbio.pipeline.config_loader import load_config
from bcbio.pipeline.main import run_main, parse_cl_args

def main(config_file, fc_dir=None, run_info_yaml=None, num_cores=None):
    config = load_config(config_file)
    work_dir = os.getcwd()
    if config.get("log_dir", None) is None:
        config["log_dir"] = os.path.join(work_dir, "log")
    if num_cores:
        config["algorithm"]["num_cores"] = int(num_cores)
    run_main(config, config_file, work_dir, fc_dir, run_info_yaml)

if __name__ == "__main__":
    config_file, kwargs = parse_cl_args(sys.argv[1:])
    main(config_file, **kwargs)
