#!/usr/bin/env python
"""Run an automated analysis pipeline in a distributed cluster architecture.

Usage:
  run_distributed_job.py <config_file> <fc_dir> [<run_info_yaml>]
"""
import sys

import yaml

from bcbio.pipeline.run_info import get_run_info
from bcbio.distributed.manage import run_and_monitor

def main(config_file, fc_dir, run_info_yaml=None):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    assert config["algorithm"]["num_cores"] == "messaging", \
           "Use this script only with configured 'messaging' parallelization"
    workers_needed = _needed_workers(get_run_info(fc_dir, config, run_info_yaml)[-1])
    task_module = "bcbio.distributed.tasks"
    args = [config_file, fc_dir]
    if run_info_yaml:
        args.append(run_info_yaml)
    run_and_monitor(config, config_file, workers_needed, args, task_module)

def _needed_workers(run_info):
    """Determine workers needed to run multiplex flowcells in parallel.
    """
    names = []
    for lane in run_info["details"]:
        for multiplex in lane.get("multiplex", [{"barcode_id": ""}]):
            names.append((lane.get("name", ""), lane["description"], multiplex["barcode_id"]))
    return len(set(names))

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    main(*sys.argv[1:])
