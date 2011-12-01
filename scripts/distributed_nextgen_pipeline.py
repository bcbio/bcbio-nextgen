#!/usr/bin/env python
"""Run an automated analysis pipeline in a distributed cluster architecture.

Usage:
  distributed_nextgen_pipeline.py <config_file> <fc_dir> [<run_info_yaml>]
     -n <number of processes to start>
"""
import sys
from optparse import OptionParser

import yaml

from bcbio.pipeline.run_info import get_run_info
from bcbio.distributed.manage import run_and_monitor
from bcbio.pipeline.config_loader import load_config

def main(config_file, fc_dir, run_info_yaml=None, num_workers=None):
    config = load_config(config_file)
    assert config["algorithm"]["num_cores"] == "messaging", \
           "Use this script only with configured 'messaging' parallelization"
    if num_workers is None:
        if config["distributed"].get("num_workers", "") == "all":
            cp = config["distributed"]["cluster_platform"]
            cluster = __import__("bcbio.distributed.{0}".format(cp), fromlist=[cp])
            num_workers = cluster.available_nodes(config["distributed"]["platform_args"]) - 1
        if num_workers is None:
            num_workers = _needed_workers(get_run_info(fc_dir, config, run_info_yaml)[-1])
    task_module = "bcbio.distributed.tasks"
    args = [config_file, fc_dir]
    if run_info_yaml:
        args.append(run_info_yaml)
    run_and_monitor(config, config_file, args, num_workers, task_module)

def _needed_workers(run_info):
    """Determine workers needed to run multiplex flowcells in parallel.
    """
    names = []
    for xs in run_info["details"]:
        for x in xs:
            names.append(x.get("name", (x["lane"], x["barcode_id"])))
    return len(set(names))

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-n", "--numworkers", dest="num_workers", default=None)
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    kwargs = {"num_workers": options.num_workers}
    main(*args, **kwargs)
