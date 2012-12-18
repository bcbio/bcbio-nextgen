#!/usr/bin/env python
"""Run an automated analysis pipeline on nextgen sequencing information.

Handles runs in local or distributed mode based on the command line or
configured parameters.

The <config file> is a global YAML configuration file specifying details
about the system. An example configuration file is in 'config/post_process.yaml'.

<fc_dir> is an optional parameter specifying a directory of Illumina output
or fastq files to process. If configured to connect to a Galaxy LIMS system,
this can retrieve run information directly from Galaxy for processing.

<YAML run information> is on optional file specifies details about the
flowcell lanes, instead of retrieving it from Galaxy. An example
configuration file is located in 'config/run_info.yaml' This allows running
on files in arbitrary locations with no connection to Galaxy required.

Usage:
  bcbio_nextgen.py <config_file> [<fc_dir>] [<run_info_yaml>]
     -t type of parallelization to use:
          - local: Non-distributed, possibly multiple if n > 1 (default)
          - ipython: IPython distributed processing
          - messaging: RabbitMQ distributed messaging queue
     -n number of processes to start
     -p profile name for ipython parallelization
"""
import os
import sys

from bcbio.pipeline.run_info import get_run_info
from bcbio.distributed import manage as messaging
from bcbio.pipeline.config_utils import load_config
from bcbio.pipeline.main import run_main, parse_cl_args

def main(config_file, fc_dir=None, run_info_yaml=None, numcores=None,
         paralleltype=None, profile="default"):
    work_dir = os.getcwd()
    config = load_config(config_file)
    if config.get("log_dir", None) is None:
        config["log_dir"] = os.path.join(work_dir, "log")
    paralleltype, numcores = _get_cores_and_type(config, fc_dir, run_info_yaml,
                                                 numcores, paralleltype)
    parallel = {"type": paralleltype, "cores": numcores,
                "profile": profile,
                "module": "bcbio.distributed"}
    if parallel["type"] in ["local", "messaging-main"]:
        if numcores is None:
            config["algorithm"]["num_cores"] = numcores
        run_main(config, config_file, work_dir, parallel,
                 fc_dir, run_info_yaml)
    elif parallel["type"] == "messaging":
        parallel["task_module"] = "bcbio.distributed.tasks"
        args = [config_file, fc_dir]
        if run_info_yaml:
            args.append(run_info_yaml)
        messaging.run_and_monitor(config, config_file, args, parallel) 
    elif parallel["type"] == "ipython":
        run_main(config, config_file, work_dir, parallel,
                 fc_dir, run_info_yaml)
    else:
        raise ValueError("Unexpected type of parallel run: %s" % parallel["type"])

def _get_cores_and_type(config, fc_dir, run_info_yaml,
                        numcores=None, paralleltype=None):
    """Return core and parallelization approach from combo of config and commandline.

    Prefers passed commandline parameters over pre-configured, defaulting
    to a local run on a single core.

    The preferred approach is to pass in values explicitly on the commandline
    and this helps maintain back compatibility.
    """
    config_cores = config["algorithm"].get("num_cores", None)
    if config_cores:
        try:
            config_cores = int(config_cores)
            if numcores is None:
                numcores = config_cores
        except ValueError:
            if paralleltype is None:
                paralleltype = config_cores
    if paralleltype is None:
        paralleltype = "local"
        
    if numcores is None:
        if config["distributed"].get("num_workers", "") == "all":
            cp = config["distributed"]["cluster_platform"]
            cluster = __import__("bcbio.distributed.{0}".format(cp), fromlist=[cp])
            numcores = cluster.available_nodes(config["distributed"]["platform_args"]) - 1
        if numcores is None:
            if paralleltype == "local":
                numcores = 1
            else:
                numcores = _needed_workers(get_run_info(fc_dir, config, run_info_yaml)[-1])
    return paralleltype, int(numcores)

def _needed_workers(run_info):
    """Determine workers needed to run multiplex flowcells in parallel.
    """
    names = []
    for xs in run_info["details"]:
        for x in xs:
            names.append(x.get("name", (x["lane"], x["barcode_id"])))
    return len(set(names))

if __name__ == "__main__":
    config_file, kwargs = parse_cl_args(sys.argv[1:])
    main(config_file, **kwargs)
