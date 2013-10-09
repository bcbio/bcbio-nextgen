#!/usr/bin/env python -E
"""Run an automated analysis pipeline for high throughput sequencing data.

Handles runs in local or distributed mode based on the command line or
configured parameters.

The <config file> is a global YAML configuration file specifying details
about the system. An example configuration file is in 'config/post_process.yaml'.
This is optional for automated installations.

<fc_dir> is an optional parameter specifying a directory of Illumina output
or fastq files to process. If configured to connect to a Galaxy LIMS system,
this can retrieve run information directly from Galaxy for processing.

<YAML run information> is an optional file specifies details about the
flowcell lanes, instead of retrieving it from Galaxy. An example
configuration file is located in 'config/bcbio_sample.yaml' This allows running
on files in arbitrary locations with no connection to Galaxy required.

Usage:
  bcbio_nextgen.py <config_file> [<fc_dir>] [<run_info_yaml>]
     -t type of parallelization to use:
          - local: Non-distributed, possibly multiple if n > 1 (default)
          - ipython: IPython distributed processing
     -n total number of processes to use
     -s scheduler for ipython parallelization (lsf, sge, slurm)
     -q queue to submit jobs for ipython parallelization
"""
import os
import sys

from bcbio import install, workflow
from bcbio.pipeline.main import run_main, parse_cl_args
from bcbio.server import main as server_main

def main(**kwargs):
    kwargs["work_dir"] = os.getcwd()
    run_main(**kwargs)

if __name__ == "__main__":
    kwargs = parse_cl_args(sys.argv[1:])
    if "upgrade" in kwargs and kwargs["upgrade"]:
        install.upgrade_bcbio(kwargs["args"])
    elif "server" in kwargs and kwargs["server"]:
        server_main.start(kwargs["args"])
    else:
        if kwargs.get("workflow"):
            setup_info = workflow.setup(kwargs["workflow"], kwargs["inputs"])
            if setup_info is None: # no automated run after setup
                sys.exit(0)
            workdir, new_kwargs = setup_info
            os.chdir(workdir)
            kwargs.update(new_kwargs)
        main(**kwargs)
