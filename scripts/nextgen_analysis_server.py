#!/usr/bin/env python
"""Start a nextgen analysis server that handles processing from a distributed task queue.

This reads configuration details and then starts a Celery (http://celeryproject.org>
server that will handle requests that are passed via an external message queue. This
allows distributed processing of analysis sections with less assumptions about the
system architecture.

Usage:
  nextgen_analysis_server.py <post_process.yaml>
"""
import os
import sys
import subprocess

import yaml

from bcbio import utils
from bcbio.distributed.messaging import create_celeryconfig

def main(config_file):
    task_module = "bcbio.distributed.tasks"
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    with utils.curdir_tmpdir() as work_dir:
        dirs = {"work": work_dir, "config": os.path.dirname(config_file)}
        with create_celeryconfig(task_module, dirs, config):
            run_celeryd(work_dir)

def run_celeryd(work_dir):
    with utils.chdir(work_dir):
        subprocess.check_call("celeryd")

if __name__ == "__main__":
    main(*sys.argv[1:])
