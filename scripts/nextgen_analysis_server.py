#!/usr/bin/env python
"""Start a nextgen analysis server that handles processing from a distributed task queue.

This reads configuration details and then starts a Celery (http://celeryproject.org>
server that will handle requests that are passed via an external message queue. This
allows distributed processing of analysis sections with less assumptions about the
system architecture.

Usage:
  nextgen_analysis_server.py <post_process.yaml>
   [--queues=list,of,queues: can specify specific queues to listen for jobs
                             on. No argument runs the default queue, which
                             handles processing alignments. 'toplevel' handles
                             managing the full work process.]
   [--tasks=task.module.import: Specify the module of tasks to make available.
                                Defaults to bcbio.distributed.tasks if not specified.]
   [--basedir=<dirname>: Base directory to work in. Defaults to current directory.]
"""
import os
import sys
import subprocess
import optparse

import yaml
from celery import signals

from bcbio import utils
from bcbio.distributed.messaging import create_celeryconfig
from bcbio.pipeline.config_loader import load_config
from bcbio.log import logger, setup_logging

def main(config_file, queues=None, task_module=None, base_dir=None):
    if base_dir is None:
        base_dir = os.getcwd()
    if task_module is None:
        task_module = "bcbio.distributed.tasks"
    config = load_config(config_file)
    if config.get("log_dir", None) is None:
        config["log_dir"] = os.path.join(base_dir, "log")
    signals.setup_logging.connect(celery_logger(config))
    setup_logging(config)
    logger.info("Starting distributed worker process: {0}".format(queues if queues else ""))
    with utils.chdir(base_dir):
        with utils.curdir_tmpdir() as work_dir:
            dirs = {"work": work_dir, "config": os.path.dirname(config_file)}
            with create_celeryconfig(task_module, dirs, config,
                                     os.path.abspath(config_file)):
                run_celeryd(work_dir, queues)

def celery_logger(config):
    def _worker(**kwds):
        setup_logging(config)
    return _worker

def run_celeryd(work_dir, queues):
    with utils.chdir(work_dir):
        cl = ["celeryd"]
        if queues:
            cl += ["-Q", queues]
        subprocess.check_call(cl)

if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option("-q", "--queues", dest="queues", action="store",
                      default=None)
    parser.add_option("-t", "--tasks", dest="task_module", action="store",
                      default=None)
    parser.add_option("-d", "--basedir", dest="basedir", action="store",
                      default=None)
    (options, args) = parser.parse_args()
    if len(args) != 1:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    main(args[0], options.queues, options.task_module, options.basedir)
