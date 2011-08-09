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
"""
import os
import sys
import subprocess
import optparse

import yaml

from bcbio import utils
from bcbio.distributed.messaging import create_celeryconfig

def main(config_file, queues=None):
    task_module = "bcbio.distributed.tasks"
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    with utils.curdir_tmpdir() as work_dir:
        dirs = {"work": work_dir, "config": os.path.dirname(config_file)}
        with create_celeryconfig(task_module, dirs, config,
                                 os.path.abspath(config_file)):
            run_celeryd(work_dir, queues)

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
    (options, args) = parser.parse_args()
    if len(args) != 1:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    main(args[0], options.queues)
