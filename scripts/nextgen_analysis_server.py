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
import multiprocessing

import yaml
from mako.template import Template

from bcbio import utils

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    amqp_config = utils.read_galaxy_amqp_config(config["galaxy_config"])
    with utils.curdir_tmpdir() as work_dir:
        write_celeryconfig(work_dir, config, amqp_config)
        run_celeryd(work_dir)

def run_celeryd(work_dir):
    with utils.chdir(work_dir):
        subprocess.check_call("celeryd")

def write_celeryconfig(work_dir, config, amqp_config):
    out_file = os.path.join(work_dir, "celeryconfig.py")
    amqp_config["rabbitmq_vhost"] = config["analysis"]["rabbitmq_vhost"]
    amqp_config["cores"] = multiprocessing.cpu_count()
    amqp_config["task_import"] = "bcbio.distributed.tasks"
    with open(out_file, "w") as out_handle:
        out_handle.write(Template(celeryconfig_tmpl).render(**amqp_config))

celeryconfig_tmpl = """
CELERY_IMPORTS = ("${task_import}", )

BROKER_HOST = "${host}"
BROKER_PORT = "${port}"
BROKER_USER = "${userid}"
BROKER_PASSWORD = "${password}"
BROKER_VHOST = "${rabbitmq_vhost}"
CELERY_RESULT_BACKEND= "amqp"
CELERY_TASK_SERIALIZER = "json"
CELERYD_CONCURRENCY = ${cores}
"""

if __name__ == "__main__":
    main(*sys.argv[1:])
