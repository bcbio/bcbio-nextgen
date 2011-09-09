"""Run distributed tasks using the Celery distributed task queue.

http://celeryproject.org/
"""

import os
import sys
import time
import contextlib
import multiprocessing

from mako.template import Template

from bcbio import utils

def runner(dirs, config, config_file, wait=True):
    """Run a set of tasks using Celery, waiting for results or asynchronously.

    Initialize with the configuration and directory information,
    used to prepare a Celery configuration file and imports. It
    returns a function which acts like standard map; provide the function
    name instead of the function itself when calling.

    After name lookup, Celery runs the function in parallel; Celery servers
    can be remote or local but must have access to a shared filesystem. The
    function polls if wait is True, returning when all results are available.
    """
    task_module = "bcbio.distributed.tasks"
    with create_celeryconfig(task_module, dirs, config, config_file):
        __import__(task_module)
        tasks = sys.modules[task_module]
        from celery.task.sets import TaskSet
        def _run(fn_name, xs):
            fn = getattr(tasks, fn_name)
            job = TaskSet(tasks=[apply(fn.subtask, (x,)) for x in xs])
            result = job.apply_async()
            out = []
            if wait:
                while not result.ready():
                    time.sleep(5)
                for x in result.join():
                    if x:
                        out.extend(x)
            return out
        return _run

# ## Utility functions

_celeryconfig_tmpl = """
CELERY_IMPORTS = ("${task_import}", )

BROKER_HOST = "${host}"
BROKER_PORT = "${port}"
BROKER_USER = "${userid}"
BROKER_PASSWORD = "${password}"
BROKER_VHOST = "${rabbitmq_vhost}"
CELERY_RESULT_BACKEND= "amqp"
CELERY_TASK_SERIALIZER = "json"
CELERYD_CONCURRENCY = ${cores}
CELERY_ACKS_LATE = False
CELERYD_PREFETCH_MULTIPLIER = 1
CELERY_ROUTES = {"bcbio.distributed.tasks.analyze_and_upload": {"queue": "toplevel"},
                 "bcbio.distributed.tasks.long_term_storage" : {"queue": "storage"}}
BCBIO_CONFIG_FILE = "${config_file}"
"""

@contextlib.contextmanager
def create_celeryconfig(task_module, dirs, config, config_file):
    amqp_config = utils.read_galaxy_amqp_config(config["galaxy_config"], dirs["config"])
    if not amqp_config.has_key("host") or not amqp_config.has_key("userid"):
        raise ValueError("universe_wsgi.ini does not have RabbitMQ messaging details set")
    out_file = os.path.join(dirs["work"], "celeryconfig.py")
    amqp_config["rabbitmq_vhost"] = config["distributed"]["rabbitmq_vhost"]
    cores = config["distributed"].get("cores_per_host", 0)
    if cores < 1:
        cores = multiprocessing.cpu_count()
    amqp_config["cores"] = cores
    amqp_config["task_import"] = task_module
    amqp_config["config_file"] = config_file
    with open(out_file, "w") as out_handle:
        out_handle.write(Template(_celeryconfig_tmpl).render(**amqp_config))
    try:
        yield out_file
    finally:
        pyc_file = "%s.pyc" % os.path.splitext(out_file)[0]
        for fname in [pyc_file, out_file]:
            if os.path.exists(fname):
                os.remove(fname)
