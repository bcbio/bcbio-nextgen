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

def runner(dirs, config):
    """Run a set of tasks asynchronously using Celery.

    This is initialized with the configuration and directory information,
    which is used to prepare a Celery configuration file and imports. It
    returns a function which acts like standard map, except that the
    function name is provided instead of the function itself.

    The name is looked up and the function is run in parallel on Celery
    servers, which can be remotely located but are assumed to have access
    to the same filesystem. We will poll and wait until all results are
    ready, returning them.
    """
    task_module = "bcbio.distributed.tasks"
    with create_celeryconfig(task_module, dirs, config):
        __import__(task_module)
        tasks = sys.modules[task_module]
        from celery.task.sets import TaskSet
        def _run(fn_name, xs):
            fn = getattr(tasks, fn_name)
            job = TaskSet(tasks=[apply(fn.subtask, (x,)) for x in xs])
            result = job.apply_async()
            while not result.ready():
                time.sleep(5)
            out = []
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
"""

@contextlib.contextmanager
def create_celeryconfig(task_module, dirs, config):
    amqp_config = utils.read_galaxy_amqp_config(config["galaxy_config"], dirs["config"])
    out_file = os.path.join(dirs["work"], "celeryconfig.py")
    amqp_config["rabbitmq_vhost"] = config["distributed"]["rabbitmq_vhost"]
    amqp_config["cores"] = multiprocessing.cpu_count()
    amqp_config["task_import"] = task_module
    with open(out_file, "w") as out_handle:
        out_handle.write(Template(_celeryconfig_tmpl).render(**amqp_config))
    try:
        yield out_file
    finally:
        pyc_file = "%s.pyc" % os.path.splitext(out_file)[0]
        for fname in [pyc_file, out_file]:
            if os.path.exists(fname):
                os.remove(fname)
