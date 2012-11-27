"""Run distributed tasks using the Celery distributed task queue.

http://celeryproject.org/
"""

import os
import sys
import time
import contextlib
import multiprocessing
import subprocess

from mako.template import Template

from bcbio import utils
from bcbio.distributed import ipython

def parallel_runner(parallel, dirs, config, config_file):
    """Process a supplied function: single, multi-processor or distributed.
    """
    def run_parallel(fn_name, items, metadata=None):
        if parallel["type"].startswith("messaging"):
            task_module = "{base}.tasks".format(base=parallel["module"])
            runner_fn = runner(task_module, dirs, config, config_file)
            return runner_fn(fn_name, items)
        elif parallel["type"] == "ipython":
            return ipython.runner(parallel, fn_name, items, dirs["work"], config)
        else:
            out = []
            fn = getattr(__import__("{base}.multitasks".format(base=parallel["module"]),
                                    fromlist=["multitasks"]),
                         fn_name)
            cores = cores_including_resources(int(parallel["cores"]), metadata, config)
            with utils.cpmap(cores) as cpmap:
                for data in cpmap(fn, filter(lambda x: x is not None, items)):
                    if data:
                        out.extend(data)
            return out
    return run_parallel

def runner(task_module, dirs, config, config_file, wait=True):
    """Run a set of tasks using Celery, waiting for results or asynchronously.

    Initialize with the configuration and directory information,
    used to prepare a Celery configuration file and imports. It
    returns a function which acts like standard map; provide the function
    name instead of the function itself when calling.

    After name lookup, Celery runs the function in parallel; Celery servers
    can be remote or local but must have access to a shared filesystem. The
    function polls if wait is True, returning when all results are available.
    """
    with create_celeryconfig(task_module, dirs, config, config_file):
        sys.path.append(dirs["work"])
        __import__(task_module)
        tasks = sys.modules[task_module]
        from celery.task.sets import TaskSet
        def _run(fn_name, xs):
            fn = getattr(tasks, fn_name)
            job = TaskSet(tasks=[apply(fn.subtask, (x,)) for x in xs])
            result = job.apply_async()
            out = []
            if wait:
                with _close_taskset(result):
                    while not result.ready():
                        time.sleep(5)
                        if result.failed():
                            raise ValueError("Failed distributed task; cleaning up")
                    for x in result.join():
                        if x:
                            out.extend(x)
            return out
        return _run

@contextlib.contextmanager
def _close_taskset(ts):
    """Revoke existing jobs if a taskset fails; raise original error.
    """
    try:
        yield None
    except:
        try:
            raise
        finally:
            try:
                ts.revoke()
            except:
                pass

# ## Handle memory bound processes on multi-core machines

def cores_including_resources(cores, metadata, config):
    """Retrieve number of cores to use, considering program resources.
    """
    if metadata is None: metadata = {}
    required_memory = -1
    for program in metadata.get("programs", []):
        presources = config.get("resources", {}).get(program, {})
        memory = presources.get("memory", None)
        if memory:
            if memory.endswith("g"):
                memory = int(memory[:-1])
            else:
                raise NotImplementedError("Unpexpected units on memory: %s", memory)
            if memory > required_memory:
                required_memory = memory
    if required_memory > 0:
        cur_memory = _machine_memory()
        cores = min(cores,
                    int(round(float(cur_memory) / float(required_memory))))
    if cores < 1:
        cores = 1
    return cores

def _machine_memory():
    """Retrieve available memory on current machine using 'free.'
    """
    with contextlib.closing(subprocess.Popen(["free", "-g"],
                                             stdout=subprocess.PIPE).stdout) as stdout:
        for line in stdout:
            if line.startswith("Mem:"):
                return int(line.split()[1])

# ## Utility functions

_celeryconfig_tmpl = """
CELERY_IMPORTS = ("${task_import}", )

BROKER_URL = "amqp://${userid}:${password}@${host}:${port}/${rabbitmq_vhost}"
CELERY_RESULT_BACKEND= "amqp"
CELERY_TASK_SERIALIZER = "json"
CELERYD_CONCURRENCY = ${cores}
CELERY_ACKS_LATE = False
CELERYD_PREFETCH_MULTIPLIER = 1
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
