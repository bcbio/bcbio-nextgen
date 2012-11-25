"""Distributed execution using an IPython cluster.

Uses IPython parallel to setup a cluster and manage execution:

http://ipython.org/ipython-doc/stable/parallel/index.html

Borrowed from Rory Kirchner's Bipy cluster implementation:

https://github.com/roryk/bipy/blob/master/bipy/cluster/__init__.py
"""
import copy
import time
import uuid
import subprocess
import contextlib

from IPython.parallel import Client

def _start(workers_needed, profile, cluster_id, delay):
    """Starts cluster from commandline.
    """
    subprocess.check_call(["ipcluster", "start",
                           "--daemonize=True",
                           "--delay=%s" % delay, 
                           "--log-level=%s" % "WARN",
                           #"--cluster-id=%s" % cluster_id,
                           "--n=%s" % workers_needed,
                           "--profile=%s" % profile])

def _stop(profile, cluster_id):
    subprocess.check_call(["ipcluster", "stop", "--profile=%s" % profile,
                           #"--cluster-id=%s" % cluster_id
                           ])

def _is_up(profile, cluster_id, n):
    try:
        #client = Client(profile=profile, cluster_id=cluster_id)
        client = Client(profile=profile)
        up = len(client.ids)
    except IOError, msg:
        return False
    else:
        return up >= n

@contextlib.contextmanager
def cluster_view(parallel):
    """Provide a view on an ipython cluster for processing.

    parallel is a dictionary with:
      - profile: The name of the ipython profile to use
      - cores: The number of cores to start for processing.
      - queue_type: Optionally, the type of parallel queue
        to start. Defaults to a standard parallel queue, can
        also specify 'multicore' for a multiple core machine
        and 'io' for an I/O intensive queue.
    """
    delay = 10
    max_delay = 300
    profile = parallel["profile"]
    if parallel.get("queue_type", None):
        profile = "%s_%s" % (profile, parallel["queue_type"])
    cluster_id = str(uuid.uuid1())
    # need at least two processes to run main and workers
    _start(parallel["cores"], profile, cluster_id, delay)
    try:
        slept = 0
        target_cores = 1 if parallel.get("queue_type", None) == "multicore" \
                       else parallel["cores"]
        while not _is_up(profile, cluster_id, target_cores):
            time.sleep(delay)
            slept += delay
            if slept > max_delay:
                raise IOError("Cluster startup timed out.")
        #client = Client(profile=profile, cluster_id=cluster_id)
        client = Client(profile=profile)
        yield client.load_balanced_view()
    finally:
        _stop(profile, cluster_id)

def dictadd(orig, k, v):
    """Imitates immutability by adding a key/value to a new dictionary.
    Works around not being able to deepcopy view objects; can remove this
    once we create views on demand.
    """
    view = orig.pop("view", None)
    new = copy.deepcopy(orig)
    new[k] = v
    if view:
        orig["view"] = view
        new["view"] = view
    return new

def _get_queue_type(fn):
    if hasattr(fn, "metadata"):
        return fn.metadata.get("queue_type", None)
    else:
        return None

def runner(parallel, fn_name, items):
    """Run a task on an ipython parallel cluster, allowing alternative queue types.

    This will spawn clusters for custom queue types like multicore and high I/O
    tasks on demand.
    TODO: spawn standard queues on demand as well.
    """
    out = []
    fn = getattr(__import__("{base}.ipythontasks".format(base=parallel["module"]),
                            fromlist=["ipythontasks"]),
                 fn_name)
    queue_type = _get_queue_type(fn)
    if queue_type:
        parallel = dictadd(parallel, "queue_type", queue_type)
    if queue_type == "multicore":
        with cluster_view(parallel) as view:
            for args in items:
                if args:
                    data = view.apply_sync(fn, args)
                    if data:
                        out.append(data)
    else:
        xs = [x for x in items if x is not None]
        if len(xs) > 0:
            for data in parallel["view"].map_sync(fn, xs):
                if data:
                    out.extend(data)
    return out
