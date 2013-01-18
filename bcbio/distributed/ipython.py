"""Distributed execution using an IPython cluster.

Uses IPython parallel to setup a cluster and manage execution:

http://ipython.org/ipython-doc/stable/parallel/index.html

Borrowed from Rory Kirchner's Bipy cluster implementation:

https://github.com/roryk/bipy/blob/master/bipy/cluster/__init__.py
"""
import os
import copy
import time
import uuid
import subprocess
import contextlib

from bcbio import utils
from bcbio.log import setup_logging, logger

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
def cluster_view(parallel, config):
    """Provide a view on an ipython cluster for processing.

    parallel is a dictionary with:
      - profile: The name of the ipython profile to use
      - cores: The number of cores to start for processing.
      - queue_type: Optionally, the type of parallel queue
        to start. Defaults to a standard parallel queue, can
        also specify 'multicore' for a multiple core machine
        and 'io' for an I/O intensive queue.
    """
    delay = 5
    max_delay = 300
    max_tries = 10
    profile = parallel["profile"]
    if parallel.get("queue_type", None):
        profile = "%s_%s" % (profile, parallel["queue_type"])
    cluster_id = str(uuid.uuid1())
    num_tries = 0
    while 1:
        try:
            _start(parallel["cores"], profile, cluster_id, delay)
            break
        except subprocess.CalledProcessError:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(delay)
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
        # push config to all engines and force them to set up logging
        client[:]['config'] = config
        client[:].execute('from bcbio.log import setup_logging')
        client[:].execute('setup_logging(config)')
        client[:].execute('from bcbio.log import logger')
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

def runner(parallel, fn_name, items, work_dir, config):
    """Run a task on an ipython parallel cluster, allowing alternative queue types.

    This will spawn clusters for parallel and custom queue types like multicore
    and high I/O tasks on demand.

    A checkpoint directory keeps track of finished tasks, avoiding spinning up clusters
    for sections that have been previous processed.
    """
    setup_logging(config)
    out = []
    checkpoint_dir = utils.safe_makedir(os.path.join(work_dir, "checkpoints_ipython"))
    checkpoint_file = os.path.join(checkpoint_dir, "%s.done" % fn_name)
    fn = getattr(__import__("{base}.ipythontasks".format(base=parallel["module"]),
                            fromlist=["ipythontasks"]),
                 fn_name)
    queue_type = _get_queue_type(fn)
    if queue_type:
        parallel = dictadd(parallel, "queue_type", queue_type)
    # already finished, run locally on current machine to collect details
    if os.path.exists(checkpoint_file):
        logger.info("ipython: %s -- local; checkpoint passed" % fn_name)
        for args in items:
            if args:
                data = fn(args)
                if data:
                    out.extend(data)
    # Run on a multicore queue with available cores on the same machine
    elif queue_type == "multicore":
        logger.info("ipython: %s -- multicore" % fn_name)
        with cluster_view(parallel, config) as view:
            for args in items:
                if args:
                    data = view.apply_sync(fn, args)
                    if data:
                        out.extend(data)
    # Run on a standard parallel queue
    else:
        logger.info("ipython: %s -- parallel" % fn_name)
        with cluster_view(parallel, config) as view:
            xs = [x for x in items if x is not None]
            if len(xs) > 0:
                for data in view.map_sync(fn, xs):
                    if data:
                        out.extend(data)
    with open(checkpoint_file, "w") as out_handle:
        out_handle.write("done\n")
    return out
