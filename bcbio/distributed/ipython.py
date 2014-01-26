"""Distributed execution using an IPython cluster.

Uses IPython parallel to setup a cluster and manage execution:

http://ipython.org/ipython-doc/stable/parallel/index.html

Cluster implementation from ipython-cluster-helper:

https://github.com/roryk/ipython-cluster-helper
"""
import os
import contextlib
import copy

from bcbio import utils
from bcbio.log import logger, get_log_dir
from bcbio.pipeline import config_utils
from bcbio.provenance import system
from bcbio.distributed import resources

from cluster_helper import cluster as ipython_cluster

def dictadd(orig, k, v):
    """Imitates immutability by adding a key/value to a new dictionary.
    Works around not being able to deepcopy view objects.
    """
    view = orig.pop("view", None)
    new = copy.deepcopy(orig)
    new[k] = v
    if view:
        orig["view"] = view
        new["view"] = view
    return new

cur_num = 0
def _get_checkpoint_file(cdir, fn_name):
    """Retrieve checkpoint file for this step, with step number and function name.
    """
    global cur_num
    fname = os.path.join(cdir, "%s-%s.done" % (cur_num, fn_name))
    cur_num += 1
    return fname

def _view_from_parallel(parallel, work_dir, config):
    """Translate parallel map into options for a cluster view.
    """
    profile_dir = utils.safe_makedir(os.path.join(work_dir, get_log_dir(config), "ipython"))
    return ipython_cluster.cluster_view(parallel["scheduler"].lower(), parallel["queue"],
                                        parallel["num_jobs"], parallel["cores_per_job"],
                                        profile=profile_dir, start_wait=parallel["timeout"],
                                        extra_params={"resources": parallel["resources"],
                                                      "mem": parallel["mem"],
                                                      "tag": parallel.get("tag")},
                                        retries=parallel.get("retries"))

def _get_ipython_fn(fn_name, parallel):
    return getattr(__import__("{base}.ipythontasks".format(base=parallel["module"]),
                              fromlist=["ipythontasks"]),
                   fn_name)

@contextlib.contextmanager
def global_parallel(parallel, name, fn_names, items, dirs, config,
                    multiplier=1, max_multicore=None):
    """Add an IPython cluster to be used for multiple remote functions.

    Allows sharing of a single cluster across multiple functions with
    identical resource requirements. Falls back into local execution for
    non-distributed clusters or completed jobs.
    """
    checkpoint_dir = utils.safe_makedir(os.path.join(dirs["work"], "checkpoints_ipython"))
    checkpoint_file = os.path.join(checkpoint_dir, "global-%s.done" % name)
    sysinfo = system.get_info(dirs, parallel)
    try:
        if parallel["type"] != "ipython":
            parallel = copy.deepcopy(parallel)
            parallel["multiplier"] = multiplier
            parallel["max_multicore"] = max_multicore
            yield parallel
        elif os.path.exists(checkpoint_file):
            parallel["checkpoint"] = True
            yield parallel
        else:
            items = [x for x in items if x is not None]
            jobr = resources.calculate([_get_ipython_fn(x, parallel) for x in fn_names],
                                       parallel, items, sysinfo, config, multiplier=multiplier,
                                       max_multicore=max_multicore)
            parallel = dictadd(parallel, "cores_per_job", jobr.cores_per_job)
            parallel = dictadd(parallel, "num_jobs", jobr.num_jobs)
            parallel = dictadd(parallel, "mem", jobr.memory_per_job)
            with _view_from_parallel(parallel, dirs["work"], config) as view:
                parallel["checkpoint"] = False
                parallel["view"] = view
                yield parallel
    except:
        raise
    else:
        parallel["view"] = None
        parallel["checkpoint"] = False
        with open(checkpoint_file, "w") as out_handle:
            out_handle.write("done\n")

def runner(parallel, fn_name, items, work_dir, sysinfo, config):
    """Run a task on an ipython parallel cluster, allowing alternative queue types.

    This will spawn clusters for parallel and custom queue types like multicore
    and high I/O tasks on demand.

    A checkpoint directory keeps track of finished tasks, avoiding spinning up clusters
    for sections that have been previous processed.

    The parallel input function can contain a pre-created view on an
    existing Ipython cluster, in which case this will be reused
    instead of creating a new cluster.
    """
    out = []
    items = [x for x in items if x is not None]
    algs = [config_utils.get_algorithm_config(x) for x in items]
    if ((len(algs) > 0 and not algs[0].get("resource_check", True))
          or parallel.get("view") or parallel.get("checkpoint")):
        checkpoint_file = None
    else:
        checkpoint_dir = utils.safe_makedir(os.path.join(work_dir, "checkpoints_ipython"))
        checkpoint_file = _get_checkpoint_file(checkpoint_dir, fn_name)
    fn = _get_ipython_fn(fn_name, parallel)
    if parallel.get("view") is None:
        jobr = resources.calculate([fn], parallel, items, sysinfo, config)
        parallel = dictadd(parallel, "cores_per_job", jobr.cores_per_job)
        parallel = dictadd(parallel, "num_jobs", jobr.num_jobs)
        parallel = dictadd(parallel, "mem", jobr.memory_per_job)
    # already finished, run locally on current machine to collect details
    if parallel.get("checkpoint") or (checkpoint_file and os.path.exists(checkpoint_file)):
        logger.info("ipython: %s -- local; checkpoint passed" % fn_name)
        for args in items:
            if args:
                data = fn(args)
                if data:
                    out.extend(data)
    # Run on a standard parallel queue
    else:
        logger.info("ipython: %s" % fn_name)
        if len(items) > 0:
            items = [config_utils.add_cores_to_config(x, parallel["cores_per_job"], parallel) for x in items]
            if parallel.get("view"):
                for data in parallel["view"].map_sync(fn, items, track=False):
                    if data:
                        out.extend(data)
            else:
                with _view_from_parallel(parallel, work_dir, config) as view:
                    for data in view.map_sync(fn, items, track=False):
                        if data:
                            out.extend(data)
    if checkpoint_file:
        with open(checkpoint_file, "w") as out_handle:
            out_handle.write("done\n")
    return out
