"""Distributed execution using an IPython cluster.

Uses IPython parallel to setup a cluster and manage execution:

http://ipython.org/ipython-doc/stable/parallel/index.html

Cluster implementation from ipython-cluster-helper:

https://github.com/roryk/ipython-cluster-helper
"""
import os
import copy

from bcbio import utils
from bcbio.log import setup_logging, logger
from bcbio.pipeline import config_utils

from cluster_helper import cluster as ipython_cluster

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

def _find_cores_per_job(fn, parallel, item_count, config):
    """Determine cores and workers to use for this stage based on function metadata.
    """
    all_cores = [1]
    for prog in (fn.metadata.get("resources", []) if hasattr(fn, "metadata") else []):
        resources = config_utils.get_resources(prog, config)
        cores = resources.get("cores")
        if cores:
            all_cores.append(cores)
    cores_per_job = max(all_cores)
    total = parallel["cores"]
    if total > cores_per_job:
        return min(total // cores_per_job, item_count), cores_per_job
    else:
        return 1, total

cur_num = 0
def _get_checkpoint_file(cdir, fn_name):
    """Retrieve checkpoint file for this step, with step number and function name.
    """
    global cur_num
    fname = os.path.join(cdir, "%s-%s.done" % (cur_num, fn_name))
    cur_num += 1
    return fname

def is_std_config_arg(x):
    return (isinstance(x, dict) and x.has_key("algorithm") and x.has_key("resources"))

def is_nested_config_arg(x):
    return (isinstance(x, dict) and x.has_key("config") and
            is_std_config_arg(x["config"]))

def add_cores_to_config(args, cores_per_job):
    """Add information about available cores for a job to configuration.
    Ugly hack to update core information in a configuration dictionary.
    """
    new_i = None
    for i, arg in enumerate(args):
        if is_std_config_arg(arg) or is_nested_config_arg(arg):
            new_i = i
            break
    if new_i is None:
        raise ValueError("Could not find configuration in args: %s" % args)

    new_arg = copy.deepcopy(args[new_i])
    if is_nested_config_arg(new_arg):
        new_arg["config"]["algorithm"]["num_cores"] = int(cores_per_job)
    elif is_std_config_arg(new_arg):
        new_arg["algorithm"]["num_cores"] = int(cores_per_job)
    else:
        raise ValueError("Unexpected configuration dictionary: %s" % new_arg)
    args = list(args)[:]
    args[new_i] = new_arg
    return args

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
    checkpoint_file = _get_checkpoint_file(checkpoint_dir, fn_name)
    fn = getattr(__import__("{base}.ipythontasks".format(base=parallel["module"]),
                            fromlist=["ipythontasks"]),
                 fn_name)
    items = [x for x in items if x is not None]
    num_jobs, cores_per_job = _find_cores_per_job(fn, parallel, len(items), config)
    parallel = dictadd(parallel, "cores_per_job", cores_per_job)
    parallel = dictadd(parallel, "num_jobs", num_jobs)
    # already finished, run locally on current machine to collect details
    if os.path.exists(checkpoint_file):
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
            items = [add_cores_to_config(x, cores_per_job) for x in items]
            with ipython_cluster.cluster_view(parallel["scheduler"].lower(), parallel["queue"],
                                              parallel["num_jobs"], parallel["cores_per_job"],
                                              profile=parallel["profile"]) as view:
                for data in view.map_sync(fn, items, track=False):
                    if data:
                        out.extend(data)
    with open(checkpoint_file, "w") as out_handle:
        out_handle.write("done\n")
    return out
