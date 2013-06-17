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
from bcbio.log import logger
from bcbio.pipeline import config_utils

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

def _get_resource_programs(fn, algs):
    """Retrieve programs used in analysis based on algorithm configurations.

    Helps avoid requiring core information from unused programs.
    """
    # standard list of programs we always use
    used_progs = set(["gatk", "gemini"])
    for alg in algs:
        # get aligners used
        aligner = alg.get("aligner")
        if aligner:
            used_progs.add(aligner)
    for prog in (fn.metadata.get("resources", []) if hasattr(fn, "metadata") else []):
        if prog in used_progs:
            yield prog

def find_cores_per_job(fns, parallel, items, config):
    """Determine cores and workers to use for this stage based on function metadata.
    """
    all_cores = [1]
    algs = [get_algorithm_config(x) for x in items]
    for fn in fns:
        for prog in _get_resource_programs(fn, algs):
            resources = config_utils.get_resources(prog, config)
            cores = resources.get("cores")
            if cores:
                all_cores.append(cores)
    cores_per_job = max(all_cores)
    total = parallel["cores"]
    if total > cores_per_job:
        return min(total // cores_per_job, len(items)), cores_per_job
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
    return isinstance(x, dict) and "algorithm" in x and "resources" in x and not "files" in x

def is_nested_config_arg(x):
    return isinstance(x, dict) and "config" in x and is_std_config_arg(x["config"])

def get_algorithm_config(xs):
    """Flexibly extract algorithm configuration for a sample from any function arguments.
    """
    for x in xs:
        if is_std_config_arg(x):
            return x["algorithm"]
        elif is_nested_config_arg(x):
            return x["config"]["algorithm"]
    raise ValueError("Did not find algorithm configuration in items: {0}"
                     .format(xs))

def add_cores_to_config(args, cores_per_job, parallel=None):
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
        if parallel:
            new_arg["config"]["parallel"] = parallel
    elif is_std_config_arg(new_arg):
        new_arg["algorithm"]["num_cores"] = int(cores_per_job)
        if parallel:
            new_arg["parallel"] = parallel
    else:
        raise ValueError("Unexpected configuration dictionary: %s" % new_arg)
    args = list(args)[:]
    args[new_i] = new_arg
    return args

def _view_from_parallel(parallel):
    """Translate parallel map into options for a cluster view.
    """
    return ipython_cluster.cluster_view(parallel["scheduler"].lower(), parallel["queue"],
                                        parallel["num_jobs"], parallel["cores_per_job"],
                                        profile=parallel["profile"], start_wait=parallel["timeout"],
                                        extra_params={"resources": parallel["resources"]},
                                        retries=parallel.get("retries"))

def _get_ipython_fn(fn_name, parallel):
    return getattr(__import__("{base}.ipythontasks".format(base=parallel["module"]),
                            fromlist=["ipythontasks"]),
                   fn_name)

@contextlib.contextmanager
def global_parallel(parallel, name, fn_names, items, work_dir, config):
    """Add an IPython cluster to be used for multiple remote functions.

    Allows sharing of a single cluster across multiple functions with
    identical resource requirements. Falls back into local execution for
    non-distributed clusters or completed jobs.
    """
    checkpoint_dir = utils.safe_makedir(os.path.join(work_dir, "checkpoints_ipython"))
    checkpoint_file = os.path.join(checkpoint_dir, "global-%s.done" % name)
    try:
        if parallel["type"] != "ipython" or os.path.exists(checkpoint_file):
            yield parallel
        else:
            items = [x for x in items if x is not None]
            num_jobs, cores_per_job = find_cores_per_job([_get_ipython_fn(x) for x in fn_names],
                                                         parallel, items, config)
            parallel = dictadd(parallel, "cores_per_job", cores_per_job)
            parallel = dictadd(parallel, "num_jobs", num_jobs)
            with _view_from_parallel(parallel) as view:
                parallel["view"] = view
                yield parallel
    except:
        raise
    else:
        with open(checkpoint_file, "w") as out_handle:
            out_handle.write("done\n")

def runner(parallel, fn_name, items, work_dir, config):
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
    checkpoint_dir = utils.safe_makedir(os.path.join(work_dir, "checkpoints_ipython"))
    checkpoint_file = _get_checkpoint_file(checkpoint_dir, fn_name)
    fn = _get_ipython_fn(fn_name, parallel)
    items = [x for x in items if x is not None]
    if parallel.get("view") is None:
        num_jobs, cores_per_job = find_cores_per_job([fn], parallel, items, config)
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
            items = [add_cores_to_config(x, parallel["cores_per_job"], parallel) for x in items]
            if parallel.get("view"):
                for data in parallel["view"].map_sync(fn, items, track=False):
                    if data:
                        out.extend(data)
            else:
                with _view_from_parallel(parallel) as view:
                    for data in view.map_sync(fn, items, track=False):
                        if data:
                            out.extend(data)
    with open(checkpoint_file, "w") as out_handle:
        out_handle.write("done\n")
    return out
