"""Distributed execution using an IPython cluster.

Uses IPython parallel to setup a cluster and manage execution:

http://ipython.org/ipython-doc/stable/parallel/index.html

Cluster implementation from ipython-cluster-helper:

https://github.com/roryk/ipython-cluster-helper
"""
import collections
import math
import os
import time

try:
    import ipyparallel
    # msgpack not working with IPython 4.0 ipyparallel
    msgpack = None
except ImportError:
    try:
        import msgpack
    except ImportError:
        msgpack = None

from bcbio import utils, setpath
from bcbio.log import logger, get_log_dir
from bcbio.pipeline import config_utils
from bcbio.provenance import diagnostics

from cluster_helper import cluster as ipython_cluster

def create(parallel, dirs, config):
    """Create a cluster based on the provided parallel arguments.

    Returns an IPython view on the cluster, enabling processing on jobs.

    Adds a mincores specification if he have machines with a larger
    number of cores to allow jobs to be batched together for shared
    memory usage.
    """
    profile_dir = utils.safe_makedir(os.path.join(dirs["work"], get_log_dir(config), "ipython"))
    has_mincores = any(x.startswith("mincores=") for x in parallel["resources"])
    cores = min(_get_common_cores(config["resources"]), parallel["system_cores"])
    if cores > 1 and not has_mincores:
        adj_cores = max(1, int(math.floor(cores * float(parallel.get("mem_pct", 1.0)))))
        # if we have less scheduled cores than per machine, use the scheduled count
        if cores > parallel["cores"]:
            cores = parallel["cores"]
        # if we have less total cores required for the entire process, use that
        elif adj_cores > parallel["num_jobs"] * parallel["cores_per_job"]:
            cores = parallel["num_jobs"] * parallel["cores_per_job"]
        else:
            cores = adj_cores
            cores = per_machine_target_cores(cores, parallel["num_jobs"])
        parallel["resources"].append("mincores=%s" % cores)
    return ipython_cluster.cluster_view(parallel["scheduler"].lower(), parallel["queue"],
                                        parallel["num_jobs"], parallel["cores_per_job"],
                                        profile=profile_dir, start_wait=parallel["timeout"],
                                        extra_params={"resources": parallel["resources"],
                                                      "mem": parallel["mem"],
                                                      "tag": parallel.get("tag"),
                                                      "run_local": parallel.get("run_local"),
                                                      "local_controller": parallel.get("local_controller")},
                                        retries=parallel.get("retries"))

def per_machine_target_cores(cores, num_jobs):
    """Select target cores on larger machines to leave room for batch script and controller.

    On resource constrained environments, we want to pack all bcbio submissions onto a specific
    number of machines. This gives up some cores to enable sharing cores with the controller
    and batch script on larger machines.
    """
    if cores >= 32 and num_jobs == 1:
        cores = cores - 2
    elif cores >= 16 and num_jobs in [1, 2]:
        cores = cores - 1
    return cores

def _get_common_cores(resources):
    """Retrieve the most common configured number of cores in the input file.
    """
    all_cores = []
    for vs in resources.values():
        cores = vs.get("cores")
        if cores:
            all_cores.append(int(vs["cores"]))
    return collections.Counter(all_cores).most_common(1)[0][0]

def stop(view):
    try:
        ipython_cluster.stop_from_view(view)
        time.sleep(10)
    except:
        logger.exception("Did not stop IPython cluster correctly")

def _get_ipython_fn(fn_name, parallel):
    import_fn_name = parallel.get("wrapper", fn_name)
    return getattr(__import__("{base}.ipythontasks".format(base=parallel["module"]),
                              fromlist=["ipythontasks"]),
                   import_fn_name)

def zip_args(args, config=None):
    """Compress arguments using msgpack.
    """
    if msgpack:
        return [msgpack.packb(x, use_single_float=True, use_bin_type=True) for x in args]
    else:
        return args

def unzip_args(args):
    """Uncompress arguments using msgpack.
    """
    if msgpack:
        return [msgpack.unpackb(x) for x in args]
    else:
        return args

def runner(view, parallel, dirs, config):
    """Run a task on an ipython parallel cluster, allowing alternative queue types.

    view provides map-style access to an existing Ipython cluster.
    """
    def run(fn_name, items):
        setpath.prepend_bcbiopath()
        out = []
        fn, fn_name = (fn_name, fn_name.__name__) if callable(fn_name) else (_get_ipython_fn(fn_name, parallel), fn_name)
        items = [x for x in items if x is not None]
        items = diagnostics.track_parallel(items, fn_name)
        logger.info("ipython: %s" % fn_name)
        if len(items) > 0:
            items = [config_utils.add_cores_to_config(x, parallel["cores_per_job"], parallel) for x in items]
            if "wrapper" in parallel:
                wrap_parallel = {k: v for k, v in parallel.items() if k in set(["fresources"])}
                items = [[fn_name] + parallel.get("wrapper_args", []) + [wrap_parallel] + list(x) for x in items]
            items = zip_args([args for args in items])
            for data in view.map_sync(fn, items, track=False):
                if data:
                    out.extend(unzip_args(data))
        return out
    return run
