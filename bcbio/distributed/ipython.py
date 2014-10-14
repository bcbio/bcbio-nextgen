"""Distributed execution using an IPython cluster.

Uses IPython parallel to setup a cluster and manage execution:

http://ipython.org/ipython-doc/stable/parallel/index.html

Cluster implementation from ipython-cluster-helper:

https://github.com/roryk/ipython-cluster-helper
"""
import os

from bcbio import utils
from bcbio.log import logger, get_log_dir
from bcbio.pipeline import config_utils
from bcbio.provenance import diagnostics

from cluster_helper import cluster as ipython_cluster

def create(parallel, dirs, config):
    """Create a cluster based on the provided parallel arguments.

    Returns an IPython view on the cluster, enabling processing on jobs.
    """
    profile_dir = utils.safe_makedir(os.path.join(dirs["work"], get_log_dir(config), "ipython"))
    return ipython_cluster.cluster_view(parallel["scheduler"].lower(), parallel["queue"],
                                        parallel["num_jobs"], parallel["cores_per_job"],
                                        profile=profile_dir, start_wait=parallel["timeout"],
                                        extra_params={"resources": parallel["resources"],
                                                      "mem": parallel["mem"],
                                                      "tag": parallel.get("tag"),
                                                      "run_local": parallel.get("run_local")},
                                        retries=parallel.get("retries"))

def _get_ipython_fn(fn_name, parallel):
    import_fn_name = parallel.get("wrapper", fn_name)
    return getattr(__import__("{base}.ipythontasks".format(base=parallel["module"]),
                              fromlist=["ipythontasks"]),
                   import_fn_name)

def runner(view, parallel, dirs, config):
    """Run a task on an ipython parallel cluster, allowing alternative queue types.

    view provides map-style access to an existing Ipython cluster.
    """
    def run(fn_name, items):
        out = []
        items = [x for x in items if x is not None]
        items = diagnostics.track_parallel(items, fn_name)
        fn = _get_ipython_fn(fn_name, parallel)
        logger.info("ipython: %s" % fn_name)
        if len(items) > 0:
            items = [config_utils.add_cores_to_config(x, parallel["cores_per_job"], parallel) for x in items]
            if "wrapper" in parallel:
                wrap_parallel = {k: v for k, v in parallel.items() if k in set(["fresources"])}
                items = [[fn_name] + parallel.get("wrapper_args", []) + [wrap_parallel] + list(x) for x in items]
            for data in view.map_sync(fn, items, track=False):
                if data:
                    out.extend(data)
        return out
    return run
