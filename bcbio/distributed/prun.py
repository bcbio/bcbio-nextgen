"""Generalized running of parallel tasks in multiple environments.
"""
import contextlib
import os

from bcbio import utils
from bcbio.log import logger
from bcbio.provenance import system
from bcbio.distributed import multi, resources

@contextlib.contextmanager
def start(parallel, items, config, dirs=None, name=None, multiplier=1,
          max_multicore=None):
    """Start a parallel cluster or machines to be used for running remote
    functions.

    Returns a function used to process, in parallel items with a given function.

    Allows sharing of a single cluster across multiple functions with
    identical resource requirements. Uses local execution for non-distributed
    clusters or completed jobs.

    A checkpoint directory keeps track of finished tasks, avoiding spinning up
    clusters for sections that have been previous processed.

    multiplier - Number of expected jobs per initial input item. Used to avoid
    underscheduling cores when an item is split during processing.
    max_multicore -- The maximum number of cores to use for each process. Can be
    used to process less multicore usage when jobs run faster on more single
    cores.
    """
    if name:
        checkpoint_dir = utils.safe_makedir(os.path.join(dirs["work"],
                                                         "checkpoints_parallel"))
        checkpoint_file = os.path.join(checkpoint_dir, "%s.done" % name)
    else:
        checkpoint_file = None
    sysinfo = system.get_info(dirs, parallel, config.get("resources", {}))
    items = [x for x in items if x is not None] if items else []
    max_multicore = int(max_multicore or sysinfo.get("cores", 1))
    parallel = resources.calculate(parallel, items, sysinfo, config,
                                   multiplier=multiplier,
                                   max_multicore=max_multicore)
    try:
        view = None
        if parallel["type"] == "ipython":
            if checkpoint_file and os.path.exists(checkpoint_file):
                logger.info("Running locally instead of distributed -- checkpoint passed: %s" % name)
                parallel["cores_per_job"] = 1
                parallel["num_jobs"] = 1
                parallel["checkpointed"] = True
                yield multi.runner(parallel, config)
            else:
                from bcbio.distributed import ipython
                with ipython.create(parallel, dirs, config) as view:
                    yield ipython.runner(view, parallel, dirs, config)
        else:
            yield multi.runner(parallel, config)
    except:
        if view is not None:
            from bcbio.distributed import ipython
            ipython.stop(view)
        raise
    else:
        for x in ["cores_per_job", "num_jobs", "mem"]:
            parallel.pop(x, None)
        if checkpoint_file:
            with open(checkpoint_file, "w") as out_handle:
                out_handle.write("done\n")
