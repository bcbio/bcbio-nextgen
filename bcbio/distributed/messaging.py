"""Run distributed tasks in parallel using IPython or joblib on multiple cores.
"""
import functools

try:
    import joblib
except ImportError:
    joblib = False

from bcbio.distributed import ipython
from bcbio.log import logger, setup_local_logging
from bcbio.provenance import diagnostics, system

def parallel_runner(parallel, dirs, config):
    """Process a supplied function: single, multi-processor or distributed.
    """
    def run_parallel(fn_name, items, metadata=None):
        items = [x for x in items if x is not None]
        if len(items) == 0:
            return []
        items = diagnostics.track_parallel(items, fn_name)
        sysinfo = system.get_info(dirs, parallel)
        if parallel["type"] == "ipython":
            return ipython.runner(parallel, fn_name, items, dirs["work"], sysinfo, config)
        else:
            imodule = parallel.get("module", "bcbio.distributed")
            logger.info("multiprocessing: %s" % fn_name)
            fn = getattr(__import__("{base}.multitasks".format(base=imodule),
                                    fromlist=["multitasks"]),
                         fn_name)
            return run_multicore(fn, items, config, parallel["cores"])
    return run_parallel

def zeromq_aware_logging(f):
    """Ensure multiprocessing logging uses ZeroMQ queues.

    ZeroMQ and local stdout/stderr do not behave nicely when intertwined. This
    ensures the local logging uses existing ZeroMQ logging queues.
    """
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        config = None
        for arg in args:
            if ipython.is_std_config_arg(arg):
                config = arg
                break
            elif ipython.is_nested_config_arg(arg):
                config = arg["config"]
                break
        assert config, "Could not find config dictionary in function arguments."
        if config.get("parallel", {}).get("log_queue"):
            handler = setup_local_logging(config, config["parallel"])
        else:
            handler = None
        try:
            out = f(*args, **kwargs)
        finally:
            if handler and hasattr(handler, "close"):
                handler.close()
        return out
    return wrapper

def run_multicore(fn, items, config, cores=None):
    """Run the function using multiple cores on the given items to process.
    """
    if cores is None:
        cores = config["algorithm"].get("num_cores", 1)
    parallel = {"type": "local", "cores": cores}
    sysinfo = system.get_info({}, parallel)
    jobr = ipython.find_job_resources([fn], parallel, items, sysinfo, config,
                                      parallel.get("multiplier", 1),
                                      max_multicore=int(sysinfo["cores"]))
    items = [ipython.add_cores_to_config(x, jobr.cores_per_job) for x in items]
    if joblib is None:
        raise ImportError("Need joblib for multiprocessing parallelization")
    out = []
    for data in joblib.Parallel(jobr.num_jobs)(joblib.delayed(fn)(x) for x in items):
        if data:
            out.extend(data)
    return out
