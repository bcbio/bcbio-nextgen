"""Run tasks in parallel on a single machine using multiple cores.
"""
import functools

try:
    import joblib
except ImportError:
    joblib = False

from bcbio.distributed import resources
from bcbio.log import logger, setup_local_logging
from bcbio.pipeline import config_utils
from bcbio.provenance import diagnostics, system

def runner(parallel, config):
    """Run functions, provided by string name, on multiple cores on the current machine.
    """
    def run_parallel(fn_name, items):
        items = [x for x in items if x is not None]
        if len(items) == 0:
            return []
        items = diagnostics.track_parallel(items, fn_name)
        fn, fn_name = (fn_name, fn_name.__name__) if callable(fn_name) else (get_fn(fn_name, parallel), fn_name)
        logger.info("multiprocessing: %s" % fn_name)
        if "wrapper" in parallel:
            wrap_parallel = {k: v for k, v in parallel.items() if k in set(["fresources", "checkpointed"])}
            items = [[fn_name] + parallel.get("wrapper_args", []) + [wrap_parallel] + list(x) for x in items]
        return run_multicore(fn, items, config, parallel=parallel)
    return run_parallel

def get_fn(fn_name, parallel):
    taskmod = "multitasks"
    imodule = parallel.get("module", "bcbio.distributed")
    import_fn_name = parallel.get("wrapper", fn_name)
    return getattr(__import__("{base}.{taskmod}".format(base=imodule, taskmod=taskmod),
                              fromlist=[taskmod]),
                   import_fn_name)

def zeromq_aware_logging(f):
    """Ensure multiprocessing logging uses ZeroMQ queues.

    ZeroMQ and local stdout/stderr do not behave nicely when intertwined. This
    ensures the local logging uses existing ZeroMQ logging queues.
    """
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        config = None
        for arg in args:
            if config_utils.is_std_config_arg(arg):
                config = arg
                break
            elif config_utils.is_nested_config_arg(arg):
                config = arg["config"]
            elif isinstance(arg, (list, tuple)) and config_utils.is_nested_config_arg(arg[0]):
                config = arg[0]["config"]
                break
        assert config, "Could not find config dictionary in function arguments."
        if config.get("parallel", {}).get("log_queue") and not config.get("parallel", {}).get("wrapper"):
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

def run_multicore(fn, items, config, parallel=None):
    """Run the function using multiple cores on the given items to process.
    """
    if len(items) == 0:
        return []
    if parallel is None or "num_jobs" not in parallel:
        if parallel is None:
            parallel = {"type": "local", "cores": config["algorithm"].get("num_cores", 1)}
        sysinfo = system.get_info({}, parallel)
        parallel = resources.calculate(parallel, items, sysinfo, config,
                                       parallel.get("multiplier", 1),
                                       max_multicore=int(parallel.get("max_multicore", sysinfo["cores"])))
    items = [config_utils.add_cores_to_config(x, parallel["cores_per_job"]) for x in items]
    if joblib is None:
        raise ImportError("Need joblib for multiprocessing parallelization")
    out = []
    for data in joblib.Parallel(parallel["num_jobs"], batch_size=1, backend="multiprocessing")(joblib.delayed(fn)(*x) for x in items):
        if data:
            out.extend(data)
    return out
