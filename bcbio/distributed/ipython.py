"""Distributed execution using an IPython cluster.

Uses IPython parallel to setup a cluster and manage execution:

http://ipython.org/ipython-doc/stable/parallel/index.html

Cluster implementation from ipython-cluster-helper:

https://github.com/roryk/ipython-cluster-helper
"""
import os
import collections
import contextlib
import copy
import math

from bcbio import utils
from bcbio.log import logger, get_log_dir
from bcbio.pipeline import config_utils
from bcbio.provenance import system

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

def dictdissoc(orig, k):
    """Imitates immutability: create a new dictionary with the key dropped.
    """
    v = orig.pop(k, None)
    new = copy.deepcopy(orig)
    orig[k] = v
    return new

def _get_resource_programs(fn, algs):
    """Retrieve programs used in analysis based on algorithm configurations.

    Helps avoid requiring core information from unused programs.
    """
    used_progs = _get_used_programs(fn, algs)
    # standard list of programs we always use
    # XXX Need to expose this in a top-level way to allow more multiprocessing
    for prog in (fn.metadata.get("resources", []) if hasattr(fn, "metadata") else []):
        if prog in used_progs:
            yield prog

def _get_ensure_functions(fn, algs):
    used_progs = _get_used_programs(fn, algs)
    for prog in (fn.metadata.get("ensure", {}).keys() if hasattr(fn, "metadata") else []):
        if prog in used_progs:
            yield fn.metadata["ensure"][prog]

def _get_used_programs(fn, algs):
    used_progs = set(["gatk", "gemini", "bcbio_coverage", "samtools",
                      "snpEff", "cufflinks", "picard"])
    for alg in algs:
        # get aligners used
        aligner = alg.get("aligner")
        if aligner:
            used_progs.add(aligner)
        vc = alg.get("variantcaller")
        if vc:
            if isinstance(vc, (list, tuple)):
                for x in vc:
                    used_progs.add(x)
            else:
                used_progs.add(vc)
    if config_utils.use_vqsr(algs):
        used_progs.add("gatk-vqsr")
    return used_progs


def _str_memory_to_gb(memory):
    val = float(memory[:-1])
    units = memory[-1]
    if units.lower() == "m":
        val = val / 1000.0
    else:
        assert units.lower() == "g", "Unexpected memory units: %s" % memory
    return val


def _get_prog_memory(resources):
    """Get expected memory usage, in Gb per core, for a program from resource specification.
    """
    out = None
    for jvm_opt in resources.get("jvm_opts", []):
        if jvm_opt.startswith("-Xmx"):
            out = _str_memory_to_gb(jvm_opt[4:])
    memory = resources.get("memory")
    if memory:
        out = _str_memory_to_gb(memory)
    return out


def _scale_cores_to_memory(cores, mem_per_core, sysinfo, system_memory):
    """Scale multicore usage to avoid excessive memory usage based on system information.
    """
    total_mem = "%.1f" % (cores * mem_per_core + system_memory)
    if "cores" not in sysinfo:
        return cores, total_mem
    cores = min(cores, int(sysinfo["cores"]))
    total_mem = min(float(total_mem), float(sysinfo["memory"]) - system_memory)
    cores = max(1, min(cores, int(math.floor(float(total_mem) / mem_per_core))))
    return cores, total_mem

def _scale_jobs_to_memory(jobs, mem_per_core, sysinfo):
    """When scheduling jobs with single cores, avoid overscheduling due to memory.
    """
    if "cores" not in sysinfo:
        return jobs
    sys_mem_per_core = float(sysinfo["memory"]) / float(sysinfo["cores"])
    if sys_mem_per_core < mem_per_core:
        pct = sys_mem_per_core / float(mem_per_core)
        target_jobs = int(math.floor(jobs * pct))
        return max(target_jobs, 1)
    else:
        return jobs

def find_job_resources(fns, parallel, items, sysinfo, config, multiplier=1,
                       max_multicore=None):
    """Determine cores and workers to use for this stage based on function metadata.
    multiplier specifies the number of regions items will be split into during
    processing.
    max_multicore specifies an optional limit on the maximum cores. Can use to
    force single core processing during specific tasks.
    sysinfo specifies cores and memory on processing nodes, allowing us to tailor
    jobs for available resources.
    """
    assert len(items) > 0, "Finding job resources but no items to process"
    all_cores = [1]
    # Use modest 2Gb memory usage per core as min baseline
    all_memory = [2]
    # Provide 250Mb of additional memory for the system
    system_memory = 0.25
    algs = [get_algorithm_config(x) for x in items]
    for fn in fns:
        for prog in _get_resource_programs(fn, algs):
            resources = config_utils.get_resources(prog, config)
            cores = resources.get("cores", 1)
            memory = _get_prog_memory(resources)
            all_cores.append(cores)
            if memory:
                all_memory.append(memory)
            logger.debug("{prog} requests {cores} cores and {memory}g "
                         "memory for each core.".format(**locals()))

    cores_per_job = max(all_cores)
    if max_multicore:
        cores_per_job = min(cores_per_job, max_multicore)
    if "cores" in sysinfo:
        cores_per_job = min(cores_per_job, int(sysinfo["cores"]))
    memory_per_core = max(all_memory)

    # these callbacks make sure the cores and memory meet minimum requirements
    for fn in fns:
        for ensure in _get_ensure_functions(fn, algs):
            cores_per_job, memory_per_core = ensure(cores_per_job, memory_per_core)

    total = parallel["cores"]
    if total > cores_per_job:
        num_jobs = total // cores_per_job
    else:
        num_jobs, cores_per_job = 1, total
    if cores_per_job == 1:
        memory_per_job = "%.1f" % (memory_per_core + system_memory)
        num_jobs = _scale_jobs_to_memory(num_jobs, memory_per_core, sysinfo)
    else:
        cores_per_job, memory_per_job = _scale_cores_to_memory(cores_per_job,
                                                               memory_per_core, sysinfo,
                                                               system_memory)
    # do not overschedule if we don't have extra items to process
    num_jobs = min(num_jobs, len(items) * multiplier)
    JobResources = collections.namedtuple("JobResources", "num_jobs cores_per_job memory_per_job")
    logger.debug("Configuring %d jobs to run, using %d cores each with %sg of "
                 "memory reserved for each job" % (num_jobs, cores_per_job,
                                                   str(memory_per_job)))
    return JobResources(num_jobs, cores_per_job, str(memory_per_job))

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
            new_arg["config"]["parallel"] = dictdissoc(parallel, "view")
    elif is_std_config_arg(new_arg):
        new_arg["algorithm"]["num_cores"] = int(cores_per_job)
        if parallel:
            new_arg["parallel"] = dictdissoc(parallel, "view")
    else:
        raise ValueError("Unexpected configuration dictionary: %s" % new_arg)
    args = list(args)[:]
    args[new_i] = new_arg
    return args

def _view_from_parallel(parallel, work_dir, config):
    """Translate parallel map into options for a cluster view.
    """
    profile_dir = utils.safe_makedir(os.path.join(work_dir, get_log_dir(config), "ipython"))
    return ipython_cluster.cluster_view(parallel["scheduler"].lower(), parallel["queue"],
                                        parallel["num_jobs"], parallel["cores_per_job"],
                                        profile=profile_dir, start_wait=parallel["timeout"],
                                        extra_params={"resources": parallel["resources"],
                                                      "mem": parallel["mem"]},
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
            jobr = find_job_resources([_get_ipython_fn(x, parallel) for x in fn_names],
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
    algs = [get_algorithm_config(x) for x in items]
    if ((len(algs) > 0 and not algs[0].get("resource_check", True))
          or parallel.get("view") or parallel.get("checkpoint")):
        checkpoint_file = None
    else:
        checkpoint_dir = utils.safe_makedir(os.path.join(work_dir, "checkpoints_ipython"))
        checkpoint_file = _get_checkpoint_file(checkpoint_dir, fn_name)
    fn = _get_ipython_fn(fn_name, parallel)
    if parallel.get("view") is None:
        jobr = find_job_resources([fn], parallel, items, sysinfo, config)
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
            items = [add_cores_to_config(x, parallel["cores_per_job"], parallel) for x in items]
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
