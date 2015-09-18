"""Estimate resources required for processing a set of tasks.

Uses annotations provided in multitasks.py for each function to identify utilized
programs, then extracts resource requirements from the input bcbio_system file.
"""
import copy
import math

from bcbio.pipeline import config_utils
from bcbio.log import logger

def _get_resource_programs(progs, algs):
    """Retrieve programs used in analysis based on algorithm configurations.
    Handles special cases like aligners and variant callers.
    """
    checks = {"gatk-vqsr": config_utils.use_vqsr,
              "snpeff": config_utils.use_snpeff,
              "bcbio-variation-recall": config_utils.use_bcbio_variation_recall}
    parent_child = {"vardict": _parent_prefix("vardict")}
    out = set([])
    for p in progs:
        if p == "aligner":
            for alg in algs:
                aligner = alg.get("aligner")
                if aligner:
                    out.add(aligner)
        elif p in ["variantcaller", "svcaller"]:
            if p == "variantcaller": 
                for key, fn in parent_child.items():
                    if fn(algs):
                        out.add(key)
            for alg in algs:
                callers = alg.get(p)
                if callers:
                    if isinstance(callers, (list, tuple)):
                        for x in callers:
                            out.add(x)
                    else:
                        out.add(callers)
        elif p in checks:
            if checks[p](algs):
                out.add(p)
        else:
            out.add(p)
    return sorted(list(out))

def _parent_prefix(prefix):
    """Identify a parent prefix we should add to resources if present in a caller name.
    """
    def run(algs):
        for alg in algs:
            vcs = alg.get("variantcaller")
            if vcs:
                if not isinstance(vcs, (list, tuple)):
                    vcs = [vcs]
                return any(vc.startswith(prefix) for vc in vcs)
    return run

def _ensure_min_resources(progs, cores, memory, min_memory):
    """Ensure setting match minimum resources required for used programs.
    """
    for p in progs:
        if p in min_memory:
            if not memory or cores * memory < min_memory[p]:
                memory = float(min_memory[p]) / cores
    return cores, memory

def _str_memory_to_gb(memory):
    val = float(memory[:-1])
    units = memory[-1]
    if units.lower() == "m":
        val = val / 1000.0
    else:
        assert units.lower() == "g", "Unexpected memory units: %s" % memory
    return val

def _get_prog_memory(resources, cores_per_job):
    """Get expected memory usage, in Gb per core, for a program from resource specification.
    """
    out = None
    for jvm_opt in resources.get("jvm_opts", []):
        if jvm_opt.startswith("-Xmx"):
            out = _str_memory_to_gb(jvm_opt[4:])
    memory = resources.get("memory")
    if memory:
        out = _str_memory_to_gb(memory)
    prog_cores = resources.get("cores")
    # if a single core with memory is requested for the job
    # and we run multiple cores, scale down to avoid overscheduling
    if out and prog_cores and int(prog_cores) == 1 and cores_per_job > int(prog_cores):
        out = out / float(cores_per_job)
    return out

def _scale_cores_to_memory(cores, mem_per_core, sysinfo, system_memory):
    """Scale multicore usage to avoid excessive memory usage based on system information.
    """
    total_mem = "%.2f" % (cores * mem_per_core + system_memory)
    if "cores" not in sysinfo:
        return cores, total_mem, 1.0

    total_mem = min(float(total_mem), float(sysinfo["memory"]) - system_memory)
    cores = min(cores, int(sysinfo["cores"]))
    mem_cores = int(math.floor(float(total_mem) / mem_per_core))  # cores based on available memory
    if mem_cores < 1:
        out_cores = 1
    elif mem_cores < cores:
        out_cores = mem_cores
    else:
        out_cores = cores
    mem_pct = float(out_cores) / float(cores)
    return out_cores, total_mem, mem_pct

def _scale_jobs_to_memory(jobs, mem_per_core, sysinfo):
    """When scheduling jobs with single cores, avoid overscheduling due to memory.
    """
    if "cores" not in sysinfo:
        return jobs, 1.0
    sys_mem_per_core = float(sysinfo["memory"]) / float(sysinfo["cores"])
    if sys_mem_per_core < mem_per_core:
        pct = sys_mem_per_core / float(mem_per_core)
        target_jobs = int(math.floor(jobs * pct))
        return max(target_jobs, 1), pct
    else:
        return jobs, 1.0

def calculate(parallel, items, sysinfo, config, multiplier=1,
              max_multicore=None):
    """Determine cores and workers to use for this stage based on used programs.
    multiplier specifies the number of regions items will be split into during
    processing.
    max_multicore specifies an optional limit on the maximum cores. Can use to
    force single core processing during specific tasks.
    sysinfo specifies cores and memory on processing nodes, allowing us to tailor
    jobs for available resources.
    """
    assert len(items) > 0, "Finding job resources but no items to process"
    all_cores = []
    all_memory = []
    # Provide 100Mb of additional memory for the system
    system_memory = 0.10
    algs = [config_utils.get_algorithm_config(x) for x in items]
    progs = _get_resource_programs(parallel.get("progs", []), algs)
    # Calculate cores
    for prog in progs:
        resources = config_utils.get_resources(prog, config)
        all_cores.append(resources.get("cores", 1))
    if len(all_cores) == 0:
        all_cores.append(1)
    cores_per_job = max(all_cores)
    if max_multicore:
        cores_per_job = min(cores_per_job, max_multicore)
    if "cores" in sysinfo:
        cores_per_job = min(cores_per_job, int(sysinfo["cores"]))
    total = parallel["cores"]
    if total > cores_per_job:
        num_jobs = total // cores_per_job
    else:
        num_jobs, cores_per_job = 1, total

    # Calculate memory. Use 1Gb memory usage per core as min baseline if not specified
    for prog in progs:
        resources = config_utils.get_resources(prog, config)
        memory = _get_prog_memory(resources, cores_per_job)
        if memory:
            all_memory.append(memory)
    if len(all_memory) == 0:
        all_memory.append(1)
    memory_per_core = max(all_memory)

    logger.debug("Resource requests: {progs}; memory: {memory}; cores: {cores}".format(
        progs=", ".join(progs), memory=", ".join("%.2f" % x for x in all_memory),
        cores=", ".join(str(x) for x in all_cores)))

    cores_per_job, memory_per_core = _ensure_min_resources(progs, cores_per_job, memory_per_core,
                                                           min_memory=parallel.get("ensure_mem", {}))
    if cores_per_job == 1:
        memory_per_job = "%.2f" % memory_per_core
        num_jobs, mem_pct = _scale_jobs_to_memory(num_jobs, memory_per_core, sysinfo)
    else:
        cores_per_job, memory_per_job, mem_pct = _scale_cores_to_memory(cores_per_job,
                                                                        memory_per_core, sysinfo,
                                                                        system_memory)
        # For local runs with multiple jobs and multiple cores, potentially scale jobs down
        if num_jobs > 1 and parallel.get("type") == "local":
            memory_per_core = float(memory_per_job) / cores_per_job
            num_jobs, _ = _scale_jobs_to_memory(num_jobs, memory_per_core, sysinfo)

    # do not overschedule if we don't have extra items to process
    num_jobs = min(num_jobs, len(items) * multiplier)
    logger.debug("Configuring %d jobs to run, using %d cores each with %sg of "
                 "memory reserved for each job" % (num_jobs, cores_per_job,
                                                   str(memory_per_job)))
    parallel = copy.deepcopy(parallel)
    parallel["cores_per_job"] = cores_per_job
    parallel["num_jobs"] = num_jobs
    parallel["mem"] = str(memory_per_job)
    parallel["mem_pct"] = "%.2f" % mem_pct
    parallel["system_cores"] = sysinfo.get("cores", 1)
    return parallel
