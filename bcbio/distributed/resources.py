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
    out = set([])
    for p in progs:
        if p == "aligner":
            for alg in algs:
                aligner = alg.get("aligner")
                if aligner:
                    out.add(aligner)
        elif p == "variantcaller":
            for alg in algs:
                vc = alg.get("variantcaller")
                if vc:
                    if isinstance(vc, (list, tuple)):
                        for x in vc:
                            out.add(x)
                    else:
                        out.add(vc)
        elif p == "gatk-vqsr":
            if config_utils.use_vqsr(algs):
                out.add("gatk-vqsr")
        else:
            out.add(p)
    return sorted(list(out))

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
    # Provide 250Mb of additional memory for the system
    system_memory = 0.25
    algs = [config_utils.get_algorithm_config(x) for x in items]
    progs = _get_resource_programs(parallel.get("progs", []), algs)
    for prog in progs:
        resources = config_utils.get_resources(prog, config)
        cores = resources.get("cores", 1)
        memory = _get_prog_memory(resources)
        all_cores.append(cores)
        if memory:
            all_memory.append(memory)
    # Use modest 1Gb memory usage per core as min baseline if not specified
    if len(all_memory) == 0:
        all_memory.append(1)
    if len(all_cores) == 0:
        all_cores.append(1)
    logger.debug("Resource requests: {progs}; memory: {memory}; cores: {cores}".format(
        progs=", ".join(progs), memory=", ".join("%.1f" % x for x in all_memory),
        cores=", ".join(str(x) for x in all_cores)))

    cores_per_job = max(all_cores)
    if max_multicore:
        cores_per_job = min(cores_per_job, max_multicore)
    if "cores" in sysinfo:
        cores_per_job = min(cores_per_job, int(sysinfo["cores"]))
    memory_per_core = max(all_memory)

    cores_per_job, memory_per_core = _ensure_min_resources(progs, cores_per_job, memory_per_core,
                                                           min_memory=parallel.get("ensure_mem", {}))

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
    logger.debug("Configuring %d jobs to run, using %d cores each with %sg of "
                 "memory reserved for each job" % (num_jobs, cores_per_job,
                                                   str(memory_per_job)))
    parallel = copy.deepcopy(parallel)
    parallel["cores_per_job"] = cores_per_job
    parallel["num_jobs"] = num_jobs
    parallel["mem"] = str(memory_per_job)
    return parallel
