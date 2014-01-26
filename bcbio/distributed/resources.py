"""Estimate resources required for processing a set of tasks.

Uses annotations provided in multitasks.py for each function to identify utilized
programs, then extracts resource requirements from the input bcbio_system file.
"""
import collections
import math

from bcbio.pipeline import config_utils
from bcbio.log import logger

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
                      "snpEff", "cufflinks", "picard", "rnaseqc"])
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

def calculate(fns, parallel, items, sysinfo, config, multiplier=1,
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
    all_memory = []
    # Provide 250Mb of additional memory for the system
    system_memory = 0.25
    algs = [config_utils.get_algorithm_config(x) for x in items]
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

    # Use modest 1Gb memory usage per core as min baseline if not specified
    if len(all_memory) == 0:
        all_memory.append(1)

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
