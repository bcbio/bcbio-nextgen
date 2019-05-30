#!/usr/bin/env python -Es

import os
import yaml
import math
from argparse import ArgumentParser
from bcbio.install import _get_data_dir
from bcbio import utils
from bcbio.distributed import clargs, resources, ipython as ip
from bcbio.pipeline.main import _pair_samples_with_pipelines


def ipc_fn(parallel):
    cores_per_job = parallel['cores_per_job']
    num_jobs = parallel['num_jobs']
    mincores = int(parallel['resources'][0].split("=")[1])
    if mincores > cores_per_job:
        if cores_per_job > 1:
            mincores = cores_per_job
        else:
            mincores = int(math.ceil(mincores / float(cores_per_job)))
            num_jobs = int(math.ceil(num_jobs / float(mincores)))
    print(f"num of jobs: {num_jobs} with {mincores} engines")


def ipython_fn(parallel, config):
    has_mincores = any(x.startswith("mincores=") for x in parallel["resources"])
    common_cores = min(ip._get_common_cores(config["resources"]), parallel["system_cores"])
    # if has_mincores, check memory usage use min(adj_cores, mincores)
    # min(common_cores, parallel[cores])
    if common_cores > 1 and not has_mincores:
        adj_cores = max(1, int(math.floor(common_cores * float(parallel.get("mem_pct", 1.0)))))
        num_engines = adj_cores
        # if we have less scheduled cores than per machine, use the scheduled count
        if num_engines > parallel["cores"]:
            num_engines = parallel["cores"]
        print(f"common {common_cores} cores_per_job {parallel['cores_per_job']} cores {num_engines}")
        if (parallel['cores_per_job'] * num_engines) > common_cores:
            num_engines = int(math.floor(float(common_cores) / parallel['cores_per_job']))

        parallel["resources"].append("mincores=%s" % num_engines)
    return parallel


def ipython_current(parallel, config):
    has_mincores = any(x.startswith("mincores=") for x in parallel["resources"])
    cores = min(ip._get_common_cores(config["resources"]), parallel["system_cores"])
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
            cores = ip.per_machine_target_cores(cores, parallel["num_jobs"] // cores)
        parallel["resources"].append("mincores=%s" % cores)
    return parallel

if __name__ == "__main__":
    parser = ArgumentParser(description="Merge fastq or bam files")
    parser.add_argument("--yaml-file", required=True,
                        help="yaml file")
    parser.add_argument("--sys-info", required=True,
                        help="system yaml file created by bcbio"
                        " or --sys-info  cores;memory for custom numbers")
    parser.add_argument("--progs", required=True, default=[], action='append', help="look for those tools")
    parser.add_argument("--galaxy", help="custom galaxy file.")
    parser.add_argument("--fixed", help="fixed resources fn", action="store_true")
    parser.add_argument("-n", "--numcores", type=int,
                        default=1, help="Number of concurrent jobs to process.")
    parser.add_argument("-c", "--cores-per-job", type=int,
                        default=1, help="Number of cores to use.")
    parser.add_argument("-m", "--memory-per-job", default=2, help="Memory in GB to reserve per job.")
    parser.add_argument("--timeout", default=15, help="Time to wait before giving up starting.")
    parser.add_argument("--retries", default=0, type=int,
                        help=("Number of retries of failed tasks during "
                              "distributed processing. Default 0 "
                              "(no retries)"))
    parser.add_argument("-s", "--scheduler", help="Type of scheduler to use.",
                        choices=["lsf", "slurm", "torque", "sge", "pbspro"])
    parser.add_argument("-r", "--resources", help="Extra scheduler resource flags.", default=[], action="append")
    parser.add_argument("-q", "--queue", help="Queue to submit jobs to.")
    parser.add_argument("-p", "--tag", help="Tag name to label jobs on the cluster", default="bcb-prep")
    parser.add_argument("-t", "--paralleltype",
                        choices=["local", "ipython"],
                        default="local", help="Run with iptyhon")

    args = parser.parse_args()
    system_config = os.path.join(_get_data_dir(), "galaxy", "bcbio_system.yaml")
    if args.galaxy:
        system_config = args.galaxy
    with open(system_config) as in_handle:
        config = yaml.safe_load(in_handle)

    parallel = clargs.to_parallel(args)
    parallel.update({'progs': args.progs})
    dirs = {'work': os.path.abspath(os.getcwd())}
    if args.sys_info.find(";") > -1:
        info = args.sys_info.split(";")
        sysinfo = {'cores': int(info[0]), 'memory': float(info[1])}
    else:
        if utils.file_exists(args.sys_info):
            sysinfo = yaml.safe_load(open(args.sys_info))[0]
    print(f"system info {sysinfo}")
    samples = []
    pipelines, config = _pair_samples_with_pipelines(args.yaml_file, config)
    for s in pipelines:
        samples = [item for item in pipelines[s]]
    print(f"number of samples {len(samples)}") % len(samples)
    print("after calculate fn")
    parallel = resources.calculate(parallel, samples, sysinfo, config)
    print(parallel)
    if args.fixed:
        print("after fixed ipython fn")
        parallel = ipython_fn(parallel, config)
    else:
        print("after ipython fn")
        parallel = ipython_current(parallel, config)
    print(parallel)
    ipc_fn(parallel)
