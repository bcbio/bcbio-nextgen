#!/usr/bin/env python -Es

import os
import yaml
from argparse import ArgumentParser
from bcbio.install import _get_data_dir
from bcbio import utils
from bcbio.distributed import clargs, resources
from bcbio.pipeline.main import _pair_samples_with_pipelines


if __name__ == "__main__":
    parser = ArgumentParser(description="Merge fastq or bam files")
    parser.add_argument("--yaml-file", required=True,
                        help="yaml file")
    parser.add_argument("--sys-info", required=True,
                        help="system yaml file created by bcbio"
                        " or --sys-info  cores;memory for custom numbers")
    parser.add_argument("--progs", required=True, default=[], action='append', help="look for those tools")
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
    with open(system_config) as in_handle:
        config = yaml.load(in_handle)

    parallel = clargs.to_parallel(args)
    parallel.update({'progs': args.progs})
    dirs = {'work': os.path.abspath(os.getcwd())}
    if args.sys_info.find(";") > -1:
        info = args.sys_info.split(";")
        sysinfo = {'cores': int(info[0]), 'memory': float(info[1])}
    else:
        if utils.file_exists(args.sys_info):
            sysinfo = yaml.load(open(args.sys_info))[0]
    print "system info %s" % sysinfo
    samples = []
    pipelines, config = _pair_samples_with_pipelines(args.yaml_file, config)
    for s in pipelines:
        samples = [item for item in pipelines[s]]
    print "number of samples %s" % len(samples)
    parallel = resources.calculate(parallel, samples, sysinfo, config)
    print parallel
