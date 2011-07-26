#!/usr/bin/env python
"""Run an automated analysis pipeline in a distributed cluster architecture.

Automate:
 - starting working nodes to process the data
 - kicking off an analysis
 - cleaning up nodes on finishing

Currently works on LSF managed clusters but is written generally to work
on other architectures.

Usage:
  run_distributed_job.py <config_file> <fc_dir> [<run_info_yaml>]
"""
import sys
import time
import subprocess

import yaml

from bcbio.distributed import lsf

def main(config_file, fc_dir, run_info_yaml=None):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    assert config["algorithm"]["num_cores"] == "messaging", \
           "Designed for use with messaging parallelization"
    cluster = globals()[config["distributed"]["cluster_platform"]]
    print "Starting cluster workers"
    jobids = start_workers(cluster, config, config_file)
    try:
        print "Running analysis"
        run_analysis(config_file, fc_dir, run_info_yaml, cluster, config)
    finally:
        print "Cleaning up cluster workers"
        stop_workers(cluster, jobids)

def start_workers(cluster, config, config_file):
    args = config["distributed"]["platform_args"].split()
    program_cl = [config["analysis"]["worker_program"], config_file]
    jobids = [cluster.submit_job(args, program_cl)
              for _ in range(config["distributed"]["num_workers"])]
    while not(cluster.are_running(jobids)):
        time.sleep(5)
    return jobids

def run_analysis(config_file, fc_dir, run_info_yaml, cluster, config):
    # get arguments, ignoring need for multiple cores in arguments
    # since analysis is single threaded runner
    args = []
    ignore = False
    for a in config["distributed"]["platform_args"].split():
        if ignore:
            ignore = False
        elif a == "-n":
            ignore = True
        else:
            args.append(a)
    program_cl = [config["analysis"]["process_program"], config_file, fc_dir]
    if run_info_yaml:
        program_cl.append(run_info_yaml)
    jobid = cluster.submit_job(args, program_cl)
    try:
        _monitor_analysis(cluster, jobid)
    except:
        stop_workers(cluster, [jobid])
        raise

def _monitor_analysis(cluster, jobid):
    # wait for job to start
    while not(cluster.are_running([jobid])):
        time.sleep(5)
    # wait for job to finish
    while cluster.are_running([jobid]):
        time.sleep(5)

def stop_workers(cluster, jobids):
    for jobid in jobids:
        cluster.stop_job(jobid)

if __name__ == "__main__":
    main(*sys.argv[1:])
