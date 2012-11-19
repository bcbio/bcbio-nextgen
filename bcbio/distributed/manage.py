"""Manage processes on a cluster

Automate:
 - starting working nodes to process the data
 - kicking off an analysis
 - cleaning up nodes on finishing

Currently works on LSF and SGE managed clusters; it's readily generalizable to
other architectures as well.
"""
import time
import math

def run_and_monitor(config, config_file, args, parallel):
    """Run a distributed analysis in s cluster environment, monitoring outputs.
    """
    args = [x for x in args if x is not None]
    cp = config["distributed"]["cluster_platform"]
    cluster = __import__("bcbio.distributed.{0}".format(cp), fromlist=[cp])
    jobids = []
    try:
        print "Starting manager"
        manager_id = start_analysis_manager(cluster, args, config)
        time.sleep(60)
        print "Starting cluster workers"
        jobids.extend(start_workers(cluster, config, config_file, parallel))
        jobids.append(manager_id)
        while not(cluster.are_running(jobids)):
            time.sleep(5)
        print "Running analysis"
        monitor_analysis(cluster, manager_id)
    finally:
        print "Cleaning up cluster workers"
        stop_workers(cluster, jobids)

def start_workers(cluster, config, config_file, parallel):
    """Initiate worker nodes on cluster, returning jobs IDs for management.
    """
    # we can manually specify workers or dynamically get as many as needed
    num_workers = config["distributed"].get("num_workers", None)
    if num_workers in [None, "all"]:
        cores_per_host = config["distributed"].get("cores_per_host", 1)
        if cores_per_host == 0:
            raise ValueError("Set num_workers or cores_per_host in YAML config")
        assert parallel["cores"] is not None, \
               "Supply workers needed if not configured in YAML"
        num_workers = int(math.ceil(float(parallel["cores"]) / cores_per_host))
    program_cl = [config["analysis"]["worker_program"], config_file]
    if parallel.get("task_module", None):
        program_cl.append("--tasks={0}".format(parallel["task_module"]))
    if parallel.get("queues", None):
        program_cl.append("--queues={0}".format(parallel["queues"]))
    args = config["distributed"]["platform_args"].split()
    return [cluster.submit_job(args, program_cl) for _ in range(num_workers)]

def start_analysis_manager(cluster, args, config):
    """Start analysis manager node on cluster.
    """
    cluster_args = config["distributed"]["platform_args"].split()
    program_cl = ["bcbio_nextgen.py"] + args + ["-t", "messaging-main"]
    job_id = cluster.submit_job(cluster_args, program_cl)
    # wait for job to start
    # Avoid this for systems where everything queues as batches
    #while not(cluster.are_running([job_id])):
    #    time.sleep(5)
    return job_id

def monitor_analysis(cluster, job_id):
    """Wait for manager cluster job to finish
    """
    while cluster.are_running([job_id]):
        time.sleep(5)

def stop_workers(cluster, jobids):
    for jobid in jobids:
        try:
            cluster.stop_job(jobid)
        except:
            pass
