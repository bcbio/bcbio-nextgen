"""Commandline interaction with SGE cluster schedulers.
"""
import re
import subprocess

_jobid_pat = re.compile('Your job (?P<jobid>\d+) \("')

def submit_job(scheduler_args, command):
    """Submit a job to the scheduler, returning the supplied job ID.
    """
    cl = ["qsub", "-cwd", "-b", "y", "-j", "y"] + scheduler_args + command
    status = subprocess.check_output(cl)
    match = _jobid_pat.search(status)
    return match.groups("jobid")[0]

def stop_job(jobid):
    cl = ["qdel", jobid]
    subprocess.check_call(cl)

def are_running(jobids):
    """Check if submitted job IDs are running.
    """
    run_info = subprocess.check_output(["qstat"])
    running = []
    for parts in (l.split() for l in run_info.split("\n") if l.strip()):
        if len(parts) >= 5:
            pid, _, _, _, status = parts[:5]
            if status.lower() in ["r"]:
                running.append(pid)
    want_running = set(running).intersection(set(jobids))
    return len(want_running) == len(jobids)
