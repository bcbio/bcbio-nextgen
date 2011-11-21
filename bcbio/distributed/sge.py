"""Commandline interaction with SGE cluster schedulers.
"""
import re
import time
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
    # handle SGE errors, retrying to get the current status
    max_retries = 10
    tried = 0
    while 1:
        try:
            run_info = subprocess.check_output(["qstat"])
            break
        except:
            tried += 1
            if tried > max_retries:
                raise
            time.sleep(5)
    running = []
    for parts in (l.split() for l in run_info.split("\n") if l.strip()):
        if len(parts) >= 5:
            pid, _, _, _, status = parts[:5]
            if status.lower() in ["r"]:
                running.append(pid)
    want_running = set(running).intersection(set(jobids))
    return len(want_running) == len(jobids)

def available_nodes(scheduler_args):
    """Retrieve a count of available nodes in the configured queue.
    """
    cl = ["qstat", "-f"]
    info = subprocess.check_output(cl)
    total = 0
    for i, line in enumerate(info.split("\n")):
        if i > 1 and not line.startswith("----") and line.startswith(tuple(scheduler_args)):
            _, _, counts = line.split()[:3]
            _, _, avail = counts.split("/")
            total += int(avail)
    return total if total > 0 else None
