"""Commandline interaction with LSF schedulers.
"""
import re
import subprocess

_jobid_pat = re.compile("Job <(?P<jobid>\d+)> is")

def submit_job(scheduler_args, command):
    """Submit a job to the scheduler, returning the supplied job ID.
    """
    cl = ["bsub"] + scheduler_args + command
    status = subprocess.check_output(cl)
    match = _jobid_pat.search(status)
    return match.groups("jobid")[0]

def stop_job(jobid):
    cl = ["bkill", jobid]
    subprocess.check_call(cl)

def are_running(jobids):
    """Check if all of the submitted job IDs are running.
    """
    run_info = subprocess.check_output(["bjobs"])
    running = []
    for parts in (l.split() for l in run_info.split("\n") if l.strip()):
        if len(parts) >= 3:
            pid, _, status = parts[:3]
            if status.lower() in ["run"]:
                running.append(pid)
    want_running = set(running).intersection(set(jobids))
    return len(want_running) == len(jobids)
