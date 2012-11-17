"""Distributed execution using an IPython cluster.

Uses IPython parallel to setup a cluster and manage execution:

http://ipython.org/ipython-doc/stable/parallel/index.html

Borrowed from Rory Kirchner's Bipy cluster implementation:

https://github.com/roryk/bipy/blob/master/bipy/cluster/__init__.py
"""
import time
import subprocess

from IPython.parallel import Client

def _start(workers_needed, profile, delay):
    """Starts cluster from commandline.

    XXX: in the future, add "--cluster-id=" + self._cluster_id to
    this, to run each new cluster with a different ID, so we can
    reuse the same profile. right now there is a bug in ipython that
    doesn't support this"""
    subprocess.check_call(["ipcluster", "start",
                           "--daemonize=True",
                           "--delay=%s" % delay, 
                           "--log-level=%s" % 30,
                           "--n=%s" % workers_needed,
                           "--profile=%s" % profile])

def _stop(profile):
    # add carg = "--cluster-id=%s" % (self._cluster_id) when
    # this gets fixed in iPython
    subprocess.check_call(["ipcluster", "stop", "--profile=%s" % profile])

def _is_up(client, n):
    try:
        up = len(client.ids)
    except IOError:
        return False
    else:
        not_up = n - up
        if not_up > 0:
            return False
        else:
            return True

def run_and_monitor(config, config_file, args, workers_needed=None):
    """Run a distributed analysis after starting an Ipython parallel environment.
    """
    delay = 10
    max_delay = 300
    profile = "default"
    _start(workers_needed, profile, delay)
    client = Client(profile=profile)
    try:
        slept = 0
        while not _is_up(client, workers_needed):
            time.sleep(delay)
            slept += delay
            if slept > max_delay:
                raise IOError("Cluster startup timed out.")
        cluster_view = client.load_balanced_view()
    finally:
        _stop(profile)
