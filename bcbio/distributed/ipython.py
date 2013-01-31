"""Distributed execution using an IPython cluster.

Uses IPython parallel to setup a cluster and manage execution:

http://ipython.org/ipython-doc/stable/parallel/index.html

Borrowed from Rory Kirchner's Bipy cluster implementation:

https://github.com/roryk/bipy/blob/master/bipy/cluster/__init__.py
"""
import os
import copy
import pipes
import time
import uuid
import subprocess
import contextlib

from bcbio import utils
from bcbio.log import setup_logging, logger

from IPython.parallel import Client
from IPython.parallel.apps import launcher
from IPython.utils import traitlets

# ## Custom launchers

class BcbioLSFEngineSetLauncher(launcher.LSFEngineSetLauncher):
    """Custom launcher handling heterogeneous clusters on LSF.
    """
    cores = traitlets.Integer()
    default_template = traitlets.Unicode("""#!/bin/sh
    #BSUB -q {queue}
    #BSUB -J bcbio-ipengine[1-{n}]
    #BSUB -oo bcbio-ipengine.bsub.%%J
    #BSUB -n {cores}
    #BSUB -R "span[hosts=1]"
    %s --profile-dir="{profile_dir}" --cluster-id="{cluster_id}"
    """%(' '.join(map(pipes.quote, launcher.ipengine_cmd_argv))))
    def start(self, n):
        return super(BcbioLSFEngineSetLauncher, self).start(n)

class BcbioSGEEngineSetLauncher(launcher.SGEEngineSetLauncher):
    """Custom launcher handling heterogeneous clusters on SGE.
    """
    cores = traitlets.Integer()
    default_template = traitlets.Unicode("""#$ -V
#$ -cwd
#$ -b y
#$ -j y
#$ -S /bin/sh
#$ -q {queue}
#$ -N bcbio-ipengine
#$ -t 1-{n}
#$ -pe threaded {cores}
%s --profile-dir="{profile_dir}" --cluster-id="{cluster_id}"
"""%(' '.join(map(pipes.quote, launcher.ipengine_cmd_argv))))

    def start(self, n):
        return super(BcbioSGEEngineSetLauncher, self).start(n)

# ## Control clusters

def _start(parallel, profile, cluster_id, delay):
    """Starts cluster from commandline.
    """
    scheduler = parallel["scheduler"].upper()
    engine_class = "bcbio.distributed.ipython.Bcbio%sEngineSetLauncher" % scheduler
    subprocess.check_call(
        ["ipcluster", "start",
         "--daemonize=True",
         "--delay=%s" % delay,
         "--log-level=%s" % "WARN",
         "--profile=%s" % profile,
         #"--cluster-id=%s" % cluster_id,
         "--n=%s" % parallel["num_jobs"],
         "--%s.cores=%s" % (engine_class, parallel["cores_per_job"]),
         "--HubFactory.ip=*",
         "--IPClusterStart.controller_launcher_class=%s" % scheduler,
         "--IPClusterStart.engine_launcher_class=%s" % engine_class,
         "--%sLauncher.queue=%s" % (scheduler, parallel["queue"]),
         ])

def _stop(profile, cluster_id):
    subprocess.check_call(["ipcluster", "stop", "--profile=%s" % profile,
                           #"--cluster-id=%s" % cluster_id
                           ])

def _is_up(profile, cluster_id, n):
    try:
        #client = Client(profile=profile, cluster_id=cluster_id)
        client = Client(profile=profile)
        up = len(client.ids)
    except IOError, msg:
        return False
    else:
        return up >= n

@contextlib.contextmanager
def cluster_view(parallel, config):
    """Provide a view on an ipython cluster for processing.

    parallel is a dictionary with:
      - scheduler: The type of cluster to start (lsf, sge).
      - num_jobs: Number of jobs to start.
      - cores_per_job: The number of cores to use for each job.
    """
    delay = 5
    max_delay = 300
    max_tries = 10
    profile = "bcbio_nextgen"
    cluster_id = str(uuid.uuid1())
    num_tries = 0
    while 1:
        try:
            _start(parallel, profile, cluster_id, delay)
            break
        except subprocess.CalledProcessError:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(delay)
    try:
        slept = 0
        while not _is_up(profile, cluster_id, parallel["num_jobs"]):
            time.sleep(delay)
            slept += delay
            if slept > max_delay:
                raise IOError("Cluster startup timed out.")
        #client = Client(profile=profile, cluster_id=cluster_id)
        client = Client(profile=profile)
        # push config to all engines and force them to set up logging
        client[:]['config'] = config
        client[:].execute('from bcbio.log import setup_logging')
        client[:].execute('setup_logging(config)')
        client[:].execute('from bcbio.log import logger')
        yield client.load_balanced_view()
    finally:
        _stop(profile, cluster_id)

def dictadd(orig, k, v):
    """Imitates immutability by adding a key/value to a new dictionary.
    Works around not being able to deepcopy view objects; can remove this
    once we create views on demand.
    """
    view = orig.pop("view", None)
    new = copy.deepcopy(orig)
    new[k] = v
    if view:
        orig["view"] = view
        new["view"] = view
    return new

def _find_cores_per_job(fn, parallel, config):
    """Determine cores and workers to use for this stage based on function metadata.

    TODO: Generalize. Currently handles single core jobs.
    """
    total = parallel["cores"]
    return total, 1

def runner(parallel, fn_name, items, work_dir, config):
    """Run a task on an ipython parallel cluster, allowing alternative queue types.

    This will spawn clusters for parallel and custom queue types like multicore
    and high I/O tasks on demand.

    A checkpoint directory keeps track of finished tasks, avoiding spinning up clusters
    for sections that have been previous processed.
    """
    setup_logging(config)
    out = []
    checkpoint_dir = utils.safe_makedir(os.path.join(work_dir, "checkpoints_ipython"))
    checkpoint_file = os.path.join(checkpoint_dir, "%s.done" % fn_name)
    fn = getattr(__import__("{base}.ipythontasks".format(base=parallel["module"]),
                            fromlist=["ipythontasks"]),
                 fn_name)
    num_jobs, cores_per_job = _find_cores_per_job(fn, parallel, config)
    parallel = dictadd(parallel, "cores_per_job", cores_per_job)
    parallel = dictadd(parallel, "num_jobs", num_jobs)
    # already finished, run locally on current machine to collect details
    if os.path.exists(checkpoint_file):
        logger.info("ipython: %s -- local; checkpoint passed" % fn_name)
        for args in items:
            if args:
                data = fn(args)
                if data:
                    out.extend(data)
    # Run on a standard parallel queue
    else:
        logger.info("ipython: %s" % fn_name)
        with cluster_view(parallel, config) as view:
            xs = [x for x in items if x is not None]
            if len(xs) > 0:
                for data in view.map_sync(fn, xs):
                    if data:
                        out.extend(data)
    with open(checkpoint_file, "w") as out_handle:
        out_handle.write("done\n")
    return out
