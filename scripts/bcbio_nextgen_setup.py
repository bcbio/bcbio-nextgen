#!/usr/bin/env python
"""Automate setup of bcbio_nextgen pipeline for common configurations.

Handles:
  Setup of ipython queues for parallel and multicore jobs.

  TODO: IO intensive jobs
"""
import os
import sys
import shutil
import pipes
import argparse

import sh

from IPython.parallel.apps import launcher

def main(args):
    print "Setting up profiles for ipython: %s" % args.scheduler
    setup_ipython(args.scheduler, _prep_queues(args.queues))

def setup_ipython(scheduler, queues):
    for ptype, queue in zip(["", "multicore"], queues):
        if ptype:
            profile = "%s_%s" % (scheduler, ptype)
        else:
            profile = scheduler
        sh.ipython("profile", "create", reset=True, parallel=True, profile=profile)
        config_dir = _find_ipython_dir(profile)
        update_ipcontroller_config(config_dir)
        update_ipcluster_config(config_dir, scheduler, ptype, queue)

def _prep_queues(queues):
    if len(queues) == 1:
        return queues * 3
    else:
        assert len(queues) >= 2, queues
        return queues

def _find_ipython_dir(profile):
    """Find ipython directory associated with the given profile.
    """
    to_test = [os.path.join(os.environ["HOME"], ".ipython"),
               os.path.join(os.environ["HOME"], ".config", "ipython")]
    for basedir in to_test:
        profile_dir = os.path.join(basedir, "profile_%s" % profile)
        if os.path.exists(profile_dir) and os.path.isdir(profile_dir):
            return profile_dir
    raise ValueError("ipython configuration directory for %s not found" % profile)

def _update_config_file(new_config, profile_dir, fname):
    config_file = os.path.join(profile_dir, fname) 
    orig_config_file = config_file + ".orig"
    shutil.move(config_file, orig_config_file)
    with open(orig_config_file) as in_handle:
        with open(config_file, "w") as out_handle:
            for line in in_handle:
                out_handle.write(line)
                if line.startswith("c = get_config"):
                    out_handle.write(new_config + "\n")
    os.remove(orig_config_file)

# ## ipcontroller

def _ipcontroller_config():
    return "\n".join(
        ["# Added by bcbio_nextgen",
         "c.HubFactory.ip = '*'"])

def update_ipcontroller_config(profile_dir):
    _update_config_file(_ipcontroller_config(),
                        profile_dir, "ipcontroller_config.py")

# ## ipcluster

def get_ipcluster_config(scheduler, parallel_type, queue):
    batch_configs = {"lsf": lsf_batch_config,
                     "sge": sge_batch_config}
    return "\n".join(
        ["# Added by bcbio_nextgen",
         "c.IPClusterStart.controller_launcher_class = '%s'" % scheduler.upper(),
         "c.IPClusterStart.engine_launcher_class = '%s'" % scheduler.upper()] +
        batch_configs[scheduler](parallel_type, queue))

def update_ipcluster_config(profile_dir, scheduler, parallel_type, queue):
    _update_config_file(get_ipcluster_config(scheduler, parallel_type, queue),
                        profile_dir, "ipcluster_config.py")

# ## LSF/bsub

def _get_std_cmd(cmd):
    launcher_cmds = {"ipengine": launcher.ipengine_cmd_argv,
                     "ipcontroller": launcher.ipcontroller_cmd_argv}
    return '%s --log-to-file --profile-dir="{profile_dir}" --cluster-id="{cluster_id}"' % \
        (" ".join(map(pipes.quote, launcher_cmds[cmd])))

LSF_BATCHES = {
"": """#!/bin/sh
#BSUB -q {queue}
#BSUB -J %s[1-{n}]
#BSUB -oo %s.bsub.%%J
%s
""",
"multicore": """#!/bin/sh
#BSUB -q {queue}
#BSUB -J %s
#BSUB -oo %s.bsub.%%J
#BSUB -n {n}
#BSUB -R "span[hosts=1]"
%s
"""}

def lsf_batch_config(parallel_type, queue):
    out = []
    for ltype, launcher in [("ipengine", "LSFEngineSetLauncher"),
                            ("ipcontroller", "LSFControllerLauncher")]:
        batch_str = LSF_BATCHES[parallel_type] % (ltype, ltype, _get_std_cmd(ltype))
        final_str = 'c.%s.batch_template = """%s\n"""' % (launcher, batch_str.rstrip())
        out.extend(final_str.split("\n"))
    out.append("c.LSFLauncher.queue = '%s'" % queue)
    return out
        
# ## SGE/qsub

SGE_BATCHES = {
"" : """#$ -V
#$ -cwd
#$ -b y
#$ -j y
#$ -S /bin/sh
#$ -q {queue}
#$ -N %s
#$ -t 1-{n}
%s
""",
"multicore": """#$ -V
#$ -cwd
#$ -b y
#$ -j y
#$ -S /bin/sh
#$ -q {queue}
#$ -N %s
#$ -pe threaded {n}
%s
"""}

def sge_batch_config(parallel_type, queue):
    out = []
    for ltype, launcher in [("ipengine", "SGEEngineSetLauncher"),
                            ("ipcontroller", "SGEControllerLauncher")]:
        batch_str = SGE_BATCHES[parallel_type] % (ltype, _get_std_cmd(ltype))
        final_str = 'c.%s.batch_template = """%s\n"""' % (launcher, batch_str.rstrip())
        out.extend(final_str.split("\n"))
    out.append("c.SGELauncher.queue = '%s'" % queue)
    return out

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--scheduler", type=lambda x: x.lower(), default="lsf",
                        help="Type of cluster scheduler: lsf or sge")
    parser.add_argument("-q", "--queues", nargs="+", required=True,
                        help="Queues to place jobs on. With a single queue "\
                            "will use that queue for all types. With multiple "\
                            "queues, assign the first to parallel, second to "\
                            "multicore jobs, and third to IO intensive jobs.")
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        main(parser.parse_args())
