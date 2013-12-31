"""Identify system information for distributed systems, used to manage resources.

This provides a background on cluster and single multicore systems allowing
jobs to be reasonably distributed in cases of higher memory usage.
"""
import copy
import multiprocessing
import os
import resource
import shlex
import socket
import subprocess

import yaml

from bcbio import utils
from bcbio.log import logger

def _get_cache_file(dirs, parallel):
    base_dir = utils.safe_makedir(os.path.join(dirs["work"], "provenance"))
    return os.path.join(base_dir, "system-%s-%s.yaml" % (parallel["type"],
                                                         parallel.get("queue", "default")))

def write_info(dirs, run_parallel, parallel, config):
    """Write cluster or local filesystem resources, spinning up cluster if not present.
    """
    if parallel["type"] in ["ipython"]:
        out_file = _get_cache_file(dirs, parallel)
        if not utils.file_exists(out_file):
            sys_config = copy.deepcopy(config)
            sys_config["algorithm"]["resource_check"] = False
            minfos = _get_machine_info(parallel, run_parallel, sys_config)
            with open(out_file, "w") as out_handle:
                yaml.dump(minfos, out_handle, default_flow_style=False, allow_unicode=False)

def _get_machine_info(parallel, run_parallel, sys_config):
    """Get machine resource information from the job scheduler via either the command-line or the queue.
    """
    if parallel.get("queue") and parallel.get("scheduler"):
        # dictionary as switch statement; can add new scheduler implementation functions as (lowercase) keys
        sched_info_dict = {
                            "slurm": _slurm_info,
                            "torque": _torque_info
                          }
        try:
            return sched_info_dict[parallel["scheduler"].lower()](parallel["queue"])
        except KeyError:
            logger.info("Resource query function not implemented for scheduler \"{0}\"; "
                         "submitting job to queue".format(parallel["scheduler"]))
        except:
            # If something goes wrong, just hit the queue
            logger.warn("Couldn't get machine information from resource query function for queue "
                        "'{0}' on scheduler \"{1}\"; "
                         "submitting job to queue".format(parallel["queue"], parallel["scheduler"]))
    return run_parallel("machine_info", [[sys_config]])

def _slurm_info(queue):
    """Returns machine information for a slurm job scheduler.
    """
    cl = "sinfo -h -p {} --format '%c %m'".format(queue)
    num_cpus, mem = subprocess.check_output(shlex.split(cl)).split()
    # if the queue contains multiple memory configurations, the minimum value is printed with a trailing '+'
    mem = mem.replace('+', '')
    return [{"cores": int(num_cpus), "memory": float(mem) / 1024.0, "name": "slurm_machine"}]

def _torque_info(queue):
    """Return machine information for a torque job scheduler using pbsnodes.
    """
    pbs_out = subprocess.check_output("pbsnodes")
    cores = int(pbs_out[pbs_out.find("np = ")+5])
    name = pbs_out.split("\n")[0]
    mem = min([float(string.split("=")[1].rstrip("kb")) / 1048576.0
               for string in pbs_out.split(",") if "availmem" in string])
    return [{"cores": cores, "memory": mem, "name": name}]

def _combine_machine_info(xs):
    if len(xs) == 1:
        return xs[0]
    else:
        raise NotImplementedError("Add logic to pick specification from non-homogeneous clusters.")

def get_info(dirs, parallel):
    """Retrieve cluster or local filesystem resources from pre-retrieved information.
    """
    if parallel["type"] in ["ipython"]:
        cache_file = _get_cache_file(dirs, parallel)
        if utils.file_exists(cache_file):
            with open(cache_file) as in_handle:
                minfo = yaml.load(in_handle)
            return _combine_machine_info(minfo)
        else:
            return {}
    else:
        return _combine_machine_info(machine_info())

def machine_info():
    """Retrieve core and memory information for the current machine.
    """
    import psutil
    BYTES_IN_GIG = 1073741824
    free_bytes = psutil.virtual_memory().available
    return [{"memory": float(free_bytes / BYTES_IN_GIG), "cores": multiprocessing.cpu_count(),
             "name": socket.gethostname()}]

def open_file_limit():
    return resource.getrlimit(resource.RLIMIT_NOFILE)[0]
