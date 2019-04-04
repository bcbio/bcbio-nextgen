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
from xml.etree import ElementTree as ET

from bcbio import utils
from bcbio.log import logger

def _get_cache_file(dirs, parallel):
    base_dir = utils.safe_makedir(os.path.join(dirs["work"], "provenance"))
    return os.path.join(base_dir, "system-%s-%s.yaml" % (parallel["type"],
                                                         parallel.get("queue", "default")))

def write_info(dirs, parallel, config):
    """Write cluster or local filesystem resources, spinning up cluster if not present.
    """
    if parallel["type"] in ["ipython"] and not parallel.get("run_local"):
        out_file = _get_cache_file(dirs, parallel)
        if not utils.file_exists(out_file):
            sys_config = copy.deepcopy(config)
            minfos = _get_machine_info(parallel, sys_config, dirs, config)
            with open(out_file, "w") as out_handle:
                yaml.safe_dump(minfos, out_handle, default_flow_style=False, allow_unicode=False)

def _get_machine_info(parallel, sys_config, dirs, config):
    """Get machine resource information from the job scheduler via either the command line or the queue.
    """
    if parallel.get("queue") and parallel.get("scheduler"):
        # dictionary as switch statement; can add new scheduler implementation functions as (lowercase) keys
        sched_info_dict = {
                            "slurm": _slurm_info,
                            "torque": _torque_info,
                            "sge": _sge_info
                          }
        if parallel["scheduler"].lower() in sched_info_dict:
            try:
                return sched_info_dict[parallel["scheduler"].lower()](parallel.get("queue", ""))
            except:
                # If something goes wrong, just hit the queue
                logger.exception("Couldn't get machine information from resource query function for queue "
                                 "'{0}' on scheduler \"{1}\"; "
                                 "submitting job to queue".format(parallel.get("queue", ""), parallel["scheduler"]))
        else:
            logger.info("Resource query function not implemented for scheduler \"{0}\"; "
                         "submitting job to queue".format(parallel["scheduler"]))
    from bcbio.distributed import prun
    with prun.start(parallel, [[sys_config]], config, dirs) as run_parallel:
        return run_parallel("machine_info", [[sys_config]])

def _slurm_info(queue):
    """Returns machine information for a slurm job scheduler.
    """
    cl = "sinfo -h -p {} --format '%c %m %D'".format(queue)
    num_cpus, mem, num_nodes = subprocess.check_output(shlex.split(cl)).decode().split()
    # if the queue contains multiple memory configurations, the minimum value is printed with a trailing '+'
    mem = float(mem.replace('+', ''))
    num_cpus = int(num_cpus.replace('+', ''))
    # handle small clusters where we need to allocate memory for bcbio and the controller
    # This will typically be on cloud AWS machines
    bcbio_mem = 2000
    controller_mem = 4000
    if int(num_nodes) < 3 and mem > (bcbio_mem + controller_mem) * 2:
        mem = mem - bcbio_mem - controller_mem
    return [{"cores": int(num_cpus), "memory": mem / 1024.0, "name": "slurm_machine"}]

def _torque_info(queue):
    """Return machine information for a torque job scheduler using pbsnodes.

    To identify which host to use it tries to parse available hosts
    from qstat -Qf `acl_hosts`. If found, it uses these and gets the
    first node from pbsnodes matching to the list. If no attached
    hosts are available, it uses the first host found from pbsnodes.
    """
    nodes = _torque_queue_nodes(queue)
    pbs_out = subprocess.check_output(["pbsnodes"]).decode()
    info = {}
    for i, line in enumerate(pbs_out.split("\n")):
        if i == 0 and len(nodes) == 0:
            info["name"] = line.strip()
        elif line.startswith(nodes):
            info["name"] = line.strip()
        elif info.get("name"):
            if line.strip().startswith("np = "):
                info["cores"] = int(line.replace("np = ", "").strip())
            elif line.strip().startswith("status = "):
                mem = [x for x in pbs_out.split(",") if x.startswith("physmem=")][0]
                info["memory"] = float(mem.split("=")[1].rstrip("kb")) / 1048576.0
                return [info]

def _torque_queue_nodes(queue):
    """Retrieve the nodes available for a queue.

    Parses out nodes from `acl_hosts` in qstat -Qf and extracts the
    initial names of nodes used in pbsnodes.
    """
    qstat_out = subprocess.check_output(["qstat", "-Qf", queue]).decode()
    hosts = []
    in_hosts = False
    for line in qstat_out.split("\n"):
        if line.strip().startswith("acl_hosts = "):
            hosts.extend(line.replace("acl_hosts = ", "").strip().split(","))
            in_hosts = True
        elif in_hosts:
            if line.find(" = ") > 0:
                break
            else:
                hosts.extend(line.strip().split(","))
    return tuple([h.split(".")[0].strip() for h in hosts if h.strip()])

def median_left(x):
    if len(x) < 1:
        return None
    sortedval = sorted(x)
    centre = int((len(sortedval)-1)/2)
    return sortedval[centre]

def _sge_info(queue):
    """Returns machine information for an sge job scheduler.
    """
    qhost_out = subprocess.check_output(["qhost", "-q", "-xml"]).decode()
    qstat_queue = ["-q", queue] if queue and "," not in queue else []
    qstat_out = subprocess.check_output(["qstat", "-f", "-xml"] + qstat_queue).decode()
    slot_info = _sge_get_slots(qstat_out)
    mem_info = _sge_get_mem(qhost_out, queue)
    machine_keys = slot_info.keys()
    #num_cpus_vec = [slot_info[x]["slots_total"] for x in machine_keys]
    #mem_vec = [mem_info[x]["mem_total"] for x in machine_keys]
    mem_per_slot = [mem_info[x]["mem_total"] / float(slot_info[x]["slots_total"]) for x in machine_keys]
    min_ratio_index = mem_per_slot.index(median_left(mem_per_slot))
    mem_info[machine_keys[min_ratio_index]]["mem_total"]
    return [{"cores": slot_info[machine_keys[min_ratio_index]]["slots_total"],
             "memory": mem_info[machine_keys[min_ratio_index]]["mem_total"],
             "name": "sge_machine"}]

def _sge_get_slots(xmlstring):
    """ Get slot information from qstat
    """
    rootxml = ET.fromstring(xmlstring)
    my_machine_dict = {}
    for queue_list in rootxml.iter("Queue-List"):
        # find all hosts supporting queues
        my_hostname = queue_list.find("name").text.rsplit("@")[-1]
        my_slots = queue_list.find("slots_total").text
        my_machine_dict[my_hostname] = {}
        my_machine_dict[my_hostname]["slots_total"] = int(my_slots)
    return my_machine_dict

def _sge_get_mem(xmlstring, queue_name):
    """ Get memory information from qhost
    """
    rootxml = ET.fromstring(xmlstring)
    my_machine_dict = {}
    # on some machines rootxml.tag looks like "{...}qhost" where the "{...}" gets prepended to all attributes
    rootTag = rootxml.tag.rstrip("qhost")
    for host in rootxml.findall(rootTag + 'host'):
        # find all hosts supporting queues
        for queues in host.findall(rootTag + 'queue'):
            # if the user specified queue matches that in the xml:
            if not queue_name or any(q in queues.attrib['name'] for q in queue_name.split(",")):
                my_machine_dict[host.attrib['name']] = {}
                # values from xml for number of processors and mem_total on each machine
                for hostvalues in host.findall(rootTag + 'hostvalue'):
                    if('mem_total' == hostvalues.attrib['name']):
                        if hostvalues.text.lower().endswith('g'):
                            multip = 1
                        elif hostvalues.text.lower().endswith('m'):
                            multip = 1 / float(1024)
                        elif hostvalues.text.lower().endswith('t'):
                            multip = 1024
                        else:
                            raise Exception("Unrecognized suffix in mem_tot from SGE")
                        my_machine_dict[host.attrib['name']]['mem_total'] = \
                                float(hostvalues.text[:-1]) * float(multip)
                break
    return my_machine_dict

def _combine_machine_info(xs):
    if len(xs) == 1:
        return xs[0]
    else:
        raise NotImplementedError("Add logic to pick specification from non-homogeneous clusters.")

def get_info(dirs, parallel, resources=None):
    """Retrieve cluster or local filesystem resources from pre-retrieved information.
    """
    # Allow custom specification of cores/memory in resources
    if resources and isinstance(resources, dict) and "machine" in resources:
        minfo = resources["machine"]
        assert "memory" in minfo, "Require memory specification (Gb) in machine resources: %s" % minfo
        assert "cores" in minfo, "Require core specification in machine resources: %s" % minfo
        return minfo
    if parallel["type"] in ["ipython"] and not parallel["queue"] == "localrun":
        cache_file = _get_cache_file(dirs, parallel)
        if utils.file_exists(cache_file):
            with open(cache_file) as in_handle:
                minfo = yaml.safe_load(in_handle)
            return _combine_machine_info(minfo)
        else:
            return {}
    else:
        return _combine_machine_info(machine_info())

def machine_info():
    """Retrieve core and memory information for the current machine.
    """
    import psutil
    BYTES_IN_GIG = 1073741824.0
    free_bytes = psutil.virtual_memory().total
    return [{"memory": float("%.1f" % (free_bytes / BYTES_IN_GIG)), "cores": multiprocessing.cpu_count(),
             "name": socket.gethostname()}]

def open_file_limit():
    return resource.getrlimit(resource.RLIMIT_NOFILE)[0]
