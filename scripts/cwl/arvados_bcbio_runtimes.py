#!/usr/bin/env python
"""Summarize runtimes for bcbio Arvados CWL runs.

Usage:
  arvados_bcbio_runtimes.py <run_UUID>
"""
import csv
import json
import math
import operator
import os
import pprint
import sys

import arrow

def main(run_uuid):
    out_file = "%s-stats.csv" % run_uuid
    writer = csv.writer(open(out_file, "w"))
    client = _get_api_client()
    from arvados.collection import Collection
    req = client.container_requests().get(uuid=run_uuid).execute()
    jobs = client.container_requests().list(
            filters=[["requesting_container_uuid", "=", req["container_uuid"]]], limit=10000).execute()["items"]
    group_runtimes = {}
    writer.writerow(["group", "run", "time", "walltime", "cores", "memory"])
    for job in jobs:
        try:
            j1, j2 = job["name"].rsplit("_", 1)
            int(j2)
            name = j1
        except ValueError:
            name = job["name"]
        sample, vc = get_sample_variantcaller(job)
        name = "%s-%s" % (sample, name)
        if job["name"].startswith("variantcall_batch") and vc:
            name += "-%s" % vc
        print(name, job["uuid"], job["log_uuid"])
        log = client.collections().get(uuid=job["log_uuid"]).execute()
        logc = Collection(log["portable_data_hash"])
        machine_info = get_machine_info(logc)
        runtime = get_runtime(logc)
        writer.writerow([name, job["name"], str(runtime).rsplit(":", 1)[0], "",
                         machine_info["cores"], machine_info["memory"]])
        if name not in group_runtimes:
            group_runtimes[name] = [runtime]
        else:
            group_runtimes[name].append(runtime)
    for k in sorted(group_runtimes.keys()):
        writer.writerow([k, "",
                         str(reduce(operator.add, group_runtimes[k])).rsplit(":", 1)[0],
                         str(max(group_runtimes[k])).rsplit(":", 1)[0],
                         "", ""])

def get_machine_info(logc):
    procs = []
    mem = None
    with logc.open("node-info.txt") as in_handle:
        for line in in_handle:
            if line.startswith("node-info  processor"):
                procs.append(int(line.strip().split()[-1]))
            elif line.startswith("node-info  MemAvailable:"):
                mem = float(line.split()[-2]) / 1024.0 / 1024.0
    return {"cores": max(procs) + 1,
            "memory": "%sG" % int(math.floor(mem))}

def get_runtime(logc):
    with logc.open("crunchstat.txt") as in_handle:
        tstart = arrow.get(in_handle.readline().split()[0])
        for line in in_handle:
            last = line
        tend = arrow.get(last.split()[0])
    return tend - tstart

def get_in_inputs(key, data):
    if isinstance(data, dict):
        for k, v in data.items():
            if k == key:
                return v
            elif isinstance(v, (list, tuple, dict)):
                out = get_in_inputs(key, v)
                if out:
                    return out
    elif isinstance(data, (list, tuple)):
        out = [get_in_inputs(key, x) for x in data]
        out = [x for x in out if x]
        if out:
            return out[0]

def get_sample_variantcaller(job):
    from arvados.collection import Collection
    needs_json = job["name"].startswith("variantcall_batch")
    sample = None
    for f in job["mounts"]:
        if f.endswith("-sort.bam"):
            sample = os.path.basename(f).replace("-sort.bam", "")
    if not sample or needs_json:
        for k, v in job["mounts"].items():
            if k.endswith("cwl.inputs.json"):
                c = Collection(v["portable_data_hash"])
                with c.open("cwl.inputs.json", "r") as in_handle:
                    inputs = json.load(in_handle)
                return [get_in_inputs("description", inputs),
                        get_in_inputs("config__algorithm__variantcaller", inputs)]
    else:
        return sample, None
    return None, None

def _get_api_client():
    import arvados
    return arvados.api("v1")

if __name__ == "__main__":
    main(*sys.argv[1:])
