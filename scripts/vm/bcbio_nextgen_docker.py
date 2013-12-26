#!/usr/bin/env python
"""Install bcbio-nextgen, using code and tools isolated in a docker container.

Work in progress script to explore the best ways to integrate docker isolated
software with external data.
"""
import argparse
import grp
import json
import os
import pwd
import subprocess
import sys
import time

import requests

# ## Installation

def install(args):
    port = 8085
    image = "chapmanb/bcbio-nextgen-devel"
    docker_biodata_dir = "/mnt/biodata"
    pull(image)
    cid = None
    try:
        cid = start(image, port, args.datadir, docker_biodata_dir)
        wait(port)
        print "Running data installation with docker container", cid
        r = install_data(args, port)
        if r.status_code != 200:
            print "Problem installing data. For detailed logs, run: "
            print "docker logs {0}".format(cid)
    finally:
        if cid:
            stop(cid)

def install_data(args, port):
    payload = json.dumps({"genomes": args.genomes, "aligners": args.aligners,
                          "install_data": args.install_data})
    return requests.get("http://localhost:{port}/install".format(port=port), params={"args": payload})

def pull(image):
    print "Retrieving bcbio-nextgen docker images with code and tools"
    subprocess.check_call(["docker", "pull", image])

# ## Start and stop bcbio-nextgen docker container

def start(image, port, datadir, docker_biodata_dir):
    data_bind = "{datadir}:{docker_biodata_dir}".format(**locals())
    cmd = ("docker run -d -p {port}:{port} -v {data_bind} {image} /bin/bash -c '" +
           user_create_cmd() +
           "bcbio_nextgen.py server --port={port} --biodata_dir={docker_biodata_dir}\"'")
    process = subprocess.Popen(cmd.format(**locals()), shell=True, stdout=subprocess.PIPE)
    cid, _ = process.communicate()
    return cid.rstrip()

def user_create_cmd():
    user = pwd.getpwuid(os.getuid())
    group = grp.getgrgid(os.getgid())
    container_bcbio_dir = "/usr/local/share"
    homedir = "/home/{user.pw_name}".format(**locals())
    cmd = ("addgroup --gid {group.gr_gid} {group.gr_name} && "
           "useradd -m -d {homedir} -g {group.gr_gid} -o -u {user.pw_uid} {user.pw_name} && "
           "chown -R {user.pw_name} {container_bcbio_dir} && "
           "cd {homedir} && "
           "su - -s /bin/bash {user.pw_name} -c \"")
    return cmd.format(**locals())

def wait(port):
    """Wait for server to start.
    """
    num_tries = 0
    max_tries = 10
    while 1:
        try:
            requests.get("http://localhost:{port}".format(**locals()))
            break
        except requests.exceptions.ConnectionError:
            if num_tries > max_tries:
                raise
            else:
                num_tries += 1
                time.sleep(2)

def stop(cid):
    subprocess.check_call(["docker", "kill", cid])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Automatic installation for bcbio-nextgen pipelines, with docker.")
    subparsers = parser.add_subparsers(title="[sub-commands]")
    # installation
    parser_i = subparsers.add_parser("install", help="Install or upgrade bcbio-nextgen docker container and data.")
    parser_i.add_argument("datadir", help="Directory to install genome data",
                          type=lambda x: (os.path.abspath(os.path.expanduser(x))))
    parser_i.add_argument("--genomes", help="Genomes to download",
                          action="append", default=["GRCh37"])
    parser_i.add_argument("--aligners", help="Aligner indexes to download",
                          action="append", default=["bwa"])
    parser_i.add_argument("--nodata", help="Do not install data dependencies",
                          dest="install_data", action="store_false", default=True)
    parser_i.set_defaults(func=install)
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        args.func(args)

"""
def docker_py_start():
    # XXX Does not appear to bind ports correctly
    # Swap to API instead of command line calls later as it stabilizes
    import docker
    ports = {8085: ("0.0.0.0", 8085)}
    binds = {"/usr/local/share/bcbio_nextgen": "/mnt/biodata"}
    client = docker.Client()
    cid = client.create_container("chapmanb/bcbio-nextgen-devel",
                                  command="bcbio_nextgen.py server --port=%s" % ports.keys()[0],
                                  volumes=binds.values(), ports=ports.keys())
    client.start(cid, port_bindings=ports, binds=binds)
"""
