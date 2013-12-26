#!/usr/bin/env python
"""Install bcbio-nextgen, using code and tools isolated in a docker container.

Work in progress script to explore the best ways to integrate docker isolated
software with external data.
"""
import argparse
import json
import os
import subprocess
import sys
import time

import requests

def main(args):
    port = 8085
    image = "chapmanb/bcbio-nextgen-devel"
    docker_biodata_dir = "/mnt/biodata"
    pull(image)
    try:
        cid = start(image, port, args.datadir, docker_biodata_dir)
        wait(port)
        print "Running data installation with docker container", cid
        install_data(args, port)
    finally:
        stop(cid)

def start(image, port, datadir, docker_biodata_dir):
    data_bind = "{datadir}:{docker_biodata_dir}".format(**locals())
    cmd = ("docker run -d -p {port}:{port} -v {data_bind} {image} "
           "bcbio_nextgen.py server --port={port} --biodata_dir={docker_biodata_dir}")
    process = subprocess.Popen(cmd.format(**locals()), shell=True, stdout=subprocess.PIPE)
    cid, _ = process.communicate()
    return cid.rstrip()

def install_data(args, port):
    payload = json.dumps({"genomes": args.genomes, "aligners": args.aligners,
                          "install_data": args.install_data})
    r = requests.get("http://localhost:{port}/install".format(port=port), params={"args": payload})
    print r

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

def pull(image):
    print "Retrieving bcbio-nextgen docker images with code and tools"
    subprocess.check_call(["docker", "pull", image])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Automatic installation for bcbio-nextgen pipelines, with docker.")
    parser.add_argument("datadir", help="Directory to install genome data",
                        type=lambda x: (os.path.abspath(os.path.expanduser(x))))
    parser.add_argument("--genomes", help="Genomes to download",
                        action="append", default=["GRCh37"])
    parser.add_argument("--aligners", help="Aligner indexes to download",
                        action="append", default=["bwa"])
    parser.add_argument("--nodata", help="Do not install data dependencies",
                        dest="install_data", action="store_false", default=True)
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        main(parser.parse_args())

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
