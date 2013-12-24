#!/usr/bin/env python
"""Install bcbio-nextgen, using code and tools isolated in a docker container.

Work in progress script to explore the best ways to integrate docker isolated
software with external data.
"""
import subprocess

def main():
    image = "chapmanb/bcbio-nextgen-devel"
    #pull(image)
    cid = start(image)
    print(cid)
    #stop(cid)

def start(image):
    port = 8085
    data_bind = "/usr/local/share/bcbio_nextgen:/mnt/biodata"
    cmd = ("docker run -d -p {port}:{port} -v {data_bind} {image} "
           "bcbio_nextgen.py server --port={port}")
    process = subprocess.Popen(cmd.format(**locals()), shell=True, stdout=subprocess.PIPE)
    cid, _ = process.communicate()
    return cid

def stop(cid):
    subprocess.check_call(["docker", "kill", cid])

def pull(image):
    subprocess.check_call(["docker", "pull", image])

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

if __name__ == "__main__":
    main()
