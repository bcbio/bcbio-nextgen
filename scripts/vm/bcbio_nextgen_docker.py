#!/usr/bin/env python
"""Run bcbio_nextgen, using code and tools isolated in a docker container.
"""
import subprocess

def main():
    port = 8085
    image = "chapmanb/bcbio-nextgen-devel"
    data_bind = "/usr/local/share/bcbio_nextgen:/mnt/biodata"
    cmd = ("docker run -d -p {port}:{port} -v {data_bind} {image} "
           "bcbio_nextgen.py server --port={port}")
    subprocess.check_call(cmd.format(**locals()), shell=True)

def main_docker_py():
    # XXX Does not appear to bind ports correctly
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
