#!/usr/bin/env python
"""Run and install bcbio-nextgen, using code and tools isolated in a docker container.

Work in progress script to explore the best ways to integrate docker isolated
software with external data.
"""
import argparse
import contextlib
import grp
import json
import os
import pwd
import subprocess
import sys
import time

import requests
import yaml

# default information about docker container
DOCKER = {"port": 8085,
          "biodata_dir": "/mnt/biodata",
          "image": "chapmanb/bcbio-nextgen-devel"}

# ## Running analysis

def run(args):
    args = add_defaults(args)
    with bcbio_docker(DOCKER, args) as cid:
        print "***", cid

# ## Installation

def install(args):
    args = add_defaults(args)
    if not args.datadir:
        print("Must specify a `--datadir` for installation or save location with `saveconfig`.\n"
              "bcbio-nextgen not upgraded.")
        return
    if args.install_tools:
        pull(DOCKER["image"])
    success = True
    with bcbio_docker(DOCKER, args) as cid:
        print("Running data installation with docker container: %s" % cid)
        r = install_data(args, DOCKER["port"])
        if r is None or r.status_code != 200:
            success = False
            print("Problem installing data. For detailed logs, run:\n"
                  "docker logs {0}".format(cid))
    if success:
        print("bcbio-nextgen successfully upgraded")

def install_data(args, port):
    payload = json.dumps({"genomes": args.genomes, "aligners": args.aligners,
                          "install_data": args.install_data})
    try:
        return requests.get("http://localhost:{port}/install".format(port=port), params={"args": payload})
    except requests.exceptions.ConnectionError:
        return None

def pull(image):
    print("Retrieving bcbio-nextgen docker images with code and tools")
    subprocess.check_call(["docker", "pull", image])

# ## Defaults

TOSAVE_DEFAULTS = {"port": 8085, "datadir": None}

def save_defaults(args):
    """Save user specific defaults to a yaml configuration file.
    """
    out = get_defaults()
    for k in TOSAVE_DEFAULTS:
        karg = getattr(args, k, None)
        if karg and karg != TOSAVE_DEFAULTS[k]:
            out[k] = karg
    if len(out) > 0:
        with open(_get_config_file(just_filename=True), "w") as out_handle:
            yaml.dump(out, out_handle, default_flow_style=False, allow_unicode=False)

def add_defaults(args):
    """Add user configured defaults to supplied command line arguments.
    """
    config_defaults = get_defaults()
    for k in TOSAVE_DEFAULTS:
        karg = getattr(args, k, None)
        if not karg or karg == TOSAVE_DEFAULTS[k]:
            if k in config_defaults:
                setattr(args, k, config_defaults[k])
    return args

def get_defaults():
    """Retrieve saved default configurations.
    """
    config_file = _get_config_file()
    if config_file:
        with open(config_file) as in_handle:
            return yaml.load(in_handle)
    else:
        return {}

def _get_config_file(just_filename=False):
    """Retrieve standard user configuration file.
    Uses location from appdirs (https://github.com/ActiveState/appdirs). Could
    pull this in as dependency for more broad platform support.
    """
    config_dir = os.path.join(os.getenv('XDG_CONFIG_HOME', os.path.expanduser("~/.config")),
                              "bcbio-nextgen")
    if not os.path.exists(config_dir):
        os.makedirs(config_dir)
    config_file = os.path.join(config_dir, "bcbio-docker-config.yaml")
    if just_filename or os.path.exists(config_file):
        return config_file
    else:
        return None


# ## Start and stop bcbio-nextgen docker container

@contextlib.contextmanager
def bcbio_docker(dconf, args):
    """Provide a running bcbio-nextgen docker server with automatic stop on completion.
    """
    cid = None
    try:
        cid = start(dconf["image"], args.port, dconf["port"], args.datadir, dconf["biodata_dir"])
        wait(dconf["port"])
        yield cid
    finally:
        if cid:
            stop(cid)

def start(image, hport, cport, datadir, docker_biodata_dir):
    mounts = " ".join("-v %s" % x for x in prepare_mounts(datadir, docker_biodata_dir))
    cmd = ("docker run -d -p {hport}:{cport} {mounts} {image} /bin/bash -c '" +
           user_create_cmd() +
           "bcbio_nextgen.py server --port={cport} --biodata_dir={docker_biodata_dir}\"'")
    process = subprocess.Popen(cmd.format(**locals()), shell=True, stdout=subprocess.PIPE)
    cid, _ = process.communicate()
    return cid.rstrip()

def prepare_mounts(datadir, docker_biodata_dir):
    """Create set of system mountpoints to link into Docker container.
    """
    mounts = []
    for d in ["genomes", "liftOver", "gemini_data", "galaxy"]:
        cur_d = os.path.normpath(os.path.realpath(os.path.join(datadir, d)))
        if not os.path.exists(cur_d):
            os.makedirs(cur_d)
        mounts.append("{cur_d}:{docker_biodata_dir}/{d}".format(**locals()))
    return mounts

def user_create_cmd():
    """Create a user on the docker container with equivalent UID/GIDs to external user.
    """
    user = pwd.getpwuid(os.getuid())
    group = grp.getgrgid(os.getgid())
    container_bcbio_dir = "/usr/local/share"
    homedir = "/home/{user.pw_name}".format(**locals())
    cmd = ("addgroup --gid {group.gr_gid} {group.gr_name} && "
           "useradd -m -d {homedir} -g {group.gr_gid} -o -u {user.pw_uid} {user.pw_name} && "
           "su - -s /bin/bash {user.pw_name} -c \"cd {homedir} && ")
    return cmd.format(**locals())

def wait(port):
    """Wait for server to start.
    """
    num_tries = 0
    max_tries = 40
    while 1:
        try:
            requests.get("http://localhost:{port}".format(**locals()))
            break
        except requests.exceptions.ConnectionError:
            if num_tries > max_tries:
                raise
            else:
                num_tries += 1
                time.sleep(1)

def stop(cid):
    subprocess.check_call(["docker", "kill", cid], stdout=subprocess.PIPE)
    #subprocess.check_call(["docker", "rm", cid], stdout=subprocess.PIPE)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Automatic installation for bcbio-nextgen pipelines, with docker.")
    parser.add_argument("--port", default=8085, help="External port to connect to docker image.")
    parser.add_argument("--datadir", help="Directory to install genome data and associated files.",
                        type=lambda x: (os.path.abspath(os.path.expanduser(x))))
    subparsers = parser.add_subparsers(title="[sub-commands]")
    # installation
    parser_i = subparsers.add_parser("install", help="Install or upgrade bcbio-nextgen docker container and data.")
    parser_i.add_argument("--genomes", help="Genomes to download",
                          action="append", default=["GRCh37"])
    parser_i.add_argument("--aligners", help="Aligner indexes to download",
                          action="append", default=["bwa"])
    parser_i.add_argument("--data", help="Install or upgrade data dependencies",
                          dest="install_data", action="store_true", default=False)
    parser_i.add_argument("--tools", help="Install or upgrade tool dependencies",
                          dest="install_tools", action="store_true", default=False)
    parser_i.set_defaults(func=install)
    # running
    parser_r = subparsers.add_parser("run", help="Run an automated analysis.")
    parser_r.add_argument("run_config", help="YAML file with details about samples to process.")
    parser_r.add_argument("--fcdir", help="A directory of Illumina output or fastq files to process")
    parser_r.add_argument("--globalconfig", help="Global YAML configuration file specifying details. "
                          "Defaults to installed bcbio_system.yaml.")
    parser_r.set_defaults(func=run)
    # configuration
    parser_c = subparsers.add_parser("saveconfig", help="Save standard configuration variables for current user. "
                                     "Avoids need to specify on the command line in future runs.")
    parser_c.set_defaults(func=save_defaults)
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
