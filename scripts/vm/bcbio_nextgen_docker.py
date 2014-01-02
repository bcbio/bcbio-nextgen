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

import docker
import requests
import yaml

# default information about docker container
DOCKER = {"port": 8085,
          "biodata_dir": "/mnt/biodata",
          "input_dir": "/mnt/inputs",
          "work_dir": "/mnt/work",
          "image": "chapmanb/bcbio-nextgen-devel"}

# ## Running analysis

def run(args):
    """Run a full analysis on a local machine, utilizing multiple cores.
    """
    args = update_check_args(args, "Could not run analysis.")
    with open(args.sample_config) as in_handle:
        sample_config, mounts = update_config_mounts(yaml.load(in_handle), DOCKER["input_dir"])
    mounts += prepare_system_mounts(args.datadir, DOCKER["biodata_dir"])
    mounts.append("%s:%s" % (os.getcwd(), DOCKER["work_dir"]))
    system_config, system_mounts = read_system_config(args, DOCKER)
    dockerc = docker.Client()
    with bcbio_docker(DOCKER, mounts + system_mounts, args) as cid:
        print("Running analysis using docker container: %s" % cid)
        payload = {"work_dir": DOCKER["work_dir"],
                   "system_config": system_config,
                   "sample_config": sample_config,
                   "numcores": args.numcores}
        r = requests.get("http://localhost:{port}/run".format(port=args.port), params={"args": json.dumps(payload)})
        run_id = r.text
        clogs = dockerc.logs(cid, stream=True)
        # monitor processing status, writing logging information
        for log_info in clogs:
            print(log_info.rstrip())
            r = requests.get("http://localhost:{port}/status".format(port=args.port), params={"run_id": run_id})
            if r.text != "running":
                break

def read_system_config(args, DOCKER):
    if args.systemconfig:
        f = args.systemconfig
    else:
        f = os.path.join(args.datadir, "galaxy", "bcbio_system.yaml")
    with open(f) as in_handle:
        config = yaml.load(in_handle)
    # Map external galaxy specifications over to docker container
    mounts = []
    for k in ["galaxy_config"]:
        if k in config:
            dirname, base = os.path.split(os.path.normpath(os.path.realpath(config[k])))
            container_dir = os.path.join(DOCKER["input_dir"], "system", "galaxy", k)
            mounts.append("%s:%s" % (dirname, container_dir))
            mounts.extend(_find_genome_directory_mounts(dirname, container_dir))
            config[k] = os.path.join(container_dir, base)
    return config, mounts

def _find_genome_directory_mounts(dirname, container_dir):
    """Handle external non-docker installed biodata located relative to config directory.

    Need a general way to handle mounting these and adjusting paths, but this handles
    the special case used in testing.
    """
    mounts = []
    sam_loc = os.path.join(dirname, "tool-data", "sam_fa_indices.loc")
    genome_dir = None
    if os.path.exists(sam_loc):
        with open(sam_loc) as in_handle:
            for line in in_handle:
                if line.startswith("index"):
                    genome_dir = line.split()[-1].strip()
                    break
    if genome_dir and not os.path.isabs(genome_dir):
        rel_genome_dir = os.path.dirname(os.path.dirname(os.path.dirname(genome_dir)))
        mounts.append("%s:%s" % (os.path.normpath(os.path.join(os.path.dirname(sam_loc), rel_genome_dir)),
                                 os.path.normpath(os.path.join(os.path.join(container_dir, "tool-data"),
                                                               rel_genome_dir))))
    return mounts

def update_config_mounts(config, input_dir):
    """Update input configuration with local docker container mounts.
    Maps input files into docker mounts and resolved relative and symlinked paths.
    """
    absdetails = []
    directories = []
    for d in config["details"]:
        d = abs_file_paths(d, base_dirs=[args.fcdir] if args.fcdir else None,
                           ignore=["description", "analysis", "resources",
                                   "genome_build", "lane"])
        d["algorithm"] = abs_file_paths(d["algorithm"], base_dirs=[args.fcdir] if args.fcdir else None,
                                        ignore=["variantcaller", "realign", "recalibrate",
                                                "phasing", "svcaller"])
        absdetails.append(d)
        directories.extend(_get_directories(d))
    mounts = {}
    for i, d in enumerate(sorted(set(directories))):
        mounts[d] = os.path.join(input_dir, str(i))
    config["details"] = [_remap_directories(d, mounts) for d in absdetails]
    return config, ["%s:%s" % (k, v) for k, v in mounts.iteritems()]

def _remap_directories(xs, mounts):
    """Remap files to point to internal docker container mounts.
    """
    if not isinstance(xs, dict):
        return xs
    out = {}
    for k, v in xs.iteritems():
        if isinstance(v, dict):
            out[k] = _remap_directories(v, mounts)
        elif v and isinstance(v, basestring) and os.path.exists(v) and os.path.isabs(v):
            dirname, basename = os.path.split(v)
            out[k] = os.path.join(mounts[dirname], basename)
        elif v and isinstance(v, (list, tuple)) and os.path.exists(v[0]):
            ready_vs = []
            for x in v:
                dirname, basename = os.path.split(x)
                ready_vs.append(os.path.join(mounts[dirname], basename))
            out[k] = ready_vs
        else:
            out[k] = v
    return out

def _get_directories(xs):
    """Retrieve all directories specified in an input file.
    """
    out = []
    if not isinstance(xs, dict):
        return out
    for k, v in xs.iteritems():
        if isinstance(v, dict):
            out.extend(_get_directories(v))
        elif v and isinstance(v, basestring) and os.path.exists(v) and os.path.isabs(v):
            out.append(os.path.dirname(v))
        elif v and isinstance(v, (list, tuple)) and os.path.exists(v[0]):
            out.extend(os.path.dirname(x) for x in v)
    return out

def _normalize_path(x, base_dirs):
    for base_dir in base_dirs:
        if os.path.exists(os.path.join(base_dir, x)):
            return os.path.normpath(os.path.realpath(os.path.join(base_dir, x)))
    return None

def abs_file_paths(xs, base_dirs=None, ignore=[]):
    """Expand files to be absolute, non-symlinked file paths.
    """
    if not isinstance(xs, dict):
        return xs
    base_dirs = base_dirs if base_dirs else []
    base_dirs.append(os.getcwd())
    ignore_keys = set(ignore if ignore else [])
    out = {}
    for k, v in xs.iteritems():
        if k not in ignore_keys and v and isinstance(v, basestring) and _normalize_path(v, base_dirs):
            out[k] = _normalize_path(v, base_dirs)
        elif k not in ignore_keys and v and isinstance(v, (list, tuple)) and _normalize_path(v[0], base_dirs):
            out[k] = [_normalize_path(x, base_dirs) for x in v]
        else:
            out[k] = v
    return out

# ## Installation

def install(args):
    args = update_check_args(args, "bcbio-nextgen not upgraded.")
    if args.install_tools:
        pull(DOCKER["image"])
    success = True
    mounts = prepare_system_mounts(args.datadir, DOCKER["biodata_dir"])
    with bcbio_docker(DOCKER, mounts, args) as cid:
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

def update_check_args(args, command_info):
    args = add_defaults(args)
    if not args.datadir:
        print("Must specify a `--datadir` or save the default location with `saveconfig`.\n" + command_info)
        sys.exit(1)
    return args

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
def bcbio_docker(dconf, mounts, args):
    """Provide a running bcbio-nextgen docker server with automatic stop on completion.
    """
    cid = None
    if args.develrepo:
        yield start_devel(dconf["image"], args.port, dconf["port"], mounts,
                          dconf["biodata_dir"], args.develrepo)
    else:
        try:
            cid = start(dconf["image"], args.port, dconf["port"], mounts, dconf["biodata_dir"])
            wait(dconf["port"])
            yield cid
        finally:
            if cid:
                stop(cid)

def start_devel(image, hport, cport, mounts, docker_biodata_dir, repo):
    """Start a docker container for development, attached to the provided code repo.
    Uses a standard name 'bcbio-develrepo' to avoid launching multiple images on
    re-run.
    """
    name = "bcbio-develrepo"
    # look for existing running processes
    process = subprocess.Popen(["docker", "ps"], stdout=subprocess.PIPE)
    containers, _ = process.communicate()
    for line in containers.split("\n"):
        if line.find(name) >= 0:
            return line.split()[0]
    # start a new container running bash for development
    mounts.append("%s:/tmp/bcbio-nextgen" % repo)
    mounts = " ".join("-v %s" % x for x in mounts)
    chown_cmd = ("chown -R {user.pw_name} /usr/local/share/bcbio-nextgen/anaconda/lib/python2.7/site-packages && "
                 "chown -R {user.pw_name} /usr/local/share/bcbio-nextgen/anaconda/bin && ")
    cmd = ("docker run -d -i -t -name {name} -p {hport}:{cport} {mounts} {image} "
           "/bin/bash -c '"
           + user_create_cmd(chown_cmd) +
           "/bin/bash"
           "\"'")
    process = subprocess.Popen(cmd.format(**locals()), shell=True, stdout=subprocess.PIPE)
    cid, _ = process.communicate()
    return cid.rstrip()

def start(image, hport, cport, mounts, docker_biodata_dir):
    mounts = " ".join("-v %s" % x for x in mounts)
    cmd = ("docker run -d -p {hport}:{cport} {mounts} {image} "
           "/bin/bash -c '" + user_create_cmd() +
           "bcbio_nextgen.py server --port={cport} --biodata_dir={docker_biodata_dir}"
           "\"'")
    process = subprocess.Popen(cmd.format(**locals()), shell=True, stdout=subprocess.PIPE)
    cid, _ = process.communicate()
    return cid.rstrip()

def prepare_system_mounts(datadir, docker_biodata_dir):
    """Create set of system mountpoints to link into Docker container.
    """
    mounts = []
    for d in ["genomes", "liftOver", "gemini_data", "galaxy"]:
        cur_d = os.path.normpath(os.path.realpath(os.path.join(datadir, d)))
        if not os.path.exists(cur_d):
            os.makedirs(cur_d)
        mounts.append("{cur_d}:{docker_biodata_dir}/{d}".format(**locals()))
    return mounts

def user_create_cmd(chown_cmd=""):
    """Create a user on the docker container with equivalent UID/GIDs to external user.
    """
    user = pwd.getpwuid(os.getuid())
    group = grp.getgrgid(os.getgid())
    container_bcbio_dir = "/usr/local/share"
    homedir = "/home/{user.pw_name}".format(**locals())
    cmd = ("addgroup --gid {group.gr_gid} {group.gr_name} && "
           "useradd -m -d {homedir} -g {group.gr_gid} -o -u {user.pw_uid} {user.pw_name} && "
           + chown_cmd +
           "su - -s /bin/bash {user.pw_name} -c \"cd {homedir} && ")
    return cmd.format(**locals())

def wait(port):
    """Wait for server to start.
    """
    num_tries = 0
    max_tries = 40
    while 1:
        try:
            requests.get("http://localhost:{port}/status".format(**locals()),
                         params={"run_id": "checkup"})
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
    parser.add_argument("--develrepo", help=("Specify a development repository to link. "
                                             "Used for debugging and development"))
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
    parser_r.add_argument("sample_config", help="YAML file with details about samples to process.")
    parser_r.add_argument("--fcdir", help="A directory of Illumina output or fastq files to process")
    parser_r.add_argument("--systemconfig", help="Global YAML configuration file specifying system details. "
                          "Defaults to installed bcbio_system.yaml.")
    parser_r.add_argument("-n", "--numcores", help="Total cores to use for processing",
                          type=int, default=1)
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
    ports = {8085: ("0.0.0.0", 8085)}
    binds = {"/usr/local/share/bcbio_nextgen": "/mnt/biodata"}
    client = docker.Client()
    cid = client.create_container("chapmanb/bcbio-nextgen-devel",
                                  command="bcbio_nextgen.py server --port=%s" % ports.keys()[0],
                                  volumes=binds.values(), ports=ports.keys())
    client.start(cid, port_bindings=ports, binds=binds)
"""
