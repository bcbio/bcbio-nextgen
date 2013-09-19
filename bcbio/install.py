"""Handle installation and updates of bcbio-nextgen, third party software and data.

Enables automated installation tool and in-place updates to install additional
data and software.
"""
import contextlib
import os
import shutil
import subprocess
import sys

import requests
import yaml

from bcbio import utils
from bcbio.pipeline import genome

REMOTES = {
    "requirements": "https://raw.github.com/chapmanb/bcbio-nextgen/master/requirements.txt",
    "gitrepo": "git://github.com/chapmanb/bcbio-nextgen.git",
    "cloudbiolinux": "https://github.com/chapmanb/cloudbiolinux.git",
    "genome_resources": "https://raw.github.com/chapmanb/bcbio-nextgen/master/config/genomes/%s-resources.yaml",
    }

def upgrade_bcbio(args):
    """Perform upgrade of bcbio to latest release, or from GitHub development version.

    Handles bcbio, third party tools and data.
    """
    args = add_install_defaults(args)
    pip_bin = os.path.join(os.path.dirname(sys.executable), "pip")
    if args.upgrade not in ["skip"]:
        _update_conda_packages()
    if args.upgrade in ["skip"]:
        pass
    elif args.upgrade in ["stable", "system"]:
        print("Upgrading bcbio-nextgen to latest stable version")
        sudo_cmd = [] if args.upgrade == "stable" else ["sudo"]
        subprocess.check_call(sudo_cmd + [pip_bin, "install", "-r", REMOTES["requirements"]])
    else:
        print("Upgrading bcbio-nextgen to latest development version")
        subprocess.check_call([pip_bin, "install", "--upgrade", "--no-deps",
                               "git+%s#egg=bcbio-nextgen" % REMOTES["gitrepo"]])
        subprocess.check_call([pip_bin, "install", "git+%s#egg=bcbio-nextgen" % REMOTES["gitrepo"]])

    if args.tooldir:
        with bcbio_tmpdir():
            print("Upgrading third party tools to latest versions")
            upgrade_thirdparty_tools(args, REMOTES)
    if args.install_data:
        with bcbio_tmpdir():
            print("Upgrading bcbio-nextgen data files")
            upgrade_bcbio_data(args, REMOTES)
    save_install_defaults(args)

def _default_deploy_args(args):
    flavors = {"minimal": "ngs_pipeline_minimal",
               "full": "ngs_pipeline"}
    return {"flavor": flavors[args.tooldist],
            "vm_provider": "novm",
            "hostname": "localhost",
            "fabricrc_overrides" : {"edition": "minimal",
                                    "use_sudo": args.sudo,
                                    "distribution": args.distribution or "__auto__",
                                    "dist_name": "__auto__"}}

def _update_conda_packages():
    """If installed in an anaconda directory, upgrade conda packages.
    """
    conda_bin = os.path.join(os.path.dirname(sys.executable), "conda")
    pkgs = ["biopython", "boto", "cython", "distribute", "ipython", "nose", "numpy",
            "pycrypto", "pip", "pysam", "pyyaml", "pyzmq", "requests"]
    if os.path.exists(conda_bin):
        subprocess.check_call([conda_bin, "install", "--yes"] + pkgs)

def _get_data_dir():
    base_dir = os.path.realpath(os.path.dirname(os.path.dirname(sys.executable)))
    if "anaconda" not in os.path.basename(base_dir) and "virtualenv" not in os.path.basename(base_dir):
        raise ValueError("Cannot update data for bcbio-nextgen not installed by installer.")
    return os.path.dirname(base_dir)

def upgrade_bcbio_data(args, remotes):
    """Upgrade required genome data files in place.
    """
    data_dir = _get_data_dir()
    s = _default_deploy_args(args)
    s["actions"] = ["setup_biodata"]
    s["fabricrc_overrides"]["data_files"] = data_dir
    s["fabricrc_overrides"]["galaxy_home"] = os.path.join(data_dir, "galaxy")
    cbl = get_cloudbiolinux(remotes)
    s["genomes"] = _get_biodata(cbl["biodata"], args)
    sys.path.insert(0, cbl["dir"])
    cbl_deploy = __import__("cloudbio.deploy", fromlist=["deploy"])
    cbl_deploy.deploy(s)
    _upgrade_genome_resources(s["fabricrc_overrides"]["galaxy_home"],
                              remotes["genome_resources"])

def _upgrade_genome_resources(galaxy_dir, base_url):
    """Retrieve latest version of genome resource YAML configuration files.
    """
    for dbkey, ref_file in genome.get_builds(galaxy_dir):
        # Check for a remote genome resources file
        remote_url = base_url % dbkey
        r = requests.get(remote_url)
        if r.status_code == requests.codes.ok:
            local_file = os.path.join(os.path.dirname(ref_file), os.path.basename(remote_url))
            if os.path.exists(local_file):
                with open(local_file) as in_handle:
                    local_config = yaml.load(in_handle)
                remote_config = yaml.load(r.text)
                needs_update = remote_config["version"] > local_config.get("version", 0)
                if needs_update:
                    shutil.move(local_file, local_file + ".old%s" % local_config.get("version", 0))
            else:
                needs_update = True
            if needs_update:
                print("Updating %s genome resources configuration" % dbkey)
                with open(local_file, "w") as out_handle:
                    out_handle.write(r.text)

def _get_biodata(base_file, args):
    with open(base_file) as in_handle:
        config = yaml.load(in_handle)
    config["install_liftover"] = False
    config["genome_indexes"] = args.aligners
    config["genomes"] = [g for g in config["genomes"] if g["dbkey"] in args.genomes]
    return config

def upgrade_thirdparty_tools(args, remotes):
    """Install and update third party tools used in the pipeline.
    """
    s = {"fabricrc_overrides": {"system_install": args.tooldir,
                                "local_install": os.path.join(args.tooldir, "local_install"),
                                "distribution": args.distribution,
                                "use_sudo": args.sudo,
                                "edition": "minimal"}}
    s = _default_deploy_args(args)
    s["actions"] = ["install_biolinux"]
    s["fabricrc_overrides"]["system_install"] = args.tooldir
    s["fabricrc_overrides"]["local_install"] = os.path.join(args.tooldir, "local_install")
    cbl = get_cloudbiolinux(remotes)
    sys.path.insert(0, cbl["dir"])
    cbl_deploy = __import__("cloudbio.deploy", fromlist=["deploy"])
    cbl_deploy.deploy(s)

# ## Store a local configuration file with upgrade details

def _get_install_config():
    """Return the YAML configuration file used to store upgrade information.
    """
    try:
        data_dir = _get_data_dir()
    except ValueError:
        return None
    config_dir = utils.safe_makedir(os.path.join(data_dir, "config"))
    return os.path.join(config_dir, "install-params.yaml")

def save_install_defaults(args):
    """Save installation information to make future upgrades easier.
    """
    install_config = _get_install_config()
    if install_config is None:
        return
    if utils.file_exists(install_config):
        with open(install_config) as in_handle:
            cur_config = yaml.load(in_handle)
    else:
        cur_config = {}
    if args.tooldist not in "minimal":
        cur_config["tooldist"] = args.tooldist
    if args.tooldir:
        cur_config["tooldir"] = args.tooldir
    cur_config["sudo"] = args.sudo
    for attr in ["genomes", "aligners"]:
        if not cur_config.get(attr):
            cur_config[attr] = []
        for x in getattr(args, attr):
            if x not in cur_config[attr]:
                cur_config[attr].append(x)
    with open(install_config, "w") as out_handle:
        yaml.dump(cur_config, out_handle, default_flow_style=False, allow_unicode=False)

def add_install_defaults(args):
    """Add any saved installation defaults to the upgrade.
    """
    install_config = _get_install_config()
    if install_config is None:
        return args
    with open(install_config) as in_handle:
        default_args = yaml.load(in_handle)
    if default_args.get("tooldist") and args.tooldist == "minimal":
        args.tooldir = default_args["tooldist"]
    for attr in ["genomes", "aligners"]:
        for x in default_args.get(attr, []):
            new_val =  getattr(args, attr)
            if x not in getattr(args, attr):
                new_val.append(x)
            setattr(args, attr, new_val)
    args.sudo = default_args["sudo"]
    return args

def add_subparser(subparsers):
    parser = subparsers.add_parser("upgrade", help="Install or upgrade bcbio-nextgen")
    parser.add_argument("--tooldir",
                        help="Directory to install 3rd party software tools. Leave unspecified for no tools",
                        type=lambda x: (os.path.abspath(os.path.expanduser(x))), default=None)
    parser.add_argument("--tooldist",
                        help="Type of tool distribution to install. Defaults to a minimum install.",
                        default="minimal",
                        choices=["minimal", "full"])
    parser.add_argument("-u", "--upgrade", help="Code version to upgrade",
                        choices = ["stable", "development", "system", "skip"], default="stable")
    parser.add_argument("--genomes", help="Genomes to download",
                        action="append", default=["GRCh37"])
    parser.add_argument("--aligners", help="Aligner indexes to download",
                        action="append", default=["bwa"])
    parser.add_argument("--nosudo", help="Specify we cannot use sudo for commands",
                        dest="sudo", action="store_false", default=True)
    parser.add_argument("--nodata", help="Do not install data dependencies",
                        dest="install_data", action="store_false", default=True)
    parser.add_argument("--distribution", help="Operating system distribution",
                        default="",
                        choices=["ubuntu", "debian", "centos", "scientificlinux", "macosx"])

def get_cloudbiolinux(remotes):
    base_dir = os.path.join(os.getcwd(), "cloudbiolinux")
    if not os.path.exists(base_dir):
        subprocess.check_call(["git", "clone", remotes["cloudbiolinux"]])
    return {"biodata": os.path.join(base_dir, "config", "biodata.yaml"),
            "dir": base_dir}

@contextlib.contextmanager
def bcbio_tmpdir():
    orig_dir = os.getcwd()
    work_dir = os.path.join(os.getcwd(), "tmpbcbio-install")
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    os.chdir(work_dir)
    yield work_dir
    os.chdir(orig_dir)
    shutil.rmtree(work_dir)
