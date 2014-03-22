"""Handle installation and updates of bcbio-nextgen, third party software and data.

Enables automated installation tool and in-place updates to install additional
data and software.
"""
import argparse
import collections
import contextlib
import datetime
from distutils.version import LooseVersion
import os
import shutil
import string
import subprocess
import sys

import requests
import yaml

from bcbio import broad, utils
from bcbio.pipeline import genome
from bcbio.variation import effects
from bcbio.provenance import programs

REMOTES = {
    "requirements": "https://raw.github.com/chapmanb/bcbio-nextgen/master/requirements.txt",
    "gitrepo": "git://github.com/chapmanb/bcbio-nextgen.git",
    "cloudbiolinux": "https://github.com/chapmanb/cloudbiolinux.git",
    "genome_resources": "https://raw.github.com/chapmanb/bcbio-nextgen/master/config/genomes/%s-resources.yaml",
    "snpeff_dl_url": ("http://downloads.sourceforge.net/project/snpeff/databases/v{snpeff_ver}/"
                      "snpEff_v{snpeff_ver}_{genome}.zip")}

Tool = collections.namedtuple("Tool", ["name", "fname"])

def upgrade_bcbio(args):
    """Perform upgrade of bcbio to latest release, or from GitHub development version.

    Handles bcbio, third party tools and data.
    """
    args = add_install_defaults(args)
    pip_bin = os.path.join(os.path.dirname(sys.executable), "pip")
    if args.upgrade in ["skip"]:
        pass
    elif args.upgrade in ["stable", "system"]:
        _update_conda_packages()
        print("Upgrading bcbio-nextgen to latest stable version")
        sudo_cmd = [] if args.upgrade == "stable" else ["sudo"]
        subprocess.check_call(sudo_cmd + [pip_bin, "install", "-r", REMOTES["requirements"]])
    else:
        _update_conda_packages()
        print("Upgrading bcbio-nextgen to latest development version")
        subprocess.check_call([pip_bin, "install", "git+%s#egg=bcbio-nextgen" % REMOTES["gitrepo"]])
        subprocess.check_call([pip_bin, "install", "--upgrade", "--no-deps",
                               "git+%s#egg=bcbio-nextgen" % REMOTES["gitrepo"]])

    if args.tooldir:
        with bcbio_tmpdir():
            print("Upgrading third party tools to latest versions")
            upgrade_thirdparty_tools(args, REMOTES)
    if args.install_data:
        with bcbio_tmpdir():
            print("Upgrading bcbio-nextgen data files")
            upgrade_bcbio_data(args, REMOTES)
    if args.isolate and args.tooldir:
        print("Installation directory not added to current PATH")
        print("  Add {t}/bin to PATH and {t}/lib to LD_LIBRARY_PATH".format(t=args.tooldir))
    save_install_defaults(args)
    args.datadir = _get_data_dir()
    _install_container_bcbio_system(args.datadir)
    return args

def _install_container_bcbio_system(datadir):
    """Install limited bcbio_system.yaml file for setting core and memory usage.

    Adds any non-specific programs to the exposed bcbio_system.yaml file, only
    when upgrade happening inside a docker container.
    """
    base_file = os.path.join(datadir, "config", "bcbio_system.yaml")
    if not os.path.exists(base_file):
        return
    expose_file = os.path.join(datadir, "galaxy", "bcbio_system.yaml")
    expose = set(["memory", "cores", "jvm_opts"])
    with open(base_file) as in_handle:
        config = yaml.load(in_handle)
    if os.path.exists(expose_file):
        with open(expose_file) as in_handle:
            expose_config = yaml.load(in_handle)
    else:
        expose_config = {"resources": {}}
    for pname, vals in config["resources"].iteritems():
        expose_vals = {}
        for k, v in vals.iteritems():
            if k in expose:
                expose_vals[k] = v
        if len(expose_vals) > 0 and pname not in expose_config["resources"]:
            expose_config["resources"][pname] = expose_vals
    with open(expose_file, "w") as out_handle:
        yaml.safe_dump(expose_config, out_handle, default_flow_style=False, allow_unicode=False)
    return expose_file

def _default_deploy_args(args):
    toolplus = {"data": {"bio_nextgen": []}}
    custom_add = collections.defaultdict(list)
    for x in args.toolplus:
        if not x.fname:
            for k, vs in toolplus.get(x.name, {}).iteritems():
                custom_add[k].extend(vs)
    return {"flavor": "ngs_pipeline_minimal",
            "custom_add": dict(custom_add),
            "vm_provider": "novm",
            "hostname": "localhost",
            "fabricrc_overrides": {"edition": "minimal",
                                   "use_sudo": args.sudo,
                                   "keep_isolated": args.isolate,
                                   "distribution": args.distribution or "__auto__",
                                   "dist_name": "__auto__"}}

def _update_conda_packages():
    """If installed in an anaconda directory, upgrade conda packages.
    """
    conda_bin = os.path.join(os.path.dirname(sys.executable), "conda")
    pkgs = ["biopython", "boto", "cython", "ipython", "lxml", "matplotlib",
            "nose", "numpy", "pandas", "patsy", "pycrypto", "pip", "pysam",
            "pyyaml", "pyzmq", "requests", "scipy", "setuptools", "sqlalchemy",
            "statsmodels", "tornado"]
    channels = ["-c", "https://conda.binstar.org/collections/chapmanb/bcbio"]
    if os.path.exists(conda_bin):
        subprocess.check_call([conda_bin, "install", "--yes", "numpy"])
        subprocess.check_call([conda_bin, "install", "--yes"] + channels + pkgs)

def _get_data_dir():
    base_dir = os.path.realpath(os.path.dirname(os.path.dirname(sys.executable)))
    if "anaconda" not in os.path.basename(base_dir) and "virtualenv" not in os.path.basename(base_dir):
        raise ValueError("Cannot update data for bcbio-nextgen not installed by installer.\n"
                         "bcbio-nextgen needs to be installed inside an anaconda environment \n"
                         "located in the same directory as `galaxy` `genomes` and `gemini_data` directories.")
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
    _upgrade_snpeff_data(s["fabricrc_overrides"]["galaxy_home"], args, remotes)
    if 'data' in set([x.name for x in args.toolplus]):
        gemini = os.path.join(os.path.dirname(sys.executable), "gemini")
        subprocess.check_call([gemini, "update", "--dataonly"])

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

def _upgrade_snpeff_data(galaxy_dir, args, remotes):
    """Install or upgrade snpEff databases, localized to reference directory.
    """
    for dbkey, ref_file in genome.get_builds(galaxy_dir):
        resource_file = os.path.join(os.path.dirname(ref_file), "%s-resources.yaml" % dbkey)
        if os.path.exists(resource_file):
            with open(resource_file) as in_handle:
                resources = yaml.load(in_handle)
            snpeff_db, snpeff_base_dir = effects.get_db({"genome_resources": resources,
                                                         "reference": {"fasta": {"base": ref_file}}})
            if snpeff_db:
                snpeff_db_dir = os.path.join(snpeff_base_dir, snpeff_db)
                if not os.path.exists(snpeff_db_dir):
                    print("Installing snpEff database %s in %s" % (snpeff_db, snpeff_base_dir))
                    tooldir = args.tooldir or get_defaults()["tooldir"]
                    config = {"resources": {"snpeff": {"jvm_opts": ["-Xms500m", "-Xmx1g"],
                                                       "dir": os.path.join(tooldir, "share", "java", "snpeff")}}}
                    raw_version = programs.java_versioner("snpeff", "snpEff",
                                                          stdout_flag="snpEff version SnpEff")(config)
                    snpeff_version = "".join([x for x in raw_version
                                              if x in set(string.digits + ".")]).replace(".", "_")
                    dl_url = remotes["snpeff_dl_url"].format(snpeff_ver=snpeff_version, genome=snpeff_db)
                    dl_file = os.path.basename(dl_url)
                    with utils.chdir(snpeff_base_dir):
                        subprocess.check_call(["wget", "-c", "-O", dl_file, dl_url])
                        subprocess.check_call(["unzip", dl_file])
                        os.remove(dl_file)
                    dl_dir = os.path.join(snpeff_base_dir, "data", snpeff_db)
                    os.rename(dl_dir, snpeff_db_dir)
                    os.rmdir(os.path.join(snpeff_base_dir, "data"))

def _get_biodata(base_file, args):
    with open(base_file) as in_handle:
        config = yaml.load(in_handle)
    config["install_liftover"] = False
    config["genome_indexes"] = args.aligners
    config["genomes"] = [g for g in config["genomes"] if g["dbkey"] in args.genomes]
    return config

def upgrade_thirdparty_tools(args, remotes):
    """Install and update third party tools used in the pipeline.

    Creates a manifest directory with installed programs on the system.
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
    manifest_dir = os.path.join(_get_data_dir(), "manifest")
    print("Creating manifest of installed packages in %s" % manifest_dir)
    cbl_manifest = __import__("cloudbio.manifest", fromlist=["manifest"])
    if os.path.exists(manifest_dir):
        for fname in os.listdir(manifest_dir):
            if not fname.startswith("toolplus"):
                os.remove(os.path.join(manifest_dir, fname))
    cbl_manifest.create(manifest_dir, args.tooldir)
    print("Installing additional tools")
    _install_toolplus(args, manifest_dir)

def _install_toolplus(args, manifest_dir):
    """Install additional tools we cannot distribute, updating local manifest.
    """
    toolplus_manifest = os.path.join(manifest_dir, "toolplus-packages.yaml")
    system_config = os.path.join(_get_data_dir(), "galaxy", "bcbio_system.yaml")
    toolplus_dir = os.path.join(_get_data_dir(), "toolplus")
    for tool in args.toolplus:
        if tool.name == "data":
            _install_gemini(args.tooldir, _get_data_dir(), args)
        elif tool.name in set(["gatk", "mutect"]):
            _install_gatk_jar(tool.name, tool.fname, toolplus_manifest, system_config, toolplus_dir)
        elif tool.name in set(["protected"]):  # back compatibility
            pass
        else:
            raise ValueError("Unexpected toolplus argument: %s %s" (tool.name, tool.fname))

def _install_gatk_jar(name, fname, manifest, system_config, toolplus_dir):
    """Install a jar for GATK or associated tools like MuTect.
    """
    if not fname.endswith(".jar"):
        raise ValueError("--toolplus argument for %s expects a jar file: %s" % (name, fname))
    if name == "gatk":
        version = broad.get_gatk_version(fname)
    elif name == "mutect":
        version = broad.get_mutect_version(fname)
    else:
        raise ValueError("Unexpected GATK input: %s" % name)
    store_dir = utils.safe_makedir(os.path.join(toolplus_dir, name, version))
    shutil.copyfile(fname, os.path.join(store_dir, os.path.basename(fname)))
    _update_system_file(system_config, name, {"dir": store_dir})
    _update_manifest(manifest, name, version)

def _update_manifest(manifest_file, name, version):
    """Update the toolplus manifest file with updated name and version
    """
    if os.path.exists(manifest_file):
        with open(manifest_file) as in_handle:
            manifest = yaml.load(in_handle)
    else:
        manifest = {}
    manifest[name] = {"name": name, "version": version}
    with open(manifest_file, "w") as out_handle:
        yaml.safe_dump(manifest, out_handle, default_flow_style=False, allow_unicode=False)

def _update_system_file(system_file, name, new_kvs):
    """Update the bcbio_system.yaml file with new resource information.
    """
    bak_file = system_file + ".bak%s" % datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    shutil.copyfile(system_file, bak_file)
    with open(system_file) as in_handle:
        config = yaml.load(in_handle)
    new_rs = {}
    for rname, r_kvs in config.get("resources", {}).iteritems():
        if rname == name:
            for k, v in new_kvs.iteritems():
                r_kvs[k] = v
        new_rs[rname] = r_kvs
    config["resources"] = new_rs
    with open(system_file, "w") as out_handle:
        yaml.safe_dump(config, out_handle, default_flow_style=False, allow_unicode=False)

def _install_gemini(tooldir, datadir, args):
    """Install gemini layered on top of bcbio-nextgen, sharing anaconda framework.
    """
    # check if we have an up to date version, upgrading if needed
    gemini = os.path.join(os.path.dirname(sys.executable), "gemini")
    if os.path.exists(gemini):
        vurl = "https://raw.github.com/arq5x/gemini/master/requirements.txt"
        r = requests.get(vurl)
        for line in r.text.split():
            if line.startswith("gemini=="):
                latest_version = line.split("==")[-1]
        cur_version = subprocess.check_output([gemini, "-v"], stderr=subprocess.STDOUT).strip().split()[-1]
        if LooseVersion(latest_version) > LooseVersion(cur_version):
            subprocess.check_call([gemini, "update"])
    # install from scratch inside existing Anaconda python
    else:
        url = "https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py"
        script = os.path.basename(url)
        subprocess.check_call(["wget", "-O", script, url, "--no-check-certificate"])
        cmd = [sys.executable, "-E", script, tooldir, datadir, "--notools", "--nodata", "--sharedpy"]
        if not args.sudo:
            cmd.append("--nosudo")
        subprocess.check_call(cmd)
        os.remove(script)

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
    if args.tooldir:
        cur_config["tooldir"] = args.tooldir
    cur_config["sudo"] = args.sudo
    cur_config["isolate"] = args.isolate
    for attr in ["genomes", "aligners"]:
        if not cur_config.get(attr):
            cur_config[attr] = []
        for x in getattr(args, attr):
            if x not in cur_config[attr]:
                cur_config[attr].append(x)
    # toolplus -- save non-filename inputs
    attr = "toolplus"
    if not cur_config.get(attr):
        cur_config[attr] = []
    for x in getattr(args, attr):
        if not x.fname:
            if x.name not in cur_config[attr]:
                cur_config[attr].append(x.name)
    with open(install_config, "w") as out_handle:
        yaml.safe_dump(cur_config, out_handle, default_flow_style=False, allow_unicode=False)

def add_install_defaults(args):
    """Add any saved installation defaults to the upgrade.
    """
    install_config = _get_install_config()
    if install_config is None or not utils.file_exists(install_config):
        return args
    with open(install_config) as in_handle:
        default_args = yaml.load(in_handle)
    if args.tools and args.tooldir is None:
        if "tooldir" in default_args:
            args.tooldir = str(default_args["tooldir"])
        else:
            raise ValueError("Default tool directory not yet saved in config defaults. "
                             "Specify the '--tooldir=/path/to/tools' to upgrade tools. "
                             "After a successful upgrade, the '--tools' parameter will "
                             "work for future upgrades.")
    for attr in ["genomes", "aligners", "toolplus"]:
        for x in default_args.get(attr, []):
            x = Tool(x, None) if attr == "toolplus" else str(x)
            new_val = getattr(args, attr)
            if x not in getattr(args, attr):
                new_val.append(x)
            setattr(args, attr, new_val)
    if "sudo" in default_args and not args.sudo is False:
        args.sudo = default_args["sudo"]
    if "isolate" in default_args and not args.isolate is True:
        args.isolate = default_args["isolate"]
    return args

def get_defaults():
    install_config = _get_install_config()
    if install_config is None or not utils.file_exists(install_config):
        return {}
    with open(install_config) as in_handle:
        return yaml.load(in_handle)

def _check_toolplus(x):
    """Parse options for adding non-standard/commercial tools like GATK and MuTecT.
    """
    std_choices = set(["data"])
    if x in std_choices:
        return Tool(x, None)
    elif "=" in x and len(x.split("=")) == 2:
        name, fname = x.split("=")
        fname = os.path.normpath(os.path.realpath(fname))
        if not os.path.exists(fname):
            raise argparse.ArgumentTypeError("Unexpected --toolplus argument for %s. File does not exist: %s"
                                             % (name, fname))
        return Tool(name, fname)
    else:
        raise argparse.ArgumentTypeError("Unexpected --toolplus argument. Expect toolname=filename.")

def add_subparser(subparsers):
    parser = subparsers.add_parser("upgrade", help="Install or upgrade bcbio-nextgen")
    parser.add_argument("--tooldir",
                        help="Directory to install 3rd party software tools. Leave unspecified for no tools",
                        type=lambda x: (os.path.abspath(os.path.expanduser(x))), default=None)
    parser.add_argument("--tools",
                        help="Boolean argument specifying upgrade of tools. Uses previously saved install directory",
                        action="store_true", default=False)
    parser.add_argument("-u", "--upgrade", help="Code version to upgrade",
                        choices=["stable", "development", "system", "skip"], default="skip")
    parser.add_argument("--toolplus", help="Specify additional tool categories to install",
                        action="append", default=[], type=_check_toolplus)
    parser.add_argument("--genomes", help="Genomes to download",
                        action="append", default=["GRCh37"],
                        choices=["GRCh37", "hg19", "mm10", "mm9", "rn5", "canFam3"])
    parser.add_argument("--aligners", help="Aligner indexes to download",
                        action="append", default=["bwa", "bowtie2"],
                        choices=["bowtie", "bowtie2", "bwa", "novoalign", "star", "ucsc"])
    parser.add_argument("--data", help="Upgrade data dependencies",
                        dest="install_data", action="store_true", default=False)
    parser.add_argument("--nosudo", help="Specify we cannot use sudo for commands",
                        dest="sudo", action="store_false", default=True)
    parser.add_argument("--isolate", help="Created an isolated installation without PATH updates",
                        dest="isolate", action="store_true", default=False)
    parser.add_argument("--distribution", help="Operating system distribution",
                        default="",
                        choices=["ubuntu", "debian", "centos", "scientificlinux", "macosx"])
    return parser

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
