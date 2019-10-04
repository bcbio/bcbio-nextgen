#!/usr/bin/env python
"""Automatically install required tools and data to run bcbio-nextgen pipelines.

This automates the steps required for installation and setup to make it
easier to get started with bcbio-nextgen. The defaults provide data files
for human variant calling.

Requires: git, wget, bgzip2, Python 3.x, Python 2.7 or argparse + Python 2.6 and earlier
"""
from __future__ import print_function
import collections
import contextlib
import datetime
import os
import platform
import shutil
import subprocess
import sys
try:
    import urllib2 as urllib_request
except ImportError:
    import urllib.request as urllib_request

REMOTES = {
    "requirements": "https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/requirements-conda.txt",
    "gitrepo": "https://github.com/bcbio/bcbio-nextgen.git",
    "system_config": "https://raw.github.com/bcbio/bcbio-nextgen/master/config/bcbio_system.yaml",
    "anaconda": "https://repo.continuum.io/miniconda/Miniconda3-latest-%s-x86_64.sh"}
TARGETPY = "python=3.6"

def main(args, sys_argv):
    check_arguments(args)
    check_dependencies()
    with bcbio_tmpdir():
        setup_data_dir(args)
        print("Installing isolated base python installation")
        anaconda = install_anaconda_python(args)
        print("Installing bcbio-nextgen")
        bcbio = install_conda_pkgs(anaconda, args)
        bootstrap_bcbionextgen(anaconda, args)
    print("Installing data and third party dependencies")
    system_config = write_system_config(REMOTES["system_config"], args.datadir,
                                        args.tooldir)
    setup_manifest(args.datadir)
    subprocess.check_call([bcbio, "upgrade"] + _clean_args(sys_argv, args))
    print("Finished: bcbio-nextgen, tools and data installed")
    print(" Genome data installed in:\n  %s" % args.datadir)
    if args.tooldir:
        print(" Tools installed in:\n  %s" % args.tooldir)
    print(" Ready to use system configuration at:\n  %s" % system_config)
    print(" Edit configuration file as needed to match your machine or cluster")

def _clean_args(sys_argv, args):
    """Remove data directory from arguments to pass to upgrade function.
    """
    base = [x for x in sys_argv if
            x.startswith("-") or not args.datadir == os.path.abspath(os.path.expanduser(x))]
    # Remove installer only options we don't pass on
    base = [x for x in base if x not in set(["--minimize-disk"])]
    if "--nodata" in base:
        base.remove("--nodata")
    else:
        base.append("--data")
    return base

def bootstrap_bcbionextgen(anaconda, args):
    if args.upgrade == "development":
        git_tag = "@%s" % args.revision if args.revision != "master" else ""
        subprocess.check_call([anaconda["pip"], "install", "--upgrade", "--no-deps",
                               "git+%s%s#egg=bcbio-nextgen" % (REMOTES["gitrepo"], git_tag)])

def _get_conda_channels(conda_bin):
    """Retrieve default conda channels, checking if they are pre-specified in config.

    This allows users to override defaults with specific mirrors in their .condarc
    """
    channels = ["bioconda", "conda-forge"]
    out = []
    try:
        import yaml
        config = yaml.safe_load(subprocess.check_output([conda_bin, "config", "--show"]))
    except ImportError:
        config = {}
    for c in channels:
        present = False
        for orig_c in config.get("channels") or []:
            if orig_c.endswith((c, "%s/" % c)):
                present = True
                break
        if not present:
            out += ["-c", c]
    return out

def install_conda_pkgs(anaconda, args):
    env = dict(os.environ)
    # Try to avoid user specific pkgs and envs directories
    # https://github.com/conda/conda/issues/6748
    env["CONDA_PKGS_DIRS"] = os.path.join(anaconda["dir"], "pkgs")
    env["CONDA_ENVS_DIRS"] = os.path.join(anaconda["dir"], "envs")
    if not os.path.exists(os.path.basename(REMOTES["requirements"])):
        subprocess.check_call(["wget", "--no-check-certificate", REMOTES["requirements"]])
    if args.minimize_disk:
        subprocess.check_call([anaconda["conda"], "install", "--yes", "nomkl"], env=env)
    channels = _get_conda_channels(anaconda["conda"])
    subprocess.check_call([anaconda["conda"], "install", "--yes"] + channels +
                          ["--only-deps", "bcbio-nextgen", TARGETPY], env=env)
    subprocess.check_call([anaconda["conda"], "install", "--yes"] + channels +
                          ["--file", os.path.basename(REMOTES["requirements"]), TARGETPY], env=env)
    return os.path.join(anaconda["dir"], "bin", "bcbio_nextgen.py")

def _guess_distribution():
    """Simple approach to identify if we are on a MacOSX or Linux system for Anaconda.
    """
    if platform.mac_ver()[0]:
        return "macosx"
    else:
        return "linux"

def install_anaconda_python(args):
    """Provide isolated installation of Anaconda python for running bcbio-nextgen.
    http://docs.continuum.io/anaconda/index.html
    """
    anaconda_dir = os.path.join(args.datadir, "anaconda")
    bindir = os.path.join(anaconda_dir, "bin")
    conda = os.path.join(bindir, "conda")
    if not os.path.exists(anaconda_dir) or not os.path.exists(conda):
        if os.path.exists(anaconda_dir):
            shutil.rmtree(anaconda_dir)
        dist = args.distribution if args.distribution else _guess_distribution()
        url = REMOTES["anaconda"] % ("MacOSX" if dist.lower() == "macosx" else "Linux")
        if not os.path.exists(os.path.basename(url)):
            subprocess.check_call(["wget", "--progress=dot:mega", "--no-check-certificate", url])
        subprocess.check_call("bash %s -b -p %s" %
                              (os.path.basename(url), anaconda_dir), shell=True)
    return {"conda": conda,
            "pip": os.path.join(bindir, "pip"),
            "dir": anaconda_dir}

def setup_manifest(datadir):
    """Create barebones manifest to be filled in during update
    """
    manifest_dir = os.path.join(datadir, "manifest")
    if not os.path.exists(manifest_dir):
        os.makedirs(manifest_dir)

def write_system_config(base_url, datadir, tooldir):
    """Write a bcbio_system.yaml configuration file with tool information.
    """
    out_file = os.path.join(datadir, "galaxy", os.path.basename(base_url))
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))
    if os.path.exists(out_file):
        # if no tool directory and exists, do not overwrite
        if tooldir is None:
            return out_file
        else:
            bak_file = out_file + ".bak%s" % (datetime.datetime.now().strftime("%Y%M%d_%H%M"))
            shutil.copy(out_file, bak_file)
    if tooldir:
        java_basedir = os.path.join(tooldir, "share", "java")
    rewrite_ignore = ("log",)
    with contextlib.closing(urllib_request.urlopen(base_url)) as in_handle:
        with open(out_file, "w") as out_handle:
            in_resources = False
            in_prog = None
            for line in (l.decode("utf-8") for l in in_handle):
                if line[0] != " ":
                    in_resources = line.startswith("resources")
                    in_prog = None
                elif (in_resources and line[:2] == "  " and line[2] != " "
                      and not line.strip().startswith(rewrite_ignore)):
                    in_prog = line.split(":")[0].strip()
                # Update java directories to point to install directory, avoid special cases
                elif line.strip().startswith("dir:") and in_prog and in_prog not in ["log", "tmp"]:
                    final_dir = os.path.basename(line.split()[-1])
                    if tooldir:
                        line = "%s: %s\n" % (line.split(":")[0],
                                             os.path.join(java_basedir, final_dir))
                    in_prog = None
                elif line.startswith("galaxy"):
                    line = "# %s" % line
                out_handle.write(line)
    return out_file

def setup_data_dir(args):
    if not os.path.exists(args.datadir):
        cmd = ["mkdir", "-p", args.datadir]
        subprocess.check_call(cmd)

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

def check_arguments(args):
    """Ensure argruments are consistent and correct.
    """
    if args.toolplus and not args.tooldir:
        raise argparse.ArgumentTypeError("Cannot specify --toolplus without --tooldir")

def check_dependencies():
    """Ensure required tools for installation are present.
    """
    print("Checking required dependencies")
    for dep, msg in [(["git", "--version"], "Git (http://git-scm.com/)"),
                     (["wget", "--version"], "wget"),
                     (["bzip2", "-h"], "bzip2")]:
        try:
            p = subprocess.Popen(dep, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
            out, code = p.communicate()
        except OSError:
            out = "Executable not found"
            code = 127
        if code == 127:
            raise OSError("bcbio-nextgen installer requires %s\n%s" % (msg, out))

def _check_toolplus(x):
    """Parse options for adding non-standard/commercial tools like GATK and MuTecT.
    """
    import argparse
    Tool = collections.namedtuple("Tool", ["name", "fname"])
    std_choices = set(["data", "dbnsfp", "ericscript"])
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

if __name__ == "__main__":
    try:
        import argparse
    except ImportError:
        raise ImportError("bcbio-nextgen installer requires `argparse`, included in Python 2.7.\n"
                          "Install for earlier versions with `pip install argparse` or "
                          "`easy_install argparse`.")
    parser = argparse.ArgumentParser(
        description="Automatic installation for bcbio-nextgen pipelines")
    parser.add_argument("datadir", help="Directory to install genome data",
                        type=lambda x: (os.path.abspath(os.path.expanduser(x))))
    parser.add_argument("--cores", default=1,
                        help="Number of cores to use if local indexing is necessary.")
    parser.add_argument("--tooldir",
                        help="Directory to install 3rd party software tools. Leave unspecified for no tools",
                        type=lambda x: (os.path.abspath(os.path.expanduser(x))), default=None)
    parser.add_argument("--toolplus", help="Specify additional tool categories to install",
                        action="append", default=[], type=_check_toolplus)
    parser.add_argument("--datatarget", help="Data to install. Allows customization or install of extra data.",
                        action="append", default=[],
                        choices=["variation", "rnaseq", "smallrna", "gemini", "vep", "dbnsfp",
                                 "battenberg", "kraken", "ericscript", "gnomad"])
    parser.add_argument("--genomes", help="Genomes to download",
                        action="append", default=[],
                        choices=["GRCh37", "hg19", "hg38", "hg38-noalt", "mm10", "mm9", "rn6", "rn5",
                                 "canFam3", "dm3", "galGal4", "phix", "pseudomonas_aeruginosa_ucbpp_pa14",
                                 "sacCer3", "TAIR10", "WBcel235", "xenTro3", "GRCz10", "GRCz11",
                                 "Sscrofa11.1", "BDGP6"])
    parser.add_argument("--aligners", help="Aligner indexes to download",
                        action="append", default=[],
                        choices=["bbmap", "bowtie", "bowtie2", "bwa", "minimap2", "novoalign", "rtg", "snap",
                                 "star", "ucsc", "hisat2"])
    parser.add_argument("--nodata", help="Do not install data dependencies",
                        dest="install_data", action="store_false", default=True)
    parser.add_argument("--isolate", help="Created an isolated installation without PATH updates",
                        dest="isolate", action="store_true", default=False)
    parser.add_argument("--minimize-disk", help="Try to minimize disk usage (no MKL extensions)",
                        dest="minimize_disk", action="store_true", default=False)
    parser.add_argument("-u", "--upgrade", help="Code version to install",
                        choices=["stable", "development"], default="stable")
    parser.add_argument("--revision", help="Specify a git commit hash or tag to install", default="master")
    parser.add_argument("--distribution", help="Operating system distribution",
                        default="",
                        choices=["ubuntu", "debian", "centos", "scientificlinux", "macosx"])
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        main(parser.parse_args(), sys.argv[1:])
