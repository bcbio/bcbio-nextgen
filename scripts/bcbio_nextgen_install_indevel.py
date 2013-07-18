#!/usr/bin/env python
"""Automatically install required tools and data to run bcbio-nextgen pipelines.

This automates the steps required for installation and setup to make it
easier to get started with bcbio-nextgen. The defaults provide data files
for human variant calling.

Requires: git

# XXX In progress development script using new approaches to improve automation and avoid
installation issues. Not yet fully functional.
"""
import argparse
import contextlib
import os
import shutil
import subprocess
import sys

remotes = {"requirements":
           "https://raw.github.com/chapmanb/bcbio-nextgen/master/requirements.txt",
           "anaconda":
           "http://repo.continuum.io/pkgs/miniconda/Miniconda-1.6.2-%s-x86_64.sh"}

def main(args, sys_argv):
    check_dependencies()
    with bcbio_tmpdir():
        setup_data_dir(args)
        print("Installing isolated base python installation")
        anaconda = install_anaconda_python(args, remotes)
        print("Installing bcbio-nextgen")
        install_conda_pkgs(anaconda)
        bcbio = bootstrap_bcbionextgen(anaconda, args, remotes)
        subprocess.check_call([bcbio["bcbio_nextgen.py"], "upgrade"] + sys_argv[1:])

def bootstrap_bcbionextgen(anaconda, args, remotes):
    """Install bcbio-nextgen to bootstrap rest of installation process.
    """
    subprocess.check_call([anaconda["pip"], "install", "fabric"])
    subprocess.check_call([anaconda["pip"], "install", "-r", remotes["requirements"]])
    out = {}
    for script in ["bcbio_nextgen.py"]:
        final_script = os.path.join(args.tooldir, "bin", script)
        ve_script = os.path.join(anaconda["dir"], "bin", script)
        if not os.path.exists(final_script):
            sudo_cmd = ["sudo"] if args.sudo else []
            subprocess.check_call(sudo_cmd + ["mkdir", "-p", os.path.dirname(final_script)])
            cmd = ["ln", "-s", ve_script, final_script]
            subprocess.check_call(sudo_cmd + cmd)
        out[script] = ve_script
    return out

def install_conda_pkgs(anaconda):
    pkgs = ["biopython", "boto", "cython", "numpy", "pycrypto", "pip", "pysam", "pyyaml"]
    subprocess.check_call([anaconda["conda"], "install", "--yes"] + pkgs)

def install_anaconda_python(args, remotes):
    """Provide isolated installation of Anaconda python for running bcbio-nextgen.
    http://docs.continuum.io/anaconda/index.html
    """
    anaconda_dir = os.path.join(args.datadir, "anaconda")
    if not os.path.exists(anaconda_dir):
        url = remotes["anaconda"] % ("MacOSX" if args.distribution.lower() == "macosx" else "Linux")
        if not os.path.exists(os.path.basename(url)):
            subprocess.check_call(["wget", url])
        subprocess.check_call("echo -e '\nyes\n%s\nno\n' | bash %s" %
                              (anaconda_dir, os.path.basename(url)), shell=True)
    bindir = os.path.join(anaconda_dir, "bin")
    return {"conda": os.path.join(bindir, "conda"),
            "pip": os.path.join(bindir, "pip"),
            "dir": anaconda_dir}

def setup_data_dir(args):
    if not os.path.exists(args.datadir):
        cmd = ["mkdir", "-p", args.datadir]
        if args.sudo:
            cmd.insert(0, "sudo")
        subprocess.check_call(cmd)
    if args.sudo:
        subprocess.check_call(["sudo", "chown", "-R", os.environ["USER"], args.datadir])

def check_dependencies():
    """Ensure required tools for installation are present.
    """
    print("Checking required dependencies")
    try:
        subprocess.check_call(["git", "--version"])
    except OSError:
        raise OSError("bcbio-nextgen installer requires Git (http://git-scm.com/)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Automatic installation for bcbio-nextgen pipelines")
    parser.add_argument("datadir", help="Directory to install genome data",
                        type=os.path.abspath)
    parser.add_argument("--distribution", help="Operating system distribution",
                        default="ubuntu",
                        choices=["ubuntu", "debian", "centos", "scientificlinux"])
    parser.add_argument("--tooldir",
                        help="Directory to install 3rd party software tools. Leave unspecified for no tools",
                        type=os.path.abspath, default=None)
    parser.add_argument("--genomes", help="Genomes to download",
                        action="append", default=["hg19", "GRCh37"])
    parser.add_argument("--aligners", help="Aligner indexes to download",
                        action="append", default=["bwa"])
    parser.add_argument("--nosudo", help="Specify we cannot use sudo for commands",
                        dest="sudo", action="store_false", default=True)
    parser.add_argument("--nodata", help="Do not install data dependencies",
                        dest="install_data", action="store_false", default=True)
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        main(parser.parse_args(), sys.argv[1:])
