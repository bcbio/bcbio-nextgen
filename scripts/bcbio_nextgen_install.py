#!/usr/bin/env python
"""Automatically install required tools and data to run bcbio-nextgen pipelines.

This automates the steps required for installation and setup to make it
easier to get started with bcbio-nextgen. The defaults provide data files
for human variant calling.

Requires: git
"""
import argparse
import contextlib
import datetime
import os
import platform
import shutil
import subprocess
import sys
import urllib2

remotes = {"requirements":
           "https://raw.github.com/chapmanb/bcbio-nextgen/master/requirements.txt",
           "system_config":
           "https://raw.github.com/chapmanb/bcbio-nextgen/master/config/bcbio_system.yaml",
           "anaconda":
           "http://repo.continuum.io/miniconda/Miniconda-1.9.1-%s-x86_64.sh"}

def main(args, sys_argv):
    check_dependencies()
    with bcbio_tmpdir():
        setup_data_dir(args)
        print("Installing isolated base python installation")
        anaconda = install_anaconda_python(args, remotes)
        print("Installing bcbio-nextgen")
        install_conda_pkgs(anaconda)
        bcbio = bootstrap_bcbionextgen(anaconda, args, remotes)
    print("Installing data and third party dependencies")
    subprocess.check_call([bcbio["bcbio_nextgen.py"], "upgrade"] + _clean_args(sys_argv, args))
    system_config = write_system_config(remotes["system_config"], args.datadir,
                                        args.tooldir)
    print("Finished: bcbio-nextgen, tools and data installed")
    print(" Genome data installed in:\n  %s" % args.datadir)
    if args.tooldir:
        print(" Tools installed in:\n  %s" % args.tooldir)
    print(" Ready to use system configuration at:\n  %s" % system_config)
    print(" Edit configuration file as needed to match your machine or cluster")

def _clean_args(sys_argv, args):
    """Remove data directory from arguments to pass to upgrade function.
    """
    return [x for x in sys_argv if
            x.startswith("--") or not args.datadir == os.path.abspath(os.path.expanduser(x))]

def bootstrap_bcbionextgen(anaconda, args, remotes):
    """Install bcbio-nextgen to bootstrap rest of installation process.
    """
    subprocess.check_call([anaconda["pip"], "install", "fabric"])
    subprocess.check_call([anaconda["pip"], "install", "-r", remotes["requirements"]])
    out = {}
    for script in ["bcbio_nextgen.py"]:
        ve_script = os.path.join(anaconda["dir"], "bin", script)
        if args.tooldir:
            final_script = os.path.join(args.tooldir, "bin", script)
            sudo_cmd = ["sudo"] if args.sudo else []
            subprocess.check_call(sudo_cmd + ["mkdir", "-p", os.path.dirname(final_script)])
            if os.path.lexists(final_script):
                cmd = ["rm", "-f", final_script]
                subprocess.check_call(sudo_cmd + cmd)
            cmd = ["ln", "-s", ve_script, final_script]
            subprocess.check_call(sudo_cmd + cmd)
        out[script] = ve_script
    return out

def install_conda_pkgs(anaconda):
    pkgs = ["biopython", "boto", "cython", "distribute", "ipython", "nose", "numpy",
            "pycrypto", "pip", "pysam", "pyyaml", "pyzmq", "requests"]
    subprocess.check_call([anaconda["conda"], "install", "--yes"] + pkgs)

def _guess_distribution():
    """Simple approach to identify if we are on a MacOSX or Linux system for Anaconda.
    """
    if platform.mac_ver()[0]:
        return "macosx"
    else:
        return "linux"

def install_anaconda_python(args, remotes):
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
        url = remotes["anaconda"] % ("MacOSX" if dist.lower() == "macosx" else "Linux")
        if not os.path.exists(os.path.basename(url)):
            subprocess.check_call(["wget", url])
        subprocess.check_call("echo -e '\nyes\n%s\nno\n' | bash %s" %
                              (anaconda_dir, os.path.basename(url)), shell=True)
    return {"conda": conda,
            "pip": os.path.join(bindir, "pip"),
            "dir": anaconda_dir}

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
    to_rewrite = ("gatk", "picard", "snpEff", "bcbio_variation")
    with contextlib.closing(urllib2.urlopen(base_url)) as in_handle:
        with open(out_file, "w") as out_handle:
            in_prog = None
            for line in in_handle:
                if line.strip().startswith(to_rewrite):
                    in_prog = line.split(":")[0].strip()
                elif line.strip().startswith("dir:") and in_prog:
                    if tooldir:
                        line = "%s: %s\n" % (line.split(":")[0],
                                             os.path.join(java_basedir, in_prog.lower()))
                    in_prog = None
                elif line.startswith("galaxy"):
                    line = "# %s" % line
                out_handle.write(line)
    return out_file

def setup_data_dir(args):
    if not os.path.exists(args.datadir):
        cmd = ["mkdir", "-p", args.datadir]
        if args.sudo:
            cmd.insert(0, "sudo")
        subprocess.check_call(cmd)
    if args.sudo:
        subprocess.check_call(["sudo", "chown", "-R", os.environ["USER"], args.datadir])

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
                        type=lambda x: (os.path.abspath(os.path.expanduser(x))))
    parser.add_argument("--tooldir",
                        help="Directory to install 3rd party software tools. Leave unspecified for no tools",
                        type=lambda x: (os.path.abspath(os.path.expanduser(x))), default=None)
    parser.add_argument("--tooldist",
                        help="Type of tool distribution to install. Defaults to a minimum install.",
                        default="minimal",
                        choices=["minimal", "full"])
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
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        main(parser.parse_args(), sys.argv[1:])
