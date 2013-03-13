#!/usr/bin/env python
"""Automatically install required tools and data to run bcbio-nextgen pipelines.

This automates the steps required for installation and setup to make it
easier to get started with bcbio-nextgen. The defaults provide data files
for human variant calling.

Requires: PyYAML, fabric
"""
import argparse
import contextlib
import datetime
import os
import urllib2
import shutil
import subprocess
import sys

import yaml

remotes = {"system_config":
           "https://raw.github.com/chapmanb/bcbio-nextgen/master/config/bcbio_system.yaml",
           "requirements":
           "https://raw.github.com/chapmanb/bcbio-nextgen/master/requirements.txt",
           "virtualenv":
           "https://raw.github.com/pypa/virtualenv/master/virtualenv.py"}

def main(args):
    check_dependencies()
    work_dir = os.path.join(os.getcwd(), "tmpbcbio-install")
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    os.chdir(work_dir)
    cbl = get_cloudbiolinux()
    fabricrc = write_fabricrc(cbl["fabricrc"], args.tooldir, args.datadir,
                              args.distribution, args.sudo)
    biodata = write_biodata(cbl["biodata"], args.genomes, args.aligners)
    if args.install_tools:
        print "Installing tools..."
        install_tools(cbl["tool_fabfile"], fabricrc)
    print "Installing data..."
    if args.install_data:
        install_data(cbl["data_fabfile"], fabricrc, biodata)
    print "Installing bcbio-nextgen..."
    install_bcbio_nextgen(remotes["requirements"], args.datadir, args.tooldir,
                          args.sudo)
    system_config = write_system_config(remotes["system_config"], args.datadir,
                                        args.tooldir)
    print "Finished: bcbio-nextgen, tools and data installed"
    print " Ready to use system configuration at:\n  %s" % system_config
    print " Tools installed in:\n  %s" % args.tooldir
    print " Genome data installed in:\n  %s" % args.datadir
    shutil.rmtree(work_dir)

def install_bcbio_nextgen(requirements, datadir, tooldir, use_sudo):
    """Install a virtualenv containing bcbio_nextgen depdencies.
    """
    virtualenv_dir = os.path.join(datadir, "bcbio-nextgen-virtualenv")
    if not os.path.exists(virtualenv_dir):
        subprocess.check_call(["wget", remotes["virtualenv"]])
        subprocess.check_call(["python", "virtualenv.py", "--no-site-packages",
                               "--distribute", virtualenv_dir])
        os.remove("virtualenv.py")
    pip_cmd = os.path.join(virtualenv_dir, "bin", "pip")
    subprocess.check_call([pip_cmd, "install", "--upgrade", "distribute"])
    subprocess.check_call([pip_cmd, "install", "-r", requirements])
    for script in ["bcbio_nextgen.py", "bam_to_wiggle.py"]:
        final_script = os.path.join(tooldir, "bin", script)
        ve_script = os.path.join(virtualenv_dir, "bin", script)
        if not os.path.exists(final_script):
            sudo_cmd = ["sudo"] if use_sudo else []
            subprocess.check_call(sudo_cmd + ["mkdir", "-p", os.path.dirname(final_script)])
            cmd = ["ln", "-s", ve_script, final_script]
            subprocess.check_call(sudo_cmd + cmd)

def install_tools(fabfile, fabricrc):
    subprocess.check_call(["fab", "-f", fabfile, "-H", "localhost",
                           "-c", fabricrc,
                           "install_biolinux:flavor=ngs_pipeline"])

def install_data(fabfile, fabricrc, biodata):
    subprocess.check_call(["fab", "-f", fabfile, "-H", "localhost",
                           "-c", fabricrc, "install_data_s3:%s" % biodata])

def write_system_config(base_url, datadir, tooldir):
    java_basedir = os.path.join(tooldir, "share", "java")
    out_file = os.path.join(datadir, "galaxy", os.path.basename(base_url))
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))
    if os.path.exists(out_file):
        bak_file = out_file + ".bak%s" % (datetime.datetime.now().strftime("%Y%M%d_%H%M"))
    to_rewrite = ("gatk", "picard", "snpEff", "bcbio_variation")
    with contextlib.closing(urllib2.urlopen(base_url)) as in_handle:
        with open(out_file, "w") as out_handle:
            in_prog = None
            for line in in_handle:
                if line.strip().startswith(to_rewrite):
                    in_prog = line.split(":")[0].strip()
                elif line.strip().startswith("dir:") and in_prog:
                    line = "%s: %s\n" % (line.split(":")[0],
                                         os.path.join(java_basedir, in_prog.lower()))
                    in_prog = None
                elif line.startswith("galaxy"):
                    line = "# %s" % line
                out_handle.write(line)
    return out_file

def write_biodata(base_file, genomes, aligners):
    out_file = os.path.join(os.getcwd(), os.path.basename(base_file))
    with open(base_file) as in_handle:
        config = yaml.load(in_handle)
    config["install_liftover"] = False
    config["genome_indexes"] = aligners
    config["genomes"] = [g for g in config["genomes"] if g["dbkey"] in genomes]
    with open(out_file, "w") as out_handle:
        yaml.dump(config, out_handle, allow_unicode=False, default_flow_style=False)
    return out_file

def write_fabricrc(base_file, tooldir, datadir, distribution, use_sudo):
    out_file = os.path.join(os.getcwd(), os.path.basename(base_file))
    with open(base_file) as in_handle:
        with open(out_file, "w") as out_handle:
            for line in in_handle:
                if line.startswith("system_install"):
                    line = "system_install = %s\n" % tooldir
                elif line.startswith("local_install"):
                    line = "local_install = %s/install\n" % tooldir
                elif line.startswith("data_files"):
                    line = "data_files = %s\n" % datadir
                elif line.startswith("distribution"):
                    line = "distribution = %s\n" % distribution
                elif line.startswith("use_sudo"):
                    line = "use_sudo = %s\n" % use_sudo
                elif line.startswith("edition"):
                    line = "edition = minimal\n"
                elif line.startswith("#galaxy_home"):
                    line = "galaxy_home = %s\n" % os.path.join(datadir, "galaxy")
                out_handle.write(line)
    return out_file

def get_cloudbiolinux():
    base_dir = os.path.join(os.getcwd(), "cloudbiolinux")
    if not os.path.exists(base_dir):
        subprocess.check_call(["git", "clone",
                               "git://github.com/chapmanb/cloudbiolinux.git"])
    return {"fabricrc": os.path.join(base_dir, "config", "fabricrc.txt"),
            "biodata": os.path.join(base_dir, "config", "biodata.yaml"),
            "tool_fabfile": os.path.join(base_dir, "fabfile.py"),
            "data_fabfile": os.path.join(base_dir, "data_fabfile.py")}

def check_dependencies():
    """Ensure required tools for installation are present.
    """
    print "Checking required dependencies"
    try:
        subprocess.check_call(["fab", "--version"])
    except OSError:
        raise OSError("bcbio-nextgen installer requires fabric (http://fabfile.org)")
    try:
        subprocess.check_call(["git", "--version"])
    except OSError:
        raise OSError("bcbio-nextgen installer requires Git (http://git-scm.com/)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Automatic installation for bcbio-nextgen pipelines")
    parser.add_argument("tooldir", help="Directory to install 3rd party software tools",
                        type=os.path.abspath)
    parser.add_argument("datadir", help="Directory to install genome data",
                        type=os.path.abspath)
    parser.add_argument("--distribution", help="Operating system distribution",
                        default="ubuntu",
                        choices=["ubuntu", "debian", "centos", "scientificlinux"])
    parser.add_argument("--genomes", help="Genomes to download",
                        action="append", default=["hg19", "GRCh37"])
    parser.add_argument("--aligners", help="Aligner indexes to download",
                        action="append", default=["bwa", "bowtie2", "novoalign", "ucsc"])
    parser.add_argument("--nosudo", help="Specify we cannot use sudo for commands",
                        dest="sudo", action="store_false", default=True)
    parser.add_argument("--notools", help="Do not install tool dependencies",
                        dest="install_tools", action="store_false", default=True)
    parser.add_argument("--nodata", help="Do not install data dependencies",
                        dest="install_data", action="store_false", default=True)
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        main(parser.parse_args())
