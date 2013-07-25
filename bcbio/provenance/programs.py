"""Identify program versions used for analysis, reporting in structured table.

Catalogs the full list of programs used in analysis, enabling reproduction of
results and tracking of provenance in output files.
"""
import os
import contextlib
import subprocess

from bcbio import broad, utils
from bcbio.pipeline import config_utils, version

_cl_progs = [{"cmd": "bamtools", "args": "--version", "stdout_flag": "bamtools"},
             {"cmd": "bedtools", "args": "--version", "stdout_flag": "bedtools"},
             {"cmd": "bowtie2", "args": "--version", "stdout_flag": "bowtie2-align"},
             {"cmd": "bwa", "stdout_flag": "Version:"},
             {"cmd": "fastqc", "args": "--version", "stdout_flag": "FastQC"},
             {"cmd": "freebayes", "stdout_flag": "version:"},
             {"cmd": "gemini", "args": "--version", "stdout_flag": "gemini"},
             {"cmd": "novosort", "paren_flag": "novosort"},
             {"cmd": "novoalign", "stdout_flag": "Novoalign"},
             {"cmd": "samtools", "stdout_flag": "Version"}]
# TODO: ogap, bamleftalign

def _broad_versioner(type):
    def get_version(config):
        runner = broad.runner_from_config(config)
        if type == "gatk":
            return runner.get_gatk_version()
        elif type == "picard":
            return runner.get_picard_version("ViewSam")
        else:
            raise NotImplementedError(type)
    return get_version

def jar_versioner(program_name, jar_name):
    """Retrieve version information based on jar file.
    """
    def get_version(config):
        try:
            pdir = config_utils.get_program(program_name, config, "dir")
        # not configured
        except ValueError:
            return ""
        jar = os.path.basename(config_utils.get_jar(jar_name, pdir))
        for to_remove in [jar_name, ".jar", "-standalone"]:
            jar = jar.replace(to_remove, "")
        if jar.startswith(("-", ".")):
            jar = jar[1:]
        return jar
    return get_version

_alt_progs = [{"name": "gatk", "version_fn": _broad_versioner("gatk")},
              {"name": "picard", "version_fn": _broad_versioner("picard")},
              {"name": "bcbio.variation",
               "version_fn": jar_versioner("bcbio_variation", "bcbio.variation")},
              {"name": "varscan",
               "version_fn": jar_versioner("varscan", "VarScan")}]
# TODO: cortex_var

def _parse_from_stdoutflag(stdout, x):
    for line in stdout:
        if line.find(x) >= 0:
            return line.split()[-1]
    return ""

def _parse_from_parenflag(stdout, x):
    for line in stdout:
        if line.find(x) >= 0:
            return line.split("(")[-1].split(")")[0]
    return ""

def _get_cl_version(p, config):
    """Retrieve version of a single commandline program.
    """
    try:
        prog = config_utils.get_program(p["cmd"], config)
    except config_utils.CmdNotFound:
        return "NA"
    args = p.get("args", "")

    cmd = "{prog} {args}"
    subp = subprocess.Popen(cmd.format(**locals()), stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            shell=True)
    with contextlib.closing(subp.stdout) as stdout:
        if p.get("stdout_flag"):
            v = _parse_from_stdoutflag(stdout, p["stdout_flag"])
        elif p.get("paren_flag"):
            v = _parse_from_parenflag(stdout, p["paren_flag"])
        else:
            print stdout.read()
            raise NotImplementedError("Don't know how to extract version")
    return v

def get_versions(config):
    """Retrieve details on all programs available on the system.
    """
    out = [{"program": "bcbio-nextgen",
            "version": version.__version__}]
    for p in _cl_progs:
        out.append({"program": p["cmd"],
                    "version": _get_cl_version(p, config)})
    for p in _alt_progs:
        out.append({"program": p["name"],
                    "version": p["version_fn"](config)})
    return out

def write_versions(dirs, config):
    """Write CSV file with versions used in analysis pipeline.
    """
    base_dir = utils.safe_makedir(os.path.join(dirs["work"], "provenance"))
    out_file = os.path.join(base_dir, "programs.txt")
    if not utils.file_exists(out_file):
        with open(out_file, "w") as out_handle:
            for p in get_versions(config):
                out_handle.write("{program},{version}\n".format(**p))
    return out_file
