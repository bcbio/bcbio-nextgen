"""Identify program versions used for analysis, reporting in structured table.

Catalogs the full list of programs used in analysis, enabling reproduction of
results and tracking of provenance in output files.
"""
import os
import contextlib
import subprocess

from bcbio import utils
from bcbio.pipeline import config_utils, version
from bcbio.log import logger

import HTSeq

_cl_progs = [{"cmd": "bamtofastq", "args": "--version", "stdout_flag": "This is biobambam version"},
             {"cmd": "bamtools", "args": "--version", "stdout_flag": "bamtools"},
             {"cmd": "bcftools", "stdout_flag": "Version:"},
             {"cmd": "bedtools", "args": "--version", "stdout_flag": "bedtools"},
             {"cmd": "bowtie2", "args": "--version", "stdout_flag": "bowtie2-align version"},
             {"cmd": "bwa", "stdout_flag": "Version:"},
             {"cmd": "cufflinks", "stdout_flag": "cufflinks"},
             {"cmd": "cutadapt", "args": "--version"},
             {"cmd": "fastqc", "args": "--version", "stdout_flag": "FastQC"},
             {"cmd": "freebayes", "stdout_flag": "version:"},
             {"cmd": "gemini", "args": "--version", "stdout_flag": "gemini "},
             {"cmd": "novosort", "paren_flag": "novosort"},
             {"cmd": "novoalign", "stdout_flag": "Novoalign"},
             {"cmd": "samtools", "stdout_flag": "Version:"},
             {"cmd": "sambamba", "stdout_flag": "sambamba"},
             {"cmd": "qualimap", "args": "-h", "stdout_flag": "QualiMap"},
             {"cmd": "tophat", "args": "--version", "stdout_flag": "TopHat"},
             {"cmd": "vcflib", "has_cl_version": False}]

def _broad_versioner(type):
    def get_version(config):
        from bcbio import broad
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
        if jar is "":
            logger.warn("Unable to determine version for program '{}' from jar file {}".format(
                program_name, config_utils.get_jar(jar_name, pdir)))
        return jar
    return get_version

def java_versioner(pname, jar_name, **kwargs):
    def get_version(config):
        try:
            pdir = config_utils.get_program(pname, config, "dir")
        except ValueError:
            return ""
        jar = config_utils.get_jar(jar_name, pdir)
        kwargs["cmd"] = "java"
        kwargs["args"] = "-Xms128m -Xmx256m -jar %s" % jar
        return _get_cl_version(kwargs, config)
    return get_version

_alt_progs = [{"name": "bcbio.variation",
               "version_fn": jar_versioner("bcbio_variation", "bcbio.variation")},
              {"name": "gatk", "version_fn": _broad_versioner("gatk")},
              {"name": "mutect",
               "version_fn": jar_versioner("mutect", "muTect")},
              {"name": "picard", "version_fn": _broad_versioner("picard")},
              {"name": "rnaseqc",
               "version_fn": jar_versioner("rnaseqc", "RNA-SeQC")},
              {"name": "snpeff",
               "version_fn": java_versioner("snpEff", "snpEff", stdout_flag="snpEff version SnpEff")},
              {"name": "varscan",
               "version_fn": jar_versioner("varscan", "VarScan")}]

def _parse_from_stdoutflag(stdout, x):
    for line in stdout:
        if line.find(x) >= 0:
            parts = [p for p in line[line.find(x) + len(x):].split() if p.strip()]
            return parts[0].strip()
    return ""

def _parse_from_parenflag(stdout, x):
    for line in stdout:
        if line.find(x) >= 0:
            return line.split("(")[-1].split(")")[0]
    return ""

def _get_cl_version(p, config):
    """Retrieve version of a single commandline program.
    """
    if not p.get("has_cl_version", True):
        return ""
    try:
        prog = config_utils.get_program(p["cmd"], config)
    except config_utils.CmdNotFound:
        return ""
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
            v = stdout.read().strip()
    if v.endswith("."):
        v = v[:-1]
    return v

def _get_brew_versions():
    """Retrieve versions of tools installed via brew.
    """
    from bcbio import install
    tooldir = install.get_defaults().get("tooldir")
    brew_cmd = os.path.join(tooldir, "bin", "brew") if tooldir else "brew"
    try:
        vout = subprocess.check_output([brew_cmd, "which"])
        uses_which = True
    except subprocess.CalledProcessError:
        vout = subprocess.check_output([brew_cmd, "list", "--versions"])
        uses_which = False
    except OSError:  # brew not installed/used
        vout = ""
    out = {}
    for vstr in vout.split("\n"):
        if vstr.strip():
            if uses_which:
                name, v = vstr.rstrip().split(": ")
            else:
                parts = vstr.rstrip().split()
                name = parts[0]
                v = parts[-1]
            out[name] = v
    return out

def _get_versions(config):
    """Retrieve details on all programs available on the system.
    """
    brew_vs = _get_brew_versions()
    out = [{"program": "bcbio-nextgen",
            "version": ("%s-%s" % (version.__version__, version.__git_revision__)
                        if version.__git_revision__ else version.__version__)},
           {"program": "htseq", "version": HTSeq.__version__}]
    for p in _cl_progs:
        out.append({"program": p["cmd"],
                    "version": (brew_vs[p["cmd"]] if p["cmd"] in brew_vs else
                                _get_cl_version(p, config))})
    for p in _alt_progs:
        out.append({"program": p["name"],
                    "version": (brew_vs[p["name"]] if p["name"] in brew_vs else
                                p["version_fn"](config))})
    return out

def _get_program_file(dirs):
    base_dir = utils.safe_makedir(os.path.join(dirs["work"], "provenance"))
    return os.path.join(base_dir, "programs.txt")

def write_versions(dirs, config):
    """Write CSV file with versions used in analysis pipeline.
    """
    out_file = _get_program_file(dirs)
    with open(out_file, "w") as out_handle:
        for p in _get_versions(config):
            out_handle.write("{program},{version}\n".format(**p))
    return out_file

def get_version(name, dirs=None, config=None):
    """Retrieve the current version of the given program from cached names.
    """
    if dirs:
        p = _get_program_file(dirs)
    else:
        p = config["resources"]["program_versions"]
    with open(p) as in_handle:
        for line in in_handle:
            prog, version = line.rstrip().split(",")
            if prog == name and version:
                return version
    raise KeyError("Version information not found for %s in %s" % (name, p))
