"""Identify program versions used for analysis, reporting in structured table.

Catalogs the full list of programs used in analysis, enabling reproduction of
results and tracking of provenance in output files.
"""
import os
import contextlib
import subprocess
import sys
import yaml
import toolz as tz

from bcbio import utils
from bcbio.pipeline import config_utils, version
from bcbio.pipeline import datadict as dd
from bcbio.log import logger

_cl_progs = [{"cmd": "bamtofastq", "name": "biobambam",
              "args": "--version", "stdout_flag": "This is biobambam version"},
             {"cmd": "bamtools", "args": "--version", "stdout_flag": "bamtools"},
             {"cmd": "bcftools", "stdout_flag": "Version:"},
             {"cmd": "bedtools", "args": "--version", "stdout_flag": "bedtools"},
             {"cmd": "bowtie2", "args": "--version", "stdout_flag": "bowtie2-align version"},
             {"cmd": "bwa", "stdout_flag": "Version:"},
             {"cmd": "chanjo"},
             {"cmd": "cutadapt", "args": "--version"},
             {"cmd": "fastqc", "args": "--version", "stdout_flag": "FastQC"},
             {"cmd": "freebayes", "stdout_flag": "version:"},
             {"cmd": "gemini", "args": "--version", "stdout_flag": "gemini "},
             {"cmd": "novosort", "paren_flag": "novosort"},
             {"cmd": "novoalign", "stdout_flag": "Novoalign"},
             {"cmd": "samtools", "stdout_flag": "Version:"},
             {"cmd": "qualimap", "args": "-h", "stdout_flag": "QualiMap"},
             {"cmd": "vcflib", "has_cl_version": False},
             {"cmd": "featurecounts", "args": "-v", "stdout_flag": "featureCounts"}]
_manifest_progs = ["BubbleTree", "cufflinks-binary", "cnvkit", "gatk-framework", "grabix", "htseq",
                   "lumpy-sv", "manta", "metasv", "phylowgs", "platypus-variant", "rna-star",
                   "rtg-tools","sambamba-binary", "samblaster", "scalpel", "vardict",
                   "vardict-java", "vep", "vt", "wham"]

def _broad_versioner(type):
    def get_version(config):
        from bcbio import broad
        try:
            runner = broad.runner_from_config(config)
        except ValueError:
            return ""
        if type == "gatk":
            return runner.get_gatk_version()
        elif type == "picard":
            return runner.get_picard_version("ViewSam")
        elif type == "mutect":
            try:
                runner = broad.runner_from_config(config, "mutect")
            except ValueError:
                return ""
            return runner.get_mutect_version()
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
        if not jar:
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

_alt_progs = [{"name": "bcbio_variation",
               "version_fn": jar_versioner("bcbio_variation", "bcbio.variation")},
              {"name": "gatk", "version_fn": _broad_versioner("gatk")},
              {"name": "mutect",
               "version_fn": _broad_versioner("mutect")},
              {"name": "picard", "version_fn": _broad_versioner("picard")},
              {"name": "rnaseqc",
               "version_fn": jar_versioner("rnaseqc", "RNA-SeQC")},
              {"name": "snpeff",
               "version_fn": java_versioner("snpeff", "snpEff", stdout_flag="snpEff version SnpEff")},
              {"name": "varscan",
               "version_fn": jar_versioner("varscan", "VarScan")},
              {"name": "oncofuse",
               "version_fn": jar_versioner("Oncofuse", "Oncofuse")},
              {"name": "alientrimmer",
               "version_fn": jar_versioner("AlienTrimmer", "AlienTrimmer")}
]

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

        localpy_cmd = os.path.join(os.path.dirname(sys.executable), p["cmd"])
        if os.path.exists(localpy_cmd):
            prog = localpy_cmd
        else:
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
            lines = [l.strip() for l in stdout.read().split("\n") if l.strip()]
            v = lines[-1]
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
        vout = subprocess.check_output([brew_cmd, "list", "--versions"])
    except OSError:  # brew not installed/used
        vout = ""
    out = {}
    for vstr in vout.split("\n"):
        if vstr.strip():
            parts = vstr.rstrip().split()
            name = parts[0]
            v = parts[-1]
            out[name] = v
    return out

def _get_versions(config=None):
    """Retrieve details on all programs available on the system.
    """
    out = [{"program": "bcbio-nextgen",
            "version": ("%s-%s" % (version.__version__, version.__git_revision__)
                        if version.__git_revision__ else version.__version__)}]
    manifest_dir = _get_manifest_dir(config)
    manifest_vs = _get_versions_manifest(manifest_dir)
    if manifest_vs:
        out += manifest_vs
    else:
        assert config is not None, "Need configuration to retrieve from non-manifest installs"
        brew_vs = _get_brew_versions()
        for p in _cl_progs:
            out.append({"program": p["cmd"],
                        "version": (brew_vs[p["cmd"]] if p["cmd"] in brew_vs else
                                    _get_cl_version(p, config))})
        for p in _alt_progs:
            out.append({"program": p["name"],
                        "version": (brew_vs[p["name"]] if p["name"] in brew_vs else
                                    p["version_fn"](config))})
    out.sort(key=lambda x: x["program"])
    return out

def _get_manifest_dir(data=None):
    """
    get manifest directory from the data dictionary, falling back on alternatives
    it prefers, in order:
    1. locating it from the bcbio_system.yaml file
    2. locating it from the galaxy directory
    3. location it from the python executable. 

    it can accept either the data or config dictionary
    """
    manifest_dir = None
    if data:
        bcbio_system = tz.get_in(["config", "bcbio_system"], data, None)
        bcbio_system = bcbio_system if bcbio_system else data.get("bcbio_system", None)
        if bcbio_system:
            sibling_dir = os.path.normpath(os.path.dirname(bcbio_system))
        else:
            sibling_dir = dd.get_galaxy_dir(data)
        if sibling_dir:
            manifest_dir = os.path.normpath(os.path.join(sibling_dir, os.pardir,
                                                         "manifest"))
    if not manifest_dir or not os.path.exists(manifest_dir):
        manifest_dir = os.path.join(config_utils.get_base_installdir(), "manifest")
    return manifest_dir

def _get_versions_manifest(manifest_dir):
    """Retrieve versions from a pre-existing manifest of installed software.
    """
    all_pkgs = _manifest_progs + [p.get("name", p["cmd"]) for p in _cl_progs] + [p["name"] for p in _alt_progs]
    if os.path.exists(manifest_dir):
        out = []
        for plist in ["toolplus", "brew", "python", "r", "debian", "custom"]:
            pkg_file = os.path.join(manifest_dir, "%s-packages.yaml" % plist)
            if os.path.exists(pkg_file):
                with open(pkg_file) as in_handle:
                    pkg_info = yaml.safe_load(in_handle)
                added = []
                for pkg in all_pkgs:
                    if pkg in pkg_info:
                        added.append(pkg)
                        out.append({"program": pkg, "version": pkg_info[pkg]["version"]})
                for x in added:
                    all_pkgs.remove(x)
        out.sort(key=lambda x: x["program"])
        for pkg in all_pkgs:
            out.append({"program": pkg, "version": ""})
        return out

def _get_program_file(dirs):
    if dirs.get("work"):
        base_dir = utils.safe_makedir(os.path.join(dirs["work"], "provenance"))
        return os.path.join(base_dir, "programs.txt")

def write_versions(dirs, config=None, is_wrapper=False):
    """Write CSV file with versions used in analysis pipeline.
    """
    out_file = _get_program_file(dirs)
    if is_wrapper:
        assert utils.file_exists(out_file), "Failed to create program versions from VM"
    elif out_file is None:
        for p in _get_versions(config):
            print("{program},{version}".format(**p))
    else:
        with open(out_file, "w") as out_handle:
            for p in _get_versions(config):
                out_handle.write("{program},{version}\n".format(**p))
    return out_file

def get_version_manifest(name, data=None, required=False):
    """Retrieve a version from the currently installed manifest.
    """
    manifest_dir = _get_manifest_dir(data)
    manifest_vs = _get_versions_manifest(manifest_dir)
    for x in manifest_vs:
        if x["program"] == name:
            v = x.get("version", "")
            if v:
                return v
    if required:
        raise ValueError("Did not find %s in install manifest. Could not check version." % name)
    return ""

def add_subparser(subparsers):
    """Add command line option for exporting version information.
    """
    parser = subparsers.add_parser("version",
                                   help="Export versions of used software to stdout or a file ")
    parser.add_argument("--workdir", help="Directory export programs to in workdir/provenance/programs.txt",
                        default=None)

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
