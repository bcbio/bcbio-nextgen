"""Work with Broad's Java libraries from Python.

  Picard -- BAM manipulation and analysis library.
  GATK -- Next-generation sequence processing.
"""
from contextlib import closing
import copy
from distutils.version import LooseVersion
import os
import subprocess

from bcbio.broad import picardrun
from bcbio.pipeline import config_utils
from bcbio.provenance import do, programs
from bcbio.utils import curdir_tmpdir

class BroadRunner:
    """Simplify running Broad commandline tools.
    """
    def __init__(self, picard_ref, gatk_dir, config):
        resources = config_utils.get_resources("gatk", config)
        self._jvm_opts = resources.get("jvm_opts", ["-Xms750m", "-Xmx2g"])
        self._picard_ref = config_utils.expand_path(picard_ref)
        self._gatk_dir = config_utils.expand_path(gatk_dir) or config_utils.expand_path(picard_ref)
        self._config = config
        self._gatk_version, self._picard_version, self._mutect_version = (
            None, None, None)
        self._gatk_resources = resources

    def _set_default_versions(self, config):
        """Retrieve pre-computed version information for expensive to retrieve versions.
        Starting up GATK takes a lot of resources so we do it once at start of analysis.
        """
        out = []
        for name in ["gatk", "picard", "mutect"]:
            try:
                v = programs.get_version(name, config=config)
            except KeyError:
                v = None
            out.append(v)
        self._gatk_version, self._picard_version, self._mutect_version = out

    def new_resources(self, program):
        """Set new resource usage for the given program.
        This allows customization of memory usage for particular sub-programs
        of GATK like HaplotypeCaller.
        """
        resources = config_utils.get_resources(program, self._config)
        if resources.get("jvm_opts"):
            self._jvm_opts = resources.get("jvm_opts")

    def run_fn(self, name, *args, **kwds):
        """Run pre-built functionality that used Broad tools by name.

        See the gatkrun, picardrun module for available functions.
        """
        fn = None
        to_check = [picardrun]
        for ns in to_check:
            try:
                fn = getattr(ns, name)
                break
            except AttributeError:
                pass
        assert fn is not None, "Could not find function %s in %s" % (name, to_check)
        return fn(self, *args, **kwds)

    def cl_picard(self, command, options):
        """Prepare a Picard commandline.
        """
        options = ["%s=%s" % (x, y) for x, y in options]
        options.append("VALIDATION_STRINGENCY=SILENT")
        return self._get_picard_cmd(command) + options

    def run(self, command, options, pipe=False, get_stdout=False):
        """Run a Picard command with the provided option pairs.
        """
        cl = self.cl_picard(command, options)
        if pipe:
            subprocess.Popen(cl)
        elif get_stdout:
            p = subprocess.Popen(cl, stdout=subprocess.PIPE)
            stdout = p.stdout.read()
            p.wait()
            p.stdout.close()
            return stdout
        else:
            do.run(cl, "Picard {0}".format(command), None)

    def get_picard_version(self, command):
        if self._picard_version is None:
            self._set_default_versions(self._config)
        if self._picard_version:
            return self._picard_version
        if os.path.isdir(self._picard_ref):
            picard_jar = self._get_jar(command)
            cl = ["java", "-Xms64m", "-Xmx128m", "-jar", picard_jar]
        else:
            cl = [self._picard_ref, command]
        cl += ["--version"]
        p = subprocess.Popen(cl, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        version = float(p.stdout.read().split("(")[0])
        self._picard_version = version
        p.wait()
        p.stdout.close()
        return version

    def cl_gatk(self, params, tmp_dir):
        support_nt = set()
        support_nct = set(["BaseRecalibrator"])
        gatk_jar = self._get_jar("GenomeAnalysisTK", ["GenomeAnalysisTKLite"])
        local_args = []
        cores = self._config["algorithm"].get("num_cores", 1)
        config = copy.deepcopy(self._config)
        if cores and int(cores) > 1:
            atype_index = params.index("-T") if params.count("-T") > 0 \
                          else params.index("--analysis_type")
            prog = params[atype_index + 1]
            if prog in support_nt:
                params.extend(["-nt", str(cores)])
            elif prog in support_nct:
                params.extend(["-nct", str(cores)])
                if config["algorithm"].get("memory_adjust") is None:
                    config["algorithm"]["memory_adjust"] = {"direction": "increase",
                                                            "magnitude": int(cores) // 2}
        if LooseVersion(self.gatk_major_version()) > LooseVersion("1.9"):
            if len([x for x in params if x.startswith(("-U", "--unsafe"))]) == 0:
                params.extend(["-U", "LENIENT_VCF_PROCESSING"])
            params.extend(["--read_filter", "BadCigar", "--read_filter", "NotPrimaryAlignment"])
        local_args.append("-Djava.io.tmpdir=%s" % tmp_dir)
        if "keyfile" in self._gatk_resources:
            params = ["-et", "NO_ET", "-K", self._gatk_resources["keyfile"]] + params
        return ["java"] + config_utils.adjust_opts(self._jvm_opts, config) + local_args + \
          ["-jar", gatk_jar] + [str(x) for x in params]

    def cl_mutect(self, params, tmp_dir):

        """Define parameters to run the mutect paired algorithm."""

        gatk_jar = self._get_jar("muTect")
        local_args = []
        config = copy.deepcopy(self._config)

        local_args.append("-Djava.io.tmpdir=%s" % tmp_dir)
        return ["java"] + config_utils.adjust_opts(self._jvm_opts, config) + local_args + \
          ["-jar", gatk_jar] + [str(x) for x in params]

    def run_gatk(self, params, tmp_dir=None, log_error=True, memory_retry=False,
                 data=None, region=None):
        with curdir_tmpdir() as local_tmp_dir:
            if tmp_dir is None:
                tmp_dir = local_tmp_dir
            cl = self.cl_gatk(params, tmp_dir)
            atype_index = cl.index("-T") if cl.count("-T") > 0 \
                          else cl.index("--analysis_type")
            prog = cl[atype_index + 1]
            if memory_retry:
                do.run_memory_retry(cl, "GATK: {0}".format(prog), data, region=region)
            else:
                do.run(cl, "GATK: {0}".format(prog), data, region=region,
                       log_error=log_error)

    def run_mutect(self, params, tmp_dir=None):
        with curdir_tmpdir() as local_tmp_dir:
            if tmp_dir is None:
                tmp_dir = local_tmp_dir
            cl = self.cl_mutect(params, tmp_dir)
            prog = "MuTect"
            do.run(cl, "MuTect: {0}".format(prog), None)

    def get_gatk_version(self):
        """Retrieve GATK version, handling locally and config cached versions.
        Calling version can be expensive due to all the startup and shutdown
        of JVMs, so we prefer cached version information.
        """
        if self._gatk_version is None:
            self._set_default_versions(self._config)
        if self._gatk_version is not None:
            return self._gatk_version
        else:
            gatk_jar = self._get_jar("GenomeAnalysisTK", ["GenomeAnalysisTKLite"])
            cl = ["java", "-Xms64m", "-Xmx128m", "-jar", gatk_jar, "-version"]
            with closing(subprocess.Popen(cl, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout) as stdout:
                out = stdout.read().strip()
                # versions earlier than 2.4 do not have explicit version command,
                # parse from error output from GATK
                if out.find("ERROR") >= 0:
                    flag = "The Genome Analysis Toolkit (GATK)"
                    for line in out.split("\n"):
                        if line.startswith(flag):
                            version = line.split(flag)[-1].split(",")[0].strip()
                else:
                    version = out
            if version.startswith("v"):
                version = version[1:]
            self._gatk_version = version
            return version

    def get_mutect_version(self):
        """Retrieve the Mutect version.
        """
        if self._mutect_version is None:
            self._set_default_versions(self._config)
        return self._mutect_version

    def gatk_type(self):
        """Retrieve type of GATK jar, allowing support for older GATK lite.
        Returns either `lite` (targeting GATK-lite 2.3.9) or `restricted`,
        the latest 2.4+ restricted version of GATK.
        """
        if LooseVersion(self.gatk_major_version()) > LooseVersion("2.3"):
            return "restricted"
        else:
            return "lite"

    def gatk_major_version(self):
        """Retrieve the GATK major version, handling multiple GATK distributions.

        Has special cases for GATK nightly builds, Appistry releases and
        GATK prior to 2.3.
        """
        full_version = self.get_gatk_version()
        # Working with a recent version if using nightlies
        if full_version.startswith("nightly-"):
            return "2.8"
        parts = full_version.split("-")
        if len(parts) == 4:
            appistry_release, version, subversion, githash = parts
        elif len(parts) == 3:
            version, subversion, githash = parts
        # version was not properly implemented in earlier GATKs
        else:
            version = "2.3"
        if version.startswith("v"):
            version = version[1:]
        return version

    def _get_picard_cmd(self, command):
        """Retrieve the base Picard command, handling both shell scripts and directory of jars.
        """
        resources = config_utils.get_resources("picard", self._config)
        if resources.get("jvm_opts"):
            jvm_opts = resources.get("jvm_opts")
        else:
            jvm_opts = self._jvm_opts
        if os.path.isdir(self._picard_ref):
            dist_file = self._get_jar(command)
            return ["java"] + jvm_opts + ["-jar", dist_file]
        else:
            # XXX Cannot currently set JVM opts with picard-tools script
            return [self._picard_ref, command]

    def _get_jar(self, command, alts=None):
        """Retrieve the jar for running the specified command.
        """
        dirs = []
        for bdir in [self._gatk_dir, self._picard_ref]:
            dirs.extend([bdir,
                         os.path.join(bdir, os.pardir, "gatk")])
        if alts is None: alts = []
        for check_cmd in [command] + alts:
            for dir_check in dirs:
                try:
                    check_file = config_utils.get_jar(command, dir_check)
                    return check_file
                except ValueError, msg:
                    if str(msg).find("multiple") > 0:
                        raise
                    else:
                        pass
        raise ValueError("Could not find jar %s in %s:%s" % (command, self._picard_ref, self._gatk_dir))

def _get_picard_ref(config):
    """Handle retrieval of Picard for running, handling multiple cases:

    - A directory of jar files corresponding to individual commands.
    - The ubuntu/debian picard-tools commandline, which provides a shell wrapper around
      a shared jar.

    For a directory, configure resources with:
      picard:
        dir: /path/to/jars

    For the debian commandline, configure with:
      picard:
        cmd: picard-tools
    """
    try:
        picard = config_utils.get_program("picard", config, default="notfound")
    except config_utils.CmdNotFound:
        picard = "notfound"
    if picard == "notfound" or os.path.isdir(picard):
        picard = config_utils.get_program("picard", config, "dir")
    return picard

def runner_from_config(config, program="gatk"):
    return BroadRunner(_get_picard_ref(config),
                       config_utils.get_program(program, config, "dir"),
                       config)
