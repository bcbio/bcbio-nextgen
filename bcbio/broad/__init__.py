"""Work with Broad's Java libraries from Python.

  Picard -- BAM manipulation and analysis library.
  GATK -- Next-generation sequence processing.
"""
from contextlib import closing
from distutils.version import LooseVersion
import re
import os
import subprocess

import toolz as tz

from bcbio import utils
from bcbio.broad import picardrun
from bcbio.distributed.transaction import tx_tmpdir
from bcbio.pipeline import config_utils
from bcbio.provenance import do, programs

def get_default_jvm_opts(tmp_dir=None):
    """Retrieve default JVM tuning options

    Avoids issues with multiple spun up Java processes running into out of memory errors.
    Parallel GC can use a lot of cores on big machines and primarily helps reduce task latency
    and responsiveness which are not needed for batch jobs.
    https://github.com/chapmanb/bcbio-nextgen/issues/532#issuecomment-50989027
    https://wiki.csiro.au/pages/viewpage.action?pageId=545034311
    http://stackoverflow.com/questions/9738911/javas-serial-garbage-collector-performing-far-better-than-other-garbage-collect
    """
    opts = ["-XX:+UseSerialGC"]
    if tmp_dir:
        opts.append("-Djava.io.tmpdir=%s" % tmp_dir)
    return opts

def _get_gatk_opts(config, names, tmp_dir=None, memscale=None, include_gatk=True):
    """Retrieve GATK memory specifications, moving down a list of potential specifications.
    """
    if include_gatk:
        opts = ["-U", "LENIENT_VCF_PROCESSING", "--read_filter",
                "BadCigar", "--read_filter", "NotPrimaryAlignment"]
    else:
        opts = []
    jvm_opts = ["-Xms750m", "-Xmx2g"]
    for n in names:
        resources = config_utils.get_resources(n, config)
        if resources and resources.get("jvm_opts"):
            jvm_opts = resources.get("jvm_opts")
            break
    if memscale:
        jvm_opts = config_utils.adjust_opts(jvm_opts, {"algorithm": {"memory_adjust": memscale}})
    jvm_opts += get_default_jvm_opts(tmp_dir)
    return jvm_opts + opts

def get_gatk_framework_opts(config, tmp_dir=None, memscale=None, include_gatk=True):
    return _get_gatk_opts(config, ["gatk-framework", "gatk"], tmp_dir, memscale, include_gatk=include_gatk)

def get_gatk_opts(config, tmp_dir=None, memscale=None, include_gatk=True):
    return _get_gatk_opts(config, ["gatk", "gatk-framework"], tmp_dir, memscale,
                          include_gatk=include_gatk)

def get_gatk_vqsr_opts(config, tmp_dir=None, memscale=None):
    return _get_gatk_opts(config, ["gatk-vqsr", "gatk", "gatk-framework"], tmp_dir, memscale)

def get_picard_opts(config, memscale=None):
    return _get_gatk_opts(config, ["picard", "gatk", "gatk-framework"], memscale=memscale, include_gatk=False)

def _clean_java_out(version_str):
    """Remove extra environmental information reported in java when querying for versions.

    Java will report information like _JAVA_OPTIONS environmental variables in the output.
    """
    out = []
    for line in version_str.split("\n"):
        if line.startswith("Picked up"):
            pass
        else:
            out.append(line)
    return "\n".join(out)

def _check_for_bad_version(version, program):
    if version.find("bad major version") >= 0 or version.find("UnsupportedClassVersionError") >= 0:
        raise ValueError("Problem getting version for %s. "
                         "This often indicates you're not running Java 1.7:\n%s" % (program, version))

def get_gatk_version(gatk_jar):
    cl = ["java", "-Xms128m", "-Xmx256m"] + get_default_jvm_opts() + ["-jar", gatk_jar, "-version"]
    with closing(subprocess.Popen(cl, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout) as stdout:
        out = _clean_java_out(stdout.read().strip())
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
    _check_for_bad_version(version, "GATK")
    return version

def get_mutect_version(mutect_jar):
    """Retrieves version from input jar name since there is not an easy way to get MuTect version.
    Check mutect jar for SomaticIndelDetector, which is an Appistry feature
    """
    cl = ["java", "-Xms128m", "-Xmx256m"] + get_default_jvm_opts() + ["-jar", mutect_jar, "-h"]
    with closing(subprocess.Popen(cl, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout) as stdout:
        if "SomaticIndelDetector" in stdout.read().strip():
            mutect_type = "-appistry"
        else:
            mutect_type = ""
    version = os.path.basename(mutect_jar).lower()
    for to_remove in [".jar", "-standalone", "mutect"]:
        version = version.replace(to_remove, "")
    if version.startswith(("-", ".")):
        version = version[1:]
    if not version:
        raise ValueError("Unable to determine MuTect version from jar file. "
                         "Need to have version contained in jar (ie. muTect-1.1.5.jar): %s" % mutect_jar)
    _check_for_bad_version(version, "MuTect")
    return version + mutect_type

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
            v = tz.get_in(["resources", name, "version"], config)
            if not v:
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

    def cl_picard(self, command, options, memscale=None):
        """Prepare a Picard commandline.
        """
        options = ["%s=%s" % (x, y) for x, y in options]
        options.append("VALIDATION_STRINGENCY=SILENT")
        return self._get_picard_cmd(command, memscale=memscale) + options

    def run(self, command, options, pipe=False, get_stdout=False, memscale=None):
        """Run a Picard command with the provided option pairs.
        """
        cl = self.cl_picard(command, options, memscale=memscale)
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
            cl = ["java", "-Xms64m", "-Xmx128m"] + get_default_jvm_opts() + ["-jar", picard_jar]
        else:
            cl = [self._picard_ref, command]
        cl += ["--version"]
        p = subprocess.Popen(cl, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        # fix for issue #494
        pat = re.compile('([\d|\.]*)(\(\d*\)$)')  # matches '1.96(1510)'
        m = pat.search(p.stdout.read())
        version = m.group(1)
        self._picard_version = version
        p.wait()
        p.stdout.close()
        return version

    def has_gatk(self):
        try:
            self._get_jar("GenomeAnalysisTK", ["GenomeAnalysisTKLite"])
            return True
        except ValueError, msg:
            if "Could not find jar" in str(msg):
                return False
            else:
                raise

    def cl_gatk(self, params, tmp_dir, memscale=None):
        support_nt = set()
        support_nct = set(["BaseRecalibrator"])
        gatk_jar = self._get_jar("GenomeAnalysisTK", ["GenomeAnalysisTKLite"])
        cores = self._config["algorithm"].get("num_cores", 1)
        config = self._config
        if cores and int(cores) > 1:
            atype_index = params.index("-T") if params.count("-T") > 0 \
                          else params.index("--analysis_type")
            prog = params[atype_index + 1]
            if prog in support_nt:
                params.extend(["-nt", str(cores)])
            elif prog in support_nct:
                params.extend(["-nct", str(cores)])
                if config["algorithm"].get("memory_adjust") is None:
                    config = utils.deepish_copy(config)
                    config["algorithm"]["memory_adjust"] = {"direction": "increase",
                                                            "magnitude": int(cores) // 2}
        if LooseVersion(self.gatk_major_version()) > LooseVersion("1.9"):
            if len([x for x in params if x.startswith(("-U", "--unsafe"))]) == 0:
                params.extend(["-U", "LENIENT_VCF_PROCESSING"])
            params.extend(["--read_filter", "BadCigar", "--read_filter", "NotPrimaryAlignment"])
        if memscale:
            jvm_opts = get_gatk_opts(config, tmp_dir=tmp_dir, memscale=memscale, include_gatk=False)
        else:
            # Decrease memory slightly from configuration to avoid memory allocation errors
            jvm_opts = config_utils.adjust_opts(self._jvm_opts,
                                                {"algorithm": {"memory_adjust":
                                                               {"magnitude": 1.1, "direction": "decrease"}}})
            jvm_opts += get_default_jvm_opts(tmp_dir)
        if "keyfile" in self._gatk_resources:
            params = ["-et", "NO_ET", "-K", self._gatk_resources["keyfile"]] + params
        return ["java"] + jvm_opts + ["-jar", gatk_jar] + [str(x) for x in params]

    def cl_mutect(self, params, tmp_dir):
        """Define parameters to run the mutect paired algorithm.
        """
        gatk_jar = self._get_jar("muTect", ["mutect"])
        # Decrease memory slightly from configuration to avoid memory allocation errors
        jvm_opts = config_utils.adjust_opts(self._jvm_opts,
                                            {"algorithm": {"memory_adjust":
                                                           {"magnitude": 1.1, "direction": "decrease"}}})
        return ["java"] + jvm_opts + get_default_jvm_opts(tmp_dir) + \
               ["-jar", gatk_jar] + [str(x) for x in params]

    def run_gatk(self, params, tmp_dir=None, log_error=True,
                 data=None, region=None, memscale=None):
        with tx_tmpdir(self._config) as local_tmp_dir:
            if tmp_dir is None:
                tmp_dir = local_tmp_dir
            cl = self.cl_gatk(params, tmp_dir, memscale=memscale)
            atype_index = cl.index("-T") if cl.count("-T") > 0 \
                          else cl.index("--analysis_type")
            prog = cl[atype_index + 1]
            do.run(cl, "GATK: {0}".format(prog), data, region=region,
                   log_error=log_error)

    def run_mutect(self, params, tmp_dir=None):
        with tx_tmpdir(self._config) as local_tmp_dir:
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
            self._gatk_version = get_gatk_version(gatk_jar)
            return self._gatk_version

    def get_mutect_version(self):
        """Retrieve the Mutect version.
        """
        if self._mutect_version is None:
            mutect_jar = self._get_jar("muTect", ["mutect"])
            self._mutect_version = get_mutect_version(mutect_jar)
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

    def _get_picard_cmd(self, command, memscale=None):
        """Retrieve the base Picard command, handling both shell scripts and directory of jars.
        """
        resources = config_utils.get_resources("picard", self._config)
        if memscale:
            jvm_opts = get_picard_opts(self._config, memscale=memscale)
        elif resources.get("jvm_opts"):
            jvm_opts = resources.get("jvm_opts")
        else:
            jvm_opts = self._jvm_opts
        if os.path.isdir(self._picard_ref):
            dist_file = self._get_jar(command)
            return ["java"] + jvm_opts + get_default_jvm_opts() + ["-jar", dist_file]
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
                    check_file = config_utils.get_jar(check_cmd, dir_check)
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

def runner_from_config_safe(config):
    """Retrieve a runner, returning None if GATK is not available.
    """
    try:
        return runner_from_config(config)
    except ValueError, msg:
        if str(msg).find("Could not find directory in config for gatk") >= 0:
            return None
        else:
            raise

class PicardCmdRunner:
    def __init__(self, cmd, config):
        self._cmd = cmd
        self._config = config

    def run(self, subcmd, opts, memscale=None):
        jvm_opts = get_picard_opts(self._config, memscale=memscale)
        cmd = [self._cmd] + jvm_opts + [subcmd] + ["%s=%s" % (x, y) for x, y in opts] + \
              ["VALIDATION_STRINGENCY=SILENT"]
        do.run(cmd, "Picard: %s" % subcmd)

    def run_fn(self, name, *args, **kwds):
        """Run pre-built functionality that used Broad tools by name.

        See the picardrun module for available functions.
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

def runner_from_path(cmd, config):
    """Simple command line runner that expects a bash cmd in the PATH.

    This makes Picard tools back compatible with new approach of a single
    jar + bash script.
    """
    if cmd.endswith("picard"):
        return PicardCmdRunner(cmd, config)
    else:
        raise ValueError("Do not support PATH running for %s" % cmd)