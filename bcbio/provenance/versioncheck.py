"""Check specific required program versions required during the pipeline.
"""
from distutils.version import LooseVersion
import subprocess

from bcbio import setpath, utils
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.log import logger

def samtools(items):
    """Ensure samtools has parallel processing required for piped analysis.
    """
    samtools = config_utils.get_program("samtools", items[0]["config"])
    p = subprocess.Popen([samtools, "sort", "-h"], stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    output, stderr = p.communicate()
    p.stdout.close()
    p.stderr.close()
    if str(output.decode()).find("-@") == -1 and str(stderr.decode()).find("-@") == -1:
        return ("Installed version of samtools sort does not have support for "
                "multithreading (-@ option) "
                "required to support bwa piped alignment and BAM merging. "
                "Please upgrade to the latest version "
                "from http://samtools.sourceforge.net/")

def _has_pipeline(items):
    """Only perform version checks when we're running an analysis pipeline.
    """
    return any(item.get("analysis", "") != "" for item in items)

def _needs_java(data):
    """Check if a caller needs external java for MuTect.

    No longer check for older GATK (<3.6) versions because of time cost; this
    won't be relevant to most runs so we skip the sanity check.
    """
    vc = dd.get_variantcaller(data)
    if isinstance(vc, dict):
        out = {}
        for k, v in vc.items():
            if not isinstance(v, (list, tuple)):
                v = [v]
            out[k] = v
        vc = out
    elif not isinstance(vc, (list, tuple)):
        vc = [vc]
    if "mutect" in vc or ("somatic" in vc and "mutect" in vc["somatic"]):
        return True
    if "gatk" in vc or "gatk-haplotype" in vc or ("germline" in vc and "gatk-haplotype" in vc["germline"]):
        pass
        # runner = broad.runner_from_config(data["config"])
        # version = runner.get_gatk_version()
        # if LooseVersion(version) < LooseVersion("3.6"):
        #     return True
    return False

def java(items):
    """Check for presence of external Java 1.7 for tools that require it.
    """
    if any([_needs_java(d) for d in items]):
        min_version = "1.7"
        max_version = "1.8"
        with setpath.orig_paths():
            java = utils.which("java")
            if not java:
                return ("java not found on PATH. Java %s required for MuTect and GATK < 3.6." % min_version)
            p = subprocess.Popen([java, "-Xms250m", "-Xmx250m", "-version"],
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            output, _ = p.communicate()
            p.stdout.close()
        version = ""
        for line in output.decode().split("\n"):
            if line.startswith(("java version", "openjdk version")):
                version = line.strip().split()[-1]
                if version.startswith('"'):
                    version = version[1:]
                if version.endswith('"'):
                    version = version[:-1]
        if (not version or LooseVersion(version) >= LooseVersion(max_version) or
            LooseVersion(version) < LooseVersion(min_version)):
            return ("java version %s required for running MuTect and GATK < 3.6.\n"
                    "It needs to be first on your PATH so running 'java -version' give the correct version.\n"
                    "Found version %s at %s" % (min_version, version, java))

def testall(items):
    logger.info("Testing minimum versions of installed programs")
    msgs = []
    if _has_pipeline(items):
        for fn in [samtools, java]:
            out = fn(items)
            if out:
                msgs.append(out)
    if msgs:
        raise OSError("Program problems found. You can upgrade dependencies with:\n" +
                      "bcbio_nextgen.py upgrade -u skip --tooldir=/usr/local\n\n" +
                      "\n".join(msgs))
