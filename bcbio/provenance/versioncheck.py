"""Check specific required program versions required during the pipeline.
"""
from distutils.version import LooseVersion
import subprocess

from bcbio.pipeline import config_utils
from bcbio.log import logger

def samtools(config, items):
    """Ensure samtools has parallel processing required for piped analysis.
    """
    samtools = config_utils.get_program("samtools", config)
    p = subprocess.Popen([samtools, "sort"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, _ = p.communicate()
    p.stdout.close()
    if output.find("-@") == -1:
        return ("Installed version of samtools sort does not have support for multithreading (-@ option) "
                "required to support bwa piped alignment and BAM merging. "
                "Please upgrade to the latest version "
                "from http://samtools.sourceforge.net/")

def _has_pipeline(items):
    """Only perform version checks when we're running an analysis pipeline.
    """
    return any(item.get("analysis", "") != "" for item in items)

def _is_variant(items):
    return any(item.get("analysis", "").lower().startswith("variant") for item in items)

def java(config, items):
    """GATK and derived tools requires Java 1.7.
    """
    want_version = "1.7" if _is_variant(items) else "1.6"
    try:
        java = config_utils.get_program("java", config)
    except config_utils.CmdNotFound:
        return ("java not found on PATH. Java %s required." % want_version)
    p = subprocess.Popen([java, "-Xms250m", "-Xmx250m", "-version"],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, _ = p.communicate()
    p.stdout.close()
    version = ""
    for line in output.split("\n"):
        if line.startswith(("java version", "openjdk version")):
            version = line.strip().split()[-1]
            if version.startswith('"'):
                version = version[1:]
            if version.endswith('"'):
                version = version[:-1]
    if not version or LooseVersion(version) < LooseVersion(want_version):
        return ("java version %s required for running GATK and other tools. "
                "Found version %s at %s" % (want_version, version, java))

def testall(items):
    logger.info("Testing minimum versions of installed programs")
    config = items[0]["config"]
    msgs = []
    if _has_pipeline(items):
        for fn in [samtools, java]:
            out = fn(config, items)
            if out:
                msgs.append(out)
    if msgs:
        raise OSError("Program problems found. You can upgrade dependencies with:\n" +
                      "bcbio_nextgen.py upgrade -u skip --tooldir=/usr/local\n\n" +
                      "\n".join(msgs))
