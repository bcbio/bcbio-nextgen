"""Check specific required program versions required during the pipeline.
"""
from distutils.version import LooseVersion
import subprocess

from bcbio.pipeline import config_utils
from bcbio.log import logger

def samtools(config):
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

def pdflatex(config):
    """Check for pdflatex if write_summary enabled.
    """
    if config["algorithm"].get("write_summary", True):
        try:
            config_utils.get_program("pdflatex", config)
        except config_utils.CmdNotFound:
            return ("pdflatex not found, required to build final pdf reports. "
                    "Install a LaTeX package like TeX Live or set `write_summary` to false "
                    "in the `algorithm` section of your sample configuration.")

def java(config):
    """GATK requires Java 1.7 or better.
    """
    try:
        java = config_utils.get_program("java", config)
    except config_utils.CmdNotFound:
        return ("java not found on PATH. Java 1.7 or better required.")
    p = subprocess.Popen([java, "-version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, _ = p.communicate()
    p.stdout.close()
    version = ""
    for line in output.split("\n"):
        if line.startswith("java version"):
            version = line.strip().split()[-1]
            if version.startswith('"'):
                version = version[1:]
            if version.endswith('"'):
                version = version[:-1]
    if not version or LooseVersion(version) < LooseVersion("1.7.0"):
        return ("java version 1.7 or better required for running GATK and other tools. "
                "Found version %s at %s" % (version, java))

def testall(items):
    logger.info("Testing minimum versions of installed programs")
    items = [x[0] for x in items]
    config = items[0]["config"]
    msgs = []
    for fn in [samtools, pdflatex, java]:
        out = fn(config)
        if out:
            msgs.append(out)
    if msgs:
        raise OSError("Program problems found. You can upgrade dependencies with:\n" +
                      "bcbio_nextgen.py upgrade -u skip --tooldir=/usr/local\n\n" +
                      "\n".join(msgs))
