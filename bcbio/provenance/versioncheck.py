"""Check specific required program versions required during the pipeline.
"""
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

def testall(items):
    logger.info("Testing minimum versions of installed programs")
    items = [x[0] for x in items]
    config = items[0]["config"]
    msgs = []
    for fn in [samtools]:
        out = fn(config)
        if out:
            msgs.append(out)
    if msgs:
        raise OSError("Program problems found. You can upgrade dependencies with:\n" +
                      "bcbio_nextgen.py upgrade -u skip --tooldir=/usr/local\n\n" +
                      "\n".join(msgs))
