"""Access tool command lines, handling back compatibility and file type issues.

Abstracts out
"""
import subprocess
import toolz as tz
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils

def get_tabix_cmd(config):
    """Retrieve tabix command, handling new bcftools tabix and older tabix.
    """
    try:
        bcftools = config_utils.get_program("bcftools", config)
        # bcftools has terrible error codes and stderr output, swallow those.
        bcftools_tabix = subprocess.check_output("{bcftools} 2>&1; echo $?".format(**locals()),
                                                 shell=True).decode().find("tabix") >= 0
    except config_utils.CmdNotFound:
        bcftools_tabix = False
    if bcftools_tabix:
        return "{0} tabix".format(bcftools)
    else:
        tabix = config_utils.get_program("tabix", config)
        return tabix

def get_bgzip_cmd(config, is_retry=False):
    """Retrieve command to use for bgzip, trying to use bgzip parallel threads.

    By default, parallel bgzip is enabled in bcbio. If it causes problems
    please report them. You can turn parallel bgzip off with `tools_off: [pbgzip]`
    """
    num_cores = tz.get_in(["algorithm", "num_cores"], config, 1)
    cmd = config_utils.get_program("bgzip", config)
    if (not is_retry and num_cores > 1 and
          "pbgzip" not in dd.get_tools_off({"config": config})):
        cmd += " --threads %s" % num_cores
    return cmd
