"""Access tool command lines, handling back compatibility and file type issues.

Abstracts out
"""
import subprocess
import toolz as tz
from bcbio.pipeline import config_utils

def get_tabix_cmd(config):
    """Retrieve tabix command, handling new bcftools tabix and older tabix.
    """
    try:
        bcftools = config_utils.get_program("bcftools", config)
        # bcftools has terrible error codes and stderr output, swallow those.
        bcftools_tabix = subprocess.check_output("{bcftools} 2>&1; echo $?".format(**locals()),
                                                 shell=True).find("tabix") >= 0
    except config_utils.CmdNotFound:
        bcftools_tabix = False
    if bcftools_tabix:
        return "{0} tabix".format(bcftools)
    else:
        tabix = config_utils.get_program("tabix", config)
        return tabix

def get_bgzip_cmd(config, is_retry=False):
    """Retrieve command to use for bgzip, trying to use parallel pbgzip if available.

    NOTE: pbgzip is experimental and may not be stable at the moment. Use it
    with this information in mind. To enable pbgzip, you need to add it into
    algorithm --> tools_on list in the configuration file.

    Avoids over committing cores to gzipping since run in pipe with other tools.
    Allows for retries which force single core bgzip mode.
    """
    num_cores = tz.get_in(["algorithm", "num_cores"], config, 1)
    if not is_retry and num_cores > 1 and \
            "pbgzip" in tz.get_in(["algorithm", "tools_on"], config, []):
        try:
            pbgzip = config_utils.get_program("pbgzip", config)
            return "%s -n %s " % (pbgzip, num_cores)
        except config_utils.CmdNotFound:
            pass
    return config_utils.get_program("bgzip", config)
