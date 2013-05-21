"""Centralize running of external commands, providing logging and tracking.
"""
import os
import subprocess

from bcbio import utils
from bcbio.log import logger
from bcbio.provenance import diagnostics

def run(cmd, descr, data, checks=None):
    """Run the provided command, logging details and checking for errors.
    """
    if data:
        descr = "{0} : {1}".format(descr, data["name"][-1])
    logger.info(descr)
    # TODO: Extract entity information from data input
    cmd_id = diagnostics.start_cmd(descr, data, cmd)
    try:
        _do_run(cmd, checks)
    except:
        diagnostics.end_cmd(cmd_id, False)
        raise
    finally:
        diagnostics.end_cmd(cmd_id)

def _do_run(cmd, checks):
    """Perform running and check results, raising errors for issues.
    """
    if isinstance(cmd, basestring):
        subprocess.check_call(cmd, shell=True)
    else:
        subprocess.check_call(cmd)
    # Check for problems not identified by shell return codes
    if checks:
        for check in checks:
            if not check():
                raise IOError("External command failed: %s" %
                              " ".join(cmd) if not isinstance(cmd, basestring) else cmd)

# checks for validating run completed successfully

def file_nonempty(target_file):
    def check():
        ok = utils.file_exists(target_file)
        if not ok:
            logger.info("Did not find non-empty output file {0}".format(target_file))
        return ok
    return check

def file_exists(target_file):
    def check():
        ok = os.path.exists(target_file)
        if not ok:
            logger.info("Did not find output file {0}".format(target_file))
        return ok
    return check
