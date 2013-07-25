"""Centralize running of external commands, providing logging and tracking.
"""
import collections
import contextlib
import os
import subprocess

from bcbio import utils
from bcbio.log import logger, logger_cl
from bcbio.provenance import diagnostics

def run(cmd, descr, data, checks=None):
    """Run the provided command, logging details and checking for errors.
    """
    if data:
        descr = "{0} : {1}".format(descr, data["name"][-1])
    logger.debug(descr)
    # TODO: Extract entity information from data input
    cmd_id = diagnostics.start_cmd(descr, data, cmd)
    try:
        logger_cl.debug(" ".join(cmd) if not isinstance(cmd, basestring) else cmd)
        _do_run(cmd, checks)
    except:
        diagnostics.end_cmd(cmd_id, False)
        logger.exception()
        raise
    finally:
        diagnostics.end_cmd(cmd_id)

def _do_run(cmd, checks):
    """Perform running and check results, raising errors for issues.
    """
    s = subprocess.Popen(cmd, shell=isinstance(cmd, basestring),
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    debug_stdout = collections.deque(maxlen=100)
    with contextlib.closing(s.stdout) as stdout:
        while 1:
            line = stdout.readline()
            exitcode = s.poll()
            if exitcode is not None:
                if exitcode is not None and exitcode != 0:
                    error_msg = " ".join(cmd) if not isinstance(cmd, basestring) else cmd
                    error_msg += "\n"
                    error_msg += "".join(debug_stdout)
                    raise subprocess.CalledProcessError(exitcode, error_msg)
                else:
                    break
            if line:
                debug_stdout.append(line)
                logger.debug(line.rstrip())
    # Check for problems not identified by shell return codes
    if checks:
        for check in checks:
            if not check():
                raise IOError("External command failed")

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
