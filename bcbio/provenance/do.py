"""Centralize running of external commands, providing logging and tracking.
"""
import collections
import contextlib
import os
import subprocess
import time

from bcbio import utils
from bcbio.log import logger, logger_cl, logger_stdout
from bcbio.provenance import diagnostics

def run(cmd, descr, data=None, checks=None, region=None, log_error=True,
        log_stdout=False):
    """Run the provided command, logging details and checking for errors.
    """
    descr = _descr_str(descr, data, region)
    logger.debug(descr)
    # TODO: Extract entity information from data input
    cmd_id = diagnostics.start_cmd(descr, data, cmd)
    try:
        logger_cl.debug(" ".join(cmd) if not isinstance(cmd, basestring) else cmd)
        _do_run(cmd, checks, log_stdout)
    except:
        diagnostics.end_cmd(cmd_id, False)
        if log_error:
            logger.exception()
        raise
    finally:
        diagnostics.end_cmd(cmd_id)

def run_memory_retry(cmd, descr, data=None, check=None, region=None):
    """Run command, retrying when detecting fail due to memory errors.

    This is useful for high throughput Java jobs which fail
    intermittently due to an inability to get system resources.
    """
    max_runs = 5
    num_runs = 0
    while 1:
        try:
            run(cmd, descr, data, check, region=region, log_error=False)
            break
        except subprocess.CalledProcessError, msg:
            if num_runs < max_runs and ("insufficient memory" in str(msg) or
                                        "did not provide enough memory" in str(msg) or
                                        "A fatal error has been detected" in str(msg) or
                                        "java.lang.OutOfMemoryError" in str(msg) or
                                        "Resource temporarily unavailable" in str(msg)):
                logger.info("Retrying job. Memory or resource issue with run: %s"
                            % _descr_str(descr, data, region))
                time.sleep(30)
                num_runs += 1
            else:
                logger.exception()
                raise

def _descr_str(descr, data, region):
    """Add additional useful information from data to description string.
    """
    if data:
        if "name" in data:
            descr = "{0} : {1}".format(descr, data["name"][-1])
        elif "work_bam" in data:
            descr = "{0} : {1}".format(descr, os.path.basename(data["work_bam"]))
    if region:
        descr = "{0} : {1}".format(descr, region)
    return descr

def _find_bash():
    try:
        which_bash = subprocess.check_output(["which", "bash"]).strip()
    except subprocess.CalledProcessError:
        which_bash = None
    for test_bash in [which_bash, "/bin/bash", "/usr/bin/bash", "/usr/local/bin/bash"]:
        if test_bash and os.path.exists(test_bash):
            return test_bash
    raise IOError("Could not find bash in any standard location. Needed for unix pipes")

def _normalize_cmd_args(cmd):
    """Normalize subprocess arguments to handle list commands, string and pipes.
    Piped commands set pipefail and require use of bash to help with debugging
    intermediate errors.
    """
    if isinstance(cmd, basestring):
        # check for standard or anonymous named pipes
        if cmd.find(" | ") > 0 or cmd.find(">(") or cmd.find("<("):
            return "set -o pipefail; " + cmd, True, _find_bash()
        else:
            return cmd, True, None
    else:
        return cmd, False, None

def _do_run(cmd, checks, log_stdout=False):
    """Perform running and check results, raising errors for issues.
    """
    cmd, shell_arg, executable_arg = _normalize_cmd_args(cmd)
    s = subprocess.Popen(cmd, shell=shell_arg, executable=executable_arg,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, close_fds=True)
    debug_stdout = collections.deque(maxlen=100)
    with contextlib.closing(s.stdout) as stdout:
        while 1:
            line = stdout.readline()
            if line:
                debug_stdout.append(line)
                if log_stdout:
                    logger_stdout.debug(line.rstrip())
                else:
                    logger.debug(line.rstrip())
            exitcode = s.poll()
            if exitcode is not None:
                if exitcode is not None and exitcode != 0:
                    error_msg = " ".join(cmd) if not isinstance(cmd, basestring) else cmd
                    error_msg += "\n"
                    error_msg += "".join(debug_stdout)
                    raise subprocess.CalledProcessError(exitcode, error_msg)
                else:
                    break
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

def file_reasonable_size(target_file, input_file):
    def check():
        # named pipes -- we can't calculate size
        if input_file.strip().startswith("<("):
            return True
        if input_file.endswith((".bam", ".gz")):
            scale = 5.0
        else:
            scale = 10.0
        orig_size = os.path.getsize(input_file) / pow(1024.0, 3)
        out_size = os.path.getsize(target_file) / pow(1024.0, 3)
        if out_size < (orig_size / scale):
            logger.info("Output file unexpectedly small. %.1fGb for output versus "
                        "%.1fGb for the input file. This often indicates a truncated "
                        "BAM file or memory errors during the run." % (out_size, orig_size))
            return False
        else:
            return True
    return check
