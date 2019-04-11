"""Centralize running of external commands, providing logging and tracking.
"""
import collections
import os
import subprocess

from bcbio import utils
from bcbio.log import logger, logger_cl, logger_stdout
from bcbio.pipeline import datadict as dd
from bcbio.provenance import diagnostics


import six


def run(cmd, descr=None, data=None, checks=None, region=None, log_error=True,
        log_stdout=False, env=None):
    """Run the provided command, logging details and checking for errors.
    """
    if descr:
      descr = _descr_str(descr, data, region)
      logger.debug(descr)
    cmd_id = diagnostics.start_cmd(cmd, descr or "", data)
    try:
        logger_cl.debug(" ".join(str(x) for x in cmd) if not isinstance(cmd, six.string_types) else cmd)
        _do_run(cmd, checks, log_stdout, env=env)
    except:
        diagnostics.end_cmd(cmd_id, False)
        if log_error:
            logger.exception()
        raise
    finally:
        diagnostics.end_cmd(cmd_id)

def _descr_str(descr, data, region):
    """Add additional useful information from data to description string.
    """
    if data:
        name = dd.get_sample_name(data)
        if name:
            descr = "{0} : {1}".format(descr, name)
        elif "work_bam" in data:
            descr = "{0} : {1}".format(descr, os.path.basename(data["work_bam"]))
    if region:
        descr = "{0} : {1}".format(descr, region)
    return descr

def find_bash():
    for test_bash in [find_cmd("bash"), "/bin/bash", "/usr/bin/bash", "/usr/local/bin/bash"]:
        if test_bash and os.path.exists(test_bash):
            return test_bash
    raise IOError("Could not find bash in any standard location. Needed for unix pipes")

def find_cmd(cmd):
    try:
        return subprocess.check_output(["which", cmd]).decode().strip()
    except subprocess.CalledProcessError:
        return None

def _normalize_cmd_args(cmd):
    """Normalize subprocess arguments to handle list commands, string and pipes.
    Piped commands set pipefail and require use of bash to help with debugging
    intermediate errors.
    """
    if isinstance(cmd, six.string_types):
        # check for standard or anonymous named pipes
        if cmd.find(" | ") > 0 or cmd.find(">(") or cmd.find("<("):
            return "set -o pipefail; " + cmd, True, find_bash()
        else:
            return cmd, True, None
    else:
        return [str(x) for x in cmd], False, None

def _do_run(cmd, checks, log_stdout=False, env=None):
    """Perform running and check results, raising errors for issues.
    """
    cmd, shell_arg, executable_arg = _normalize_cmd_args(cmd)
    s = subprocess.Popen(
        cmd,
        shell=shell_arg,
        executable=executable_arg,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        close_fds=True,
        env=env,
    )
    debug_stdout = collections.deque(maxlen=100)
    while 1:
        line = s.stdout.readline().decode("utf-8", errors="replace")
        if line.rstrip():
            debug_stdout.append(line)
            if log_stdout:
                logger_stdout.debug(line.rstrip())
            else:
                logger.debug(line.rstrip())
        exitcode = s.poll()
        if exitcode is not None:
            for line in s.stdout:
                debug_stdout.append(line.decode("utf-8", errors="replace"))
            if exitcode is not None and exitcode != 0:
                error_msg = " ".join(cmd) if not isinstance(cmd, six.string_types) else cmd
                error_msg += "\n"
                error_msg += "".join(debug_stdout)
                s.communicate()
                s.stdout.close()
                raise subprocess.CalledProcessError(exitcode, error_msg)
            else:
                break
    s.communicate()
    s.stdout.close()
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
        if input_file.endswith((".gz")):
            scale = 10.0
        # bams can be compressed at different levels so use larger scale factor
        # to account for that potential
        elif input_file.endswith((".bam")):
            scale = 25.0
        else:
            scale = 20.0
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
