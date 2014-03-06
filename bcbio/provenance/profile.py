"""Profiling of system resources (CPU, memory, disk, filesystem IO) during pipeline runs.
"""
import contextlib
import os

from bcbio import utils
from bcbio.log import logger

@contextlib.contextmanager
def report(label, dirs):
    """Run reporting metrics to prepare reports of resource usage.
    """
    logger.info("Timing: %s" % label)
    profile_dir = utils.safe_makedir(os.path.join(dirs["work"], "profile"))
    # Prepare and start profiling scripts
    try:
        yield None
    finally:
        # Cleanup and stop profiling
        pass
