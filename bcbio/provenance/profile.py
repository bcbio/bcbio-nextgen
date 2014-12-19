"""Profiling of system resources (CPU, memory, disk, filesystem IO) during
pipeline runs.
"""
import contextlib

from bcbio.log import logger

@contextlib.contextmanager
def report(label, dirs):
    """Log timing information for later graphing of resource usage."""
    logger.info("Timing: %s" % label)
    yield None
