"""Task definitions for the Celery message queue (http://celeryproject.org/).
"""
import time

from celery.task import task

from bcbio.pipeline import sample
from bcbio.pipeline import lane

# Global configuration for tasks in the main celeryconfig module
import celeryconfig

@task
def analyze_and_upload(*args):
    """Run full analysis and upload results to Galaxy instance.

    Workers need to run on the machine with Galaxy, but can be
    distributed to multiple nodes if needed.
    """
    print celeryconfig.BCBIO_CONFIG_FILE

@task
def process_lane(*args):
    return lane.process_lane(*args)

@task
def process_alignment(*args):
    return lane.process_alignment(*args)

@task
def process_sample(*args):
    return sample.process_sample(*args)

@task
def test(x):
    print x
    time.sleep(5)
    return x
