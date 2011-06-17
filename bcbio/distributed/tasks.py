"""Task definitions for the Celery message queue (http://celeryproject.org/).
"""
import time

from celery.task import task

from bcbio.pipeline import sample
from bcbio.pipeline import lane

@task
def test(x):
    print x
    time.sleep(5)
    return x

@task
def process_lane(*args):
    return lane.process_lane(*args)

@task
def process_alignment(*args):
    return lane.process_alignment(*args)

@task
def process_sample(*args):
    return sample.process_sample(*args)
