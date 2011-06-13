"""Task definitions for the Celery message queue (http://celeryproject.org/).
"""
import time

from celery.task import task

@task
def test(x):
    print x
    time.sleep(5)
    return x
