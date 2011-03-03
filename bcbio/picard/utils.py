"""Various utilities to handle working with Picard and GATK.
"""
import os
import tempfile
import shutil
import contextlib

@contextlib.contextmanager
def curdir_tmpdir():
    """Context manager to create and remove a temporary directory.
    """
    tmp_dir_base = os.path.join(os.getcwd(), "tmp")
    if not os.path.exists(tmp_dir_base):
        # we could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(tmp_dir_base)
        except OSError:
            pass
    tmp_dir = tempfile.mkdtemp(dir=tmp_dir_base)
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    try :
        yield tmp_dir
    finally :
        shutil.rmtree(tmp_dir)

@contextlib.contextmanager
def chdir(new_dir):
    """Context manager to temporarily change to a new directory.

    http://lucentbeing.com/blog/context-managers-and-the-with-statement-in-python/
    """
    cur_dir = os.getcwd()
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    os.chdir(new_dir)
    try :
        yield
    finally :
        os.chdir(cur_dir)

@contextlib.contextmanager
def tmpfile(*args, **kwargs):
    """Make a tempfile, safely cleaning up file descriptors on completion.
    """
    (fd, fname) = tempfile.mkstemp(*args, **kwargs)
    try:
        yield fname
    finally:
        os.close(fd)
        if os.path.exists(fname):
            os.remove(fname)
