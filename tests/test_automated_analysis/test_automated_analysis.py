"""This directory is setup with configurations to run the main functional test.

It exercises a full analysis pipeline on a smaller subset of data.
"""
import os
import subprocess
import unittest
import shutil
import contextlib

@contextlib.contextmanager
def test_workdir():
    dirname = os.path.join(os.path.dirname(__file__), "test_output")
    if os.path.exists(dirname):
        shutil.rmtree(dirname)
    os.makedirs(dirname)
    orig_dir = os.getcwd()
    try:
        os.chdir(dirname)
        yield
    finally:
        os.chdir(orig_dir)

class AutomatedAnalysisTest(unittest.TestCase):
    """Setup a full automated analysis and run the pipeline.
    """
    def setUp(self):
        """Download required sequence and reference files.

        XXX Setup repository for these with automated download if not present.
        """
        pass

    def test_run_full_pipeline(self):
        """Run full automated analysis pipeline.
        """
        with test_workdir():
            cl = ["automated_initial_analysis.py",
                  os.path.join(os.pardir, "post_process.yaml"),
                  os.path.join(os.pardir, os.pardir, "data", "110106_FC70BUKAAXX"),
                  os.path.join(os.pardir, "run_info.yaml")]
            subprocess.check_call(cl)

