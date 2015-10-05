"""Tests for the bcbio integration with the Common Workflow Language
"""
import os
import unittest
import subprocess

from nose.plugins.attrib import attr
import yaml

class CWLTest(unittest.TestCase):
    """ Run a simple CWL workflow.
    """
    def setUp(self):
        pass

    def test_1_cwl_nocontainer(self):
        cl = ["cwl-runner", "--verbose", "--no-container",
              "../cwl/bcbio2cwl.cwl", "../cwl/testinput-args.json"]

        subprocess.check_call(cl)
