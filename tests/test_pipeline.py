"""Test individual components of the analysis pipeline.
"""
import os
import unittest

from nose.plugins.attrib import attr

from bcbio.pipeline.config_utils import load_config
from bcbio.pipeline.run_info import get_run_info
from bcbio.provenance import programs

class RunInfoTest(unittest.TestCase):
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), "data")

    @attr(speed=1)
    def test_run_info_combine(self):
        """Combine multiple lanes in a test run into a single combined lane.
        """
        run_info_yaml = os.path.join(self.data_dir, "run_info-alternatives.yaml")
        _, _, run_info = get_run_info("", {}, run_info_yaml)
        assert len(run_info["details"]) == 2
        assert len(run_info["details"][0]) == 3
        x1, x2, x3 = run_info["details"][0]
        assert x1["description"] == "1: BC1"
        assert x2["description"] == "1: BC2"
        assert x3["genome_build"] == "mm9"
        x1 = run_info["details"][1][0]
        assert x1["barcode_id"] is None

    @attr(speed=1)
    def test_programs(self):
        """Identify programs and versions used in analysis.
        """
        config = load_config(os.path.join(self.data_dir, "automated",
                                          "post_process-sample.yaml"))
        print programs.get_versions(config)
