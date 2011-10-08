"""Test individual components of the analysis pipeline.
"""
import os
import unittest

from bcbio.pipeline.run_info import _generate_lane, get_run_info

class RunInfoTest(unittest.TestCase):
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), "data")

    def test_lanename(self):
        """Generate lane names from supplied filenames.
        """
        assert _generate_lane(["s_1_sequence.txt"], 2) == "1"
        assert _generate_lane(["aname-sequence.fastq"], 2) == "aname"
        assert _generate_lane(["s_1_1-sequence.txt", "s_1_2-sequence.txt"], 2) == "1"
        assert _generate_lane(["one.txt", "two.txt"], 2) == "3"

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
