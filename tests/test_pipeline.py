"""Test individual components of the analysis pipeline.
"""
import os
import unittest

from bcbio.pipeline.run_info import _generate_lane

class RunInfoTest(unittest.TestCase):
    def test_lanename(self):
        """Generate lane names from supplied filenames.
        """
        assert _generate_lane(["s_1_sequence.txt"], 2) == "1"
        assert _generate_lane(["aname-sequence.fastq"], 2) == "aname"
        assert _generate_lane(["s_1_1-sequence.txt", "s_1_2-sequence.txt"], 2) == "1"
        assert _generate_lane(["one.txt", "two.txt"], 2) == "3"
