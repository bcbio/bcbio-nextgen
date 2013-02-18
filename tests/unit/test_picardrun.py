import yaml
import unittest
from bcbio.broad.picardrun import bed2interval
from bcbio.utils import replace_suffix
import os
from tempfile import NamedTemporaryFile
import filecmp


class TestBed2interval(unittest.TestCase):

    def setUp(self):
        self.config_file = "data/bed2interval/test_bed2interval.yaml"
        with open(self.config_file) as in_handle:
            self.config = yaml.load(in_handle)
        self.in_file = self.config["input"]
        self.bed = self.config["annotation"]["bed"]
        self.correct_file = self.config["correct"]

    def test_bed2interval(self):
        tmpfile = NamedTemporaryFile()
        out_file = bed2interval(self.in_file, self.bed,
                                out_file=tmpfile.name)
        self.assertTrue(filecmp.cmp(self.correct_file, out_file))
