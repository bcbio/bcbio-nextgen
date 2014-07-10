from bcbio.bam.trim import trim_read_through
from bcbio.utils import append_stem
import yaml
import unittest
import filecmp
import os
import shutil

CONFIG_FILE = "data/trim_read_through/test_trim_read_through.yaml"
CORRECT_DIR = "data/trim_read_through/correct"

from nose.plugins.attrib import attr

class TestTrimReadThrough(unittest.TestCase):

    def setUp(self):
        with open(CONFIG_FILE) as in_handle:
            self.config = yaml.load(in_handle)
        self.root_work = os.path.dirname(self.config["dir"]["work"])

    def _find_length_filter_correct(self, out_file):
        correct_dir = os.path.join(os.path.dirname(out_file), "correct",
                                   "length_filter")
        correct_file = os.path.join(correct_dir, os.path.basename(out_file))
        return correct_file

    def _trim_single_correct(self, out_file):
        correct_dir = os.path.join(CORRECT_DIR, "trim_single")
        return os.path.join(correct_dir, os.path.basename(out_file))

    def _trim_paired_correct(self, out_file):
        correct_dir = os.path.join(CORRECT_DIR, "trim_paired")
        return os.path.join(correct_dir, os.path.basename(out_file))

    @attr("unit-broken")
    def test_pairedend(self):
        paired = self.config["input_paired"]
        out_files = trim_read_through(paired, self.config["dir"], self.config)
        correct_files = map(self._trim_paired_correct, out_files)
        self.assertTrue(all(map(filecmp.cmp, correct_files, out_files)))
        shutil.rmtree(self.root_work)

    @attr("unit-broken")
    def test_single(self):
        single = self.config["input_single"]
        out_file = trim_read_through(single, self.config["dir"], self.config)
        correct_file = self._trim_single_correct(out_file)
        self.assertTrue(filecmp.cmp(correct_file, out_file))
        shutil.rmtree(self.root_work)

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestTrimReadThrough)
    unittest.TextTestRunner(verbosity=2).run(suite)
