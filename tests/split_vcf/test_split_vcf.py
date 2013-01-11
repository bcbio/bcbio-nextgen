import yaml
import unittest
from bcbio.variation.split import split_vcf
import filecmp
import shutil
import os


class TestVcf(unittest.TestCase):

    def setUp(self):
        self.config_file = "tests/split_vcf/test_split_vcf.yaml"
        with open(self.config_file) as in_handle:
            self.config = yaml.load(in_handle)

    def test_splitvcf(self):
        in_file = self.config["input_split"]
        correct = self.config["correct_split"]
        out_files = split_vcf(in_file, self.config)
        self.assertTrue(all(map(filecmp.cmp, out_files, correct)))

        # cleanup
        data_dir = os.path.dirname(in_file)
        shutil.rmtree(os.path.join(data_dir, "split"))
        os.remove(in_file + ".gz")
        os.remove(in_file + ".gz.tbi")
