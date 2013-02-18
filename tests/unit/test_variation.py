import yaml
import unittest
from bcbio.variation.split import split_vcf
import filecmp
import shutil
import os


class TestVcf(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), "data")

    def test_splitvcf(self):
        data_subdir = os.path.join(self.data_dir, "split_vcf")
        config_file = os.path.join(data_subdir, "test_split_vcf.yaml")
        with open(config_file) as in_handle:
            config = yaml.load(in_handle)
        in_file = config["input_split"]
        correct_files = config["correct_split"]
        out_files = split_vcf(in_file, config)
        self.assertTrue(all(map(filecmp.cmp, out_files, correct_files)))

        # cleanup
        data_dir = os.path.dirname(in_file)
        shutil.rmtree(os.path.join(data_dir, "split"))
        os.remove(in_file + ".gz")
        os.remove(in_file + ".gz.tbi")
