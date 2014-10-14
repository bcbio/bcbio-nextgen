import unittest
from bcbio.pipeline import run_info
from nose.plugins.attrib import attr

@attr("unit")
class TestSanityCheck(unittest.TestCase):
    sample = {"description": "Dummy"}

    def test_multiple_file_types(self):
        files = ["a.bam", "b.fastq"]
        self.assertRaises(ValueError, run_info._sanity_check_files, self.sample, files)

    def test_multiple_bam_files(self):
        files = ["a.bam", "b.bam"]
        self.assertRaises(ValueError, run_info._sanity_check_files, self.sample, files)

    def test_too_many_files(self):
        files = ["a.fastq", "b.fastq", "c.fastq"]
        self.assertRaises(ValueError, run_info._sanity_check_files, self.sample, files)

    def test_fastq_files_are_same(self):
        files = ["a.fastq", "a.fastq"]
        self.assertRaises(ValueError, run_info._sanity_check_files, self.sample, files)

    def test_bam_works(self):
        files = ["a.bam"]
        try:
            run_info._sanity_check_files(self.sample, files)
        except:
            self.fail()

    def test_fastq_works(self):
        files = ["a.fastq"]
        try:
            run_info._sanity_check_files(self.sample, files)
        except:
            self.fail()

    def test_paired_fastq_works(self):
        files = ["a.fastq", "b.fastq"]
        try:
            run_info._sanity_check_files(self.sample, files)
        except:
            self.fail()
