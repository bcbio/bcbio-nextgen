"""This directory is setup with configurations to run the main functional test.

It exercises a full analysis pipeline on a smaller subset of data.
"""
import os
import subprocess
import unittest
import shutil
import contextlib

@contextlib.contextmanager
def make_workdir():
    dirname = os.path.join(os.path.dirname(__file__), "test_automated_output")
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
        self.data_dir = os.path.join(os.pardir, "data", "automated")

    def _install_test_files(self, data_dir):
        """Download required sequence and reference files.
        """
        read_url = "http://chapmanb.s3.amazonaws.com/110106_FC70BUKAAXX.tar.gz"
        read_dir = os.path.join(data_dir, os.pardir, "110106_FC70BUKAAXX")
        if not os.path.exists(read_dir):
            self._download_to_dir(read_url, read_dir)

        genome_url = "http://chapmanb.s3.amazonaws.com/genomes_automated_test.tar.gz"
        genome_dir = os.path.join(data_dir, os.pardir, "genomes")
        if not os.path.exists(genome_dir):
            self._download_to_dir(genome_url, genome_dir)

        rnaseq_url = "http://chapmanb.s3.amazonaws.com/110907_ERP000591.tar.gz"
        rnaseq_dir = os.path.join(data_dir, os.pardir, "110907_ERP000591")
        if not os.path.exists(rnaseq_dir):
            self._download_to_dir(rnaseq_url, rnaseq_dir)

    def _download_to_dir(self, url, dirname):
        cl = ["wget", url]
        subprocess.check_call(cl)
        cl = ["tar", "-xzvpf", os.path.basename(url)]
        subprocess.check_call(cl)
        os.rename(os.path.basename(dirname), dirname)

    def test_run_full_pipeline(self):
        """Run full automated analysis pipeline.
        """
        with make_workdir():
            self._install_test_files(self.data_dir)
            cl = ["automated_initial_analysis.py",
                  os.path.join(self.data_dir, "post_process.yaml"),
                  os.path.join(self.data_dir, os.pardir, "110106_FC70BUKAAXX"),
                  os.path.join(self.data_dir, "run_info.yaml")]
            subprocess.check_call(cl)

    def test_empty_fastq(self):
        """Handle analysis of empty fastq inputs from failed runs.
        """
        with make_workdir():
            cl = ["automated_initial_analysis.py",
                  os.path.join(self.data_dir, "post_process.yaml"),
                  os.path.join(self.data_dir, os.pardir, "110221_empty_FC12345AAXX"),
                  os.path.join(self.data_dir, "run_info-empty.yaml")]
            subprocess.check_call(cl)

    def test_rnaseq(self):
        """Run an RNA-seq analysis with TopHat and Cufflinks.
        """
        with make_workdir():
            self._install_test_files(self.data_dir)
            cl = ["automated_initial_analysis.py",
                  os.path.join(self.data_dir, "post_process.yaml"),
                  os.path.join(self.data_dir, os.pardir, "110907_ERP000591"),
                  os.path.join(self.data_dir, "run_info-rnaseq.yaml")]
            subprocess.check_call(cl)


