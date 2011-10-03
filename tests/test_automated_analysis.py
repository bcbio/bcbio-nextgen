"""This directory is setup with configurations to run the main functional test.

It exercises a full analysis pipeline on a smaller subset of data.
"""
import os
import subprocess
import unittest
import shutil
import contextlib
import collections

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
        self.data_dir = os.path.join(os.path.dirname(__file__), "data", "automated")

    def _install_test_files(self, data_dir):
        """Download required sequence and reference files.
        """
        DlInfo = collections.namedtuple("DlInfo", "fname dirname version")
        download_data = [DlInfo("110106_FC70BUKAAXX.tar.gz", None, None),
                         DlInfo("genomes_automated_test.tar.gz", "genomes", 1),
                         DlInfo("110907_ERP000591.tar.gz", None, None),
                         DlInfo("100326_FC6107FAAXX.tar.gz", None, 1)]
        for dl in download_data:
            url = "http://chapmanb.s3.amazonaws.com/{fname}".format(fname=dl.fname)
            dirname = os.path.join(data_dir, os.pardir,
                                   dl.fname.replace(".tar.gz", "") if dl.dirname is None
                                   else dl.dirname)
            if os.path.exists(dirname) and dl.version is not None:
                version_file = os.path.join(dirname, "VERSION")
                is_old = True
                if os.path.exists(version_file):
                    with open(version_file) as in_handle:
                        version = int(in_handle.read())
                    is_old = version < dl.version
                if is_old:
                    shutil.rmtree(dirname)
            if not os.path.exists(dirname):
                self._download_to_dir(url, dirname)

    def _download_to_dir(self, url, dirname):
        print dirname
        cl = ["wget", url]
        subprocess.check_call(cl)
        cl = ["tar", "-xzvpf", os.path.basename(url)]
        subprocess.check_call(cl)
        os.rename(os.path.basename(dirname), dirname)
        os.remove(os.path.basename(url))

    def test_3_full_pipeline(self):
        """Run full automated analysis pipeline with multiplexing.
        """
        self._install_test_files(self.data_dir)
        with make_workdir():
            cl = ["automated_initial_analysis.py",
                  os.path.join(self.data_dir, "post_process.yaml"),
                  os.path.join(self.data_dir, os.pardir, "110106_FC70BUKAAXX"),
                  os.path.join(self.data_dir, "run_info.yaml")]
            subprocess.check_call(cl)

    def test_4_empty_fastq(self):
        """Handle analysis of empty fastq inputs from failed runs.
        """
        with make_workdir():
            cl = ["automated_initial_analysis.py",
                  os.path.join(self.data_dir, "post_process.yaml"),
                  os.path.join(self.data_dir, os.pardir, "110221_empty_FC12345AAXX"),
                  os.path.join(self.data_dir, "run_info-empty.yaml")]
            subprocess.check_call(cl)

    def test_2_rnaseq(self):
        """Run an RNA-seq analysis with TopHat and Cufflinks.
        """
        self._install_test_files(self.data_dir)
        with make_workdir():
            cl = ["automated_initial_analysis.py",
                  os.path.join(self.data_dir, "post_process.yaml"),
                  os.path.join(self.data_dir, os.pardir, "110907_ERP000591"),
                  os.path.join(self.data_dir, "run_info-rnaseq.yaml")]
            subprocess.check_call(cl)

    def test_1_variantcall(self):
        """Test variant calling with GATK pipeline.
        """
        self._install_test_files(self.data_dir)
        with make_workdir():
            cl = ["automated_initial_analysis.py",
                  os.path.join(self.data_dir, "post_process.yaml"),
                  os.path.join(self.data_dir, os.pardir, "100326_FC6107FAAXX"),
                  os.path.join(self.data_dir, "run_info-variantcall.yaml")]
            subprocess.check_call(cl)


