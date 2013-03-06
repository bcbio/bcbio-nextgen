import os
import unittest
from bcbio.rnaseq import count
from bcbio.utils import safe_makedir
import tempfile
import stat
import shutil


class TestHtseqCount(unittest.TestCase):
    cur_dir = os.path.dirname(__file__)
    organism_dir = os.path.join(cur_dir, "data", "organisms", "mouse")
    data_dir = os.path.join(cur_dir, "data", "count", "test_data")
    correct_dir = os.path.join(cur_dir, "data", "count", "correct")
    out_dir = os.path.join(cur_dir, "htseq-test")

    def setUp(self):
        self.in_bam = os.path.join(self.data_dir, "test.bam")
        self.in_gtf = os.path.join(self.data_dir, "test.gtf")
        self.count_file = os.path.join(self.data_dir, "test.count")
        self.correct_file = os.path.join(self.correct_dir, "correct.count")
        safe_makedir(self.out_dir)


    def test_is_countfile_correct(self):
        test_file = os.path.join(self.data_dir, "test.count")
        self.assertTrue(count.is_countfile(test_file))

    def test_is_countfile_not_correct(self):
        test_file = os.path.join(self.organism_dir, "mouse.gtf")
        self.assertFalse(count.is_countfile(test_file))

    def test_htseq_is_installed_in_path(self):
        self.assertTrue(count.htseq_is_installed({"config": {}}))

    def test_htseq_is_installed_in_resource(self):
        orig_path = os.environ['PATH']
        os.environ['PATH'] = ""
        faux_htseq_count = tempfile.NamedTemporaryFile()
        os.chmod(faux_htseq_count.name, stat.S_IEXEC)
        config = {"config": {"resources": {"htseq-count":
                                           {"cmd": faux_htseq_count.name}}}}
        is_installed = count.htseq_is_installed(config)
        os.environ['PATH'] = orig_path
        self.assertTrue(is_installed)

    def test_htseq_count(self):
        data = {"work_bam": self.in_bam,
                "config": {"algorithm": {"transcripts": self.in_gtf},
                           "dirs": {"work": self.out_dir}}}
        out_file = count.htseq_count(data)
        self.assertTrue(count.is_countfile(out_file))

    def tearDown(self):
        shutil.rmtree(self.out_dir)
