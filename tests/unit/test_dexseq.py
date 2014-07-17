import os
import unittest
import shutil
import test_data
from bcbio.rnaseq import dexseq
from bcbio.rnaseq import count
from bcbio.utils import file_exists, safe_makedir
from nose.plugins.attrib import attr

@attr("unit")
class TestDEXSeqCount(unittest.TestCase):

    def setUp(self):
        self.out_dir = "dexseq-test"
        safe_makedir(self.out_dir)

    def test_dexseq_count(self):
        out_file = os.path.join(self.out_dir, "dexseq-counts.txt")
        bam_file = test_data.BAM_FILE
        dexseq_gff = test_data.DEXSEQ_GFF
        stranded = "unstranded"
        out_file = dexseq.run_count(bam_file, dexseq_gff, stranded, out_file)
        self.assertTrue(file_exists(out_file))

    def test_dexseq_combine(self):
        count_files = test_data.DEXSEQ_COUNT_FILES
        test_file = os.path.join(self.out_dir, "dexseq-combined.txt")
        out_file = count.combine_count_files(count_files, out_file=test_file,
                                             ext=".txt")
        self.assertTrue(file_exists(out_file))

    def tearDown(self):
        pass
        shutil.rmtree(self.out_dir)
