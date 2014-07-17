import os
import unittest
import shutil
import test_data
from bcbio.rnaseq import dexseq
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

    def tearDown(self):
        shutil.rmtree(self.out_dir)

