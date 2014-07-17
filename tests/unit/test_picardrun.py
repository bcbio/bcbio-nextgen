import unittest
import shutil
import os
from bcbio.broad.picardrun import bed2interval
from bcbio.utils import safe_makedir, file_exists
import test_data

from nose.plugins.attrib import attr

class TestBed2interval(unittest.TestCase):

    def setUp(self):
        self.out_dir = "test_picard"
        safe_makedir(self.out_dir)

    @attr("unit")
    def test_bed2interval(self):
        test_file = os.path.join(self.out_dir, "test.interval")
        bam_file = test_data.BAM_FILE
        bed_file = test_data.BED_FILE
        out_file = bed2interval(bam_file, bed_file, out_file=test_file)
        self.assertTrue(file_exists(out_file))

    def tearDown(self):
        shutil.rmtree(self.out_dir)
