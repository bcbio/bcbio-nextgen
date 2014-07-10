import unittest
import tempfile
from bcbio.bam.fastq import groom
from bcbio.utils import locate, file_exists
import os
import tempfile

from nose.plugins.attrib import attr

class Fastq(unittest.TestCase):

    def setUp(self):
        self.root_dir = os.path.join(os.path.dirname(__file__), "data/fastq/")

    @attr("unit")
    def test_groom(self):
        illumina_dir = os.path.join(self.root_dir, "illumina")
        test_data = locate("*.fastq", illumina_dir)
        self.assertTrue(not test_data == [])
        sanger_dir = tempfile.mkdtemp()
        out_files = [groom(x, in_qual="fastq-illumina", out_dir=sanger_dir) for
                     x in test_data]
        self.assertTrue(all(map(file_exists, out_files)))
