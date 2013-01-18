import unittest
import tempfile
from bcbio.fastq import groom
from bcbio.utils import locate, file_exists
import os
import tempfile


class Fastq(unittest.TestCase):

    def setUp(self):
        self.root_dir = os.path.join(__file__, "data/fastq/")

    def test_groom(self):
        illumina_dir = os.path.join(self.root_dir, "illumina")
        test_data = locate("*.fastq", illumina_dir)
        sanger_dir = tempfile.mkdtemp()
        out_files = [groom(x, in_qual="fastq-illumina", out_dir=sanger_dir) for
                     x in test_data]
        self.assertTrue(all(map(file_exists, out_files)))
