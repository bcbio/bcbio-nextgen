import os
import unittest
import shutil
from bcbio.rnaseq import cufflinks
from bcbio.utils import file_exists, safe_makedir
from nose.plugins.attrib import attr

DATA_DIR = os.path.join(os.path.dirname(__file__), "bcbio-nextgen-test-data", "data")

class TestCufflinks(unittest.TestCase):
    merged_gtf = os.path.join(DATA_DIR, "cufflinks", "merged.gtf")
    ref_gtf = os.path.join(DATA_DIR, "cufflinks", "ref-transcripts.gtf")
    out_dir = "cufflinks-test"

    def setUp(self):
        safe_makedir(self.out_dir)

    @attr("unit")
    def test_cufflinks_clean(self):
        clean_fn = os.path.join(self.out_dir, "clean.gtf")
        dirty_fn = os.path.join(self.out_dir, "dirty.gtf")
        clean, dirty = cufflinks.clean_assembly(self.merged_gtf, clean_fn,
                                                dirty_fn)
        assert(file_exists(clean))
        assert(os.path.exists(dirty))

    def tearDown(self):
        shutil.rmtree(self.out_dir)
