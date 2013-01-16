import unittest
import tempfile
from bcbio.utils import transform_to, filter_to, memoize_outfile
import os


class Utils(unittest.TestCase):

    def setUp(self):
        self.out_dir = tempfile.mkdtemp()

    def test_transform_to(self):
        in_file = "test.sam"
        out_file = "test.bam"

        @transform_to(".bam")
        def f(in_file, out_dir=None, out_file=None):
            return out_file

        self.assertTrue(f(in_file) == out_file)
        self.assertTrue(f(in_file, out_dir=self.out_dir) ==
                        os.path.join(self.out_dir, out_file))

    def test_filter_to(self):
        in_file = "test.sam"
        out_file = "test_sorted.sam"

        @filter_to("sorted")
        def f(in_file, out_dir=None, out_file=None):
            return out_file

        self.assertTrue(f(in_file) == out_file)
        self.assertTrue(f(in_file, out_dir=self.out_dir) ==
                        os.path.join(self.out_dir, out_file))

    def test_memoize_outfile(self):
        out_file = tempfile.NamedTemporaryFile(dir=self.out_dir)

        @memoize_outfile
        def f(out_file):
            return "not_memoized"

        # file is empty so the function should run
        print f(out_file=out_file.name)
        self.assertTrue(f(out_file=out_file.name) == "not_memoized")
        # file is non-empty to function should not run
        out_file.write("foo")
        out_file.flush()
        self.assertTrue(f(out_file=out_file.name) == out_file.name)
