import unittest
import tempfile
from bcbio.utils import (transform_to, filter_to, memoize_outfile,
                         replace_suffix, file_exists, append_stem)
from bcbio import utils
import os
import time


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

    def test_memoize_outfile_ext(self):
        temp_file = tempfile.NamedTemporaryFile(dir=self.out_dir,
                                                suffix=".sam")
        ext = ".bam"
        word = __name__

        # test by changing the extension
        @utils.memoize_outfile(ext=ext)
        def f(in_file, word, out_file=None):
            return self._write_word(out_file, word)
        self._run_calls_no_dir(f, temp_file, word)

    def test_memoize_outfile_stem(self):
        temp_file = tempfile.NamedTemporaryFile(dir=self.out_dir,
                                                suffix=".sam")
        stem = "stem"
        word = __name__

        @utils.memoize_outfile(stem=stem)
        def f(in_file, word, out_file=None):
            return self._write_word(out_file, word)
        self._run_calls_no_dir(f, temp_file, word)

    def test_memoize_outfile_stem_with_dir(self):
        temp_file = tempfile.NamedTemporaryFile(dir=self.out_dir,
                                                suffix=".sam")
        temp_dir = tempfile.mkdtemp()
        stem = "stem"
        word = __name__

        @utils.memoize_outfile(stem=stem)
        def f(in_file, word, out_dir=None, out_file=None):
            return self._write_word(out_file, word)

        self._run_calls_with_dir(f, temp_file, word, out_dir=temp_dir)

    def test_memoize_outfile_ext_with_dir(self):
        temp_file = tempfile.NamedTemporaryFile(dir=self.out_dir,
                                                suffix=".sam")
        temp_dir = tempfile.mkdtemp()
        ext = ".bam"
        word = __name__

        @utils.memoize_outfile(ext=ext)
        def f(in_file, word, out_dir=None, out_file=None):
            return self._write_word(out_file, word)

        self._run_calls_with_dir(f, temp_file, word, out_dir=temp_dir)

    def _run_calls_with_dir(self, f, temp_file, word, out_dir=None):
        first_call = f(temp_file.name, word, out_dir=out_dir)
        second_call = f(temp_file.name, word, out_dir=out_dir)
        self.assertTrue(self._read_word(first_call) == word)
        self.assertTrue(self._read_word(second_call) == word)

    def _run_calls_no_dir(self, f, temp_file, word):
        first_call = f(temp_file.name, word)
        second_call = f(temp_file.name, word)
        self.assertTrue(self._read_word(first_call) == word)
        self.assertTrue(self._read_word(second_call) == word)

    def _read_word(self, in_file):
        with open(in_file) as in_handle:
            return in_handle.read()

    def _write_word(self, out_file, word):
        with open(out_file, "w") as out_handle:
            out_handle.write(word)
            out_handle.flush()
        return out_file
