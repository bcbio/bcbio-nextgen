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

        @filter_to("_sorted")
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
        stem = "_stem"
        word = __name__

        @utils.memoize_outfile(stem=stem)
        def f(in_file, word, out_file=None):
            return self._write_word(out_file, word)
        self._run_calls_no_dir(f, temp_file, word)

    def test_memoize_outfile_stem_with_dir(self):
        temp_file = tempfile.NamedTemporaryFile(dir=self.out_dir,
                                                suffix=".sam")
        temp_dir = tempfile.mkdtemp()
        stem = "_stem"
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

    def test_replace_suffix_of_string(self):
        test_string = "/string/test/foo.txt"
        correct = "/string/test/foo.bar"
        out_string = utils.replace_suffix(test_string, ".bar")
        self.assertEquals(correct, out_string)

    def test_replace_suffix_of_list(self):
        test_list = ["/list/test/foo.txt", "/list/test/foobar.txt"]
        correct = ["/list/test/foo.bar", "/list/test/foobar.bar"]
        out_list = utils.replace_suffix(test_list, ".bar")
        for c, o in zip(correct, out_list):
            self.assertEquals(c, o)

    def test_append_stem_of_string(self):
        test_string = "/string/test/foo.txt"
        correct = "/string/test/foo_bar.txt"
        out_string = utils.append_stem(test_string, "_bar")
        self.assertEquals(correct, out_string)

    def test_append_stem_of_list(self):
        test_list = ["/list/test/foo.txt", "/list/test/foobar.txt"]
        correct = ["/list/test/foo_bar.txt", "/list/test/foobar_bar.txt"]
        out_list = utils.append_stem(test_list, "_bar")
        for c, o in zip(correct, out_list):
            self.assertEquals(c, o)

    def test_replace_directory_of_string(self):
        test_string = "/string/test/foo.txt"
        correct = "/new/dir/foo.txt"
        out_string = utils.replace_directory(test_string, "/new/dir")
        self.assertEquals(correct, out_string)

    def test_replace_directory_of_list(self):
        test_list = ["/list/test/bar.txt", "/list/test/foobar.txt"]
        correct = ["/new/dir/bar.txt", "/new/dir/foobar.txt"]
        out_list = utils.replace_directory(test_list, "/new/dir")
        for c, o in zip(correct, out_list):
            self.assertEquals(c, o)

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
