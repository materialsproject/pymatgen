#!/usr/bin/env python

"""
Created on Nov 14, 2012
"""


__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Nov 14, 2012"

import unittest
import os

from pymatgen.util.io_utils import reverse_readline, FileLock, FileLockException

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class BackwardsReaderTest(unittest.TestCase):
    NUMLINES = 3000

    def test_reverse_readline(self):
        """
        We are making sure a file containing line numbers is read in reverse
        order, i.e. the first line that is read corresponds to the last line.
        number
        """
        with open(os.path.join(test_dir, "three_thousand_lines.txt")) as f:
            for idx, line in enumerate(reverse_readline(f)):
                self.assertEqual(int(line), self.NUMLINES - idx,
                                 "read_backwards read {} whereas it should "
                                 "have read {}".format(
                                     int(line), self.NUMLINES - idx))

    def test_empty_file(self):
        """
        make sure an empty file does not throw an error when reverse_readline
        is called this was a problem with an earlier implementation
        """
        with open(os.path.join(test_dir, "empty_file.txt")) as f:
            for idx, line in enumerate(reverse_readline(f)):
                raise ValueError("an empty file is being read!")


class FileLockTest(unittest.TestCase):

    def setUp(self):
        self.file_name = "__lock__"
        self.lock = FileLock(self.file_name, timeout=1)
        self.lock.acquire()

    def test_raise(self):
        with self.assertRaises(FileLockException):
            new_lock = FileLock(self.file_name, timeout=1)
            new_lock.acquire()

    def tearDown(self):
        self.lock.release()

if __name__ == "__main__":
    unittest.main()
