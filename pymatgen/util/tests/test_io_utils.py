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

from pymatgen.util.io_utils import FileLock, FileLockException

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')

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
