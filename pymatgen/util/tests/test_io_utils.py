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
import json

from pymatgen.util.testing import PymatgenTest
from pymatgen.util.io_utils import FileLock, FileLockException, clean_json

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class FuncTest(PymatgenTest):

    def test_clean_json(self):
        #clean_json should have no effect on None types.
        d = {"hello": 1, "world": None}
        clean = clean_json(d)
        self.assertIsNone(clean["world"])
        self.assertEqual(json.dumps(d), json.dumps(clean))

        d = {"hello": self.get_test_structure()}
        self.assertRaises(TypeError, json.dumps, d)
        clean = clean_json(d)
        self.assertIsInstance(clean["hello"], basestring)
        clean_strict = clean_json(d, strict=True)
        self.assertIsInstance(clean_strict["hello"], dict)


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
