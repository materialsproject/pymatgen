#!/usr/bin/env python
import unittest

from pymatgen.util.filelock import FileLock

class FileLockTest(unittest.TestCase):

    def setUp(self):
        self.file_name = "__lock__"
        self.lock = FileLock(self.file_name, timeout=1) 
        self.lock.acquire()

    def test_raise(self):
        with self.assertRaises(FileLock.Exception):
            new_lock = FileLock(self.file_name, timeout=1) 
            new_lock.acquire()

    def tearDown(self):
        self.lock.release()

if __name__ == "__main__":
    unittest.main()


