# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import os
import unittest

from pymatgen.util.io_utils import micro_pyawk
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "test_files")


class FuncTest(PymatgenTest):
    def test_micro_pyawk(self):
        filename = os.path.join(test_dir, "OUTCAR")
        data = []

        def f(x, y):
            data.append(y.group(1).strip())

        def f2(x, y):
            return y

        micro_pyawk(filename, [["POTCAR:(.*)", f2, f]])
        self.assertEqual(len(data), 6)


if __name__ == "__main__":
    unittest.main()
