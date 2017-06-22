# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import unittest
import os

from pymatgen.util.testing import PymatgenTest
from pymatgen.ext.wmm import get_from_wmm
from pymatgen.io.vasp.inputs import Incar

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class WmmTest(PymatgenTest):
    def test_get_wmm_kpoints(self):
        si = PymatgenTest.get_structure("Si")
        file_name = os.path.join(test_dir, 'INCAR')
        incar = Incar.from_file(file_name)
        kpoints = get_from_wmm(si, incar=incar)

if __name__ == '__main__':
    unittest.main()
