# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest

import requests

from pymatgen.ext.jhu import get_kpoints
from pymatgen.io.vasp.inputs import Incar
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.util.testing import PymatgenTest

__author__ = "Joseph Montoya"
__copyright__ = "Copyright 2017, The Materials Project"
__maintainer__ = "Joseph Montoya"
__email__ = "montoyjh@lbl.gov"
__date__ = "June 22, 2017"


website_is_up = requests.get("http://muellergroup.jhu.edu:8080").status_code == 200


@unittest.skipIf(True, "This code is way too buggy to be tested.")
class JhuTest(PymatgenTest):
    _multiprocess_shared_ = True

    def test_get_kpoints(self):
        si = PymatgenTest.get_structure("Si")
        input_set = MPRelaxSet(si)
        kpoints = get_kpoints(si, incar=input_set.incar)


if __name__ == "__main__":
    unittest.main()
