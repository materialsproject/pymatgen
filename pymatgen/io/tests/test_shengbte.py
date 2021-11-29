# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
import unittest

from pymatgen.io.shengbte import Control
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "shengbte")

this_dir = os.path.dirname(os.path.abspath(__file__))

try:
    import f90nml
except ImportError:
    f90nml = None


class TestShengBTE(PymatgenTest):
    def setUp(self):
        self.filename = os.path.join(test_dir, "CONTROL-CSLD_Si")
        self.test_dict = {
            "nelements": 1,
            "natoms": 2,
            "ngrid": [25, 25, 25],
            "norientations": 0,
            "lfactor": 0.1,
            "lattvec": [
                [0.0, 2.734363999, 2.734363999],
                [2.734363999, 0.0, 2.734363999],
                [2.734363999, 2.734363999, 0.0],
            ],
            "elements": "Si",
            "types": [1, 1],
            "positions": [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
            "scell": [5, 5, 5],
            "t": 500,
            "scalebroad": 0.5,
            "isotopes": False,
            "onlyharmonic": False,
            "nonanalytic": False,
            "nanowires": False,
        }

    @unittest.skipIf(f90nml is None, "No f90nml")
    def test_from_file(self):
        io = Control.from_file(self.filename)
        self.assertEqual(io["nelements"], 1)
        self.assertEqual(io["natoms"], 2)
        self.assertArrayEqual(io["ngrid"], [25, 25, 25])
        self.assertEqual(io["norientations"], 0)
        self.assertEqual(io["lfactor"], 0.1)
        self.assertEqual(io["lattvec"][0], [0.0, 2.734363999, 2.734363999])
        self.assertEqual(io["lattvec"][1], [2.734363999, 0.0, 2.734363999])
        self.assertEqual(io["lattvec"][2], [2.734363999, 2.734363999, 0.0])
        self.assertIsInstance(io["elements"], (list, str))
        if isinstance(io["elements"], list):
            all_strings = all(isinstance(item, str) for item in io["elements"])
            self.assertTrue(all_strings)
        self.assertIsInstance(io["types"], (list, int))
        if isinstance(io["types"], list):
            all_ints = all(isinstance(item, int) for item in io["types"])
            self.assertTrue(all_ints)
        self.assertArrayEqual(io["positions"], [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]])
        self.assertArrayEqual(io["scell"], [5, 5, 5])
        self.assertEqual(io["t"], 500)
        self.assertEqual(io["scalebroad"], 0.5)
        self.assertFalse(io["isotopes"])
        self.assertFalse(io["onlyharmonic"])
        self.assertFalse(io["nonanalytic"])
        self.assertFalse(io["nanowires"])

        if os.path.exists(os.path.join(test_dir, "test_control")):
            os.remove(os.path.join(test_dir, "test_control"))
        io.to_file(filename=os.path.join(test_dir, "test_control"))

        with open(os.path.join(test_dir, "test_control")) as file:
            test_string = file.read()
        with open(os.path.join(test_dir, "CONTROL-CSLD_Si")) as reference_file:
            reference_string = reference_file.read()
        self.assertMultiLineEqual(test_string, reference_string)
        os.remove(os.path.join(test_dir, "test_control"))

    @unittest.skipIf(f90nml is None, "No f90nml")
    def test_from_dict(self):
        io = Control.from_dict(self.test_dict)
        if os.path.exists(os.path.join(test_dir, "test_control")):
            os.remove(os.path.join(test_dir, "test_control"))
        io.to_file(filename=os.path.join(test_dir, "test_control"))
        with open(os.path.join(test_dir, "test_control")) as file:
            test_string = file.read()
        with open(os.path.join(test_dir, "CONTROL-CSLD_Si")) as reference_file:
            reference_string = reference_file.read()
        self.assertMultiLineEqual(test_string, reference_string)
        os.remove(os.path.join(test_dir, "test_control"))

    @unittest.skipIf(f90nml is None, "No f90nml")
    def test_MSONable_implementation(self):
        # tests as dict and from dict methods
        Controlinfromfile = Control.from_file(self.filename)
        newControlin = Control.from_dict(Controlinfromfile.as_dict())
        self.assertDictEqual(newControlin, Controlinfromfile)
        newControlin.to_json()
