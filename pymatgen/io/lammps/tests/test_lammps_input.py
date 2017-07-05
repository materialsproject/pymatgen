# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
import unittest

from pymatgen.io.lammps.input import DictLammpsInput

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files", "lammps")


class TestLammpsInput(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.json_input = DictLammpsInput.from_file(
            "NVT", os.path.join(test_dir, "NVT.json"),
            data_filename=os.path.join(test_dir, "nvt.data"),
            is_forcefield=True)
        cls.lammps_input = DictLammpsInput.from_file("peptide", os.path.join(test_dir, "in.peptide"),
                                               lammps_data=os.path.join(test_dir, "data.peptide"),
                                               data_filename="data.peptide",
                                               is_forcefield=True)

    def test_string_rep(self):
        self.json_input.config_dict["read_data"] = "nvt.data"
        with open(os.path.join(test_dir, "nvt.inp")) as f:
            for l1, l2 in zip(str(self.json_input).split("\n"),
                              f.readlines()):
                self.assertEqual(l1.strip(), l2.strip())

    def test_from_file(self):
        # test lammps style input
        self.assertEqual(self.lammps_input.config_dict["pair_style"], "lj/charmm/coul/long 8.0 10.0 10.0")
        self.assertEqual(self.lammps_input.config_dict["fix"], ["1 all nvt temp 275.0 275.0 100.0 tchain 1",
                                                               "2 all shake 0.0001 10 100 b 4 6 8 10 12 14 18 a 31"])
        test_box = self.lammps_input.as_dict()["lammps_data"]["box_size"][0]
        box = [36.840194, 64.211560]
        self.assertEqual(box, test_box)
        # test JSON style input
        self.assertEqual(self.json_input.config_dict["kspace_style"], "ewald 0.0001")

if __name__ == "__main__":
    unittest.main()
