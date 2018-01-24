# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
import unittest

from pymatgen.io.lammps.sets import LammpsInputSet

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files", "lammps")


class TestLammpsInputSet(unittest.TestCase):

    def setUp(self):
        template_file = os.path.join(test_dir, "in.peptide.template")
        data_file = os.path.join(test_dir, "data.peptide")
        self.data_filename = "test_data.peptide"
        self.input_filename = "test_input.peptide"
        self.settings = {
            "pair_style": "lj/charmm/coul/long 8.0 10.0 10.0",
            "kspace_style": "pppm 0.0001",
            "fix_1": "1 all nvt temp 275.0 275.0 100.0 tchain 1",
            "fix_2": "2 all shake 0.0001 10 100 b 4 6 8 10 12 14 18 a 31"
        }
        self.lammps_input_set = LammpsInputSet.from_file(
            "test", template_file, self.settings, lammps_data=data_file,
            data_filename=self.data_filename)

    def test_input(self):
        self.assertEqual(self.lammps_input_set.lammps_input.settings["data_file"],
                         self.data_filename)
        for k, v in self.settings.items():
            self.assertEqual(self.lammps_input_set.lammps_input.settings[k], v)

    def test_write_input_set(self):
        self.lammps_input_set.write_input(self.input_filename)
        self.assertTrue(os.path.exists(self.input_filename))
        self.assertTrue(os.path.exists(self.data_filename))
        os.remove(self.input_filename)
        os.remove(self.data_filename)
        # now change both input and data filenames
        self.lammps_input_set.write_input("xxxx.input", "yyy.data")
        self.assertTrue(os.path.exists("xxxx.input"))
        self.assertTrue(os.path.exists("yyy.data"))
        os.remove("xxxx.input")
        os.remove("yyy.data")


if __name__ == "__main__":
    unittest.main()
