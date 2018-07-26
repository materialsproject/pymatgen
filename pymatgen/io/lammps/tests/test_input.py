# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import unittest
import os
import re
import filecmp
import shutil

import pandas as pd
from pymatgen.io.lammps.data import LammpsData

from pymatgen.io.lammps.input import write_lammps_inputs


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files", "lammps")


class FuncTest(unittest.TestCase):

    def test_write_lammps_inputs(self):
        # script template
        with open(os.path.join(test_dir, "kappa.txt")) as f:
            kappa_template = f.read()
        kappa_settings = {"method": "heat"}
        write_lammps_inputs(output_dir="heat", script_template=kappa_template,
                            settings=kappa_settings)
        with open(os.path.join("heat", "in.lammps")) as f:
            kappa_script = f.read()
        fix_hot = re.search(r"fix\s+hot\s+all\s+([^\s]+)\s+", kappa_script)
        # placeholders supposed to be filled
        self.assertEqual(fix_hot.group(1), "heat")
        fix_cold = re.search(r"fix\s+cold\s+all\s+([^\s]+)\s+", kappa_script)
        self.assertEqual(fix_cold.group(1), "heat")
        lattice = re.search(r"lattice\s+fcc\s+(.*)\n", kappa_script)
        # parentheses not supposed to be filled
        self.assertEqual(lattice.group(1), "${rho}")
        pair_style = re.search(r"pair_style\slj/cut\s+(.*)\n", kappa_script)
        self.assertEqual(pair_style.group(1), "${rc}")

        with open(os.path.join(test_dir, "in.peptide")) as f:
            peptide_script = f.read()
        # copy data file
        src = os.path.join(test_dir, "data.quartz")
        write_lammps_inputs(output_dir="path", script_template=peptide_script,
                            data=src)
        dst = os.path.join("path", "data.peptide")
        self.assertTrue(filecmp.cmp(src, dst, shallow=False))
        # write data file from obj
        obj = LammpsData.from_file(src, atom_style="atomic")
        write_lammps_inputs(output_dir="obj", script_template=peptide_script,
                            data=obj)
        obj_read = LammpsData.from_file(os.path.join("obj", "data.peptide"),
                                        atom_style="atomic")
        pd.testing.assert_frame_equal(obj_read.masses, obj.masses)
        pd.testing.assert_frame_equal(obj_read.atoms, obj.atoms)

    @classmethod
    def tearDownClass(cls):
        temp_dirs = ["heat", "path", "obj"]
        for td in temp_dirs:
            if os.path.exists(td):
                shutil.rmtree(td)


if __name__ == "__main__":
    unittest.main()
