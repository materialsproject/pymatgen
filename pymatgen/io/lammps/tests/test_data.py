# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import unittest
import os
import random

import numpy as np

from pymatgen.io.lammps.data import LammpsData


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files", "lammps")


class LammpsDataTest(unittest.TestCase):

    def test_from_file(self):
        peptide = LammpsData.from_file(filename=os.path.join(test_dir,
                                                             "data.peptide"),
                                       atom_style="full")

        # header stats
        self.assertEqual(len(peptide.atoms), 2004)
        topology = peptide.topology
        self.assertEqual(len(topology["Bonds"]), 1365)
        self.assertEqual(len(topology["Angles"]), 786)
        self.assertEqual(len(topology["Dihedrals"]), 207)
        self.assertEqual(len(topology["Impropers"]), 12)
        self.assertEqual(len(peptide.masses), 14)
        ff_coeffs = peptide.ff_coeffs
        self.assertEqual(len(ff_coeffs["Pair Coeffs"]), 14)
        self.assertEqual(len(ff_coeffs["Bond Coeffs"]), 18)
        self.assertEqual(len(ff_coeffs["Angle Coeffs"]), 31)
        self.assertEqual(len(ff_coeffs["Dihedral Coeffs"]), 21)
        self.assertEqual(len(ff_coeffs["Improper Coeffs"]), 2)
        # header box
        np.testing.assert_array_equal(peptide.box_bounds,
                                      [[36.840194, 64.211560],
                                       [41.013691, 68.385058],
                                       [29.768095, 57.139462]])
        # body  last line of each section
        self.assertDictEqual(peptide.masses[-1], {"id": 14, "mass": 1.0100})
        self.assertDictEqual(ff_coeffs["Pair Coeffs"][-1],
                             {"id": 14, "coeffs": [0.046000, 0.400014,
                                                   0.046000, 0.400014]})
        self.assertDictEqual(ff_coeffs["Bond Coeffs"][-1],
                             {"id": 18, "coeffs": [450.000000, 0.957200]})
        self.assertDictEqual(ff_coeffs["Angle Coeffs"][-1],
                             {"id": 31, "coeffs": [55.000000, 104.520000,
                                                   0.000000, 0.000000]})
        self.assertDictEqual(ff_coeffs["Dihedral Coeffs"][-1],
                             {"id": 21, "coeffs": [0.010000, 3, 0, 1.000000]})
        for c in ff_coeffs["Dihedral Coeffs"][-1]["coeffs"][1:2]:
            self.assertIsInstance(c, int)
        self.assertDictEqual(ff_coeffs["Improper Coeffs"][-1],
                             {"id": 2, "coeffs": [20.000000, 0.000000]})
        self.assertDictEqual(peptide.atoms[-1],
                             {"id": 2004, "molecule-ID": 641, "type": 14,
                              "q": 0.417, "x": 56.55074, "y": 49.75049,
                              "z": 48.61854, "nx": 1, "ny": 1, "nz": 1})
        self.assertDictEqual(peptide.velocities[-1],
                             {"id": 2004, "velocity": [-0.010076,
                                                       -0.005729,
                                                       -0.026032]})
        self.assertDictEqual(topology["Bonds"][-1],
                             {"id": 1365, "type": 18, "bond": [2002, 2003]})
        self.assertDictEqual(topology["Angles"][-1],
                             {"id": 786, "type": 31,
                              "angle": [2003, 2002, 2004]})
        self.assertDictEqual(topology["Dihedrals"][-1],
                             {"id": 207, "type": 3,
                              "dihedral": [64, 79, 80, 82]})
        self.assertDictEqual(topology["Impropers"][-1],
                             {"id": 12, "type": 2,
                              "improper": [79, 64, 80, 81]})
        # box_tilt and another atom_style
        quartz = LammpsData.from_file(filename=os.path.join(test_dir,
                                                            "data.quartz"),
                                      atom_style="atomic")
        np.testing.assert_array_equal(quartz.box_tilt, [-2.456700, 0.0, 0.0])
        self.assertDictEqual(quartz.atoms[-1],
                             {"id": 9, "type": 2, "x": 1.375998,
                              "y": -1.140800, "z": -2.443511})
        # sort_id
        nvt = LammpsData.from_file(filename=os.path.join(test_dir,
                                                         "nvt.data"),
                                   sort_id=True)
        atom_id = random.randint(1, 648)
        self.assertEqual(nvt.atoms[atom_id - 1]["id"], atom_id)
        # PairIJ Coeffs section
        virus = LammpsData.from_file(filename=os.path.join(test_dir,
                                                           "virus.data"),
                                     atom_style="angle")
        n = len(virus.masses)
        pairij = virus.ff_coeffs["PairIJ Coeffs"]
        self.assertEqual(len(pairij), n * (n + 1) / 2)
        self.assertDictEqual(pairij[-1],
                             {"id1": 4, "id2": 4, "coeffs": [1, 1, 1.1225]})


if __name__ == "__main__":
    unittest.main()
