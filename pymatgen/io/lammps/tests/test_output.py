# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
import unittest

import numpy as np

from pymatgen.io.lammps.output import LammpsLog, LammpsDump


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files", "lammps")


class LammpsLogTest(unittest.TestCase):

    def test_init(self):
        comb_file = "log.5Oct16.comb.Si.elastic.g++.1"
        comb = LammpsLog(filename=os.path.join(test_dir, comb_file))
        self.assertEqual(len(comb.runs), 6)
        # first comb run
        comb0 = comb.runs[0]
        np.testing.assert_array_equal(["Step", "Temp", "TotEng", "PotEng",
                                       "E_vdwl", "E_coul"], comb0.columns)
        self.assertEqual(len(comb0), 6)
        comb0_data = [[0, 1, -4.6295947, -4.6297237, -4.6297237, 0],
                      [5, 1, -4.6295965, -4.6297255, -4.6297255, 0]]
        np.testing.assert_array_almost_equal(comb0.iloc[[0, -1]], comb0_data)
        # final comb run
        comb_1 = comb.runs[-1]
        np.testing.assert_array_equal(["Step", "Lx", "Ly", "Lz",
                                       "Xy", "Xz", "Yz",
                                       "c_fxy[1]", "c_fxy[2]", "c_fxy[3]",
                                       "c_fxy[4]", "c_fxy[5]", "c_fxy[6]"],
                                      comb_1.columns)
        self.assertEqual(len(comb_1), 11)
        comb_1_data = [[36, 5.1293854e-06], [46, 2192.8256]]
        np.testing.assert_array_almost_equal(comb_1.iloc[[0, -1], [0, -3]],
                                             comb_1_data)

        ehex_file = "log.13Oct16.ehex.g++.8"
        ehex = LammpsLog(filename=os.path.join(test_dir, ehex_file))
        self.assertEqual(len(ehex.runs), 3)
        ehex0, ehex1, ehex2 = ehex.runs
        # ehex run #1
        np.testing.assert_array_equal(["Step", "Temp", "E_pair", "E_mol",
                                       "TotEng", "Press"], ehex0.columns)
        self.assertEqual(len(ehex0), 11)
        ehex0_data = [[0, 1.35, -4.1241917, 0, -2.0994448, -3.1961612],
                      [1000, 1.3732017, -3.7100044, 0,
                       -1.6504594, 0.83982701]]
        np.testing.assert_array_almost_equal(ehex0.iloc[[0, -1]], ehex0_data)
        # ehex run #2
        np.testing.assert_array_equal(["Step", "Temp", "c_Thot", "c_Tcold"],
                                      ehex1.columns)
        self.assertEqual(len(ehex1), 11)
        ehex1_data = [[1000, 1.35, 1.431295, 1.2955644],
                      [11000, 1.3794051, 1.692299, 1.0515688]]
        np.testing.assert_array_almost_equal(ehex1.iloc[[0, -1]], ehex1_data)
        # ehex run #3
        np.testing.assert_array_equal(["Step", "Temp", "c_Thot", "c_Tcold",
                                       "v_tdiff", "f_ave"], ehex2.columns)
        self.assertEqual(len(ehex2), 21)
        ehex2_data = [[11000, 1.3794051, 1.6903393, 1.0515688, 0, 0],
                      [31000, 1.3822489, 1.8220413, 1.0322271, -0.7550338,
                       -0.76999077]]
        np.testing.assert_array_almost_equal(ehex2.iloc[[0, -1]], ehex2_data)

        peptide_file = "log.5Oct16.peptide.g++.1"
        peptide = LammpsLog(filename=os.path.join(test_dir, peptide_file))
        peptide0 = peptide.runs[0]
        np.testing.assert_array_equal(["Step", "TotEng", "KinEng", "Temp",
                                       "PotEng", "E_bond", "E_angle",
                                       "E_dihed", "E_impro", "E_vdwl",
                                       "E_coul", "E_long", "Press"],
                                      peptide0.columns)
        self.assertEqual(len(peptide0), 7)
        peptide0_select = peptide0.loc[[0, 6], ["Step", "TotEng", "Press"]]
        peptide0_data = [[0, -5237.4580, -837.0112],
                         [300, -5251.3637, -471.5505]]
        np.testing.assert_array_almost_equal(peptide0_select, peptide0_data)


class LammpsDumpTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.rdx_10 = LammpsDump(filename=os.path.join(test_dir,
                                                      "dump.rdx.gz"))
        cls.rdx_25 = LammpsDump(filename=os.path.join(test_dir,
                                                      "dump.rdx_wc.*"),
                                parse_box=False)
        cls.tatb = LammpsDump(filename=os.path.join(test_dir, "dump.tatb"))

    def test_init(self):
        files = ["dump.rdx_wc.0", "dump.rdx_wc.25", "dump.rdx_wc.50",
                 "dump.rdx_wc.75", "dump.rdx_wc.100"]
        self.assertListEqual([os.path.join(test_dir, f) for f in files],
                             self.rdx_25.all_files)

    def test_read(self):
        # general tests + gzipped
        rdx_10_data = list(self.rdx_10.read())
        timesteps_10 = [d["timestep"] for d in rdx_10_data]
        np.testing.assert_array_equal(timesteps_10, np.arange(0, 101, 10))
        obox = rdx_10_data[0]["box"]
        np.testing.assert_array_equal(obox.bounds, np.array([(35, 48)] * 3))
        atom = rdx_10_data[-1]["data"][-1]
        np.testing.assert_array_equal(atom,
                                      [19, 2, 0.42369, 0.47347, 0.555425])
        # timestep wildcard
        rdx_25_data = list(self.rdx_25.read())
        timesteps_25 = [d["timestep"] for d in rdx_25_data]
        np.testing.assert_array_equal(timesteps_25, np.arange(0, 101, 25))
        self.assertNotIn("box", rdx_25_data[0])
        # tilted box
        tbox = list(self.tatb.read())[0]["box"]
        bounds = [[0, 13.624], [0, 17.1149153805], [0, 15.1826391451]]
        tilt = [-5.75315630927, -6.325466, 7.4257288]
        np.testing.assert_array_almost_equal(tbox.bounds, bounds)
        np.testing.assert_array_almost_equal(tbox.tilt, tilt)


if __name__ == "__main__":
    unittest.main()
