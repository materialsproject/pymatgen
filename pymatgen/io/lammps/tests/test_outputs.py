# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
import unittest
import os
import json

import numpy as np
import pandas as pd

from pymatgen.io.lammps.outputs import LammpsDump, parse_lammps_dumps,\
    parse_lammps_log


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files", "lammps")


class LammpsDumpTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(os.path.join(test_dir, "dump.rdx_wc.100")) as f:
            rdx_str = f.read()
        cls.rdx = LammpsDump.from_string(string=rdx_str)
        with open(os.path.join(test_dir, "dump.tatb")) as f:
            tatb_str = f.read()
        cls.tatb = LammpsDump.from_string(string=tatb_str)

    def test_from_string(self):
        self.assertEqual(self.rdx.timestep, 100)
        self.assertEqual(self.rdx.natoms, 21)
        np.testing.assert_array_equal(self.rdx.box.bounds,
                                      np.array([(35, 48)] * 3))
        np.testing.assert_array_equal(self.rdx.data.columns,
                                      ["id", "type", "xs", "ys", "zs"])
        rdx_data = self.rdx.data.iloc[-1]
        rdx_data_target = [19, 2, 0.42369, 0.47347, 0.555425]
        np.testing.assert_array_almost_equal(rdx_data, rdx_data_target)

        self.assertEqual(self.tatb.timestep, 0)
        self.assertEqual(self.tatb.natoms, 384)
        bounds = [[0, 13.624], [0, 17.1149153805], [0, 15.1826391451]]
        np.testing.assert_array_almost_equal(self.tatb.box.bounds, bounds)
        tilt = [-5.75315630927, -6.325466, 7.4257288]
        np.testing.assert_array_almost_equal(self.tatb.box.tilt, tilt)
        np.testing.assert_array_equal(self.tatb.data.columns,
                                      ["id", "type", "q", "x", "y", "z"])
        tatb_data = self.tatb.data.iloc[-1]
        tatb_data_target = [356, 3, -0.482096, 2.58647, 12.9577, 14.3143]
        np.testing.assert_array_almost_equal(tatb_data, tatb_data_target)

    def test_json_dict(self):
        encoded = json.dumps(self.rdx.as_dict())
        decoded = json.loads(encoded)
        rdx = LammpsDump.from_dict(decoded)
        self.assertEqual(rdx.timestep, 100)
        self.assertEqual(rdx.natoms, 21)
        np.testing.assert_array_equal(rdx.box.bounds,
                                      np.array([(35, 48)] * 3))
        pd.testing.assert_frame_equal(rdx.data, self.rdx.data)


class FuncTest(unittest.TestCase):

    def test_parse_lammps_dumps(self):
        # gzipped
        rdx_10_pattern = os.path.join(test_dir, "dump.rdx.gz")
        rdx_10 = list(parse_lammps_dumps(file_pattern=rdx_10_pattern))
        timesteps_10 = [d.timestep for d in rdx_10]
        np.testing.assert_array_equal(timesteps_10, np.arange(0, 101, 10))
        self.assertTupleEqual(rdx_10[-1].data.shape, (21, 5))
        # wildcard
        rdx_25_pattern = os.path.join(test_dir, "dump.rdx_wc.*")
        rdx_25 = list(parse_lammps_dumps(file_pattern=rdx_25_pattern))
        timesteps_25 = [d.timestep for d in rdx_25]
        np.testing.assert_array_equal(timesteps_25, np.arange(0, 101, 25))
        self.assertTupleEqual(rdx_25[-1].data.shape, (21, 5))

    def test_parse_lammps_log(self):
        comb_file = "log.5Oct16.comb.Si.elastic.g++.1"
        comb = parse_lammps_log(filename=os.path.join(test_dir, comb_file))
        self.assertEqual(len(comb), 6)
        # first comb run
        comb0 = comb[0]
        np.testing.assert_array_equal(["Step", "Temp", "TotEng", "PotEng",
                                       "E_vdwl", "E_coul"], comb0.columns)
        self.assertEqual(len(comb0), 6)
        comb0_data = [[0, 1, -4.6295947, -4.6297237, -4.6297237, 0],
                      [5, 1, -4.6295965, -4.6297255, -4.6297255, 0]]
        np.testing.assert_array_almost_equal(comb0.iloc[[0, -1]], comb0_data)
        # final comb run
        comb_1 = comb[-1]
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
        ehex = parse_lammps_log(filename=os.path.join(test_dir, ehex_file))
        self.assertEqual(len(ehex), 3)
        ehex0, ehex1, ehex2 = ehex
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
        peptide = parse_lammps_log(filename=os.path.join(test_dir,
                                                         peptide_file))
        peptide0 = peptide[0]
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


if __name__ == "__main__":
    unittest.main()
