# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import unittest
import numpy as np
import os
from pymatgen.core.structure import Structure
from pymatgen.core.units import Ry_to_eV
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.lmto import LMTOCtrl, LMTOCopl
from pymatgen.util.num import round_to_sigfigs
from pymatgen.util.testing import PymatgenTest

__author__ = "Marco Esters"
__copyright__ = "Copyright 2017, The Materials Project"
__version__ = "0.1"
__email__ = "esters@uoregon.edu"
__date__ = "Nov 30, 2017"

test_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "..", "..", "test_files", "cohp")
this_dir = os.path.dirname(os.path.abspath(__file__))


class CtrlTest(unittest.TestCase):
    def setUp(self):
        os.chdir(test_dir)
        self.ctrl_bise = LMTOCtrl.from_file(filename="CTRL.BiSe")
        self.ctrl_fe = LMTOCtrl.from_file()

    def tearDown(self):
        os.chdir(this_dir)

    def test_dict(self):
        self.assertEqual(self.ctrl_bise,
                         LMTOCtrl.from_dict(self.ctrl_bise.as_dict()))

    def test_structure(self):
        bise_poscar = Structure.from_file("POSCAR.BiSe")
        self.assertTrue(bise_poscar.matches(self.ctrl_bise.structure))
        self.assertEqual(self.ctrl_bise,
                         LMTOCtrl(self.ctrl_bise.structure,
                                  header="Bi6Se6, hexagonal"))

    def test_read_write(self):
        self.ctrl_bise.write_file(filename="CTRL.tmp")
        ctrl_tmp = LMTOCtrl.from_file(filename="CTRL.tmp")
        self.assertTrue(self.ctrl_bise.structure.matches(ctrl_tmp.structure))
        os.remove("CTRL.tmp")


class CoplTest(PymatgenTest):
    def setUp(self):
        os.chdir(test_dir)
        self.copl_bise = LMTOCopl("COPL.BiSe")
        self.copl_bise_eV = LMTOCopl(filename="COPL.BiSe", to_eV=True)
        self.copl_fe = LMTOCopl()

    def tearDown(self):
        os.chdir(this_dir)

    def test_attributes(self):
        self.assertFalse(self.copl_bise.is_spin_polarized)
        self.assertTrue(self.copl_fe.is_spin_polarized)
        self.assertEqual(len(self.copl_bise.energies), 801)
        self.assertEqual(len(self.copl_fe.energies), 801)
        self.assertEqual(len(self.copl_bise.cohp_data), 7)
        self.assertEqual(len(self.copl_fe.cohp_data), 8)

    def test_cohp_data(self):
        lengths_sites_bise = {"Bi1-Se7": (2.882, (0, 6)),
                              "Bi1-Se9": (3.102, (0, 8)),
                              "Bi3-Se11": (2.917, (2, 10)),
                              "Bi3-Se9": (3.050, (2, 8)),
                              "Bi5-Bi6": (3.073, (4, 5)),
                              "Bi5-Se11": (3.375, (4, 10)),
                              "Se7-Se8": (3.364, (6, 7))}
        for bond in self.copl_bise.cohp_data:
            self.assertEqual(self.copl_bise.cohp_data[bond]["length"],
                             lengths_sites_bise[bond][0])
            self.assertEqual(self.copl_bise.cohp_data[bond]["sites"],
                             lengths_sites_bise[bond][1])
        labels_fe = ["Fe1-Fe1"] + ["Fe1-Fe1-%d" % i for i in range(1, 8)]
        self.assertEqual(sorted(self.copl_fe.cohp_data.keys()), labels_fe)
        for bond in labels_fe:
            self.assertEqual(self.copl_fe.cohp_data[bond]["length"], 2.482)
            self.assertEqual(self.copl_fe.cohp_data[bond]["sites"], (0, 0))

    def test_energies(self):
        self.assertEqual(self.copl_bise.efermi, -0.17223)
        self.assertEqual(self.copl_bise_eV.efermi, -2.3433)
        self.assertEqual(self.copl_fe.efermi, -0.085683)
        ener_eV = np.array([round_to_sigfigs(energy, 5)
                            for energy in self.copl_bise.energies * Ry_to_eV],
                           dtype=float)
        self.assertArrayEqual(ener_eV, self.copl_bise_eV.energies)
        copl_icohp = self.copl_bise.cohp_data["Bi1-Se7"]["ICOHP"][Spin.up]
        icohp = np.array([round_to_sigfigs(i, 5)
                          for i in copl_icohp * Ry_to_eV],
                         dtype=float)
        icohp_eV = self.copl_bise_eV.cohp_data["Bi1-Se7"]["ICOHP"][Spin.up]
        self.assertArrayEqual(icohp, icohp_eV)

if __name__ == "__main__":
    unittest.main()
