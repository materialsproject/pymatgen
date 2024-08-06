from __future__ import annotations

import os

import numpy as np
from numpy.testing import assert_array_equal

from pymatgen.core.structure import Structure
from pymatgen.core.units import Ry_to_eV
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.lmto import LMTOCopl, LMTOCtrl
from pymatgen.util.num import round_to_sigfigs
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

__author__ = "Marco Esters"
__copyright__ = "Copyright 2017, The Materials Project"
__version__ = "0.1"
__email__ = "esters@uoregon.edu"
__date__ = "Nov 30, 2017"


TEST_DIR = f"{TEST_FILES_DIR}/electronic_structure/cohp"
module_dir = os.path.dirname(os.path.abspath(__file__))


class TestCtrl(PymatgenTest):
    def setUp(self):
        os.chdir(TEST_DIR)
        self.ref_bise = LMTOCtrl.from_file(filename="CTRL.BiSe")
        self.ref_fe = LMTOCtrl.from_file()

    def tearDown(self):
        os.chdir(module_dir)

    def test_dict(self):
        assert self.ref_bise == LMTOCtrl.from_dict(self.ref_bise.as_dict())

    def test_structure(self):
        bise_poscar = Structure.from_file("POSCAR.BiSe")
        assert bise_poscar.matches(self.ref_bise.structure)
        assert self.ref_bise == LMTOCtrl(self.ref_bise.structure, header="Bi6Se6, hexagonal")

    def test_read_write(self):
        ctrl_path = f"{self.tmp_path}/CTRL.tmp"
        self.ref_bise.write_file(filename=ctrl_path)
        ctrl_file = LMTOCtrl.from_file(filename=ctrl_path)
        assert self.ref_bise.structure.matches(ctrl_file.structure)


class TestCopl(PymatgenTest):
    def setUp(self):
        os.chdir(TEST_DIR)
        self.copl_bise = LMTOCopl("COPL.BiSe")
        self.copl_bise_eV = LMTOCopl(filename="COPL.BiSe", to_eV=True)
        self.copl_fe = LMTOCopl()

    def tearDown(self):
        os.chdir(module_dir)

    def test_attributes(self):
        assert not self.copl_bise.is_spin_polarized
        assert self.copl_fe.is_spin_polarized
        assert len(self.copl_bise.energies) == 801
        assert len(self.copl_fe.energies) == 801
        assert len(self.copl_bise.cohp_data) == 7
        assert len(self.copl_fe.cohp_data) == 8

    def test_cohp_data(self):
        lengths_sites_bise = {
            "Bi1-Se7": (2.882, (0, 6)),
            "Bi1-Se9": (3.102, (0, 8)),
            "Bi3-Se11": (2.917, (2, 10)),
            "Bi3-Se9": (3.050, (2, 8)),
            "Bi5-Bi6": (3.073, (4, 5)),
            "Bi5-Se11": (3.375, (4, 10)),
            "Se7-Se8": (3.364, (6, 7)),
        }
        for bond in self.copl_bise.cohp_data:
            assert self.copl_bise.cohp_data[bond]["length"] == lengths_sites_bise[bond][0]
            assert self.copl_bise.cohp_data[bond]["sites"] == lengths_sites_bise[bond][1]
        labels_fe = ["Fe1-Fe1"] + [f"Fe1-Fe1-{i}" for i in range(1, 8)]
        assert sorted(self.copl_fe.cohp_data) == labels_fe
        for bond in labels_fe:
            assert self.copl_fe.cohp_data[bond]["length"] == 2.482
            assert self.copl_fe.cohp_data[bond]["sites"] == (0, 0)

    def test_energies(self):
        assert self.copl_bise.efermi == -0.17223
        assert self.copl_bise_eV.efermi == -2.3433
        assert self.copl_fe.efermi == -0.085683
        ener_eV = np.array(
            [round_to_sigfigs(energy, 5) for energy in self.copl_bise.energies * Ry_to_eV],
            dtype=float,
        )
        assert_array_equal(ener_eV, self.copl_bise_eV.energies)
        copl_icohp = self.copl_bise.cohp_data["Bi1-Se7"]["ICOHP"][Spin.up]
        icohp = np.array([round_to_sigfigs(i, 5) for i in copl_icohp * Ry_to_eV], dtype=float)
        icohp_eV = self.copl_bise_eV.cohp_data["Bi1-Se7"]["ICOHP"][Spin.up]
        assert_array_equal(icohp, icohp_eV)
