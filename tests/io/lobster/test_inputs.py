from __future__ import annotations

import json
import os
import tempfile
import unittest

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from pytest import approx

from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Orbital, Spin
from pymatgen.io.lobster import (
    Bandoverlaps,
    Charge,
    Cohpcar,
    Doscar,
    Fatband,
    Grosspop,
    Icohplist,
    Lobsterin,
    LobsterMatrices,
    Lobsterout,
    MadelungEnergies,
    NciCobiList,
    SitePotential,
    Wavefunction,
)
from pymatgen.io.lobster.inputs import get_all_possible_basis_combinations
from pymatgen.io.vasp import Vasprun
from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

__author__ = "Janine George, Marco Esters"
__copyright__ = "Copyright 2017, The Materials Project"
__version__ = "0.2"
__email__ = "janine.george@uclouvain.be, esters@uoregon.edu"
__date__ = "Dec 10, 2017"


module_dir = os.path.dirname(os.path.abspath(__file__))


class TestCohpcar(PymatgenTest):
    def setUp(self):
        self.cohp_bise = Cohpcar(filename=f"{TEST_FILES_DIR}/cohp/COHPCAR.lobster.BiSe.gz")
        self.coop_bise = Cohpcar(
            filename=f"{TEST_FILES_DIR}/cohp/COOPCAR.lobster.BiSe.gz",
            are_coops=True,
        )
        self.cohp_fe = Cohpcar(filename=f"{TEST_FILES_DIR}/cohp/COOPCAR.lobster.gz")
        self.coop_fe = Cohpcar(
            filename=f"{TEST_FILES_DIR}/cohp/COOPCAR.lobster.gz",
            are_coops=True,
        )
        self.orb = Cohpcar(filename=f"{TEST_FILES_DIR}/cohp/COHPCAR.lobster.orbitalwise.gz")
        self.orb_notot = Cohpcar(filename=f"{TEST_FILES_DIR}/cohp/COHPCAR.lobster.notot.orbitalwise.gz")

        # Lobster 3.1 (Test data is from prerelease of Lobster 3.1)
        self.cohp_KF = Cohpcar(filename=f"{TEST_FILES_DIR}/cohp/COHPCAR.lobster.KF.gz")
        self.coop_KF = Cohpcar(
            filename=f"{TEST_FILES_DIR}/cohp/COHPCAR.lobster.KF.gz",
            are_coops=True,
        )

        # example with f electrons
        self.cohp_Na2UO4 = Cohpcar(filename=f"{TEST_FILES_DIR}/cohp/COHPCAR.lobster.Na2UO4.gz")
        self.coop_Na2UO4 = Cohpcar(
            filename=f"{TEST_FILES_DIR}/cohp/COOPCAR.lobster.Na2UO4.gz",
            are_coops=True,
        )
        self.cobi = Cohpcar(
            filename=f"{TEST_FILES_DIR}/cohp/COBICAR.lobster.gz",
            are_cobis=True,
        )

    def test_attributes(self):
        assert not self.cohp_bise.are_coops
        assert self.coop_bise.are_coops
        assert not self.cohp_bise.is_spin_polarized
        assert not self.coop_bise.is_spin_polarized
        assert not self.cohp_fe.are_coops
        assert self.coop_fe.are_coops
        assert self.cohp_fe.is_spin_polarized
        assert self.coop_fe.is_spin_polarized
        assert len(self.cohp_bise.energies) == 241
        assert len(self.coop_bise.energies) == 241
        assert len(self.cohp_fe.energies) == 301
        assert len(self.coop_fe.energies) == 301
        assert len(self.cohp_bise.cohp_data) == 12
        assert len(self.coop_bise.cohp_data) == 12
        assert len(self.cohp_fe.cohp_data) == 3
        assert len(self.coop_fe.cohp_data) == 3

        # Lobster 3.1
        assert not self.cohp_KF.are_coops
        assert self.coop_KF.are_coops
        assert not self.cohp_KF.is_spin_polarized
        assert not self.coop_KF.is_spin_polarized
        assert len(self.cohp_KF.energies) == 6
        assert len(self.coop_KF.energies) == 6
        assert len(self.cohp_KF.cohp_data) == 7
        assert len(self.coop_KF.cohp_data) == 7

        # Lobster 4.1.0
        assert not self.cohp_KF.are_cobis
        assert not self.coop_KF.are_cobis
        assert not self.cobi.are_coops
        assert self.cobi.are_cobis
        assert not self.cobi.is_spin_polarized

    def test_energies(self):
        efermi_bise = 5.90043
        elim_bise = (-0.124679, 11.9255)
        efermi_fe = 9.75576
        elim_fe = (-0.277681, 14.7725)
        efermi_KF = -2.87475
        elim_KF = (-11.25000 + efermi_KF, 7.5000 + efermi_KF)

        assert self.cohp_bise.efermi == efermi_bise
        assert self.coop_bise.efermi == efermi_bise
        assert self.cohp_fe.efermi == efermi_fe
        assert self.coop_fe.efermi == efermi_fe
        # Lobster 3.1
        assert self.cohp_KF.efermi == efermi_KF
        assert self.coop_KF.efermi == efermi_KF

        assert self.cohp_bise.energies[0] + self.cohp_bise.efermi == approx(elim_bise[0], abs=1e-4)
        assert self.cohp_bise.energies[-1] + self.cohp_bise.efermi == approx(elim_bise[1], abs=1e-4)
        assert self.coop_bise.energies[0] + self.coop_bise.efermi == approx(elim_bise[0], abs=1e-4)
        assert self.coop_bise.energies[-1] + self.coop_bise.efermi == approx(elim_bise[1], abs=1e-4)

        assert self.cohp_fe.energies[0] + self.cohp_fe.efermi == approx(elim_fe[0], abs=1e-4)
        assert self.cohp_fe.energies[-1] + self.cohp_fe.efermi == approx(elim_fe[1], abs=1e-4)
        assert self.coop_fe.energies[0] + self.coop_fe.efermi == approx(elim_fe[0], abs=1e-4)
        assert self.coop_fe.energies[-1] + self.coop_fe.efermi == approx(elim_fe[1], abs=1e-4)

        # Lobster 3.1
        assert self.cohp_KF.energies[0] + self.cohp_KF.efermi == approx(elim_KF[0], abs=1e-4)
        assert self.cohp_KF.energies[-1] + self.cohp_KF.efermi == approx(elim_KF[1], abs=1e-4)
        assert self.coop_KF.energies[0] + self.coop_KF.efermi == approx(elim_KF[0], abs=1e-4)
        assert self.coop_KF.energies[-1] + self.coop_KF.efermi == approx(elim_KF[1], abs=1e-4)

    def test_cohp_data(self):
        lengths_sites_bise = {
            "1": (2.882308829886294, (0, 6)),
            "2": (3.1014396233274444, (0, 9)),
            "3": (2.8823088298862083, (1, 7)),
            "4": (3.1014396233275434, (1, 8)),
            "5": (3.0500070394403904, (2, 9)),
            "6": (2.9167594580335807, (2, 10)),
            "7": (3.05000703944039, (3, 8)),
            "8": (2.9167594580335803, (3, 11)),
            "9": (3.3752173204052101, (4, 11)),
            "10": (3.0729354518345948, (4, 5)),
            "11": (3.3752173204052101, (5, 10)),
        }
        lengths_sites_fe = {
            "1": (2.8318907764979082, (7, 6)),
            "2": (2.4524893531900283, (7, 8)),
        }
        # Lobster 3.1
        lengths_sites_KF = {
            "1": (2.7119923200622269, (0, 1)),
            "2": (2.7119923200622269, (0, 1)),
            "3": (2.7119923576010501, (0, 1)),
            "4": (2.7119923576010501, (0, 1)),
            "5": (2.7119923200622269, (0, 1)),
            "6": (2.7119923200622269, (0, 1)),
        }

        for data in [self.cohp_bise.cohp_data, self.coop_bise.cohp_data]:
            for bond, val in data.items():
                if bond != "average":
                    assert val["length"] == lengths_sites_bise[bond][0]
                    assert val["sites"] == lengths_sites_bise[bond][1]
                    assert len(val["COHP"][Spin.up]) == 241
                    assert len(val["ICOHP"][Spin.up]) == 241
        for data in [self.cohp_fe.cohp_data, self.coop_fe.cohp_data]:
            for bond, val in data.items():
                if bond != "average":
                    assert val["length"] == lengths_sites_fe[bond][0]
                    assert val["sites"] == lengths_sites_fe[bond][1]
                    assert len(val["COHP"][Spin.up]) == 301
                    assert len(val["ICOHP"][Spin.up]) == 301

        # Lobster 3.1
        for data in [self.cohp_KF.cohp_data, self.coop_KF.cohp_data]:
            for bond, val in data.items():
                if bond != "average":
                    assert val["length"] == lengths_sites_KF[bond][0]
                    assert val["sites"] == lengths_sites_KF[bond][1]
                    assert len(val["COHP"][Spin.up]) == 6
                    assert len(val["ICOHP"][Spin.up]) == 6

    def test_orbital_resolved_cohp(self):
        orbitals = [(Orbital(i), Orbital(j)) for j in range(4) for i in range(4)]
        assert self.cohp_bise.orb_res_cohp is None
        assert self.coop_bise.orb_res_cohp is None
        assert self.cohp_fe.orb_res_cohp is None
        assert self.coop_fe.orb_res_cohp is None
        assert self.orb_notot.cohp_data["1"]["COHP"] is None
        assert self.orb_notot.cohp_data["1"]["ICOHP"] is None
        for orbs in self.orb.orb_res_cohp["1"]:
            orb_set = self.orb.orb_res_cohp["1"][orbs]["orbitals"]
            assert orb_set[0][0] == 4
            assert orb_set[1][0] == 4
            assert (orb_set[0][1], orb_set[1][1]) in orbitals

        # test d and f orbitals
        ref_list1 = [*[5] * 28, *[6] * 36, *[7] * 4]
        ref_list2 = [
            *["f0"] * 4,
            *["f1"] * 4,
            *["f2"] * 4,
            *["f3"] * 4,
            *["f_1"] * 4,
            *["f_2"] * 4,
            *["f_3"] * 4,
            *["dx2"] * 4,
            *["dxy"] * 4,
            *["dxz"] * 4,
            *["dyz"] * 4,
            *["dz2"] * 4,
            *["px"] * 4,
            *["py"] * 4,
            *["pz"] * 4,
            *["s"] * 8,
        ]
        for iorb, orbs in enumerate(sorted(self.cohp_Na2UO4.orb_res_cohp["49"])):
            orb_set = self.cohp_Na2UO4.orb_res_cohp["49"][orbs]["orbitals"]
            assert orb_set[0][0] == ref_list1[iorb]
            assert str(orb_set[0][1]) == ref_list2[iorb]

        # The sum of the orbital-resolved COHPs should be approximately
        # the total COHP. Due to small deviations in the LOBSTER calculation,
        # the precision is not very high though.
        cohp = self.orb.cohp_data["1"]["COHP"][Spin.up]
        icohp = self.orb.cohp_data["1"]["ICOHP"][Spin.up]
        tot = np.sum(
            [self.orb.orb_res_cohp["1"][orbs]["COHP"][Spin.up] for orbs in self.orb.orb_res_cohp["1"]],
            axis=0,
        )
        assert_allclose(tot, cohp, atol=1e-3)
        tot = np.sum(
            [self.orb.orb_res_cohp["1"][orbs]["ICOHP"][Spin.up] for orbs in self.orb.orb_res_cohp["1"]],
            axis=0,
        )
        assert_allclose(tot, icohp, atol=1e-3)

        # Lobster 3.1
        cohp_KF = self.cohp_KF.cohp_data["1"]["COHP"][Spin.up]
        icohp_KF = self.cohp_KF.cohp_data["1"]["ICOHP"][Spin.up]
        tot_KF = np.sum(
            [self.cohp_KF.orb_res_cohp["1"][orbs]["COHP"][Spin.up] for orbs in self.cohp_KF.orb_res_cohp["1"]],
            axis=0,
        )
        assert_allclose(tot_KF, cohp_KF, atol=1e-3)
        tot_KF = np.sum(
            [self.cohp_KF.orb_res_cohp["1"][orbs]["ICOHP"][Spin.up] for orbs in self.cohp_KF.orb_res_cohp["1"]],
            axis=0,
        )
        assert_allclose(tot_KF, icohp_KF, atol=1e-3)

        # d and f orbitals
        cohp_Na2UO4 = self.cohp_Na2UO4.cohp_data["49"]["COHP"][Spin.up]
        icohp_Na2UO4 = self.cohp_Na2UO4.cohp_data["49"]["ICOHP"][Spin.up]
        tot_Na2UO4 = np.sum(
            [
                self.cohp_Na2UO4.orb_res_cohp["49"][orbs]["COHP"][Spin.up]
                for orbs in self.cohp_Na2UO4.orb_res_cohp["49"]
            ],
            axis=0,
        )
        assert_allclose(tot_Na2UO4, cohp_Na2UO4, atol=1e-3)
        tot_Na2UO4 = np.sum(
            [
                self.cohp_Na2UO4.orb_res_cohp["49"][orbs]["ICOHP"][Spin.up]
                for orbs in self.cohp_Na2UO4.orb_res_cohp["49"]
            ],
            axis=0,
        )
        assert_allclose(tot_Na2UO4, icohp_Na2UO4, atol=1e-3)


class TestIcohplist(unittest.TestCase):
    def setUp(self):
        self.icohp_bise = Icohplist(filename=f"{TEST_FILES_DIR}/cohp/ICOHPLIST.lobster.BiSe")
        self.icoop_bise = Icohplist(
            filename=f"{TEST_FILES_DIR}/cohp/ICOOPLIST.lobster.BiSe",
            are_coops=True,
        )
        self.icohp_fe = Icohplist(filename=f"{TEST_FILES_DIR}/cohp/ICOHPLIST.lobster")
        # allow gzipped files
        self.icohp_gzipped = Icohplist(filename=f"{TEST_FILES_DIR}/cohp/ICOHPLIST.lobster.gz")
        self.icoop_fe = Icohplist(
            filename=f"{TEST_FILES_DIR}/cohp/ICOHPLIST.lobster",
            are_coops=True,
        )
        # ICOBIs and orbitalwise ICOBILIST.lobster
        self.icobi_orbitalwise = Icohplist(
            filename=f"{TEST_FILES_DIR}/cohp/ICOBILIST.lobster",
            are_cobis=True,
        )
        # TODO: test orbitalwise ICOHPs with and without spin polarization

        self.icobi = Icohplist(
            filename=f"{TEST_FILES_DIR}/cohp/ICOBILIST.lobster.withoutorbitals",
            are_cobis=True,
        )
        self.icobi_orbitalwise_spinpolarized = Icohplist(
            filename=f"{TEST_FILES_DIR}/cohp/ICOBILIST.lobster.spinpolarized",
            are_cobis=True,
        )
        # make sure the correct line is read to check if this is a orbitalwise ICOBILIST
        self.icobi_orbitalwise_add = Icohplist(
            filename=f"{TEST_FILES_DIR}/cohp/ICOBILIST.lobster.additional_case",
            are_cobis=True,
        )
        self.icobi_orbitalwise_spinpolarized_add = Icohplist(
            filename=os.path.join(
                TEST_FILES_DIR,
                "cohp",
                "ICOBILIST.lobster.spinpolarized.additional_case",
            ),
            are_cobis=True,
        )

    def test_attributes(self):
        assert not self.icohp_bise.are_coops
        assert self.icoop_bise.are_coops
        assert not self.icohp_bise.is_spin_polarized
        assert not self.icoop_bise.is_spin_polarized
        assert len(self.icohp_bise.icohplist) == 11
        assert len(self.icoop_bise.icohplist) == 11
        assert not self.icohp_fe.are_coops
        assert self.icoop_fe.are_coops
        assert self.icohp_fe.is_spin_polarized
        assert self.icoop_fe.is_spin_polarized
        assert len(self.icohp_fe.icohplist) == 2
        assert len(self.icoop_fe.icohplist) == 2
        # test are_cobis
        assert not self.icohp_fe.are_coops
        assert not self.icohp_fe.are_cobis
        assert self.icoop_fe.are_coops
        assert not self.icoop_fe.are_cobis
        assert self.icobi.are_cobis
        assert not self.icobi.are_coops

        # orbitalwise
        assert self.icobi_orbitalwise.orbitalwise
        assert not self.icobi.orbitalwise

        assert self.icobi_orbitalwise_spinpolarized.orbitalwise

        assert self.icobi_orbitalwise_add.orbitalwise
        assert self.icobi_orbitalwise_spinpolarized_add.orbitalwise

    def test_values(self):
        icohplist_bise = {
            "1": {
                "length": 2.88231,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -2.18042},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "2": {
                "length": 3.10144,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -1.14347},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "3": {
                "length": 2.88231,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -2.18042},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "4": {
                "length": 3.10144,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -1.14348},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "5": {
                "length": 3.05001,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -1.30006},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "6": {
                "length": 2.91676,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -1.96843},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "7": {
                "length": 3.05001,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -1.30006},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "8": {
                "length": 2.91676,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -1.96843},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "9": {
                "length": 3.37522,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.47531},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "10": {
                "length": 3.07294,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -2.38796},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "11": {
                "length": 3.37522,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.47531},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
        }
        icooplist_bise = {
            "1": {
                "length": 2.88231,
                "number_of_bonds": 3,
                "icohp": {Spin.up: 0.14245},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "2": {
                "length": 3.10144,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.04118},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "3": {
                "length": 2.88231,
                "number_of_bonds": 3,
                "icohp": {Spin.up: 0.14245},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "4": {
                "length": 3.10144,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.04118},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "5": {
                "length": 3.05001,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.03516},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "6": {
                "length": 2.91676,
                "number_of_bonds": 3,
                "icohp": {Spin.up: 0.10745},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "7": {
                "length": 3.05001,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.03516},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "8": {
                "length": 2.91676,
                "number_of_bonds": 3,
                "icohp": {Spin.up: 0.10745},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "9": {
                "length": 3.37522,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.12395},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "10": {
                "length": 3.07294,
                "number_of_bonds": 3,
                "icohp": {Spin.up: 0.24714},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "11": {
                "length": 3.37522,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.12395},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
        }
        icooplist_fe = {
            "1": {
                "length": 2.83189,
                "number_of_bonds": 2,
                "icohp": {Spin.up: -0.10218, Spin.down: -0.19701},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
            "2": {
                "length": 2.45249,
                "number_of_bonds": 1,
                "icohp": {Spin.up: -0.28485, Spin.down: -0.58279},
                "translation": [0, 0, 0],
                "orbitals": None,
            },
        }

        assert icohplist_bise == self.icohp_bise.icohplist
        assert self.icohp_bise.icohpcollection.extremum_icohpvalue() == -2.38796
        assert icooplist_fe == self.icoop_fe.icohplist
        assert self.icoop_fe.icohpcollection.extremum_icohpvalue() == -0.29919
        assert icooplist_bise == self.icoop_bise.icohplist
        assert self.icoop_bise.icohpcollection.extremum_icohpvalue() == 0.24714
        assert self.icobi.icohplist["1"]["icohp"][Spin.up] == approx(0.58649)
        assert self.icobi_orbitalwise.icohplist["2"]["icohp"][Spin.up] == approx(0.58649)
        assert self.icobi_orbitalwise.icohplist["1"]["icohp"][Spin.up] == approx(0.58649)
        assert self.icobi_orbitalwise_spinpolarized.icohplist["1"]["icohp"][Spin.up] == approx(0.58649 / 2, abs=1e-3)
        assert self.icobi_orbitalwise_spinpolarized.icohplist["1"]["icohp"][Spin.down] == approx(0.58649 / 2, abs=1e-3)
        assert self.icobi_orbitalwise_spinpolarized.icohplist["2"]["icohp"][Spin.down] == approx(0.58649 / 2, abs=1e-3)
        assert self.icobi.icohpcollection.extremum_icohpvalue() == 0.58649
        assert self.icobi_orbitalwise_spinpolarized.icohplist["2"]["orbitals"]["2s-6s"]["icohp"][Spin.up] == 0.0247


class TestNciCobiList(unittest.TestCase):
    def setUp(self):
        self.ncicobi = NciCobiList(filename=f"{TEST_FILES_DIR}/cohp/NcICOBILIST.lobster")
        self.ncicobi_gz = NciCobiList(filename=f"{TEST_FILES_DIR}/cohp/NcICOBILIST.lobster.gz")
        self.ncicobi_no_spin = NciCobiList(filename=f"{TEST_FILES_DIR}/cohp/NcICOBILIST.lobster.nospin")
        self.ncicobi_no_spin_wo = NciCobiList(
            filename=f"{TEST_FILES_DIR}/cohp/NcICOBILIST.lobster.nospin.withoutorbitals"
        )
        self.ncicobi_wo = NciCobiList(filename=f"{TEST_FILES_DIR}/cohp/NcICOBILIST.lobster.withoutorbitals")

    def test_ncicobilist(self):
        assert self.ncicobi.is_spin_polarized
        assert not self.ncicobi_no_spin.is_spin_polarized
        assert self.ncicobi_wo.is_spin_polarized
        assert not self.ncicobi_no_spin_wo.is_spin_polarized
        assert self.ncicobi.orbital_wise
        assert self.ncicobi_no_spin.orbital_wise
        assert not self.ncicobi_wo.orbital_wise
        assert not self.ncicobi_no_spin_wo.orbital_wise
        assert len(self.ncicobi.ncicobi_list) == 2
        assert self.ncicobi.ncicobi_list["2"]["number_of_atoms"] == 3
        assert self.ncicobi.ncicobi_list["2"]["ncicobi"][Spin.up] == approx(0.00009)
        assert self.ncicobi.ncicobi_list["2"]["ncicobi"][Spin.down] == approx(0.00009)
        assert self.ncicobi.ncicobi_list["2"]["interaction_type"] == "[X22[0,0,0]->Xs42[0,0,0]->X31[0,0,0]]"
        assert (
            self.ncicobi.ncicobi_list["2"]["ncicobi"][Spin.up] == self.ncicobi_wo.ncicobi_list["2"]["ncicobi"][Spin.up]
        )
        assert (
            self.ncicobi.ncicobi_list["2"]["ncicobi"][Spin.up] == self.ncicobi_gz.ncicobi_list["2"]["ncicobi"][Spin.up]
        )
        assert (
            self.ncicobi.ncicobi_list["2"]["interaction_type"] == self.ncicobi_gz.ncicobi_list["2"]["interaction_type"]
        )
        assert sum(self.ncicobi.ncicobi_list["2"]["ncicobi"].values()) == approx(
            self.ncicobi_no_spin.ncicobi_list["2"]["ncicobi"][Spin.up]
        )


class TestDoscar(unittest.TestCase):
    def setUp(self):
        # first for spin polarized version
        doscar = f"{TEST_FILES_DIR}/DOSCAR.lobster.spin"
        poscar = f"{TEST_FILES_DIR}/POSCAR.lobster.spin_DOS"
        # not spin polarized
        doscar2 = f"{TEST_FILES_DIR}/DOSCAR.lobster.nonspin"
        poscar2 = f"{TEST_FILES_DIR}/POSCAR.lobster.nonspin_DOS"
        f"{TEST_FILES_DIR}/DOSCAR.lobster.nonspin_zip.gz"
        f"{TEST_FILES_DIR}/POSCAR.lobster.nonspin_DOS_zip.gz"
        self.DOSCAR_spin_pol = Doscar(doscar=doscar, structure_file=poscar)
        self.DOSCAR_nonspin_pol = Doscar(doscar=doscar2, structure_file=poscar2)

        self.DOSCAR_spin_pol = Doscar(doscar=doscar, structure_file=poscar)
        self.DOSCAR_nonspin_pol = Doscar(doscar=doscar2, structure_file=poscar2)

        with open(f"{TEST_FILES_DIR}/structure_KF.json") as f:
            data = json.load(f)

        self.structure = Structure.from_dict(data)

        # test structure argument
        self.DOSCAR_spin_pol2 = Doscar(doscar=doscar, structure_file=None, structure=Structure.from_file(poscar))

    def test_complete_dos(self):
        # first for spin polarized version
        energies_spin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        tdos_up = [0.00000, 0.79999, 0.00000, 0.79999, 0.00000, 0.02577]
        tdos_down = [0.00000, 0.79999, 0.00000, 0.79999, 0.00000, 0.02586]
        fermi = 0.0

        pdos_f_2s_up = [0.00000, 0.00159, 0.00000, 0.00011, 0.00000, 0.00069]
        pdos_f_2s_down = [0.00000, 0.00159, 0.00000, 0.00011, 0.00000, 0.00069]
        pdos_f_2py_up = [0.00000, 0.00160, 0.00000, 0.25801, 0.00000, 0.00029]
        pdos_f_2py_down = [0.00000, 0.00161, 0.00000, 0.25819, 0.00000, 0.00029]
        pdos_f_2pz_up = [0.00000, 0.00161, 0.00000, 0.25823, 0.00000, 0.00029]
        pdos_f_2pz_down = [0.00000, 0.00160, 0.00000, 0.25795, 0.00000, 0.00029]
        pdos_f_2px_up = [0.00000, 0.00160, 0.00000, 0.25805, 0.00000, 0.00029]
        pdos_f_2px_down = [0.00000, 0.00161, 0.00000, 0.25814, 0.00000, 0.00029]

        assert energies_spin == self.DOSCAR_spin_pol.completedos.energies.tolist()
        assert tdos_up == self.DOSCAR_spin_pol.completedos.densities[Spin.up].tolist()
        assert tdos_down == self.DOSCAR_spin_pol.completedos.densities[Spin.down].tolist()
        assert fermi == approx(self.DOSCAR_spin_pol.completedos.efermi)

        assert_allclose(
            self.DOSCAR_spin_pol.completedos.structure.frac_coords,
            self.structure.frac_coords,
        )
        assert_allclose(
            self.DOSCAR_spin_pol2.completedos.structure.frac_coords,
            self.structure.frac_coords,
        )
        assert self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2s"][Spin.up].tolist() == pdos_f_2s_up
        assert self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2s"][Spin.down].tolist() == pdos_f_2s_down
        assert self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2p_y"][Spin.up].tolist() == pdos_f_2py_up
        assert self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2p_y"][Spin.down].tolist() == pdos_f_2py_down
        assert self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2p_z"][Spin.up].tolist() == pdos_f_2pz_up
        assert self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2p_z"][Spin.down].tolist() == pdos_f_2pz_down
        assert self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2p_x"][Spin.up].tolist() == pdos_f_2px_up
        assert self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2p_x"][Spin.down].tolist() == pdos_f_2px_down

        energies_nonspin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        tdos_nonspin = [0.00000, 1.60000, 0.00000, 1.60000, 0.00000, 0.02418]
        pdos_f_2s = [0.00000, 0.00320, 0.00000, 0.00017, 0.00000, 0.00060]
        pdos_f_2py = [0.00000, 0.00322, 0.00000, 0.51635, 0.00000, 0.00037]
        pdos_f_2pz = [0.00000, 0.00322, 0.00000, 0.51636, 0.00000, 0.00037]
        pdos_f_2px = [0.00000, 0.00322, 0.00000, 0.51634, 0.00000, 0.00037]

        assert energies_nonspin == self.DOSCAR_nonspin_pol.completedos.energies.tolist()

        assert tdos_nonspin == self.DOSCAR_nonspin_pol.completedos.densities[Spin.up].tolist()

        assert fermi == approx(self.DOSCAR_nonspin_pol.completedos.efermi)

        assert self.DOSCAR_nonspin_pol.completedos.structure == self.structure

        assert self.DOSCAR_nonspin_pol.completedos.pdos[self.structure[0]]["2s"][Spin.up].tolist() == pdos_f_2s
        assert self.DOSCAR_nonspin_pol.completedos.pdos[self.structure[0]]["2p_y"][Spin.up].tolist() == pdos_f_2py
        assert self.DOSCAR_nonspin_pol.completedos.pdos[self.structure[0]]["2p_z"][Spin.up].tolist() == pdos_f_2pz
        assert self.DOSCAR_nonspin_pol.completedos.pdos[self.structure[0]]["2p_x"][Spin.up].tolist() == pdos_f_2px

    def test_pdos(self):
        # first for spin polarized version

        pdos_f_2s_up = [0.00000, 0.00159, 0.00000, 0.00011, 0.00000, 0.00069]
        pdos_f_2s_down = [0.00000, 0.00159, 0.00000, 0.00011, 0.00000, 0.00069]
        pdos_f_2py_up = [0.00000, 0.00160, 0.00000, 0.25801, 0.00000, 0.00029]
        pdos_f_2py_down = [0.00000, 0.00161, 0.00000, 0.25819, 0.00000, 0.00029]
        pdos_f_2pz_up = [0.00000, 0.00161, 0.00000, 0.25823, 0.00000, 0.00029]
        pdos_f_2pz_down = [0.00000, 0.00160, 0.00000, 0.25795, 0.00000, 0.00029]
        pdos_f_2px_up = [0.00000, 0.00160, 0.00000, 0.25805, 0.00000, 0.00029]
        pdos_f_2px_down = [0.00000, 0.00161, 0.00000, 0.25814, 0.00000, 0.00029]

        assert self.DOSCAR_spin_pol.pdos[0]["2s"][Spin.up].tolist() == pdos_f_2s_up
        assert self.DOSCAR_spin_pol.pdos[0]["2s"][Spin.down].tolist() == pdos_f_2s_down
        assert self.DOSCAR_spin_pol.pdos[0]["2p_y"][Spin.up].tolist() == pdos_f_2py_up
        assert self.DOSCAR_spin_pol.pdos[0]["2p_y"][Spin.down].tolist() == pdos_f_2py_down
        assert self.DOSCAR_spin_pol.pdos[0]["2p_z"][Spin.up].tolist() == pdos_f_2pz_up
        assert self.DOSCAR_spin_pol.pdos[0]["2p_z"][Spin.down].tolist() == pdos_f_2pz_down
        assert self.DOSCAR_spin_pol.pdos[0]["2p_x"][Spin.up].tolist() == pdos_f_2px_up
        assert self.DOSCAR_spin_pol.pdos[0]["2p_x"][Spin.down].tolist() == pdos_f_2px_down

        # non spin
        pdos_f_2s = [0.00000, 0.00320, 0.00000, 0.00017, 0.00000, 0.00060]
        pdos_f_2py = [0.00000, 0.00322, 0.00000, 0.51635, 0.00000, 0.00037]
        pdos_f_2pz = [0.00000, 0.00322, 0.00000, 0.51636, 0.00000, 0.00037]
        pdos_f_2px = [0.00000, 0.00322, 0.00000, 0.51634, 0.00000, 0.00037]

        assert self.DOSCAR_nonspin_pol.pdos[0]["2s"][Spin.up].tolist() == pdos_f_2s
        assert self.DOSCAR_nonspin_pol.pdos[0]["2p_y"][Spin.up].tolist() == pdos_f_2py
        assert self.DOSCAR_nonspin_pol.pdos[0]["2p_z"][Spin.up].tolist() == pdos_f_2pz
        assert self.DOSCAR_nonspin_pol.pdos[0]["2p_x"][Spin.up].tolist() == pdos_f_2px

    def test_tdos(self):
        # first for spin polarized version
        energies_spin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        tdos_up = [0.00000, 0.79999, 0.00000, 0.79999, 0.00000, 0.02577]
        tdos_down = [0.00000, 0.79999, 0.00000, 0.79999, 0.00000, 0.02586]
        fermi = 0.0

        assert energies_spin == self.DOSCAR_spin_pol.tdos.energies.tolist()
        assert tdos_up == self.DOSCAR_spin_pol.tdos.densities[Spin.up].tolist()
        assert tdos_down == self.DOSCAR_spin_pol.tdos.densities[Spin.down].tolist()
        assert fermi == approx(self.DOSCAR_spin_pol.tdos.efermi)

        energies_nonspin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        tdos_nonspin = [0.00000, 1.60000, 0.00000, 1.60000, 0.00000, 0.02418]
        fermi = 0.0

        assert energies_nonspin == self.DOSCAR_nonspin_pol.tdos.energies.tolist()
        assert tdos_nonspin == self.DOSCAR_nonspin_pol.tdos.densities[Spin.up].tolist()
        assert fermi == approx(self.DOSCAR_nonspin_pol.tdos.efermi)

    def test_energies(self):
        # first for spin polarized version
        energies_spin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]

        assert energies_spin == self.DOSCAR_spin_pol.energies.tolist()

        energies_nonspin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        assert energies_nonspin == self.DOSCAR_nonspin_pol.energies.tolist()

    def test_tdensities(self):
        # first for spin polarized version
        tdos_up = [0.00000, 0.79999, 0.00000, 0.79999, 0.00000, 0.02577]
        tdos_down = [0.00000, 0.79999, 0.00000, 0.79999, 0.00000, 0.02586]

        assert tdos_up == self.DOSCAR_spin_pol.tdensities[Spin.up].tolist()
        assert tdos_down == self.DOSCAR_spin_pol.tdensities[Spin.down].tolist()

        tdos_nonspin = [0.00000, 1.60000, 0.00000, 1.60000, 0.00000, 0.02418]
        assert tdos_nonspin == self.DOSCAR_nonspin_pol.tdensities[Spin.up].tolist()

    def test_itdensities(self):
        itdos_up = [1.99997, 4.99992, 4.99992, 7.99987, 7.99987, 8.09650]
        itdos_down = [1.99997, 4.99992, 4.99992, 7.99987, 7.99987, 8.09685]
        assert itdos_up == self.DOSCAR_spin_pol.itdensities[Spin.up].tolist()
        assert itdos_down == self.DOSCAR_spin_pol.itdensities[Spin.down].tolist()

        itdos_nonspin = [4.00000, 10.00000, 10.00000, 16.00000, 16.00000, 16.09067]
        assert itdos_nonspin == self.DOSCAR_nonspin_pol.itdensities[Spin.up].tolist()

    def test_is_spin_polarized(self):
        # first for spin polarized version
        assert self.DOSCAR_spin_pol.is_spin_polarized

        assert not self.DOSCAR_nonspin_pol.is_spin_polarized


class TestCharge(PymatgenTest):
    def setUp(self):
        self.charge2 = Charge(filename=f"{TEST_FILES_DIR}/cohp/CHARGE.lobster.MnO")
        # gzipped file
        self.charge = Charge(filename=f"{TEST_FILES_DIR}/cohp/CHARGE.lobster.MnO2.gz")

    def test_attributes(self):
        charge_Loewdin = [-1.25, 1.25]
        charge_Mulliken = [-1.30, 1.30]
        atomlist = ["O1", "Mn2"]
        types = ["O", "Mn"]
        num_atoms = 2
        assert charge_Mulliken == self.charge2.Mulliken
        assert charge_Loewdin == self.charge2.Loewdin
        assert atomlist == self.charge2.atomlist
        assert types == self.charge2.types
        assert num_atoms == self.charge2.num_atoms

    def test_get_structure_with_charges(self):
        structure_dict2 = {
            "lattice": {
                "c": 3.198244,
                "volume": 23.132361565928807,
                "b": 3.1982447183003364,
                "gamma": 60.00000011873414,
                "beta": 60.00000401737447,
                "alpha": 60.00000742944491,
                "matrix": [
                    [2.769761, 0.0, 1.599122],
                    [0.923254, 2.611356, 1.599122],
                    [0.0, 0.0, 3.198244],
                ],
                "a": 3.1982443884113985,
            },
            "@class": "Structure",
            "sites": [
                {
                    "xyz": [1.846502883732, 1.305680611356, 3.198248797366],
                    "properties": {"Loewdin Charges": -1.25, "Mulliken Charges": -1.3},
                    "abc": [0.499998, 0.500001, 0.500002],
                    "species": [{"occu": 1, "element": "O"}],
                    "label": "O",
                },
                {
                    "xyz": [0.0, 0.0, 0.0],
                    "properties": {"Loewdin Charges": 1.25, "Mulliken Charges": 1.3},
                    "abc": [0.0, 0.0, 0.0],
                    "species": [{"occu": 1, "element": "Mn"}],
                    "label": "Mn",
                },
            ],
            "charge": None,
            "@module": "pymatgen.core.structure",
        }
        s2 = Structure.from_dict(structure_dict2)
        assert s2 == self.charge2.get_structure_with_charges(f"{TEST_FILES_DIR}/POSCAR.MnO")


class TestLobsterout(PymatgenTest):
    def setUp(self):
        self.lobsterout_normal = Lobsterout(filename=f"{TEST_FILES_DIR}/cohp/lobsterout.normal")
        # make sure .gz files are also read correctly
        self.lobsterout_normal = Lobsterout(filename=f"{TEST_FILES_DIR}/cohp/lobsterout.normal2.gz")
        self.lobsterout_fatband_grosspop_densityofenergies = Lobsterout(
            filename=os.path.join(
                TEST_FILES_DIR,
                "cohp",
                "lobsterout.fatband_grosspop_densityofenergy",
            )
        )
        self.lobsterout_saveprojection = Lobsterout(filename=f"{TEST_FILES_DIR}/cohp/lobsterout.saveprojection")
        self.lobsterout_skipping_all = Lobsterout(filename=f"{TEST_FILES_DIR}/cohp/lobsterout.skipping_all")
        self.lobsterout_twospins = Lobsterout(filename=f"{TEST_FILES_DIR}/cohp/lobsterout.twospins")
        self.lobsterout_GaAs = Lobsterout(filename=f"{TEST_FILES_DIR}/cohp/lobsterout.GaAs")
        self.lobsterout_from_projection = Lobsterout(filename=f"{TEST_FILES_DIR}/cohp/lobsterout_from_projection")
        self.lobsterout_onethread = Lobsterout(filename=f"{TEST_FILES_DIR}/cohp/lobsterout.onethread")
        self.lobsterout_cobi_madelung = Lobsterout(filename=f"{TEST_FILES_DIR}/cohp/lobsterout_cobi_madelung")
        self.lobsterout_doscar_lso = Lobsterout(filename=f"{TEST_FILES_DIR}/cohp/lobsterout_doscar_lso")

        # TODO: implement skipping madelung/cobi
        self.lobsterout_skipping_cobi_madelung = Lobsterout(
            filename=f"{TEST_FILES_DIR}/cohp/lobsterout.skip_cobi_madelung"
        )

    def test_attributes(self):
        assert self.lobsterout_normal.basis_functions == [
            ["3s", "4s", "3p_y", "3p_z", "3p_x", "3d_xy", "3d_yz", "3d_z^2", "3d_xz", "3d_x^2-y^2"]
        ]
        assert self.lobsterout_normal.basis_type == ["pbeVaspFit2015"]
        assert self.lobsterout_normal.charge_spilling == [0.0268]
        assert self.lobsterout_normal.dft_program == "VASP"
        assert self.lobsterout_normal.elements == ["Ti"]
        assert self.lobsterout_normal.has_charge
        assert self.lobsterout_normal.has_cohpcar
        assert self.lobsterout_normal.has_coopcar
        assert self.lobsterout_normal.has_doscar
        assert not self.lobsterout_normal.has_projection
        assert self.lobsterout_normal.has_bandoverlaps
        assert not self.lobsterout_normal.has_density_of_energies
        assert not self.lobsterout_normal.has_fatbands
        assert not self.lobsterout_normal.has_grosspopulation
        assert self.lobsterout_normal.info_lines == [
            "There are more PAW bands than local basis functions available.",
            "To prevent trouble in orthonormalization and Hamiltonian reconstruction",
            "the PAW bands from 21 and upwards will be ignored.",
        ]
        assert self.lobsterout_normal.info_orthonormalization == [
            "3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5."
        ]
        assert not self.lobsterout_normal.is_restart_from_projection
        assert self.lobsterout_normal.lobster_version == "v3.1.0"
        assert self.lobsterout_normal.number_of_spins == 1
        assert self.lobsterout_normal.number_of_threads == 8
        assert self.lobsterout_normal.timing == {
            "wall_time": {"h": "0", "min": "0", "s": "2", "ms": "702"},
            "user_time": {"h": "0", "min": "0", "s": "20", "ms": "330"},
            "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "310"},
        }
        assert self.lobsterout_normal.total_spilling[0] == approx([0.044000000000000004][0])
        assert self.lobsterout_normal.warning_lines == [
            "3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5.",
            "Generally, this is not a critical error. But to help you analyze it,",
            "I dumped the band overlap matrices to the file bandOverlaps.lobster.",
            "Please check how much they deviate from the identity matrix and decide to",
            "use your results only, if you are sure that this is ok.",
        ]

        assert self.lobsterout_fatband_grosspop_densityofenergies.basis_functions == [
            ["3s", "4s", "3p_y", "3p_z", "3p_x", "3d_xy", "3d_yz", "3d_z^2", "3d_xz", "3d_x^2-y^2"]
        ]
        assert self.lobsterout_fatband_grosspop_densityofenergies.basis_type == ["pbeVaspFit2015"]
        assert self.lobsterout_fatband_grosspop_densityofenergies.charge_spilling == [0.0268]
        assert self.lobsterout_fatband_grosspop_densityofenergies.dft_program == "VASP"
        assert self.lobsterout_fatband_grosspop_densityofenergies.elements == ["Ti"]
        assert self.lobsterout_fatband_grosspop_densityofenergies.has_charge
        assert not self.lobsterout_fatband_grosspop_densityofenergies.has_cohpcar
        assert not self.lobsterout_fatband_grosspop_densityofenergies.has_coopcar
        assert not self.lobsterout_fatband_grosspop_densityofenergies.has_doscar
        assert not self.lobsterout_fatband_grosspop_densityofenergies.has_projection
        assert self.lobsterout_fatband_grosspop_densityofenergies.has_bandoverlaps
        assert self.lobsterout_fatband_grosspop_densityofenergies.has_density_of_energies
        assert self.lobsterout_fatband_grosspop_densityofenergies.has_fatbands
        assert self.lobsterout_fatband_grosspop_densityofenergies.has_grosspopulation
        assert self.lobsterout_fatband_grosspop_densityofenergies.info_lines == [
            "There are more PAW bands than local basis functions available.",
            "To prevent trouble in orthonormalization and Hamiltonian reconstruction",
            "the PAW bands from 21 and upwards will be ignored.",
        ]
        assert self.lobsterout_fatband_grosspop_densityofenergies.info_orthonormalization == [
            "3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5."
        ]
        assert not self.lobsterout_fatband_grosspop_densityofenergies.is_restart_from_projection
        assert self.lobsterout_fatband_grosspop_densityofenergies.lobster_version == "v3.1.0"
        assert self.lobsterout_fatband_grosspop_densityofenergies.number_of_spins == 1
        assert self.lobsterout_fatband_grosspop_densityofenergies.number_of_threads == 8
        assert self.lobsterout_fatband_grosspop_densityofenergies.timing == {
            "wall_time": {"h": "0", "min": "0", "s": "4", "ms": "136"},
            "user_time": {"h": "0", "min": "0", "s": "18", "ms": "280"},
            "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "290"},
        }
        assert self.lobsterout_fatband_grosspop_densityofenergies.total_spilling[0] == approx([0.044000000000000004][0])
        assert self.lobsterout_fatband_grosspop_densityofenergies.warning_lines == [
            "3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5.",
            "Generally, this is not a critical error. But to help you analyze it,",
            "I dumped the band overlap matrices to the file bandOverlaps.lobster.",
            "Please check how much they deviate from the identity matrix and decide to",
            "use your results only, if you are sure that this is ok.",
        ]

        assert self.lobsterout_saveprojection.basis_functions == [
            [
                "3s",
                "4s",
                "3p_y",
                "3p_z",
                "3p_x",
                "3d_xy",
                "3d_yz",
                "3d_z^2",
                "3d_xz",
                "3d_x^2-y^2",
            ]
        ]
        assert self.lobsterout_saveprojection.basis_type == ["pbeVaspFit2015"]
        assert self.lobsterout_saveprojection.charge_spilling == [0.0268]
        assert self.lobsterout_saveprojection.dft_program == "VASP"
        assert self.lobsterout_saveprojection.elements == ["Ti"]
        assert self.lobsterout_saveprojection.has_charge
        assert not self.lobsterout_saveprojection.has_cohpcar
        assert not self.lobsterout_saveprojection.has_coopcar
        assert not self.lobsterout_saveprojection.has_doscar
        assert self.lobsterout_saveprojection.has_projection
        assert self.lobsterout_saveprojection.has_bandoverlaps
        assert self.lobsterout_saveprojection.has_density_of_energies
        assert not self.lobsterout_saveprojection.has_fatbands
        assert not self.lobsterout_saveprojection.has_grosspopulation
        assert self.lobsterout_saveprojection.info_lines == [
            "There are more PAW bands than local basis functions available.",
            "To prevent trouble in orthonormalization and Hamiltonian reconstruction",
            "the PAW bands from 21 and upwards will be ignored.",
        ]
        assert self.lobsterout_saveprojection.info_orthonormalization == [
            "3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5."
        ]
        assert not self.lobsterout_saveprojection.is_restart_from_projection
        assert self.lobsterout_saveprojection.lobster_version == "v3.1.0"
        assert self.lobsterout_saveprojection.number_of_spins == 1
        assert self.lobsterout_saveprojection.number_of_threads == 8
        assert self.lobsterout_saveprojection.timing == {
            "wall_time": {"h": "0", "min": "0", "s": "2", "ms": "574"},
            "user_time": {"h": "0", "min": "0", "s": "18", "ms": "250"},
            "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "320"},
        }
        assert self.lobsterout_saveprojection.total_spilling[0] == approx([0.044000000000000004][0])
        assert self.lobsterout_saveprojection.warning_lines == [
            "3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5.",
            "Generally, this is not a critical error. But to help you analyze it,",
            "I dumped the band overlap matrices to the file bandOverlaps.lobster.",
            "Please check how much they deviate from the identity matrix and decide to",
            "use your results only, if you are sure that this is ok.",
        ]

        assert self.lobsterout_skipping_all.basis_functions == [
            [
                "3s",
                "4s",
                "3p_y",
                "3p_z",
                "3p_x",
                "3d_xy",
                "3d_yz",
                "3d_z^2",
                "3d_xz",
                "3d_x^2-y^2",
            ]
        ]
        assert self.lobsterout_skipping_all.basis_type == ["pbeVaspFit2015"]
        assert self.lobsterout_skipping_all.charge_spilling == [0.0268]
        assert self.lobsterout_skipping_all.dft_program == "VASP"
        assert self.lobsterout_skipping_all.elements == ["Ti"]
        assert not self.lobsterout_skipping_all.has_charge
        assert not self.lobsterout_skipping_all.has_cohpcar
        assert not self.lobsterout_skipping_all.has_coopcar
        assert not self.lobsterout_skipping_all.has_doscar
        assert not self.lobsterout_skipping_all.has_projection
        assert self.lobsterout_skipping_all.has_bandoverlaps
        assert not self.lobsterout_skipping_all.has_density_of_energies
        assert not self.lobsterout_skipping_all.has_fatbands
        assert not self.lobsterout_skipping_all.has_grosspopulation
        assert not self.lobsterout_skipping_all.has_cobicar
        assert not self.lobsterout_skipping_all.has_madelung
        assert self.lobsterout_skipping_all.info_lines == [
            "There are more PAW bands than local basis functions available.",
            "To prevent trouble in orthonormalization and Hamiltonian reconstruction",
            "the PAW bands from 21 and upwards will be ignored.",
        ]
        assert self.lobsterout_skipping_all.info_orthonormalization == [
            "3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5."
        ]
        assert not self.lobsterout_skipping_all.is_restart_from_projection
        assert self.lobsterout_skipping_all.lobster_version == "v3.1.0"
        assert self.lobsterout_skipping_all.number_of_spins == 1
        assert self.lobsterout_skipping_all.number_of_threads == 8
        assert self.lobsterout_skipping_all.timing == {
            "wall_time": {"h": "0", "min": "0", "s": "2", "ms": "117"},
            "user_time": {"h": "0", "min": "0", "s": "16", "ms": "79"},
            "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "320"},
        }
        assert self.lobsterout_skipping_all.total_spilling[0] == approx([0.044000000000000004][0])
        assert self.lobsterout_skipping_all.warning_lines == [
            "3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5.",
            "Generally, this is not a critical error. But to help you analyze it,",
            "I dumped the band overlap matrices to the file bandOverlaps.lobster.",
            "Please check how much they deviate from the identity matrix and decide to",
            "use your results only, if you are sure that this is ok.",
        ]

        assert self.lobsterout_twospins.basis_functions == [
            [
                "4s",
                "4p_y",
                "4p_z",
                "4p_x",
                "3d_xy",
                "3d_yz",
                "3d_z^2",
                "3d_xz",
                "3d_x^2-y^2",
            ]
        ]
        assert self.lobsterout_twospins.basis_type == ["pbeVaspFit2015"]
        assert self.lobsterout_twospins.charge_spilling[0] == approx(0.36619999999999997)
        assert self.lobsterout_twospins.charge_spilling[1] == approx(0.36619999999999997)
        assert self.lobsterout_twospins.dft_program == "VASP"
        assert self.lobsterout_twospins.elements == ["Ti"]
        assert self.lobsterout_twospins.has_charge
        assert self.lobsterout_twospins.has_cohpcar
        assert self.lobsterout_twospins.has_coopcar
        assert self.lobsterout_twospins.has_doscar
        assert not self.lobsterout_twospins.has_projection
        assert self.lobsterout_twospins.has_bandoverlaps
        assert not self.lobsterout_twospins.has_density_of_energies
        assert not self.lobsterout_twospins.has_fatbands
        assert not self.lobsterout_twospins.has_grosspopulation
        assert self.lobsterout_twospins.info_lines == [
            "There are more PAW bands than local basis functions available.",
            "To prevent trouble in orthonormalization and Hamiltonian reconstruction",
            "the PAW bands from 19 and upwards will be ignored.",
        ]
        assert self.lobsterout_twospins.info_orthonormalization == [
            "60 of 294 k-points could not be orthonormalized with an accuracy of 1.0E-5."
        ]
        assert not self.lobsterout_twospins.is_restart_from_projection
        assert self.lobsterout_twospins.lobster_version == "v3.1.0"
        assert self.lobsterout_twospins.number_of_spins == 2
        assert self.lobsterout_twospins.number_of_threads == 8
        assert self.lobsterout_twospins.timing == {
            "wall_time": {"h": "0", "min": "0", "s": "3", "ms": "71"},
            "user_time": {"h": "0", "min": "0", "s": "22", "ms": "660"},
            "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "310"},
        }
        assert self.lobsterout_twospins.total_spilling[0] == approx([0.2567][0])
        assert self.lobsterout_twospins.total_spilling[1] == approx([0.2567][0])
        assert self.lobsterout_twospins.warning_lines == [
            "60 of 294 k-points could not be orthonormalized with an accuracy of 1.0E-5.",
            "Generally, this is not a critical error. But to help you analyze it,",
            "I dumped the band overlap matrices to the file bandOverlaps.lobster.",
            "Please check how much they deviate from the identity matrix and decide to",
            "use your results only, if you are sure that this is ok.",
        ]

        assert self.lobsterout_from_projection.basis_functions == []
        assert self.lobsterout_from_projection.basis_type == []
        assert self.lobsterout_from_projection.charge_spilling[0] == approx(0.0177)
        assert self.lobsterout_from_projection.dft_program is None
        assert self.lobsterout_from_projection.elements == []
        assert self.lobsterout_from_projection.has_charge
        assert self.lobsterout_from_projection.has_cohpcar
        assert self.lobsterout_from_projection.has_coopcar
        assert self.lobsterout_from_projection.has_doscar
        assert not self.lobsterout_from_projection.has_projection
        assert not self.lobsterout_from_projection.has_bandoverlaps
        assert not self.lobsterout_from_projection.has_density_of_energies
        assert not self.lobsterout_from_projection.has_fatbands
        assert not self.lobsterout_from_projection.has_grosspopulation
        assert self.lobsterout_from_projection.info_lines == []
        assert self.lobsterout_from_projection.info_orthonormalization == []
        assert self.lobsterout_from_projection.is_restart_from_projection
        assert self.lobsterout_from_projection.lobster_version == "v3.1.0"
        assert self.lobsterout_from_projection.number_of_spins == 1
        assert self.lobsterout_from_projection.number_of_threads == 8
        assert self.lobsterout_from_projection.timing == {
            "wall_time": {"h": "0", "min": "2", "s": "1", "ms": "890"},
            "user_time": {"h": "0", "min": "15", "s": "10", "ms": "530"},
            "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "400"},
        }
        assert self.lobsterout_from_projection.total_spilling[0] == approx([0.1543][0])
        assert self.lobsterout_from_projection.warning_lines == []

        assert self.lobsterout_GaAs.basis_functions == [
            ["4s", "4p_y", "4p_z", "4p_x"],
            [
                "4s",
                "4p_y",
                "4p_z",
                "4p_x",
                "3d_xy",
                "3d_yz",
                "3d_z^2",
                "3d_xz",
                "3d_x^2-y^2",
            ],
        ]
        assert self.lobsterout_GaAs.basis_type == ["Bunge", "Bunge"]
        assert self.lobsterout_GaAs.charge_spilling[0] == approx(0.0089)
        assert self.lobsterout_GaAs.dft_program == "VASP"
        assert self.lobsterout_GaAs.elements == ["As", "Ga"]
        assert self.lobsterout_GaAs.has_charge
        assert self.lobsterout_GaAs.has_cohpcar
        assert self.lobsterout_GaAs.has_coopcar
        assert self.lobsterout_GaAs.has_doscar
        assert not self.lobsterout_GaAs.has_projection
        assert not self.lobsterout_GaAs.has_bandoverlaps
        assert not self.lobsterout_GaAs.has_density_of_energies
        assert not self.lobsterout_GaAs.has_fatbands
        assert not self.lobsterout_GaAs.has_grosspopulation
        assert self.lobsterout_GaAs.info_lines == [
            "There are more PAW bands than local basis functions available.",
            "To prevent trouble in orthonormalization and Hamiltonian reconstruction",
            "the PAW bands from 14 and upwards will be ignored.",
        ]
        assert self.lobsterout_GaAs.info_orthonormalization == []
        assert not self.lobsterout_GaAs.is_restart_from_projection
        assert self.lobsterout_GaAs.lobster_version == "v3.1.0"
        assert self.lobsterout_GaAs.number_of_spins == 1
        assert self.lobsterout_GaAs.number_of_threads == 8
        assert self.lobsterout_GaAs.timing == {
            "wall_time": {"h": "0", "min": "0", "s": "2", "ms": "726"},
            "user_time": {"h": "0", "min": "0", "s": "12", "ms": "370"},
            "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "180"},
        }
        assert self.lobsterout_GaAs.total_spilling[0] == approx([0.0859][0])

        assert self.lobsterout_onethread.number_of_threads == 1
        # Test lobsterout of lobster-4.1.0
        assert self.lobsterout_cobi_madelung.has_cobicar
        assert self.lobsterout_cobi_madelung.has_cohpcar
        assert self.lobsterout_cobi_madelung.has_madelung
        assert not self.lobsterout_cobi_madelung.has_doscar_lso

        assert self.lobsterout_doscar_lso.has_doscar_lso

        assert self.lobsterout_skipping_cobi_madelung.has_cobicar is False
        assert self.lobsterout_skipping_cobi_madelung.has_madelung is False

    def test_get_doc(self):
        comparedict = {
            "restart_from_projection": False,
            "lobster_version": "v3.1.0",
            "threads": 8,
            "dft_program": "VASP",
            "charge_spilling": [0.0268],
            "total_spilling": [0.044000000000000004],
            "elements": ["Ti"],
            "basis_type": ["pbeVaspFit2015"],
            "basis_functions": [
                [
                    "3s",
                    "4s",
                    "3p_y",
                    "3p_z",
                    "3p_x",
                    "3d_xy",
                    "3d_yz",
                    "3d_z^2",
                    "3d_xz",
                    "3d_x^2-y^2",
                ]
            ],
            "timing": {
                "wall_time": {"h": "0", "min": "0", "s": "2", "ms": "702"},
                "user_time": {"h": "0", "min": "0", "s": "20", "ms": "330"},
                "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "310"},
            },
            "warning_lines": [
                "3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5.",
                "Generally, this is not a critical error. But to help you analyze it,",
                "I dumped the band overlap matrices to the file bandOverlaps.lobster.",
                "Please check how much they deviate from the identity matrix and decide to",
                "use your results only, if you are sure that this is ok.",
            ],
            "info_orthonormalization": ["3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5."],
            "info_lines": [
                "There are more PAW bands than local basis functions available.",
                "To prevent trouble in orthonormalization and Hamiltonian reconstruction",
                "the PAW bands from 21 and upwards will be ignored.",
            ],
            "has_doscar": True,
            "has_doscar_lso": False,
            "has_cohpcar": True,
            "has_coopcar": True,
            "has_charge": True,
            "has_projection": False,
            "has_bandoverlaps": True,
            "has_fatbands": False,
            "has_grosspopulation": False,
            "has_density_of_energies": False,
        }
        for key, item in self.lobsterout_normal.get_doc().items():
            if key not in ["has_cobicar", "has_madelung"]:
                if isinstance(item, str):
                    assert comparedict[key], item
                elif isinstance(item, int):
                    assert comparedict[key] == item
                elif key in ("charge_spilling", "total_spilling"):
                    assert item[0] == approx(comparedict[key][0])
                elif isinstance(item, (list, dict)):
                    assert item == comparedict[key]


class TestFatband(PymatgenTest):
    def setUp(self):
        self.fatband_SiO2_p_x = Fatband(
            filenames=f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p_x",
            Kpointsfile=f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p_x/KPOINTS",
            vasprun=f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p_x/vasprun.xml",
        )
        self.vasprun_SiO2_p_x = Vasprun(filename=f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p_x/vasprun.xml")
        self.bs_symmline = self.vasprun_SiO2_p_x.get_band_structure(line_mode=True, force_hybrid_mode=True)
        self.fatband_SiO2_p = Fatband(
            filenames=f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p",
            Kpointsfile=f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p/KPOINTS",
            vasprun=f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p/vasprun.xml",
        )
        self.vasprun_SiO2_p = Vasprun(filename=f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p/vasprun.xml")
        self.bs_symmline2 = self.vasprun_SiO2_p.get_band_structure(line_mode=True, force_hybrid_mode=True)
        self.fatband_SiO2_spin = Fatband(
            filenames=f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_Spin",
            Kpointsfile=f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_Spin/KPOINTS",
            vasprun=os.path.join(
                TEST_FILES_DIR,
                "cohp",
                "Fatband_SiO2/Test_Spin/vasprun.xml",
            ),
        )
        self.vasprun_SiO2_spin = Vasprun(
            filename=os.path.join(
                TEST_FILES_DIR,
                "cohp",
                "Fatband_SiO2/Test_Spin/vasprun.xml",
            )
        )
        self.bs_symmline_spin = self.vasprun_SiO2_p.get_band_structure(line_mode=True, force_hybrid_mode=True)

    def test_attributes(self):
        assert list(self.fatband_SiO2_p_x.label_dict["M"]) == approx([0.5, 0.0, 0.0])
        assert self.fatband_SiO2_p_x.efermi == self.vasprun_SiO2_p_x.efermi
        lattice1 = self.bs_symmline.lattice_rec.as_dict()
        lattice2 = self.fatband_SiO2_p_x.lattice.as_dict()
        for idx in range(3):
            assert lattice1["matrix"][idx] == approx(lattice2["matrix"][idx])
        assert self.fatband_SiO2_p_x.eigenvals[Spin.up][1][1] - self.fatband_SiO2_p_x.efermi == -18.245
        assert self.fatband_SiO2_p_x.is_spinpolarized is False
        assert self.fatband_SiO2_p_x.kpoints_array[3] == approx([0.03409091, 0, 0])
        assert self.fatband_SiO2_p_x.nbands == 36
        assert self.fatband_SiO2_p_x.p_eigenvals[Spin.up][2][1]["Si1"]["3p_x"] == 0.002
        assert self.fatband_SiO2_p_x.structure[0].frac_coords == approx([0.0, 0.47634315, 0.666667])
        assert self.fatband_SiO2_p_x.structure[0].species_string == "Si"
        assert self.fatband_SiO2_p_x.structure[0].coords == approx([-1.19607309, 2.0716597, 3.67462144])

        assert list(self.fatband_SiO2_p.label_dict["M"]) == approx([0.5, 0.0, 0.0])
        assert self.fatband_SiO2_p.efermi == self.vasprun_SiO2_p.efermi
        lattice1 = self.bs_symmline2.lattice_rec.as_dict()
        lattice2 = self.fatband_SiO2_p.lattice.as_dict()
        for idx in range(3):
            assert lattice1["matrix"][idx] == approx(lattice2["matrix"][idx])
        assert self.fatband_SiO2_p.eigenvals[Spin.up][1][1] - self.fatband_SiO2_p.efermi == -18.245
        assert self.fatband_SiO2_p.is_spinpolarized is False
        assert self.fatband_SiO2_p.kpoints_array[3] == approx([0.03409091, 0, 0])
        assert self.fatband_SiO2_p.nbands == 36
        assert self.fatband_SiO2_p.p_eigenvals[Spin.up][2][1]["Si1"]["3p"] == 0.042
        assert self.fatband_SiO2_p.structure[0].frac_coords == approx([0.0, 0.47634315, 0.666667])
        assert self.fatband_SiO2_p.structure[0].species_string == "Si"
        assert self.fatband_SiO2_p.structure[0].coords == approx([-1.19607309, 2.0716597, 3.67462144])

        assert list(self.fatband_SiO2_spin.label_dict["M"]) == approx([0.5, 0.0, 0.0])
        assert self.fatband_SiO2_spin.efermi == self.vasprun_SiO2_spin.efermi
        lattice1 = self.bs_symmline_spin.lattice_rec.as_dict()
        lattice2 = self.fatband_SiO2_spin.lattice.as_dict()
        for idx in range(3):
            assert lattice1["matrix"][idx] == approx(lattice2["matrix"][idx])
        assert self.fatband_SiO2_spin.eigenvals[Spin.up][1][1] - self.fatband_SiO2_spin.efermi == -18.245
        assert self.fatband_SiO2_spin.eigenvals[Spin.down][1][1] - self.fatband_SiO2_spin.efermi == -18.245
        assert self.fatband_SiO2_spin.is_spinpolarized
        assert self.fatband_SiO2_spin.kpoints_array[3] == approx([0.03409091, 0, 0])
        assert self.fatband_SiO2_spin.nbands == 36

        assert self.fatband_SiO2_spin.p_eigenvals[Spin.up][2][1]["Si1"]["3p"] == 0.042
        assert self.fatband_SiO2_spin.structure[0].frac_coords == approx([0.0, 0.47634315, 0.666667])
        assert self.fatband_SiO2_spin.structure[0].species_string == "Si"
        assert self.fatband_SiO2_spin.structure[0].coords == approx([-1.19607309, 2.0716597, 3.67462144])

    def test_raises(self):
        with pytest.raises(
            ValueError, match="The are two FATBAND files for the same atom and orbital. The program will stop"
        ):
            self.fatband_SiO2_p_x = Fatband(
                filenames=[
                    f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p_x/FATBAND_si1_3p_x.lobster",
                    f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p_x/FATBAND_si1_3p_x.lobster",
                ],
                Kpointsfile=f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p_x/KPOINTS",
                vasprun=f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p_x/vasprun.xml",
            )

        with pytest.raises(
            ValueError,
            match=r"Make sure all relevant orbitals were generated and that no duplicates \(2p and 2p_x\) are present",
        ):
            self.fatband_SiO2_p_x = Fatband(
                filenames=[
                    f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p_x/FATBAND_si1_3p_x.lobster",
                    f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p/FATBAND_si1_3p.lobster",
                ],
                Kpointsfile=f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p_x/KPOINTS",
                vasprun=f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p_x/vasprun.xml",
            )

        with pytest.raises(ValueError, match="No FATBAND files in folder or given"):
            self.fatband_SiO2_p_x = Fatband(
                filenames=".",
                Kpointsfile=f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p_x/KPOINTS",
                vasprun=f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p_x/vasprun.xml",
            )

    def test_get_bandstructure(self):
        bs_p = self.fatband_SiO2_p.get_bandstructure()
        atom1 = bs_p.structure[0]
        atom2 = self.bs_symmline2.structure[0]
        assert atom1.frac_coords[0] == approx(atom2.frac_coords[0])
        assert atom1.frac_coords[1] == approx(atom2.frac_coords[1])
        assert atom1.frac_coords[2] == approx(atom2.frac_coords[2])
        assert atom1.coords[0] == approx(atom2.coords[0])
        assert atom1.coords[1] == approx(atom2.coords[1])
        assert atom1.coords[2] == approx(atom2.coords[2])
        assert atom1.species_string == atom2.species_string
        assert bs_p.efermi == self.bs_symmline2.efermi
        branch1 = bs_p.branches[0]
        branch2 = self.bs_symmline2.branches[0]
        assert branch2["name"] == branch1["name"]
        assert branch2["start_index"] == branch1["start_index"]
        assert branch2["end_index"] == branch1["end_index"]

        assert bs_p.distance[30] == approx(self.bs_symmline2.distance[30])
        lattice1 = bs_p.lattice_rec.as_dict()
        lattice2 = self.bs_symmline2.lattice_rec.as_dict()
        assert lattice1["matrix"][0] == approx(lattice2["matrix"][0])
        assert lattice1["matrix"][1] == approx(lattice2["matrix"][1])
        assert lattice1["matrix"][2] == approx(lattice2["matrix"][2])

        assert bs_p.kpoints[8].frac_coords[0] == approx(self.bs_symmline2.kpoints[8].frac_coords[0])
        assert bs_p.kpoints[8].frac_coords[1] == approx(self.bs_symmline2.kpoints[8].frac_coords[1])
        assert bs_p.kpoints[8].frac_coords[2] == approx(self.bs_symmline2.kpoints[8].frac_coords[2])
        assert bs_p.kpoints[8].cart_coords[0] == approx(self.bs_symmline2.kpoints[8].cart_coords[0])
        assert bs_p.kpoints[8].cart_coords[1] == approx(self.bs_symmline2.kpoints[8].cart_coords[1])
        assert bs_p.kpoints[8].cart_coords[2] == approx(self.bs_symmline2.kpoints[8].cart_coords[2])
        assert bs_p.kpoints[50].frac_coords[0] == approx(self.bs_symmline2.kpoints[50].frac_coords[0])
        assert bs_p.kpoints[50].frac_coords[1] == approx(self.bs_symmline2.kpoints[50].frac_coords[1])
        assert bs_p.kpoints[50].frac_coords[2] == approx(self.bs_symmline2.kpoints[50].frac_coords[2])
        assert bs_p.kpoints[50].cart_coords[0] == approx(self.bs_symmline2.kpoints[50].cart_coords[0])
        assert bs_p.kpoints[50].cart_coords[1] == approx(self.bs_symmline2.kpoints[50].cart_coords[1])
        assert bs_p.kpoints[50].cart_coords[2] == approx(self.bs_symmline2.kpoints[50].cart_coords[2])
        assert bs_p.get_band_gap()["energy"] == approx(self.bs_symmline2.get_band_gap()["energy"], abs=1e-2)
        assert bs_p.get_projection_on_elements()[Spin.up][0][0]["Si"] == approx(3 * (0.001 + 0.064))
        assert bs_p.get_projections_on_elements_and_orbitals({"Si": ["3p"]})[Spin.up][0][0]["Si"]["3p"] == approx(0.003)
        assert bs_p.get_projections_on_elements_and_orbitals({"O": ["2p"]})[Spin.up][0][0]["O"]["2p"] == approx(
            0.002 * 3 + 0.003 * 3
        )
        dict_here = bs_p.get_projections_on_elements_and_orbitals({"Si": ["3s", "3p"], "O": ["2s", "2p"]})[Spin.up][0][
            0
        ]
        assert dict_here["Si"]["3s"] == approx(0.192)
        assert dict_here["Si"]["3p"] == approx(0.003)
        assert dict_here["O"]["2s"] == approx(0.792)
        assert dict_here["O"]["2p"] == approx(0.015)

        bs_spin = self.fatband_SiO2_spin.get_bandstructure()
        assert bs_spin.get_projection_on_elements()[Spin.up][0][0]["Si"] == approx(3 * (0.001 + 0.064))
        assert bs_spin.get_projections_on_elements_and_orbitals({"Si": ["3p"]})[Spin.up][0][0]["Si"]["3p"] == approx(
            0.003
        )
        assert bs_spin.get_projections_on_elements_and_orbitals({"O": ["2p"]})[Spin.up][0][0]["O"]["2p"] == approx(
            0.002 * 3 + 0.003 * 3
        )
        dict_here = bs_spin.get_projections_on_elements_and_orbitals({"Si": ["3s", "3p"], "O": ["2s", "2p"]})[Spin.up][
            0
        ][0]
        assert dict_here["Si"]["3s"] == approx(0.192)
        assert dict_here["Si"]["3p"] == approx(0.003)
        assert dict_here["O"]["2s"] == approx(0.792)
        assert dict_here["O"]["2p"] == approx(0.015)

        assert bs_spin.get_projection_on_elements()[Spin.up][0][0]["Si"] == approx(3 * (0.001 + 0.064))
        assert bs_spin.get_projections_on_elements_and_orbitals({"Si": ["3p"]})[Spin.down][0][0]["Si"]["3p"] == approx(
            0.003
        )
        assert bs_spin.get_projections_on_elements_and_orbitals({"O": ["2p"]})[Spin.down][0][0]["O"]["2p"] == approx(
            0.002 * 3 + 0.003 * 3
        )
        dict_here = bs_spin.get_projections_on_elements_and_orbitals({"Si": ["3s", "3p"], "O": ["2s", "2p"]})[
            Spin.down
        ][0][0]
        assert dict_here["Si"]["3s"] == approx(0.192)
        assert dict_here["Si"]["3p"] == approx(0.003)
        assert dict_here["O"]["2s"] == approx(0.792)
        assert dict_here["O"]["2p"] == approx(0.015)
        bs_p_x = self.fatband_SiO2_p_x.get_bandstructure()
        assert bs_p_x.get_projection_on_elements()[Spin.up][0][0]["Si"] == approx(3 * (0.001 + 0.064), abs=1e-2)


class TestLobsterin(unittest.TestCase):
    def setUp(self):
        self.Lobsterinfromfile = Lobsterin.from_file(f"{TEST_FILES_DIR}/cohp/lobsterin.1")
        self.Lobsterinfromfile2 = Lobsterin.from_file(f"{TEST_FILES_DIR}/cohp/lobsterin.2")
        self.Lobsterinfromfile3 = Lobsterin.from_file(f"{TEST_FILES_DIR}/cohp/lobsterin.3")
        self.Lobsterinfromfile4 = Lobsterin.from_file(f"{TEST_FILES_DIR}/cohp/lobsterin.4.gz")

    def test_from_file(self):
        # test read from file
        assert self.Lobsterinfromfile["cohpstartenergy"] == approx(-15.0)
        assert self.Lobsterinfromfile["cohpendenergy"] == approx(5.0)
        assert self.Lobsterinfromfile["basisset"] == "pbeVaspFit2015"
        assert self.Lobsterinfromfile["gaussiansmearingwidth"] == approx(0.1)
        assert self.Lobsterinfromfile["basisfunctions"][0] == "Fe 3d 4p 4s"
        assert self.Lobsterinfromfile["basisfunctions"][1] == "Co 3d 4p 4s"
        assert self.Lobsterinfromfile["skipdos"]
        assert self.Lobsterinfromfile["skipcohp"]
        assert self.Lobsterinfromfile["skipcoop"]
        assert self.Lobsterinfromfile["skippopulationanalysis"]
        assert self.Lobsterinfromfile["skipgrosspopulation"]

        # test if comments are correctly removed
        assert self.Lobsterinfromfile == self.Lobsterinfromfile2

    def test_getitem(self):
        # tests implementation of getitem, should be case independent
        assert self.Lobsterinfromfile["COHPSTARTENERGY"] == approx(-15.0)

    def test_setitem(self):
        # test implementation of setitem
        self.Lobsterinfromfile["skipCOHP"] = False
        assert self.Lobsterinfromfile["skipcohp"] is False

    def test_initialize_from_dict(self):
        # initialize from dict
        lobsterin = Lobsterin(
            {
                "cohpstartenergy": -15.0,
                "cohpendenergy": 5.0,
                "basisset": "pbeVaspFit2015",
                "gaussiansmearingwidth": 0.1,
                "basisfunctions": ["Fe 3d 4p 4s", "Co 3d 4p 4s"],
                "skipdos": True,
                "skipcohp": True,
                "skipcoop": True,
                "skippopulationanalysis": True,
                "skipgrosspopulation": True,
            }
        )
        assert lobsterin["cohpstartenergy"] == approx(-15.0)
        assert lobsterin["cohpendenergy"] == approx(5.0)
        assert lobsterin["basisset"] == "pbeVaspFit2015"
        assert lobsterin["gaussiansmearingwidth"] == approx(0.1)
        assert lobsterin["basisfunctions"][0] == "Fe 3d 4p 4s"
        assert lobsterin["basisfunctions"][1] == "Co 3d 4p 4s"
        assert {*lobsterin} >= {"skipdos", "skipcohp", "skipcoop", "skippopulationanalysis", "skipgrosspopulation"}
        with pytest.raises(IOError, match="There are duplicates for the keywords! The program will stop here."):
            lobsterin2 = Lobsterin({"cohpstartenergy": -15.0, "cohpstartEnergy": -20.0})
        lobsterin2 = Lobsterin({"cohpstartenergy": -15.0})
        # can only calculate nbands if basis functions are provided
        with pytest.raises(IOError, match="No basis functions are provided. The program cannot calculate nbands"):
            lobsterin2._get_nbands(structure=Structure.from_file(f"{TEST_FILES_DIR}/POSCAR.Fe3O4"))

    def test_standard_settings(self):
        # test standard settings
        for option in [
            "standard",
            "standard_from_projection",
            "standard_with_fatband",
            "onlyprojection",
            "onlydos",
            "onlycohp",
            "onlycoop",
            "onlycobi",
            "onlycohpcoop",
            "onlycohpcoopcobi",
        ]:
            lobsterin1 = Lobsterin.standard_calculations_from_vasp_files(
                f"{TEST_FILES_DIR}/POSCAR.Fe3O4",
                f"{TEST_FILES_DIR}/INCAR.lobster",
                f"{TEST_FILES_DIR}/POTCAR.Fe3O4",
                option=option,
            )
            assert lobsterin1["cohpstartenergy"] == approx(-35.0)
            assert lobsterin1["cohpendenergy"] == approx(5.0)
            assert lobsterin1["basisset"] == "pbeVaspFit2015"
            assert lobsterin1["gaussiansmearingwidth"] == approx(0.1)
            assert lobsterin1["basisfunctions"][0] == "Fe 3d 4p 4s "
            assert lobsterin1["basisfunctions"][1] == "O 2p 2s "

            if option in [
                "standard",
                "standard_with_fatband",
                "onlyprojection",
                "onlycohp",
                "onlycoop",
                "onlycohpcoop",
            ]:
                assert lobsterin1["saveProjectiontoFile"]
            if option in [
                "standard",
                "standard_with_fatband",
                "onlycohp",
                "onlycoop",
                "onlycohpcoop",
            ]:
                assert lobsterin1["cohpGenerator"] == "from 0.1 to 6.0 orbitalwise"
            if option in ["standard"]:
                assert "skipdos" not in lobsterin1
                assert "skipcohp" not in lobsterin1
                assert "skipcoop" not in lobsterin1
            if option in ["standard_with_fatband"]:
                assert lobsterin1["createFatband"] == ["Fe 3d 4p 4s ", "O 2p 2s "]
                assert "skipdos" not in lobsterin1
                assert "skipcohp" not in lobsterin1
                assert "skipcoop" not in lobsterin1
            if option in ["standard_from_projection"]:
                assert lobsterin1["loadProjectionFromFile"], True
            if option in [
                "onlyprojection",
                "onlycohp",
                "onlycoop",
                "onlycobi",
                "onlycohpcoop",
                "onlycohpcoopcobi",
            ]:
                assert lobsterin1["skipdos"], True
                assert lobsterin1["skipPopulationAnalysis"], True
                assert lobsterin1["skipGrossPopulation"], True
                assert lobsterin1["skipMadelungEnergy"], True

            if option in ["onlydos"]:
                assert lobsterin1["skipPopulationAnalysis"], True
                assert lobsterin1["skipGrossPopulation"], True
                assert lobsterin1["skipcohp"], True
                assert lobsterin1["skipcoop"], True
                assert lobsterin1["skipcobi"], True
                assert lobsterin1["skipMadelungEnergy"], True
            if option in ["onlycohp"]:
                assert lobsterin1["skipcoop"], True
                assert lobsterin1["skipcobi"], True
            if option in ["onlycoop"]:
                assert lobsterin1["skipcohp"], True
                assert lobsterin1["skipcobi"], True
            if option in ["onlyprojection"]:
                assert lobsterin1["skipdos"], True
            if option in ["onlymadelung"]:
                assert lobsterin1["skipPopulationAnalysis"], True
                assert lobsterin1["skipGrossPopulation"], True
                assert lobsterin1["skipcohp"], True
                assert lobsterin1["skipcoop"], True
                assert lobsterin1["skipcobi"], True
                assert lobsterin1["skipdos"], True
        # test basis functions by dict
        lobsterin_new = Lobsterin.standard_calculations_from_vasp_files(
            f"{TEST_FILES_DIR}/POSCAR.Fe3O4",
            f"{TEST_FILES_DIR}/INCAR.lobster",
            dict_for_basis={"Fe": "3d 4p 4s", "O": "2s 2p"},
            option="standard",
        )
        assert lobsterin_new["basisfunctions"] == ["Fe 3d 4p 4s", "O 2s 2p"]

        # test gaussian smearing
        lobsterin_new = Lobsterin.standard_calculations_from_vasp_files(
            f"{TEST_FILES_DIR}/POSCAR.Fe3O4",
            f"{TEST_FILES_DIR}/INCAR.lobster2",
            dict_for_basis={"Fe": "3d 4p 4s", "O": "2s 2p"},
            option="standard",
        )
        assert "gaussiansmearingwidth" not in lobsterin_new

        # fatband and ISMEAR=-5 does not work together
        with pytest.raises(ValueError, match="ISMEAR has to be 0 for a fatband calculation with Lobster"):
            lobsterin_new = Lobsterin.standard_calculations_from_vasp_files(
                f"{TEST_FILES_DIR}/POSCAR.Fe3O4",
                f"{TEST_FILES_DIR}/INCAR.lobster2",
                dict_for_basis={"Fe": "3d 4p 4s", "O": "2s 2p"},
                option="standard_with_fatband",
            )

    def test_standard_with_energy_range_from_vasprun(self):
        # test standard_with_energy_range_from_vasprun
        lobsterin_comp = Lobsterin.standard_calculations_from_vasp_files(
            f"{TEST_FILES_DIR}/POSCAR.C2.gz",
            f"{TEST_FILES_DIR}/INCAR.C2.gz",
            f"{TEST_FILES_DIR}/POTCAR.C2.gz",
            f"{TEST_FILES_DIR}/vasprun.xml.C2.gz",
            option="standard_with_energy_range_from_vasprun",
        )
        assert lobsterin_comp["COHPstartEnergy"] == -28.3679
        assert lobsterin_comp["COHPendEnergy"] == 32.8968
        assert lobsterin_comp["COHPSteps"] == 301

    def test_diff(self):
        # test diff
        assert self.Lobsterinfromfile.diff(self.Lobsterinfromfile2)["Different"] == {}
        assert self.Lobsterinfromfile.diff(self.Lobsterinfromfile2)["Same"]["COHPSTARTENERGY"] == approx(-15.0)

        # test diff in both directions
        for entry in self.Lobsterinfromfile.diff(self.Lobsterinfromfile3)["Same"]:
            assert entry in self.Lobsterinfromfile3.diff(self.Lobsterinfromfile)["Same"]
        for entry in self.Lobsterinfromfile3.diff(self.Lobsterinfromfile)["Same"]:
            assert entry in self.Lobsterinfromfile.diff(self.Lobsterinfromfile3)["Same"]
        for entry in self.Lobsterinfromfile.diff(self.Lobsterinfromfile3)["Different"]:
            assert entry in self.Lobsterinfromfile3.diff(self.Lobsterinfromfile)["Different"]
        for entry in self.Lobsterinfromfile3.diff(self.Lobsterinfromfile)["Different"]:
            assert entry in self.Lobsterinfromfile.diff(self.Lobsterinfromfile3)["Different"]

        assert (
            self.Lobsterinfromfile.diff(self.Lobsterinfromfile3)["Different"]["SKIPCOHP"]["lobsterin1"]
            == self.Lobsterinfromfile3.diff(self.Lobsterinfromfile)["Different"]["SKIPCOHP"]["lobsterin2"]
        )

    def test_dict_functionality(self):
        assert self.Lobsterinfromfile.get("COHPstartEnergy") == -15.0
        assert self.Lobsterinfromfile.get("COHPstartEnergy") == -15.0
        assert self.Lobsterinfromfile.get("COhPstartenergy") == -15.0
        lobsterincopy = self.Lobsterinfromfile.copy()
        lobsterincopy.update({"cohpstarteNergy": -10.00})
        assert lobsterincopy["cohpstartenergy"] == -10.0
        lobsterincopy.pop("cohpstarteNergy")
        assert "cohpstartenergy" not in lobsterincopy
        lobsterincopy.pop("cohpendenergY")
        lobsterincopy["cohpsteps"] = 100
        assert lobsterincopy["cohpsteps"] == 100
        before = len(lobsterincopy.items())
        lobsterincopy.popitem()
        after = len(lobsterincopy.items())
        assert before != after

    def test_read_write_lobsterin(self):
        outfile_path = tempfile.mkstemp()[1]
        lobsterin1 = Lobsterin.from_file(f"{TEST_FILES_DIR}/cohp/lobsterin.1")
        lobsterin1.write_lobsterin(outfile_path)
        lobsterin2 = Lobsterin.from_file(outfile_path)
        assert lobsterin1.diff(lobsterin2)["Different"] == {}

        # TODO: will integer vs float break cohpsteps?

    def test_get_basis(self):
        # get basis functions
        lobsterin1 = Lobsterin({})
        potcar = Potcar.from_file(f"{TEST_FILES_DIR}/POTCAR.Fe3O4")
        potcar_names = [name["symbol"] for name in potcar.spec]

        assert lobsterin1.get_basis(
            Structure.from_file(f"{TEST_FILES_DIR}/Fe3O4.cif"),
            potcar_symbols=potcar_names,
        ) == ["Fe 3d 4p 4s ", "O 2p 2s "]
        potcar = Potcar.from_file(f"{TEST_FILES_DIR}/cohp/POTCAR.GaAs")
        potcar_names = [name["symbol"] for name in potcar.spec]
        assert lobsterin1.get_basis(
            Structure.from_file(f"{TEST_FILES_DIR}/cohp/POSCAR.GaAs"),
            potcar_symbols=potcar_names,
        ) == ["Ga 3d 4p 4s ", "As 4p 4s "]

    def test_get_all_possible_basis_functions(self):
        potcar = Potcar.from_file(f"{TEST_FILES_DIR}/POTCAR.Fe3O4")
        potcar_names = [name["symbol"] for name in potcar.spec]
        result = Lobsterin.get_all_possible_basis_functions(
            Structure.from_file(f"{TEST_FILES_DIR}/Fe3O4.cif"),
            potcar_symbols=potcar_names,
        )
        assert result[0] == {"Fe": "3d 4s", "O": "2p 2s"}
        assert result[1] == {"Fe": "3d 4s 4p", "O": "2p 2s"}

        potcar2 = Potcar.from_file(f"{TEST_FILES_DIR}/POT_GGA_PAW_PBE_54/POTCAR.Fe.gz")
        Potcar_names2 = [name["symbol"] for name in potcar2.spec]
        result2 = Lobsterin.get_all_possible_basis_functions(
            Structure.from_file(f"{TEST_FILES_DIR}/Fe.cif"),
            potcar_symbols=Potcar_names2,
        )
        assert result2[0] == {"Fe": "3d 4s"}

    def test_get_potcar_symbols(self):
        lobsterin1 = Lobsterin({})
        assert lobsterin1._get_potcar_symbols(f"{TEST_FILES_DIR}/POTCAR.Fe3O4") == ["Fe", "O"]
        assert lobsterin1._get_potcar_symbols(f"{TEST_FILES_DIR}/cohp/POTCAR.GaAs") == ["Ga_d", "As"]

    def test_write_lobsterin(self):
        # write lobsterin, read it and compare it
        outfile_path = tempfile.mkstemp()[1]
        lobsterin1 = Lobsterin.standard_calculations_from_vasp_files(
            f"{TEST_FILES_DIR}/POSCAR.Fe3O4",
            f"{TEST_FILES_DIR}/INCAR.lobster",
            f"{TEST_FILES_DIR}/POTCAR.Fe3O4",
            option="standard",
        )
        lobsterin1.write_lobsterin(outfile_path)
        lobsterin2 = Lobsterin.from_file(outfile_path)
        assert lobsterin1.diff(lobsterin2)["Different"] == {}

    def test_write_incar(self):
        # write INCAR and compare
        outfile_path = tempfile.mkstemp()[1]
        lobsterin1 = Lobsterin.standard_calculations_from_vasp_files(
            f"{TEST_FILES_DIR}/POSCAR.Fe3O4",
            f"{TEST_FILES_DIR}/INCAR.lobster",
            f"{TEST_FILES_DIR}/POTCAR.Fe3O4",
            option="standard",
        )
        lobsterin1.write_INCAR(
            f"{TEST_FILES_DIR}/INCAR.lobster3",
            outfile_path,
            f"{TEST_FILES_DIR}/POSCAR.Fe3O4",
        )

        incar1 = Incar.from_file(f"{TEST_FILES_DIR}/INCAR.lobster3")
        incar2 = Incar.from_file(outfile_path)

        assert incar1.diff(incar2)["Different"] == {
            "ISYM": {"INCAR1": 2, "INCAR2": -1},
            "NBANDS": {"INCAR1": None, "INCAR2": 86},
            "NSW": {"INCAR1": 500, "INCAR2": 0},
            "LWAVE": {"INCAR1": False, "INCAR2": True},
        }

    def test_write_kpoints(self):
        # line mode
        outfile_path = tempfile.mkstemp()[1]
        outfile_path2 = tempfile.mkstemp(prefix="POSCAR")[1]
        lobsterin1 = Lobsterin({})
        # test writing primitive cell
        lobsterin1.write_POSCAR_with_standard_primitive(
            POSCAR_input=f"{TEST_FILES_DIR}/POSCAR.Fe3O4", POSCAR_output=outfile_path2
        )

        lobsterin1.write_KPOINTS(
            POSCAR_input=outfile_path2,
            KPOINTS_output=outfile_path,
            kpoints_line_density=58,
        )
        kpoint = Kpoints.from_file(outfile_path)
        assert kpoint.num_kpts == 562
        assert kpoint.kpts[-1][0] == approx(-0.5)
        assert kpoint.kpts[-1][1] == approx(0.5)
        assert kpoint.kpts[-1][2] == approx(0.5)
        assert kpoint.labels[-1] == "T"
        kpoint2 = Kpoints.from_file(f"{TEST_FILES_DIR}/KPOINTS_band.lobster")

        labels = []
        number = 0
        for label in kpoint.labels:
            if label is not None:
                if number != 0:
                    if label != labels[number - 1]:
                        labels.append(label)
                        number += 1
                else:
                    labels.append(label)
                    number += 1

        labels2 = []
        number2 = 0
        for label in kpoint2.labels:
            if label is not None:
                if number2 != 0:
                    if label != labels2[number2 - 1]:
                        labels2.append(label)
                        number2 += 1
                else:
                    labels2.append(label)
                    number2 += 1
        assert labels == labels2

        # without line mode
        lobsterin1.write_KPOINTS(POSCAR_input=outfile_path2, KPOINTS_output=outfile_path, line_mode=False)
        kpoint = Kpoints.from_file(outfile_path)
        kpoint2 = Kpoints.from_file(f"{TEST_FILES_DIR}/IBZKPT.lobster")

        for num_kpt, list_kpoint in enumerate(kpoint.kpts):
            assert list_kpoint[0] == approx(kpoint2.kpts[num_kpt][0])
            assert list_kpoint[1] == approx(kpoint2.kpts[num_kpt][1])
            assert list_kpoint[2] == approx(kpoint2.kpts[num_kpt][2])

        assert kpoint.num_kpts == 108

        # without line mode, use grid instead of reciprocal density
        lobsterin1.write_KPOINTS(
            POSCAR_input=outfile_path2,
            KPOINTS_output=outfile_path,
            line_mode=False,
            from_grid=True,
            input_grid=[6, 6, 3],
        )
        kpoint = Kpoints.from_file(outfile_path)
        kpoint2 = Kpoints.from_file(f"{TEST_FILES_DIR}/IBZKPT.lobster")

        for num_kpt, list_kpoint in enumerate(kpoint.kpts):
            assert list_kpoint[0] == approx(kpoint2.kpts[num_kpt][0])
            assert list_kpoint[1] == approx(kpoint2.kpts[num_kpt][1])
            assert list_kpoint[2] == approx(kpoint2.kpts[num_kpt][2])

        assert kpoint.num_kpts == 108

        #
        # #without line mode, using a certain grid, isym=0 instead of -1
        lobsterin1.write_KPOINTS(
            POSCAR_input=f"{TEST_FILES_DIR}/cohp/POSCAR.Li",
            KPOINTS_output=outfile_path,
            line_mode=False,
            from_grid=True,
            input_grid=[3, 3, 3],
            isym=0,
        )

        kpoint1 = Kpoints.from_file(outfile_path)
        kpoint2 = Kpoints.from_file(f"{TEST_FILES_DIR}/cohp/IBZKPT_3_3_3_Li")
        for ikpoint, kpoint in enumerate(kpoint1.kpts):
            assert self.is_kpoint_in_list(
                kpoint,
                kpoint2.kpts,
                kpoint1.kpts_weights[ikpoint],
                kpoint2.kpts_weights,
            )
        for ikpoint, kpoint in enumerate(kpoint2.kpts):
            assert self.is_kpoint_in_list(
                kpoint,
                kpoint1.kpts,
                kpoint2.kpts_weights[ikpoint],
                kpoint1.kpts_weights,
            )

        lobsterin1.write_KPOINTS(
            POSCAR_input=f"{TEST_FILES_DIR}/cohp/POSCAR.Li",
            KPOINTS_output=outfile_path,
            line_mode=False,
            from_grid=True,
            input_grid=[2, 2, 2],
            isym=0,
        )

        kpoint1 = Kpoints.from_file(outfile_path)
        kpoint2 = Kpoints.from_file(f"{TEST_FILES_DIR}/cohp/IBZKPT_2_2_2_Li")
        for ikpoint, kpoint in enumerate(kpoint1.kpts):
            assert self.is_kpoint_in_list(
                kpoint,
                kpoint2.kpts,
                kpoint1.kpts_weights[ikpoint],
                kpoint2.kpts_weights,
            )
        for ikpoint, kpoint in enumerate(kpoint2.kpts):
            assert self.is_kpoint_in_list(
                kpoint,
                kpoint1.kpts,
                kpoint2.kpts_weights[ikpoint],
                kpoint1.kpts_weights,
            )

    def is_kpoint_in_list(self, kpoint, kpointlist, weight, weightlist) -> bool:
        found = 0
        for ikpoint2, kpoint2 in enumerate(kpointlist):
            if (
                np.isclose(kpoint[0], kpoint2[0])
                and np.isclose(kpoint[1], kpoint2[1])
                and np.isclose(kpoint[2], kpoint2[2])
            ):
                if weight == weightlist[ikpoint2]:
                    found += 1
            elif (
                np.isclose(-kpoint[0], kpoint2[0])
                and np.isclose(-kpoint[1], kpoint2[1])
                and np.isclose(-kpoint[2], kpoint2[2])
            ) and weight == weightlist[ikpoint2]:
                found += 1
        return found == 1

    def test_msonable_implementation(self):
        # tests as dict and from dict methods
        new_lobsterin = Lobsterin.from_dict(self.Lobsterinfromfile.as_dict())
        assert new_lobsterin == self.Lobsterinfromfile
        new_lobsterin.to_json()


class TestBandoverlaps(unittest.TestCase):
    def setUp(self):
        # test spin polarlized calc and non spinpolarized calc

        self.bandoverlaps1 = Bandoverlaps(f"{TEST_FILES_DIR}/cohp/bandOverlaps.lobster.1")
        self.bandoverlaps2 = Bandoverlaps(f"{TEST_FILES_DIR}/cohp/bandOverlaps.lobster.2")

        self.bandoverlaps1_new = Bandoverlaps(f"{TEST_FILES_DIR}/cohp/bandOverlaps.lobster.new.1")
        self.bandoverlaps2_new = Bandoverlaps(f"{TEST_FILES_DIR}/cohp/bandOverlaps.lobster.new.2")

    def test_attributes(self):
        # bandoverlapsdict
        assert self.bandoverlaps1.bandoverlapsdict[Spin.up]["0.5 0 0"]["maxDeviation"] == approx(0.000278953)
        assert self.bandoverlaps1_new.bandoverlapsdict[Spin.up]["0 0 0"]["maxDeviation"] == approx(0.0640933)
        assert self.bandoverlaps1.bandoverlapsdict[Spin.up]["0.5 0 0"]["matrix"][-1][-1] == approx(0.0188058)
        assert self.bandoverlaps1_new.bandoverlapsdict[Spin.up]["0 0 0"]["matrix"][-1][-1] == approx(1.0)
        assert self.bandoverlaps1.bandoverlapsdict[Spin.up]["0.5 0 0"]["matrix"][0][0] == approx(1)
        assert self.bandoverlaps1_new.bandoverlapsdict[Spin.up]["0 0 0"]["matrix"][0][0] == approx(0.995849)

        assert self.bandoverlaps1.bandoverlapsdict[Spin.down]["0.0261194 0.0261194 0.473881"]["maxDeviation"] == approx(
            4.31567e-05
        )
        assert self.bandoverlaps1_new.bandoverlapsdict[Spin.down]["0 0 0"]["maxDeviation"] == approx(0.064369)
        assert self.bandoverlaps1.bandoverlapsdict[Spin.down]["0.0261194 0.0261194 0.473881"]["matrix"][0][
            -1
        ] == approx(4.0066e-07)
        assert self.bandoverlaps1_new.bandoverlapsdict[Spin.down]["0 0 0"]["matrix"][0][-1] == approx(1.37447e-09)

        # maxDeviation
        assert self.bandoverlaps1.max_deviation[0] == approx(0.000278953)
        assert self.bandoverlaps1_new.max_deviation[0] == approx(0.39824)
        assert self.bandoverlaps1.max_deviation[-1] == approx(4.31567e-05)
        assert self.bandoverlaps1_new.max_deviation[-1] == approx(0.324898)

        assert self.bandoverlaps2.max_deviation[0] == approx(0.000473319)
        assert self.bandoverlaps2_new.max_deviation[0] == approx(0.403249)
        assert self.bandoverlaps2.max_deviation[-1] == approx(1.48451e-05)
        assert self.bandoverlaps2_new.max_deviation[-1] == approx(0.45154)

    def test_has_good_quality(self):
        assert not self.bandoverlaps1.has_good_quality_maxDeviation(limit_maxDeviation=0.1)
        assert not self.bandoverlaps1_new.has_good_quality_maxDeviation(limit_maxDeviation=0.1)
        assert not self.bandoverlaps1.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=9,
            number_occ_bands_spin_down=5,
            limit_deviation=0.1,
            spin_polarized=True,
        )
        assert not self.bandoverlaps1_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=9,
            number_occ_bands_spin_down=5,
            limit_deviation=0.1,
            spin_polarized=True,
        )
        assert self.bandoverlaps1.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=3,
            number_occ_bands_spin_down=0,
            limit_deviation=0.001,
            spin_polarized=True,
        )
        assert self.bandoverlaps1_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=3,
            number_occ_bands_spin_down=0,
            limit_deviation=0.01,
            spin_polarized=True,
        )
        assert not self.bandoverlaps1.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=1,
            number_occ_bands_spin_down=1,
            limit_deviation=0.000001,
            spin_polarized=True,
        )
        assert not self.bandoverlaps1_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=1,
            number_occ_bands_spin_down=1,
            limit_deviation=0.000001,
            spin_polarized=True,
        )
        assert not self.bandoverlaps1.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=1,
            number_occ_bands_spin_down=0,
            limit_deviation=0.000001,
            spin_polarized=True,
        )
        assert not self.bandoverlaps1_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=1,
            number_occ_bands_spin_down=0,
            limit_deviation=0.000001,
            spin_polarized=True,
        )
        assert not self.bandoverlaps1.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=0,
            number_occ_bands_spin_down=1,
            limit_deviation=0.000001,
            spin_polarized=True,
        )
        assert not self.bandoverlaps1_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=0,
            number_occ_bands_spin_down=1,
            limit_deviation=0.000001,
            spin_polarized=True,
        )
        assert not self.bandoverlaps1.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=4,
            number_occ_bands_spin_down=4,
            limit_deviation=0.001,
            spin_polarized=True,
        )
        assert not self.bandoverlaps1_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=4,
            number_occ_bands_spin_down=4,
            limit_deviation=0.001,
            spin_polarized=True,
        )

        assert self.bandoverlaps1.has_good_quality_maxDeviation(limit_maxDeviation=100)
        assert self.bandoverlaps1_new.has_good_quality_maxDeviation(limit_maxDeviation=100)
        assert self.bandoverlaps2.has_good_quality_maxDeviation()
        assert not self.bandoverlaps2_new.has_good_quality_maxDeviation()
        assert not self.bandoverlaps2.has_good_quality_maxDeviation(limit_maxDeviation=0.0000001)
        assert not self.bandoverlaps2_new.has_good_quality_maxDeviation(limit_maxDeviation=0.0000001)
        assert not self.bandoverlaps2.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=10, limit_deviation=0.0000001
        )
        assert not self.bandoverlaps2_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=10, limit_deviation=0.0000001
        )
        assert self.bandoverlaps2.has_good_quality_check_occupied_bands(number_occ_bands_spin_up=1, limit_deviation=0.1)
        assert self.bandoverlaps2_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=1, limit_deviation=0.1
        )

        assert not self.bandoverlaps2.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=1, limit_deviation=1e-8
        )
        assert not self.bandoverlaps2_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=1, limit_deviation=1e-8
        )
        assert self.bandoverlaps2.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=10, limit_deviation=0.1
        )
        assert self.bandoverlaps2_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=2, limit_deviation=0.1
        )

        assert self.bandoverlaps2.has_good_quality_check_occupied_bands(number_occ_bands_spin_up=1, limit_deviation=0.1)
        assert self.bandoverlaps2_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=1, limit_deviation=0.1
        )


class TestGrosspop(unittest.TestCase):
    def setUp(self):
        self.grosspop1 = Grosspop(f"{TEST_FILES_DIR}/cohp/GROSSPOP.lobster")

    def test_attributes(self):
        assert self.grosspop1.list_dict_grosspop[0]["Mulliken GP"]["3s"] == approx(0.52)
        assert self.grosspop1.list_dict_grosspop[0]["Mulliken GP"]["3p_y"] == approx(0.38)
        assert self.grosspop1.list_dict_grosspop[0]["Mulliken GP"]["3p_z"] == approx(0.37)
        assert self.grosspop1.list_dict_grosspop[0]["Mulliken GP"]["3p_x"] == approx(0.37)
        assert self.grosspop1.list_dict_grosspop[0]["Mulliken GP"]["total"] == approx(1.64)
        assert self.grosspop1.list_dict_grosspop[0]["element"] == "Si"
        assert self.grosspop1.list_dict_grosspop[0]["Loewdin GP"]["3s"] == approx(0.61)
        assert self.grosspop1.list_dict_grosspop[0]["Loewdin GP"]["3p_y"] == approx(0.52)
        assert self.grosspop1.list_dict_grosspop[0]["Loewdin GP"]["3p_z"] == approx(0.52)
        assert self.grosspop1.list_dict_grosspop[0]["Loewdin GP"]["3p_x"] == approx(0.52)
        assert self.grosspop1.list_dict_grosspop[0]["Loewdin GP"]["total"] == approx(2.16)
        assert self.grosspop1.list_dict_grosspop[5]["Mulliken GP"]["2s"] == approx(1.80)
        assert self.grosspop1.list_dict_grosspop[5]["Loewdin GP"]["2s"] == approx(1.60)
        assert self.grosspop1.list_dict_grosspop[5]["element"] == "O"
        assert self.grosspop1.list_dict_grosspop[8]["Mulliken GP"]["2s"] == approx(1.80)
        assert self.grosspop1.list_dict_grosspop[8]["Loewdin GP"]["2s"] == approx(1.60)
        assert self.grosspop1.list_dict_grosspop[8]["element"] == "O"

    def test_structure_with_grosspop(self):
        struct_dict = {
            "@module": "pymatgen.core.structure",
            "@class": "Structure",
            "charge": None,
            "lattice": {
                "matrix": [
                    [5.021897888834907, 4.53806e-11, 0.0],
                    [-2.5109484443388332, 4.349090983701526, 0.0],
                    [0.0, 0.0, 5.511929408565514],
                ],
                "a": 5.021897888834907,
                "b": 5.0218974974248045,
                "c": 5.511929408565514,
                "alpha": 90.0,
                "beta": 90.0,
                "gamma": 119.99999598960493,
                "volume": 120.38434608659402,
            },
            "sites": [
                {
                    "species": [{"element": "Si", "occu": 1}],
                    "abc": [-3e-16, 0.4763431475490085, 0.6666669999999968],
                    "xyz": [-1.1960730853096477, 2.0716596881533986, 3.674621443020128],
                    "label": "Si",
                    "properties": {"Total Mulliken GP": 1.64, "Total Loewdin GP": 2.16},
                },
                {
                    "species": [{"element": "Si", "occu": 1}],
                    "abc": [0.5236568524509936, 0.5236568524509926, 0.0],
                    "xyz": [1.3148758827683875, 2.277431295571896, 0.0],
                    "label": "Si",
                    "properties": {"Total Mulliken GP": 1.64, "Total Loewdin GP": 2.16},
                },
                {
                    "species": [{"element": "Si", "occu": 1}],
                    "abc": [0.4763431475490066, -1.2e-15, 0.3333330000000032],
                    "xyz": [
                        2.392146647037334,
                        2.1611518932482004e-11,
                        1.8373079655453863,
                    ],
                    "label": "Si",
                    "properties": {"Total Mulliken GP": 1.64, "Total Loewdin GP": 2.16},
                },
                {
                    "species": [{"element": "O", "occu": 1}],
                    "abc": [0.1589037798059321, 0.7440031622164922, 0.4613477252144715],
                    "xyz": [-1.0701550264153763, 3.235737444648381, 2.5429160941844473],
                    "label": "O",
                    "properties": {"Total Mulliken GP": 7.18, "Total Loewdin GP": 6.92},
                },
                {
                    "species": [{"element": "O", "occu": 1}],
                    "abc": [0.2559968377835071, 0.4149006175894398, 0.7946807252144676],
                    "xyz": [0.2437959189219816, 1.8044405351020447, 4.380224059729795],
                    "label": "O",
                    "properties": {"Total Mulliken GP": 7.18, "Total Loewdin GP": 6.92},
                },
                {
                    "species": [{"element": "O", "occu": 1}],
                    "abc": [0.5850993824105679, 0.8410962201940679, 0.1280147252144683],
                    "xyz": [0.8263601076506712, 3.6580039876980064, 0.7056081286390611],
                    "label": "O",
                    "properties": {"Total Mulliken GP": 7.18, "Total Loewdin GP": 6.92},
                },
                {
                    "species": [{"element": "O", "occu": 1}],
                    "abc": [0.7440031622164928, 0.1589037798059326, 0.5386522747855285],
                    "xyz": [3.337308710918233, 0.6910869960638374, 2.969013314381067],
                    "label": "O",
                    "properties": {"Total Mulliken GP": 7.18, "Total Loewdin GP": 6.92},
                },
                {
                    "species": [{"element": "O", "occu": 1}],
                    "abc": [0.4149006175894392, 0.2559968377835, 0.2053192747855324],
                    "xyz": [1.4407936739605638, 1.1133535390791505, 1.13170534883572],
                    "label": "O",
                    "properties": {"Total Mulliken GP": 7.18, "Total Loewdin GP": 6.92},
                },
                {
                    "species": [{"element": "O", "occu": 1}],
                    "abc": [0.841096220194068, 0.5850993824105675, 0.8719852747855317],
                    "xyz": [2.754744948452184, 2.5446504486493, 4.806321279926453],
                    "label": "O",
                    "properties": {"Total Mulliken GP": 7.18, "Total Loewdin GP": 6.92},
                },
            ],
        }

        new_structure = self.grosspop1.get_structure_with_total_grosspop(f"{TEST_FILES_DIR}/cohp/POSCAR.SiO2")
        assert_allclose(new_structure.frac_coords, Structure.from_dict(struct_dict).frac_coords)


class TestUtils(PymatgenTest):
    def test_get_all_possible_basis_combinations(self):
        # this basis is just for testing (not correct)
        min_basis = ["Li 1s 2s ", "Na 1s 2s", "Si 1s 2s"]
        max_basis = ["Li 1s 2p 2s ", "Na 1s 2p 2s", "Si 1s 2s"]
        combinations_basis = get_all_possible_basis_combinations(min_basis, max_basis)
        assert combinations_basis == [
            ["Li 1s 2s", "Na 1s 2s", "Si 1s 2s"],
            ["Li 1s 2s", "Na 1s 2s 2p", "Si 1s 2s"],
            ["Li 1s 2s 2p", "Na 1s 2s", "Si 1s 2s"],
            ["Li 1s 2s 2p", "Na 1s 2s 2p", "Si 1s 2s"],
        ]

        min_basis = ["Li 1s 2s"]
        max_basis = ["Li 1s 2s 2p 3s"]
        combinations_basis = get_all_possible_basis_combinations(min_basis, max_basis)
        assert combinations_basis == [["Li 1s 2s"], ["Li 1s 2s 2p"], ["Li 1s 2s 3s"], ["Li 1s 2s 2p 3s"]]

        min_basis = ["Li 1s 2s", "Na 1s 2s"]
        max_basis = ["Li 1s 2s 2p 3s", "Na 1s 2s 2p 3s"]
        combinations_basis = get_all_possible_basis_combinations(min_basis, max_basis)
        assert combinations_basis == [
            ["Li 1s 2s", "Na 1s 2s"],
            ["Li 1s 2s", "Na 1s 2s 2p"],
            ["Li 1s 2s", "Na 1s 2s 3s"],
            ["Li 1s 2s", "Na 1s 2s 2p 3s"],
            ["Li 1s 2s 2p", "Na 1s 2s"],
            ["Li 1s 2s 2p", "Na 1s 2s 2p"],
            ["Li 1s 2s 2p", "Na 1s 2s 3s"],
            ["Li 1s 2s 2p", "Na 1s 2s 2p 3s"],
            ["Li 1s 2s 3s", "Na 1s 2s"],
            ["Li 1s 2s 3s", "Na 1s 2s 2p"],
            ["Li 1s 2s 3s", "Na 1s 2s 3s"],
            ["Li 1s 2s 3s", "Na 1s 2s 2p 3s"],
            ["Li 1s 2s 2p 3s", "Na 1s 2s"],
            ["Li 1s 2s 2p 3s", "Na 1s 2s 2p"],
            ["Li 1s 2s 2p 3s", "Na 1s 2s 3s"],
            ["Li 1s 2s 2p 3s", "Na 1s 2s 2p 3s"],
        ]

        min_basis = ["Si 1s 2s 2p", "Na 1s 2s"]
        max_basis = ["Si 1s 2s 2p 3s", "Na 1s 2s 2p 3s"]
        combinations_basis = get_all_possible_basis_combinations(min_basis, max_basis)
        assert combinations_basis == [
            ["Si 1s 2s 2p", "Na 1s 2s"],
            ["Si 1s 2s 2p", "Na 1s 2s 2p"],
            ["Si 1s 2s 2p", "Na 1s 2s 3s"],
            ["Si 1s 2s 2p", "Na 1s 2s 2p 3s"],
            ["Si 1s 2s 2p 3s", "Na 1s 2s"],
            ["Si 1s 2s 2p 3s", "Na 1s 2s 2p"],
            ["Si 1s 2s 2p 3s", "Na 1s 2s 3s"],
            ["Si 1s 2s 2p 3s", "Na 1s 2s 2p 3s"],
        ]


class TestWavefunction(PymatgenTest):
    def test_parse_file(self):
        grid, points, real, imaginary, distance = Wavefunction._parse_file(
            os.path.join(
                TEST_FILES_DIR,
                "cohp",
                "LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz",
            )
        )
        assert_array_equal([41, 41, 41], grid)
        assert points[4][0] == approx(0.0000)
        assert points[4][1] == approx(0.0000)
        assert points[4][2] == approx(0.4000)
        assert real[8] == approx(1.38863e-01)
        assert imaginary[8] == approx(2.89645e-01)
        assert len(imaginary) == 41 * 41 * 41
        assert len(real) == 41 * 41 * 41
        assert len(points) == 41 * 41 * 41
        assert distance[0] == approx(0.0000)

    def test_set_volumetric_data(self):
        wave1 = Wavefunction(
            filename=f"{TEST_FILES_DIR}/cohp/LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz",
            structure=Structure.from_file(f"{TEST_FILES_DIR}/cohp/POSCAR_O.gz"),
        )

        wave1.set_volumetric_data(grid=wave1.grid, structure=wave1.structure)
        assert hasattr(wave1, "volumetricdata_real")
        assert hasattr(wave1, "volumetricdata_imaginary")

    def test_get_volumetricdata_real(self):
        wave1 = Wavefunction(
            filename=os.path.join(
                TEST_FILES_DIR,
                "cohp",
                "LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz",
            ),
            structure=Structure.from_file(f"{TEST_FILES_DIR}/cohp/POSCAR_O.gz"),
        )
        volumetricdata_real = wave1.get_volumetricdata_real()
        assert volumetricdata_real.data["total"][0, 0, 0] == approx(-3.0966)

    def test_get_volumetricdata_imaginary(self):
        wave1 = Wavefunction(
            filename=os.path.join(
                TEST_FILES_DIR,
                "cohp",
                "LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz",
            ),
            structure=Structure.from_file(f"{TEST_FILES_DIR}/cohp/POSCAR_O.gz"),
        )
        volumetricdata_imaginary = wave1.get_volumetricdata_imaginary()
        assert volumetricdata_imaginary.data["total"][0, 0, 0] == approx(-6.45895e00)

    def test_get_volumetricdata_density(self):
        wave1 = Wavefunction(
            filename=os.path.join(TEST_FILES_DIR, "cohp", "LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz"),
            structure=Structure.from_file(f"{TEST_FILES_DIR}/cohp/POSCAR_O.gz"),
        )
        volumetricdata_density = wave1.get_volumetricdata_density()
        assert volumetricdata_density.data["total"][0, 0, 0] == approx((-3.0966 * -3.0966) + (-6.45895 * -6.45895))

    def test_write_file(self):
        wave1 = Wavefunction(
            filename=f"{TEST_FILES_DIR}/cohp/LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz",
            structure=Structure.from_file(f"{TEST_FILES_DIR}/cohp/POSCAR_O.gz"),
        )
        wave1.write_file(filename=f"{self.tmp_path}/wavecar_test.vasp", part="real")
        assert os.path.isfile("wavecar_test.vasp")

        wave1.write_file(filename=f"{self.tmp_path}/wavecar_test.vasp", part="imaginary")
        assert os.path.isfile("wavecar_test.vasp")
        wave1.write_file(filename=f"{self.tmp_path}/density.vasp", part="density")
        assert os.path.isfile("density.vasp")


class TestSitePotentials(PymatgenTest):
    def setUp(self) -> None:
        self.sitepotential = SitePotential(filename=f"{TEST_FILES_DIR}/cohp/SitePotentials.lobster.perovskite")

    def test_attributes(self):
        assert self.sitepotential.sitepotentials_Loewdin == [-8.77, -17.08, 9.57, 9.57, 8.45]
        assert self.sitepotential.sitepotentials_Mulliken == [-11.38, -19.62, 11.18, 11.18, 10.09]
        assert self.sitepotential.madelungenergies_Loewdin == approx(-28.64)
        assert self.sitepotential.madelungenergies_Mulliken == approx(-40.02)
        assert self.sitepotential.atomlist == ["La1", "Ta2", "N3", "N4", "O5"]
        assert self.sitepotential.types == ["La", "Ta", "N", "N", "O"]
        assert self.sitepotential.num_atoms == 5
        assert self.sitepotential.ewald_splitting == approx(3.14)

    def test_get_structure(self):
        structure = self.sitepotential.get_structure_with_site_potentials(f"{TEST_FILES_DIR}/cohp/POSCAR.perovskite")
        assert structure.site_properties["Loewdin Site Potentials (eV)"] == [-8.77, -17.08, 9.57, 9.57, 8.45]
        assert structure.site_properties["Mulliken Site Potentials (eV)"] == [-11.38, -19.62, 11.18, 11.18, 10.09]


class TestMadelungEnergies(PymatgenTest):
    def setUp(self) -> None:
        self.madelungenergies = MadelungEnergies(filename=f"{TEST_FILES_DIR}/cohp/MadelungEnergies.lobster.perovskite")

    def test_attributes(self):
        assert self.madelungenergies.madelungenergies_Loewdin == approx(-28.64)
        assert self.madelungenergies.madelungenergies_Mulliken == approx(-40.02)
        assert self.madelungenergies.ewald_splitting == approx(3.14)


class TestLobsterMatrices(PymatgenTest):
    def setUp(self) -> None:
        self.hamilton_matrices = LobsterMatrices(
            filename=f"{TEST_FILES_DIR}/cohp/Na_hamiltonMatrices.lobster.gz", e_fermi=-2.79650354
        )
        self.transfer_matrices = LobsterMatrices(filename=f"{TEST_FILES_DIR}/cohp/C_transferMatrices.lobster.gz")
        self.overlap_matrices = LobsterMatrices(filename=f"{TEST_FILES_DIR}/cohp/Si_overlapMatrices.lobster.gz")
        self.coeff_matrices = LobsterMatrices(filename=f"{TEST_FILES_DIR}/cohp/Si_coefficientMatricesLSO1.lobster.gz")

    def test_attributes(self):
        # hamilton matrices
        assert self.hamilton_matrices.average_onsite_energies == pytest.approx(
            {"Na1_3s": 0.58855353, "Na1_2p_y": -25.72719646, "Na1_2p_z": -25.72719646, "Na1_2p_x": -25.72719646}
        )
        ref_onsite_energies = [
            [-0.22519646, -25.76989646, -25.76989646, -25.76989646],
            [1.40230354, -25.68449646, -25.68449646, -25.68449646],
        ]
        assert_allclose(self.hamilton_matrices.onsite_energies, ref_onsite_energies)

        ref_imag_mat_spin_up = np.zeros((4, 4))

        assert_allclose(self.hamilton_matrices.hamilton_matrices["1"][Spin.up].imag, ref_imag_mat_spin_up)

        ref_real_mat_spin_up = [
            [-3.0217, 0.0, 0.0, 0.0],
            [0.0, -28.5664, 0.0, 0.0],
            [0.0, 0.0, -28.5664, 0.0],
            [0.0, 0.0, 0.0, -28.5664],
        ]
        assert_allclose(self.hamilton_matrices.hamilton_matrices["1"][Spin.up].real, ref_real_mat_spin_up)

        # overlap matrices
        assert self.overlap_matrices.average_onsite_overlaps == pytest.approx(
            {"Si1_3s": 1.00000009, "Si1_3p_y": 0.99999995, "Si1_3p_z": 0.99999995, "Si1_3p_x": 0.99999995}
        )
        ref_onsite_ovelaps = [[1.00000009, 0.99999995, 0.99999995, 0.99999995]]

        assert_allclose(self.overlap_matrices.onsite_overlaps, ref_onsite_ovelaps)

        ref_imag_mat = np.zeros((4, 4))

        assert_allclose(self.overlap_matrices.overlap_matrices["1"].imag, ref_imag_mat)

        ref_real_mat = [
            [1.00000009, 0.0, 0.0, 0.0],
            [0.0, 0.99999995, 0.0, 0.0],
            [0.0, 0.0, 0.99999995, 0.0],
            [0.0, 0.0, 0.0, 0.99999995],
        ]

        assert_allclose(self.overlap_matrices.overlap_matrices["1"].real, ref_real_mat)

        assert len(self.overlap_matrices.overlap_matrices) == 1
        # transfer matrices
        ref_onsite_transfer = [
            [-0.70523233, -0.07099237, -0.65987499, -0.07090411],
            [-0.03735031, -0.66865552, 0.69253776, 0.80648063],
        ]
        assert_allclose(self.transfer_matrices.onsite_transfer, ref_onsite_transfer)

        ref_imag_mat_spin_down = [
            [-0.99920553, 0.0, 0.0, 0.0],
            [0.0, 0.71219607, -0.06090336, -0.08690835],
            [0.0, -0.04539545, -0.69302453, 0.08323944],
            [0.0, -0.12220894, -0.09749622, -0.53739499],
        ]

        assert_allclose(self.transfer_matrices.transfer_matrices["1"][Spin.down].imag, ref_imag_mat_spin_down)

        ref_real_mat_spin_down = [
            [-0.03735031, 0.0, 0.0, 0.0],
            [0.0, -0.66865552, 0.06086057, 0.13042529],
            [-0.0, 0.04262018, 0.69253776, -0.12491928],
            [0.0, 0.11473763, 0.09742773, 0.80648063],
        ]

        assert_allclose(self.transfer_matrices.transfer_matrices["1"][Spin.down].real, ref_real_mat_spin_down)

        # coefficient matrices
        assert list(self.coeff_matrices.coefficient_matrices["1"]) == [Spin.up, Spin.down]
        assert self.coeff_matrices.average_onsite_coefficient == pytest.approx(
            {
                "Si1_3s": 0.6232626450000001,
                "Si1_3p_y": -0.029367565000000012,
                "Si1_3p_z": -0.50003867,
                "Si1_3p_x": 0.13529422,
            }
        )

        ref_imag_mat_spin_up = [
            [-0.59697342, 0.0, 0.0, 0.0],
            [0.0, 0.50603774, 0.50538255, -0.26664607],
            [0.0, -0.45269894, 0.56996771, 0.23223275],
            [0.0, 0.47836456, 0.00476861, 0.50184424],
        ]

        assert_allclose(self.coeff_matrices.coefficient_matrices["1"][Spin.up].imag, ref_imag_mat_spin_up)

        ref_real_mat_spin_up = [
            [0.80226096, 0.0, 0.0, 0.0],
            [0.0, -0.33931137, -0.42979933, -0.34286226],
            [0.0, 0.30354633, -0.48472536, 0.29861248],
            [0.0, -0.32075579, -0.00405544, 0.64528776],
        ]

        assert_allclose(self.coeff_matrices.coefficient_matrices["1"][Spin.up].real, ref_real_mat_spin_up)

    def test_raises(self):
        with pytest.raises(ValueError, match="Please provide the fermi energy in eV"):
            self.hamilton_matrices = LobsterMatrices(filename=f"{TEST_FILES_DIR}/cohp/Na_hamiltonMatrices.lobster.gz")

        with pytest.raises(
            OSError,
            match=r"Please check provided input file, it seems to be empty",
        ):
            self.hamilton_matrices = LobsterMatrices(filename=f"{TEST_FILES_DIR}/cohp/hamiltonMatrices.lobster")
