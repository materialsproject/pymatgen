from __future__ import annotations

import copy
import gzip
import os

import numpy as np
import orjson
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from pytest import approx

from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.cohp import IcohpCollection
from pymatgen.electronic_structure.core import Orbital, Spin
from pymatgen.io.lobster import (
    Bandoverlaps,
    Bwdf,
    Charge,
    Cohpcar,
    Doscar,
    Fatband,
    Grosspop,
    Icohplist,
    LobsterMatrices,
    Lobsterout,
    MadelungEnergies,
    NciCobiList,
    Polarization,
    SitePotential,
    Wavefunction,
)
from pymatgen.io.lobster.outputs import _get_lines
from pymatgen.io.vasp import Vasprun
from pymatgen.util.testing import TEST_FILES_DIR, VASP_IN_DIR, VASP_OUT_DIR, MatSciTest

TEST_DIR = f"{TEST_FILES_DIR}/electronic_structure/cohp"

__author__ = "Janine George, Marco Esters"
__copyright__ = "Copyright 2017, The Materials Project"
__version__ = "0.2"
__email__ = "janine.george@uclouvain.be, esters@uoregon.edu"
__date__ = "Dec 10, 2017"


class TestBwdf(MatSciTest):
    def setup_method(self):
        self.bwdf_coop = Bwdf(filename=f"{TEST_DIR}/BWDF.lobster.AlN.gz")
        self.bwdf_cohp = Bwdf(filename=f"{TEST_DIR}/BWDFCOHP.lobster.NaCl.gz")

    def test_attributes(self):
        assert self.bwdf_coop.bin_width == approx(0.02005000, abs=1e-4)
        assert self.bwdf_cohp.bin_width == approx(0.02005000, abs=1e-4)
        assert len(self.bwdf_coop.centers) == 201
        assert len(self.bwdf_cohp.centers) == 143
        assert self.bwdf_coop.bwdf[Spin.down][-2] == approx(0.00082, abs=1e-4)
        assert self.bwdf_coop.bwdf[Spin.up][0] == approx(0.81161, abs=1e-4)
        assert self.bwdf_cohp.bwdf[Spin.up][103] == approx(-0.01392, abs=1e-4)


class TestCohpcar(MatSciTest):
    def setup_method(self):
        self.cohp_bise = Cohpcar(filename=f"{TEST_DIR}/COHPCAR.lobster.BiSe.gz")
        self.coop_bise = Cohpcar(
            filename=f"{TEST_DIR}/COOPCAR.lobster.BiSe.gz",
            are_coops=True,
        )

        # Make sure Cohpcar also works with terminating line ending char
        gz_path = f"{TEST_DIR}/COOPCAR.lobster.gz"
        with gzip.open(gz_path, "rt", encoding="utf-8") as f:
            content = f.read() + "\n"

        # Test default filename (None should be redirected to "COHPCAR.lobster")
        with open("COHPCAR.lobster", "w", encoding="utf-8") as f:
            f.write(content)

        self.cohp_fe = Cohpcar(filename=None)
        self.coop_fe = Cohpcar(
            filename=f"{TEST_DIR}/COOPCAR.lobster.gz",
            are_coops=True,
        )
        self.orb = Cohpcar(filename=f"{TEST_DIR}/COHPCAR.lobster.orbitalwise.gz")
        self.orb_notot = Cohpcar(filename=f"{TEST_DIR}/COHPCAR.lobster.notot.orbitalwise.gz")

        # Lobster 3.1 (Test data is from prerelease of Lobster 3.1)
        self.cohp_KF = Cohpcar(filename=f"{TEST_DIR}/COHPCAR.lobster.KF.gz")
        self.coop_KF = Cohpcar(
            filename=f"{TEST_DIR}/COHPCAR.lobster.KF.gz",
            are_coops=True,
        )

        # example with f electrons
        self.cohp_Na2UO4 = Cohpcar(filename=f"{TEST_DIR}/COHPCAR.lobster.Na2UO4.gz")
        self.coop_Na2UO4 = Cohpcar(
            filename=f"{TEST_DIR}/COOPCAR.lobster.Na2UO4.gz",
            are_coops=True,
        )
        self.cobi = Cohpcar(
            filename=f"{TEST_DIR}/COBICAR.lobster.gz",
            are_cobis=True,
        )
        # 3 center
        self.cobi2 = Cohpcar(
            filename=f"{TEST_DIR}/COBICAR.lobster.GeTe",
            are_cobis=False,
            are_multi_center_cobis=True,
        )
        # 4 center
        self.cobi3 = Cohpcar(
            filename=f"{TEST_DIR}/COBICAR.lobster.GeTe_4center",
            are_cobis=False,
            are_multi_center_cobis=True,
        )
        # partially orbital-resolved
        self.cobi4 = Cohpcar(
            filename=f"{TEST_DIR}/COBICAR.lobster.GeTe.multi.orbitalwise",
            are_cobis=False,
            are_multi_center_cobis=True,
        )
        # fully orbital-resolved
        self.cobi5 = Cohpcar(
            filename=f"{TEST_DIR}/COBICAR.lobster.GeTe.multi.orbitalwise.full",
            are_cobis=False,
            are_multi_center_cobis=True,
        )
        # spin polarized
        # fully orbital-resolved
        self.cobi6 = Cohpcar(
            filename=f"{TEST_DIR}/COBICAR.lobster.B2H6.spin",
            are_cobis=False,
            are_multi_center_cobis=True,
        )

        # COHPCAR.LCFO.lobster from v5.1.1
        self.cohp_nacl_lcfo = Cohpcar(filename=f"{TEST_DIR}/COHPCAR.LCFO.lobster.NaCl.gz", is_lcfo=True)

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

        # test multi-center cobis
        assert not self.cobi2.are_cobis
        assert not self.cobi2.are_coops
        assert self.cobi2.are_multi_center_cobis

        # test COHPCAR.LCFO.lobster v 5.1.1
        assert len(self.cohp_nacl_lcfo.orb_res_cohp) == 16
        assert self.cohp_nacl_lcfo.is_lcfo
        assert not self.cohp_nacl_lcfo.is_spin_polarized
        assert not self.cohp_nacl_lcfo.are_coops
        assert not self.cohp_nacl_lcfo.are_cobis
        assert len(self.cohp_nacl_lcfo.energies) == 11

    def test_energies(self):
        efermi_bise = 5.90043
        elim_bise = (-0.124679, 11.9255)
        efermi_fe = 9.75576
        elim_fe = (-0.277681, 14.7725)
        efermi_KF = -2.87475
        elim_KF = (-11.25000 + efermi_KF, 7.5000 + efermi_KF)

        assert self.cohp_bise.efermi == approx(efermi_bise)
        assert self.coop_bise.efermi == approx(efermi_bise)
        assert self.cohp_fe.efermi == approx(efermi_fe)
        assert self.coop_fe.efermi == approx(efermi_fe)
        # Lobster 3.1
        assert self.cohp_KF.efermi == approx(efermi_KF)
        assert self.coop_KF.efermi == approx(efermi_KF)

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

        for data in [self.cobi2.cohp_data]:
            for bond, val in data.items():
                if bond != "average":
                    if int(bond) >= 13:
                        assert len(val["COHP"][Spin.up]) == 11
                        assert len(val["cells"]) == 3
                    else:
                        assert len(val["COHP"][Spin.up]) == 11
                        assert len(val["cells"]) == 2

        for data in [self.cobi3.cohp_data, self.cobi4.cohp_data]:
            for bond, val in data.items():
                if bond != "average":
                    if int(bond) >= 13:
                        assert len(val["cells"]) == 4
                    else:
                        assert len(val["cells"]) == 2
        for data in [self.cobi5.cohp_data]:
            for bond, val in data.items():
                if bond != "average":
                    if int(bond) >= 25:
                        assert len(val["cells"]) == 4
                    else:
                        assert len(val["cells"]) == 2
        for data in [self.cobi6.cohp_data]:
            for bond, val in data.items():
                if bond != "average":
                    if int(bond) >= 21:
                        assert len(val["cells"]) == 3
                        assert len(val["COHP"][Spin.up]) == 12
                        assert len(val["COHP"][Spin.down]) == 12
                        for cohp1, cohp2 in zip(val["COHP"][Spin.up], val["COHP"][Spin.down], strict=False):
                            assert cohp1 == approx(cohp2, abs=1e-4)
                    else:
                        assert len(val["cells"]) == 2
                        assert len(val["COHP"][Spin.up]) == 12
                        assert len(val["COHP"][Spin.down]) == 12
                        for cohp1, cohp2 in zip(val["COHP"][Spin.up], val["COHP"][Spin.down], strict=False):
                            assert cohp1 == approx(cohp2, abs=1e-3)

    def test_orbital_resolved_cohp(self):
        orbitals = [(Orbital(jj), Orbital(ii)) for ii in range(4) for jj in range(4)]
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

        assert "5s-4s-5s-4s" in self.cobi4.orb_res_cohp["13"]
        assert "5px-4px-5px-4px" in self.cobi4.orb_res_cohp["13"]
        assert len(self.cobi4.orb_res_cohp["13"]["5px-4px-5px-4px"]["COHP"][Spin.up]) == 11

        assert "5s-4s-5s-4s" in self.cobi5.orb_res_cohp["25"]
        assert "5px-4px-5px-4px" in self.cobi5.orb_res_cohp["25"]
        assert len(self.cobi5.orb_res_cohp["25"]["5px-4px-5px-4px"]["COHP"][Spin.up]) == 11

        assert len(self.cobi6.orb_res_cohp["21"]["2py-1s-2s"]["COHP"][Spin.up]) == 12
        assert len(self.cobi6.orb_res_cohp["21"]["2py-1s-2s"]["COHP"][Spin.down]) == 12


class TestDoscar:
    def setup_method(self):
        # first for spin polarized version
        doscar = f"{VASP_OUT_DIR}/DOSCAR.lobster.spin"
        poscar = f"{VASP_IN_DIR}/POSCAR.lobster.spin_DOS"

        # not spin polarized
        doscar2 = f"{VASP_OUT_DIR}/DOSCAR.lobster.nonspin"
        poscar2 = f"{VASP_IN_DIR}/POSCAR.lobster.nonspin_DOS"

        # DOSCAR.LCFO.lobster
        doscar3 = f"{VASP_OUT_DIR}/DOSCAR.LCFO.lobster.AlN"
        poscar3 = f"{VASP_IN_DIR}/POSCAR.AlN"

        self.DOSCAR_spin_pol = Doscar(doscar=doscar, structure_file=poscar)
        self.DOSCAR_nonspin_pol = Doscar(doscar=doscar2, structure_file=poscar2)

        self.DOSCAR_spin_pol = Doscar(doscar=doscar, structure_file=poscar)
        self.DOSCAR_nonspin_pol = Doscar(doscar=doscar2, structure_file=poscar2)

        self.DOSCAR_lcfo = Doscar(doscar=doscar3, structure_file=poscar3, is_lcfo=True)

        with open(f"{TEST_FILES_DIR}/electronic_structure/dos/structure_KF.json", "rb") as file:
            data = orjson.loads(file.read())

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

        assert_allclose(energies_spin, self.DOSCAR_spin_pol.completedos.energies)
        assert_allclose(tdos_up, self.DOSCAR_spin_pol.completedos.densities[Spin.up])
        assert_allclose(tdos_down, self.DOSCAR_spin_pol.completedos.densities[Spin.down])
        assert fermi == approx(self.DOSCAR_spin_pol.completedos.efermi)

        assert_allclose(
            self.DOSCAR_spin_pol.completedos.structure.frac_coords,
            self.structure.frac_coords,
        )
        assert_allclose(
            self.DOSCAR_spin_pol2.completedos.structure.frac_coords,
            self.structure.frac_coords,
        )
        assert_allclose(self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2s"][Spin.up], pdos_f_2s_up)
        assert_allclose(self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2s"][Spin.down], pdos_f_2s_down)
        assert_allclose(self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2p_y"][Spin.up], pdos_f_2py_up)
        assert_allclose(self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2p_y"][Spin.down], pdos_f_2py_down)
        assert_allclose(self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2p_z"][Spin.up], pdos_f_2pz_up)
        assert_allclose(self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2p_z"][Spin.down], pdos_f_2pz_down)
        assert_allclose(self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2p_x"][Spin.up], pdos_f_2px_up)
        assert_allclose(self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2p_x"][Spin.down], pdos_f_2px_down)

        energies_nonspin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        tdos_nonspin = [0.00000, 1.60000, 0.00000, 1.60000, 0.00000, 0.02418]
        pdos_f_2s = [0.00000, 0.00320, 0.00000, 0.00017, 0.00000, 0.00060]
        pdos_f_2py = [0.00000, 0.00322, 0.00000, 0.51635, 0.00000, 0.00037]
        pdos_f_2pz = [0.00000, 0.00322, 0.00000, 0.51636, 0.00000, 0.00037]
        pdos_f_2px = [0.00000, 0.00322, 0.00000, 0.51634, 0.00000, 0.00037]

        assert_allclose(energies_nonspin, self.DOSCAR_nonspin_pol.completedos.energies)

        assert_allclose(tdos_nonspin, self.DOSCAR_nonspin_pol.completedos.densities[Spin.up])

        assert fermi == approx(self.DOSCAR_nonspin_pol.completedos.efermi)

        assert self.DOSCAR_nonspin_pol.completedos.structure == self.structure

        assert_allclose(self.DOSCAR_nonspin_pol.completedos.pdos[self.structure[0]]["2s"][Spin.up], pdos_f_2s)
        assert_allclose(self.DOSCAR_nonspin_pol.completedos.pdos[self.structure[0]]["2p_y"][Spin.up], pdos_f_2py)
        assert_allclose(self.DOSCAR_nonspin_pol.completedos.pdos[self.structure[0]]["2p_z"][Spin.up], pdos_f_2pz)
        assert_allclose(self.DOSCAR_nonspin_pol.completedos.pdos[self.structure[0]]["2p_x"][Spin.up], pdos_f_2px)

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

        assert_allclose(self.DOSCAR_spin_pol.pdos[0]["2s"][Spin.up], pdos_f_2s_up)
        assert_allclose(self.DOSCAR_spin_pol.pdos[0]["2s"][Spin.down], pdos_f_2s_down)
        assert_allclose(self.DOSCAR_spin_pol.pdos[0]["2p_y"][Spin.up], pdos_f_2py_up)
        assert_allclose(self.DOSCAR_spin_pol.pdos[0]["2p_y"][Spin.down], pdos_f_2py_down)
        assert_allclose(self.DOSCAR_spin_pol.pdos[0]["2p_z"][Spin.up], pdos_f_2pz_up)
        assert_allclose(self.DOSCAR_spin_pol.pdos[0]["2p_z"][Spin.down], pdos_f_2pz_down)
        assert_allclose(self.DOSCAR_spin_pol.pdos[0]["2p_x"][Spin.up], pdos_f_2px_up)
        assert_allclose(self.DOSCAR_spin_pol.pdos[0]["2p_x"][Spin.down], pdos_f_2px_down)

        # non spin
        pdos_f_2s = [0.00000, 0.00320, 0.00000, 0.00017, 0.00000, 0.00060]
        pdos_f_2py = [0.00000, 0.00322, 0.00000, 0.51635, 0.00000, 0.00037]
        pdos_f_2pz = [0.00000, 0.00322, 0.00000, 0.51636, 0.00000, 0.00037]
        pdos_f_2px = [0.00000, 0.00322, 0.00000, 0.51634, 0.00000, 0.00037]

        assert_allclose(self.DOSCAR_nonspin_pol.pdos[0]["2s"][Spin.up], pdos_f_2s)
        assert_allclose(self.DOSCAR_nonspin_pol.pdos[0]["2p_y"][Spin.up], pdos_f_2py)
        assert_allclose(self.DOSCAR_nonspin_pol.pdos[0]["2p_z"][Spin.up], pdos_f_2pz)
        assert_allclose(self.DOSCAR_nonspin_pol.pdos[0]["2p_x"][Spin.up], pdos_f_2px)

        # test with DOSCAR.LCFO.lobster file
        pdos_1a1_AlN = [
            0.0,
            0.22594,
            0.01335,
            0.0,
            0.0,
            0.0,
            0.00228,
            0.02836,
            0.03053,
            0.01612,
            0.02379,
        ]
        pdos_3py_Al = [
            0.0,
            0.02794,
            0.00069,
            0.0,
            0.0,
            0.0,
            0.00216,
            0.0682,
            0.06966,
            0.04402,
            0.16579,
        ]
        pdos_2s_N = [
            0.0,
            0.25324,
            0.0157,
            0.0,
            0.0,
            0.0,
            0.0006,
            0.01747,
            0.02247,
            0.01589,
            0.03565,
        ]

        assert self.DOSCAR_lcfo._is_lcfo
        assert_allclose(self.DOSCAR_lcfo.pdos[0]["1a1"][Spin.down], pdos_1a1_AlN)
        assert_allclose(self.DOSCAR_lcfo.pdos[1]["3p_y"][Spin.down], pdos_3py_Al)
        assert_allclose(self.DOSCAR_lcfo.pdos[2]["2s"][Spin.down], pdos_2s_N)

    def test_tdos(self):
        # first for spin polarized version
        energies_spin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        tdos_up = [0.00000, 0.79999, 0.00000, 0.79999, 0.00000, 0.02577]
        tdos_down = [0.00000, 0.79999, 0.00000, 0.79999, 0.00000, 0.02586]
        fermi = 0.0

        assert_allclose(energies_spin, self.DOSCAR_spin_pol.tdos.energies)
        assert_allclose(tdos_up, self.DOSCAR_spin_pol.tdos.densities[Spin.up])
        assert_allclose(tdos_down, self.DOSCAR_spin_pol.tdos.densities[Spin.down])
        assert fermi == approx(self.DOSCAR_spin_pol.tdos.efermi)

        energies_nonspin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        tdos_nonspin = [0.00000, 1.60000, 0.00000, 1.60000, 0.00000, 0.02418]
        fermi = 0.0

        assert_allclose(energies_nonspin, self.DOSCAR_nonspin_pol.tdos.energies)
        assert_allclose(tdos_nonspin, self.DOSCAR_nonspin_pol.tdos.densities[Spin.up])
        assert fermi == approx(self.DOSCAR_nonspin_pol.tdos.efermi)

    def test_energies(self):
        # first for spin polarized version
        energies_spin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]

        assert_allclose(energies_spin, self.DOSCAR_spin_pol.energies)

        energies_nonspin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        assert_allclose(energies_nonspin, self.DOSCAR_nonspin_pol.energies)

    def test_tdensities(self):
        # first for spin polarized version
        tdos_up = [0.00000, 0.79999, 0.00000, 0.79999, 0.00000, 0.02577]
        tdos_down = [0.00000, 0.79999, 0.00000, 0.79999, 0.00000, 0.02586]

        assert_allclose(tdos_up, self.DOSCAR_spin_pol.tdensities[Spin.up])
        assert_allclose(tdos_down, self.DOSCAR_spin_pol.tdensities[Spin.down])

        tdos_nonspin = [0.00000, 1.60000, 0.00000, 1.60000, 0.00000, 0.02418]
        assert_allclose(tdos_nonspin, self.DOSCAR_nonspin_pol.tdensities[Spin.up])

        # test with DOSCAR.LCFO.lobster file
        tdos_up = [
            0.0,
            1.75477,
            0.11803,
            0.0,
            0.0,
            0.0,
            0.04156,
            0.82291,
            0.74449,
            0.42481,
            1.04535,
        ]

        assert_allclose(tdos_up, self.DOSCAR_lcfo.tdensities[Spin.up])

    def test_itdensities(self):
        itdos_up = [1.99997, 4.99992, 4.99992, 7.99987, 7.99987, 8.09650]
        itdos_down = [1.99997, 4.99992, 4.99992, 7.99987, 7.99987, 8.09685]
        assert_allclose(itdos_up, self.DOSCAR_spin_pol.itdensities[Spin.up])
        assert_allclose(itdos_down, self.DOSCAR_spin_pol.itdensities[Spin.down])

        itdos_nonspin = [4.00000, 10.00000, 10.00000, 16.00000, 16.00000, 16.09067]
        assert_allclose(itdos_nonspin, self.DOSCAR_nonspin_pol.itdensities[Spin.up])

    def test_is_spin_polarized(self):
        # first for spin polarized version
        assert self.DOSCAR_spin_pol.is_spin_polarized

        assert not self.DOSCAR_nonspin_pol.is_spin_polarized


class TestCharge(MatSciTest):
    def setup_method(self):
        self.charge2 = Charge(filename=f"{TEST_DIR}/CHARGE.lobster.MnO")
        # gzipped file
        self.charge = Charge(filename=f"{TEST_DIR}/CHARGE.lobster.MnO2.gz")
        self.charge_lcfo = Charge(filename=f"{TEST_DIR}/CHARGE.LCFO.lobster.ALN.gz", is_lcfo=True)

    def test_attributes(self):
        assert self.charge2.mulliken == approx([-1.30, 1.30])
        assert self.charge2.loewdin == approx([-1.25, 1.25])
        assert self.charge2.atomlist == ["O1", "Mn2"]
        assert self.charge2.types == ["O", "Mn"]
        assert self.charge2.num_atoms == 2

        # test with CHARG.LCFO.lobster file
        assert self.charge_lcfo.is_lcfo
        assert self.charge_lcfo.num_atoms == 3
        assert self.charge_lcfo.types == ["AlN", "Al", "N"]
        assert self.charge_lcfo.atomlist == ["AlN1", "Al2", "N3"]
        assert_allclose(self.charge_lcfo.loewdin, [0.0, 1.02, -1.02])
        assert not self.charge_lcfo.mulliken

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
        assert s2 == self.charge2.get_structure_with_charges(f"{VASP_IN_DIR}/POSCAR_MnO")

    def test_exception(self):
        structure_file = f"{VASP_IN_DIR}/POSCAR.AlN"
        with pytest.raises(ValueError, match="CHARGE.LCFO.lobster charges are not sorted site wise"):
            self.charge_lcfo.get_structure_with_charges(structure_filename=structure_file)

    def test_msonable(self):
        dict_data = self.charge2.as_dict()
        charge_from_dict = Charge.from_dict(dict_data)
        all_attributes = vars(self.charge2)
        for attr_name, attr_value in all_attributes.items():
            assert getattr(charge_from_dict, attr_name) == attr_value


class TestLobsterout(MatSciTest):
    def setup_method(self):
        self.lobsterout_normal = Lobsterout(filename=f"{TEST_DIR}/lobsterout.normal")
        # make sure .gz files are also read correctly
        self.lobsterout_normal = Lobsterout(filename=f"{TEST_DIR}/lobsterout.normal2.gz")
        self.lobsterout_fatband_grosspop_densityofenergies = Lobsterout(
            filename=f"{TEST_DIR}/lobsterout.fatband_grosspop_densityofenergy"
        )
        self.lobsterout_saveprojection = Lobsterout(filename=f"{TEST_DIR}/lobsterout.saveprojection")
        self.lobsterout_skipping_all = Lobsterout(filename=f"{TEST_DIR}/lobsterout.skipping_all")
        self.lobsterout_twospins = Lobsterout(filename=f"{TEST_DIR}/lobsterout.twospins")
        self.lobsterout_GaAs = Lobsterout(filename=f"{TEST_DIR}/lobsterout.GaAs")
        self.lobsterout_from_projection = Lobsterout(filename=f"{TEST_DIR}/lobsterout_from_projection")
        self.lobsterout_onethread = Lobsterout(filename=f"{TEST_DIR}/lobsterout.onethread")
        self.lobsterout_cobi_madelung = Lobsterout(filename=f"{TEST_DIR}/lobsterout_cobi_madelung")
        self.lobsterout_doscar_lso = Lobsterout(filename=f"{TEST_DIR}/lobsterout_doscar_lso")

        # TODO: implement skipping madelung/cobi
        self.lobsterout_skipping_cobi_madelung = Lobsterout(filename=f"{TEST_DIR}/lobsterout.skip_cobi_madelung")

        # Lobster v5.1.1
        self.lobsterout_v511 = Lobsterout(filename=f"{TEST_DIR}/lobsterout_v511.gz")

    def test_attributes(self):
        assert self.lobsterout_normal.basis_functions == [
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
        assert self.lobsterout_normal.basis_type == ["pbeVaspFit2015"]
        assert_allclose(self.lobsterout_normal.charge_spilling, [0.0268])
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
        assert self.lobsterout_fatband_grosspop_densityofenergies.basis_type == ["pbeVaspFit2015"]
        assert_allclose(self.lobsterout_fatband_grosspop_densityofenergies.charge_spilling, [0.0268])
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
        assert_allclose(self.lobsterout_saveprojection.charge_spilling, [0.0268])
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
        assert_allclose(self.lobsterout_skipping_all.charge_spilling, [0.0268])
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
        assert self.lobsterout_GaAs.total_spilling[0] == approx(0.0859)

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
        ref_data = {
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
            "has_doscar_lcfo": False,
            "has_cohpcar": True,
            "has_cohpcar_lcfo": False,
            "has_coopcar": True,
            "has_coopcar_lcfo": False,
            "has_cobicar_lcfo": False,
            "has_charge": True,
            "has_projection": False,
            "has_bandoverlaps": True,
            "has_fatbands": False,
            "has_grosspopulation": False,
            "has_polarization": False,
            "has_density_of_energies": False,
            "has_mofecar": False,
        }
        for key, item in self.lobsterout_normal.get_doc().items():
            if key not in ["has_cobicar", "has_madelung"]:
                if isinstance(item, str):
                    assert ref_data[key], item
                elif isinstance(item, int):
                    assert ref_data[key] == item
                elif key in ("charge_spilling", "total_spilling"):
                    assert item[0] == approx(ref_data[key][0])
                elif isinstance(item, list | dict):
                    assert item == ref_data[key]
        ref_data_v511 = {
            "restart_from_projection": False,
            "lobster_version": "v5.1.1",
            "threads": 8,
            "dft_program": "VASP",
            "charge_spilling": [0.0111, 0.0111],
            "total_spilling": [],
            "elements": ["N", "Al"],
            "basis_type": ["pbevaspfit2015", "pbevaspfit2015"],
            "basis_functions": [
                ["2s", "2p_y", "2p_z", "2p_x"],
                ["3s", "3p_y", "3p_z", "3p_x"],
            ],
            "timing": {
                "wall_time": {"h": "0", "min": "4", "s": "8", "ms": "368"},
                "user_time": {"h": "0", "min": "22", "s": "34", "ms": "960"},
                "sys_time": {"h": "0", "min": "0", "s": "4", "ms": "100"},
            },
            "warning_lines": [
                "2 of 1354 k-points could not be orthonormalized with an accuracy of 1.0E-5.",
                "Generally, this is not a critical error. But to help you analyze it,",
                "I dumped the band overlap matrices to the file bandOverlaps.lobster.",
                "Please check how much they deviate from the identity matrix and decide to",
                "use your results only, if you are sure that this is ok.",
            ],
            "info_orthonormalization": ["2 of 1354 k-points could not be orthonormalized with an accuracy of 1.0E-5."],
            "info_lines": [
                "No errors detected in LCFO definition...",
                "Some atoms were missing in the LCFO basis you defined.",
                "In order to avoid problems with matrix handling later,",
                "I added each remaining atom as a new fragment.",
                "This is NOT an error.",
                "You did not specify a value for the point density",
                "using the gridDensityForPrinting keyword.",
                "Therefore, I will use the default value of 0.05.",
                "You did not specify a value for the point density",
                "using the gridBufferForPrinting keyword.",
                "Therefore, I will use the default value of 3.04 Angstroms.",
            ],
            "has_doscar": True,
            "has_doscar_lso": False,
            "has_doscar_lcfo": True,
            "has_cohpcar": True,
            "has_cohpcar_lcfo": True,
            "has_coopcar": True,
            "has_coopcar_lcfo": False,
            "has_cobicar": True,
            "has_cobicar_lcfo": True,
            "has_charge": True,
            "has_madelung": True,
            "has_mofecar": True,
            "has_projection": True,
            "has_bandoverlaps": True,
            "has_fatbands": False,
            "has_grosspopulation": True,
            "has_polarization": True,
            "has_density_of_energies": False,
        }

        for key, item in self.lobsterout_v511.get_doc().items():
            if isinstance(item, str):
                assert ref_data_v511[key], item
            elif isinstance(item, int):
                assert ref_data_v511[key] == item
            elif key == "charge_spilling":
                assert item[0] == approx(ref_data_v511[key][0])
            elif isinstance(item, list | dict):
                assert item == ref_data_v511[key]

    def test_msonable(self):
        dict_data = self.lobsterout_normal.as_dict()
        lobsterout_from_dict = Lobsterout.from_dict(dict_data)
        assert dict_data == lobsterout_from_dict.as_dict()
        # test initialization with empty attributes (ensure file is not read again)
        dict_data_empty = dict.fromkeys(self.lobsterout_doscar_lso._ATTRIBUTES, None)
        lobsterout_empty_init_dict = Lobsterout.from_dict(dict_data_empty).as_dict()
        for attribute in lobsterout_empty_init_dict:
            if "@" not in attribute:
                assert lobsterout_empty_init_dict[attribute] is None

        with pytest.raises(ValueError, match="invalid=val is not a valid attribute for Lobsterout"):
            Lobsterout(filename=None, invalid="val")


class TestFatband(MatSciTest):
    def setup_method(self):
        self.structure = Vasprun(
            filename=f"{TEST_DIR}/Fatband_SiO2/Test_p_x/vasprun.xml",
            ionic_step_skip=None,
            ionic_step_offset=0,
            parse_dos=True,
            parse_eigen=False,
            parse_projected_eigen=False,
            parse_potcar_file=False,
            occu_tol=1e-8,
            exception_on_bad_xml=True,
        ).final_structure
        self.fatband_SiO2_p_x = Fatband(
            filenames=f"{TEST_DIR}/Fatband_SiO2/Test_p_x",
            kpoints_file=f"{TEST_DIR}/Fatband_SiO2/Test_p_x/KPOINTS",
            structure=self.structure,
            vasprun_file=f"{TEST_DIR}/Fatband_SiO2/Test_p_x/vasprun.xml",
        )
        self.vasprun_SiO2_p_x = Vasprun(filename=f"{TEST_DIR}/Fatband_SiO2/Test_p_x/vasprun.xml")
        self.bs_symmline = self.vasprun_SiO2_p_x.get_band_structure(line_mode=True, force_hybrid_mode=True)
        self.fatband_SiO2_p = Fatband(
            filenames=f"{TEST_DIR}/Fatband_SiO2/Test_p",
            kpoints_file=f"{TEST_DIR}/Fatband_SiO2/Test_p/KPOINTS",
            vasprun_file=f"{TEST_DIR}/Fatband_SiO2/Test_p/vasprun.xml",
            structure=self.structure,
        )
        self.fatband_SiO2_p2 = Fatband(
            filenames=f"{TEST_DIR}/Fatband_SiO2/Test_p",
            kpoints_file=f"{TEST_DIR}/Fatband_SiO2/Test_p/KPOINTS",
            structure=self.structure,
            vasprun_file=None,
            efermi=1.0647039,
        )
        self.vasprun_SiO2_p = Vasprun(filename=f"{TEST_DIR}/Fatband_SiO2/Test_p/vasprun.xml")
        self.bs_symmline2 = self.vasprun_SiO2_p.get_band_structure(line_mode=True, force_hybrid_mode=True)
        self.fatband_SiO2_spin = Fatband(
            filenames=f"{TEST_DIR}/Fatband_SiO2/Test_Spin",
            kpoints_file=f"{TEST_DIR}/Fatband_SiO2/Test_Spin/KPOINTS",
            vasprun_file=f"{TEST_DIR}/Fatband_SiO2/Test_Spin/vasprun.xml",
            structure=self.structure,
        )

        self.vasprun_SiO2_spin = Vasprun(filename=f"{TEST_DIR}/Fatband_SiO2/Test_Spin/vasprun.xml")
        self.bs_symmline_spin = self.vasprun_SiO2_p.get_band_structure(line_mode=True, force_hybrid_mode=True)

    def test_attributes(self):
        assert_allclose(list(self.fatband_SiO2_p_x.label_dict["M"]), [0.5, 0.0, 0.0])
        assert self.fatband_SiO2_p_x.efermi == self.vasprun_SiO2_p_x.efermi
        lattice1 = self.bs_symmline.lattice_rec.as_dict()
        lattice2 = self.fatband_SiO2_p_x.lattice.as_dict()
        for idx in range(3):
            assert lattice1["matrix"][idx] == approx(lattice2["matrix"][idx])
        assert self.fatband_SiO2_p_x.eigenvals[Spin.up][1][1] - self.fatband_SiO2_p_x.efermi == approx(-18.245)
        assert self.fatband_SiO2_p_x.is_spinpolarized is False
        assert_allclose(self.fatband_SiO2_p_x.kpoints_array[3], [0.03409091, 0, 0])
        assert self.fatband_SiO2_p_x.nbands == 36
        assert self.fatband_SiO2_p_x.p_eigenvals[Spin.up][2][1]["Si1"]["3p_x"] == approx(0.002)
        assert_allclose(self.fatband_SiO2_p_x.structure[0].frac_coords, [0.0, 0.47634315, 0.666667])
        assert self.fatband_SiO2_p_x.structure[0].species_string == "Si"
        assert_allclose(self.fatband_SiO2_p_x.structure[0].coords, [-1.19607309, 2.0716597, 3.67462144])

        assert_allclose(list(self.fatband_SiO2_p.label_dict["M"]), [0.5, 0.0, 0.0])
        assert self.fatband_SiO2_p.efermi == self.vasprun_SiO2_p.efermi
        lattice1 = self.bs_symmline2.lattice_rec.as_dict()
        lattice2 = self.fatband_SiO2_p.lattice.as_dict()
        for idx in range(3):
            assert lattice1["matrix"][idx] == approx(lattice2["matrix"][idx])
        assert self.fatband_SiO2_p.eigenvals[Spin.up][1][1] - self.fatband_SiO2_p.efermi == approx(-18.245)
        assert self.fatband_SiO2_p.is_spinpolarized is False
        assert_allclose(self.fatband_SiO2_p.kpoints_array[3], [0.03409091, 0, 0])
        assert self.fatband_SiO2_p.nbands == 36
        assert self.fatband_SiO2_p.p_eigenvals[Spin.up][2][1]["Si1"]["3p"] == approx(0.042)
        assert_allclose(self.fatband_SiO2_p.structure[0].frac_coords, [0.0, 0.47634315, 0.666667])
        assert self.fatband_SiO2_p.structure[0].species_string == "Si"
        assert_allclose(self.fatband_SiO2_p.structure[0].coords, [-1.19607309, 2.0716597, 3.67462144])
        assert self.fatband_SiO2_p.efermi == approx(1.0647039)

        assert list(self.fatband_SiO2_spin.label_dict["M"]) == approx([0.5, 0.0, 0.0])
        assert self.fatband_SiO2_spin.efermi == self.vasprun_SiO2_spin.efermi
        lattice1 = self.bs_symmline_spin.lattice_rec.as_dict()
        lattice2 = self.fatband_SiO2_spin.lattice.as_dict()
        for idx in range(3):
            assert lattice1["matrix"][idx] == approx(lattice2["matrix"][idx])
        assert self.fatband_SiO2_spin.eigenvals[Spin.up][1][1] - self.fatband_SiO2_spin.efermi == approx(-18.245)
        assert self.fatband_SiO2_spin.eigenvals[Spin.down][1][1] - self.fatband_SiO2_spin.efermi == approx(-18.245)
        assert self.fatband_SiO2_spin.is_spinpolarized
        assert_allclose(self.fatband_SiO2_spin.kpoints_array[3], [0.03409091, 0, 0])
        assert self.fatband_SiO2_spin.nbands == 36

        assert self.fatband_SiO2_spin.p_eigenvals[Spin.up][2][1]["Si1"]["3p"] == approx(0.042)
        assert_allclose(self.fatband_SiO2_spin.structure[0].frac_coords, [0.0, 0.47634315, 0.666667])
        assert self.fatband_SiO2_spin.structure[0].species_string == "Si"
        assert_allclose(self.fatband_SiO2_spin.structure[0].coords, [-1.19607309, 2.0716597, 3.67462144])

    def test_raises(self):
        with pytest.raises(ValueError, match="vasprun_file or efermi have to be provided"):
            Fatband(
                filenames=f"{TEST_DIR}/Fatband_SiO2/Test_Spin",
                kpoints_file=f"{TEST_DIR}/Fatband_SiO2/Test_Spin/KPOINTS",
                vasprun_file=None,
                structure=self.structure,
            )
        with pytest.raises(
            ValueError,
            match="The are two FATBAND files for the same atom and orbital. The program will stop",
        ):
            self.fatband_SiO2_p_x = Fatband(
                filenames=[
                    f"{TEST_DIR}/Fatband_SiO2/Test_p_x/FATBAND_si1_3p_x.lobster",
                    f"{TEST_DIR}/Fatband_SiO2/Test_p_x/FATBAND_si1_3p_x.lobster",
                ],
                kpoints_file=f"{TEST_DIR}/Fatband_SiO2/Test_p_x/KPOINTS",
                vasprun_file=f"{TEST_DIR}/Fatband_SiO2/Test_p_x/vasprun.xml",
                structure=self.structure,
            )

        with pytest.raises(ValueError, match="A structure object has to be provided"):
            self.fatband_SiO2_p_x = Fatband(
                filenames=[
                    f"{TEST_DIR}/Fatband_SiO2/Test_p_x/FATBAND_si1_3p_x.lobster",
                    f"{TEST_DIR}/Fatband_SiO2/Test_p_x/FATBAND_si1_3p_x.lobster",
                ],
                kpoints_file=f"{TEST_DIR}/Fatband_SiO2/Test_p_x/KPOINTS",
                vasprun_file=f"{TEST_DIR}/Fatband_SiO2/Test_p_x/vasprun.xml",
                structure=None,
            )

        with pytest.raises(
            ValueError,
            match=r"Make sure all relevant orbitals were generated and that no duplicates \(2p and 2p_x\) are present",
        ):
            self.fatband_SiO2_p_x = Fatband(
                filenames=[
                    f"{TEST_DIR}/Fatband_SiO2/Test_p_x/FATBAND_si1_3p_x.lobster",
                    f"{TEST_DIR}/Fatband_SiO2/Test_p/FATBAND_si1_3p.lobster",
                ],
                kpoints_file=f"{TEST_DIR}/Fatband_SiO2/Test_p_x/KPOINTS",
                vasprun_file=f"{TEST_DIR}/Fatband_SiO2/Test_p_x/vasprun.xml",
                structure=self.structure,
            )

        with pytest.raises(ValueError, match="No FATBAND files in folder or given"):
            self.fatband_SiO2_p_x = Fatband(
                filenames=".",
                kpoints_file=f"{TEST_DIR}/Fatband_SiO2/Test_p_x/KPOINTS",
                vasprun_file=f"{TEST_DIR}/Fatband_SiO2/Test_p_x/vasprun.xml",
                structure=self.structure,
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


class TestBandoverlaps:
    def setup_method(self):
        # test spin-polarized calc and non spin-polarized calc

        self.band_overlaps1 = Bandoverlaps(f"{TEST_DIR}/bandOverlaps.lobster.1")
        self.band_overlaps2 = Bandoverlaps(f"{TEST_DIR}/bandOverlaps.lobster.2")

        self.band_overlaps1_new = Bandoverlaps(f"{TEST_DIR}/bandOverlaps.lobster.new.1")
        self.band_overlaps2_new = Bandoverlaps(f"{TEST_DIR}/bandOverlaps.lobster.new.2")

    def test_attributes(self):
        # band_overlaps_dict
        bo_dict = self.band_overlaps1.band_overlaps_dict
        assert bo_dict[Spin.up]["max_deviations"][0] == approx(0.000278953)
        assert self.band_overlaps1_new.band_overlaps_dict[Spin.up]["max_deviations"][10] == approx(0.0640933)
        assert bo_dict[Spin.up]["matrices"][0].item(-1, -1) == approx(0.0188058)
        assert self.band_overlaps1_new.band_overlaps_dict[Spin.up]["matrices"][10].item(-1, -1) == approx(1.0)
        assert bo_dict[Spin.up]["matrices"][0].item(0, 0) == approx(1)
        assert self.band_overlaps1_new.band_overlaps_dict[Spin.up]["matrices"][10].item(0, 0) == approx(0.995849)

        assert bo_dict[Spin.down]["max_deviations"][-1] == approx(4.31567e-05)
        assert self.band_overlaps1_new.band_overlaps_dict[Spin.down]["max_deviations"][9] == approx(0.064369)
        assert bo_dict[Spin.down]["matrices"][-1].item(0, -1) == approx(4.0066e-07)
        assert self.band_overlaps1_new.band_overlaps_dict[Spin.down]["matrices"][9].item(0, -1) == approx(1.37447e-09)

        # maxDeviation
        assert self.band_overlaps1.max_deviation[0] == approx(0.000278953)
        assert self.band_overlaps1_new.max_deviation[0] == approx(0.39824)
        assert self.band_overlaps1.max_deviation[-1] == approx(4.31567e-05)
        assert self.band_overlaps1_new.max_deviation[-1] == approx(0.324898)

        assert self.band_overlaps2.max_deviation[0] == approx(0.000473319)
        assert self.band_overlaps2_new.max_deviation[0] == approx(0.403249)
        assert self.band_overlaps2.max_deviation[-1] == approx(1.48451e-05)
        assert self.band_overlaps2_new.max_deviation[-1] == approx(0.45154)

    def test_has_good_quality_maxDeviation(self):
        assert not self.band_overlaps1.has_good_quality_maxDeviation(limit_maxDeviation=0.1)
        assert not self.band_overlaps1_new.has_good_quality_maxDeviation(limit_maxDeviation=0.1)

        assert self.band_overlaps1.has_good_quality_maxDeviation(limit_maxDeviation=100)
        assert self.band_overlaps1_new.has_good_quality_maxDeviation(limit_maxDeviation=100)
        assert self.band_overlaps2.has_good_quality_maxDeviation()
        assert not self.band_overlaps2_new.has_good_quality_maxDeviation()
        assert not self.band_overlaps2.has_good_quality_maxDeviation(limit_maxDeviation=0.0000001)
        assert not self.band_overlaps2_new.has_good_quality_maxDeviation(limit_maxDeviation=0.0000001)

    def test_has_good_quality_check_occupied_bands(self):
        assert not self.band_overlaps1.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=9,
            number_occ_bands_spin_down=5,
            limit_deviation=0.1,
            spin_polarized=True,
        )
        assert not self.band_overlaps1_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=9,
            number_occ_bands_spin_down=5,
            limit_deviation=0.1,
            spin_polarized=True,
        )
        assert self.band_overlaps1.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=3,
            number_occ_bands_spin_down=0,
            limit_deviation=1,
            spin_polarized=True,
        )
        assert self.band_overlaps1_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=3,
            number_occ_bands_spin_down=0,
            limit_deviation=1,
            spin_polarized=True,
        )
        assert not self.band_overlaps1.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=1,
            number_occ_bands_spin_down=1,
            limit_deviation=1e-6,
            spin_polarized=True,
        )
        assert not self.band_overlaps1_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=1,
            number_occ_bands_spin_down=1,
            limit_deviation=1e-6,
            spin_polarized=True,
        )
        assert not self.band_overlaps1.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=1,
            number_occ_bands_spin_down=0,
            limit_deviation=1e-6,
            spin_polarized=True,
        )
        assert not self.band_overlaps1_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=1,
            number_occ_bands_spin_down=0,
            limit_deviation=1e-6,
            spin_polarized=True,
        )
        assert not self.band_overlaps1.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=0,
            number_occ_bands_spin_down=1,
            limit_deviation=1e-6,
            spin_polarized=True,
        )
        assert not self.band_overlaps1_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=0,
            number_occ_bands_spin_down=1,
            limit_deviation=1e-6,
            spin_polarized=True,
        )
        assert not self.band_overlaps1.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=4,
            number_occ_bands_spin_down=4,
            limit_deviation=1e-3,
            spin_polarized=True,
        )
        assert not self.band_overlaps1_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=4,
            number_occ_bands_spin_down=4,
            limit_deviation=1e-3,
            spin_polarized=True,
        )
        assert not self.band_overlaps2.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=10, limit_deviation=1e-7
        )
        assert not self.band_overlaps2_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=10, limit_deviation=1e-7
        )
        assert self.band_overlaps2.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=1, limit_deviation=0.1
        )

        assert not self.band_overlaps2.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=1, limit_deviation=1e-8
        )
        assert not self.band_overlaps2_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=1, limit_deviation=1e-8
        )
        assert self.band_overlaps2.has_good_quality_check_occupied_bands(number_occ_bands_spin_up=10, limit_deviation=1)
        assert self.band_overlaps2_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=2, limit_deviation=0.1
        )
        assert self.band_overlaps2.has_good_quality_check_occupied_bands(number_occ_bands_spin_up=1, limit_deviation=1)
        assert self.band_overlaps2_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=1, limit_deviation=2
        )

    def test_has_good_quality_check_occupied_bands_patched(self):
        """Test with patched data."""

        limit_deviation = 0.1

        rng = np.random.default_rng(42)  # set seed for reproducibility

        band_overlaps = copy.deepcopy(self.band_overlaps1_new)

        number_occ_bands_spin_up_all = list(range(band_overlaps.band_overlaps_dict[Spin.up]["matrices"][0].shape[0]))
        number_occ_bands_spin_down_all = list(
            range(band_overlaps.band_overlaps_dict[Spin.down]["matrices"][0].shape[0])
        )

        for actual_deviation in [0.05, 0.1, 0.2, 0.5, 1.0]:
            for spin in (Spin.up, Spin.down):
                for number_occ_bands_spin_up, number_occ_bands_spin_down in zip(
                    number_occ_bands_spin_up_all, number_occ_bands_spin_down_all, strict=False
                ):
                    for i_arr, array in enumerate(band_overlaps.band_overlaps_dict[spin]["matrices"]):
                        number_occ_bands = number_occ_bands_spin_up if spin is Spin.up else number_occ_bands_spin_down

                        shape = array.shape
                        assert np.all(np.array(shape) >= number_occ_bands)
                        assert len(shape) == 2
                        assert shape[0] == shape[1]

                        # Generate a noisy background array
                        patch_array = rng.uniform(0, 10, shape)

                        # Patch the top-left sub-array (the part that would be checked)
                        patch_array[:number_occ_bands, :number_occ_bands] = np.identity(number_occ_bands) + rng.uniform(
                            0, actual_deviation, (number_occ_bands, number_occ_bands)
                        )

                        band_overlaps.band_overlaps_dict[spin]["matrices"][i_arr] = patch_array

                    result = band_overlaps.has_good_quality_check_occupied_bands(
                        number_occ_bands_spin_up=number_occ_bands_spin_up,
                        number_occ_bands_spin_down=number_occ_bands_spin_down,
                        spin_polarized=True,
                        limit_deviation=limit_deviation,
                    )
                    # Assert for expected results
                    if (
                        (
                            actual_deviation == 0.05
                            and number_occ_bands_spin_up <= 7
                            and number_occ_bands_spin_down <= 7
                            and spin is Spin.up
                        )
                        or (actual_deviation == 0.05 and spin is Spin.down)
                        or actual_deviation == 0.1
                        or (
                            actual_deviation in [0.2, 0.5, 1.0]
                            and number_occ_bands_spin_up == 0
                            and number_occ_bands_spin_down == 0
                        )
                    ):
                        assert result
                    else:
                        assert not result

    def test_exceptions(self):
        with pytest.raises(ValueError, match="number_occ_bands_spin_down has to be specified"):
            self.band_overlaps1.has_good_quality_check_occupied_bands(
                number_occ_bands_spin_up=4,
                spin_polarized=True,
            )
        with pytest.raises(ValueError, match="number_occ_bands_spin_down has to be specified"):
            self.band_overlaps1_new.has_good_quality_check_occupied_bands(
                number_occ_bands_spin_up=4,
                spin_polarized=True,
            )

    def test_msonable(self):
        dict_data = self.band_overlaps2_new.as_dict()
        bandoverlaps_from_dict = Bandoverlaps.from_dict(dict_data)
        all_attributes = vars(self.band_overlaps2_new)
        for attr_name, attr_value in all_attributes.items():
            assert getattr(bandoverlaps_from_dict, attr_name) == attr_value

    def test_keys(self):
        bo_dict = self.band_overlaps1.band_overlaps_dict
        bo_dict_new = self.band_overlaps1_new.band_overlaps_dict
        bo_dict_2 = self.band_overlaps2.band_overlaps_dict
        assert len(bo_dict[Spin.up]["k_points"]) == 408
        assert len(bo_dict_2[Spin.up]["max_deviations"]) == 2
        assert len(bo_dict_new[Spin.down]["matrices"]) == 73


class TestGrosspop:
    def setup_method(self):
        self.grosspop1 = Grosspop(f"{TEST_DIR}/GROSSPOP.lobster")
        self.grosspop_511_sp = Grosspop(f"{TEST_DIR}/GROSSPOP_511_sp.lobster.AlN.gz")
        self.grosspop_511_nsp = Grosspop(f"{TEST_DIR}/GROSSPOP_511_nsp.lobster.NaCl.gz")
        self.grosspop_511_lcfo = Grosspop(f"{TEST_DIR}/GROSSPOP.LCFO.lobster.AlN.gz", is_lcfo=True)

    def test_attributes(self):
        gross_pop_list = self.grosspop1.list_dict_grosspop
        gross_pop_list_511_sp = self.grosspop_511_sp.list_dict_grosspop
        gross_pop_list_511_nsp = self.grosspop_511_nsp.list_dict_grosspop
        gross_pop_list_lcfo = self.grosspop_511_lcfo.list_dict_grosspop

        assert gross_pop_list[0]["Mulliken GP"]["3s"] == approx(0.52)
        assert gross_pop_list[0]["Mulliken GP"]["3p_y"] == approx(0.38)
        assert gross_pop_list[0]["Mulliken GP"]["3p_z"] == approx(0.37)
        assert gross_pop_list[0]["Mulliken GP"]["3p_x"] == approx(0.37)
        assert gross_pop_list[0]["Mulliken GP"]["total"] == approx(1.64)
        assert gross_pop_list[0]["element"] == "Si"
        assert gross_pop_list[0]["Loewdin GP"]["3s"] == approx(0.61)
        assert gross_pop_list[0]["Loewdin GP"]["3p_y"] == approx(0.52)
        assert gross_pop_list[0]["Loewdin GP"]["3p_z"] == approx(0.52)
        assert gross_pop_list[0]["Loewdin GP"]["3p_x"] == approx(0.52)
        assert gross_pop_list[0]["Loewdin GP"]["total"] == approx(2.16)
        assert gross_pop_list[5]["Mulliken GP"]["2s"] == approx(1.80)
        assert gross_pop_list[5]["Loewdin GP"]["2s"] == approx(1.60)
        assert gross_pop_list[5]["element"] == "O"
        assert gross_pop_list[8]["Mulliken GP"]["2s"] == approx(1.80)
        assert gross_pop_list[8]["Loewdin GP"]["2s"] == approx(1.60)
        assert gross_pop_list[8]["element"] == "O"

        # v5.1 spin polarized
        assert gross_pop_list_511_sp[1]["Mulliken GP"]["3p_y"][Spin.up] == approx(0.19)
        assert gross_pop_list_511_sp[3]["Loewdin GP"]["2s"][Spin.down] == approx(0.7)

        # v5.1 non spin polarized
        assert gross_pop_list_511_nsp[0]["Mulliken GP"]["3s"] == approx(0.22)
        assert gross_pop_list_511_nsp[-1]["Loewdin GP"]["total"] == approx(7.67)

        # v.5.1.1 LCFO
        assert self.grosspop_511_lcfo.is_lcfo
        assert gross_pop_list_lcfo[0]["mol"] == "AlN"
        assert gross_pop_list_lcfo[0]["Loewdin GP"]["4a1"][Spin.up] == approx(0.07)

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

        new_structure = self.grosspop1.get_structure_with_total_grosspop(f"{TEST_DIR}/POSCAR.SiO2")
        assert_allclose(new_structure.frac_coords, Structure.from_dict(struct_dict).frac_coords)

    def test_exception(self):
        structure_file = f"{VASP_IN_DIR}/POSCAR.AlN"
        with pytest.raises(ValueError, match="The GROSSPOP.LCFO.lobster data is not site wise"):
            self.grosspop_511_lcfo.get_structure_with_total_grosspop(structure_filename=structure_file)

    def test_msonable(self):
        dict_data = self.grosspop1.as_dict()
        grosspop_from_dict = Grosspop.from_dict(dict_data)
        all_attributes = vars(self.grosspop1)
        for attr_name, attr_value in all_attributes.items():
            assert getattr(grosspop_from_dict, attr_name) == attr_value


class TestIcohplist(MatSciTest):
    def setup_method(self):
        self.icohp_bise = Icohplist(filename=f"{TEST_DIR}/ICOHPLIST.lobster.BiSe")
        self.icoop_bise = Icohplist(
            filename=f"{TEST_DIR}/ICOOPLIST.lobster.BiSe",
            are_coops=True,
        )
        self.icohp_fe = Icohplist(filename=f"{TEST_DIR}/ICOHPLIST.lobster")
        # allow gzipped files
        self.icohp_gzipped = Icohplist(filename=f"{TEST_DIR}/ICOHPLIST.lobster.gz")
        self.icoop_fe = Icohplist(
            filename=f"{TEST_DIR}/ICOHPLIST.lobster",
            are_coops=True,
        )
        # ICOHPLIST.lobster from Lobster >v5
        self.icohp_aln_511_sp = Icohplist(filename=f"{TEST_DIR}/ICOHPLIST_511_sp.lobster.AlN.gz")
        self.icohp_nacl_511_nsp = Icohplist(filename=f"{TEST_DIR}/ICOHPLIST_511_nsp.lobster.NaCl.gz")

        # ICOHPLIST.LCFO.lobster from Lobster v5.1.1
        self.icohp_lcfo = Icohplist(filename=f"{TEST_DIR}/ICOHPLIST.LCFO.lobster.AlN.gz", is_lcfo=True)
        self.icohp_lcfo_non_orbitalwise = Icohplist(
            filename=f"{TEST_DIR}/ICOHPLIST_non_orbitalwise.LCFO.lobster.AlN.gz",
            is_lcfo=True,
        )

        # ICOBIs and orbitalwise ICOBILIST.lobster
        self.icobi_orbitalwise = Icohplist(
            filename=f"{TEST_DIR}/ICOBILIST.lobster",
            are_cobis=True,
        )

        self.icobi = Icohplist(
            filename=f"{TEST_DIR}/ICOBILIST.lobster.withoutorbitals",
            are_cobis=True,
        )
        self.icobi_orbitalwise_spinpolarized = Icohplist(
            filename=f"{TEST_DIR}/ICOBILIST.lobster.spinpolarized",
            are_cobis=True,
        )
        # make sure the correct line is read to check if this is a orbitalwise ICOBILIST
        self.icobi_orbitalwise_add = Icohplist(
            filename=f"{TEST_DIR}/ICOBILIST.lobster.additional_case",
            are_cobis=True,
        )
        self.icobi_orbitalwise_spinpolarized_add = Icohplist(
            filename=f"{TEST_DIR}/ICOBILIST.lobster.spinpolarized.additional_case",
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

        # >v5 ICOHPLIST
        assert self.icohp_aln_511_sp.is_spin_polarized
        assert len(self.icohp_aln_511_sp.icohplist) == 64
        assert not self.icohp_nacl_511_nsp.is_spin_polarized
        assert len(self.icohp_nacl_511_nsp.icohplist) == 152

        # v5.1.1 LCFO
        assert self.icohp_lcfo.is_lcfo
        assert self.icohp_lcfo.is_spin_polarized
        assert len(self.icohp_lcfo.icohplist) == 28
        assert not self.icohp_lcfo_non_orbitalwise.orbitalwise
        assert len(self.icohp_lcfo_non_orbitalwise.icohplist) == 28

    def test_values(self):
        icohplist_bise = {
            "1": {
                "length": 2.88231,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -2.18042},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "2": {
                "length": 3.10144,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -1.14347},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "3": {
                "length": 2.88231,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -2.18042},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "4": {
                "length": 3.10144,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -1.14348},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "5": {
                "length": 3.05001,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -1.30006},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "6": {
                "length": 2.91676,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -1.96843},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "7": {
                "length": 3.05001,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -1.30006},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "8": {
                "length": 2.91676,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -1.96843},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "9": {
                "length": 3.37522,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.47531},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "10": {
                "length": 3.07294,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -2.38796},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "11": {
                "length": 3.37522,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.47531},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
        }
        icooplist_bise = {
            "1": {
                "length": 2.88231,
                "number_of_bonds": 3,
                "icohp": {Spin.up: 0.14245},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "2": {
                "length": 3.10144,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.04118},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "3": {
                "length": 2.88231,
                "number_of_bonds": 3,
                "icohp": {Spin.up: 0.14245},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "4": {
                "length": 3.10144,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.04118},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "5": {
                "length": 3.05001,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.03516},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "6": {
                "length": 2.91676,
                "number_of_bonds": 3,
                "icohp": {Spin.up: 0.10745},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "7": {
                "length": 3.05001,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.03516},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "8": {
                "length": 2.91676,
                "number_of_bonds": 3,
                "icohp": {Spin.up: 0.10745},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "9": {
                "length": 3.37522,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.12395},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "10": {
                "length": 3.07294,
                "number_of_bonds": 3,
                "icohp": {Spin.up: 0.24714},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "11": {
                "length": 3.37522,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.12395},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
        }
        icooplist_fe = {
            "1": {
                "length": 2.83189,
                "number_of_bonds": 2,
                "icohp": {Spin.up: -0.10218, Spin.down: -0.19701},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
            "2": {
                "length": 2.45249,
                "number_of_bonds": 1,
                "icohp": {Spin.up: -0.28485, Spin.down: -0.58279},
                "translation": (0, 0, 0),
                "orbitals": None,
            },
        }

        assert icohplist_bise == self.icohp_bise.icohplist
        assert self.icohp_bise.icohpcollection.extremum_icohpvalue() == approx(-2.38796)
        assert icooplist_fe == self.icoop_fe.icohplist
        assert self.icoop_fe.icohpcollection.extremum_icohpvalue() == approx(-0.29919)
        assert icooplist_bise == self.icoop_bise.icohplist
        assert self.icoop_bise.icohpcollection.extremum_icohpvalue() == approx(0.24714)
        assert self.icobi.icohplist["1"]["icohp"][Spin.up] == approx(0.58649)
        assert self.icobi_orbitalwise.icohplist["2"]["icohp"][Spin.up] == approx(0.58649)
        assert self.icobi_orbitalwise.icohplist["1"]["icohp"][Spin.up] == approx(0.58649)
        assert self.icobi_orbitalwise_spinpolarized.icohplist["1"]["icohp"][Spin.up] == approx(0.58649 / 2, abs=1e-3)
        assert self.icobi_orbitalwise_spinpolarized.icohplist["1"]["icohp"][Spin.down] == approx(0.58649 / 2, abs=1e-3)
        assert self.icobi_orbitalwise_spinpolarized.icohplist["2"]["icohp"][Spin.down] == approx(0.58649 / 2, abs=1e-3)
        assert self.icobi.icohpcollection.extremum_icohpvalue() == approx(0.58649)
        assert self.icobi_orbitalwise_spinpolarized.icohplist["2"]["orbitals"]["2s-6s"]["icohp"][Spin.up] == approx(
            0.0247
        )

        # >v5 ICOHPLIST
        assert self.icohp_aln_511_sp.icohplist["2"]["icohp"][Spin.up] == approx(-0.21482)
        assert self.icohp_aln_511_sp.icohplist["2"]["icohp"][Spin.down] == approx(-0.21493)
        assert self.icohp_nacl_511_nsp.icohplist["13"]["icohp"][Spin.up] == approx(-0.03021)
        assert self.icohp_nacl_511_nsp.icohplist["10"]["orbitals"]["3s-2py"]["icohp"][Spin.up] == approx(-0.00113)

        # v5.1.1 LCFO
        assert self.icohp_lcfo.icohplist["15"]["orbitals"]["3a1-3p"]["icohp"][Spin.up] == approx(-0.00586)
        assert self.icohp_lcfo_non_orbitalwise.icohplist["16"]["icohp"][Spin.up] == approx(-0.2983)
        assert self.icohp_lcfo_non_orbitalwise.icohplist["16"]["icohp"][Spin.down] == approx(-0.29842)

    def test_msonable(self):
        dict_data = self.icobi_orbitalwise_spinpolarized.as_dict()
        icohplist_from_dict = Icohplist.from_dict(dict_data)
        all_attributes = vars(self.icobi_orbitalwise_spinpolarized)
        for attr_name, attr_value in all_attributes.items():
            if isinstance(attr_value, IcohpCollection):
                assert getattr(icohplist_from_dict, attr_name).as_dict() == attr_value.as_dict()
            else:
                assert getattr(icohplist_from_dict, attr_name) == attr_value

    def test_missing_trailing_newline(self):
        fname = f"{self.tmp_path}/icohplist"
        with open(fname, mode="w", encoding="utf-8") as f:
            f.write(
                "1   Co1   O1   1.00000   0   0   0   -0.50000   -1.00000\n"
                "2   Co2   O2   1.10000   0   0   0   -0.60000   -1.10000"
            )

        ip = Icohplist(filename=fname)
        assert len(ip.icohplist) == 2
        assert ip.icohplist["1"]["icohp"][Spin.up] == approx(-0.5)


class TestNciCobiList:
    def setup_method(self):
        self.ncicobi = NciCobiList(filename=f"{TEST_DIR}/NcICOBILIST.lobster")
        self.ncicobi_gz = NciCobiList(filename=f"{TEST_DIR}/NcICOBILIST.lobster.gz")
        self.ncicobi_no_spin = NciCobiList(filename=f"{TEST_DIR}/NcICOBILIST.lobster.nospin")
        self.ncicobi_no_spin_wo = NciCobiList(filename=f"{TEST_DIR}/NcICOBILIST.lobster.nospin.withoutorbitals")
        self.ncicobi_wo = NciCobiList(filename=f"{TEST_DIR}/NcICOBILIST.lobster.withoutorbitals")

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


class TestWavefunction(MatSciTest):
    def test_parse_file(self):
        grid, points, real, imaginary, distance = Wavefunction._parse_file(
            f"{TEST_DIR}/LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz"
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
            filename=f"{TEST_DIR}/LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR_O.gz"),
        )

        wave1.set_volumetric_data(grid=wave1.grid, structure=wave1.structure)
        assert wave1.volumetricdata_real.data["total"][0, 0, 0] == approx(-3.0966)
        assert wave1.volumetricdata_imaginary.data["total"][0, 0, 0] == approx(-6.45895e00)

    def test_get_volumetricdata_real(self):
        wave1 = Wavefunction(
            filename=f"{TEST_DIR}/LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR_O.gz"),
        )
        volumetricdata_real = wave1.get_volumetricdata_real()
        assert volumetricdata_real.data["total"][0, 0, 0] == approx(-3.0966)

    def test_get_volumetricdata_imaginary(self):
        wave1 = Wavefunction(
            filename=f"{TEST_DIR}/LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR_O.gz"),
        )
        volumetricdata_imaginary = wave1.get_volumetricdata_imaginary()
        assert volumetricdata_imaginary.data["total"][0, 0, 0] == approx(-6.45895e00)

    def test_get_volumetricdata_density(self):
        wave1 = Wavefunction(
            filename=f"{TEST_DIR}/LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR_O.gz"),
        )
        volumetricdata_density = wave1.get_volumetricdata_density()
        assert volumetricdata_density.data["total"][0, 0, 0] == approx((-3.0966 * -3.0966) + (-6.45895 * -6.45895))

    def test_write_file(self):
        wave1 = Wavefunction(
            filename=f"{TEST_DIR}/LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR_O.gz"),
        )
        real_wavecar_path = f"{self.tmp_path}/real-wavecar.vasp"
        wave1.write_file(filename=real_wavecar_path, part="real")
        assert os.path.isfile(real_wavecar_path)

        imag_wavecar_path = f"{self.tmp_path}/imaginary-wavecar.vasp"
        wave1.write_file(filename=imag_wavecar_path, part="imaginary")
        assert os.path.isfile(imag_wavecar_path)

        density_wavecar_path = f"{self.tmp_path}/density-wavecar.vasp"
        wave1.write_file(filename=density_wavecar_path, part="density")
        assert os.path.isfile(density_wavecar_path)


class TestSitePotentials(MatSciTest):
    def setup_method(self) -> None:
        self.sitepotential = SitePotential(filename=f"{TEST_DIR}/SitePotentials.lobster.perovskite")

    def test_attributes(self):
        assert self.sitepotential.sitepotentials_loewdin == [
            -8.77,
            -17.08,
            9.57,
            9.57,
            8.45,
        ]
        assert self.sitepotential.sitepotentials_mulliken == [
            -11.38,
            -19.62,
            11.18,
            11.18,
            10.09,
        ]
        assert self.sitepotential.madelungenergies_loewdin == approx(-28.64)
        assert self.sitepotential.madelungenergies_mulliken == approx(-40.02)
        assert self.sitepotential.atomlist == ["La1", "Ta2", "N3", "N4", "O5"]
        assert self.sitepotential.types == ["La", "Ta", "N", "N", "O"]
        assert self.sitepotential.num_atoms == 5
        assert self.sitepotential.ewald_splitting == approx(3.14)

    def test_get_structure(self):
        structure = self.sitepotential.get_structure_with_site_potentials(f"{TEST_DIR}/POSCAR.perovskite")
        assert structure.site_properties["Loewdin Site Potentials (eV)"] == [
            -8.77,
            -17.08,
            9.57,
            9.57,
            8.45,
        ]
        assert structure.site_properties["Mulliken Site Potentials (eV)"] == [
            -11.38,
            -19.62,
            11.18,
            11.18,
            10.09,
        ]

    def test_msonable(self):
        dict_data = self.sitepotential.as_dict()
        sitepotential_from_dict = SitePotential.from_dict(dict_data)
        all_attributes = vars(self.sitepotential)
        for attr_name, attr_value in all_attributes.items():
            assert getattr(sitepotential_from_dict, attr_name) == attr_value


class TestMadelungEnergies(MatSciTest):
    def setup_method(self) -> None:
        self.madelungenergies = MadelungEnergies(filename=f"{TEST_DIR}/MadelungEnergies.lobster.perovskite")

    def test_attributes(self):
        assert self.madelungenergies.madelungenergies_loewdin == approx(-28.64)
        assert self.madelungenergies.madelungenergies_mulliken == approx(-40.02)
        assert self.madelungenergies.ewald_splitting == approx(3.14)

    def test_msonable(self):
        dict_data = self.madelungenergies.as_dict()
        madelung_from_dict = MadelungEnergies.from_dict(dict_data)
        all_attributes = vars(self.madelungenergies)
        for attr_name, attr_value in all_attributes.items():
            assert getattr(madelung_from_dict, attr_name) == attr_value


class TestLobsterMatrices(MatSciTest):
    def setup_method(self) -> None:
        self.hamilton_matrices = LobsterMatrices(
            filename=f"{TEST_DIR}/Na_hamiltonMatrices.lobster.gz", e_fermi=-2.79650354
        )
        self.transfer_matrices = LobsterMatrices(filename=f"{TEST_DIR}/C_transferMatrices.lobster.gz")
        self.overlap_matrices = LobsterMatrices(filename=f"{TEST_DIR}/Si_overlapMatrices.lobster.gz")
        self.coeff_matrices = LobsterMatrices(filename=f"{TEST_DIR}/Si_coefficientMatricesLSO1.lobster.gz")

    def test_attributes(self):
        # hamilton matrices
        assert self.hamilton_matrices.average_onsite_energies == approx(
            {
                "Na1_3s": 0.58855353,
                "Na1_2p_y": -25.72719646,
                "Na1_2p_z": -25.72719646,
                "Na1_2p_x": -25.72719646,
            }
        )
        ref_onsite_energies = [
            [-0.22519646, -25.76989646, -25.76989646, -25.76989646],
            [1.40230354, -25.68449646, -25.68449646, -25.68449646],
        ]
        assert_allclose(self.hamilton_matrices.onsite_energies, ref_onsite_energies)

        ref_imag_mat_spin_up = np.zeros((4, 4))

        assert_allclose(
            self.hamilton_matrices.hamilton_matrices["1"][Spin.up].imag,
            ref_imag_mat_spin_up,
        )

        ref_real_mat_spin_up = [
            [-3.0217, 0.0, 0.0, 0.0],
            [0.0, -28.5664, 0.0, 0.0],
            [0.0, 0.0, -28.5664, 0.0],
            [0.0, 0.0, 0.0, -28.5664],
        ]
        assert_allclose(
            self.hamilton_matrices.hamilton_matrices["1"][Spin.up].real,
            ref_real_mat_spin_up,
        )

        # overlap matrices
        assert self.overlap_matrices.average_onsite_overlaps == approx(
            {
                "Si1_3s": 1.00000009,
                "Si1_3p_y": 0.99999995,
                "Si1_3p_z": 0.99999995,
                "Si1_3p_x": 0.99999995,
            }
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

        assert_allclose(
            self.transfer_matrices.transfer_matrices["1"][Spin.down].imag,
            ref_imag_mat_spin_down,
        )

        ref_real_mat_spin_down = [
            [-0.03735031, 0.0, 0.0, 0.0],
            [0.0, -0.66865552, 0.06086057, 0.13042529],
            [-0.0, 0.04262018, 0.69253776, -0.12491928],
            [0.0, 0.11473763, 0.09742773, 0.80648063],
        ]

        assert_allclose(
            self.transfer_matrices.transfer_matrices["1"][Spin.down].real,
            ref_real_mat_spin_down,
        )

        # coefficient matrices
        assert list(self.coeff_matrices.coefficient_matrices["1"]) == [
            Spin.up,
            Spin.down,
        ]
        assert self.coeff_matrices.average_onsite_coefficient == approx(
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

        assert_allclose(
            self.coeff_matrices.coefficient_matrices["1"][Spin.up].imag,
            ref_imag_mat_spin_up,
        )

        ref_real_mat_spin_up = [
            [0.80226096, 0.0, 0.0, 0.0],
            [0.0, -0.33931137, -0.42979933, -0.34286226],
            [0.0, 0.30354633, -0.48472536, 0.29861248],
            [0.0, -0.32075579, -0.00405544, 0.64528776],
        ]

        assert_allclose(
            self.coeff_matrices.coefficient_matrices["1"][Spin.up].real,
            ref_real_mat_spin_up,
        )

    def test_raises(self):
        with pytest.raises(ValueError, match="Please provide the fermi energy in eV"):
            self.hamilton_matrices = LobsterMatrices(filename=f"{TEST_DIR}/Na_hamiltonMatrices.lobster.gz")

        with pytest.raises(
            RuntimeError,
            match="Please check provided input file, it seems to be empty",
        ):
            self.hamilton_matrices = LobsterMatrices(filename=f"{TEST_DIR}/hamiltonMatrices.lobster")


class TestPolarization(MatSciTest):
    def setup_method(self) -> None:
        self.polarization = Polarization(filename=f"{TEST_DIR}/POLARIZATION.lobster.AlN.gz")

    def test_attributes(self):
        assert self.polarization.rel_loewdin_pol_vector == {
            "x": -0.0,
            "y": -0.01,
            "z": 45.62,
            "abs": 45.62,
            "unit": "uC/cm2",
        }
        assert self.polarization.rel_mulliken_pol_vector == {
            "x": -0.0,
            "y": -0.02,
            "z": 56.14,
            "abs": 56.14,
            "unit": "uC/cm2",
        }


def test_get_lines():
    """Ensure `_get_lines` is not trailing end char sensitive."""
    with open("without-end-char", mode="wb") as f:
        f.write(b"first line\nsecond line")

    with open("with-end-char", mode="wb") as f:
        f.write(b"first line\nsecond line\n")

    without_end_char = _get_lines("without-end-char")
    with_end_char = _get_lines("with-end-char")

    assert len(with_end_char) == len(without_end_char) == 2
