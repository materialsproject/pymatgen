from __future__ import annotations

import copy
import json
import os
from typing import TYPE_CHECKING

import numpy as np
import pytest
from monty.json import MontyEncoder
from numpy.testing import assert_allclose, assert_array_equal
from pytest import approx

from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.lobster.future import (
    BWDF,
    CHARGE,
    CHARGE_LCFO,
    COBICAR,
    COBICAR_LCFO,
    COOPCAR,
    DOSCAR,
    DOSCAR_LCFO,
    GROSSPOP,
    GROSSPOP_LCFO,
    ICOBILIST,
    ICOHPLIST,
    ICOHPLIST_LCFO,
    ICOOPLIST,
    POLARIZATION,
    BandOverlaps,
    Fatband,
    Fatbands,
    LobsterMatrices,
    LobsterOut,
    MadelungEnergies,
    NcICOBILIST,
    SitePotentials,
    Wavefunction,
)
from pymatgen.io.vasp import Vasprun
from pymatgen.util.testing import TEST_FILES_DIR, VASP_OUT_DIR, MatSciTest

if TYPE_CHECKING:
    from monty.json import MSONable

TEST_DIR = f"{TEST_FILES_DIR}/electronic_structure/cohp"

__author__ = "Janine George, Marco Esters"
__copyright__ = "Copyright 2017, The Materials Project"
__version__ = "0.2"
__email__ = "janine.george@uclouvain.be, esters@uoregon.edu"
__date__ = "Dec 10, 2017"


class TestBWDF(MatSciTest):
    """Test BWDF and BWDFCOHP classes."""

    def setup_method(self):
        self.bwdf_coop = BWDF(filename=f"{TEST_DIR}/BWDF.lobster.AlN.gz")
        self.bwdf_cohp = BWDF(filename=f"{TEST_DIR}/BWDFCOHP.lobster.NaCl.gz")

    def test_attributes(self):
        assert len(self.bwdf_coop.centers) == 201
        assert len(self.bwdf_cohp.centers) == 143
        assert self.bwdf_coop.bwdf[Spin.down][-2] == approx(0.00082, abs=1e-4)
        assert self.bwdf_coop.bwdf[Spin.up][0] == approx(0.81161, abs=1e-4)
        assert self.bwdf_cohp.bwdf[Spin.up][103] == approx(-0.01392, abs=1e-4)


class TestDOSCAR(MatSciTest):
    def setup_method(self):
        doscar = f"{VASP_OUT_DIR}/DOSCAR.lobster.spin"
        doscar2 = f"{VASP_OUT_DIR}/DOSCAR.lobster.nonspin"
        doscar3 = f"{VASP_OUT_DIR}/DOSCAR.LCFO.lobster.AlN"

        self.doscar_spin_pol = DOSCAR(filename=doscar)
        self.doscar_nonspin_pol = DOSCAR(filename=doscar2)
        self.doscar_lcfo = DOSCAR_LCFO(filename=doscar3)
        self.doscar_spin_pol2 = DOSCAR(filename=doscar)

    def test_pdos(self):
        """Test projected densities of states (PDOS) from DOSCAR files."""
        expected_pdos_spin = {
            "2s": {
                Spin.up: [0.00000, 0.00159, 0.00000, 0.00011, 0.00000, 0.00069],
                Spin.down: [0.00000, 0.00159, 0.00000, 0.00011, 0.00000, 0.00069],
            },
            "2p_y": {
                Spin.up: [0.00000, 0.00160, 0.00000, 0.25801, 0.00000, 0.00029],
                Spin.down: [0.00000, 0.00161, 0.00000, 0.25819, 0.00000, 0.00029],
            },
            "2p_z": {
                Spin.up: [0.00000, 0.00161, 0.00000, 0.25823, 0.00000, 0.00029],
                Spin.down: [0.00000, 0.00160, 0.00000, 0.25795, 0.00000, 0.00029],
            },
            "2p_x": {
                Spin.up: [0.00000, 0.00160, 0.00000, 0.25805, 0.00000, 0.00029],
                Spin.down: [0.00000, 0.00161, 0.00000, 0.25814, 0.00000, 0.00029],
            },
        }

        expected_pdos_nonspin = {
            "2s": [0.00000, 0.00320, 0.00000, 0.00017, 0.00000, 0.00060],
            "2p_y": [0.00000, 0.00322, 0.00000, 0.51635, 0.00000, 0.00037],
            "2p_z": [0.00000, 0.00322, 0.00000, 0.51636, 0.00000, 0.00037],
            "2p_x": [0.00000, 0.00322, 0.00000, 0.51634, 0.00000, 0.00037],
        }

        for orbital, spin_data in expected_pdos_spin.items():
            for spin, expected_values in spin_data.items():
                assert_allclose(
                    self.doscar_spin_pol.projected_dos["F1"][orbital].densities[spin],
                    expected_values,
                )

        for orbital, expected_values in expected_pdos_nonspin.items():
            assert_allclose(
                self.doscar_nonspin_pol.projected_dos["F1"][orbital].densities[Spin.up],
                expected_values,
            )

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

        assert self.doscar_lcfo.is_lcfo
        assert_allclose(
            self.doscar_lcfo.projected_dos["AlN_1"]["1a1"].densities[Spin.down], pdos_1a1_AlN
        )
        assert_allclose(
            self.doscar_lcfo.projected_dos["Al_1"]["3p_y"].densities[Spin.down], pdos_3py_Al
        )
        assert_allclose(
            self.doscar_lcfo.projected_dos["N_1"]["2s"].densities[Spin.down], pdos_2s_N
        )

    def test_tdos(self):
        """Test total densities of states (TDOS) from DOSCAR files."""
        energies_spin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        tdos_up = [0.00000, 0.79999, 0.00000, 0.79999, 0.00000, 0.02577]
        tdos_down = [0.00000, 0.79999, 0.00000, 0.79999, 0.00000, 0.02586]

        assert_allclose(energies_spin, self.doscar_spin_pol.total_dos.energies)
        assert_allclose(tdos_up, self.doscar_spin_pol.total_dos.densities[Spin.up])
        assert_allclose(tdos_down, self.doscar_spin_pol.total_dos.densities[Spin.down])
        # assert fermi == approx(self.doscar_spin_pol.total_dos.efermi)

        energies_nonspin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        tdos_nonspin = [0.00000, 1.60000, 0.00000, 1.60000, 0.00000, 0.02418]

        assert_allclose(energies_nonspin, self.doscar_nonspin_pol.total_dos.energies)
        assert_allclose(tdos_nonspin, self.doscar_nonspin_pol.total_dos.densities[Spin.up])
        # assert fermi == approx(self.doscar_nonspin_pol.total_dos.efermi)

    def test_energies(self):
        """Test energies from DOSCAR files."""
        energies_spin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]

        assert_allclose(energies_spin, self.doscar_spin_pol.energies)

        energies_nonspin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        assert_allclose(energies_nonspin, self.doscar_nonspin_pol.energies)

    def test_itdensities(self):
        """Test integrated total densities from DOSCAR files."""
        itdos_up = [1.99997, 4.99992, 4.99992, 7.99987, 7.99987, 8.09650]
        itdos_down = [1.99997, 4.99992, 4.99992, 7.99987, 7.99987, 8.09685]
        assert_allclose(itdos_up, self.doscar_spin_pol.integrated_total_dos.densities[Spin.up])
        assert_allclose(itdos_down, self.doscar_spin_pol.integrated_total_dos.densities[Spin.down])

        itdos_nonspin = [4.00000, 10.00000, 10.00000, 16.00000, 16.00000, 16.09067]
        assert_allclose(itdos_nonspin, self.doscar_nonspin_pol.integrated_total_dos.densities[Spin.up])

    def test_is_spin_polarized(self):
        """Test is_spin_polarized attribute from DOSCAR files."""
        assert self.doscar_spin_pol.is_spin_polarized
        assert not self.doscar_nonspin_pol.is_spin_polarized


class TestCHARGE(MatSciTest):
    def setup_method(self):
        """Setup for CHARGE and CHARGE_LCFO tests."""
        self.charge2 = CHARGE(filename=f"{TEST_DIR}/CHARGE.lobster.MnO")
        self.charge = CHARGE(filename=f"{TEST_DIR}/CHARGE.lobster.MnO2.gz")
        self.charge_lcfo = CHARGE_LCFO(
            filename=f"{TEST_DIR}/CHARGE.LCFO.lobster.ALN.gz"
        )

    def test_attributes(self):
        """Test attributes of CHARGE and CHARGE_LCFO classes."""
        charge_loewdin = [-1.25, 1.25]
        charge_mulliken = [-1.30, 1.30]

        assert charge_mulliken == self.charge2.mulliken
        assert charge_loewdin == self.charge2.loewdin

        assert self.charge_lcfo.is_lcfo

        assert_allclose(self.charge_lcfo.loewdin, [0.0, 1.02, -1.02])
        assert not self.charge_lcfo.mulliken

    def test_msonable(self):
        """Test MSONable functionality of CHARGE class."""
        dict_data = self.charge2.as_dict()
        charge_from_dict = CHARGE.from_dict(dict_data)
        all_attributes = vars(self.charge2)
        for attr_name, attr_value in all_attributes.items():
            assert getattr(charge_from_dict, attr_name) == attr_value


class TestLobsterOut(MatSciTest):
    def setup_method(self):
        """Setup for LobsterOut tests."""
        self.lobsterout_normal = LobsterOut(filename=f"{TEST_DIR}/lobsterout.normal")
        # make sure .gz files are also read correctly
        self.lobsterout_normal = LobsterOut(
            filename=f"{TEST_DIR}/lobsterout.normal2.gz"
        )
        self.lobsterout_fatband_grosspop_densityofenergies = LobsterOut(
            filename=f"{TEST_DIR}/lobsterout.fatband_grosspop_densityofenergy"
        )
        self.lobsterout_saveprojection = LobsterOut(
            filename=f"{TEST_DIR}/lobsterout.saveprojection"
        )
        self.lobsterout_skipping_all = LobsterOut(
            filename=f"{TEST_DIR}/lobsterout.skipping_all"
        )
        self.lobsterout_twospins = LobsterOut(
            filename=f"{TEST_DIR}/lobsterout.twospins"
        )
        self.lobsterout_GaAs = LobsterOut(filename=f"{TEST_DIR}/lobsterout.GaAs")
        self.lobsterout_from_projection = LobsterOut(
            filename=f"{TEST_DIR}/lobsterout_from_projection"
        )
        self.lobsterout_onethread = LobsterOut(
            filename=f"{TEST_DIR}/lobsterout.onethread"
        )
        self.lobsterout_cobi_madelung = LobsterOut(
            filename=f"{TEST_DIR}/lobsterout_cobi_madelung"
        )
        self.lobsterout_doscar_lso = LobsterOut(
            filename=f"{TEST_DIR}/lobsterout_doscar_lso"
        )

        # TODO: implement skipping madelung/cobi
        self.lobsterout_skipping_cobi_madelung = LobsterOut(
            filename=f"{TEST_DIR}/lobsterout.skip_cobi_madelung"
        )

        self.lobsterout_v511 = LobsterOut(filename=f"{TEST_DIR}/lobsterout_v511.gz")

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
        assert self.lobsterout_normal.lobster_version == "3.1.0"
        assert self.lobsterout_normal.number_of_spins == 1
        assert self.lobsterout_normal.number_of_threads == 8
        assert self.lobsterout_normal.timing == {
            "wall_time": {"h": "0", "min": "0", "s": "2", "ms": "702"},
            "user_time": {"h": "0", "min": "0", "s": "20", "ms": "330"},
            "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "310"},
        }
        assert self.lobsterout_normal.total_spilling[0] == approx(
            [0.044000000000000004][0]
        )
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
        assert self.lobsterout_fatband_grosspop_densityofenergies.basis_type == [
            "pbeVaspFit2015"
        ]
        assert_allclose(
            self.lobsterout_fatband_grosspop_densityofenergies.charge_spilling, [0.0268]
        )
        assert self.lobsterout_fatband_grosspop_densityofenergies.dft_program == "VASP"
        assert self.lobsterout_fatband_grosspop_densityofenergies.elements == ["Ti"]
        assert self.lobsterout_fatband_grosspop_densityofenergies.has_charge
        assert not self.lobsterout_fatband_grosspop_densityofenergies.has_cohpcar
        assert not self.lobsterout_fatband_grosspop_densityofenergies.has_coopcar
        assert not self.lobsterout_fatband_grosspop_densityofenergies.has_doscar
        assert not self.lobsterout_fatband_grosspop_densityofenergies.has_projection
        assert self.lobsterout_fatband_grosspop_densityofenergies.has_bandoverlaps
        assert (
            self.lobsterout_fatband_grosspop_densityofenergies.has_density_of_energies
        )
        assert self.lobsterout_fatband_grosspop_densityofenergies.has_fatbands
        assert self.lobsterout_fatband_grosspop_densityofenergies.has_grosspopulation
        assert self.lobsterout_fatband_grosspop_densityofenergies.info_lines == [
            "There are more PAW bands than local basis functions available.",
            "To prevent trouble in orthonormalization and Hamiltonian reconstruction",
            "the PAW bands from 21 and upwards will be ignored.",
        ]
        assert (
            self.lobsterout_fatband_grosspop_densityofenergies.info_orthonormalization
            == [
                "3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5."
            ]
        )
        assert (
            not self.lobsterout_fatband_grosspop_densityofenergies.is_restart_from_projection
        )
        assert (
            self.lobsterout_fatband_grosspop_densityofenergies.lobster_version
            == "3.1.0"
        )
        assert self.lobsterout_fatband_grosspop_densityofenergies.number_of_spins == 1
        assert self.lobsterout_fatband_grosspop_densityofenergies.number_of_threads == 8
        assert self.lobsterout_fatband_grosspop_densityofenergies.timing == {
            "wall_time": {"h": "0", "min": "0", "s": "4", "ms": "136"},
            "user_time": {"h": "0", "min": "0", "s": "18", "ms": "280"},
            "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "290"},
        }
        assert self.lobsterout_fatband_grosspop_densityofenergies.total_spilling[
            0
        ] == approx([0.044000000000000004][0])
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
        assert self.lobsterout_saveprojection.lobster_version == "3.1.0"
        assert self.lobsterout_saveprojection.number_of_spins == 1
        assert self.lobsterout_saveprojection.number_of_threads == 8
        assert self.lobsterout_saveprojection.timing == {
            "wall_time": {"h": "0", "min": "0", "s": "2", "ms": "574"},
            "user_time": {"h": "0", "min": "0", "s": "18", "ms": "250"},
            "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "320"},
        }
        assert self.lobsterout_saveprojection.total_spilling[0] == approx(
            [0.044000000000000004][0]
        )
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
        assert self.lobsterout_skipping_all.lobster_version == "3.1.0"
        assert self.lobsterout_skipping_all.number_of_spins == 1
        assert self.lobsterout_skipping_all.number_of_threads == 8
        assert self.lobsterout_skipping_all.timing == {
            "wall_time": {"h": "0", "min": "0", "s": "2", "ms": "117"},
            "user_time": {"h": "0", "min": "0", "s": "16", "ms": "79"},
            "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "320"},
        }
        assert self.lobsterout_skipping_all.total_spilling[0] == approx(
            [0.044000000000000004][0]
        )
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
        assert self.lobsterout_twospins.charge_spilling[0] == approx(
            0.36619999999999997
        )
        assert self.lobsterout_twospins.charge_spilling[1] == approx(
            0.36619999999999997
        )
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
        assert self.lobsterout_twospins.lobster_version == "3.1.0"
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
        assert self.lobsterout_from_projection.lobster_version == "3.1.0"
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
        assert self.lobsterout_GaAs.lobster_version == "3.1.0"
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

    def test_msonable(self):
        """Test the as_dict and from_dict methods for Lobsterout."""
        dict_data = self.lobsterout_normal.as_dict()

        lobsterout_from_dict = LobsterOut.from_dict(dict_data)
        assert dict_data == lobsterout_from_dict.as_dict()

        with pytest.raises(
            TypeError, match="got an unexpected keyword argument 'invalid'"
        ):
            LobsterOut(filename=None, invalid="invalid")  # type: ignore[testing]


class TestFatbands(MatSciTest):
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
        self.fatbands_sio2_p_x = Fatbands(
            directory=f"{TEST_DIR}/Fatband_SiO2/Test_p_x",
            structure=self.structure,
        )
        self.vasprun_sio2_p_x = Vasprun(
            filename=f"{TEST_DIR}/Fatband_SiO2/Test_p_x/vasprun.xml"
        )
        self.bs_symmline = self.vasprun_sio2_p_x.get_band_structure(
            line_mode=True, force_hybrid_mode=True
        )
        self.fatbands_sio2_p = Fatbands(
            directory=f"{TEST_DIR}/Fatband_SiO2/Test_p",
            structure=self.structure,
        )
        self.single_fatband = Fatband(
            filename=f"{TEST_DIR}/Fatband_SiO2/Test_p/FATBAND_o4_2p.lobster"
        )
        self.vasprun_sio2_p = Vasprun(
            filename=f"{TEST_DIR}/Fatband_SiO2/Test_p/vasprun.xml"
        )
        self.bs_symmline2 = self.vasprun_sio2_p.get_band_structure(
            line_mode=True, force_hybrid_mode=True
        )
        self.fatbands_sio2_spin = Fatbands(
            directory=f"{TEST_DIR}/Fatband_SiO2/Test_Spin",
            structure=self.structure,
        )

        self.vasprun_sio2_spin = Vasprun(
            filename=f"{TEST_DIR}/Fatband_SiO2/Test_Spin/vasprun.xml"
        )
        self.bs_symmline_spin = self.vasprun_sio2_p.get_band_structure(
            line_mode=True, force_hybrid_mode=True
        )

    def test_attributes(self):
        assert self.fatbands_sio2_p_x.efermi == self.vasprun_sio2_p_x.efermi
        lattice1 = self.bs_symmline.lattice_rec.as_dict()
        lattice2 = self.fatbands_sio2_p_x.reciprocal_lattice.as_dict()
        for idx in range(3):
            assert lattice1["matrix"][idx] == approx(lattice2["matrix"][idx])

        assert self.fatbands_sio2_p_x.fatbands[1]["energies"][Spin.up][1][1] == approx(
            -18.245
        )
        assert len(self.fatbands_sio2_p_x.spins) == 1
        assert_allclose(self.fatbands_sio2_p_x.kpoints.kpts[3], [0.03409091, 0, 0])

        assert len(self.fatbands_sio2_p_x.fatbands[0]["projections"][Spin.up]) == len(
            self.fatbands_sio2_p_x.kpoints.kpts
        )

        assert self.single_fatband.nbands == 36
        assert self.single_fatband.center == "O4"
        assert self.single_fatband.orbital == "2p"
        assert len(self.single_fatband.spins) == 1
        assert len(self.single_fatband.fatband["energies"][Spin.up][0]) == 36

        assert self.fatbands_sio2_p_x.fatbands[-1]["center"] == "Si3"
        assert self.fatbands_sio2_p_x.fatbands[-1]["orbital"] == "3s"
        assert self.fatbands_sio2_p_x.fatbands[-1]["projections"][Spin.up][2][
            1
        ] == approx(0.013)
        assert self.fatbands_sio2_p_x.fatbands[-1]["energies"][Spin.up][2][2] == approx(
            -18.245
        )
        assert_allclose(
            self.fatbands_sio2_p_x.structure[0].frac_coords, [0.0, 0.47634315, 0.666667]
        )
        assert self.fatbands_sio2_p_x.structure[0].species_string == "Si"
        assert_allclose(
            self.fatbands_sio2_p_x.structure[0].coords,
            [-1.19607309, 2.0716597, 3.67462144],
        )

        assert self.fatbands_sio2_p.efermi == self.vasprun_sio2_p.efermi
        lattice1 = self.bs_symmline2.lattice_rec.as_dict()
        lattice2 = self.fatbands_sio2_p.reciprocal_lattice.as_dict()
        for idx in range(3):
            assert lattice1["matrix"][idx] == approx(lattice2["matrix"][idx])
        assert self.fatbands_sio2_p.fatbands[0]["energies"][Spin.up][1][1] == approx(
            -18.245
        )
        assert len(self.fatbands_sio2_p.spins) == 1
        assert_allclose(self.fatbands_sio2_p.kpoints.kpts[3], [0.03409091, 0, 0])

        assert self.fatbands_sio2_p.fatbands[2]["projections"][Spin.up][2][1] == approx(
            0.002
        )
        assert_allclose(
            self.fatbands_sio2_p.structure[0].frac_coords, [0.0, 0.47634315, 0.666667]
        )
        assert self.fatbands_sio2_p.structure[0].species_string == "Si"
        assert_allclose(
            self.fatbands_sio2_p.structure[0].coords,
            [-1.19607309, 2.0716597, 3.67462144],
        )
        assert self.fatbands_sio2_p.efermi == approx(1.0647039)

        assert self.fatbands_sio2_spin.efermi == self.vasprun_sio2_spin.efermi
        lattice1 = self.bs_symmline_spin.lattice_rec.as_dict()
        lattice2 = self.fatbands_sio2_spin.reciprocal_lattice.as_dict()
        for idx in range(3):
            assert lattice1["matrix"][idx] == approx(lattice2["matrix"][idx])
        assert self.fatbands_sio2_spin.fatbands[1]["energies"][Spin.up][1][1] == approx(
            -18.245
        )
        assert self.fatbands_sio2_spin.fatbands[0]["energies"][Spin.down][1][
            1
        ] == approx(-18.245)
        assert len(self.fatbands_sio2_spin.spins) == 2
        assert_allclose(self.fatbands_sio2_spin.kpoints.kpts[3], [0.03409091, 0, 0])
        assert len(self.fatbands_sio2_spin.fatbands[0]["energies"][Spin.up][0]) == 36
        assert len(self.fatbands_sio2_spin.fatbands[0]["energies"][Spin.down][0]) == 36

        assert self.fatbands_sio2_spin.fatbands[0]["projections"][Spin.up][2][
            1
        ] == approx(0.003)
        assert_allclose(
            self.fatbands_sio2_spin.structure[0].frac_coords,
            [0.0, 0.47634315, 0.666667],
        )
        assert self.fatbands_sio2_spin.structure[0].species_string == "Si"
        assert_allclose(
            self.fatbands_sio2_spin.structure[0].coords,
            [-1.19607309, 2.0716597, 3.67462144],
        )


class TestBandOverlaps(MatSciTest):
    def setup_method(self):
        with pytest.raises(RuntimeError, match="could not convert string"):
            self.band_overlaps1 = BandOverlaps(f"{TEST_DIR}/bandOverlaps.lobster.1")

        self.band_overlaps1 = BandOverlaps(
            f"{TEST_DIR}/bandOverlaps.lobster.1", process_immediately=False
        )
        self.band_overlaps2 = BandOverlaps(
            f"{TEST_DIR}/bandOverlaps.lobster.2", process_immediately=False
        )

        self.band_overlaps1.lobster_version = "3.1.0"
        self.band_overlaps2.lobster_version = "2.7.0"

        self.band_overlaps1.process()
        self.band_overlaps2.process()

        self.band_overlaps1_new = BandOverlaps(f"{TEST_DIR}/bandOverlaps.lobster.new.1")
        self.band_overlaps2_new = BandOverlaps(f"{TEST_DIR}/bandOverlaps.lobster.new.2")

    def test_attributes(self):
        # band_overlaps_dict
        bo_dict = self.band_overlaps1.band_overlaps
        assert bo_dict["max_deviations"][Spin.up][0] == approx(0.000278953)
        assert self.band_overlaps1_new.band_overlaps["max_deviations"][Spin.up][
            10
        ] == approx(0.0640933)
        assert bo_dict["matrices"][Spin.up][0].item(-1, -1) == approx(0.0188058)
        assert self.band_overlaps1_new.band_overlaps["matrices"][Spin.up][10].item(
            -1, -1
        ) == approx(1.0)
        assert bo_dict["matrices"][Spin.up][0].item(0, 0) == approx(1)
        assert self.band_overlaps1_new.band_overlaps["matrices"][Spin.up][10].item(
            0, 0
        ) == approx(0.995849)

        assert bo_dict["max_deviations"][Spin.down][-1] == approx(4.31567e-05)
        assert self.band_overlaps1_new.band_overlaps["max_deviations"][Spin.down][
            9
        ] == approx(0.064369)
        assert bo_dict["matrices"][Spin.down][-1].item(0, -1) == approx(4.0066e-07)
        assert self.band_overlaps1_new.band_overlaps["matrices"][Spin.down][9].item(
            0, -1
        ) == approx(1.37447e-09)

    def test_has_good_quality_maxDeviation(self):
        assert not self.band_overlaps1.has_good_quality_max_deviation(
            limit_max_deviation=0.1
        )
        assert not self.band_overlaps1_new.has_good_quality_max_deviation(
            limit_max_deviation=0.1
        )

        assert self.band_overlaps1.has_good_quality_max_deviation(
            limit_max_deviation=100
        )
        assert self.band_overlaps1_new.has_good_quality_max_deviation(
            limit_max_deviation=100
        )
        assert self.band_overlaps2.has_good_quality_max_deviation()
        assert not self.band_overlaps2_new.has_good_quality_max_deviation()
        assert not self.band_overlaps2.has_good_quality_max_deviation(
            limit_max_deviation=0.0000001
        )
        assert not self.band_overlaps2_new.has_good_quality_max_deviation(
            limit_max_deviation=0.0000001
        )

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
        assert self.band_overlaps2.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=10, limit_deviation=1
        )
        assert self.band_overlaps2_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=2, limit_deviation=0.1
        )
        assert self.band_overlaps2.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=1, limit_deviation=1
        )
        assert self.band_overlaps2_new.has_good_quality_check_occupied_bands(
            number_occ_bands_spin_up=1, limit_deviation=2
        )

    def test_has_good_quality_check_occupied_bands_patched(self):
        """Test with patched data."""

        limit_deviation = 0.1

        rng = np.random.default_rng(42)  # set seed for reproducibility

        band_overlaps = copy.deepcopy(self.band_overlaps1_new)

        number_occ_bands_spin_up_all = list(
            range(band_overlaps.band_overlaps["matrices"][Spin.up][0].shape[0])
        )
        number_occ_bands_spin_down_all = list(
            range(band_overlaps.band_overlaps["matrices"][Spin.down][0].shape[0])
        )

        for actual_deviation in [0.05, 0.1, 0.2, 0.5, 1.0]:
            for spin in (Spin.up, Spin.down):
                for number_occ_bands_spin_up, number_occ_bands_spin_down in zip(
                    number_occ_bands_spin_up_all,
                    number_occ_bands_spin_down_all,
                    strict=False,
                ):
                    for i_arr, array in enumerate(
                        band_overlaps.band_overlaps["matrices"][spin]
                    ):
                        number_occ_bands = (
                            number_occ_bands_spin_up
                            if spin is Spin.up
                            else number_occ_bands_spin_down
                        )

                        shape = array.shape
                        assert np.all(np.array(shape) >= number_occ_bands)
                        assert len(shape) == 2
                        assert shape[0] == shape[1]

                        # Generate a noisy background array
                        patch_array = rng.uniform(0, 10, shape)

                        # Patch the top-left sub-array (the part that would be checked)
                        patch_array[:number_occ_bands, :number_occ_bands] = np.identity(
                            number_occ_bands
                        ) + rng.uniform(
                            0, actual_deviation, (number_occ_bands, number_occ_bands)
                        )

                        band_overlaps.band_overlaps["matrices"][spin][
                            i_arr
                        ] = patch_array

                    result = band_overlaps.has_good_quality_check_occupied_bands(
                        number_occ_bands_spin_up=number_occ_bands_spin_up,
                        number_occ_bands_spin_down=number_occ_bands_spin_down,
                        spin_polarized=True,
                        limit_deviation=limit_deviation,
                    )

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
        with pytest.raises(
            ValueError, match="number_occ_bands_spin_down has to be specified"
        ):
            self.band_overlaps1.has_good_quality_check_occupied_bands(
                number_occ_bands_spin_up=4,
                spin_polarized=True,
            )
        with pytest.raises(
            ValueError, match="number_occ_bands_spin_down has to be specified"
        ):
            self.band_overlaps1_new.has_good_quality_check_occupied_bands(
                number_occ_bands_spin_up=4,
                spin_polarized=True,
            )

    def test_keys(self):
        bo_dict = self.band_overlaps1.band_overlaps
        bo_dict_new = self.band_overlaps1_new.band_overlaps
        bo_dict_2 = self.band_overlaps2.band_overlaps
        assert len(bo_dict["k_points"][Spin.up]) == 408
        assert len(bo_dict_2["max_deviations"][Spin.up]) == 2
        assert len(bo_dict_new["matrices"][Spin.down]) == 73


class TestGROSSPOP(MatSciTest):

    def setup_method(self):
        self.grosspop1 = GROSSPOP(f"{TEST_DIR}/GROSSPOP.lobster")
        self.grosspop_511_sp = GROSSPOP(f"{TEST_DIR}/GROSSPOP_511_sp.lobster.AlN.gz")
        self.grosspop_511_nsp = GROSSPOP(f"{TEST_DIR}/GROSSPOP_511_nsp.lobster.NaCl.gz")
        self.grosspop_511_lcfo = GROSSPOP_LCFO(
            f"{TEST_DIR}/GROSSPOP.LCFO.lobster.AlN.gz"
        )

    def test_attributes(self):
        gross_pop_list = self.grosspop1.populations
        gross_pop_list_511_sp = self.grosspop_511_sp.populations
        gross_pop_list_511_nsp = self.grosspop_511_nsp.populations
        gross_pop_list_lcfo = self.grosspop_511_lcfo.populations

        assert gross_pop_list["Si1"]["3p_y"][Spin.up]["mulliken"] == approx(0.38)
        assert gross_pop_list["Si1"]["3p_z"][Spin.up]["mulliken"] == approx(0.37)
        assert gross_pop_list["Si1"]["3p_x"][Spin.up]["mulliken"] == approx(0.37)
        assert gross_pop_list["Si1"]["3p_y"][Spin.up]["loewdin"] == approx(0.52)
        assert gross_pop_list["Si1"]["3p_z"][Spin.up]["loewdin"] == approx(0.52)
        assert gross_pop_list["Si1"]["3p_x"][Spin.up]["loewdin"] == approx(0.52)
        assert gross_pop_list["O5"]["2s"][Spin.up]["mulliken"] == approx(1.80)
        assert gross_pop_list["O5"]["2s"][Spin.up]["loewdin"] == approx(1.60)
        assert gross_pop_list["O8"]["2s"][Spin.up]["mulliken"] == approx(1.80)
        assert gross_pop_list["O8"]["2s"][Spin.up]["loewdin"] == approx(1.60)
        assert len(gross_pop_list) == 9

        # v5.1 spin polarized
        assert len(self.grosspop_511_sp.spins) == 2
        assert gross_pop_list_511_sp["Al1"]["3p_x"][Spin.up]["mulliken"] == approx(0.19)
        assert gross_pop_list_511_sp["N3"]["2s"][Spin.down]["loewdin"] == approx(0.7)

        # v5.1 non spin polarized
        assert len(self.grosspop_511_nsp.spins) == 1
        assert self.grosspop_511_lcfo.is_lcfo
        assert gross_pop_list_511_nsp["Na1"]["3s"][Spin.up]["mulliken"] == approx(0.22)
        assert gross_pop_list_511_nsp["Na1"]["3s"][Spin.up]["loewdin"] == approx(0.33)

        # v.5.1.1 LCFO
        assert self.grosspop_511_lcfo.is_lcfo
        assert gross_pop_list_lcfo["AlN1"]["1a1"][Spin.up]["loewdin"] == approx(0.81)
        assert gross_pop_list_lcfo["AlN1"]["1a1"][Spin.down]["loewdin"] == approx(0.81)

        with pytest.raises(KeyError):
            assert gross_pop_list_lcfo["AlN1"]["1a1"][Spin.up]["mulliken"]

    def test_msonable(self):
        dict_data = self.grosspop1.as_dict()
        grosspop_from_dict = GROSSPOP.from_dict(dict_data)
        all_attributes = vars(self.grosspop1)
        for attr_name, attr_value in all_attributes.items():
            assert getattr(grosspop_from_dict, attr_name) == attr_value


class TestICOXXLIST(MatSciTest):

    def setup_method(self):
        self.icohp_bise = ICOHPLIST(
            filename=f"{TEST_DIR}/ICOHPLIST.lobster.BiSe", process_immediately=False
        )
        self.icohp_bise.lobster_version = "3.1.0"
        self.icohp_bise.process()

        self.icoop_bise = ICOOPLIST(
            filename=f"{TEST_DIR}/ICOOPLIST.lobster.BiSe", process_immediately=False
        )
        self.icoop_bise.lobster_version = "3.2.0"
        self.icoop_bise.process()

        self.icohp_fe = ICOHPLIST(
            filename=f"{TEST_DIR}/ICOHPLIST.lobster", process_immediately=False
        )
        self.icohp_fe.lobster_version = "2.7.0"
        self.icohp_fe.process()

        self.icohp_gzipped = ICOHPLIST(
            filename=f"{TEST_DIR}/ICOHPLIST.lobster.gz", process_immediately=False
        )
        self.icohp_gzipped.lobster_version = "3.1.0"
        self.icohp_gzipped.process()

        self.icoop_fe = ICOOPLIST(
            filename=f"{TEST_DIR}/ICOOPLIST.lobster", process_immediately=False
        )
        self.icoop_fe.lobster_version = "2.7.0"
        self.icoop_fe.process()

        self.icohp_aln_511_sp = ICOHPLIST(
            filename=f"{TEST_DIR}/ICOHPLIST_511_sp.lobster.AlN.gz"
        )
        self.icohp_nacl_511_nsp = ICOHPLIST(
            filename=f"{TEST_DIR}/ICOHPLIST_511_nsp.lobster.NaCl.gz",
            process_immediately=False,
        )
        self.icohp_nacl_511_nsp.lobster_version = "5.0.5"
        self.icohp_nacl_511_nsp.process()

        # ICOHPLIST.LCFO.lobster from Lobster v5.1.1
        self.icohp_lcfo = ICOHPLIST_LCFO(
            filename=f"{TEST_DIR}/ICOHPLIST.LCFO.lobster.AlN.gz"
        )
        self.icohp_lcfo_non_orbitalwise = ICOHPLIST_LCFO(
            filename=f"{TEST_DIR}/ICOHPLIST_non_orbitalwise.LCFO.lobster.AlN.gz",
        )

        self.icobi_orbitalwise = ICOBILIST(
            filename=f"{TEST_DIR}/ICOBILIST.lobster", process_immediately=False
        )
        self.icobi_orbitalwise.lobster_version = "3.1.0"
        self.icobi_orbitalwise.process()

        self.icobi = ICOBILIST(
            filename=f"{TEST_DIR}/ICOBILIST.lobster.withoutorbitals",
            process_immediately=False,
        )
        self.icobi.lobster_version = "3.1.0"
        self.icobi.process()

        self.icobi_orbitalwise_spinpolarized = ICOBILIST(
            filename=f"{TEST_DIR}/ICOBILIST.lobster.spinpolarized",
            process_immediately=False,
        )
        self.icobi_orbitalwise_spinpolarized.lobster_version = "4.5.0"
        self.icobi_orbitalwise_spinpolarized.process()
        # make sure the correct line is read to check if this is a orbitalwise ICOBILIST
        self.icobi_orbitalwise_add = ICOBILIST(
            filename=f"{TEST_DIR}/ICOBILIST.lobster.additional_case",
            process_immediately=False,
        )
        self.icobi_orbitalwise_add.lobster_version = "3.1.0"
        self.icobi_orbitalwise_add.process()

        self.icobi_orbitalwise_spinpolarized_add = ICOBILIST(
            filename=f"{TEST_DIR}/ICOBILIST.lobster.spinpolarized.additional_case",
            process_immediately=False,
        )
        self.icobi_orbitalwise_spinpolarized_add.lobster_version = "4.5.0"
        self.icobi_orbitalwise_spinpolarized_add.process()

    def test_attributes(self):

        assert len(self.icohp_bise.spins) == 1
        assert len(self.icohp_bise.interactions) == 11
        assert len(self.icohp_fe.spins) == 2
        assert len(self.icohp_fe.interactions) == 2
        assert len(self.icoop_fe.spins) == 2
        assert len(self.icoop_fe.interactions) == 2

        # >v5 ICOHPLIST
        assert len(self.icohp_aln_511_sp.spins) == 2

        assert len(self.icohp_aln_511_sp.interactions) == 1088
        assert len(self.icohp_nacl_511_nsp.spins) == 1
        assert len(self.icohp_nacl_511_nsp.interactions) == 2584

        # v5.1.1 LCFO
        assert self.icohp_lcfo.is_lcfo
        assert len(self.icohp_lcfo.spins) == 2
        assert len(self.icohp_lcfo.interactions) == 1180
        assert len(self.icohp_lcfo_non_orbitalwise.interactions) == 28

    def test_values(self):
        icohplist_bise = [
            {
                "index": 1,
                "centers": ["Bi1", "Se7"],
                "length": 2.88231,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: -2.18042},
            },
            {
                "index": 2,
                "centers": ["Bi1", "Se10"],
                "length": 3.10144,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: -1.14347},
            },
            {
                "index": 3,
                "centers": ["Bi2", "Se8"],
                "length": 2.88231,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: -2.18042},
            },
            {
                "index": 4,
                "centers": ["Bi2", "Se9"],
                "length": 3.10144,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: -1.14348},
            },
            {
                "index": 5,
                "centers": ["Bi3", "Se10"],
                "length": 3.05001,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: -1.30006},
            },
            {
                "index": 6,
                "centers": ["Bi3", "Se11"],
                "length": 2.91676,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: -1.96843},
            },
            {
                "index": 7,
                "centers": ["Bi4", "Se9"],
                "length": 3.05001,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: -1.30006},
            },
            {
                "index": 8,
                "centers": ["Bi4", "Se12"],
                "length": 2.91676,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: -1.96843},
            },
            {
                "index": 9,
                "centers": ["Bi5", "Se12"],
                "length": 3.37522,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: -0.47531},
            },
            {
                "index": 10,
                "centers": ["Bi5", "Bi6"],
                "length": 3.07294,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: -2.38796},
            },
            {
                "index": 11,
                "centers": ["Bi6", "Se11"],
                "length": 3.37522,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: -0.47531},
            },
        ]
        icooplist_bise = [
            {
                "index": 1,
                "centers": ["Bi1", "Se7"],
                "length": 2.88231,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: 0.14245},
            },
            {
                "index": 2,
                "centers": ["Bi1", "Se10"],
                "length": 3.10144,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: -0.04118},
            },
            {
                "index": 3,
                "centers": ["Bi2", "Se8"],
                "length": 2.88231,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: 0.14245},
            },
            {
                "index": 4,
                "centers": ["Bi2", "Se9"],
                "length": 3.10144,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: -0.04118},
            },
            {
                "index": 5,
                "centers": ["Bi3", "Se10"],
                "length": 3.05001,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: -0.03516},
            },
            {
                "index": 6,
                "centers": ["Bi3", "Se11"],
                "length": 2.91676,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: 0.10745},
            },
            {
                "index": 7,
                "centers": ["Bi4", "Se9"],
                "length": 3.05001,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: -0.03516},
            },
            {
                "index": 8,
                "centers": ["Bi4", "Se12"],
                "length": 2.91676,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: 0.10745},
            },
            {
                "index": 9,
                "centers": ["Bi5", "Se12"],
                "length": 3.37522,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: -0.12395},
            },
            {
                "index": 10,
                "centers": ["Bi5", "Bi6"],
                "length": 3.07294,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: 0.24714},
            },
            {
                "index": 11,
                "centers": ["Bi6", "Se11"],
                "length": 3.37522,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: -0.12395},
            },
        ]
        icooplist_fe = [
            {
                "index": 1,
                "centers": ["Fe8", "Fe7"],
                "length": 2.83189,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: -0.11389, Spin.down: -0.20828},
            },
            {
                "index": 2,
                "centers": ["Fe8", "Fe9"],
                "length": 2.45249,
                "cells": [[], []],
                "orbitals": [None, None],
                "icoxx": {Spin.up: -0.04087, Spin.down: -0.05756},
            },
        ]

        assert icohplist_bise == self.icohp_bise.interactions
        assert icooplist_fe == self.icoop_fe.interactions
        assert icooplist_bise == self.icoop_bise.interactions

        assert self.icobi.interactions[1]["icoxx"][Spin.up] == approx(0.58649)
        assert self.icobi_orbitalwise.interactions[2]["icoxx"][Spin.up] == approx(
            0.02559
        )
        assert self.icobi_orbitalwise.interactions[1]["icoxx"][Spin.up] == approx(
            0.04940
        )
        assert self.icobi_orbitalwise_spinpolarized.interactions[1]["icoxx"][
            Spin.up
        ] == approx(0.04940 / 2, abs=1e-3)
        assert self.icobi_orbitalwise_spinpolarized.interactions[1]["icoxx"][
            Spin.down
        ] == approx(0.04940 / 2, abs=1e-3)
        assert self.icobi_orbitalwise_spinpolarized.interactions[2]["icoxx"][
            Spin.down
        ] == approx(0.01279, abs=1e-3)
        assert self.icobi_orbitalwise_spinpolarized.interactions[2]["orbitals"] == [
            "2p_y",
            "6s",
        ]

        # >v5 ICOHPLIST
        assert self.icohp_aln_511_sp.interactions[2]["icoxx"][Spin.up] == approx(
            0.00102
        )
        assert self.icohp_aln_511_sp.interactions[2]["icoxx"][Spin.down] == approx(
            0.00104
        )
        assert self.icohp_nacl_511_nsp.interactions[13]["icoxx"][Spin.up] == approx(0.0)
        assert self.icohp_nacl_511_nsp.interactions[10]["orbitals"] == ["2p_y", "2p_z"]

        # v5.1.1 LCFO
        assert self.icohp_lcfo.interactions[15]["orbitals"] == ["2a1", "4e"]
        assert self.icohp_lcfo_non_orbitalwise.interactions[16]["icoxx"][
            Spin.up
        ] == approx(-0.21495)
        assert self.icohp_lcfo_non_orbitalwise.interactions[16]["icoxx"][
            Spin.down
        ] == approx(-0.21498)

    def test_msonable(self):
        dict_data = self.icobi_orbitalwise_spinpolarized.as_dict()
        icohplist_from_dict = ICOHPLIST.from_dict(dict_data)
        all_attributes = vars(self.icobi_orbitalwise_spinpolarized)
        for attr_name, attr_value in all_attributes.items():
            if isinstance(attr_value, np.ndarray):
                assert_array_equal(getattr(icohplist_from_dict, attr_name), attr_value)
            else:
                assert getattr(icohplist_from_dict, attr_name) == attr_value


class TestWavefunction(MatSciTest):
    def test_parse_file(self):
        wf = Wavefunction(
            filename=f"{TEST_DIR}/LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR_O.gz"),
        )

        assert_array_equal([41, 41, 41], wf.grid)
        assert wf.points[4][0] == approx(0.0000)
        assert wf.points[4][1] == approx(0.0000)
        assert wf.points[4][2] == approx(0.4000)
        assert wf.reals[8] == approx(1.38863e-01)
        assert wf.imaginaries[8] == approx(2.89645e-01)
        assert len(wf.imaginaries) == 41 * 41 * 41
        assert len(wf.reals) == 41 * 41 * 41
        assert len(wf.points) == 41 * 41 * 41
        assert wf.distances[0] == approx(0.0000)

    def test_set_volumetric_data(self):
        wave1 = Wavefunction(
            filename=f"{TEST_DIR}/LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR_O.gz"),
        )

        wave1.set_volumetric_data(grid=wave1.grid, structure=wave1.structure)
        assert wave1.volumetricdata_real.data["total"][0, 0, 0] == approx(-3.0966)
        assert wave1.volumetricdata_imaginary.data["total"][0, 0, 0] == approx(
            -6.45895e00
        )

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
        assert volumetricdata_density.data["total"][0, 0, 0] == approx(
            (-3.0966 * -3.0966) + (-6.45895 * -6.45895)
        )

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
        self.sitepotentials = SitePotentials(
            filename=f"{TEST_DIR}/SitePotentials.lobster.perovskite"
        )

    def test_attributes(self):
        assert self.sitepotentials.site_potentials_loewdin == [
            -8.77,
            -17.08,
            9.57,
            9.57,
            8.45,
        ]
        assert self.sitepotentials.site_potentials_mulliken == [
            -11.38,
            -19.62,
            11.18,
            11.18,
            10.09,
        ]
        assert self.sitepotentials.madelung_energies_loewdin == approx(-28.64)
        assert self.sitepotentials.madelung_energies_mulliken == approx(-40.02)
        assert self.sitepotentials.centers == ["La1", "Ta2", "N3", "N4", "O5"]
        assert len(self.sitepotentials.centers) == 5
        assert self.sitepotentials.ewald_splitting == approx(3.14)

    def test_msonable(self):
        dict_data = self.sitepotentials.as_dict()
        sitepotential_from_dict = SitePotentials.from_dict(dict_data)
        all_attributes = vars(self.sitepotentials)
        for attr_name, attr_value in all_attributes.items():
            assert getattr(sitepotential_from_dict, attr_name) == attr_value


class TestMadelungEnergies(MatSciTest):
    def setup_method(self) -> None:
        self.madelungenergies = MadelungEnergies(
            filename=f"{TEST_DIR}/MadelungEnergies.lobster.perovskite"
        )

    def test_attributes(self):
        assert self.madelungenergies.madelung_energies_loewdin == approx(-28.64)
        assert self.madelungenergies.madelung_energies_mulliken == approx(-40.02)
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
            filename=f"{TEST_DIR}/Na_hamiltonMatrices.lobster.gz", efermi=-2.79650354
        )
        self.transfer_matrices = LobsterMatrices(
            filename=f"{TEST_DIR}/C_transferMatrices.lobster.gz"
        )
        self.overlap_matrices = LobsterMatrices(
            filename=f"{TEST_DIR}/Si_overlapMatrices.lobster.gz"
        )
        self.coeff_matrices = LobsterMatrices(
            filename=f"{TEST_DIR}/Si_coefficientMatricesLSO1.lobster.gz"
        )

    def test_attributes(self):
        assert self.hamilton_matrices.efermi == -2.79650354
        assert self.hamilton_matrices.centers == ["Na1", "Na1", "Na1", "Na1"]

        assert 1 in self.hamilton_matrices.matrices
        assert len(self.hamilton_matrices.matrices) == 1

        assert isinstance(self.hamilton_matrices.matrices[1], dict)
        for spin in self.hamilton_matrices.matrices[1]:
            assert spin in [Spin.up, Spin.down]
            assert isinstance(self.hamilton_matrices.matrices[1][spin], np.ndarray)
            assert self.hamilton_matrices.matrices[1][spin].shape == (4, 4)

        assert self.hamilton_matrices.orbitals == ["3s", "2p_y", "2p_z", "2p_x"]

        with pytest.raises(KeyError):
            assert self.hamilton_matrices.matrices[0][Spin.down]

        assert self.hamilton_matrices.matrices[1][Spin.up][0, 0].real == approx(
            -3.02170000
        )
        assert self.hamilton_matrices.matrices[1][Spin.up][0, 0].imag == approx(0.0)

        assert self.hamilton_matrices.get_onsite_values("Na1", "3s") == approx(
            (-3.0217 + 2.79650354 - 1.39420000 + 2.79650354) / 2
        )

        assert self.hamilton_matrices.get_onsite_values("Na1", "2p_x") == approx(
            (-28.56640000 + 2.79650354 - 28.48100000 + 2.79650354) / 2
        )

        onsite_values = self.hamilton_matrices.get_onsite_values()
        assert isinstance(onsite_values, dict)

        for key in onsite_values:
            assert key in ["Na1_3s", "Na1_2p_y", "Na1_2p_z", "Na1_2p_x"]
            assert isinstance(onsite_values[key], float)

        assert self.overlap_matrices.efermi is None
        assert self.overlap_matrices.centers == ["Si1", "Si1", "Si1", "Si1"]
        assert self.overlap_matrices.orbitals == ["3s", "3p_y", "3p_z", "3p_x"]

        assert 1 in self.overlap_matrices.matrices
        assert len(self.overlap_matrices.matrices) == 1

        assert isinstance(self.overlap_matrices.matrices[1], dict)
        for spin in self.overlap_matrices.matrices[1]:
            assert spin in [None]
            assert isinstance(self.overlap_matrices.matrices[1][spin], np.ndarray)
            assert self.overlap_matrices.matrices[1][spin].shape == (4, 4)

        for m in range(4):
            assert self.overlap_matrices.matrices[1][None][m, m].real == approx(
                1.00000000
            )

        assert self.transfer_matrices.efermi is None
        assert self.transfer_matrices.centers == ["C1", "C1", "C1", "C1"]
        assert self.transfer_matrices.orbitals == ["2s", "2p_y", "2p_z", "2p_x"]

        assert isinstance(self.transfer_matrices.matrices[1], dict)

        assert 1 in self.transfer_matrices.matrices
        assert len(self.transfer_matrices.matrices) == 1

        for spin in self.transfer_matrices.matrices[1]:
            assert spin in [Spin.up, Spin.down]
            assert isinstance(self.transfer_matrices.matrices[1][spin], np.ndarray)
            assert self.transfer_matrices.matrices[1][spin].shape == (4, 4)

        assert self.coeff_matrices.efermi is None
        assert self.coeff_matrices.centers == ["Si1", "Si1", "Si1", "Si1"]
        assert self.coeff_matrices.orbitals == ["3s", "3p_y", "3p_z", "3p_x"]

        assert isinstance(self.coeff_matrices.matrices[1], dict)
        assert 1 in self.coeff_matrices.matrices
        assert len(self.coeff_matrices.matrices) == 1
        for spin in self.coeff_matrices.matrices[1]:
            assert spin in [Spin.up, Spin.down]
            assert isinstance(self.coeff_matrices.matrices[1][spin], np.ndarray)
            assert self.coeff_matrices.matrices[1][spin].shape == (4, 4)


class TestPOLARIZATION(MatSciTest):
    def setup_method(self) -> None:
        self.polarization = POLARIZATION(
            filename=f"{TEST_DIR}/POLARIZATION.lobster.AlN.gz"
        )

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


class TestCOBICAR(MatSciTest):
    """Tests for COBICAR class."""

    def test_read_cobicar_spin(self):
        cobicar = COBICAR(TEST_DIR + "/COBICAR.lobster.B2H6.spin")

        assert cobicar.is_spin_polarized
        assert cobicar.data.shape == (
            cobicar.num_data,
            cobicar.num_bonds * 2 * (len(cobicar.spins)) + 1,
        )

        for interaction in cobicar.interactions:
            assert "coxx" in interaction
            assert "icoxx" in interaction

            assert len(interaction["coxx"]) == 2

            assert Spin.up in interaction["icoxx"]
            assert Spin.down in interaction["icoxx"]

        assert len(cobicar.energies) == cobicar.num_data

        assert (
            len(cobicar.get_interactions_by_properties(centers=["H4", "B1", "H7"])) == 5
        )
        assert (
            len(
                cobicar.get_interactions_by_properties(
                    centers=["H4", "B1", "H7"], orbitals=["2s"]
                )
            )
            == 1
        )

        assert cobicar.get_data_by_properties(
            centers=["H4", "B1", "H7"],
            orbitals=["2s"],
            spins=[Spin.up, Spin.down],
        ).shape == (cobicar.num_data, 4)

        assert cobicar.get_data_by_properties(
            centers=["H4", "B1", "H7"],
            orbitals=["2s"],
            spins=[Spin.down],
            data_type="icoxx",
        ).shape == (cobicar.num_data, 1)

    def test_read_cobicar_4_centers(self):
        cobicar = COBICAR(TEST_DIR + "/COBICAR.lobster.GeTe_4center")

        assert (
            len(cobicar.get_interactions_by_properties(centers=["Ge", "Ge", "T", "Te"]))
            == 1
        )

    def test_read_cobicar_4_centers_orbital_resolved(self, tmp_path):
        cobicar = COBICAR(TEST_DIR + "/COBICAR.lobster.GeTe.multi.orbitalwise")

        assert (
            len(
                cobicar.get_interactions_by_properties(
                    centers=["Ge", "Ge", "T", "Te"],
                    orbitals=["5p_z", "4p_z", "5p_z", "p_x"],
                )
            )
            == 2
        )

        assert len(cobicar.get_interactions_by_properties(indices=[13])) == 257
        assert len(cobicar.get_interactions_by_properties(orbitals=["2p_x"])) == 0
        assert len(cobicar.get_interactions_by_properties(cells=[[1, 0, 0]])) == 257

        interactions = cobicar.get_interaction_indices_by_properties(
            centers=["Ge", "Ge", "T", "Te"],
            orbitals=["5p_z", "4p_z", "5p_z", "p_x"],
        )

        assert len(interactions) == 2
        assert cobicar.interactions[interactions[0]]["centers"] == [
            "Te8",
            "Ge1",
            "Te8",
            "Ge1",
        ]
        assert cobicar.interactions[interactions[0]]["orbitals"] == [
            "5p_z",
            "4p_z",
            "5p_z",
            "4p_x",
        ]

        data_indices = cobicar.interaction_indices_to_data_indices_mapping(interactions)

        assert len(data_indices) == 4
        assert data_indices[0] == interactions[0] * 2 + 1
        assert data_indices[1] == interactions[0] * 2 + 2
        assert data_indices[2] == interactions[1] * 2 + 1
        assert data_indices[3] == interactions[1] * 2 + 2

        cobicar.save(f"{tmp_path}/cobicar.json")
        cobicar_from_json = COBICAR.load(f"{tmp_path}/cobicar.json")

        assert cobicar_from_json.is_spin_polarized == cobicar.is_spin_polarized
        assert cobicar_from_json.num_data == cobicar.num_data
        assert cobicar_from_json.num_bonds == cobicar.num_bonds
        assert len(cobicar_from_json.interactions) == len(cobicar.interactions)

        for interaction1, interaction2 in zip(
            cobicar_from_json.interactions, cobicar.interactions, strict=True
        ):
            assert interaction1["centers"] == interaction2["centers"]
            assert interaction1["orbitals"] == interaction2["orbitals"]
            assert_allclose(
                interaction1["coxx"][Spin.up], interaction2["coxx"][Spin.up]
            )
            assert_allclose(
                interaction1["icoxx"][Spin.up], interaction2["icoxx"][Spin.up]
            )

        assert_allclose(cobicar_from_json.data, cobicar.data)

        ya_cobicar = COBICAR.from_dict(cobicar.as_dict())

        for interaction1, interaction2 in zip(
            ya_cobicar.interactions, cobicar_from_json.interactions, strict=True
        ):
            assert interaction1["centers"] == interaction2["centers"]
            assert interaction1["orbitals"] == interaction2["orbitals"]
            assert_allclose(
                interaction1["coxx"][Spin.up], interaction2["coxx"][Spin.up]
            )
            assert_allclose(
                interaction1["icoxx"][Spin.up], interaction2["icoxx"][Spin.up]
            )

        assert_allclose(ya_cobicar.data, cobicar_from_json.data)


class TestCOHPCAR(MatSciTest):
    """Tests for COHPCAR class."""

    def test_read_cobicar_lcfo(self):
        cohpcar = COBICAR_LCFO(TEST_DIR + "/COHPCAR.LCFO.lobster.NaCl.gz")

        assert cohpcar.is_lcfo

        assert cohpcar.interactions[2]["centers"] == ["NaCl_1", "NaCl_1"]
        assert cohpcar.interactions[2]["orbitals"] == ["1a1", "1a1"]
        assert cohpcar.interactions[2]["length"] == approx(2.8473125412856759)

        assert cohpcar.data.shape == (cohpcar.num_data, cohpcar.num_bonds * 2 + 1)

        assert cohpcar.interactions[-1]["centers"] == ["Na2_2", "Cl_6"]
        assert cohpcar.interactions[-1]["orbitals"] == ["2a1u", "3p_x"]

        assert (
            len(
                cohpcar.get_interactions_by_properties(
                    centers=["NaCl_1", "Na2_2"],
                    orbitals=["a1", "a1"],
                )
            )
            == 48
        )


class TestCOOPCAR(MatSciTest):
    """Tests for COOPCAR class."""

    def test_coopcar(self):
        coopcar = COOPCAR(TEST_DIR + "/COOPCAR.lobster.gz")

        assert coopcar.is_spin_polarized
        assert coopcar.data.shape == (
            coopcar.num_data,
            coopcar.num_bonds * 2 * (len(coopcar.spins)) + 1,
        )

        for interaction in coopcar.interactions:
            assert "coxx" in interaction
            assert "icoxx" in interaction

            assert len(interaction["coxx"]) == 2

            assert Spin.up in interaction["icoxx"]
            assert Spin.down in interaction["icoxx"]

        assert len(coopcar.energies) == coopcar.num_data

        assert len(coopcar.get_interactions_by_properties(centers=["Fe8", "Fe7"])) == 1
        assert coopcar.get_data_by_properties(centers=["Fe8", "Fe9"])[0, -1] == approx(
            -0.00099
        )

    def test_coopcar_2(self):
        coopcar = COOPCAR(TEST_DIR + "/COOPCAR.lobster.BiSe.gz")

        assert not coopcar.is_spin_polarized
        assert coopcar.data.shape == (coopcar.num_data, coopcar.num_bonds * 2 + 1)

    def test_coopcar_3(self):
        coopcar = COOPCAR(TEST_DIR + "/COOPCAR.lobster.KF.gz")

        assert not coopcar.is_spin_polarized
        assert coopcar.energies.shape == (coopcar.num_data,)

        assert coopcar.data.shape == (coopcar.num_data, coopcar.num_bonds * 2 + 1)

        assert coopcar.interactions[0]["centers"] == ["Average"]
        assert coopcar.interactions[1]["centers"] == ["F1", "K2"]


class TestNcICOBILIST(MatSciTest):
    """Tests for NcICOBILIST class."""

    def test_ncicobilist(self):
        ncicobi = NcICOBILIST(
            filename=f"{TEST_DIR}/NcICOBILIST.lobster.nospin.withoutorbitals"
        )

        assert len(ncicobi.spins) == 1
        assert ncicobi.interactions[0]["centers"] == ["X1", "X20"]
        assert ncicobi.interactions[0]["orbitals"] == [None, None]
        assert ncicobi.interactions[0]["icoxx"][Spin.up] == approx(0)

        with pytest.raises(KeyError):
            ncicobi.interactions[0]["icoxx"][Spin.down]

        assert ncicobi.data.shape == (2, 1)

        assert ncicobi.get_interactions_by_properties(
            centers=["X22"],
        )[
            0
        ]["icoxx"][
            Spin.up
        ] == approx(0.00018)

        assert ncicobi.get_data_by_properties(
            centers=["X22"],
            spins=[Spin.up],
        ) == approx(0.00018)

    def test_ncicobilist_spin(self):
        ncicobi = NcICOBILIST(filename=f"{TEST_DIR}/NcICOBILIST.lobster")

        assert len(ncicobi.spins) == 2

        assert ncicobi.interactions[0]["centers"] == ["X1", "X20"]
        assert ncicobi.interactions[0]["orbitals"] == [None, None]
        assert ncicobi.interactions[1]["icoxx"][Spin.up] == approx(0.00009)
        assert ncicobi.interactions[1]["icoxx"][Spin.down] == approx(0.00009)

        assert ncicobi.data.shape == (24, 2)

        interaction = ncicobi.get_interactions_by_properties(
            orbitals=["4d_x^2-y^2", "4d_x^2-y^2"]
        )

        assert len(interaction) == 1
        assert interaction[0]["index"] == 2
        assert interaction[0]["centers"] == ["X22", "Xs42", "X31"]
        assert interaction[0]["orbitals"] == ["4d_x^2-y^2", "3p", "4d_x^2-y^2"]

        for interaction in ncicobi.interactions:
            assert "icoxx" in interaction
            assert Spin.up in interaction["icoxx"]
            assert Spin.down in interaction["icoxx"]

            assert "centers" in interaction
            assert "orbitals" in interaction
            assert "length" in interaction


class TestMsonable(MatSciTest):

    def setup_method(self) -> None:
        self.objects_to_test: dict[type[MSONable], str] = {
            BWDF: f"{TEST_DIR}/BWDF.lobster.AlN.gz",
            CHARGE: f"{TEST_DIR}/CHARGE.lobster.MnO2.gz",
            CHARGE_LCFO: f"{TEST_DIR}/CHARGE.LCFO.lobster.ALN.gz",
            COBICAR: f"{TEST_DIR}/COBICAR.lobster.B2H6.spin",
            COBICAR_LCFO: f"{TEST_DIR}/COHPCAR.LCFO.lobster.NaCl.gz",
            COOPCAR: f"{TEST_DIR}/COOPCAR.lobster.gz",
            GROSSPOP: f"{TEST_DIR}/GROSSPOP.lobster",
            GROSSPOP_LCFO: f"{TEST_DIR}/GROSSPOP.LCFO.lobster.AlN.gz",
            ICOHPLIST: f"{TEST_DIR}/ICOHPLIST_511_sp.lobster.AlN.gz",
            POLARIZATION: f"{TEST_DIR}/POLARIZATION.lobster.AlN.gz",
            BandOverlaps: f"{TEST_DIR}/bandOverlaps.lobster.new.1",
            Fatband: f"{TEST_DIR}/Fatband_SiO2/Test_p/FATBAND_o4_2p.lobster",
            LobsterOut: f"{TEST_DIR}/lobsterout.normal",
            MadelungEnergies: f"{TEST_DIR}/MadelungEnergies.lobster.perovskite",
            NcICOBILIST: f"{TEST_DIR}/NcICOBILIST.lobster.nospin.withoutorbitals",
            SitePotentials: f"{TEST_DIR}/SitePotentials.lobster.perovskite",
        }
        self.instances_to_test = [
            DOSCAR(
                filename=f"{VASP_OUT_DIR}/DOSCAR.lobster.spin",
            ),
            DOSCAR_LCFO(
                filename=f"{VASP_OUT_DIR}/DOSCAR.LCFO.lobster.AlN",
            ),
            LobsterMatrices(
                filename=f"{TEST_DIR}/Na_hamiltonMatrices.lobster.gz",
                efermi=-2.79650354,
            ),
            Wavefunction(
                filename=f"{TEST_DIR}/LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz",
                structure=Structure.from_file(f"{TEST_DIR}/POSCAR_O.gz"),
            ),
        ]

    def test_json_save_load(self, tmp_path):
        """Tests saving and loading of all MSONable classes in this package."""

        def check_msonability(instance: MSONable) -> None:
            instance.save(f"{tmp_path}/{instance.__class__.__name__.lower()}.json")
            instance_from_json = instance.load(
                f"{tmp_path}/{instance.__class__.__name__.lower()}.json"
            )

            json1 = json.dumps(instance.as_dict(), cls=MontyEncoder, sort_keys=True)
            json2 = json.dumps(
                instance_from_json.as_dict(), cls=MontyEncoder, sort_keys=True
            )
            assert json1 == json2

            for attr_name in vars(instance):
                assert hasattr(instance_from_json, attr_name)

        for obj, file_path in self.objects_to_test.items():
            instance: MSONable = obj(filename=file_path)

            check_msonability(instance)

        for instance in self.instances_to_test:
            check_msonability(instance)
