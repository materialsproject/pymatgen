from __future__ import annotations

import gzip
import json
import os
import sys
import xml
from io import StringIO
from pathlib import Path
from shutil import copyfile, copyfileobj

import numpy as np
import pytest
from monty.io import zopen
from monty.shutil import decompress_file
from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.core import Element
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.bandstructure import BandStructure, BandStructureSymmLine
from pymatgen.electronic_structure.core import Magmom, Orbital, OrbitalType, Spin
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, Potcar
from pymatgen.io.vasp.outputs import (
    WSWQ,
    BandgapProps,
    BSVasprun,
    Chgcar,
    Dynmat,
    Eigenval,
    Elfcar,
    KpointOptProps,
    Locpot,
    Oszicar,
    Outcar,
    Procar,
    UnconvergedVASPWarning,
    Vaspout,
    VaspParseError,
    Vasprun,
    Wavecar,
    Waveder,
    Xdatcar,
    get_band_structure_from_vasp_multiple_branches,
)
from pymatgen.io.wannier90 import Unk
from pymatgen.util.testing import FAKE_POTCAR_DIR, TEST_FILES_DIR, VASP_IN_DIR, VASP_OUT_DIR, MatSciTest

try:
    import h5py
except ImportError:
    h5py = None

TEST_DIR = f"{TEST_FILES_DIR}/io/vasp"


class TestVasprun(MatSciTest):
    def test_vasprun_soc(self):
        # Test that SOC vaspruns are parsed appropriately, giving just Spin.Up tdos, idos and pdos
        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.int_Te_SOC.xml.gz")
        dos_density_dicts_to_check = [vasp_run.complete_dos.densities, vasp_run.tdos.densities, vasp_run.idos.densities]
        dos_density_dicts_to_check += [
            densities for orbital_dict in vasp_run.complete_dos.pdos.values() for densities in orbital_dict.values()
        ]
        for i, dos_density_dict in enumerate(dos_density_dicts_to_check):
            assert set(dos_density_dict.keys()) == {Spin.up}, f"Failed spin keys check for {i}th dos obj!"

        assert vasp_run.complete_dos.spin_polarization is None

    def test_vasprun_ml(self):
        # Test for ML MD simulation
        # The trajectory data is stored in md_data
        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.ml_md.xml.gz")
        assert len(vasp_run.md_data) == 100
        for frame in vasp_run.md_data:
            assert "structure" in frame
            assert "forces" in frame
            assert "energy" in frame
        assert vasp_run.md_data[-1]["energy"]["total"] == approx(-491.51831988)
        assert vasp_run.md_n_steps == 100
        assert vasp_run.converged_ionic

    def test_vasprun_md(self):
        # Test for simple MD simulation (no ML).
        # Does not generate the `md_data` attribute in Vasprun. Data based on `ionic_steps`
        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.md.xml.gz")
        assert len(vasp_run.ionic_steps) == 10
        assert vasp_run.final_energy == approx(-327.73014059)
        assert vasp_run.md_n_steps == 10
        assert vasp_run.converged_ionic

    def test_vasprun_ediffg_set_to_0(self):
        # Test for case where EDIFFG is set to 0. This should pass if all ionic steps
        # complete and are electronically converged.
        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.ediffg_set_to_0.xml.gz")
        assert len(vasp_run.ionic_steps) == 3
        assert vasp_run.final_energy == approx(-34.60164204)
        assert vasp_run.converged_ionic is True
        assert vasp_run.converged_electronic is True
        assert vasp_run.converged is True
        assert vasp_run.parameters["EDIFFG"] == 0
        assert vasp_run.parameters["EDIFF"] == approx(1e-5)

    def test_bad_random_seed(self):
        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.bad_random_seed.xml.gz")
        assert vasp_run.incar["ISMEAR"] == 0
        assert vasp_run.incar["RANDOM_SEED"] is None

    def test_multiple_dielectric(self):
        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.GW0.xml.gz")
        assert len(vasp_run.dielectric_data) == 4
        assert "HEAD OF MICROSCOPIC DIELECTRIC TENSOR (INDEPENDENT PARTICLE)" in vasp_run.dielectric_data

    def test_charge_charge_dielectric(self):
        """
        VASP 5.4.4 writes out two dielectric functions to vasprun.xml
        These are the "density-density" and "velocity-velocity" linear response functions.
        See the comments in `linear_optics.F` for details.
        """
        vasp_run = Vasprun(
            f"{VASP_OUT_DIR}/vasprun.dielectric_5.4.4.xml.gz",
            parse_potcar_file=False,
        )
        assert vasp_run.dielectric is not None
        assert "density" in vasp_run.dielectric_data
        assert "velocity" in vasp_run.dielectric_data

    def test_BSE(self):
        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.BSE.xml.gz")
        absorption_coeff = vasp_run.optical_absorption_coeff
        assert absorption_coeff[1] == approx(0.8327903762077188)
        assert vasp_run.final_structure == vasp_run.initial_structure
        assert "freq_dependent" in vasp_run.dielectric_data

    def test_vasprun_with_more_than_two_unlabelled_dielectric_functions(self):
        with pytest.warns(
            UserWarning,
            match="Additional unlabelled dielectric data in vasprun.xml are stored as unlabelled.",
        ):
            vr = Vasprun(f"{VASP_OUT_DIR}/vasprun.dielectric_bad.xml.gz")
        assert "unlabelled" in vr.dielectric_data

    def test_bad_vasprun(self):
        with pytest.raises(xml.etree.ElementTree.ParseError):
            Vasprun(f"{VASP_OUT_DIR}/vasprun.bad.xml.gz")

        with pytest.warns(
            UserWarning,
            match="XML is malformed. Parsing has stopped but partial data is available",
        ):
            vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.bad.xml.gz", exception_on_bad_xml=False)
        assert len(vasp_run.ionic_steps) == 1
        assert vasp_run.final_energy == approx(-269.00551374)

    @pytest.mark.filterwarnings("ignore::pymatgen.io.vasp.outputs.UnconvergedVASPWarning")
    def test_runtype(self):
        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.GW0.xml.gz")
        assert vasp_run.run_type in "HF"

        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.pbesol_vdw.xml.gz")
        assert vasp_run.run_type in "PBEsol+vdW-DFT-D3-BJ"

        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.hse06.xml.gz")
        assert vasp_run.run_type in "HSE06"

        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.scan_rvv10.xml.gz")
        assert vasp_run.run_type in "SCAN+rVV10"

        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.dfpt.ionic.xml.gz")
        assert vasp_run.run_type in "GGA"

        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.dfpt.xml.gz")
        assert vasp_run.run_type in "GGA+U"

        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.r2scan.xml.gz")
        assert vasp_run.run_type in "R2SCAN"

        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.scan.xml.gz")
        assert vasp_run.run_type in "SCAN"

        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.pbesol.xml.gz")
        assert vasp_run.run_type in "PBEsol"

        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.rscan.xml.gz")
        assert vasp_run.run_type in "RSCAN"

        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.random.xml.gz")
        assert vasp_run.run_type in "RANDOMFUNCTIONAL"

        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.unknown.xml.gz")
        with pytest.warns(UserWarning, match="Unknown run type!"):
            assert vasp_run.run_type in "unknown"

    def test_vdw(self):
        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.vdw.xml.gz")
        assert vasp_run.final_energy == approx(-9.78310677)

    def test_energies(self):
        # VASP 5.4.1
        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.etest1.xml.gz")
        assert vasp_run.final_energy == approx(-11.18981538)

        # VASP 6.2.1
        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.etest2.xml.gz")
        assert vasp_run.final_energy == approx(-11.18986774)

        # VASP 5.4.1
        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.etest3.xml.gz")
        assert vasp_run.final_energy == approx(-15.89355325)

        # VASP 6.2.1
        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.etest4.xml.gz")
        assert vasp_run.final_energy == approx(-15.89364691)

    def test_nonlmn(self):
        filepath = f"{VASP_OUT_DIR}/vasprun.nonlm.xml.gz"
        vasp_run = Vasprun(filepath, parse_potcar_file=False)
        orbs = list(vasp_run.complete_dos.pdos[vasp_run.final_structure[0]])
        assert OrbitalType.s in orbs

    @pytest.mark.filterwarnings("ignore:invalid value encountered in divide")
    def test_standard(self):
        filepath = f"{VASP_OUT_DIR}/vasprun.xml.gz"
        vasp_run = Vasprun(filepath, parse_potcar_file=False)

        # Test NELM parsing
        assert vasp_run.parameters["NELM"] == 60

        # test pDOS parsing
        assert vasp_run.complete_dos.spin_polarization == approx(1.0)
        assert Vasprun(f"{VASP_OUT_DIR}/vasprun.etest1.xml.gz").complete_dos.spin_polarization is None

        pdos0 = vasp_run.complete_dos.pdos[vasp_run.final_structure[0]]
        assert pdos0[Orbital.s][Spin.up][16] == approx(0.0026)
        assert pdos0[Orbital.pz][Spin.down][16] == approx(0.0012)
        assert pdos0[Orbital.s][Spin.up].shape == (301,)

        pdos0_norm = vasp_run.complete_dos_normalized.pdos[vasp_run.final_structure[0]]
        assert pdos0_norm[Orbital.s][Spin.up][16] == approx(0.0026)  # the site data should not change
        assert pdos0_norm[Orbital.s][Spin.up].shape == (301,)

        cdos_norm, cdos = vasp_run.complete_dos_normalized, vasp_run.complete_dos
        ratio = np.nanmax(cdos.densities[Spin.up] / cdos_norm.densities[Spin.up])
        assert ratio == approx(vasp_run.final_structure.volume)  # the site data should not change

        # check you can normalize an existing DOS
        cdos_norm2 = cdos.get_normalized()
        ratio = np.nanmax(cdos.densities[Spin.up] / cdos_norm2.densities[Spin.up])
        assert ratio == approx(vasp_run.final_structure.volume)  # the site data should not change

        # but doing so twice should not change the data
        cdos_norm3 = cdos_norm2.get_normalized()
        ratio = np.nanmax(cdos.densities[Spin.up] / cdos_norm3.densities[Spin.up])
        assert ratio == approx(vasp_run.final_structure.volume)  # the site data should not change

        pdos0_norm = vasp_run.complete_dos_normalized.pdos[vasp_run.final_structure[0]]
        assert pdos0_norm[Orbital.s][Spin.up][16] == approx(0.0026)  # the site data should not change
        assert pdos0_norm[Orbital.s][Spin.up].shape == (301,)

        cdos_norm, cdos = vasp_run.complete_dos_normalized, vasp_run.complete_dos
        ratio = np.nanmax(cdos.densities[Spin.up] / cdos_norm.densities[Spin.up])
        assert ratio == approx(vasp_run.final_structure.volume)  # the site data should not change

        filepath2 = f"{VASP_OUT_DIR}/vasprun.lifepo4.xml.gz"
        vasprun_ggau = Vasprun(filepath2, parse_projected_eigen=True, parse_potcar_file=False)
        total_sc_steps = sum(len(i["electronic_steps"]) for i in vasp_run.ionic_steps)
        assert len(vasp_run.ionic_steps) == 29
        assert len(vasp_run.structures) == len(vasp_run.ionic_steps)

        trajectory = vasp_run.get_trajectory()
        assert len(trajectory) == len(vasp_run.ionic_steps)
        assert "forces" in trajectory[0].site_properties

        for idx, step in enumerate(vasp_run.ionic_steps):
            assert vasp_run.structures[idx] == step["structure"]

        assert all(
            vasp_run.structures[idx] == vasp_run.ionic_steps[idx]["structure"]
            for idx in range(len(vasp_run.ionic_steps))
        )

        assert total_sc_steps == 308, "Incorrect number of energies read from vasprun.xml"

        assert vasp_run.atomic_symbols == ["Li"] + 4 * ["Fe"] + 4 * ["P"] + 16 * ["O"]
        assert vasp_run.final_structure.reduced_formula == "LiFe4(PO4)4"
        assert isinstance(vasp_run.incar, Incar), f"{vasp_run.incar=}"
        assert isinstance(vasp_run.kpoints, Kpoints), f"{vasp_run.kpoints=}"
        assert isinstance(vasp_run.eigenvalues, dict), f"{vasp_run.eigenvalues=}"
        assert vasp_run.final_energy == approx(-269.38319884, abs=1e-7)
        assert vasp_run.tdos.get_gap() == approx(2.0698, abs=1e-4)
        expected = (2.539, 4.0906, 1.5516, False)
        assert vasp_run.eigenvalue_band_properties == approx(expected)
        assert vasp_run.is_hubbard is False
        assert vasp_run.potcar_symbols == [
            "PAW_PBE Li 17Jan2003",
            "PAW_PBE Fe 06Sep2000",
            "PAW_PBE Fe 06Sep2000",
            "PAW_PBE P 17Jan2003",
            "PAW_PBE O 08Apr2002",
        ]
        assert isinstance(vasp_run.kpoints, Kpoints), f"{vasp_run.kpoints=}"
        assert isinstance(vasp_run.actual_kpoints, list), f"{vasp_run.actual_kpoints=}"
        assert isinstance(vasp_run.actual_kpoints_weights, list), f"{vasp_run.actual_kpoints_weights=}"
        assert isinstance(vasp_run.actual_kpoints_weights[0], float), f"{vasp_run.actual_kpoints_weights[0]=}"
        for atom_doses in vasp_run.pdos:
            for orbital_dos in atom_doses:
                assert isinstance(orbital_dos, Orbital), f"{orbital_dos=}"

        # test skipping ionic steps.
        vasprun_skip = Vasprun(filepath, 3, parse_potcar_file=False)
        assert vasprun_skip.nionic_steps == 29
        assert len(vasprun_skip.ionic_steps) == int(vasp_run.nionic_steps / 3) + 1
        assert len(vasprun_skip.ionic_steps) == len(vasprun_skip.structures)
        assert len(vasprun_skip.ionic_steps) == int(vasp_run.nionic_steps / 3) + 1
        # Check that nionic_steps is preserved no matter what.
        assert vasprun_skip.nionic_steps == vasp_run.nionic_steps

        assert vasprun_skip.final_energy != approx(vasp_run.final_energy)

        # Test with ionic_step_offset
        vasprun_offset = Vasprun(filepath, 3, 6, parse_potcar_file=False)
        assert len(vasprun_offset.ionic_steps) == len(vasp_run.ionic_steps) // 3 - 1
        assert vasprun_offset.structures[0] == vasprun_skip.structures[2]

        assert vasprun_ggau.is_hubbard
        assert vasprun_ggau.hubbards["Fe"] == approx(4.3)
        assert vasprun_ggau.projected_eigenvalues[Spin.up][0][0][96][0] == approx(0.0032)
        dct = vasprun_ggau.as_dict()
        assert dct["elements"] == ["Fe", "Li", "O", "P"]
        assert dct["nelements"] == 4

        entry = vasp_run.get_computed_entry(inc_structure=True)
        entry_id_toks = entry.entry_id.split("-")
        assert entry_id_toks[0] == "vasprun"
        assert entry_id_toks[1] == "20100729"
        assert entry_id_toks[2] == "15.0"
        assert entry_id_toks[3] == "da7b01a471dc249323505c0676ae7350"

        assert entry.parameters["run_type"] == "PBEO or other Hybrid Functional"

    def test_unconverged(self):
        filepath = f"{VASP_OUT_DIR}/vasprun.unconverged.xml.gz"
        with pytest.warns(
            UnconvergedVASPWarning,
            match="vasprun.unconverged.xml.gz is an unconverged VASP run",
        ) as warns:
            vasprun_unconverged = Vasprun(filepath, parse_potcar_file=False)
        assert len(warns) >= 1

        assert vasprun_unconverged.converged_ionic
        assert not vasprun_unconverged.converged_electronic
        assert not vasprun_unconverged.converged

    @pytest.mark.filterwarnings("ignore:MaterialsProjectCompatibility is deprecated")
    def test_dfpt(self):
        filepath = f"{VASP_OUT_DIR}/vasprun.dfpt.xml.gz"
        vasprun_dfpt = Vasprun(filepath, parse_potcar_file=False)
        assert vasprun_dfpt.epsilon_static[0][0] == approx(3.26105533)
        assert vasprun_dfpt.epsilon_static[0][1] == approx(-0.00459066)
        assert vasprun_dfpt.epsilon_static[2][2] == approx(3.24330517)
        assert vasprun_dfpt.epsilon_static_wolfe[0][0] == approx(3.33402531)
        assert vasprun_dfpt.epsilon_static_wolfe[0][1] == approx(-0.00559998)
        assert vasprun_dfpt.epsilon_static_wolfe[2][2] == approx(3.31237357)
        assert vasprun_dfpt.converged

        entry = vasprun_dfpt.get_computed_entry()
        entry = MaterialsProjectCompatibility(check_potcar_hash=False).process_entry(entry)
        assert entry.uncorrected_energy + entry.correction == approx(entry.energy)

    @pytest.mark.filterwarnings("ignore::pymatgen.io.vasp.outputs.UnconvergedVASPWarning")
    def test_dfpt_ionic(self):
        filepath = f"{VASP_OUT_DIR}/vasprun.dfpt.ionic.xml.gz"
        vasprun_dfpt_ionic = Vasprun(filepath, parse_potcar_file=False)
        assert vasprun_dfpt_ionic.epsilon_ionic[0][0] == approx(515.73485838)
        assert vasprun_dfpt_ionic.epsilon_ionic[0][1] == approx(-0.00263523)
        assert vasprun_dfpt_ionic.epsilon_ionic[2][2] == approx(19.02110169)

    def test_dfpt_unconverged(self):
        filepath = f"{VASP_OUT_DIR}/vasprun.dfpt.unconverged.xml.gz"
        with pytest.warns(UnconvergedVASPWarning):
            vasprun_dfpt_unconverged = Vasprun(filepath, parse_potcar_file=False)
        assert not vasprun_dfpt_unconverged.converged_electronic
        assert vasprun_dfpt_unconverged.converged_ionic
        assert not vasprun_dfpt_unconverged.converged

    def test_chi(self):
        filepath = f"{VASP_OUT_DIR}/vasprun.chi.xml.gz"
        vasprun_chi = Vasprun(filepath, parse_potcar_file=False)
        assert vasprun_chi.incar.get("ALGO") == "Chi"

    def test_uniform(self):
        vasprun_uniform = Vasprun(f"{VASP_OUT_DIR}/vasprun.uniform.xml.gz", parse_potcar_file=False)
        assert vasprun_uniform.kpoints.style == Kpoints.supported_modes.Reciprocal

    def test_no_projected(self):
        vasprun_no_pdos = Vasprun(f"{VASP_OUT_DIR}/vasprun_Li_no_projected.xml.gz", parse_potcar_file=False)
        assert vasprun_no_pdos.complete_dos is not None
        assert not vasprun_no_pdos.dos_has_errors

    def test_dielectric(self):
        vasprun_diel = Vasprun(f"{VASP_OUT_DIR}/vasprun.dielectric.xml.gz", parse_potcar_file=False)
        assert vasprun_diel.dielectric[0][10] == approx(0.4294)
        assert vasprun_diel.dielectric[1][51][0] == approx(19.941)
        assert vasprun_diel.dielectric[1][51][1] == approx(19.941)
        assert vasprun_diel.dielectric[1][51][2] == approx(19.941)
        assert vasprun_diel.dielectric[1][51][3] == approx(0.0)
        assert vasprun_diel.dielectric[2][85][0] == approx(34.186)
        assert vasprun_diel.dielectric[2][85][1] == approx(34.186)
        assert vasprun_diel.dielectric[2][85][2] == approx(34.186)
        assert vasprun_diel.dielectric[2][85][3] == approx(0.0)

    def test_dielectric_vasp608(self):
        # test reading dielectric constant in vasp 6.0.8
        vasp_xml_path = f"{VASP_OUT_DIR}/vasprun.dielectric_6.0.8.xml.gz"
        vasprun_diel = Vasprun(vasp_xml_path, parse_potcar_file=False)
        assert vasprun_diel.dielectric[0][10] == approx(0.4338)
        assert vasprun_diel.dielectric[1][51][0] == approx(5.267)
        assert vasprun_diel.dielectric_data["density"][0][10] == approx(0.4338)
        assert vasprun_diel.dielectric_data["density"][1][51][0] == approx(5.267)
        assert vasprun_diel.dielectric_data["velocity"][0][10] == approx(0.4338)
        assert vasprun_diel.dielectric_data["velocity"][1][51][0] == approx(1.0741)
        assert len(vasprun_diel.dielectric_data) == 2

    def test_indirect_vasprun(self):
        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.indirect.xml.gz")
        gap, cbm, vbm, direct = vasp_run.eigenvalue_band_properties
        assert gap == approx(0.6119)
        assert cbm == approx(6.2231)
        assert vbm == approx(5.6112)
        assert not direct

    def test_optical_vasprun(self):
        vasp_xml_path = f"{VASP_OUT_DIR}/vasprun.optical_transitions.xml.gz"
        vasprun_optical = Vasprun(vasp_xml_path, parse_potcar_file=False)
        optical_trans = vasprun_optical.optical_transition
        assert optical_trans[0][0] == approx(3.084)
        assert optical_trans[3][0] == approx(3.087)
        assert optical_trans[0][1] == approx(0.001)
        assert optical_trans[1][1] == approx(0.001)
        assert optical_trans[7][1] == approx(0.001)
        assert optical_trans[19][1] == approx(0.001)
        assert optical_trans[54][0] == approx(3.3799999999)
        assert optical_trans[55][0] == approx(3.381)
        assert optical_trans[56][0] == approx(3.381)
        assert optical_trans[54][1] == approx(10554.9860)
        assert optical_trans[55][1] == approx(0.0)
        assert optical_trans[56][1] == approx(0.001)

    def test_force_constants(self):
        vasprun_fc = Vasprun(f"{VASP_OUT_DIR}/vasprun.dfpt.phonon.xml.gz", parse_potcar_file=False)
        assert vasprun_fc.force_constants.shape == (16, 16, 3, 3)
        assert_allclose(
            vasprun_fc.force_constants[8, 9],
            [
                [-0.00184451, 0, 0],
                [0, -0.00933824, -0.03021279],
                [0, -0.03021279, 0.01202547],
            ],
        )
        assert vasprun_fc.normalmode_eigenvals.size == 48
        assert_allclose(
            vasprun_fc.normalmode_eigenvals[17:29],
            [
                -0.59067079,
                -0.59067079,
                -0.59067003,
                -0.59067003,
                -0.59067003,
                -0.59067003,
                -0.585009,
                -0.585009,
                -0.58500895,
                -0.58500883,
                -0.5062956,
                -0.5062956,
            ],
        )
        assert vasprun_fc.normalmode_eigenvecs.shape == (48, 16, 3)
        assert_allclose(
            vasprun_fc.normalmode_eigenvecs[33],
            [
                [0.0884346, -0.08837289, -0.24995639],
                [-0.0884346, 0.08837289, 0.24995639],
                [0.15306645, -0.05105771, -0.14441306],
                [-0.15306645, 0.05105771, 0.14441306],
                [-0.0884346, 0.08837289, 0.24995639],
                [0.0884346, -0.08837289, -0.24995639],
                [-0.15306645, 0.05105771, 0.14441306],
                [0.15306645, -0.05105771, -0.14441306],
                [-0.0884346, 0.08837289, 0.24995639],
                [0.0884346, -0.08837289, -0.24995639],
                [-0.15306645, 0.05105771, 0.14441306],
                [0.15306645, -0.05105771, -0.14441306],
                [0.0884346, -0.08837289, -0.24995639],
                [-0.0884346, 0.08837289, 0.24995639],
                [0.15306645, -0.05105771, -0.14441306],
                [-0.15306645, 0.05105771, 0.14441306],
            ],
        )

    def test_xe(self):
        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.xe.xml.gz", parse_potcar_file=False)
        assert vasp_run.atomic_symbols == ["Xe"]

    def test_invalid_element(self):
        with pytest.raises(ValueError, match="'Z' is not a valid Element"):
            Vasprun(f"{VASP_OUT_DIR}/vasprun.wrong_sp.xml.gz")

    def test_selective_dynamics(self):
        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun.indirect.xml.gz")
        assert list(vasp_run.final_structure.site_properties.get("selective_dynamics")) == [[True] * 3, [False] * 3]

    def test_as_dict(self):
        filepath = f"{VASP_OUT_DIR}/vasprun.xml.gz"
        vasp_run = Vasprun(filepath, parse_potcar_file=False)
        # Test that as_dict() is json-serializable
        assert json.dumps(vasp_run.as_dict()) is not None
        assert vasp_run.as_dict()["input"]["potcar_type"] == [
            "PAW_PBE",
            "PAW_PBE",
            "PAW_PBE",
            "PAW_PBE",
            "PAW_PBE",
        ]
        assert vasp_run.as_dict()["input"]["nkpoints"] == 24

    def test_get_band_structure(self):
        filepath = f"{VASP_OUT_DIR}/vasprun_Si_bands.xml.gz"
        vasp_run = Vasprun(filepath, parse_projected_eigen=True, parse_potcar_file=False)
        band_struct = vasp_run.get_band_structure(kpoints_filename=f"{VASP_IN_DIR}/KPOINTS_Si_bands")
        cbm = band_struct.get_cbm()
        vbm = band_struct.get_vbm()
        assert cbm["kpoint_index"] == [13], "wrong cbm kpoint index"
        assert cbm["energy"] == approx(6.2301), "wrong cbm energy"
        assert cbm["band_index"] == {Spin.up: [4], Spin.down: [4]}, "wrong cbm bands"
        assert vbm["kpoint_index"] == [0, 63, 64]
        assert vbm["energy"] == approx(5.6158), "wrong vbm energy"
        assert vbm["band_index"] == {
            Spin.up: [1, 2, 3],
            Spin.down: [1, 2, 3],
        }, "wrong vbm bands"
        assert vbm["kpoint"].label == "\\Gamma", "wrong vbm label"
        assert cbm["kpoint"].label is None, "wrong cbm label"

        projected = band_struct.get_projection_on_elements()
        assert projected[Spin.up][0][0]["Si"] == approx(0.4238)
        projected = band_struct.get_projections_on_elements_and_orbitals({"Si": ["s"]})
        assert projected[Spin.up][0][0]["Si"]["s"] == approx(0.4238)

        # Test compressed files case 1: compressed KPOINTS in current dir
        copyfile(f"{VASP_OUT_DIR}/vasprun_Si_bands.xml.gz", "vasprun.xml.gz")

        # Check for error if no KPOINTS file
        vasp_run = Vasprun("vasprun.xml.gz", parse_projected_eigen=True, parse_potcar_file=False)
        with pytest.raises(
            VaspParseError,
            match="KPOINTS not found but needed to obtain band structure along symmetry lines",
        ):
            _ = vasp_run.get_band_structure(line_mode=True)

        # Check KPOINTS.gz successfully inferred and used if present
        with (
            open(f"{VASP_IN_DIR}/KPOINTS_Si_bands", "rb") as f_in,
            gzip.open("KPOINTS.gz", "wb") as f_out,
        ):
            copyfileobj(f_in, f_out)
        bs_kpts_gzip = vasp_run.get_band_structure()
        assert band_struct.efermi == bs_kpts_gzip.efermi
        assert band_struct.as_dict() == bs_kpts_gzip.as_dict()

        # Test compressed files case 2: compressed vasprun in another dir
        os.mkdir("deeper")
        copyfile(f"{VASP_IN_DIR}/KPOINTS_Si_bands", Path("deeper") / "KPOINTS")
        copyfile(f"{VASP_OUT_DIR}/vasprun_Si_bands.xml.gz", Path("deeper") / "vasprun.xml.gz")
        vasp_run = Vasprun(
            f"{'deeper'}/vasprun.xml.gz",
            parse_projected_eigen=True,
            parse_potcar_file=False,
        )
        bs_vasprun_gzip = vasp_run.get_band_structure(line_mode=True)
        assert band_struct.efermi == bs_vasprun_gzip.efermi
        assert band_struct.as_dict() == bs_vasprun_gzip.as_dict()

        # test hybrid band structures
        vasp_run.actual_kpoints_weights[-1] = 0.0
        band_struct = vasp_run.get_band_structure(kpoints_filename=f"{VASP_IN_DIR}/KPOINTS_Si_bands")
        cbm = band_struct.get_cbm()
        vbm = band_struct.get_vbm()
        assert cbm["kpoint_index"] == [0]
        assert cbm["energy"] == approx(6.3676)
        assert cbm["kpoint"].label is None
        assert vbm["kpoint_index"] == [0]
        assert vbm["energy"] == approx(2.8218)
        assert vbm["kpoint"].label is None

        # test self-consistent band structure calculation for non-hybrid functionals
        vasp_run = Vasprun(
            f"{VASP_OUT_DIR}/vasprun.force_hybrid_like_calc.xml.gz",
            parse_projected_eigen=True,
            parse_potcar_file=False,
        )
        band_struct = vasp_run.get_band_structure(
            kpoints_filename=f"{VASP_IN_DIR}/KPOINTS_force_hybrid_like_calc",
            force_hybrid_mode=True,
            line_mode=True,
        )

        dict_to_test = band_struct.get_band_gap()

        assert dict_to_test["direct"]
        assert dict_to_test["energy"] == approx(6.007899999999999)
        assert dict_to_test["transition"] == "\\Gamma-\\Gamma"
        assert band_struct.get_branch(0)[0]["start_index"] == 0
        assert band_struct.get_branch(0)[0]["end_index"] == 0

    def test_projected_magnetization(self):
        filepath = f"{VASP_OUT_DIR}/vasprun.lvel.Si2H.xml.gz"
        vasp_run = Vasprun(filepath, parse_projected_eigen=True)
        assert vasp_run.projected_magnetization is not None
        assert vasp_run.projected_magnetization.shape == (76, 240, 4, 9, 3)
        assert vasp_run.projected_magnetization[0, 0, 0, 0, 0] == approx(-0.0712)

    def test_smart_efermi(self):
        # branch 1 - E_fermi does not cross a band
        vrun = Vasprun(f"{VASP_OUT_DIR}/vasprun.LiF.xml.gz")
        smart_fermi = vrun.calculate_efermi()
        assert smart_fermi == approx(vrun.efermi, abs=1e-4)
        eigen_gap = vrun.eigenvalue_band_properties[0]
        bs_gap = vrun.get_band_structure(efermi=smart_fermi).get_band_gap()["energy"]
        assert bs_gap == approx(eigen_gap, abs=1e-3)

        # branch 2 - E_fermi crosses a band but bandgap=0
        vrun = Vasprun(f"{VASP_OUT_DIR}/vasprun.Al.xml.gz")
        smart_fermi = vrun.calculate_efermi()
        assert smart_fermi == approx(vrun.efermi, abs=1e-4)
        eigen_gap = vrun.eigenvalue_band_properties[0]
        bs_gap = vrun.get_band_structure(efermi=smart_fermi).get_band_gap()["energy"]
        assert bs_gap == approx(eigen_gap, abs=1e-3)

        # branch 3 - E_fermi crosses a band in an insulator
        vrun = Vasprun(f"{VASP_OUT_DIR}/vasprun.LiH_bad_efermi.xml.gz")
        smart_fermi = vrun.calculate_efermi()
        assert smart_fermi != approx(vrun.efermi, abs=1e-4)
        eigen_gap = vrun.eigenvalue_band_properties[0]
        bs_gap = vrun.get_band_structure(efermi="smart").get_band_gap()["energy"]
        assert bs_gap == approx(eigen_gap, abs=1e-3)
        assert vrun.get_band_structure(efermi=None).get_band_gap()["energy"] != approx(eigen_gap, abs=1e-3)
        assert bs_gap != 0

        # branch 4 - E_fermi incorrectly placed inside a band
        vrun = Vasprun(f"{VASP_OUT_DIR}/vasprun.bad_fermi.xml.gz")
        smart_fermi = vrun.calculate_efermi()
        assert smart_fermi == approx(6.0165)

    def test_float_overflow(self):
        # test we interpret VASP's *********** for overflowed values as NaNs
        # https://github.com/materialsproject/pymatgen/pull/3452
        filepath = f"{VASP_OUT_DIR}/vasprun.sc_overflow.xml.gz"
        with pytest.warns(UserWarning, match="Float overflow .* encountered in vasprun"):
            vasp_run = Vasprun(filepath)
        first_ionic_step = vasp_run.ionic_steps[0]
        elec_step = first_ionic_step["electronic_steps"][29]
        assert np.isnan(elec_step["e_wo_entrp"])
        assert np.isnan(elec_step["e_fr_energy"])
        assert np.isnan(first_ionic_step["forces"]).any()

    def test_update_potcar(self):
        filepath = f"{VASP_OUT_DIR}/vasprun.xml.gz"
        potcar_path = f"{VASP_IN_DIR}/POTCAR_LiFePO4.gz"
        potcar_path2 = f"{VASP_IN_DIR}/POTCAR_2_LiFePO4.gz"

        vasp_run = Vasprun(filepath, parse_potcar_file=False)
        potcars = Potcar.from_file(potcar_path)
        expected_spec = [{"titel": titel, "hash": None, "summary_stats": {}} for titel in vasp_run.potcar_symbols]
        assert vasp_run.potcar_spec == expected_spec

        vasp_run.update_potcar_spec(potcar_path)
        potcars = Potcar.from_file(potcar_path)
        expected_spec = []
        for titel in vasp_run.potcar_symbols:
            for potcar in potcars:
                if titel == potcar.TITEL:
                    break
            expected_spec += [
                {
                    "titel": titel,
                    "hash": potcar.md5_header_hash,
                    "summary_stats": potcar._summary_stats,
                }
            ]
        assert vasp_run.potcar_spec == expected_spec

        with pytest.warns(UserWarning, match="No POTCAR file with matching TITEL fields was found in"):
            Vasprun(filepath, parse_potcar_file=potcar_path2)

        vasp_run = Vasprun(filepath, parse_potcar_file=potcar_path)
        assert vasp_run.potcar_spec == expected_spec

        with pytest.warns(UserWarning, match="No POTCAR file with matching TITEL fields was found in"):
            Vasprun(filepath, parse_potcar_file=potcar_path2)

    def test_search_for_potcar(self):
        filepath = f"{VASP_OUT_DIR}/vasprun.xml.gz"
        vasp_run = Vasprun(filepath, parse_potcar_file=True)
        assert [spec["titel"] for spec in vasp_run.potcar_spec] == [
            "PAW_PBE Li 17Jan2003",
            "PAW_PBE Fe 06Sep2000",
            "PAW_PBE Fe 06Sep2000",
            "PAW_PBE P 17Jan2003",
            "PAW_PBE O 08Apr2002",
        ]

    def test_potcar_not_found(self):
        filepath = f"{VASP_OUT_DIR}/vasprun.xml.gz"
        # Ensure no potcar is found and nothing is updated
        with pytest.warns(UserWarning, match="No POTCAR file with matching TITEL fields was found in") as warns:
            vasp_run = Vasprun(filepath, parse_potcar_file=".")
        assert len(warns) >= 2, f"{len(warns)=}"
        assert vasp_run.potcar_spec == [
            {"titel": "PAW_PBE Li 17Jan2003", "hash": None, "summary_stats": {}},
            {"titel": "PAW_PBE Fe 06Sep2000", "hash": None, "summary_stats": {}},
            {"titel": "PAW_PBE Fe 06Sep2000", "hash": None, "summary_stats": {}},
            {"titel": "PAW_PBE P 17Jan2003", "hash": None, "summary_stats": {}},
            {"titel": "PAW_PBE O 08Apr2002", "hash": None, "summary_stats": {}},
        ]

    def test_parsing_chemical_shift_calculations(self):
        filepath = f"{TEST_DIR}/fixtures/nmr/cs/basic/vasprun.xml.chemical_shift.scstep"
        vasp_run = Vasprun(filepath)

        n_estep = len(vasp_run.ionic_steps[-1]["electronic_steps"])
        assert n_estep == 10
        assert vasp_run.converged

    def test_parsing_efg_calcs(self):
        filepath = f"{TEST_DIR}/fixtures/nmr/efg/AlPO4/vasprun.xml"
        vasp_run = Vasprun(filepath)
        n_elec_steps = len(vasp_run.ionic_steps[-1]["electronic_steps"])
        assert n_elec_steps == 18
        assert vasp_run.converged

    def test_charged_structure(self):
        vpath = f"{VASP_OUT_DIR}/vasprun.charged.xml.gz"
        potcar_path = f"{FAKE_POTCAR_DIR}/POT_GGA_PAW_PBE/POTCAR.Si.gz"
        vasp_run = Vasprun(vpath, parse_potcar_file=False)
        vasp_run.update_charge_from_potcar(potcar_path)
        assert vasp_run.parameters.get("NELECT", 8) == 9
        assert vasp_run.structures[0].charge == -1
        assert vasp_run.initial_structure.charge == -1
        assert vasp_run.final_structure.charge == -1

        vpath = f"{VASP_OUT_DIR}/vasprun.split.charged.xml.gz"
        potcar_path = f"{VASP_IN_DIR}/POTCAR_split_charged.gz"
        vasp_run = Vasprun(vpath, parse_potcar_file=False)
        vasp_run.update_charge_from_potcar(potcar_path)
        assert vasp_run.parameters.get("NELECT", 0) == 7
        assert vasp_run.structures[-1].charge == -1
        assert vasp_run.initial_structure.charge == -1
        assert vasp_run.final_structure.charge == -1

    def test_kpointset_electronvelocities(self):
        vpath = f"{VASP_OUT_DIR}/vasprun.lvel.Si2H.xml.gz"
        vasp_run = Vasprun(vpath, parse_potcar_file=False)
        assert vasp_run.eigenvalues[Spin.up].shape[0] == len(vasp_run.actual_kpoints)

    def test_eigenvalue_band_properties_separate_spins(self):
        eig = Vasprun(f"{VASP_OUT_DIR}/vasprun_eig_separate_spins.xml.gz", separate_spins=True)
        props = eig.eigenvalue_band_properties
        eig2 = Vasprun(f"{VASP_OUT_DIR}/vasprun_eig_separate_spins.xml.gz", separate_spins=False)
        props2 = eig2.eigenvalue_band_properties
        assert props[0][0] == approx(2.8772, abs=1e-4)
        assert props[0][1] == approx(1.2810, abs=1e-4)
        assert props[1][0] == approx(3.6741, abs=1e-4)
        assert props[1][1] == approx(1.6225, abs=1e-4)
        assert props[2][0] == approx(0.7969, abs=1e-4)
        assert props[2][1] == approx(0.3415, abs=1e-4)
        assert props2[0] == approx(np.min(props[1]) - np.max(props[2]), abs=1e-4)
        assert props[3][0]
        assert props[3][1]

    def test_kpoints_opt(self):
        vasp_run = Vasprun(f"{TEST_DIR}/fixtures/kpoints_opt/vasprun.xml.gz", parse_projected_eigen=True)
        # This calculation was run using KPOINTS_OPT
        kpt_opt_props = vasp_run.kpoints_opt_props
        # Check the k-points were read correctly.
        assert len(vasp_run.actual_kpoints) == 10
        assert len(vasp_run.actual_kpoints_weights) == 10
        assert len(kpt_opt_props.actual_kpoints) == 100
        assert len(kpt_opt_props.actual_kpoints_weights) == 100
        # Check the eigenvalues were read correctly.
        assert vasp_run.eigenvalues[Spin.up].shape == (10, 24, 2)
        assert kpt_opt_props.eigenvalues[Spin.up].shape == (100, 24, 2)
        assert vasp_run.eigenvalues[Spin.up][0, 0, 0] == approx(-6.1471)
        assert kpt_opt_props.eigenvalues[Spin.up][0, 0, 0] == approx(-6.1536)
        # Check the projected eigenvalues were read correctly
        assert vasp_run.projected_eigenvalues[Spin.up].shape == (10, 24, 8, 9)
        assert kpt_opt_props.projected_eigenvalues[Spin.up].shape == (100, 24, 8, 9)
        assert vasp_run.projected_eigenvalues[Spin.up][0, 1, 0, 0] == approx(0.0492)
        # I think these zeroes are a bug in VASP (maybe my VASP) transcribing from PROCAR_OPT to vasprun.xml
        # No matter. The point of the parser is to read what's in the file.
        assert kpt_opt_props.projected_eigenvalues[Spin.up][0, 1, 0, 0] == approx(0.0000)
        # Test as_dict
        vasp_run_dct = vasp_run.as_dict()
        assert vasp_run_dct["input"]["nkpoints_opt"] == 100
        assert vasp_run_dct["input"]["nkpoints"] == 10
        assert vasp_run_dct["output"]["eigenvalues_kpoints_opt"]["1"][0][0][0] == approx(-6.1536)

    def test_kpoints_opt_band_structure(self):
        vasp_run = Vasprun(
            f"{TEST_DIR}/fixtures/kpoints_opt/vasprun.xml.gz", parse_potcar_file=False, parse_projected_eigen=True
        )
        bs = vasp_run.get_band_structure(f"{TEST_DIR}/fixtures/kpoints_opt/KPOINTS_OPT")
        assert isinstance(bs, BandStructureSymmLine)
        cbm = bs.get_cbm()
        vbm = bs.get_vbm()
        assert cbm["kpoint_index"] == [38], "wrong cbm kpoint index"
        assert cbm["energy"] == approx(6.4394), "wrong cbm energy"
        assert cbm["band_index"] == {Spin.up: [16], Spin.down: [16]}, "wrong cbm bands"
        assert vbm["kpoint_index"] == [0, 39, 40]
        assert vbm["energy"] == approx(5.7562), "wrong vbm energy"
        assert vbm["band_index"] == {
            Spin.down: [13, 14, 15],
            Spin.up: [13, 14, 15],
        }, "wrong vbm bands"
        vbm_kp_label = vbm["kpoint"].label
        assert vbm["kpoint"].label == "\\Gamma", f"Unpexpected {vbm_kp_label=}"
        cmb_kp_label = cbm["kpoint"].label
        assert cmb_kp_label is None, f"Unpexpected {cmb_kp_label=}"
        # Test projection
        projected = bs.get_projection_on_elements()
        assert np.isnan(projected[Spin.up][0][0]["Si"])
        # Due to some error in my VASP, the transcription of PROCAR_OPT into
        # vasprun.xml is filled to the brim with errors in the projections.
        # At some point we might get a healthier vasprun.xml, but the point here
        # is to test the parser, not VASP.
        projected = bs.get_projections_on_elements_and_orbitals({"Si": ["s"]})
        assert projected[Spin.up][0][58]["Si"]["s"] == approx(-0.0271)

    def test_parse_potcar_cwd_relative(self):
        # Test to ensure that common use cases of vasprun parsing work
        # in the current working directory using relative paths,
        # either when leading ./ is specified or not for vasprun.xml
        # See gh-3586
        copyfile(f"{VASP_OUT_DIR}/vasprun.Al.xml.gz", "vasprun.xml.gz")

        potcar_path = f"{VASP_IN_DIR}/fake_potcars/POTPAW_PBE_54/POTCAR.Al.gz"
        copyfile(potcar_path, "POTCAR.gz")

        potcar = Potcar.from_file(potcar_path)
        for leading_path in ("", "./"):
            vrun = Vasprun(f"{leading_path}vasprun.xml.gz", parse_potcar_file=True)
            # Note that the TITEL is not updated in Vasprun.potcar_spec
            # Since the fake POTCARs modify the TITEL (to indicate fakeness), can't compare
            for ipot in range(len(potcar)):
                assert vrun.potcar_spec[ipot]["hash"] == potcar[ipot].md5_header_hash
                assert vrun.potcar_spec[ipot]["summary_stats"] == potcar[ipot]._summary_stats


class TestOutcar:
    @pytest.mark.parametrize("compressed", [True, False])
    def test_init(self, compressed: bool):
        """Test from both compressed and uncompressed versions,
        as there was a bug in monty causing different behaviours.
        """
        if compressed:
            outcar = Outcar(f"{VASP_OUT_DIR}/OUTCAR.gz")
        else:
            copyfile(f"{VASP_OUT_DIR}/OUTCAR.gz", "./OUTCAR.gz")
            decompress_file("./OUTCAR.gz")
            outcar = Outcar("./OUTCAR")

        expected_mag = (
            {"d": 0.0, "p": 0.003, "s": 0.002, "tot": 0.005},
            {"d": 0.798, "p": 0.008, "s": 0.007, "tot": 0.813},
            {"d": 0.798, "p": 0.008, "s": 0.007, "tot": 0.813},
            {"d": 0.0, "p": -0.117, "s": 0.005, "tot": -0.112},
            {"d": 0.0, "p": -0.165, "s": 0.004, "tot": -0.162},
            {"d": 0.0, "p": -0.117, "s": 0.005, "tot": -0.112},
            {"d": 0.0, "p": -0.165, "s": 0.004, "tot": -0.162},
        )
        expected_chg = (
            {"p": 0.154, "s": 0.078, "d": 0.0, "tot": 0.232},
            {"p": 0.707, "s": 0.463, "d": 8.316, "tot": 9.486},
            {"p": 0.707, "s": 0.463, "d": 8.316, "tot": 9.486},
            {"p": 3.388, "s": 1.576, "d": 0.0, "tot": 4.964},
            {"p": 3.365, "s": 1.582, "d": 0.0, "tot": 4.947},
            {"p": 3.388, "s": 1.576, "d": 0.0, "tot": 4.964},
            {"p": 3.365, "s": 1.582, "d": 0.0, "tot": 4.947},
        )

        assert outcar.magnetization == approx(expected_mag, abs=1e-5), "Wrong magnetization read from Outcar"
        assert outcar.charge == approx(expected_chg, abs=1e-5), "Wrong charge read from Outcar"
        assert not outcar.is_stopped
        assert outcar.run_stats == {
            "System time (sec)": 0.938,
            "Total CPU time used (sec)": 545.142,
            "Elapsed time (sec)": 546.709,
            "Maximum memory used (kb)": 0.0,
            "Average memory used (kb)": 0.0,
            "User time (sec)": 544.204,
            "cores": 8,
        }
        assert outcar.efermi == approx(2.0112)
        assert outcar.nelect == approx(44.9999991)
        assert outcar.total_mag == approx(0.9999998)

        assert outcar.as_dict() is not None

        assert not outcar.lepsilon

        toten = 0
        for k in outcar.final_energy_contribs:
            toten += outcar.final_energy_contribs[k]
        assert toten == approx(outcar.final_energy, abs=1e-6)

    def test_stopped_old(self):
        filepath = f"{VASP_OUT_DIR}/OUTCAR.stopped.gz"
        outcar = Outcar(filepath)
        assert outcar.is_stopped

        outcar = Outcar(f"{VASP_OUT_DIR}/OUTCAR.lepsilon_old_born.gz")

        assert outcar.lepsilon
        assert outcar.dielectric_tensor[0][0] == approx(3.716432)
        assert outcar.dielectric_tensor[0][1] == approx(-0.20464)
        assert outcar.dielectric_tensor[1][2] == approx(-0.20464)
        assert outcar.dielectric_ionic_tensor[0][0] == approx(0.001419)
        assert outcar.dielectric_ionic_tensor[0][2] == approx(0.001419)
        assert outcar.dielectric_ionic_tensor[2][2] == approx(0.001419)
        assert outcar.piezo_tensor[0][0] == approx(0.52799)
        assert outcar.piezo_tensor[1][3] == approx(0.35998)
        assert outcar.piezo_tensor[2][5] == approx(0.35997)
        assert outcar.piezo_ionic_tensor[0][0] == approx(0.05868)
        assert outcar.piezo_ionic_tensor[1][3] == approx(0.06241)
        assert outcar.piezo_ionic_tensor[2][5] == approx(0.06242)
        assert outcar.born[0][1][2] == approx(-0.385)
        assert outcar.born[1][2][0] == approx(0.36465)
        assert outcar.internal_strain_tensor[0][0][0] == approx(-572.5437, abs=1e-4)
        assert outcar.internal_strain_tensor[0][1][0] == approx(683.2985, abs=1e-4)
        assert outcar.internal_strain_tensor[0][1][3] == approx(73.07059, abs=1e-4)
        assert outcar.internal_strain_tensor[1][0][0] == approx(570.98927, abs=1e-4)
        assert outcar.internal_strain_tensor[1][1][0] == approx(-683.68519, abs=1e-4)
        assert outcar.internal_strain_tensor[1][2][2] == approx(570.98927, abs=1e-4)

    def test_stopped(self):
        filepath = f"{VASP_OUT_DIR}/OUTCAR.stopped.gz"
        outcar = Outcar(filepath)
        assert outcar.is_stopped

        outcar = Outcar(f"{VASP_OUT_DIR}/OUTCAR.lepsilon.gz")

        assert outcar.lepsilon
        assert outcar.dielectric_tensor[0][0] == approx(3.716432)
        assert outcar.dielectric_tensor[0][1] == approx(-0.20464)
        assert outcar.dielectric_tensor[1][2] == approx(-0.20464)
        assert outcar.dielectric_ionic_tensor[0][0] == approx(0.001419)
        assert outcar.dielectric_ionic_tensor[0][2] == approx(0.001419)
        assert outcar.dielectric_ionic_tensor[2][2] == approx(0.001419)
        assert outcar.piezo_tensor[0][0] == approx(0.52799)
        assert outcar.piezo_tensor[1][3] == approx(0.35998)
        assert outcar.piezo_tensor[2][5] == approx(0.35997)
        assert outcar.piezo_ionic_tensor[0][0] == approx(0.05868)
        assert outcar.piezo_ionic_tensor[1][3] == approx(0.06241)
        assert outcar.piezo_ionic_tensor[2][5] == approx(0.06242)
        assert outcar.born[0][1][2] == approx(-0.385)
        assert outcar.born[1][2][0] == approx(0.36465)
        assert outcar.internal_strain_tensor[0][0][0] == approx(-572.5437, abs=1e-4)
        assert outcar.internal_strain_tensor[0][1][0] == approx(683.2985, abs=1e-4)
        assert outcar.internal_strain_tensor[0][1][3] == approx(73.07059, abs=1e-4)
        assert outcar.internal_strain_tensor[1][0][0] == approx(570.98927, abs=1e-4)
        assert outcar.internal_strain_tensor[1][1][0] == approx(-683.68519, abs=1e-4)
        assert outcar.internal_strain_tensor[1][2][2] == approx(570.98927, abs=1e-4)

    def test_soc(self):
        filepath = f"{VASP_OUT_DIR}/OUTCAR.NiO_SOC.gz"
        outcar = Outcar(filepath)
        expected_mag = (
            {
                "s": Magmom([0.0, 0.0, -0.001]),
                "p": Magmom([0.0, 0.0, -0.003]),
                "d": Magmom([0.0, 0.0, 1.674]),
                "tot": Magmom([0.0, 0.0, 1.671]),
            },
            {
                "s": Magmom([0.0, 0.0, 0.001]),
                "p": Magmom([0.0, 0.0, 0.003]),
                "d": Magmom([0.0, 0.0, -1.674]),
                "tot": Magmom([0.0, 0.0, -1.671]),
            },
            {
                "s": Magmom([0.0, 0.0, 0.0]),
                "p": Magmom([0.0, 0.0, 0.0]),
                "d": Magmom([0.0, 0.0, 0.0]),
                "tot": Magmom([0.0, 0.0, 0.0]),
            },
            {
                "s": Magmom([0.0, 0.0, 0.0]),
                "p": Magmom([0.0, 0.0, 0.0]),
                "d": Magmom([0.0, 0.0, 0.0]),
                "tot": Magmom([0.0, 0.0, 0.0]),
            },
        )
        # test note: Magmom class uses np.allclose() when testing for equality
        # so fine to use == operator here
        assert outcar.magnetization == expected_mag, "Wrong vector magnetization read from Outcar for SOC calculation"

        assert outcar.noncollinear is True

    def test_polarization(self):
        filepath = f"{VASP_OUT_DIR}/OUTCAR.BaTiO3.polar"
        outcar = Outcar(filepath)
        assert outcar.spin
        assert outcar.noncollinear is False
        assert outcar.p_ion == approx([0.0, 0.0, -5.56684])
        assert outcar.p_elec == approx([0.00024, 0.00019, 3.61674])

    def test_pseudo_zval(self):
        filepath = f"{VASP_OUT_DIR}/OUTCAR.BaTiO3.polar"
        outcar = Outcar(filepath)
        assert outcar.zval_dict == approx({"Ba": 10.0, "Ti": 10.0, "O": 6.0})

        filepath = f"{VASP_OUT_DIR}/OUTCAR.LaSnNO2.polar"
        outcar = Outcar(filepath)
        assert outcar.zval_dict == approx({"La": 11.0, "N": 5.0, "O": 6.0, "Sn": 14.0})

    def test_dielectric(self):
        filepath = f"{VASP_OUT_DIR}/OUTCAR.dielectric"
        outcar = Outcar(filepath)
        outcar.read_corrections()
        assert outcar.data["dipol_quadrupol_correction"] == approx(0.03565)
        assert outcar.final_energy == approx(-797.46294064)

    def test_freq_dielectric(self):
        filepath = f"{VASP_OUT_DIR}/OUTCAR.LOPTICS"
        outcar = Outcar(filepath)
        outcar.read_freq_dielectric()
        assert outcar.dielectric_energies[0] == approx(0)
        assert outcar.dielectric_energies[-1] == approx(39.826101)
        assert outcar.dielectric_tensor_function[0][0, 0] == approx(8.96938800)
        assert outcar.dielectric_tensor_function[-1][0, 0] == approx(7.36167000e-01 + 1.53800000e-03j)
        assert len(outcar.dielectric_energies) == len(outcar.dielectric_tensor_function)
        assert_allclose(
            outcar.dielectric_tensor_function[0],
            outcar.dielectric_tensor_function[0].transpose(),
        )

        plasma_freq = outcar.plasma_frequencies
        assert_allclose(plasma_freq["intraband"], np.zeros((3, 3)))
        assert_allclose(
            plasma_freq["interband"],
            [
                [367.49, 63.939, 11.976],
                [63.939, 381.155, -24.461],
                [11.976, -24.461, 297.844],
            ],
        )

    def test_freq_dielectric_vasp544(self):
        filepath = f"{VASP_OUT_DIR}/OUTCAR.LOPTICS.vasp544"
        outcar = Outcar(filepath)
        outcar.read_freq_dielectric()
        assert outcar.dielectric_energies[0] == approx(0)
        assert outcar.dielectric_energies[-1] == approx(39.63964)
        assert outcar.dielectric_tensor_function[0][0, 0] == approx(12.769435 + 0j)
        assert outcar.dielectric_tensor_function[-1][0, 0] == approx(0.828615 + 0.016594j)
        assert len(outcar.dielectric_energies) == len(outcar.dielectric_tensor_function)
        assert_allclose(
            outcar.dielectric_tensor_function[0],
            outcar.dielectric_tensor_function[0].transpose(),
        )

    def test_parse_sci_notation(self):
        invalid_pattern = "23535.35 35235.34 325325.3"
        valid_pattern1 = " 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00"
        valid_pattern2 = " 0.62963E+00 0.15467E+02 0.15467E+02 0.15467E+02-0.30654E-16-0.91612E-16 0.52388E-16"

        assert Outcar._parse_sci_notation(invalid_pattern) == []
        assert Outcar._parse_sci_notation(valid_pattern1) == [0, 0, 0, 0, 0, 0, 0]
        assert Outcar._parse_sci_notation(valid_pattern2) == [
            0.62963,
            0.15467e02,
            0.15467e02,
            0.15467e02,
            -0.30654e-16,
            -0.91612e-16,
            0.52388e-16,
        ]

    def test_read_elastic_tensor(self):
        filepath = f"{VASP_OUT_DIR}/OUTCAR.total_tensor.Li2O.gz"
        outcar = Outcar(filepath)

        outcar.read_elastic_tensor()

        assert outcar.data["elastic_tensor"][0][0] == approx(1986.3391)
        assert outcar.data["elastic_tensor"][0][1] == approx(187.8324)
        assert outcar.data["elastic_tensor"][3][3] == approx(586.3034)

    def test_read_lcalcpol(self):
        # outcar with electrons Angst units
        folder = "io/vasp/fixtures/BTO_221_99_polarization/interpolation_6_polarization/"
        filepath = TEST_FILES_DIR / folder / "OUTCAR"
        outcar = Outcar(filepath)

        outcar.read_lcalcpol()

        p_ion = [0.0, 0.0, 19.70306]
        p_elec = [4.02264, 4.02263, -4.08851]
        p_sp1 = [2.01124, 2.01124, -2.04426]
        p_sp2 = [2.01139, 2.01139, -2.04426]

        assert outcar.p_ion == approx(p_ion)
        assert outcar.p_elec == approx(p_elec)
        assert outcar.p_sp1 == approx(p_sp1)
        assert outcar.p_sp2 == approx(p_sp2)

        # outcar with |e| Angst units
        filepath = f"{VASP_OUT_DIR}/OUTCAR_vasp_6.3.gz"
        outcar = Outcar(filepath)

        outcar.read_lcalcpol()

        p_ion = [-79.03374, -0.0, -28.44354]
        p_elec = [9.01127e00, -1.00000e-05, 3.24308e00]
        p_sp1 = [4.50564, 0.0, 1.62154]
        p_sp2 = [4.50563e00, -1.00000e-05, 1.62154e00]

        assert outcar.p_ion == approx(p_ion)
        assert outcar.p_elec == approx(p_elec)
        assert outcar.p_sp1 == approx(p_sp1)
        assert outcar.p_sp2 == approx(p_sp2)

    def test_read_piezo_tensor(self):
        filepath = f"{VASP_OUT_DIR}/OUTCAR.lepsilon.gz"
        outcar = Outcar(filepath)

        outcar.read_piezo_tensor()
        assert outcar.data["piezo_tensor"][0][0] == approx(0.52799)
        assert outcar.data["piezo_tensor"][1][3] == approx(0.35998)
        assert outcar.data["piezo_tensor"][2][5] == approx(0.35997)

    def test_core_state_eigen(self):
        filepath = f"{VASP_OUT_DIR}/OUTCAR.CL.gz"
        cl = Outcar(filepath).read_core_state_eigen()
        assert cl[6]["2s"][-1] == approx(-174.4779)
        filepath = f"{VASP_OUT_DIR}/OUTCAR.icorelevel"
        outcar = Outcar(filepath)
        cl = outcar.read_core_state_eigen()
        assert cl[4]["3d"][-1] == approx(-31.4522)

        # test serialization
        outcar.as_dict()

    def test_avg_core_poten(self):
        filepath = f"{VASP_OUT_DIR}/OUTCAR.lepsilon.gz"
        cp = Outcar(filepath).read_avg_core_poten()
        assert cp[-1][1] == approx(-90.0487)

        filepath = f"{VASP_OUT_DIR}/OUTCAR.gz"
        cp = Outcar(filepath).read_avg_core_poten()
        assert cp[0][6] == approx(-73.1068)

        filepath = f"{VASP_OUT_DIR}/OUTCAR.bad_core_poten.gz"
        cp = Outcar(filepath).read_avg_core_poten()
        assert cp[0][1] == approx(-101.5055)

    def test_single_atom(self):
        filepath = f"{VASP_OUT_DIR}/OUTCAR.Al"
        outcar = Outcar(filepath)
        expected_mag = ({"p": 0.0, "s": 0.0, "d": 0.0, "tot": 0.0},)
        expected_chg = ({"p": 0.343, "s": 0.425, "d": 0.0, "tot": 0.768},)

        assert outcar.magnetization == approx(expected_mag)
        assert outcar.charge == approx(expected_chg)
        assert not outcar.is_stopped
        assert outcar.run_stats == {
            "System time (sec)": 0.592,
            "Total CPU time used (sec)": 50.194,
            "Elapsed time (sec)": 52.337,
            "Maximum memory used (kb)": 62900.0,
            "Average memory used (kb)": 0.0,
            "User time (sec)": 49.602,
            "cores": 32,
        }
        assert outcar.efermi == approx(8.0942)
        assert outcar.nelect == 3
        assert outcar.total_mag == approx(8.2e-06)

        assert outcar.as_dict() is not None

    def test_chemical_shielding(self):
        filename = f"{TEST_DIR}/fixtures/nmr/cs/core.diff/hydromagnesite/OUTCAR"
        outcar = Outcar(filename)
        expected_chemical_shielding = [
            [191.9974, 69.5232, 0.6342],
            [195.0808, 68.183, 0.833],
            [192.0389, 69.5762, 0.6329],
            [195.0844, 68.1756, 0.8336],
            [192.005, 69.5289, 0.6339],
            [195.0913, 68.1859, 0.833],
            [192.0237, 69.565, 0.6333],
            [195.0788, 68.1733, 0.8337],
        ]

        assert len(outcar.data["chemical_shielding"]["valence_only"][20:28]) == approx(len(expected_chemical_shielding))

        assert_allclose(
            outcar.data["chemical_shielding"]["valence_and_core"][20:28],
            expected_chemical_shielding,
            atol=1e-5,
        )

    def test_chemical_shielding_with_different_core_contribution(self):
        filename = f"{TEST_DIR}/fixtures/nmr/cs/core.diff/core.diff.chemical.shifts.OUTCAR"
        outcar = Outcar(filename)
        c_vo = outcar.data["chemical_shielding"]["valence_only"][7]
        assert list(c_vo) == approx([198.7009, 73.7484, 1])
        c_vc = outcar.data["chemical_shielding"]["valence_and_core"][7]
        assert list(c_vc) == approx([-1.9406, 73.7484, 1])

    def test_cs_raw_tensors(self):
        filename = f"{TEST_DIR}/fixtures/nmr/cs/core.diff/core.diff.chemical.shifts.OUTCAR"
        outcar = Outcar(filename)
        unsym_tensors = outcar.data["unsym_cs_tensor"]
        assert unsym_tensors[0] == [
            [-145.814605, -4.263425, 0.000301],
            [4.263434, -145.812238, -8.7e-05],
            [0.000136, -0.000189, -142.794068],
        ]
        assert unsym_tensors[29] == [
            [287.789318, -53.799325, 30.900024],
            [-53.799571, 225.668117, -17.839598],
            [3.801103, -2.195218, 88.896756],
        ]

    def test_cs_g0_contribution(self):
        filename = f"{TEST_DIR}/fixtures/nmr/cs/core.diff/core.diff.chemical.shifts.OUTCAR"
        outcar = Outcar(filename)
        g0_contrib = outcar.data["cs_g0_contribution"]
        assert g0_contrib == [
            [-8.773535, 9e-06, 1e-06],
            [1.7e-05, -8.773536, -0.0792],
            [-6e-06, -0.008328, -9.320237],
        ]

    def test_cs_core_contribution(self):
        filename = f"{TEST_DIR}/fixtures/nmr/cs/core.diff/core.diff.chemical.shifts.OUTCAR"
        outcar = Outcar(filename)
        core_contrib = outcar.data["cs_core_contribution"]
        assert core_contrib == {
            "Mg": -412.8248405,
            "C": -200.5098812,
            "O": -271.0766979,
        }

    def test_nmr_efg(self):
        filename = f"{TEST_DIR}/fixtures/nmr/efg/AlPO4/OUTCAR"
        outcar = Outcar(filename)
        expected_efg = [
            {"eta": 0.465, "nuclear_quadrupole_moment": 146.6, "cq": -5.573},
            {"eta": 0.465, "nuclear_quadrupole_moment": 146.6, "cq": -5.573},
            {"eta": 0.137, "nuclear_quadrupole_moment": 146.6, "cq": 6.327},
            {"eta": 0.137, "nuclear_quadrupole_moment": 146.6, "cq": 6.327},
            {"eta": 0.112, "nuclear_quadrupole_moment": 146.6, "cq": -7.453},
            {"eta": 0.112, "nuclear_quadrupole_moment": 146.6, "cq": -7.453},
            {"eta": 0.42, "nuclear_quadrupole_moment": 146.6, "cq": -5.58},
            {"eta": 0.42, "nuclear_quadrupole_moment": 146.6, "cq": -5.58},
        ]
        assert len(outcar.data["efg"][2:10]) == len(expected_efg)
        for e1, e2 in zip(outcar.data["efg"][2:10], expected_efg, strict=True):
            for k in e1:
                assert e1[k] == approx(e2[k], abs=1e-5)

        exepected_tensors = [
            [[11.11, 1.371, 2.652], [1.371, 3.635, -3.572], [2.652, -3.572, -14.746]],
            [[11.11, -1.371, 2.652], [-1.371, 3.635, 3.572], [2.652, 3.572, -14.746]],
            [[-3.098, 6.511, 7.732], [6.511, 1.419, 11.445], [7.732, 11.445, 1.678]],
            [
                [-3.098, -6.511, 7.732],
                [-6.511, 1.419, -11.445],
                [7.732, -11.445, 1.678],
            ],
            [
                [2.344, -10.775, -7.006],
                [-10.775, -7.152, -11.309],
                [-7.006, -11.309, 4.808],
            ],
            [
                [2.344, 10.775, -7.006],
                [10.775, -7.152, 11.309],
                [-7.006, 11.309, 4.808],
            ],
            [[2.404, -0.588, -6.83], [-0.588, 10.435, 3.159], [-6.83, 3.159, -12.839]],
            [[2.404, 0.588, -6.83], [0.588, 10.435, -3.159], [-6.83, -3.159, -12.839]],
        ]

        assert len(outcar.data["unsym_efg_tensor"][2:10]) == len(exepected_tensors)
        for e1, e2 in zip(outcar.data["unsym_efg_tensor"][2:10], exepected_tensors, strict=True):
            assert_allclose(e1, e2)

    def test_read_fermi_contact_shift(self):
        filepath = f"{VASP_OUT_DIR}/OUTCAR_fc"
        outcar = Outcar(filepath)
        outcar.read_fermi_contact_shift()
        assert outcar.data["fermi_contact_shift"]["fch"][0][0] == approx(-0.002)
        assert outcar.data["fermi_contact_shift"]["th"][0][0] == approx(-0.052)
        assert outcar.data["fermi_contact_shift"]["dh"][0][0] == approx(0.0)

    def test_drift(self):
        outcar = Outcar(f"{VASP_OUT_DIR}/OUTCAR.gz")
        assert len(outcar.drift) == 5
        assert np.sum(outcar.drift) == approx(0)

        outcar = Outcar(f"{VASP_OUT_DIR}/OUTCAR.CL.gz")
        assert len(outcar.drift) == 79
        assert np.sum(outcar.drift) == approx(0.448010)

    def test_electrostatic_potential(self):
        outcar = Outcar(f"{VASP_OUT_DIR}/OUTCAR.gz")
        assert outcar.ngf == [54, 30, 54]
        assert_allclose(outcar.sampling_radii, [0.9748, 0.9791, 0.7215])
        assert_allclose(
            outcar.electrostatic_potential,
            [-26.0704, -45.5046, -45.5046, -72.9539, -73.0621, -72.9539, -73.0621],
        )

    def test_mag_electrostatic_error(self):
        outcar = Outcar(f"{VASP_OUT_DIR}/OUTCAR.electrostaticerror.gz")
        # fmt: off
        assert outcar.electrostatic_potential == [
            -21.1667, -19.6865, -22.3983, -22.3307, -20.5213, -20.9292, -21.5063, -21.3554, -21.74,
            -21.7018, -20.3422, -20.6128, -21.4405, -21.0022, -21.975, -21.915, -21.0156, -21.9027,
            -22.3712, -21.5816, -21.8535, -20.5061, -22.2474, -22.1904, -22.2203, -20.1727, -21.1068,
            -20.1669, -22.1272, -21.3446, -82.4717, -83.035, -81.8289, -82.5957, -81.7813, -82.5011,
            -82.6098, -82.2885, -81.606, -99.1621, -99.3146, -99.1742, -99.4728, -100.2139, -99.852,
            -99.3575, -99.4135, -98.9092, -99.8867, -99.3707, -99.0794, -98.8376, -99.3656, -98.6474,
            -99.3264, -98.844, -99.074, -98.9354, -99.1643, -99.2412, -68.7667, -68.2528, -66.7326,
            -67.7113, -69.2228, -67.014, -69.1456, -67.3151, -68.2625, -67.6156, -69.8112, -68.9266,
            -67.8286, -69.3289, -68.7017, -67.2834, -68.4665, -68.0188, -67.7083, -69.7195, -67.4078,
            -67.9646, -68.584, -69.2387, -69.7822, -67.0701, -67.8236, -68.2468, -68.6533, -68.3218,
            -67.5923, -69.1266, -68.4615, -68.302, -67.999, -68.6709, -68.9973, -67.4147, -68.4463,
            -68.0899, -67.665, -69.6705, -68.6433, -68.4288, -66.9027, -67.3211, -68.604, -69.1299,
            -67.5565, -69.0845, -67.4289, -66.6864, -67.6484, -67.9783, -67.7661, -66.9797, -67.8007,
            -68.3194, -69.3671, -67.2708,
        ]
        # fmt: on

    def test_onsite_density_matrix(self):
        outcar = Outcar(f"{VASP_OUT_DIR}/OUTCAR.LinearResponseU.gz")
        matrices = outcar.data["onsite_density_matrices"]
        assert matrices[0][Spin.up][0][0] == approx(1.0227)
        assert len(matrices[0][Spin.up]) == 5
        assert len(matrices[0][Spin.up][0]) == 5
        assert "onsite_density_matrices" in outcar.as_dict()
        outcar = Outcar(f"{VASP_OUT_DIR}/OUTCAR_merged_numbers")
        matrices = outcar.data["onsite_density_matrices"]
        assert matrices[0][Spin.up][0][-1] == approx(0.0)
        assert len(matrices[0][Spin.up]) == 7
        assert len(matrices[0][Spin.up][0]) == 7
        assert "onsite_density_matrices" in outcar.as_dict()
        outcar = Outcar(f"{VASP_OUT_DIR}/OUTCAR_merged_numbers2")
        assert "onsite_density_matrices" in outcar.as_dict()

    def test_nbands(self):
        # Test VASP 5.2.11
        nbands = Outcar(f"{VASP_OUT_DIR}/OUTCAR.gz").data["nbands"]
        assert nbands == 33
        assert isinstance(nbands, int)

        # Test VASP 5.4.4
        assert Outcar(f"{VASP_OUT_DIR}/OUTCAR.LOPTICS.vasp544").data["nbands"] == 128

        # Test VASP 6.3.0
        assert Outcar(f"{VASP_OUT_DIR}/OUTCAR_vasp_6.3.gz").data["nbands"] == 64

        # Test NBANDS set by user but overridden by VASP
        # VASP 6.3.2
        assert Outcar(f"{VASP_OUT_DIR}/OUTCAR.nbands_overridden.gz").data["nbands"] == 32

    def test_nplwvs(self):
        outcar = Outcar(f"{VASP_OUT_DIR}/OUTCAR.gz")
        assert outcar.data["nplwv"] == [[34560]]
        # fmt: off
        assert outcar.data["nplwvs_at_kpoints"] == [
            1719, 1714, 1722, 1728, 1722, 1726, 1722, 1720, 1717, 1724, 1715, 1724, 1726,
            1724, 1728, 1715, 1722, 1715, 1726, 1730, 1730, 1715, 1716, 1729, 1727, 1723,
            1721, 1712, 1723, 1719, 1717, 1717, 1724, 1719, 1719, 1727, 1726, 1730, 1719,
            1720, 1718, 1717, 1722, 1719, 1709, 1714, 1724, 1726, 1718, 1713, 1720, 1713,
            1711, 1713, 1715, 1717, 1728, 1726, 1712, 1722, 1714, 1713, 1717, 1714, 1714,
            1717, 1712, 1710, 1721, 1722, 1724, 1720, 1726, 1719, 1722, 1714,
        ]
        # fmt: on
        outcar = Outcar(f"{VASP_OUT_DIR}/OUTCAR.CL.gz")
        assert outcar.data["nplwv"] == [[None]]
        assert outcar.data["nplwvs_at_kpoints"] == [85687]

    def test_serial_compilation(self):
        outcar = Outcar(f"{VASP_OUT_DIR}/OUTCAR.serial.gz")
        assert outcar.data["nplwv"] == [[74088]]
        assert outcar.data["nplwvs_at_kpoints"] == [4418, 4390, 4421, 4404]

    def test_vasp620_format(self):
        filepath = f"{VASP_OUT_DIR}/OUTCAR.vasp.6.2.0.gz"
        outcar = Outcar(filepath)
        assert outcar.run_stats["Average memory used (kb)"] is None

        filepath = f"{VASP_OUT_DIR}/OUTCAR.vasp.6.2.1.mpi.gz"
        outcar = Outcar(filepath)
        assert outcar.run_stats["cores"] == 64

    def test_energies(self):
        # VASP 5.2.1
        outcar = Outcar(f"{VASP_OUT_DIR}/OUTCAR.etest1.gz")
        assert outcar.final_energy == approx(-11.18981538)
        assert outcar.final_energy_wo_entrp == approx(-11.13480014)
        assert outcar.final_fr_energy == approx(-11.21732300)

        # VASP 6.2.1
        outcar = Outcar(f"{VASP_OUT_DIR}/OUTCAR.etest2.gz")
        assert outcar.final_energy == approx(-11.18986774)
        assert outcar.final_energy_wo_entrp == approx(-11.13485250)
        assert outcar.final_fr_energy == approx(-11.21737536)

        # VASP 5.2.1
        outcar = Outcar(f"{VASP_OUT_DIR}/OUTCAR.etest3.gz")
        assert outcar.final_energy == approx(-15.89355325)
        assert outcar.final_energy_wo_entrp == approx(-15.83853800)
        assert outcar.final_fr_energy == approx(-15.92106087)

        # VASP 6.2.1
        outcar = Outcar(f"{VASP_OUT_DIR}/OUTCAR.etest4.gz")
        assert outcar.final_energy == approx(-15.89364691)
        assert outcar.final_energy_wo_entrp == approx(-15.83863167)
        assert outcar.final_fr_energy == approx(-15.92115453)

    def test_read_table_pattern(self):
        outcar = Outcar(f"{VASP_OUT_DIR}/OUTCAR.gz")

        header_pattern = r"\(the norm of the test charge is\s+[\.\-\d]+\)"
        table_pattern = r"((?:\s+\d+\s*[\.\-\d]+)+)"
        footer_pattern = r"\s+E-fermi :"

        pots = outcar.read_table_pattern(header_pattern, table_pattern, footer_pattern, last_one_only=True)
        ref_last = [
            ["       1 -26.0704       2 -45.5046       3 -45.5046       4 -72.9539       5 -73.0621"],
            ["       6 -72.9539       7 -73.0621"],
        ]
        assert pots == ref_last

        pots = outcar.read_table_pattern(
            header_pattern,
            table_pattern,
            footer_pattern,
            last_one_only=False,
            first_one_only=True,
        )
        ref_first = [
            ["       1 -26.1149       2 -45.5359       3 -45.5359       4 -72.9831       5 -73.1068"],
            ["       6 -72.9831       7 -73.1068"],
        ]
        assert pots == ref_first

        with pytest.raises(
            ValueError,
            match="last_one_only and first_one_only options are incompatible",
        ):
            outcar.read_table_pattern(
                header_pattern,
                table_pattern,
                footer_pattern,
                last_one_only=True,
                first_one_only=True,
            )


class TestBSVasprun(MatSciTest):
    def test_get_band_structure(self):
        filepath = f"{VASP_OUT_DIR}/vasprun_Si_bands.xml.gz"
        vasprun = BSVasprun(filepath, parse_potcar_file=False)
        bs = vasprun.get_band_structure(kpoints_filename=f"{VASP_IN_DIR}/KPOINTS_Si_bands")
        cbm = bs.get_cbm()
        vbm = bs.get_vbm()
        assert cbm["kpoint_index"] == [13], "wrong cbm kpoint index"
        assert cbm["energy"] == approx(6.2301), "wrong cbm energy"
        assert cbm["band_index"] == {Spin.up: [4], Spin.down: [4]}, "wrong cbm bands"
        assert vbm["kpoint_index"] == [0, 63, 64]
        assert vbm["energy"] == approx(5.6158), "wrong vbm energy"
        assert vbm["band_index"] == {
            Spin.up: [1, 2, 3],
            Spin.down: [1, 2, 3],
        }, "wrong vbm bands"
        assert vbm["kpoint"].label == "\\Gamma", "wrong vbm label"
        assert cbm["kpoint"].label is None, "wrong cbm label"
        dct = vasprun.as_dict()
        assert "eigenvalues" in dct["output"]

    def test_kpoints_opt(self):
        vasp_run = BSVasprun(
            f"{TEST_DIR}/fixtures/kpoints_opt/vasprun.xml.gz", parse_potcar_file=False, parse_projected_eigen=True
        )
        bs = vasp_run.get_band_structure(f"{TEST_DIR}/fixtures/kpoints_opt/KPOINTS_OPT")
        assert isinstance(bs, BandStructureSymmLine)
        cbm = bs.get_cbm()
        vbm = bs.get_vbm()
        assert cbm["kpoint_index"] == [38], "wrong cbm kpoint index"
        assert cbm["energy"] == approx(6.4394), "wrong cbm energy"
        assert cbm["band_index"] == {Spin.up: [16], Spin.down: [16]}, "wrong cbm bands"
        # Strangely, when I call with parse_projected_eigen, it gives empty Spin.down,
        # but without parse_projected_eigen it does not give it.
        # So at one point it called the empty key.
        assert vbm["kpoint_index"] == [0, 39, 40]
        assert vbm["energy"] == approx(5.7562), "wrong vbm energy"
        assert vbm["band_index"] == {
            Spin.down: [13, 14, 15],
            Spin.up: [13, 14, 15],
        }, "wrong vbm bands"
        assert vbm["kpoint"].label == "\\Gamma", "wrong vbm label"
        assert cbm["kpoint"].label is None, "wrong cbm label"
        # Test projection
        projected = bs.get_projection_on_elements()
        assert np.isnan(projected[Spin.up][0][0]["Si"])
        # Due to some error in my VASP, the transcription of PROCAR_OPT into
        # vasprun.xml is filled to the brim with errors in the projections.
        # At some point we might get a healthier vasprun.xml, but the point here
        # is to test the parser, not VASP.
        projected = bs.get_projections_on_elements_and_orbitals({"Si": ["s"]})
        assert projected[Spin.up][0][58]["Si"]["s"] == approx(-0.0271)
        vrun_dct = vasp_run.as_dict()
        assert {*vrun_dct["output"]} >= {"eigenvalues", "eigenvalues_kpoints_opt"}


class TestOszicar(MatSciTest):
    def test_init(self):
        fpath = f"{VASP_OUT_DIR}/OSZICAR"
        oszicar = Oszicar(fpath)
        assert len(oszicar.electronic_steps) == len(oszicar.ionic_steps)
        assert len(oszicar.all_energies) == 60
        assert oszicar.final_energy == approx(-526.63928)
        assert set(oszicar.ionic_steps[-1]) == {"F", "E0", "dE", "mag"}

    def test_static(self):
        fpath = f"{TEST_DIR}/fixtures/static_silicon/OSZICAR"
        oszicar = Oszicar(fpath)
        assert oszicar.final_energy == approx(-10.645278)
        assert set(oszicar.ionic_steps[-1]) == {"F", "E0", "dE", "mag"}


class TestGetBandStructureFromVaspMultipleBranches:
    def test_read_multi_branches(self):
        """TODO: This functionality still needs a test."""

    def test_missing_vasprun_in_branch_dir(self):
        """Test vasprun.xml missing from branch_*."""
        os.makedirs("no_vasp/branch_0", exist_ok=False)

        with pytest.raises(FileNotFoundError, match="cannot find vasprun.xml in directory"):
            get_band_structure_from_vasp_multiple_branches("no_vasp")

    def test_no_branch_head(self):
        """Test branch_0 is missing and read dir_name/vasprun.xml directly."""

        copyfile(f"{VASP_OUT_DIR}/vasprun.force_hybrid_like_calc.xml.gz", "./vasprun.xml.gz")
        decompress_file("./vasprun.xml.gz")

        with pytest.warns(DeprecationWarning, match="no branch dir found, reading directly from"):
            bs = get_band_structure_from_vasp_multiple_branches(".")
        assert isinstance(bs, BandStructure)

    def test_cannot_read_anything(self):
        """Test no branch_0/, no dir_name/vasprun.xml, no vasprun.xml at all."""
        with pytest.raises(FileNotFoundError, match="failed to find any vasprun.xml in selected"):
            get_band_structure_from_vasp_multiple_branches(".")


class TestLocpot(MatSciTest):
    def test_init(self):
        filepath = f"{VASP_OUT_DIR}/LOCPOT.gz"
        locpot = Locpot.from_file(filepath)
        assert sum(locpot.get_average_along_axis(0)) == approx(-217.05226954)
        assert locpot.get_axis_grid(0)[-1] == approx(2.87629, abs=1e-2)
        assert locpot.get_axis_grid(1)[-1] == approx(2.87629, abs=1e-2)
        assert locpot.get_axis_grid(2)[-1] == approx(2.87629, abs=1e-2)

        # make sure locpot constructor works with data_aug=None
        poscar, data, _data_aug = Locpot.parse_file(filepath)
        l2 = Locpot(poscar=poscar, data=data, data_aug=None)
        assert l2.data_aug == {}

    def test_vasp_6x_style(self):
        filepath = f"{VASP_OUT_DIR}/LOCPOT.vasp642.gz"
        locpot = Locpot.from_file(filepath)
        assert locpot.dim == (2, 2, 5)
        assert {str(ele) for ele in locpot.structure.composition} == {"Mg", "Si"}


class TestChgcar(MatSciTest):
    @classmethod
    def setup_class(cls):
        filepath = f"{VASP_OUT_DIR}/CHGCAR.nospin.gz"
        cls.chgcar_no_spin = Chgcar.from_file(filepath)

        filepath = f"{VASP_OUT_DIR}/CHGCAR.spin.gz"
        cls.chgcar_spin = Chgcar.from_file(filepath)

        filepath = f"{VASP_OUT_DIR}/CHGCAR.Fe3O4.gz"
        cls.chgcar_fe3o4 = Chgcar.from_file(filepath)

        filepath = f"{VASP_OUT_DIR}/CHGCAR.NiO_SOC.gz"
        cls.chgcar_NiO_soc = Chgcar.from_file(filepath)

    def test_init(self):
        assert self.chgcar_no_spin.get_integrated_diff(0, 2)[0, 1] == approx(0)
        assert self.chgcar_spin.get_integrated_diff(0, 1)[0, 1] == approx(-0.0043896932237534022)
        # test sum
        chgcar = self.chgcar_spin + self.chgcar_spin
        assert chgcar.get_integrated_diff(0, 1)[0, 1] == approx(-0.0043896932237534022 * 2)

        chgcar = sum([self.chgcar_spin, self.chgcar_spin])
        assert chgcar.get_integrated_diff(0, 1)[0, 1] == approx(-0.0043896932237534022 * 2)

        chgcar = self.chgcar_spin - self.chgcar_spin
        assert chgcar.get_integrated_diff(0, 1)[0, 1] == approx(0)

        expected = [
            1.56472768,
            3.25985108,
            3.49205728,
            3.66275028,
            3.8045896,
            5.10813352,
        ]
        actual = self.chgcar_fe3o4.get_integrated_diff(0, 3, 6)
        assert_allclose(actual[:, 1], expected)

    def test_write(self):
        self.chgcar_spin.write_file(out_path := f"{self.tmp_path}/CHGCAR_pmg")
        with open(out_path, encoding="utf-8") as file:
            for idx, line in enumerate(file):
                if idx in (22130, 44255):
                    assert line == "augmentation occupancies   1  15\n"

    def test_soc_chgcar(self):
        assert set(self.chgcar_NiO_soc.data) == {
            "total",
            "diff_x",
            "diff_y",
            "diff_z",
            "diff",
        }
        assert self.chgcar_NiO_soc.is_soc
        assert self.chgcar_NiO_soc.data["diff"].shape == self.chgcar_NiO_soc.data["diff_y"].shape

        # check our construction of chg.data['diff'] makes sense
        # this has been checked visually too and seems reasonable
        assert abs(self.chgcar_NiO_soc.data["diff"][0][0][0]) == np.linalg.norm(
            [self.chgcar_NiO_soc.data[f"diff_{key}"][0][0][0] for key in "xyz"]
        )

        # and that the net magnetization is about zero
        # note: we get ~ 0.08 here, seems a little high compared to
        # vasp output, but might be due to chgcar limitations?
        assert self.chgcar_NiO_soc.net_magnetization == approx(0.0, abs=1e-0)

        self.chgcar_NiO_soc.write_file(out_path := f"{self.tmp_path}/CHGCAR_pmg_soc")
        chg_from_file = Chgcar.from_file(out_path)
        assert chg_from_file.is_soc

    @pytest.mark.skipif(h5py is None, reason="h5py required for HDF5 support.")
    def test_hdf5(self):
        chgcar = Chgcar.from_file(f"{VASP_OUT_DIR}/CHGCAR.NiO_SOC.gz")
        chgcar.to_hdf5(out_path := f"{self.tmp_path}/chgcar_test.hdf5")

        with h5py.File(out_path, mode="r") as dct:
            assert_allclose(dct["vdata"]["total"], chgcar.data["total"])
            assert_allclose(dct["vdata"]["diff"], chgcar.data["diff"])
            assert_allclose(dct["lattice"], chgcar.structure.lattice.matrix)
            assert_allclose(dct["fcoords"], chgcar.structure.frac_coords)
            for z in dct["Z"]:
                assert z in [Element.Ni.Z, Element.O.Z]

            for sp in dct["species"]:
                assert sp in [b"Ni", b"O"]

        chgcar2 = Chgcar.from_hdf5(out_path)
        assert_allclose(chgcar2.data["total"], chgcar.data["total"])

    def test_spin_data(self):
        for v in self.chgcar_spin.spin_data.values():
            assert v.shape == (48, 48, 48)

    def test_add(self):
        chgcar_sum = self.chgcar_spin + self.chgcar_spin
        assert_allclose(chgcar_sum.data["total"], self.chgcar_spin.data["total"] * 2)
        chgcar_copy = self.chgcar_spin.copy()
        chgcar_copy.structure = self.get_structure("Li2O")
        with pytest.warns(
            UserWarning,
            match="Structures are different. Make sure you know what you are doing...",
        ) as warns:
            chgcar_sum = chgcar_copy + self.chgcar_spin
        assert len(warns) == 1
        with pytest.raises(
            ValueError,
            match=r"operands could not be broadcast together with shapes \(48,48,48\) \(72,72,72\)",
        ):
            _ = self.chgcar_spin + self.chgcar_fe3o4
        with pytest.raises(
            ValueError,
            match="Data have different keys! Maybe one is spin-polarized and the other is not",
        ):
            _ = self.chgcar_spin + self.chgcar_no_spin

    def test_as_dict_and_from_dict(self):
        dct = self.chgcar_NiO_soc.as_dict()
        chgcar_from_dict = Chgcar.from_dict(dct)
        assert_allclose(self.chgcar_NiO_soc.data["total"], chgcar_from_dict.data["total"])
        assert_allclose(
            self.chgcar_NiO_soc.structure.lattice.matrix,
            chgcar_from_dict.structure.lattice.matrix,
        )


class TestAeccars(MatSciTest):
    # https://github.com/materialsproject/pymatgen/pull/3343
    def test_read_write_file(self):
        aeccar0_test = Chgcar.from_file(f"{TEST_FILES_DIR}/command_line/bader/AECCAR0.gz")
        aeccar0_outpath = f"{self.tmp_path}/AECCAR0_test"
        aeccar0_test.write_file(aeccar0_outpath)
        aeccar0_read = Chgcar.from_file(aeccar0_outpath)
        assert_allclose(aeccar0_test.data["total"], aeccar0_read.data["total"])

        aeccar2 = Chgcar.from_file(f"{TEST_FILES_DIR}/command_line/bader/AECCAR2.gz")
        aeccar2_outpath = f"{self.tmp_path}/AECCAR2_test"
        aeccar2.write_file(aeccar2_outpath)
        aeccar2_read = Chgcar.from_file(aeccar2_outpath)
        assert_allclose(aeccar2.data["total"], aeccar2_read.data["total"])


class TestElfcar(MatSciTest):
    def test_init(self):
        elfcar = Elfcar.from_file(f"{VASP_OUT_DIR}/ELFCAR.gz")
        assert np.mean(elfcar.data["total"]) == approx(0.19076207645194002)
        assert np.mean(elfcar.data["diff"]) == approx(0.19076046677910055)
        reconstituted = Elfcar.from_dict(elfcar.as_dict())
        assert elfcar.data == reconstituted.data
        assert elfcar.poscar.structure == reconstituted.poscar.structure

    def test_alpha(self):
        elfcar = Elfcar.from_file(f"{VASP_OUT_DIR}/ELFCAR.gz")
        alpha = elfcar.get_alpha()
        assert np.median(alpha.data["total"]) == approx(2.936678808979031)

    def test_interpolation(self):
        elfcar = Elfcar.from_file(f"{VASP_OUT_DIR}/ELFCAR.gz")
        assert elfcar.value_at(0.4, 0.5, 0.6) == approx(0.0918471)
        assert len(elfcar.linear_slice([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])) == 100


class TestProcar(MatSciTest):
    def test_init(self):
        filepath = f"{VASP_OUT_DIR}/PROCAR.simple"
        procar = Procar(filepath)
        assert procar.get_occupation(0, "d")[Spin.up] == approx(0)
        assert procar.get_occupation(0, "s")[Spin.up] == approx(0.35381249999999997)
        assert procar.get_occupation(0, "p")[Spin.up] == approx(1.19540625)
        with pytest.raises(ValueError, match="'m' is not in list"):
            procar.get_occupation(1, "m")
        assert procar.nbands == 10
        assert procar.nkpoints == 10
        assert procar.nions == 3
        filepath = f"{VASP_OUT_DIR}/PROCAR.gz"
        procar = Procar(filepath)
        assert procar.get_occupation(0, "dxy")[Spin.up] == approx(0.96214813853000025)
        assert procar.get_occupation(0, "dxy")[Spin.down] == approx(0.85796295426000124)

    def test_soc_procar(self):
        filepath = f"{VASP_OUT_DIR}/PROCAR.SOC.gz"
        procar = Procar(filepath)
        assert procar.nions == 4
        assert procar.nkpoints == 25
        assert procar.nspins == 1
        assert procar.is_soc  # SOC PROCAR
        nb = procar.nbands
        nk = procar.nkpoints
        assert procar.eigenvalues[Spin.up].shape == (nk, nb)
        assert procar.kpoints.shape == (nk, 3)
        assert len(procar.weights) == nk
        assert_allclose(procar.kpoints[0][0], 0.0)
        assert procar.occupancies[Spin.up].shape == (nk, nb)

        # spot check some values:
        assert procar.data[Spin.up][0, 1, 1, 0] == approx(0.095)
        assert procar.data[Spin.up][0, 1, 1, 1] == approx(0)

        assert procar.xyz_data["x"][0, 1, 1, 0] == approx(-0.063)
        assert procar.xyz_data["z"][0, 1, 1, 1] == approx(0)

    def test_multiple_procars(self):
        filepaths = [
            f"{VASP_OUT_DIR}/PROCAR.split1.gz",
            f"{VASP_OUT_DIR}/PROCAR.split2.gz",
        ]
        procar = Procar(filepaths)
        assert procar.nions == 4
        assert procar.nkpoints == 96  # 96 overall, 48 in first PROCAR, 96 in second (48 duplicates)
        assert procar.nspins == 1  # SOC PROCAR, also with LORBIT = 14
        assert procar.is_soc  # SOC PROCAR
        nb = procar.nbands
        nk = procar.nkpoints
        assert procar.eigenvalues[Spin.up].shape == (nk, nb)
        assert procar.kpoints.shape == (nk, 3)
        assert len(procar.weights) == nk
        assert procar.occupancies[Spin.up].shape == (nk, nb)

        # spot check some values:
        assert procar.data[Spin.up][0, 1, 1, 0] == approx(0.094)
        assert procar.data[Spin.up][0, 1, 1, 1] == approx(0)

        assert procar.xyz_data["x"][0, 1, 1, 0] == approx(0)
        assert procar.xyz_data["z"][0, 1, 1, 1] == approx(0)

        assert procar.phase_factors[Spin.up][0, 1, 0, 0] == approx(-0.159 + 0.295j)

    def test_phase_factors(self):
        filepath = f"{VASP_OUT_DIR}/PROCAR.phase.gz"
        procar = Procar(filepath)
        assert procar.phase_factors[Spin.up][0, 0, 0, 0] == approx(-0.746 + 0.099j)
        assert procar.phase_factors[Spin.down][0, 0, 0, 0] == approx(0.372 - 0.654j)

        # Two Li should have same phase factor.
        assert procar.phase_factors[Spin.up][0, 0, 0, 0] == approx(procar.phase_factors[Spin.up][0, 0, 1, 0])
        assert procar.phase_factors[Spin.up][0, 0, 2, 0] == approx(-0.053 + 0.007j)
        assert procar.phase_factors[Spin.down][0, 0, 2, 0] == approx(0.027 - 0.047j)

        # new style phase factors (VASP 5.4.4+)
        filepath = f"{VASP_OUT_DIR}/PROCAR.new_format_5.4.4.gz"
        procar = Procar(filepath)
        assert procar.phase_factors[Spin.up][0, 0, 0, 0] == approx(-0.13 + 0.199j)

    def test_get_projection_on_elements(self):
        filepath = f"{VASP_OUT_DIR}/PROCAR.simple"
        procar = Procar(filepath)
        struct = Structure(
            Lattice.cubic(3.0),
            ["Li", "Na", "K"],
            [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25], [0.75, 0.75, 0.75]],
        )
        dct = procar.get_projection_on_elements(struct)
        assert dct[Spin.up][2][2] == approx({"Na": 0.042, "K": 0.646, "Li": 0.042})
        # https://github.com/materialsproject/pymatgen/pull/3261
        struct.replace_species({"K": "Na"})
        d2 = procar.get_projection_on_elements(struct)
        assert d2[Spin.up][2][2] == approx({"Na": 0.688, "Li": 0.042})


class TestXdatcar:
    def test_init(self):
        filepath = f"{VASP_OUT_DIR}/XDATCAR_4"
        xdatcar = Xdatcar(filepath)
        structures = xdatcar.structures
        assert len(structures) == 4
        for struct in structures:
            assert struct.formula == "Li2 O1"

        filepath = f"{VASP_OUT_DIR}/XDATCAR_5"
        xdatcar = Xdatcar(filepath)
        structures = xdatcar.structures
        assert len(structures) == 4
        for struct in structures:
            assert struct.formula == "Li2 O1"

        xdatcar.concatenate(f"{VASP_OUT_DIR}/XDATCAR_4")
        assert len(xdatcar.structures) == 8
        assert xdatcar.get_str() is not None

        filepath = f"{VASP_OUT_DIR}/XDATCAR_6"
        xdatcar = Xdatcar(filepath)
        structures = xdatcar.structures

        assert structures[0].lattice != structures[-1].lattice

        xdatcar = Xdatcar(f"{VASP_OUT_DIR}/XDATCAR_monatomic.gz")
        assert len(xdatcar.structures) == 10
        assert len(xdatcar) == len(xdatcar.structures)

        idxs = [1, 2, 3]
        assert xdatcar[idxs] == [xdatcar.structures[i] for i in idxs]
        assert all(isinstance(struct, Structure) for struct in xdatcar)

        assert all(len(structure.composition) == 1 for structure in xdatcar.structures)

    def test_bad_format(self):
        # ensure XDATCAR can be read even when formatting is poor
        xdatcar = Xdatcar(f"{VASP_OUT_DIR}/XDATCAR.bad_fmt.gz")
        assert isinstance(xdatcar, Xdatcar)


class TestDynmat:
    def test_init(self):
        filepath = f"{VASP_OUT_DIR}/DYNMAT"
        dct = Dynmat(filepath)
        assert dct.nspecs == 2
        assert dct.natoms == 6
        assert dct.ndisps == 3
        assert_allclose(dct.masses, [63.546, 196.966])
        assert 4 in dct.data
        assert 2 in dct.data[4]
        assert_allclose(dct.data[4][2]["dispvec"], [0.0, 0.05, 0.0])
        assert_allclose(dct.data[4][2]["dynmat"][3], [0.055046, -0.298080, 0.0])
        # TODO: test get_phonon_frequencies once cross-checked


class TestWavecar(MatSciTest):
    def setup_method(self):
        latt_mat = np.array(np.eye(3) * 10, dtype=float)  # lattice vectors
        self.vol = np.dot(latt_mat[0, :], np.cross(latt_mat[1, :], latt_mat[2, :]))  # unit cell volume
        # reciprocal lattice vectors
        recip_latt_mat = [
            np.cross(latt_mat[1, :], latt_mat[2, :]),
            np.cross(latt_mat[2, :], latt_mat[0, :]),
            np.cross(latt_mat[0, :], latt_mat[1, :]),
        ]
        self.recip_latt_mat = 2 * np.pi * np.array(recip_latt_mat) / self.vol
        self.latt_mat = latt_mat
        self.wavecar = Wavecar(f"{VASP_OUT_DIR}/WAVECAR.N2")
        self.wH2 = Wavecar(f"{VASP_OUT_DIR}/WAVECAR.H2_low_symm")
        self.wH2_gamma = Wavecar(f"{VASP_OUT_DIR}/WAVECAR.H2_low_symm.gamma")
        self.w_ncl = Wavecar(f"{VASP_OUT_DIR}/WAVECAR.H2.ncl")
        self.w_frac_encut = Wavecar(f"{VASP_OUT_DIR}/WAVECAR.frac_encut")

    def test_standard(self):
        wavecar = self.wavecar
        a = np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]])
        vol = np.dot(a[0, :], np.cross(a[1, :], a[2, :]))
        b = np.array(
            [
                np.cross(a[1, :], a[2, :]),
                np.cross(a[2, :], a[0, :]),
                np.cross(a[0, :], a[1, :]),
            ]
        )
        b = 2 * np.pi * b / vol

        assert wavecar.filename == f"{VASP_OUT_DIR}/WAVECAR.N2"
        assert wavecar.efermi == approx(-5.7232, abs=1e-4)
        assert wavecar.encut == approx(25.0)
        assert wavecar.nb == 9
        assert wavecar.nk == 1
        assert_allclose(wavecar.a, a)
        assert_allclose(wavecar.b, b)
        assert wavecar.vol == approx(vol)
        assert len(wavecar.kpoints) == wavecar.nk
        assert len(wavecar.coeffs) == wavecar.nk
        assert len(wavecar.coeffs[0]) == wavecar.nb
        assert len(wavecar.band_energy) == wavecar.nk
        assert wavecar.band_energy[0].shape == (wavecar.nb, 3)
        assert len(wavecar.Gpoints[0]) <= 257
        for k in range(wavecar.nk):
            for b in range(wavecar.nb):
                assert len(wavecar.coeffs[k][b]) == len(wavecar.Gpoints[k])

        # Test WAVECAR with fractional ENCUT
        assert self.w_frac_encut.encut == approx(100.5)

        # Test malformed WAVECARs
        with pytest.raises(ValueError, match="Invalid rtag=.+, must be one of"):
            Wavecar(f"{VASP_OUT_DIR}/WAVECAR.N2.malformed")

        with pytest.raises(ValueError, match="invalid vasp_type='poop'"):
            Wavecar(f"{VASP_OUT_DIR}/WAVECAR.N2", vasp_type="poop")

        with pytest.raises(
            ValueError,
            match=r"Incorrect vasp_type='g'. Please open an issue if you are certain",
        ):
            Wavecar(f"{VASP_OUT_DIR}/WAVECAR.N2", vasp_type="g")

        with pytest.raises(ValueError, match=r"cannot reshape array of size 257 into shape \(2,128\)"):
            Wavecar(f"{VASP_OUT_DIR}/WAVECAR.N2", vasp_type="n")

        saved_stdout = sys.stdout
        try:
            out = StringIO()
            sys.stdout = out
            Wavecar(f"{VASP_OUT_DIR}/WAVECAR.N2", verbose=True)
            assert out.getvalue().strip() != ""
        finally:
            sys.stdout = saved_stdout

    def test_n2_45210(self):
        wavecar = Wavecar(f"{VASP_OUT_DIR}/WAVECAR.N2.45210")
        assert wavecar.filename == f"{VASP_OUT_DIR}/WAVECAR.N2.45210"
        assert wavecar.efermi == approx(-5.7232, abs=1e-4)
        assert wavecar.encut == approx(25.0)
        assert wavecar.nb == 9
        assert wavecar.nk == 1
        assert_allclose(wavecar.a, self.latt_mat)
        assert_allclose(wavecar.b, self.recip_latt_mat)
        assert wavecar.vol == approx(self.vol)
        assert len(wavecar.kpoints) == wavecar.nk
        assert len(wavecar.coeffs) == wavecar.nk
        assert len(wavecar.coeffs[0]) == wavecar.nb
        assert len(wavecar.band_energy) == wavecar.nk
        assert wavecar.band_energy[0].shape == (wavecar.nb, 3)
        assert len(wavecar.Gpoints[0]) <= 257

    def test_n2_spin(self):
        w = Wavecar(f"{VASP_OUT_DIR}/WAVECAR.N2.spin")
        assert len(w.coeffs) == 2
        assert len(w.band_energy) == 2
        assert len(w.kpoints) == w.nk
        assert len(w.Gpoints) == w.nk
        assert len(w.coeffs[0][0]) == w.nb
        assert len(w.band_energy[0]) == w.nk

        orig_gen_g_points = Wavecar._generate_G_points
        try:
            Wavecar._generate_G_points = lambda _x, _y, gamma: []
            with pytest.raises(ValueError, match=r"not enough values to unpack \(expected 3, got 0\)"):
                Wavecar(f"{VASP_OUT_DIR}/WAVECAR.N2")
        finally:
            Wavecar._generate_G_points = orig_gen_g_points

    def test_generate_nbmax(self):
        self.wavecar._generate_nbmax()
        assert self.wavecar._nbmax.tolist() == [5, 5, 5]

    def test_generate_g_points(self):
        for k in range(self.wavecar.nk):
            kp = self.wavecar.kpoints[k]
            assert len(self.wavecar._generate_G_points(kp)) <= 257

    def test_evaluate_wavefunc(self):
        self.wavecar.Gpoints.append(np.array([0, 0, 0]))
        self.wavecar.kpoints.append(np.array([0, 0, 0]))
        self.wavecar.coeffs.append([[1 + 1j]])
        assert self.wavecar.evaluate_wavefunc(-1, -1, [0, 0, 0]) == approx((1 + 1j) / np.sqrt(self.vol), abs=1e-4)
        assert self.wavecar.evaluate_wavefunc(0, 0, [0, 0, 0]) == approx(
            np.sum(self.wavecar.coeffs[0][0]) / np.sqrt(self.vol), abs=1e-4
        )
        w = Wavecar(f"{VASP_OUT_DIR}/WAVECAR.N2.spin")
        w.Gpoints.append(np.array([0, 0, 0]))
        w.kpoints.append(np.array([0, 0, 0]))
        w.coeffs[0].append([[1 + 1j]])
        assert w.evaluate_wavefunc(-1, -1, [0, 0, 0]) == approx((1 + 1j) / np.sqrt(self.vol), abs=1e-4)

    def test_fft_mesh_basic(self):
        mesh = self.wavecar.fft_mesh(0, 5)
        ind = np.argmax(np.abs(mesh))
        assert np.unravel_index(ind, mesh.shape) == (14, 1, 1)
        assert mesh[tuple((self.wavecar.ng / 2).astype(int))] == 0j
        mesh = self.wavecar.fft_mesh(0, 5, shift=False)
        ind = np.argmax(np.abs(mesh))
        assert np.unravel_index(ind, mesh.shape) == (6, 8, 8)
        assert mesh[0, 0, 0] == 0j

    def test_fft_mesh_advanced(self):
        ik = ib = 0
        mesh = self.wH2.fft_mesh(ik, ib)
        mesh_gamma = self.wH2_gamma.fft_mesh(ik, ib)
        mesh_ncl = self.w_ncl.fft_mesh(ik, ib)

        # check equality of plane-wave coefficients
        ind_max = np.unravel_index(np.argmax(np.abs(mesh)), mesh.shape)
        phase = mesh[ind_max] / mesh_gamma[ind_max]
        assert np.max(np.abs(mesh - phase * mesh_gamma)) <= 1.0e-6

        # transform to real space for further checking
        mesh = np.fft.ifftn(mesh)
        mesh_gamma = np.fft.ifftn(mesh_gamma)
        mesh_ncl = np.fft.ifftn(mesh_ncl)

        # check equality in real space for regular vs. gamma only
        ind_max = np.unravel_index(np.argmax(np.abs(mesh)), mesh.shape)
        phase = mesh[ind_max] / mesh_gamma[ind_max]
        assert np.max(np.abs(mesh - phase * mesh_gamma)) <= 1.0e-6

        # spot check some points in real space
        p1 = (
            int(mesh.shape[0] / 2),
            int(mesh.shape[1] / 2) - 1,
            int(mesh.shape[2] / 2) - 2,
        )
        p2 = (p1[0] + 1, p1[1], p1[2])
        c = np.diag((5, 4, 6))  # this needs to match POSCAR, which we don't have
        r1 = np.dot(np.array(p1) / mesh.shape, c)
        r2 = np.dot(np.array(p2) / mesh.shape, c)

        # check equality of FFT and slow FT for regular mesh (ratio, to account for normalization)
        v1 = self.wH2.evaluate_wavefunc(ik, ib, r1)
        v2 = self.wH2.evaluate_wavefunc(ik, ib, r2)
        assert np.abs(mesh[p1]) / np.abs(mesh[p2]) == approx(np.abs(v1) / np.abs(v2), abs=1e-6)

        # spot check one value that we happen to know from reference run
        assert v1 == approx(-0.01947068011502887 + 0.23340228099620275j, abs=1e-8)

        # check equality of FFT and slow FT for gamma-only mesh (ratio again)
        v1_gamma = self.wH2_gamma.evaluate_wavefunc(ik, ib, r1)
        v2_gamma = self.wH2_gamma.evaluate_wavefunc(ik, ib, r2)
        assert np.abs(mesh_gamma[p1]) / np.abs(mesh_gamma[p2]) == approx(np.abs(v1_gamma) / np.abs(v2_gamma), abs=1e-6)

        # check equality of FFT and slow FT for ncl mesh (ratio again)
        v1_ncl = self.w_ncl.evaluate_wavefunc(ik, ib, r1)
        v2_ncl = self.w_ncl.evaluate_wavefunc(ik, ib, r2)
        assert np.abs(mesh_ncl[p1]) / np.abs(mesh_ncl[p2]) == approx(np.abs(v1_ncl) / np.abs(v2_ncl), abs=1e-6)

    def test_get_parchg(self):
        poscar = Poscar.from_file(f"{VASP_IN_DIR}/POSCAR")

        wavecar = self.wavecar
        chgcar = wavecar.get_parchg(poscar, 0, 0, spin=0, phase=False)
        assert "total" in chgcar.data
        assert "diff" not in chgcar.data
        assert chgcar.data["total"].size == np.prod(wavecar.ng * 2)
        assert np.all(chgcar.data["total"] > 0.0)

        chgcar = wavecar.get_parchg(poscar, 0, 0, spin=0, phase=True)
        assert "total" in chgcar.data
        assert "diff" not in chgcar.data
        assert chgcar.data["total"].size == np.prod(wavecar.ng * 2)
        assert not np.all(chgcar.data["total"] > 0.0)

        wavecar = Wavecar(f"{VASP_OUT_DIR}/WAVECAR.N2.spin")
        chgcar = wavecar.get_parchg(poscar, 0, 0, phase=False, scale=1)
        assert "total" in chgcar.data
        assert "diff" in chgcar.data
        assert chgcar.data["total"].size == np.prod(wavecar.ng)
        assert np.all(chgcar.data["total"] > 0.0)
        assert not np.all(chgcar.data["diff"] > 0.0)

        chgcar = wavecar.get_parchg(poscar, 0, 0, spin=0, phase=False)
        assert "total" in chgcar.data
        assert "diff" not in chgcar.data
        assert chgcar.data["total"].size == np.prod(wavecar.ng * 2)
        assert np.all(chgcar.data["total"] > 0.0)

        chgcar = wavecar.get_parchg(poscar, 0, 0, spin=0, phase=True)
        assert "total" in chgcar.data
        assert "diff" not in chgcar.data
        assert chgcar.data["total"].size == np.prod(wavecar.ng * 2)
        assert not np.all(chgcar.data["total"] > 0.0)

        wavecar = self.w_ncl
        wavecar.coeffs.append([np.ones((2, 100))])
        chgcar = wavecar.get_parchg(poscar, -1, 0, phase=False, spinor=None)
        assert "total" in chgcar.data
        assert "diff" not in chgcar.data
        assert chgcar.data["total"].size == np.prod(wavecar.ng * 2)
        # this assert was disabled as it started failing during the numpy v2 migration
        # on 2024-08-06. unclear what it was testing in the first place
        # assert not np.all(chgcar.data["total"] > 0.0)

        chgcar = wavecar.get_parchg(poscar, -1, 0, phase=True, spinor=0)
        assert "total" in chgcar.data
        assert "diff" not in chgcar.data
        assert chgcar.data["total"].size == np.prod(wavecar.ng * 2)
        assert not np.all(chgcar.data["total"] > 0.0)

        wavecar.coeffs[-1] = [np.zeros((2, 100))]
        chgcar = wavecar.get_parchg(poscar, -1, 0, phase=False, spinor=1)
        assert "total" in chgcar.data
        assert "diff" not in chgcar.data
        assert chgcar.data["total"].size == np.prod(wavecar.ng * 2)
        assert_allclose(chgcar.data["total"], 0.0)

    def test_write_unks(self):
        unk_std = Unk.from_file(f"{TEST_FILES_DIR}/io/wannier90/UNK.N2.std")
        unk_ncl = Unk.from_file(f"{TEST_FILES_DIR}/io/wannier90/UNK.H2.ncl")

        with pytest.raises(ValueError, match="invalid directory"):
            self.wavecar.write_unks(f"{TEST_FILES_DIR}/io/wannier90/UNK.N2.std")

        # different grids
        self.wavecar.write_unks("./unk_dir")
        assert len(list(Path("./unk_dir").glob("UNK*"))) == 1
        unk = Unk.from_file("./unk_dir/UNK00001.1")
        assert unk != unk_std

        # correct grid
        self.wavecar.ng = np.array([12, 12, 12])
        self.wavecar.write_unks(".")
        unk = Unk.from_file("UNK00001.1")
        assert unk == unk_std

        # ncl test
        self.w_ncl.write_unks(".")
        unk = Unk.from_file("UNK00001.NC")
        assert unk == unk_ncl


class TestEigenval(MatSciTest):
    def test_init(self):
        eig = Eigenval(f"{VASP_OUT_DIR}/EIGENVAL.gz")
        assert eig.ispin == 1
        assert eig.nkpt == len(eig.kpoints)
        assert eig.nkpt == len(eig.kpoints_weights)
        assert eig.nkpt == eig.eigenvalues[Spin.up].shape[0]
        assert eig.nelect == 16
        assert eig.nbands == eig.eigenvalues[Spin.up].shape[1]
        assert np.max(eig.eigenvalues[Spin.up]) > 0
        assert np.min(eig.eigenvalues[Spin.up]) < 0

    def test_ispin2(self):
        eig = Eigenval(f"{VASP_OUT_DIR}/EIGENVAL.ispin2.gz")
        assert eig.ispin == 2
        assert eig.nkpt == eig.eigenvalues[Spin.up].shape[0]
        assert eig.nbands == eig.eigenvalues[Spin.up].shape[1]
        assert eig.nkpt == eig.eigenvalues[Spin.down].shape[0]
        assert eig.nbands == eig.eigenvalues[Spin.down].shape[1]

    def test_eigenvalue_band_properties(self):
        eig = Eigenval(f"{VASP_OUT_DIR}/EIGENVAL.gz")
        props = eig.eigenvalue_band_properties
        assert props[0] == approx(6.4153, abs=1e-4)
        assert props[1] == approx(7.5587, abs=1e-4)
        assert props[2] == approx(1.1434, abs=1e-4)
        assert props[3] is False

    def test_eigenvalue_band_properties_separate_spins(self):
        eig = Eigenval(f"{VASP_OUT_DIR}/EIGENVAL_separate_spins.gz", separate_spins=True)
        props = eig.eigenvalue_band_properties
        eig2 = Eigenval(f"{VASP_OUT_DIR}/EIGENVAL_separate_spins.gz", separate_spins=False)
        props2 = eig2.eigenvalue_band_properties

        assert np.array(props)[:3, :2].flat == approx([2.8772, 1.2810, 3.6741, 1.6225, 0.7969, 0.3415], abs=1e-4)
        assert props2[0] == approx(np.min(props[1]) - np.max(props[2]), abs=1e-4)
        assert props[3][0]
        assert props[3][1]


class TestWaveder(MatSciTest):
    def setup_method(self):
        wder = Waveder.from_binary(f"{VASP_OUT_DIR}/WAVEDER", "float64")
        assert wder.nbands == 36
        assert wder.nkpoints == 56
        band_i = band_j = kp_index = spin_index = cart_dir_index = 0
        cder = wder.get_orbital_derivative_between_states(band_i, band_j, kp_index, spin_index, cart_dir_index)
        assert cder == approx(-1.33639226092e-103, abs=1e-114)

    def test_consistency(self):
        wder_ref = np.loadtxt(f"{VASP_OUT_DIR}/WAVEDERF.Si.gz", skiprows=1)

        def _check(wder):
            with zopen(f"{VASP_OUT_DIR}/WAVEDERF.Si.gz", mode="rt", encoding="utf-8") as file:
                first_line = [int(a) for a in file.readline().split()]
            assert wder.nkpoints == first_line[1]
            assert wder.nbands == first_line[2]
            assert [wder.get_orbital_derivative_between_states(0, idx, 0, 0, 0).real for idx in range(10)] == approx(
                wder_ref[:10, 6], abs=1e-10
            )
            assert wder.cder[0, :10, 0, 0, 0].real == approx(wder_ref[:10, 6], abs=1e-10)
            assert wder.cder[0, :10, 0, 0, 0].imag == approx(wder_ref[:10, 7], abs=1e-10)
            assert wder.cder[0, :10, 0, 0, 1].real == approx(wder_ref[:10, 8], abs=1e-10)
            assert wder.cder[0, :10, 0, 0, 1].imag == approx(wder_ref[:10, 9], abs=1e-10)
            assert wder.cder[0, :10, 0, 0, 2].real == approx(wder_ref[:10, 10], abs=1e-10)
            assert wder.cder[0, :10, 0, 0, 2].imag == approx(wder_ref[:10, 11], abs=1e-10)

        wder = Waveder.from_binary(f"{VASP_OUT_DIR}/WAVEDER.Si")
        _check(wder)
        wderf = Waveder.from_formatted(f"{VASP_OUT_DIR}/WAVEDERF.Si.gz")
        _check(wderf)


class TestWSWQ(MatSciTest):
    def setup_method(self):
        self.wswq = WSWQ.from_file(f"{VASP_OUT_DIR}/WSWQ.gz")

    def test_consistency(self):
        assert self.wswq.nbands == 18
        assert self.wswq.nkpoints == 20
        assert self.wswq.nspin == 2
        assert self.wswq.me_real.shape == (2, 20, 18, 18)
        assert self.wswq.me_imag.shape == (2, 20, 18, 18)
        for itr, (r, i) in enumerate(zip(self.wswq.me_real[0][0][4], self.wswq.me_imag[0][0][4], strict=True)):
            if itr == 4:
                assert np.linalg.norm([r, i]) > 0.999
            else:
                assert np.linalg.norm([r, i]) < 0.001


try:
    import h5py
except ImportError:
    h5py = None


@pytest.mark.skipif(condition=h5py is None, reason="h5py must be installed to use the .Vaspout class.")
class TestVaspout(MatSciTest):
    def setup_method(self):
        self.vaspout = Vaspout(f"{VASP_OUT_DIR}/vaspout.line_mode_band_structure.h5.gz")
        self.vaspout_kpoints_opt = Vaspout(f"{VASP_OUT_DIR}/vaspout.kpoints_opt.h5.gz")

    def test_parse(self):
        from pymatgen.io.vasp.inputs import Incar

        assert self.vaspout.final_energy == approx(-8.953035077096956)
        assert self.vaspout.kpoints.num_kpts == 163
        assert all(not attr for attr in [self.vaspout.is_spin, self.vaspout.is_hubbard])

        input_docs = [(self.vaspout.incar, Incar), (self.vaspout.kpoints, Kpoints), (self.vaspout.potcar, Potcar)]
        assert all(isinstance(*doc) for doc in input_docs)

        assert len(self.vaspout.ionic_steps) == 1

        # double check that these POTCARs have been scrambled
        assert all("FAKE" in psingle.data for psingle in self.vaspout.potcar)

    def test_as_dict(self):
        vout_dict = self.vaspout.as_dict()
        assert isinstance(vout_dict, dict)
        assert all(
            key in vout_dict
            for key in [
                "vasp_version",
                "has_vasp_completed",
                "nsites",
                "unit_cell_formula",
                "reduced_cell_formula",
                "pretty_formula",
                "is_hubbard",
                "hubbards",
                "elements",
                "nelements",
                "run_type",
                "input",
                "output",
            ]
        )

    def test_remove_potcar(self):
        new_vaspout_file = f"{self.tmp_path}/vaspout.h5.gz"
        self.vaspout.remove_potcar_and_write_file(filename=new_vaspout_file)
        cleansed_vout = Vaspout(new_vaspout_file)

        cleansed_vout_d = cleansed_vout.as_dict()
        assert all(cleansed_vout_d[k] == v for k, v in self.vaspout.as_dict().items() if k != "potcar")
        assert cleansed_vout.potcar is None

    def test_kpoints_opt(self):
        assert isinstance(self.vaspout_kpoints_opt.kpoints_opt_props, KpointOptProps)
        assert isinstance(self.vaspout_kpoints_opt.kpoints_opt_props.kpoints, Kpoints)
        assert isinstance(self.vaspout_kpoints_opt.kpoints_opt_props.actual_kpoints, list)
        assert all(
            getattr(self.vaspout_kpoints_opt.kpoints_opt_props, k, None) is not None
            for k in (
                "tdos",
                "idos",
                "pdos",
                "efermi",
                "eigenvalues",
            )
        )

        assert all(
            isinstance(bg_prop, BandgapProps)
            for v in self.vaspout_kpoints_opt.bandgap_props.values()
            for bg_prop in v.values()
        )
