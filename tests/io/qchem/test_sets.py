from __future__ import annotations

import os

import pytest

from pymatgen.io.qchem.sets import (
    ForceSet,
    FreqSet,
    OptSet,
    PESScanSet,
    QChemDictSet,
    QCInput,
    SinglePointSet,
    TransitionStateSet,
)
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

__author__ = "Samuel Blau, Brandon Wood, Shyam Dwaraknath, Evan Spotte-Smith, Ryan Kingsbury"
__copyright__ = "Copyright 2018-2022, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"


test_dir = f"{TEST_FILES_DIR}/molecules"


class TestQChemDictSet(PymatgenTest):
    def test_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_DictSet = QChemDictSet(
            molecule=test_molecule,
            job_type="opt",
            basis_set="6-31G*",
            scf_algorithm="diis",
        )
        ref_dict = {
            "job_type": "opt",
            "gen_scfman": "true",
            "basis": "6-31g*",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "geom_opt_max_cycles": "200",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        assert test_DictSet.rem == ref_dict
        assert test_DictSet.pcm == {}
        assert test_DictSet.solvent == {}
        assert test_DictSet.smx == {}
        assert test_DictSet.molecule == test_molecule

    def test_full_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule

        test_DictSet = QChemDictSet(
            molecule=test_molecule,
            job_type="opt",
            basis_set="6-31g*",
            scf_algorithm="diis",
            dft_rung=2,
            pcm_dielectric=10.0,
            max_scf_cycles=35,
        )
        ref_dict = {
            "job_type": "opt",
            "gen_scfman": "true",
            "basis": "6-31g*",
            "max_scf_cycles": "35",
            "method": "b97-d3",
            "dft_d": "d3_bj",
            "geom_opt_max_cycles": "200",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "pcm",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        assert test_DictSet.rem == ref_dict
        ref_dict = {
            "heavypoints": "194",
            "hpoints": "194",
            "radii": "uff",
            "theory": "cpcm",
            "vdwscale": "1.1",
        }
        assert test_DictSet.pcm == ref_dict
        assert test_DictSet.solvent == {"dielectric": "10.0"}
        assert test_DictSet.molecule == test_molecule

        test_DictSet = QChemDictSet(
            molecule=test_molecule,
            job_type="opt",
            basis_set="6-31g*",
            scf_algorithm="diis",
            dft_rung=2,
            smd_solvent="water",
            max_scf_cycles=35,
        )
        ref_dict = {
            "job_type": "opt",
            "gen_scfman": "true",
            "basis": "6-31g*",
            "max_scf_cycles": "35",
            "method": "b97-d3",
            "dft_d": "d3_bj",
            "geom_opt_max_cycles": "200",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "smd",
            "ideriv": "1",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        assert test_DictSet.rem == ref_dict
        assert test_DictSet.smx == {"solvent": "water"}

    def test_overwrite_input(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        overwrite_inputs = {
            "rem": {
                "method": "b3lyp",
                "basis": "6-31g*",
                "thresh": "10",
                "s2thresh": "12",
                "xc_grid": "000150000302",
            }
        }
        test_OptSet = OptSet(molecule=test_molecule, overwrite_inputs=overwrite_inputs)
        act_rem = {
            "job_type": "opt",
            "gen_scfman": "true",
            "basis": "6-31g*",
            "max_scf_cycles": "100",
            "method": "b3lyp",
            "scf_algorithm": "diis",
            "xc_grid": "000150000302",
            "geom_opt_max_cycles": "200",
            "thresh": "10",
            "s2thresh": "12",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        assert act_rem == test_OptSet.rem

    def test_double_solvation(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        raised_error = False
        dict_set = None
        try:
            dict_set = QChemDictSet(
                molecule=test_molecule,
                job_type="opt",
                basis_set="6-31g*",
                scf_algorithm="diis",
                dft_rung=1,
                pcm_dielectric=10.0,
                smd_solvent="water",
                max_scf_cycles=35,
            )
        except ValueError:
            raised_error = True

        assert raised_error
        assert dict_set is None

    def test_pcm_write(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        dict_set = QChemDictSet(
            molecule=test_molecule,
            job_type="opt",
            basis_set="6-31g*",
            scf_algorithm="diis",
            dft_rung=5,
            pcm_dielectric=10.0,
            max_scf_cycles=35,
        )
        dict_set.write("mol.qin")
        test_dict = QCInput.from_file("mol.qin").as_dict()
        rem = {
            "job_type": "opt",
            "basis": "6-31G*",
            "max_scf_cycles": "35",
            "method": "wb97m(2)",
            "geom_opt_max_cycles": "200",
            "gen_scfman": "true",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "pcm",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        pcm = {
            "heavypoints": "194",
            "hpoints": "194",
            "radii": "uff",
            "theory": "cpcm",
            "vdwscale": "1.1",
        }
        qc_input = QCInput(molecule=test_molecule, rem=rem, pcm=pcm, solvent={"dielectric": 10.0})
        for k, v in qc_input.as_dict().items():
            assert v == test_dict[k]
        os.remove("mol.qin")

    def test_isosvp_write(self):
        """Also tests overwrite_inputs with a RHOISO value."""
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        dict_set = QChemDictSet(
            molecule=test_molecule,
            job_type="opt",
            basis_set="def2-SVPD",
            scf_algorithm="diis",
            dft_rung=4,
            isosvp_dielectric=62,
            max_scf_cycles=35,
            overwrite_inputs={"svp": {"RHOISO": 0.0009}},
        )
        dict_set.write("mol.qin")
        test_dict = QCInput.from_file("mol.qin").as_dict()
        rem = {
            "job_type": "opt",
            "basis": "def2-SVPD",
            "max_scf_cycles": "35",
            "method": "wb97mv",
            "geom_opt_max_cycles": "200",
            "gen_scfman": "false",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "isosvp",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }

        qc_input = QCInput(
            molecule=test_molecule,
            rem=rem,
            pcm_nonels=None,
            svp={"RHOISO": 0.0009, "DIELST": 62, "NPTLEB": 1202, "ITRNGR": 2, "IROTGR": 2},
        )
        for k, v in qc_input.as_dict().items():
            assert v == test_dict[k]
        os.remove("mol.qin")

    def test_smd_write(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        dict_set = QChemDictSet(
            molecule=test_molecule,
            job_type="opt",
            basis_set="6-31g*",
            scf_algorithm="diis",
            dft_rung=4,
            smd_solvent="water",
            max_scf_cycles=35,
        )
        dict_set.write("mol.qin")
        test_dict = QCInput.from_file("mol.qin").as_dict()
        rem = {
            "job_type": "opt",
            "basis": "6-31G*",
            "max_scf_cycles": "35",
            "method": "wb97mv",
            "geom_opt_max_cycles": "200",
            "gen_scfman": "true",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "smd",
            "ideriv": "1",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        qc_input = QCInput(molecule=test_molecule, rem=rem, smx={"solvent": "water"})
        for k, v in qc_input.as_dict().items():
            assert v == test_dict[k]
        os.remove("mol.qin")

    def test_cmirs_write(self):
        """Also tests overwrite_inputs with a RHOISO value."""
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        dict_set = QChemDictSet(
            molecule=test_molecule,
            job_type="opt",
            basis_set="def2-SVPD",
            scf_algorithm="diis",
            dft_rung=4,
            cmirs_solvent="water",
            max_scf_cycles=35,
            overwrite_inputs={"svp": {"RHOISO": 0.0005}},
        )
        dict_set.write("mol.qin")
        test_dict = QCInput.from_file("mol.qin").as_dict()
        rem = {
            "job_type": "opt",
            "basis": "def2-SVPD",
            "max_scf_cycles": "35",
            "method": "wb97mv",
            "geom_opt_max_cycles": "200",
            "gen_scfman": "false",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "isosvp",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        pcm_nonels = {
            "A": -0.006496,
            "B": 0.050833,
            "C": -566.7,
            "D": -30.503,
            "gamma": 3.2,
            "solvrho": 0.05,
            "Delta": 7,
            "GauLag_N": 40,
        }

        qc_input = QCInput(
            molecule=test_molecule,
            rem=rem,
            pcm_nonels=pcm_nonels,
            svp={"RHOISO": 0.0005, "DIELST": 78.39, "NPTLEB": 1202, "ITRNGR": 2, "IROTGR": 2, "IPNRF": 1, "IDEFESR": 1},
        )
        for k, v in qc_input.as_dict().items():
            assert v == test_dict[k]
        os.remove("mol.qin")

    def test_custom_smd_write(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        dict_set = QChemDictSet(
            molecule=test_molecule,
            job_type="opt",
            basis_set="6-31g*",
            scf_algorithm="diis",
            dft_rung=4,
            smd_solvent="custom",
            custom_smd="90.00,1.415,0.00,0.735,20.2,0.00,0.00",
            max_scf_cycles=35,
        )
        dict_set.write("mol.qin")
        test_dict = QCInput.from_file("mol.qin").as_dict()
        rem = {
            "job_type": "opt",
            "basis": "6-31G*",
            "max_scf_cycles": "35",
            "method": "wb97mv",
            "geom_opt_max_cycles": "200",
            "gen_scfman": "true",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "smd",
            "ideriv": "1",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        qc_input = QCInput(molecule=test_molecule, rem=rem, smx={"solvent": "other"})
        for k, v in qc_input.as_dict().items():
            assert v == test_dict[k]
        os.remove("mol.qin")
        with open("solvent_data") as sd:
            lines = sd.readlines()
            assert lines[0] == "90.00,1.415,0.00,0.735,20.2,0.00,0.00"
        os.remove("solvent_data")

    def test_solvation_warnings(self):
        """Tests warnings / errors resulting from nonsensical overwrite_inputs."""
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        with pytest.raises(RuntimeError, match="CMIRS is only parameterized"):
            QChemDictSet(
                molecule=test_molecule,
                job_type="opt",
                basis_set="def2-SVPD",
                scf_algorithm="diis",
                dft_rung=5,
                cmirs_solvent="water",
                max_scf_cycles=35,
                overwrite_inputs={"svp": {"RHOISO": 0.0007}},
            )
        with pytest.warns(UserWarning, match="Setting IDEFESR=0"):
            QChemDictSet(
                molecule=test_molecule,
                job_type="opt",
                basis_set="def2-SVPD",
                scf_algorithm="diis",
                dft_rung=5,
                cmirs_solvent="water",
                max_scf_cycles=35,
                overwrite_inputs={"svp": {"IDEFESR": 0}},
            )
        with pytest.warns(UserWarning, match="Setting IDEFESR=1"):
            QChemDictSet(
                molecule=test_molecule,
                job_type="opt",
                basis_set="def2-SVPD",
                scf_algorithm="diis",
                dft_rung=5,
                isosvp_dielectric=78,
                max_scf_cycles=35,
                overwrite_inputs={"svp": {"IDEFESR": 1}},
            )
        with pytest.warns(UserWarning, match="Setting DIELST"):
            QChemDictSet(
                molecule=test_molecule,
                job_type="opt",
                basis_set="def2-SVPD",
                scf_algorithm="diis",
                dft_rung=5,
                pcm_dielectric=78,
                max_scf_cycles=35,
                overwrite_inputs={"svp": {"DIELST": 67}},
            )
        with pytest.warns(UserWarning, match="The solvent section will be ignored"):
            QChemDictSet(
                molecule=test_molecule,
                job_type="opt",
                basis_set="def2-SVPD",
                scf_algorithm="diis",
                dft_rung=5,
                isosvp_dielectric=78,
                max_scf_cycles=35,
                overwrite_inputs={"solvent": {"dielectric": 67}},
            )


class TestSinglePointSet(PymatgenTest):
    def test_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_SPSet = SinglePointSet(molecule=test_molecule)
        ref_dict = {
            "job_type": "sp",
            "gen_scfman": "true",
            "basis": "def2-tzvpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        assert test_SPSet.rem == ref_dict
        assert test_SPSet.pcm == {}
        assert test_SPSet.solvent == {}
        assert test_SPSet.molecule == test_molecule

    def test_scf_extra_print(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_SPSet = SinglePointSet(molecule=test_molecule, extra_scf_print=True)
        ref_dict = {
            "job_type": "sp",
            "gen_scfman": "true",
            "basis": "def2-tzvpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
            "scf_convergence": "8",
            "scf_final_print": "3",
        }
        assert test_SPSet.rem == ref_dict
        assert test_SPSet.pcm == {}
        assert test_SPSet.solvent == {}
        assert test_SPSet.molecule == test_molecule

    def test_pcm_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_SPSet = SinglePointSet(molecule=test_molecule, pcm_dielectric=10.0)
        ref_dict = {
            "job_type": "sp",
            "gen_scfman": "true",
            "basis": "def2-tzvpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "pcm",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        assert test_SPSet.rem == ref_dict
        ref_dict = {
            "heavypoints": "194",
            "hpoints": "194",
            "radii": "uff",
            "theory": "cpcm",
            "vdwscale": "1.1",
        }
        assert test_SPSet.pcm == ref_dict
        assert test_SPSet.solvent == {"dielectric": "10.0"}
        assert test_SPSet.molecule == test_molecule

    def test_isosvp_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_SPSet = SinglePointSet(molecule=test_molecule, isosvp_dielectric=10.0)
        ref_dict = {
            "job_type": "sp",
            "gen_scfman": "false",
            "basis": "def2-tzvpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "isosvp",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        assert test_SPSet.rem == ref_dict
        ref_dict = {"dielst": "10.0", "rhoiso": "0.001", "nptleb": "1202", "itrngr": "2", "irotgr": "2"}
        assert test_SPSet.svp == ref_dict
        assert test_SPSet.molecule == test_molecule

    def test_smd_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_SPSet = SinglePointSet(molecule=test_molecule, smd_solvent="water")
        ref_dict = {
            "job_type": "sp",
            "gen_scfman": "true",
            "basis": "def2-tzvpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "smd",
            "ideriv": "1",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        assert test_SPSet.rem == ref_dict
        assert test_SPSet.smx == {"solvent": "water"}
        assert test_SPSet.molecule == test_molecule

    def test_cmirs_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_SPSet = SinglePointSet(
            molecule=test_molecule, cmirs_solvent="benzene", overwrite_inputs={"svp": {"RHOISO": 0.0005}}
        )
        ref_dict = {
            "job_type": "sp",
            "gen_scfman": "false",
            "basis": "def2-tzvpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "isosvp",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        assert test_SPSet.rem == ref_dict
        ref_dict = {
            "dielst": "2.28",
            "rhoiso": "0.0005",
            "nptleb": "1202",
            "itrngr": "2",
            "irotgr": "2",
            "ipnrf": "1",
            "idefesr": "1",
        }
        assert test_SPSet.svp == ref_dict
        ref_dict = {
            "a": "-0.00572",
            "b": "0.01116",
            "c": None,
            "d": None,
            "gamma": None,
            "solvrho": "0.0421",
            "gaulag_n": "40",
            "delta": "7",
        }
        assert test_SPSet.pcm_nonels == ref_dict
        assert test_SPSet.molecule == test_molecule

    def test_plots_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_SPSet = SinglePointSet(molecule=test_molecule, smd_solvent="water", plot_cubes=True)
        ref_dict = {
            "job_type": "sp",
            "gen_scfman": "true",
            "basis": "def2-tzvpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "smd",
            "ideriv": "1",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
            "plots": "true",
            "make_cube_files": "true",
        }
        assert test_SPSet.rem == ref_dict
        assert test_SPSet.plots == {"grid_spacing": "0.05", "total_density": "0"}
        assert test_SPSet.smx == {"solvent": "water"}
        assert test_SPSet.molecule == test_molecule


class TestOptSet(PymatgenTest):
    def test_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_OptSet = OptSet(molecule=test_molecule)
        ref_dict = {
            "job_type": "opt",
            "gen_scfman": "true",
            "basis": "def2-svpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "geom_opt_max_cycles": "200",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        assert test_OptSet.rem == ref_dict
        assert test_OptSet.pcm == {}
        assert test_OptSet.solvent == {}
        assert test_OptSet.smx == {}
        assert test_OptSet.molecule == test_molecule

    def test_pcm_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_OptSet = OptSet(molecule=test_molecule, pcm_dielectric=10.0)
        ref_dict = {
            "job_type": "opt",
            "gen_scfman": "true",
            "basis": "def2-svpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "geom_opt_max_cycles": "200",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "pcm",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        assert test_OptSet.rem == ref_dict
        ref_dict = {
            "heavypoints": "194",
            "hpoints": "194",
            "radii": "uff",
            "theory": "cpcm",
            "vdwscale": "1.1",
        }
        assert test_OptSet.pcm == ref_dict
        assert test_OptSet.solvent == {"dielectric": "10.0"}
        assert test_OptSet.molecule == test_molecule

    def test_smd_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_OptSet = OptSet(molecule=test_molecule, smd_solvent="water")
        ref_dict = {
            "job_type": "opt",
            "gen_scfman": "true",
            "basis": "def2-svpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "geom_opt_max_cycles": "200",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "smd",
            "ideriv": "1",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        assert test_OptSet.rem == ref_dict
        assert test_OptSet.smx == {"solvent": "water"}
        assert test_OptSet.molecule == test_molecule

    def test_overwrite_opt_input(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        overwrite_inputs = {"opt": {"FIXED": ["1 XYZ", "2 XY"]}}
        test_OptSet = OptSet(molecule=test_molecule, overwrite_inputs=overwrite_inputs)
        act_opt = {"fixed": ["1 XYZ", "2 XY"]}
        assert act_opt == test_OptSet.opt

    def test_nbo_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_OptSet = OptSet(molecule=test_molecule, nbo_params={})
        ref_dict = {
            "job_type": "opt",
            "gen_scfman": "true",
            "geom_opt_max_cycles": "200",
            "basis": "def2-svpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
            "nbo": "true",
        }
        assert test_OptSet.rem == ref_dict
        assert test_OptSet.nbo == {}
        assert test_OptSet.molecule == test_molecule

    def test_v5_vs_v6(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        v5_OptSet = OptSet(molecule=test_molecule, qchem_version=5, basis_set="def2-tzvpd", geom_opt={})
        ref_dict = {
            "job_type": "opt",
            "gen_scfman": "true",
            "geom_opt_max_cycles": "200",
            "basis": "def2-tzvpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
            "geom_opt2": "3",
        }
        assert v5_OptSet.rem == ref_dict
        assert v5_OptSet.geom_opt == {"maxiter": "200"}
        assert v5_OptSet.molecule == test_molecule

        v6_OptSet = OptSet(molecule=test_molecule, qchem_version=6, basis_set="def2-tzvpd", geom_opt={})
        ref_dict = {
            "job_type": "opt",
            "gen_scfman": "true",
            "geom_opt_max_cycles": "200",
            "basis": "def2-tzvpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        assert v6_OptSet.rem == ref_dict
        ref_dict = {
            "maxiter": "200",
            "coordinates": "redundant",
            "max_displacement": "0.1",
            "optimization_restart": "false",
        }
        assert v6_OptSet.geom_opt == ref_dict
        assert v6_OptSet.molecule == test_molecule

        v6_OptSet_modified = OptSet(
            molecule=test_molecule,
            qchem_version=6,
            basis_set="def2-tzvpd",
            geom_opt={"coordinates": "delocalized", "initial_hessian": "read"},
        )
        ref_dict = {
            "maxiter": "200",
            "coordinates": "delocalized",
            "max_displacement": "0.1",
            "initial_hessian": "read",
            "optimization_restart": "false",
        }
        assert v6_OptSet_modified.geom_opt == ref_dict


class TestTransitionStateSet(PymatgenTest):
    def test_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_TSSet = TransitionStateSet(molecule=test_molecule)
        ref_dict = {
            "job_type": "ts",
            "gen_scfman": "true",
            "basis": "def2-svpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "geom_opt_max_cycles": "200",
            "resp_charges": "true",
            "sym_ignore": "true",
            "symmetry": "false",
        }
        assert test_TSSet.rem == ref_dict
        assert test_TSSet.pcm == {}
        assert test_TSSet.solvent == {}
        assert test_TSSet.smx == {}
        assert test_TSSet.molecule == test_molecule

    def test_pcm_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_TSSet = TransitionStateSet(molecule=test_molecule, pcm_dielectric=10.0)
        ref_dict = {
            "job_type": "ts",
            "gen_scfman": "true",
            "basis": "def2-svpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "geom_opt_max_cycles": "200",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "pcm",
            "resp_charges": "true",
            "sym_ignore": "true",
            "symmetry": "false",
        }
        assert test_TSSet.rem == ref_dict
        ref_dict = {"heavypoints": "194", "hpoints": "194", "radii": "uff", "theory": "cpcm", "vdwscale": "1.1"}
        assert test_TSSet.pcm == ref_dict
        assert test_TSSet.solvent == {"dielectric": "10.0"}
        assert test_TSSet.molecule == test_molecule

    def test_smd_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_TSSet = TransitionStateSet(molecule=test_molecule, smd_solvent="water")
        ref_dict = {
            "job_type": "ts",
            "gen_scfman": "true",
            "basis": "def2-svpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "geom_opt_max_cycles": "200",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "smd",
            "ideriv": "1",
            "resp_charges": "true",
            "sym_ignore": "true",
            "symmetry": "false",
        }
        assert test_TSSet.rem == ref_dict
        assert test_TSSet.smx == {"solvent": "water"}
        assert test_TSSet.molecule == test_molecule


class TestForceSet(PymatgenTest):
    def test_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_forceset = ForceSet(molecule=test_molecule)
        ref_dict = {
            "job_type": "force",
            "gen_scfman": "true",
            "basis": "def2-tzvpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        assert test_forceset.rem == ref_dict
        assert test_forceset.pcm == {}
        assert test_forceset.solvent == {}
        assert test_forceset.molecule == test_molecule

    def test_pcm_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_forceset = ForceSet(molecule=test_molecule, pcm_dielectric=10.0)
        ref_dict = {
            "job_type": "force",
            "gen_scfman": "true",
            "basis": "def2-tzvpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "pcm",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        assert test_forceset.rem == ref_dict
        ref_dict = {"heavypoints": "194", "hpoints": "194", "radii": "uff", "theory": "cpcm", "vdwscale": "1.1"}
        assert test_forceset.pcm == ref_dict
        assert test_forceset.solvent == {"dielectric": "10.0"}
        assert test_forceset.molecule == test_molecule

    def test_smd_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_forceset = ForceSet(molecule=test_molecule, smd_solvent="water")
        ref_dict = {
            "job_type": "force",
            "gen_scfman": "true",
            "basis": "def2-tzvpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "smd",
            "ideriv": "1",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        assert test_forceset.rem == ref_dict
        assert test_forceset.smx, {"solvent": "water"}
        assert test_forceset.molecule == test_molecule


class TestPESScanSet(PymatgenTest):
    def test_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pes_scan.qin").molecule

        test_pes_scan = PESScanSet(molecule=test_molecule, scan_variables={"stre": ["3 6 1.5 1.9 0.01"]})
        ref_dict = {
            "job_type": "pes_scan",
            "gen_scfman": "true",
            "basis": "def2-svpd",
            "max_scf_cycles": "100",
            "geom_opt_max_cycles": "200",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "resp_charges": "true",
            "sym_ignore": "true",
            "symmetry": "false",
        }
        assert test_pes_scan.rem == ref_dict
        assert test_pes_scan.pcm == {}
        assert test_pes_scan.solvent == {}
        assert test_pes_scan.smx == {}
        assert test_pes_scan.scan == {"stre": ["3 6 1.5 1.9 0.01"]}
        assert test_pes_scan.molecule == test_molecule

    def test_pcm_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pes_scan.qin").molecule
        test_pes_scan = PESScanSet(
            molecule=test_molecule, pcm_dielectric=10.0, scan_variables={"stre": ["3 6 1.5 1.9 0.01"]}
        )
        ref_dict = {
            "job_type": "pes_scan",
            "gen_scfman": "true",
            "basis": "def2-svpd",
            "max_scf_cycles": "100",
            "geom_opt_max_cycles": "200",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "pcm",
            "resp_charges": "true",
            "sym_ignore": "true",
            "symmetry": "false",
        }
        assert test_pes_scan.rem == ref_dict
        ref_dict = {"heavypoints": "194", "hpoints": "194", "radii": "uff", "theory": "cpcm", "vdwscale": "1.1"}
        assert test_pes_scan.pcm == ref_dict
        assert test_pes_scan.solvent == {"dielectric": "10.0"}
        assert test_pes_scan.scan == {"stre": ["3 6 1.5 1.9 0.01"]}
        assert test_pes_scan.molecule == test_molecule

    def test_smd_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pes_scan.qin").molecule
        test_pes_scan = PESScanSet(
            molecule=test_molecule, smd_solvent="water", scan_variables={"stre": ["3 6 1.5 1.9 0.01"]}
        )
        ref_dict = {
            "job_type": "pes_scan",
            "gen_scfman": "true",
            "basis": "def2-svpd",
            "max_scf_cycles": "100",
            "geom_opt_max_cycles": "200",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "smd",
            "ideriv": "1",
            "resp_charges": "true",
            "sym_ignore": "true",
            "symmetry": "false",
        }
        assert test_pes_scan.rem == ref_dict
        assert test_pes_scan.smx == {"solvent": "water"}
        assert test_pes_scan.scan == {"stre": ["3 6 1.5 1.9 0.01"]}
        assert test_pes_scan.molecule == test_molecule


class TestFreqSet(PymatgenTest):
    def test_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_FreqSet = FreqSet(molecule=test_molecule)
        ref_dict = {
            "job_type": "freq",
            "gen_scfman": "true",
            "basis": "def2-svpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        assert test_FreqSet.rem == ref_dict
        assert test_FreqSet.pcm == {}
        assert test_FreqSet.solvent == {}
        assert test_FreqSet.molecule == test_molecule

    def test_pcm_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_FreqSet = FreqSet(molecule=test_molecule, pcm_dielectric=10.0)
        ref_dict = {
            "job_type": "freq",
            "gen_scfman": "true",
            "basis": "def2-svpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "pcm",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        assert test_FreqSet.rem == ref_dict
        ref_dict = {
            "heavypoints": "194",
            "hpoints": "194",
            "radii": "uff",
            "theory": "cpcm",
            "vdwscale": "1.1",
        }
        assert test_FreqSet.pcm == ref_dict
        assert test_FreqSet.solvent == {"dielectric": "10.0"}
        assert test_FreqSet.molecule == test_molecule

    def test_smd_init(self):
        test_molecule = QCInput.from_file(f"{test_dir}/new_qchem_files/pcm.qin").molecule
        test_FreqSet = FreqSet(molecule=test_molecule, smd_solvent="water")
        ref_dict = {
            "job_type": "freq",
            "gen_scfman": "true",
            "basis": "def2-svpd",
            "max_scf_cycles": "100",
            "method": "wb97mv",
            "scf_algorithm": "diis",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "solvent_method": "smd",
            "ideriv": "1",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        assert test_FreqSet.rem == ref_dict
        assert test_FreqSet.smx == {"solvent": "water"}
        assert test_FreqSet.molecule == test_molecule
