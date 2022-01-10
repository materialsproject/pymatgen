# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import os
import unittest

from pymatgen.io.qchem.sets import *
from pymatgen.util.testing import PymatgenTest

__author__ = "Samuel Blau, Brandon Wood, Shyam Dwaraknath, Evan Spotte-Smith"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"


test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules")


class QChemDictSetTest(PymatgenTest):
    def test_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_DictSet = QChemDictSet(
            molecule=test_molecule,
            job_type="opt",
            basis_set="6-31G*",
            scf_algorithm="diis",
        )
        self.assertEqual(
            test_DictSet.rem,
            {
                "job_type": "opt",
                "gen_scfman": "true",
                "basis": "6-31g*",
                "max_scf_cycles": "200",
                "method": "wb97xv",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "geom_opt_max_cycles": "200",
                "symmetry": "false",
                "sym_ignore": "true",
                "resp_charges": "true",
            },
        )
        self.assertEqual(test_DictSet.pcm, {})
        self.assertEqual(test_DictSet.solvent, {})
        self.assertEqual(test_DictSet.smx, {})
        self.assertEqual(test_DictSet.molecule, test_molecule)

    def test_full_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule

        test_DictSet = QChemDictSet(
            molecule=test_molecule,
            job_type="opt",
            basis_set="6-31g*",
            scf_algorithm="diis",
            dft_rung=1,
            pcm_dielectric=10.0,
            max_scf_cycles=35,
        )
        self.assertEqual(
            test_DictSet.rem,
            {
                "job_type": "opt",
                "gen_scfman": "true",
                "basis": "6-31g*",
                "max_scf_cycles": "35",
                "method": "b3lyp",
                "geom_opt_max_cycles": "200",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "solvent_method": "pcm",
                "symmetry": "false",
                "sym_ignore": "true",
                "resp_charges": "true",
            },
        )
        self.assertEqual(
            test_DictSet.pcm,
            {
                "heavypoints": "194",
                "hpoints": "194",
                "radii": "uff",
                "theory": "cpcm",
                "vdwscale": "1.1",
            },
        )
        self.assertEqual(test_DictSet.solvent, {"dielectric": 10.0})
        self.assertEqual(test_DictSet.molecule, test_molecule)

        test_DictSet = QChemDictSet(
            molecule=test_molecule,
            job_type="opt",
            basis_set="6-31g*",
            scf_algorithm="diis",
            dft_rung=1,
            smd_solvent="water",
            max_scf_cycles=35,
        )
        self.assertEqual(
            test_DictSet.rem,
            {
                "job_type": "opt",
                "gen_scfman": "true",
                "basis": "6-31g*",
                "max_scf_cycles": "35",
                "method": "b3lyp",
                "geom_opt_max_cycles": "200",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "solvent_method": "smd",
                "ideriv": "1",
                "symmetry": "false",
                "sym_ignore": "true",
                "resp_charges": "true",
            },
        )
        self.assertEqual(test_DictSet.smx, {"solvent": "water"})

    def test_overwrite_input(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        overwrite_inputs = {
            "rem": {
                "method": "b3lyp",
                "basis": "6-31g*",
                "thresh": 10,
                "xc_grid": "000150000302",
            }
        }
        test_OptSet = OptSet(molecule=test_molecule, overwrite_inputs=overwrite_inputs)
        act_rem = {
            "job_type": "opt",
            "gen_scfman": "true",
            "basis": "6-31g*",
            "max_scf_cycles": "200",
            "method": "b3lyp",
            "scf_algorithm": "diis",
            "xc_grid": "000150000302",
            "geom_opt_max_cycles": "200",
            "thresh": 10,
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        self.assertDictEqual(act_rem, test_OptSet.rem)

    def test_double_solvation(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
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

        self.assertTrue(raised_error)
        self.assertEqual(dict_set, None)

    def test_pcm_write(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
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
            "method": "wb97mv",
            "geom_opt_max_cycles": "200",
            "gen_scfman": "true",
            "scf_algorithm": "diis",
            "xc_grid": "3",
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
        qc_input = QCInput(molecule=test_molecule, rem=rem, pcm=pcm, solvent={"dielectric": "10.0"})
        for k, v in qc_input.as_dict().items():
            self.assertEqual(v, test_dict[k])
        os.remove("mol.qin")

    def test_smd_write(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        dict_set = QChemDictSet(
            molecule=test_molecule,
            job_type="opt",
            basis_set="6-31g*",
            scf_algorithm="diis",
            dft_rung=5,
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
            "solvent_method": "smd",
            "ideriv": "1",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        qc_input = QCInput(molecule=test_molecule, rem=rem, smx={"solvent": "water"})
        for k, v in qc_input.as_dict().items():
            self.assertEqual(v, test_dict[k])
        os.remove("mol.qin")

    def test_custom_smd_write(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        dict_set = QChemDictSet(
            molecule=test_molecule,
            job_type="opt",
            basis_set="6-31g*",
            scf_algorithm="diis",
            dft_rung=5,
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
            "solvent_method": "smd",
            "ideriv": "1",
            "symmetry": "false",
            "sym_ignore": "true",
            "resp_charges": "true",
        }
        qc_input = QCInput(molecule=test_molecule, rem=rem, smx={"solvent": "other"})
        for k, v in qc_input.as_dict().items():
            self.assertEqual(v, test_dict[k])
        os.remove("mol.qin")
        with open("solvent_data") as sd:
            lines = sd.readlines()
            self.assertEqual(lines[0], "90.00,1.415,0.00,0.735,20.2,0.00,0.00")
        os.remove("solvent_data")


class SinglePointSetTest(PymatgenTest):
    def test_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_SPSet = SinglePointSet(molecule=test_molecule)
        self.assertEqual(
            test_SPSet.rem,
            {
                "job_type": "sp",
                "gen_scfman": "true",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "method": "wb97xd",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "symmetry": "false",
                "sym_ignore": "true",
                "resp_charges": "true",
            },
        )
        self.assertEqual(test_SPSet.pcm, {})
        self.assertEqual(test_SPSet.solvent, {})
        self.assertEqual(test_SPSet.molecule, test_molecule)

    def test_pcm_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_SPSet = SinglePointSet(molecule=test_molecule, pcm_dielectric=10.0)
        self.assertEqual(
            test_SPSet.rem,
            {
                "job_type": "sp",
                "gen_scfman": "true",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "method": "wb97xd",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "solvent_method": "pcm",
                "symmetry": "false",
                "sym_ignore": "true",
                "resp_charges": "true",
            },
        )
        self.assertEqual(
            test_SPSet.pcm,
            {
                "heavypoints": "194",
                "hpoints": "194",
                "radii": "uff",
                "theory": "cpcm",
                "vdwscale": "1.1",
            },
        )
        self.assertEqual(test_SPSet.solvent, {"dielectric": 10.0})
        self.assertEqual(test_SPSet.molecule, test_molecule)

    def test_smd_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_SPSet = SinglePointSet(molecule=test_molecule, smd_solvent="water")
        self.assertEqual(
            test_SPSet.rem,
            {
                "job_type": "sp",
                "gen_scfman": "true",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "method": "wb97xd",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "solvent_method": "smd",
                "ideriv": "1",
                "symmetry": "false",
                "sym_ignore": "true",
                "resp_charges": "true",
            },
        )
        self.assertEqual(test_SPSet.smx, {"solvent": "water"})
        self.assertEqual(test_SPSet.molecule, test_molecule)

    def test_plots_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_SPSet = SinglePointSet(molecule=test_molecule, smd_solvent="water", plot_cubes=True)
        self.assertEqual(
            test_SPSet.rem,
            {
                "job_type": "sp",
                "gen_scfman": "true",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "method": "wb97xd",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "solvent_method": "smd",
                "ideriv": "1",
                "symmetry": "false",
                "sym_ignore": "true",
                "resp_charges": "true",
                "plots": "true",
                "make_cube_files": "true",
            },
        )
        self.assertEqual(test_SPSet.plots, {"grid_spacing": "0.05", "total_density": "0"})
        self.assertEqual(test_SPSet.smx, {"solvent": "water"})
        self.assertEqual(test_SPSet.molecule, test_molecule)


class OptSetTest(PymatgenTest):
    def test_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_OptSet = OptSet(molecule=test_molecule)
        self.assertEqual(
            test_OptSet.rem,
            {
                "job_type": "opt",
                "gen_scfman": "true",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "method": "wb97xd",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "geom_opt_max_cycles": "200",
                "symmetry": "false",
                "sym_ignore": "true",
                "resp_charges": "true",
            },
        )
        self.assertEqual(test_OptSet.pcm, {})
        self.assertEqual(test_OptSet.solvent, {})
        self.assertEqual(test_OptSet.smx, {})
        self.assertEqual(test_OptSet.molecule, test_molecule)

    def test_pcm_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_OptSet = OptSet(molecule=test_molecule, pcm_dielectric=10.0)
        self.assertEqual(
            test_OptSet.rem,
            {
                "job_type": "opt",
                "gen_scfman": "true",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "method": "wb97xd",
                "geom_opt_max_cycles": "200",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "solvent_method": "pcm",
                "symmetry": "false",
                "sym_ignore": "true",
                "resp_charges": "true",
            },
        )
        self.assertEqual(
            test_OptSet.pcm,
            {
                "heavypoints": "194",
                "hpoints": "194",
                "radii": "uff",
                "theory": "cpcm",
                "vdwscale": "1.1",
            },
        )
        self.assertEqual(test_OptSet.solvent, {"dielectric": 10.0})
        self.assertEqual(test_OptSet.molecule, test_molecule)

    def test_smd_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_OptSet = OptSet(molecule=test_molecule, smd_solvent="water")
        self.assertEqual(
            test_OptSet.rem,
            {
                "job_type": "opt",
                "gen_scfman": "true",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "method": "wb97xd",
                "geom_opt_max_cycles": "200",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "solvent_method": "smd",
                "ideriv": "1",
                "symmetry": "false",
                "sym_ignore": "true",
                "resp_charges": "true",
            },
        )
        self.assertEqual(test_OptSet.smx, {"solvent": "water"})
        self.assertEqual(test_OptSet.molecule, test_molecule)

    def test_overwrite_opt_input(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        overwrite_inputs = {"opt": {"FIXED": ["1 XYZ", "2 XY"]}}
        test_OptSet = OptSet(molecule=test_molecule, overwrite_inputs=overwrite_inputs)
        act_opt = {"fixed": ["1 XYZ", "2 XY"]}
        self.assertDictEqual(act_opt, test_OptSet.opt)

    def test_nbo_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_OptSet = OptSet(molecule=test_molecule, nbo_params={})
        self.assertEqual(
            test_OptSet.rem,
            {
                "job_type": "opt",
                "gen_scfman": "true",
                "geom_opt_max_cycles": "200",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "method": "wb97xd",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "symmetry": "false",
                "sym_ignore": "true",
                "resp_charges": "true",
                "nbo": "true",
            },
        )
        self.assertEqual(test_OptSet.nbo, {})
        self.assertEqual(test_OptSet.molecule, test_molecule)


class TransitionStateSetTest(PymatgenTest):
    def test_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_TSSet = TransitionStateSet(molecule=test_molecule)
        self.assertEqual(
            test_TSSet.rem,
            {
                "job_type": "ts",
                "gen_scfman": "true",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "method": "wb97xd",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "geom_opt_max_cycles": "200",
                "resp_charges": "true",
                "sym_ignore": "true",
                "symmetry": "false",
            },
        )
        self.assertEqual(test_TSSet.pcm, {})
        self.assertEqual(test_TSSet.solvent, {})
        self.assertEqual(test_TSSet.smx, {})
        self.assertEqual(test_TSSet.molecule, test_molecule)

    def test_pcm_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_TSSet = TransitionStateSet(molecule=test_molecule, pcm_dielectric=10.0)
        self.assertEqual(
            test_TSSet.rem,
            {
                "job_type": "ts",
                "gen_scfman": "true",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "method": "wb97xd",
                "geom_opt_max_cycles": "200",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "solvent_method": "pcm",
                "resp_charges": "true",
                "sym_ignore": "true",
                "symmetry": "false",
            },
        )
        self.assertEqual(
            test_TSSet.pcm,
            {"heavypoints": "194", "hpoints": "194", "radii": "uff", "theory": "cpcm", "vdwscale": "1.1"},
        )
        self.assertEqual(test_TSSet.solvent, {"dielectric": 10.0})
        self.assertEqual(test_TSSet.molecule, test_molecule)

    def test_smd_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_TSSet = TransitionStateSet(molecule=test_molecule, smd_solvent="water")
        self.assertEqual(
            test_TSSet.rem,
            {
                "job_type": "ts",
                "gen_scfman": "true",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "method": "wb97xd",
                "geom_opt_max_cycles": "200",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "solvent_method": "smd",
                "ideriv": "1",
                "resp_charges": "true",
                "sym_ignore": "true",
                "symmetry": "false",
            },
        )
        self.assertEqual(test_TSSet.smx, {"solvent": "water"})
        self.assertEqual(test_TSSet.molecule, test_molecule)


class ForceSetTest(PymatgenTest):
    def test_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_forceset = ForceSet(molecule=test_molecule)
        self.assertEqual(
            test_forceset.rem,
            {
                "job_type": "force",
                "gen_scfman": "true",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "method": "wb97xd",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "symmetry": "false",
                "sym_ignore": "true",
                "resp_charges": "true",
            },
        )
        self.assertEqual(test_forceset.pcm, {})
        self.assertEqual(test_forceset.solvent, {})
        self.assertEqual(test_forceset.molecule, test_molecule)

    def test_pcm_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_forceset = ForceSet(molecule=test_molecule, pcm_dielectric=10.0)
        self.assertEqual(
            test_forceset.rem,
            {
                "job_type": "force",
                "gen_scfman": "true",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "method": "wb97xd",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "solvent_method": "pcm",
                "symmetry": "false",
                "sym_ignore": "true",
                "resp_charges": "true",
            },
        )
        self.assertEqual(
            test_forceset.pcm,
            {"heavypoints": "194", "hpoints": "194", "radii": "uff", "theory": "cpcm", "vdwscale": "1.1"},
        )
        self.assertEqual(test_forceset.solvent, {"dielectric": 10.0})
        self.assertEqual(test_forceset.molecule, test_molecule)

    def test_smd_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_forceset = ForceSet(molecule=test_molecule, smd_solvent="water")
        self.assertEqual(
            test_forceset.rem,
            {
                "job_type": "force",
                "gen_scfman": "true",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "method": "wb97xd",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "solvent_method": "smd",
                "ideriv": "1",
                "symmetry": "false",
                "sym_ignore": "true",
                "resp_charges": "true",
            },
        )
        self.assertEqual(test_forceset.smx, {"solvent": "water"})
        self.assertEqual(test_forceset.molecule, test_molecule)


class PESScanSetTest(PymatgenTest):
    def test_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pes_scan.qin")).molecule

        test_pes_scan = PESScanSet(molecule=test_molecule, scan_variables={"stre": ["3 6 1.5 1.9 0.01"]})
        self.assertEqual(
            test_pes_scan.rem,
            {
                "job_type": "pes_scan",
                "gen_scfman": "true",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "geom_opt_max_cycles": "200",
                "method": "wb97xd",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "resp_charges": "true",
                "sym_ignore": "true",
                "symmetry": "false",
            },
        )
        self.assertEqual(test_pes_scan.pcm, dict())
        self.assertEqual(test_pes_scan.solvent, dict())
        self.assertEqual(test_pes_scan.smx, dict())
        self.assertEqual(test_pes_scan.scan, {"stre": ["3 6 1.5 1.9 0.01"]})
        self.assertEqual(test_pes_scan.molecule, test_molecule)

    def test_pcm_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pes_scan.qin")).molecule
        test_pes_scan = PESScanSet(
            molecule=test_molecule, pcm_dielectric=10.0, scan_variables={"stre": ["3 6 1.5 1.9 0.01"]}
        )
        self.assertEqual(
            test_pes_scan.rem,
            {
                "job_type": "pes_scan",
                "gen_scfman": "true",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "geom_opt_max_cycles": "200",
                "method": "wb97xd",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "solvent_method": "pcm",
                "resp_charges": "true",
                "sym_ignore": "true",
                "symmetry": "false",
            },
        )
        self.assertEqual(
            test_pes_scan.pcm,
            {"heavypoints": "194", "hpoints": "194", "radii": "uff", "theory": "cpcm", "vdwscale": "1.1"},
        )
        self.assertEqual(test_pes_scan.solvent, {"dielectric": 10.0})
        self.assertEqual(test_pes_scan.scan, {"stre": ["3 6 1.5 1.9 0.01"]})
        self.assertEqual(test_pes_scan.molecule, test_molecule)

    def test_smd_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pes_scan.qin")).molecule
        test_pes_scan = PESScanSet(
            molecule=test_molecule, smd_solvent="water", scan_variables={"stre": ["3 6 1.5 1.9 0.01"]}
        )
        self.assertEqual(
            test_pes_scan.rem,
            {
                "job_type": "pes_scan",
                "gen_scfman": "true",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "geom_opt_max_cycles": "200",
                "method": "wb97xd",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "solvent_method": "smd",
                "ideriv": "1",
                "resp_charges": "true",
                "sym_ignore": "true",
                "symmetry": "false",
            },
        )
        self.assertEqual(test_pes_scan.smx, {"solvent": "water"})
        self.assertEqual(test_pes_scan.scan, {"stre": ["3 6 1.5 1.9 0.01"]})
        self.assertEqual(test_pes_scan.molecule, test_molecule)


class FreqSetTest(PymatgenTest):
    def test_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_FreqSet = FreqSet(molecule=test_molecule)
        self.assertEqual(
            test_FreqSet.rem,
            {
                "job_type": "freq",
                "gen_scfman": "true",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "method": "wb97xd",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "symmetry": "false",
                "sym_ignore": "true",
                "resp_charges": "true",
            },
        )
        self.assertEqual(test_FreqSet.pcm, {})
        self.assertEqual(test_FreqSet.solvent, {})
        self.assertEqual(test_FreqSet.molecule, test_molecule)

    def test_pcm_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_FreqSet = FreqSet(molecule=test_molecule, pcm_dielectric=10.0)
        self.assertEqual(
            test_FreqSet.rem,
            {
                "job_type": "freq",
                "gen_scfman": "true",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "method": "wb97xd",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "solvent_method": "pcm",
                "symmetry": "false",
                "sym_ignore": "true",
                "resp_charges": "true",
            },
        )
        self.assertEqual(
            test_FreqSet.pcm,
            {
                "heavypoints": "194",
                "hpoints": "194",
                "radii": "uff",
                "theory": "cpcm",
                "vdwscale": "1.1",
            },
        )
        self.assertEqual(test_FreqSet.solvent, {"dielectric": 10.0})
        self.assertEqual(test_FreqSet.molecule, test_molecule)

    def test_smd_init(self):
        test_molecule = QCInput.from_file(os.path.join(test_dir, "new_qchem_files/pcm.qin")).molecule
        test_FreqSet = FreqSet(molecule=test_molecule, smd_solvent="water")
        self.assertEqual(
            test_FreqSet.rem,
            {
                "job_type": "freq",
                "gen_scfman": "true",
                "basis": "def2-tzvppd",
                "max_scf_cycles": "200",
                "method": "wb97xd",
                "scf_algorithm": "diis",
                "xc_grid": "3",
                "solvent_method": "smd",
                "ideriv": "1",
                "symmetry": "false",
                "sym_ignore": "true",
                "resp_charges": "true",
            },
        )
        self.assertEqual(test_FreqSet.smx, {"solvent": "water"})
        self.assertEqual(test_FreqSet.molecule, test_molecule)


if __name__ == "__main__":
    unittest.main()
