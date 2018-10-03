# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import unittest
import os
import json
import gzip
import numpy as np
import warnings

from shutil import copyfile, copyfileobj
from monty.tempfile import ScratchDir

import xml.etree.cElementTree as ET

from pymatgen.core.periodic_table import Element
from pymatgen.electronic_structure.core import OrbitalType
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.outputs import Chgcar, Locpot, Oszicar, Outcar, \
    Vasprun, Procar, Xdatcar, Dynmat, BSVasprun, UnconvergedVASPWarning, \
    VaspParserError, Wavecar
from pymatgen import Spin, Orbital, Lattice, Structure
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.electronic_structure.core import Magmom
from pymatgen.util.testing import PymatgenTest

"""
Created on Jul 16, 2012
"""

__author__ = "Shyue Ping Ong, Stephen Dacek, Mark Turiansky"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Jul 16, 2012"

test_dir = os.path.abspath(
    os.path.join(os.path.dirname(__file__),
                 "..", "..", "..", "..",
                 'test_files'))


class VasprunTest(unittest.TestCase):
    _multiprocess_shared_ = True

    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.resetwarnings()

    def test_multiple_dielectric(self):
        v = Vasprun(os.path.join(test_dir, "vasprun.GW0.xml"))
        self.assertEqual(len(v.other_dielectric), 3)

    def test_charge_charge_dielectric(self):
        """
        VASP 5.4.4 writes out two dielectric functions to vasprun.xml
        These are the "density-density" and "velocity-velocity" linear response functions.
        See the comments in `linear_optics.F` for details.
        """
        v = Vasprun(os.path.join(test_dir, "vasprun.xml.dielectric_5.4.4"),
                    parse_potcar_file=False)
        self.assertEqual(v.dielectric is not None, True)
        self.assertEqual('density' in v.dielectric_data, True)
        self.assertEqual('velocity' in v.dielectric_data, True)

    def test_optical_absorption_coeff(self):
        v = Vasprun(os.path.join(test_dir, "vasprun.BSE.xml.gz"))
        absorption_coeff = v.optical_absorption_coeff
        self.assertEqual(absorption_coeff[1], 24966408728.917931)

    def test_vasprun_with_more_than_two_unlabelled_dielectric_functions(self):
        with self.assertRaises(NotImplementedError):
            Vasprun(os.path.join(test_dir, "vasprun.xml.dielectric_bad"),
                    parse_potcar_file=False)

    def test_bad_vasprun(self):
        self.assertRaises(ET.ParseError,
                          Vasprun, os.path.join(test_dir, "bad_vasprun.xml"))

        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # Trigger a warning.
            v = Vasprun(os.path.join(test_dir, "bad_vasprun.xml"),
                        exception_on_bad_xml=False)
            # Verify some things
            self.assertEqual(len(v.ionic_steps), 1)
            self.assertAlmostEqual(v.final_energy, -269.00551374)
            self.assertTrue(issubclass(w[-1].category,
                                       UserWarning))

    def test_vdw(self):
        v = Vasprun(os.path.join(test_dir, "vasprun.xml.vdw"))
        self.assertAlmostEqual(v.final_energy, -9.78310677)

    def test_nonlmn(self):

        filepath = os.path.join(test_dir, 'vasprun.xml.nonlm')
        vasprun = Vasprun(filepath, parse_potcar_file=False)
        orbs = list(vasprun.complete_dos.pdos[vasprun.final_structure[
            0]].keys())
        self.assertIn(OrbitalType.s, orbs)

    def test_standard(self):
        filepath = os.path.join(test_dir, 'vasprun.xml')
        vasprun = Vasprun(filepath, parse_potcar_file=False)

        # Test NELM parsing.
        self.assertEqual(vasprun.parameters["NELM"], 60)
        # test pdos parsing

        pdos0 = vasprun.complete_dos.pdos[vasprun.final_structure[0]]
        self.assertAlmostEqual(pdos0[Orbital.s][Spin.up][16], 0.0026)
        self.assertAlmostEqual(pdos0[Orbital.pz][Spin.down][16], 0.0012)
        self.assertEqual(pdos0[Orbital.s][Spin.up].shape, (301,))

        filepath2 = os.path.join(test_dir, 'lifepo4.xml')
        vasprun_ggau = Vasprun(filepath2, parse_projected_eigen=True,
                               parse_potcar_file=False)
        totalscsteps = sum([len(i['electronic_steps'])
                            for i in vasprun.ionic_steps])
        self.assertEqual(29, len(vasprun.ionic_steps))
        self.assertEqual(len(vasprun.structures), len(vasprun.ionic_steps))
        self.assertEqual(vasprun.lattice,
                         vasprun.lattice_rec.reciprocal_lattice)

        for i, step in enumerate(vasprun.ionic_steps):
            self.assertEqual(vasprun.structures[i], step["structure"])

        self.assertTrue(all([vasprun.structures[i] == vasprun.ionic_steps[i][
            "structure"] for i in range(len(vasprun.ionic_steps))]))

        self.assertEqual(308, totalscsteps,
                         "Incorrect number of energies read from vasprun.xml")

        self.assertEqual(['Li'] + 4 * ['Fe'] + 4 * ['P'] + 16 * ["O"],
                         vasprun.atomic_symbols)
        self.assertEqual(vasprun.final_structure.composition.reduced_formula,
                         "LiFe4(PO4)4")
        self.assertIsNotNone(vasprun.incar, "Incar cannot be read")
        self.assertIsNotNone(vasprun.kpoints, "Kpoints cannot be read")
        self.assertIsNotNone(vasprun.eigenvalues, "Eigenvalues cannot be read")
        self.assertAlmostEqual(vasprun.final_energy, -269.38319884, 7)
        self.assertAlmostEqual(vasprun.tdos.get_gap(), 2.0589, 4)
        expectedans = (2.539, 4.0906, 1.5516, False)
        (gap, cbm, vbm, direct) = vasprun.eigenvalue_band_properties
        self.assertAlmostEqual(gap, expectedans[0])
        self.assertAlmostEqual(cbm, expectedans[1])
        self.assertAlmostEqual(vbm, expectedans[2])
        self.assertEqual(direct, expectedans[3])
        self.assertFalse(vasprun.is_hubbard)
        self.assertEqual(vasprun.potcar_symbols,
                         ['PAW_PBE Li 17Jan2003', 'PAW_PBE Fe 06Sep2000',
                          'PAW_PBE Fe 06Sep2000', 'PAW_PBE P 17Jan2003',
                          'PAW_PBE O 08Apr2002'])
        self.assertIsNotNone(vasprun.kpoints, "Kpoints cannot be read")
        self.assertIsNotNone(vasprun.actual_kpoints,
                             "Actual kpoints cannot be read")
        self.assertIsNotNone(vasprun.actual_kpoints_weights,
                             "Actual kpoints weights cannot be read")
        for atomdoses in vasprun.pdos:
            for orbitaldos in atomdoses:
                self.assertIsNotNone(orbitaldos, "Partial Dos cannot be read")

        # test skipping ionic steps.
        vasprun_skip = Vasprun(filepath, 3, parse_potcar_file=False)
        self.assertEqual(vasprun_skip.nionic_steps, 29)
        self.assertEqual(len(vasprun_skip.ionic_steps),
                         int(vasprun.nionic_steps / 3) + 1)
        self.assertEqual(len(vasprun_skip.ionic_steps),
                         len(vasprun_skip.structures))
        self.assertEqual(len(vasprun_skip.ionic_steps),
                         int(vasprun.nionic_steps / 3) + 1)
        # Check that nionic_steps is preserved no matter what.
        self.assertEqual(vasprun_skip.nionic_steps,
                         vasprun.nionic_steps)

        self.assertNotAlmostEqual(vasprun_skip.final_energy,
                                  vasprun.final_energy)

        # Test with ionic_step_offset
        vasprun_offset = Vasprun(filepath, 3, 6, parse_potcar_file=False)
        self.assertEqual(len(vasprun_offset.ionic_steps),
                         int(len(vasprun.ionic_steps) / 3) - 1)
        self.assertEqual(vasprun_offset.structures[0],
                         vasprun_skip.structures[2])

        self.assertTrue(vasprun_ggau.is_hubbard)
        self.assertEqual(vasprun_ggau.hubbards["Fe"], 4.3)
        self.assertAlmostEqual(vasprun_ggau.projected_eigenvalues[Spin.up][
            0][0][96][0], 0.0032)
        d = vasprun_ggau.as_dict()
        self.assertEqual(d["elements"], ["Fe", "Li", "O", "P"])
        self.assertEqual(d["nelements"], 4)

    def test_unconverged(self):
        filepath = os.path.join(test_dir, 'vasprun.xml.unconverged')
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # Trigger a warning.
            vasprun_unconverged = Vasprun(filepath, parse_potcar_file=False)
            # Verify some things
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[-1].category,
                                       UnconvergedVASPWarning))

            self.assertTrue(vasprun_unconverged.converged_ionic)
            self.assertFalse(vasprun_unconverged.converged_electronic)
            self.assertFalse(vasprun_unconverged.converged)

    def test_dfpt(self):
        filepath = os.path.join(test_dir, 'vasprun.xml.dfpt')
        vasprun_dfpt = Vasprun(filepath, parse_potcar_file=False)
        self.assertAlmostEqual(vasprun_dfpt.epsilon_static[0][0], 3.26105533)
        self.assertAlmostEqual(vasprun_dfpt.epsilon_static[0][1], -0.00459066)
        self.assertAlmostEqual(vasprun_dfpt.epsilon_static[2][2], 3.24330517)
        self.assertAlmostEqual(vasprun_dfpt.epsilon_static_wolfe[0][0],
                               3.33402531)
        self.assertAlmostEqual(vasprun_dfpt.epsilon_static_wolfe[0][1],
                               -0.00559998)
        self.assertAlmostEqual(vasprun_dfpt.epsilon_static_wolfe[2][2],
                               3.31237357)
        self.assertTrue(vasprun_dfpt.converged)

        entry = vasprun_dfpt.get_computed_entry()
        entry = MaterialsProjectCompatibility(
            check_potcar_hash=False).process_entry(entry)
        self.assertAlmostEqual(entry.uncorrected_energy + entry.correction,
                               entry.energy)

    def test_dfpt_ionic(self):
        filepath = os.path.join(test_dir, 'vasprun.xml.dfpt.ionic')
        vasprun_dfpt_ionic = Vasprun(filepath, parse_potcar_file=False)
        self.assertAlmostEqual(vasprun_dfpt_ionic.epsilon_ionic[0][0],
                               515.73485838)
        self.assertAlmostEqual(vasprun_dfpt_ionic.epsilon_ionic[0][1],
                               -0.00263523)
        self.assertAlmostEqual(vasprun_dfpt_ionic.epsilon_ionic[2][2],
                               19.02110169)

    def test_dfpt_unconverged(self):
        filepath = os.path.join(test_dir, 'vasprun.xml.dfpt.unconverged')
        vasprun_dfpt_unconv = Vasprun(filepath, parse_potcar_file=False)
        self.assertFalse(vasprun_dfpt_unconv.converged_electronic)
        self.assertTrue(vasprun_dfpt_unconv.converged_ionic)
        self.assertFalse(vasprun_dfpt_unconv.converged)

    def test_uniform(self):
        vasprun_uniform = Vasprun(os.path.join(test_dir, "vasprun.xml.uniform"),
                                  parse_potcar_file=False)
        self.assertEqual(vasprun_uniform.kpoints.style,
                         Kpoints.supported_modes.Reciprocal)

    def test_no_projected(self):
        vasprun_no_pdos = Vasprun(os.path.join(test_dir, "Li_no_projected.xml"),
                                  parse_potcar_file=False)
        self.assertIsNotNone(vasprun_no_pdos.complete_dos)
        self.assertFalse(vasprun_no_pdos.dos_has_errors)

    def test_dielectric(self):
        vasprun_diel = Vasprun(os.path.join(test_dir, "vasprun.xml.dielectric"),
                               parse_potcar_file=False)
        self.assertAlmostEqual(0.4294, vasprun_diel.dielectric[0][10])
        self.assertAlmostEqual(19.941, vasprun_diel.dielectric[1][51][0])
        self.assertAlmostEqual(19.941, vasprun_diel.dielectric[1][51][1])
        self.assertAlmostEqual(19.941, vasprun_diel.dielectric[1][51][2])
        self.assertAlmostEqual(0.0, vasprun_diel.dielectric[1][51][3])
        self.assertAlmostEqual(34.186, vasprun_diel.dielectric[2][85][0])
        self.assertAlmostEqual(34.186, vasprun_diel.dielectric[2][85][1])
        self.assertAlmostEqual(34.186, vasprun_diel.dielectric[2][85][2])
        self.assertAlmostEqual(0.0, vasprun_diel.dielectric[2][85][3])

    def test_indirect_vasprun(self):
        v = Vasprun(os.path.join(test_dir, "vasprun.xml.indirect.gz"))
        (gap, cbm, vbm, direct) = v.eigenvalue_band_properties
        self.assertFalse(direct)

    def test_optical_vasprun(self):
        vasprun_optical = Vasprun(
            os.path.join(test_dir, "vasprun.xml.opticaltransitions"),
            parse_potcar_file=False)
        self.assertAlmostEqual(3.084, vasprun_optical.optical_transition[0][0])
        self.assertAlmostEqual(3.087, vasprun_optical.optical_transition[3][0])
        self.assertAlmostEqual(0.001, vasprun_optical.optical_transition[0][1])
        self.assertAlmostEqual(0.001, vasprun_optical.optical_transition[1][1])
        self.assertAlmostEqual(0.001, vasprun_optical.optical_transition[7][1])
        self.assertAlmostEqual(0.001, vasprun_optical.optical_transition[19][1])
        self.assertAlmostEqual(3.3799999999,
                               vasprun_optical.optical_transition[54][0])
        self.assertAlmostEqual(3.381, vasprun_optical.optical_transition[55][0])
        self.assertAlmostEqual(3.381, vasprun_optical.optical_transition[56][0])
        self.assertAlmostEqual(10554.9860,
                               vasprun_optical.optical_transition[54][1])
        self.assertAlmostEqual(0.0, vasprun_optical.optical_transition[55][1])
        self.assertAlmostEqual(0.001, vasprun_optical.optical_transition[56][1])

    def test_force_constants(self):
        vasprun_fc = Vasprun(os.path.join(test_dir, "vasprun.xml.dfpt.phonon"),
                             parse_potcar_file=False)
        fc_ans = [[-0.00184451, -0., -0.],
                  [-0., -0.00933824, -0.03021279],
                  [-0., -0.03021279, 0.01202547]]
        nm_ans = [[0.0884346, -0.08837289, -0.24995639],
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
                  [-0.15306645, 0.05105771, 0.14441306]]
        nm_eigenval_ans = [-0.59067079, -0.59067079, -0.59067003, -0.59067003,
                           -0.59067003, -0.59067003, -0.585009, -0.585009,
                           -0.58500895, -0.58500883, -0.5062956, -0.5062956]
        self.assertEqual(vasprun_fc.force_constants.shape, (16, 16, 3, 3))
        self.assertTrue(np.allclose(vasprun_fc.force_constants[8, 9], fc_ans))
        self.assertEqual(vasprun_fc.normalmode_eigenvals.size, 48)
        self.assertTrue(np.allclose(vasprun_fc.normalmode_eigenvals[17:29],
                                    nm_eigenval_ans))
        self.assertEqual(vasprun_fc.normalmode_eigenvecs.shape, (48, 16, 3))
        self.assertTrue(
            np.allclose(vasprun_fc.normalmode_eigenvecs[33], nm_ans))

    def test_Xe(self):
        vr = Vasprun(os.path.join(test_dir, 'vasprun.xml.xe'),
                     parse_potcar_file=False)
        self.assertEqual(vr.atomic_symbols, ['Xe'])

    def test_invalid_element(self):
        self.assertRaises(ValueError, Vasprun,
                          os.path.join(test_dir, 'vasprun.xml.wrong_sp'))

    def test_selective_dynamics(self):
        vsd = Vasprun(os.path.join(test_dir, 'vasprun.xml.indirect.gz'))
        np.testing.assert_array_equal(
            vsd.final_structure.site_properties.get('selective_dynamics'),
            [[True] * 3, [False] * 3], "Selective dynamics parsing error")

    def test_as_dict(self):
        filepath = os.path.join(test_dir, 'vasprun.xml')
        vasprun = Vasprun(filepath,
                          parse_potcar_file=False)
        # Test that as_dict() is json-serializable
        self.assertIsNotNone(json.dumps(vasprun.as_dict()))
        self.assertEqual(
            vasprun.as_dict()["input"]["potcar_type"],
            ['PAW_PBE', 'PAW_PBE', 'PAW_PBE', 'PAW_PBE', 'PAW_PBE'])
        self.assertEqual(vasprun.as_dict()['input']['nkpoints'], 24)

    def test_get_band_structure(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            filepath = os.path.join(test_dir, 'vasprun_Si_bands.xml')
            vasprun = Vasprun(filepath,
                              parse_projected_eigen=True,
                              parse_potcar_file=False)
            bs = vasprun.get_band_structure(kpoints_filename=os.path.join(test_dir,
                                                                          'KPOINTS_Si_bands'))
            cbm = bs.get_cbm()
            vbm = bs.get_vbm()
            self.assertEqual(cbm['kpoint_index'], [13],
                             "wrong cbm kpoint index")
            self.assertAlmostEqual(cbm['energy'], 6.2301, "wrong cbm energy")
            self.assertEqual(cbm['band_index'], {Spin.up: [4], Spin.down: [4]},
                             "wrong cbm bands")
            self.assertEqual(vbm['kpoint_index'], [0, 63, 64])
            self.assertAlmostEqual(vbm['energy'], 5.6158, "wrong vbm energy")
            self.assertEqual(vbm['band_index'], {Spin.up: [1, 2, 3],
                                                 Spin.down: [1, 2, 3]},
                             "wrong vbm bands")
            self.assertEqual(vbm['kpoint'].label, "\\Gamma", "wrong vbm label")
            self.assertEqual(cbm['kpoint'].label, None, "wrong cbm label")

            projected = bs.get_projection_on_elements()
            self.assertAlmostEqual(projected[Spin.up][0][0]["Si"], 0.4238)
            projected = bs.get_projections_on_elements_and_orbitals(
                {"Si": ["s"]})
            self.assertAlmostEqual(projected[Spin.up][0][0]["Si"]["s"], 0.4238)

            # Test compressed files case 1: compressed KPOINTS in current dir
            with ScratchDir("./"):
                copyfile(os.path.join(test_dir, 'vasprun_Si_bands.xml'),
                         'vasprun.xml')

                # Check for error if no KPOINTS file
                vasprun = Vasprun('vasprun.xml',
                                  parse_projected_eigen=True,
                                  parse_potcar_file=False)
                with self.assertRaises(VaspParserError):
                    _ = vasprun.get_band_structure(line_mode=True)

                # Check KPOINTS.gz succesfully inferred and used if present
                with open(os.path.join(test_dir, 'KPOINTS_Si_bands'),
                          'rb') as f_in:
                    with gzip.open('KPOINTS.gz', 'wb') as f_out:
                        copyfileobj(f_in, f_out)
                bs_kpts_gzip = vasprun.get_band_structure()
                self.assertEqual(bs.efermi, bs_kpts_gzip.efermi)
                self.assertEqual(bs.as_dict(), bs_kpts_gzip.as_dict())

            # Test compressed files case 2: compressed vasprun in another dir
            with ScratchDir("./"):
                os.mkdir('deeper')
                copyfile(os.path.join(test_dir, 'KPOINTS_Si_bands'),
                         os.path.join('deeper', 'KPOINTS'))
                with open(os.path.join(test_dir, 'vasprun_Si_bands.xml'),
                          'rb') as f_in:
                    with gzip.open(os.path.join('deeper', 'vasprun.xml.gz'),
                                   'wb') as f_out:
                        copyfileobj(f_in, f_out)
                vasprun = Vasprun(os.path.join('deeper', 'vasprun.xml.gz'),
                                  parse_projected_eigen=True,
                                  parse_potcar_file=False)
                bs_vasprun_gzip = vasprun.get_band_structure(line_mode=True)
                self.assertEqual(bs.efermi, bs_vasprun_gzip.efermi)
                self.assertEqual(bs.as_dict(), bs_vasprun_gzip.as_dict())


            # test hybrid band structures
            vasprun.actual_kpoints_weights[-1] = 0.
            bs = vasprun.get_band_structure(kpoints_filename=os.path.join(test_dir,
                                                                          'KPOINTS_Si_bands'))
            cbm = bs.get_cbm()
            vbm = bs.get_vbm()
            self.assertEqual(cbm['kpoint_index'], [0])
            self.assertAlmostEqual(cbm['energy'], 6.3676)
            self.assertEqual(cbm['kpoint'].label, None)
            self.assertEqual(vbm['kpoint_index'], [0])
            self.assertAlmostEqual(vbm['energy'], 2.8218)
            self.assertEqual(vbm['kpoint'].label, None)

    def test_sc_step_overflow(self):
        filepath = os.path.join(test_dir, 'vasprun.xml.sc_overflow')
        # with warnings.catch_warnings(record=True) as w:
        #     warnings.simplefilter("always")
        #     vasprun = Vasprun(filepath)
        #     self.assertEqual(len(w), 3)
        vasprun = Vasprun(filepath)
        estep = vasprun.ionic_steps[0]['electronic_steps'][29]
        self.assertTrue(np.isnan(estep['e_wo_entrp']))

    def test_update_potcar(self):
        filepath = os.path.join(test_dir, 'vasprun.xml')
        potcar_path = os.path.join(test_dir, 'POTCAR.LiFePO4.gz')
        potcar_path2 = os.path.join(test_dir, 'POTCAR2.LiFePO4.gz')
        vasprun = Vasprun(filepath, parse_potcar_file=False)
        self.assertEqual(vasprun.potcar_spec,
                         [{"titel": "PAW_PBE Li 17Jan2003", "hash": None},
                          {"titel": "PAW_PBE Fe 06Sep2000", "hash": None},
                          {"titel": "PAW_PBE Fe 06Sep2000", "hash": None},
                          {"titel": "PAW_PBE P 17Jan2003", "hash": None},
                          {"titel": "PAW_PBE O 08Apr2002", "hash": None}])

        vasprun.update_potcar_spec(potcar_path)
        self.assertEqual(vasprun.potcar_spec, [{"titel": "PAW_PBE Li 17Jan2003",
                                                "hash": "65e83282d1707ec078c1012afbd05be8"},
                                               {"titel": "PAW_PBE Fe 06Sep2000",
                                                "hash": "9530da8244e4dac17580869b4adab115"},
                                               {"titel": "PAW_PBE Fe 06Sep2000",
                                                "hash": "9530da8244e4dac17580869b4adab115"},
                                               {"titel": "PAW_PBE P 17Jan2003",
                                                "hash": "7dc3393307131ae67785a0cdacb61d5f"},
                                               {"titel": "PAW_PBE O 08Apr2002",
                                                "hash": "7a25bc5b9a5393f46600a4939d357982"}])

        vasprun2 = Vasprun(filepath, parse_potcar_file=False)
        self.assertRaises(ValueError, vasprun2.update_potcar_spec, potcar_path2)
        vasprun = Vasprun(filepath, parse_potcar_file=potcar_path)

        self.assertEqual(vasprun.potcar_spec, [{"titel": "PAW_PBE Li 17Jan2003",
                                                "hash": "65e83282d1707ec078c1012afbd05be8"},
                                               {"titel": "PAW_PBE Fe 06Sep2000",
                                                "hash": "9530da8244e4dac17580869b4adab115"},
                                               {"titel": "PAW_PBE Fe 06Sep2000",
                                                "hash": "9530da8244e4dac17580869b4adab115"},
                                               {"titel": "PAW_PBE P 17Jan2003",
                                                "hash": "7dc3393307131ae67785a0cdacb61d5f"},
                                               {"titel": "PAW_PBE O 08Apr2002",
                                                "hash": "7a25bc5b9a5393f46600a4939d357982"}])

        self.assertRaises(ValueError, Vasprun, filepath,
                          parse_potcar_file=potcar_path2)

    def test_search_for_potcar(self):
        filepath = os.path.join(test_dir, 'vasprun.xml')
        vasprun = Vasprun(filepath, parse_potcar_file=True)
        self.assertEqual(vasprun.potcar_spec, [{"titel": "PAW_PBE Li 17Jan2003",
                                                "hash": "65e83282d1707ec078c1012afbd05be8"},
                                               {"titel": "PAW_PBE Fe 06Sep2000",
                                                "hash": "9530da8244e4dac17580869b4adab115"},
                                               {"titel": "PAW_PBE Fe 06Sep2000",
                                                "hash": "9530da8244e4dac17580869b4adab115"},
                                               {"titel": "PAW_PBE P 17Jan2003",
                                                "hash": "7dc3393307131ae67785a0cdacb61d5f"},
                                               {"titel": "PAW_PBE O 08Apr2002",
                                                "hash": "7a25bc5b9a5393f46600a4939d357982"}])

    def test_potcar_not_found(self):
        filepath = os.path.join(test_dir, 'vasprun.xml')
        # Ensure no potcar is found and nothing is updated
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            vasprun = Vasprun(filepath, parse_potcar_file='.')
            self.assertEqual(len(w), 2)
            self.assertEqual(vasprun.potcar_spec, [{"titel": "PAW_PBE Li 17Jan2003", "hash": None},
                                                   {"titel": "PAW_PBE Fe 06Sep2000", "hash": None},
                                                   {"titel": "PAW_PBE Fe 06Sep2000", "hash": None},
                                                   {"titel": "PAW_PBE P 17Jan2003", "hash": None},
                                                   {"titel": "PAW_PBE O 08Apr2002", "hash": None}])

    def test_parsing_chemical_shift_calculations(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            filepath = os.path.join(test_dir, "nmr", "cs", "basic",
                                    'vasprun.xml.chemical_shift.scstep')
            vasprun = Vasprun(filepath)
            nestep = len(vasprun.ionic_steps[-1]['electronic_steps'])
            self.assertEqual(nestep, 10)
            self.assertTrue(vasprun.converged)

    def test_parsing_efg_calcs(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            filepath = os.path.join(test_dir,  "nmr", "efg", "AlPO4",
                                    'vasprun.xml')
            vasprun = Vasprun(filepath)
            nestep = len(vasprun.ionic_steps[-1]['electronic_steps'])
            self.assertEqual(nestep, 18)
            self.assertTrue(vasprun.converged)

    def test_charged_structure(self):
        vpath = os.path.join(test_dir, 'vasprun.charged.xml')
        potcar_path = os.path.join(test_dir, 'POT_GGA_PAW_PBE', 'POTCAR.Si.gz')
        vasprun = Vasprun(vpath, parse_potcar_file=False)
        vasprun.update_charge_from_potcar(potcar_path)
        self.assertEqual(vasprun.parameters.get("NELECT", 8), 9)
        self.assertEqual(vasprun.structures[0].charge, 1)

        vpath = os.path.join(test_dir, 'vasprun.split.charged.xml')
        potcar_path = os.path.join(test_dir, 'POTCAR.split.charged.gz')
        vasprun = Vasprun(vpath, parse_potcar_file=False)
        vasprun.update_charge_from_potcar(potcar_path)
        self.assertEqual(vasprun.parameters.get('NELECT', 0), 7)
        self.assertEqual(vasprun.structures[-1].charge, 1)


class OutcarTest(PymatgenTest):

    _multiprocess_shared_ = True

    def test_init(self):
        for f in ['OUTCAR', 'OUTCAR.gz']:
            filepath = os.path.join(test_dir, f)
            outcar = Outcar(filepath)
            expected_mag = ({'d': 0.0, 'p': 0.003, 's': 0.002, 'tot': 0.005},
                            {'d': 0.798, 'p': 0.008, 's': 0.007, 'tot': 0.813},
                            {'d': 0.798, 'p': 0.008, 's': 0.007, 'tot': 0.813},
                            {'d': 0.0, 'p': -0.117, 's': 0.005, 'tot': -0.112},
                            {'d': 0.0, 'p': -0.165, 's': 0.004, 'tot': -0.162},
                            {'d': 0.0, 'p': -0.117, 's': 0.005, 'tot': -0.112},
                            {'d': 0.0, 'p': -0.165, 's': 0.004, 'tot': -0.162})
            expected_chg = ({'p': 0.154, 's': 0.078, 'd': 0.0, 'tot': 0.232},
                            {'p': 0.707, 's': 0.463, 'd': 8.316, 'tot': 9.486},
                            {'p': 0.707, 's': 0.463, 'd': 8.316, 'tot': 9.486},
                            {'p': 3.388, 's': 1.576, 'd': 0.0, 'tot': 4.964},
                            {'p': 3.365, 's': 1.582, 'd': 0.0, 'tot': 4.947},
                            {'p': 3.388, 's': 1.576, 'd': 0.0, 'tot': 4.964},
                            {'p': 3.365, 's': 1.582, 'd': 0.0, 'tot': 4.947})

            self.assertAlmostEqual(outcar.magnetization, expected_mag, 5,
                                   "Wrong magnetization read from Outcar")
            self.assertAlmostEqual(outcar.charge, expected_chg, 5,
                                   "Wrong charge read from Outcar")
            self.assertFalse(outcar.is_stopped)
            self.assertEqual(outcar.run_stats, {'System time (sec)': 0.938,
                                                'Total CPU time used (sec)': 545.142,
                                                'Elapsed time (sec)': 546.709,
                                                'Maximum memory used (kb)': 0.0,
                                                'Average memory used (kb)': 0.0,
                                                'User time (sec)': 544.204,
                                                'cores': '8'})
            self.assertAlmostEqual(outcar.efermi, 2.0112)
            self.assertAlmostEqual(outcar.nelect, 44.9999991)
            self.assertAlmostEqual(outcar.total_mag, 0.9999998)

            self.assertIsNotNone(outcar.as_dict())

            self.assertFalse(outcar.lepsilon)

    def test_stopped(self):
        filepath = os.path.join(test_dir, 'OUTCAR.stopped')
        outcar = Outcar(filepath)
        self.assertTrue(outcar.is_stopped)
        for f in ['OUTCAR.lepsilon', 'OUTCAR.lepsilon.gz']:
            filepath = os.path.join(test_dir, f)
            outcar = Outcar(filepath)

            self.assertTrue(outcar.lepsilon)
            self.assertAlmostEqual(outcar.dielectric_tensor[0][0], 3.716432)
            self.assertAlmostEqual(outcar.dielectric_tensor[0][1], -0.20464)
            self.assertAlmostEqual(outcar.dielectric_tensor[1][2], -0.20464)
            self.assertAlmostEqual(outcar.dielectric_ionic_tensor[0][0],
                                   0.001419)
            self.assertAlmostEqual(outcar.dielectric_ionic_tensor[0][2],
                                   0.001419)
            self.assertAlmostEqual(outcar.dielectric_ionic_tensor[2][2],
                                   0.001419)
            self.assertAlmostEqual(outcar.piezo_tensor[0][0], 0.52799)
            self.assertAlmostEqual(outcar.piezo_tensor[1][3], 0.35998)
            self.assertAlmostEqual(outcar.piezo_tensor[2][5], 0.35997)
            self.assertAlmostEqual(outcar.piezo_ionic_tensor[0][0], 0.05868)
            self.assertAlmostEqual(outcar.piezo_ionic_tensor[1][3], 0.06241)
            self.assertAlmostEqual(outcar.piezo_ionic_tensor[2][5], 0.06242)
            self.assertAlmostEqual(outcar.born[0][1][2], -0.385)
            self.assertAlmostEqual(outcar.born[1][2][0], 0.36465)
            self.assertAlmostEqual(outcar.internal_strain_tensor[0][0][0], -572.5437,places=4)
            self.assertAlmostEqual(outcar.internal_strain_tensor[0][1][0], 683.2985,places=4)
            self.assertAlmostEqual(outcar.internal_strain_tensor[0][1][3], 73.07059,places=4)
            self.assertAlmostEqual(outcar.internal_strain_tensor[1][0][0], 570.98927,places=4)
            self.assertAlmostEqual(outcar.internal_strain_tensor[1][1][0], -683.68519,places=4)
            self.assertAlmostEqual(outcar.internal_strain_tensor[1][2][2], 570.98927,places=4)

    def test_soc(self):
        filepath = os.path.join(test_dir, 'OUTCAR.NiO_SOC.gz')
        outcar = Outcar(filepath)
        expected_mag = (
            {'s': Magmom([0.0, 0.0, -0.001]), 'p': Magmom([0.0, 0.0, -0.003]),
             'd': Magmom([0.0, 0.0, 1.674]), 'tot': Magmom([0.0, 0.0, 1.671])},
            {'s': Magmom([0.0, 0.0, 0.001]), 'p': Magmom([0.0, 0.0, 0.003]),
             'd': Magmom([0.0, 0.0, -1.674]),
             'tot': Magmom([0.0, 0.0, -1.671])},
            {'s': Magmom([0.0, 0.0, 0.0]), 'p': Magmom([0.0, 0.0, 0.0]),
             'd': Magmom([0.0, 0.0, 0.0]), 'tot': Magmom([0.0, 0.0, 0.0])},
            {'s': Magmom([0.0, 0.0, 0.0]), 'p': Magmom([0.0, 0.0, 0.0]),
             'd': Magmom([0.0, 0.0, 0.0]), 'tot': Magmom([0.0, 0.0, 0.0])}
        )
        # test note: Magmom class uses np.allclose() when testing for equality
        # so fine to use assertEqual here
        self.assertEqual(outcar.magnetization, expected_mag,
                         "Wrong vector magnetization read from Outcar for SOC calculation")

    def test_polarization(self):
        filepath = os.path.join(test_dir, "OUTCAR.BaTiO3.polar")
        outcar = Outcar(filepath)
        self.assertEqual(outcar.spin, True)
        self.assertEqual(outcar.noncollinear, False)
        self.assertAlmostEqual(outcar.p_ion[0], 0.0)
        self.assertAlmostEqual(outcar.p_ion[1], 0.0)
        self.assertAlmostEqual(outcar.p_ion[2], -5.56684)
        self.assertAlmostEqual(outcar.p_sp1[0], 2.00068)
        self.assertAlmostEqual(outcar.p_sp2[0], -2.00044)
        self.assertAlmostEqual(outcar.p_elec[0], 0.00024)
        self.assertAlmostEqual(outcar.p_elec[1], 0.00019)
        self.assertAlmostEqual(outcar.p_elec[2], 3.61674)

    def test_pseudo_zval(self):
        filepath = os.path.join(test_dir, "OUTCAR.BaTiO3.polar")
        outcar = Outcar(filepath)
        self.assertDictEqual({'Ba': 10.00, 'Ti': 10.00, 'O': 6.00},
                             outcar.zval_dict)

    def test_dielectric(self):
        filepath = os.path.join(test_dir, "OUTCAR.dielectric")
        outcar = Outcar(filepath)
        outcar.read_corrections()
        self.assertAlmostEqual(outcar.data["dipol_quadrupol_correction"],
                               0.03565)
        self.assertAlmostEqual(outcar.final_energy, -797.46760559)

    def test_freq_dielectric(self):
        filepath = os.path.join(test_dir, "OUTCAR.LOPTICS")
        outcar = Outcar(filepath)
        outcar.read_freq_dielectric()
        self.assertAlmostEqual(outcar.frequencies[0], 0)
        self.assertAlmostEqual(outcar.frequencies[-1], 39.826101)
        self.assertAlmostEqual(outcar.dielectric_tensor_function[0][0, 0],
                               8.96938800)
        self.assertAlmostEqual(outcar.dielectric_tensor_function[-1][0, 0],
                               7.36167000e-01 + 1.53800000e-03j)
        self.assertEqual(len(outcar.frequencies),
                         len(outcar.dielectric_tensor_function))
        np.testing.assert_array_equal(outcar.dielectric_tensor_function[0],
                                      outcar.dielectric_tensor_function[
                                          0].transpose())

    def test_freq_dielectric_vasp544(self):
        filepath = os.path.join(test_dir, "OUTCAR.LOPTICS.vasp544")
        outcar = Outcar(filepath)
        outcar.read_freq_dielectric()
        self.assertAlmostEqual(outcar.frequencies[0], 0)
        self.assertAlmostEqual(outcar.frequencies[-1], 39.63964)
        self.assertAlmostEqual(outcar.dielectric_tensor_function[0][0, 0],
                               12.769435 + 0j)
        self.assertAlmostEqual(outcar.dielectric_tensor_function[-1][0, 0],
                               0.828615 + 0.016594j)
        self.assertEqual(len(outcar.frequencies),
                         len(outcar.dielectric_tensor_function))
        np.testing.assert_array_equal(outcar.dielectric_tensor_function[0],
                                      outcar.dielectric_tensor_function[
                                          0].transpose())

    def test_read_elastic_tensor(self):
        filepath = os.path.join(test_dir, "OUTCAR.total_tensor.Li2O.gz")
        outcar = Outcar(filepath)

        outcar.read_elastic_tensor()

        self.assertAlmostEqual(outcar.data["elastic_tensor"][0][0], 1986.3391)
        self.assertAlmostEqual(outcar.data["elastic_tensor"][0][1], 187.8324)
        self.assertAlmostEqual(outcar.data["elastic_tensor"][3][3], 586.3034)

    def test_read_piezo_tensor(self):
        filepath = os.path.join(test_dir, "OUTCAR.lepsilon.gz")
        outcar = Outcar(filepath)

        outcar.read_piezo_tensor()
        self.assertAlmostEqual(outcar.data["piezo_tensor"][0][0], 0.52799)
        self.assertAlmostEqual(outcar.data["piezo_tensor"][1][3], 0.35998)
        self.assertAlmostEqual(outcar.data["piezo_tensor"][2][5], 0.35997)

    def test_core_state_eigen(self):
        filepath = os.path.join(test_dir, "OUTCAR.CL")
        cl = Outcar(filepath).read_core_state_eigen()
        self.assertAlmostEqual(cl[6]["2s"][-1], -174.4779)
        filepath = os.path.join(test_dir, "OUTCAR.icorelevel")
        cl = Outcar(filepath).read_core_state_eigen()
        self.assertAlmostEqual(cl[4]["3d"][-1], -31.4522)

    def test_avg_core_poten(self):
        filepath = os.path.join(test_dir, "OUTCAR.lepsilon")
        cp = Outcar(filepath).read_avg_core_poten()
        self.assertAlmostEqual(cp[-1][1], -90.0487)
        filepath = os.path.join(test_dir, "OUTCAR")
        cp = Outcar(filepath).read_avg_core_poten()
        self.assertAlmostEqual(cp[0][6], -73.1068)

    def test_single_atom(self):
        filepath = os.path.join(test_dir, "OUTCAR.Al")
        outcar = Outcar(filepath)
        expected_mag = ({u'p': 0.0, u's': 0.0, u'd': 0.0, u'tot': 0.0},)
        expected_chg = ({u'p': 0.343, u's': 0.425, u'd': 0.0, u'tot': 0.768},)

        self.assertAlmostEqual(outcar.magnetization, expected_mag)
        self.assertAlmostEqual(outcar.charge, expected_chg)
        self.assertFalse(outcar.is_stopped)
        self.assertEqual(outcar.run_stats, {'System time (sec)': 0.592,
                                            'Total CPU time used (sec)': 50.194,
                                            'Elapsed time (sec)': 52.337,
                                            'Maximum memory used (kb)': 62900.0,
                                            'Average memory used (kb)': 0.0,
                                            'User time (sec)': 49.602,
                                            'cores': '32'})
        self.assertAlmostEqual(outcar.efermi, 8.0942)
        self.assertAlmostEqual(outcar.nelect, 3)
        self.assertAlmostEqual(outcar.total_mag, 8.2e-06)

        self.assertIsNotNone(outcar.as_dict())

    def test_chemical_shielding(self):
        filename = os.path.join(test_dir, "nmr", "cs", "core.diff",
                                "hydromagnesite", "OUTCAR")
        outcar = Outcar(filename)
        expected_chemical_shielding = [[191.9974, 69.5232, 0.6342],
                                       [195.0808, 68.183, 0.833],
                                       [192.0389, 69.5762, 0.6329],
                                       [195.0844, 68.1756, 0.8336],
                                       [192.005, 69.5289, 0.6339],
                                       [195.0913, 68.1859, 0.833],
                                       [192.0237, 69.565, 0.6333],
                                       [195.0788, 68.1733, 0.8337]]

        self.assertAlmostEqual(
            len(outcar.data["chemical_shielding"]["valence_only"][20: 28]),
            len(expected_chemical_shielding))

        self.assertArrayAlmostEqual(outcar.data["chemical_shielding"]["valence_and_core"][20:28],
                                    expected_chemical_shielding, decimal=5)

    def test_chemical_shielding_with_different_core_contribution(self):
        filename = os.path.join(test_dir, "nmr", "cs", "core.diff",
                                "core.diff.chemical.shifts.OUTCAR")
        outcar = Outcar(filename)
        c_vo = outcar.data["chemical_shielding"]["valence_only"][7]
        for x1, x2 in zip(list(c_vo),
                          [198.7009, 73.7484, 1.0000]):
            self.assertAlmostEqual(x1, x2)
        c_vc = outcar.data["chemical_shielding"]["valence_and_core"][7]
        for x1, x2 in zip(list(c_vc),
                          [-1.9406, 73.7484, 1.0000]):
            self.assertAlmostEqual(x1, x2)

    def test_cs_raw_tensors(self):
        filename = os.path.join(test_dir, "nmr", "cs", "core.diff",
                                "core.diff.chemical.shifts.OUTCAR")
        outcar = Outcar(filename)
        unsym_tensors = outcar.data["unsym_cs_tensor"]
        self.assertEqual(unsym_tensors[0],
                         [[-145.814605, -4.263425, 0.000301],
                          [4.263434, -145.812238, -8.7e-05],
                          [0.000136, -0.000189, -142.794068]])
        self.assertEqual(unsym_tensors[29],
                         [[287.789318, -53.799325, 30.900024],
                          [-53.799571, 225.668117, -17.839598],
                          [3.801103, -2.195218, 88.896756]])

    def test_cs_g0_contribution(self):
        filename = os.path.join(test_dir, "nmr", "cs", "core.diff",
                                "core.diff.chemical.shifts.OUTCAR")
        outcar = Outcar(filename)
        g0_contrib = outcar.data["cs_g0_contribution"]
        self.assertEqual(g0_contrib,
                         [[-8.773535, 9e-06, 1e-06],
                          [1.7e-05, -8.773536, -0.0792],
                          [-6e-06, -0.008328, -9.320237]])

    def test_cs_core_contribution(self):
        filename = os.path.join(test_dir, "nmr", "cs", "core.diff",
                                "core.diff.chemical.shifts.OUTCAR")
        outcar = Outcar(filename)
        core_contrib = outcar.data["cs_core_contribution"]
        self.assertEqual(core_contrib,
                         {'Mg': -412.8248405,
                          'C': -200.5098812,
                          'O': -271.0766979})

    def test_nmr_efg(self):
        filename = os.path.join(test_dir, "nmr", "efg", "AlPO4", "OUTCAR")
        outcar = Outcar(filename)
        expected_efg = [
            {'eta': 0.465, 'nuclear_quadrupole_moment': 146.6, 'cq': -5.573},
            {'eta': 0.465, 'nuclear_quadrupole_moment': 146.6, 'cq': -5.573},
            {'eta': 0.137, 'nuclear_quadrupole_moment': 146.6, 'cq': 6.327},
            {'eta': 0.137, 'nuclear_quadrupole_moment': 146.6, 'cq': 6.327},
            {'eta': 0.112, 'nuclear_quadrupole_moment': 146.6, 'cq': -7.453},
            {'eta': 0.112, 'nuclear_quadrupole_moment': 146.6, 'cq': -7.453},
            {'eta': 0.42, 'nuclear_quadrupole_moment': 146.6, 'cq': -5.58},
            {'eta': 0.42, 'nuclear_quadrupole_moment': 146.6, 'cq': -5.58}]
        self.assertEqual(len(outcar.data["efg"][2:10]), len(expected_efg))
        for e1, e2 in zip(outcar.data["efg"][2:10], expected_efg):
            for k in e1.keys():
                self.assertAlmostEqual(e1[k], e2[k], places=5)

        exepected_tensors = [[[11.11, 1.371, 2.652], [1.371, 3.635, -3.572], [2.652, -3.572, -14.746]],
                             [[11.11, -1.371, 2.652], [-1.371, 3.635, 3.572], [2.652, 3.572, -14.746]],
                             [[-3.098, 6.511, 7.732], [6.511, 1.419, 11.445], [7.732, 11.445, 1.678]],
                             [[-3.098, -6.511, 7.732], [-6.511, 1.419, -11.445], [7.732, -11.445, 1.678]],
                             [[2.344, -10.775, -7.006], [-10.775, -7.152, -11.309], [-7.006, -11.309, 4.808]],
                             [[2.344, 10.775, -7.006], [10.775, -7.152, 11.309], [-7.006, 11.309, 4.808]],
                             [[2.404, -0.588, -6.83], [-0.588, 10.435, 3.159], [-6.83, 3.159, -12.839]],
                             [[2.404, 0.588, -6.83], [0.588, 10.435, -3.159], [-6.83, -3.159, -12.839]]]

        self.assertEqual(len(outcar.data["unsym_efg_tensor"][2:10]), len(exepected_tensors))
        for e1, e2 in zip(outcar.data["unsym_efg_tensor"][2:10], exepected_tensors):
            self.assertArrayAlmostEqual(e1, e2)

    def test_read_fermi_contact_shift(self):
        filepath = os.path.join(test_dir, "OUTCAR_fc")
        outcar = Outcar(filepath)
        outcar.read_fermi_contact_shift()
        self.assertAlmostEqual(outcar.data["fermi_contact_shift"][u'fch'][0][0],
                               -0.002)
        self.assertAlmostEqual(outcar.data["fermi_contact_shift"][u'th'][0][0],
                               -0.052)
        self.assertAlmostEqual(outcar.data["fermi_contact_shift"][u'dh'][0][0],
                               0.0)

    def test_drift(self):
        outcar = Outcar(os.path.join(test_dir, "OUTCAR"))
        self.assertEqual(len(outcar.drift), 5)
        self.assertAlmostEqual(np.sum(outcar.drift), 0)

        outcar = Outcar(os.path.join(test_dir, "OUTCAR.CL"))
        self.assertEqual(len(outcar.drift), 79)
        self.assertAlmostEqual(np.sum(outcar.drift), 0.448010)

    def test_electrostatic_potential(self):

        outcar = Outcar(os.path.join(test_dir, "OUTCAR"))
        self.assertEqual(outcar.ngf, [54, 30, 54])
        self.assertTrue(
            np.allclose(outcar.sampling_radii, [0.9748, 0.9791, 0.7215]))
        self.assertTrue(np.allclose(outcar.electrostatic_potential,
                                    [-26.0704, -45.5046, -45.5046, -72.9539,
                                     -73.0621, -72.9539, -73.0621]))

    def test_mag_electrostatic_error(self):
        outcar = Outcar(os.path.join(test_dir, "OUTCAR.electrostaticerror.gz"))
        self.assertEqual(outcar.electrostatic_potential,
                         [-21.1667, -19.6865, -22.3983, -22.3307, -20.5213, -20.9292, -21.5063, -21.3554, -21.74,
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
                          -68.3194, -69.3671, -67.2708])


class BSVasprunTest(unittest.TestCase):
    _multiprocess_shared_ = True

    def test_get_band_structure(self):
        filepath = os.path.join(test_dir, 'vasprun_Si_bands.xml')
        vasprun = BSVasprun(filepath, parse_potcar_file=False)
        bs = vasprun.get_band_structure(kpoints_filename=os.path.join(test_dir,
                                                                      'KPOINTS_Si_bands'))
        cbm = bs.get_cbm()
        vbm = bs.get_vbm()
        self.assertEqual(cbm['kpoint_index'], [13], "wrong cbm kpoint index")
        self.assertAlmostEqual(cbm['energy'], 6.2301, "wrong cbm energy")
        self.assertEqual(cbm['band_index'], {Spin.up: [4], Spin.down: [4]},
                         "wrong cbm bands")
        self.assertEqual(vbm['kpoint_index'], [0, 63, 64])
        self.assertAlmostEqual(vbm['energy'], 5.6158, "wrong vbm energy")
        self.assertEqual(vbm['band_index'], {Spin.up: [1, 2, 3],
                                             Spin.down: [1, 2, 3]},
                         "wrong vbm bands")
        self.assertEqual(vbm['kpoint'].label, "\\Gamma", "wrong vbm label")
        self.assertEqual(cbm['kpoint'].label, None, "wrong cbm label")
        d = vasprun.as_dict()
        self.assertIn("eigenvalues", d["output"])


class OszicarTest(unittest.TestCase):

    def test_init(self):
        filepath = os.path.join(test_dir, 'OSZICAR')
        oszicar = Oszicar(filepath)
        self.assertEqual(len(oszicar.electronic_steps),
                         len(oszicar.ionic_steps))
        self.assertEqual(len(oszicar.all_energies), 60)
        self.assertAlmostEqual(oszicar.final_energy, -526.63928)


class LocpotTest(unittest.TestCase):

    def test_init(self):
        filepath = os.path.join(test_dir, 'LOCPOT')
        locpot = Locpot.from_file(filepath)
        self.assertAlmostEqual(-217.05226954,
                               sum(locpot.get_average_along_axis(0)))
        self.assertAlmostEqual(locpot.get_axis_grid(0)[-1], 2.87629, 2)
        self.assertAlmostEqual(locpot.get_axis_grid(1)[-1], 2.87629, 2)
        self.assertAlmostEqual(locpot.get_axis_grid(2)[-1], 2.87629, 2)


class ChgcarTest(PymatgenTest):
    _multiprocess_shared_ = True

    def test_init(self):
        filepath = os.path.join(test_dir, 'CHGCAR.nospin')
        chg = Chgcar.from_file(filepath)
        self.assertAlmostEqual(chg.get_integrated_diff(0, 2)[0, 1], 0)
        filepath = os.path.join(test_dir, 'CHGCAR.spin')
        chg = Chgcar.from_file(filepath)
        self.assertAlmostEqual(chg.get_integrated_diff(0, 1)[0, 1],
                               -0.0043896932237534022)
        # test sum
        chg += chg
        self.assertAlmostEqual(chg.get_integrated_diff(0, 1)[0, 1],
                               -0.0043896932237534022 * 2)

        filepath = os.path.join(test_dir, 'CHGCAR.Fe3O4')
        chg = Chgcar.from_file(filepath)
        ans = [1.56472768, 3.25985108, 3.49205728, 3.66275028, 3.8045896,
               5.10813352]
        myans = chg.get_integrated_diff(0, 3, 6)
        self.assertTrue(np.allclose(myans[:, 1], ans))

    def test_write(self):
        filepath = os.path.join(test_dir, 'CHGCAR.spin')
        chg = Chgcar.from_file(filepath)
        chg.write_file("CHGCAR_pmg")
        with open("CHGCAR_pmg") as f:
            for i, line in enumerate(f):
                if i == 22130:
                    self.assertEqual("augmentation occupancies   1  15\n", line)
                if i == 44255:
                    self.assertEqual("augmentation occupancies   1  15\n", line)
        os.remove("CHGCAR_pmg")

    def test_soc_chgcar(self):

        filepath = os.path.join(test_dir, "CHGCAR.NiO_SOC.gz")
        chg = Chgcar.from_file(filepath)
        self.assertEqual(set(chg.data.keys()),
                         {'total', 'diff_x', 'diff_y', 'diff_z', 'diff'})
        self.assertTrue(chg.is_soc)
        self.assertEqual(chg.data['diff'].shape, chg.data['diff_y'].shape)

        # check our construction of chg.data['diff'] makes sense
        # this has been checked visually too and seems reasonable
        self.assertEqual(abs(chg.data['diff'][0][0][0]),
                         np.linalg.norm([chg.data['diff_x'][0][0][0],
                                         chg.data['diff_y'][0][0][0],
                                         chg.data['diff_z'][0][0][0]]))

        # and that the net magnetization is about zero
        # note: we get ~ 0.08 here, seems a little high compared to
        # vasp output, but might be due to chgcar limitations?
        self.assertAlmostEqual(chg.net_magnetization, 0.0, places=0)

        chg.write_file("CHGCAR_pmg_soc")
        chg_from_file = Chgcar.from_file("CHGCAR_pmg_soc")
        self.assertTrue(chg_from_file.is_soc)
        os.remove("CHGCAR_pmg_soc")

    def test_hdf5(self):
        chgcar = Chgcar.from_file(os.path.join(test_dir, "CHGCAR.NiO_SOC.gz"))
        chgcar.to_hdf5("chgcar_test.hdf5")
        import h5py
        with h5py.File("chgcar_test.hdf5", "r") as f:
            self.assertArrayAlmostEqual(np.array(f["vdata"]["total"]),
                                        chgcar.data["total"])
            self.assertArrayAlmostEqual(np.array(f["vdata"]["diff"]),
                                        chgcar.data["diff"])
            self.assertArrayAlmostEqual(np.array(f["lattice"]),
                                        chgcar.structure.lattice.matrix)
            self.assertArrayAlmostEqual(np.array(f["fcoords"]),
                                        chgcar.structure.frac_coords)
            for z in f["Z"]:
                self.assertIn(z, [Element.Ni.Z, Element.O.Z])

            for sp in f["species"]:
                self.assertIn(sp, ["Ni", "O"])

        chgcar2 = Chgcar.from_hdf5("chgcar_test.hdf5")
        self.assertArrayAlmostEqual(chgcar2.data["total"],
                                    chgcar.data["total"])

        os.remove("chgcar_test.hdf5")

    def test_as_dict_and_from_dict(self):
        chgcar = Chgcar.from_file(os.path.join(test_dir, "CHGCAR.NiO_SOC.gz"))
        d = chgcar.as_dict()
        chgcar_from_dict = Chgcar.from_dict(d)
        self.assertArrayAlmostEqual(chgcar.data['total'], chgcar_from_dict.data['total'])
        self.assertArrayAlmostEqual(chgcar.structure.lattice.matrix,
                                    chgcar_from_dict.structure.lattice.matrix)


class ProcarTest(unittest.TestCase):
    _multiprocess_shared_ = True
    def test_init(self):
        filepath = os.path.join(test_dir, 'PROCAR.simple')
        p = Procar(filepath)
        self.assertAlmostEqual(p.get_occupation(0, 'd')[Spin.up], 0)
        self.assertAlmostEqual(p.get_occupation(0, 's')[Spin.up],
                               0.35381249999999997)
        self.assertAlmostEqual(p.get_occupation(0, 'p')[Spin.up], 1.19540625)
        self.assertRaises(ValueError, p.get_occupation, 1, 'm')
        self.assertEqual(p.nbands, 10)
        self.assertEqual(p.nkpoints, 10)
        self.assertEqual(p.nions, 3)
        lat = Lattice.cubic(3.)
        s = Structure(lat, ["Li", "Na", "K"], [[0., 0., 0.],
                                               [0.25, 0.25, 0.25],
                                               [0.75, 0.75, 0.75]])
        d = p.get_projection_on_elements(s)
        self.assertAlmostEqual(d[Spin.up][2][2],
                               {'Na': 0.042, 'K': 0.646, 'Li': 0.042})
        filepath = os.path.join(test_dir, 'PROCAR')
        p = Procar(filepath)
        self.assertAlmostEqual(p.get_occupation(0, 'dxy')[Spin.up],
                               0.96214813853000025)
        self.assertAlmostEqual(p.get_occupation(0, 'dxy')[Spin.down],
                               0.85796295426000124)

    def test_phase_factors(self):
        filepath = os.path.join(test_dir, 'PROCAR.phase')
        p = Procar(filepath)
        self.assertAlmostEqual(p.phase_factors[Spin.up][0, 0, 0, 0],
                               -0.746 + 0.099j)
        self.assertAlmostEqual(p.phase_factors[Spin.down][0, 0, 0, 0],
                               0.372 - 0.654j)

        # Two Li should have same phase factor.
        self.assertAlmostEqual(p.phase_factors[Spin.up][0, 0, 0, 0],
                               p.phase_factors[Spin.up][0, 0, 1, 0])
        self.assertAlmostEqual(p.phase_factors[Spin.up][0, 0, 2, 0],
                               -0.053 + 0.007j)
        self.assertAlmostEqual(p.phase_factors[Spin.down][0, 0, 2, 0],
                               0.027 - 0.047j)


class XdatcarTest(unittest.TestCase):

    def test_init(self):
        filepath = os.path.join(test_dir, 'XDATCAR_4')
        x = Xdatcar(filepath)
        structures = x.structures
        self.assertEqual(len(structures), 4)
        for s in structures:
            self.assertEqual(s.formula, "Li2 O1")

        filepath = os.path.join(test_dir, 'XDATCAR_5')
        x = Xdatcar(filepath)
        structures = x.structures
        self.assertEqual(len(structures), 4)
        for s in structures:
            self.assertEqual(s.formula, "Li2 O1")

        x.concatenate(os.path.join(test_dir, 'XDATCAR_4'))
        self.assertEqual(len(x.structures), 8)
        self.assertIsNotNone(x.get_string())


class DynmatTest(unittest.TestCase):

    def test_init(self):
        # nosetests pymatgen/io/vasp/tests/test_outputs.py:DynmatTest.test_init
        filepath = os.path.join(test_dir, 'DYNMAT')
        d = Dynmat(filepath)
        self.assertEqual(d.nspecs, 2)
        self.assertEqual(d.natoms, 6)
        self.assertEqual(d.ndisps, 3)
        self.assertTrue(np.allclose(d.masses, [63.546, 196.966]))
        self.assertTrue(4 in d.data)
        self.assertTrue(2 in d.data[4])
        self.assertTrue(np.allclose(
            d.data[4][2]['dispvec'], [0., 0.05, 0.]
        ))
        self.assertTrue(np.allclose(
            d.data[4][2]['dynmat'][3], [0.055046, -0.298080, 0.]
        ))
        # TODO: test get_phonon_frequencies once cross-checked


class WavecarTest(unittest.TestCase):
    _multiprocess_shared_ = True

    def setUp(self):
        a = np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0],
                          [0.0, 0.0, 10.0]])
        self.vol = np.dot(a[0, :], np.cross(a[1, :], a[2, :]))
        b = np.array([np.cross(a[1, :], a[2, :]),
                      np.cross(a[2, :], a[0, :]),
                      np.cross(a[0, :], a[1, :])])
        self.b = 2 * np.pi * b / self.vol
        self.a = a
        self.w = Wavecar(os.path.join(test_dir, 'WAVECAR.N2'))

    def test_standard(self):
        w = self.w
        a = np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0],
                     [0.0, 0.0, 10.0]])
        vol = np.dot(a[0, :], np.cross(a[1, :], a[2, :]))
        b = np.array([np.cross(a[1, :], a[2, :]),
                           np.cross(a[2, :], a[0, :]),
                           np.cross(a[0, :], a[1, :])])
        b = 2 * np.pi * b / vol

        self.assertEqual(w.filename, os.path.join(test_dir, 'WAVECAR.N2'))
        self.assertAlmostEqual(w.efermi, -5.7232, places=4)
        self.assertEqual(w.encut, 25)
        self.assertEqual(w.nb, 9)
        self.assertEqual(w.nk, 1)
        self.assertTrue(np.allclose(w.a, a))
        self.assertTrue(np.allclose(w.b, b))
        self.assertAlmostEqual(w.vol, vol)
        self.assertEqual(len(w.kpoints), w.nk)
        self.assertEqual(len(w.coeffs), w.nk)
        self.assertEqual(len(w.coeffs[0]), w.nb)
        self.assertEqual(len(w.band_energy), w.nk)
        self.assertEqual(w.band_energy[0].shape, (w.nb, 3))
        self.assertLessEqual(len(w.Gpoints[0]), 257)
        for k in range(w.nk):
            for b in range(w.nb):
                self.assertEqual(len(w.coeffs[k][b]),
                                 len(w.Gpoints[k]))

        with self.assertRaises(ValueError):
            Wavecar(os.path.join(test_dir, 'WAVECAR.N2.malformed'))

        import sys
        from io import StringIO
        saved_stdout = sys.stdout
        try:
            out = StringIO()
            sys.stdout = out
            Wavecar(os.path.join(test_dir, 'WAVECAR.N2'), verbose=True)
            self.assertNotEqual(out.getvalue().strip(), '')
        finally:
            sys.stdout = saved_stdout

    def test_n2_45210(self):
        w = Wavecar(os.path.join(test_dir, 'WAVECAR.N2.45210'))
        self.assertEqual(w.filename, os.path.join(test_dir,
                                                       'WAVECAR.N2.45210'))
        self.assertAlmostEqual(w.efermi, -5.7232, places=4)
        self.assertEqual(w.encut, 25)
        self.assertEqual(w.nb, 9)
        self.assertEqual(w.nk, 1)
        self.assertTrue(np.allclose(w.a, self.a))
        self.assertTrue(np.allclose(w.b, self.b))
        self.assertAlmostEqual(w.vol, self.vol)
        self.assertEqual(len(w.kpoints), w.nk)
        self.assertEqual(len(w.coeffs), w.nk)
        self.assertEqual(len(w.coeffs[0]), w.nb)
        self.assertEqual(len(w.band_energy), w.nk)
        self.assertEqual(w.band_energy[0].shape, (w.nb, 3))
        self.assertLessEqual(len(w.Gpoints[0]), 257)

    def test_n2_spin(self):
        w = Wavecar(os.path.join(test_dir, 'WAVECAR.N2.spin'))
        self.assertEqual(len(w.coeffs), 2)
        self.assertEqual(len(w.band_energy), 2)
        self.assertEqual(len(w.kpoints), w.nk)
        self.assertEqual(len(w.Gpoints), w.nk)
        self.assertEqual(len(w.coeffs[0][0]), w.nb)
        self.assertEqual(len(w.band_energy[0]), w.nk)

        temp_ggp = Wavecar._generate_G_points
        try:
            Wavecar._generate_G_points = lambda x, y: []
            with self.assertRaises(ValueError):
                Wavecar(os.path.join(test_dir, 'WAVECAR.N2'))
        finally:
            Wavecar._generate_G_points = temp_ggp

    def test__generate_nbmax(self):
        self.w._generate_nbmax()
        self.assertEqual(self.w._nbmax.tolist(), [5, 5, 5])

    def test__generate_G_points(self):
        for k in range(self.w.nk):
            kp = self.w.kpoints[k]
            self.assertLessEqual(len(self.w._generate_G_points(kp)), 257)

    def test_evaluate_wavefunc(self):
        self.w.Gpoints.append(np.array([0, 0, 0]))
        self.w.kpoints.append(np.array([0, 0, 0]))
        self.w.coeffs.append([[1 + 1j]])
        self.assertAlmostEqual(self.w.evaluate_wavefunc(-1, -1, [0, 0, 0]),
                               (1 + 1j) / np.sqrt(self.vol), places=4)
        self.assertAlmostEqual(self.w.evaluate_wavefunc(0, 0, [0, 0, 0]),
                               np.sum(self.w.coeffs[0][0]) / np.sqrt(self.vol),
                               places=4)

    def test_fft_mesh(self):
        mesh = self.w.fft_mesh(0, 5)
        ind = np.argmax(np.abs(mesh))
        self.assertEqual(np.unravel_index(ind, mesh.shape), (14, 1, 1))
        self.assertEqual(mesh[tuple((self.w.ng / 2).astype(np.int))], 0j)
        mesh = self.w.fft_mesh(0, 5, shift=False)
        ind = np.argmax(np.abs(mesh))
        self.assertEqual(np.unravel_index(ind, mesh.shape), (6, 8, 8))
        self.assertEqual(mesh[0, 0, 0], 0j)


if __name__ == "__main__":
    unittest.main()
