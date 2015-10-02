# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
Created on Jul 16, 2012
"""


__author__ = "Shyue Ping Ong, Stephen Dacek"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Jul 16, 2012"

import unittest
import os
import json
import numpy as np
import warnings

from pymatgen.io.vasp.outputs import Chgcar, Locpot, Oszicar, Outcar, \
    Vasprun, Procar, Xdatcar, Dynmat
from pymatgen import Spin, Orbital, Lattice, Structure
from pymatgen.entries.compatibility import MaterialsProjectCompatibility

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')


class VasprunTest(unittest.TestCase):

    def test_properties(self):

        filepath = os.path.join(test_dir, 'vasprun.xml.nonlm')
        vasprun = Vasprun(filepath, parse_potcar_file=False)
        orbs = list(vasprun.complete_dos.pdos[vasprun.final_structure[
            0]].keys())
        self.assertIn("S", orbs)
        filepath = os.path.join(test_dir, 'vasprun.xml')
        vasprun = Vasprun(filepath, parse_potcar_file=False)

        #Test NELM parsing.
        self.assertEqual(vasprun.parameters["NELM"], 60)
        #test pdos parsing

        pdos0 = vasprun.complete_dos.pdos[vasprun.final_structure[0]]
        self.assertAlmostEqual(pdos0[Orbital.s][1][16], 0.0026)
        self.assertAlmostEqual(pdos0[Orbital.pz][-1][16], 0.0012)
        self.assertEqual(pdos0[Orbital.s][1].shape, (301, ))


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

        #test skipping ionic steps.
        vasprun_skip = Vasprun(filepath, 3, parse_potcar_file=False)
        self.assertEqual(vasprun_skip.nionic_steps, 29)
        self.assertEqual(len(vasprun_skip.ionic_steps),
                         int(vasprun.nionic_steps / 3) + 1)
        self.assertEqual(len(vasprun_skip.ionic_steps),
                         len(vasprun_skip.structures))
        self.assertEqual(len(vasprun_skip.ionic_steps),
                         int(vasprun.nionic_steps / 3) + 1)
        #Check that nionic_steps is preserved no matter what.
        self.assertEqual(vasprun_skip.nionic_steps,
                         vasprun.nionic_steps)

        self.assertNotAlmostEqual(vasprun_skip.final_energy,
                                  vasprun.final_energy)

        #Test with ionic_step_offset
        vasprun_offset = Vasprun(filepath, 3, 6, parse_potcar_file=False)
        self.assertEqual(len(vasprun_offset.ionic_steps),
                         int(len(vasprun.ionic_steps) / 3) - 1)
        self.assertEqual(vasprun_offset.structures[0],
                         vasprun_skip.structures[2])

        self.assertTrue(vasprun_ggau.is_hubbard)
        self.assertEqual(vasprun_ggau.hubbards["Fe"], 4.3)
        self.assertAlmostEqual(vasprun_ggau.projected_eigenvalues[(Spin.up, 0,
                                                                   0, 96,
                                                                   Orbital.s)],
                               0.0032)
        d = vasprun_ggau.as_dict()
        self.assertEqual(d["elements"], ["Fe", "Li", "O", "P"])
        self.assertEqual(d["nelements"], 4)

        filepath = os.path.join(test_dir, 'vasprun.xml.unconverged')
        vasprun_unconverged = Vasprun(filepath, parse_potcar_file=False)
        self.assertTrue(vasprun_unconverged.converged_ionic)
        self.assertFalse(vasprun_unconverged.converged_electronic)
        self.assertFalse(vasprun_unconverged.converged)

        filepath = os.path.join(test_dir, 'vasprun.xml.dfpt')
        vasprun_dfpt = Vasprun(filepath, parse_potcar_file=False)
        self.assertAlmostEqual(vasprun_dfpt.epsilon_static[0][0], 3.26105533)
        self.assertAlmostEqual(vasprun_dfpt.epsilon_static[0][1], -0.00459066)
        self.assertAlmostEqual(vasprun_dfpt.epsilon_static[2][2], 3.24330517)
        self.assertAlmostEqual(vasprun_dfpt.epsilon_static_wolfe[0][0], 3.33402531)
        self.assertAlmostEqual(vasprun_dfpt.epsilon_static_wolfe[0][1], -0.00559998)
        self.assertAlmostEqual(vasprun_dfpt.epsilon_static_wolfe[2][2], 3.31237357)
        self.assertTrue(vasprun_dfpt.converged)

        entry = vasprun_dfpt.get_computed_entry()
        entry = MaterialsProjectCompatibility(check_potcar_hash=False).process_entry(entry)
        self.assertAlmostEqual(entry.uncorrected_energy + entry.correction,
                               entry.energy)


        filepath = os.path.join(test_dir, 'vasprun.xml.dfpt.ionic')
        vasprun_dfpt_ionic = Vasprun(filepath, parse_potcar_file=False)
        self.assertAlmostEqual(vasprun_dfpt_ionic.epsilon_ionic[0][0], 515.73485838)
        self.assertAlmostEqual(vasprun_dfpt_ionic.epsilon_ionic[0][1], -0.00263523)
        self.assertAlmostEqual(vasprun_dfpt_ionic.epsilon_ionic[2][2], 19.02110169)


        filepath = os.path.join(test_dir, 'vasprun.xml.dfpt.unconverged')
        vasprun_dfpt_unconv = Vasprun(filepath, parse_potcar_file=False)
        self.assertFalse(vasprun_dfpt_unconv.converged_electronic)
        self.assertTrue(vasprun_dfpt_unconv.converged_ionic)
        self.assertFalse(vasprun_dfpt_unconv.converged)

        vasprun_uniform = Vasprun(os.path.join(test_dir, "vasprun.xml.uniform"),
                                  parse_potcar_file=False)
        self.assertEqual(vasprun_uniform.kpoints.style, "Reciprocal")


        vasprun_no_pdos = Vasprun(os.path.join(test_dir, "Li_no_projected.xml"),
                                  parse_potcar_file=False)
        self.assertIsNotNone(vasprun_no_pdos.complete_dos)
        self.assertFalse(vasprun_no_pdos.dos_has_errors)

        vasprun_diel = Vasprun(os.path.join(test_dir, "vasprun.xml.dielectric"),
                               parse_potcar_file=False)
        self.assertAlmostEqual(0.4294,vasprun_diel.dielectric[0][10])
        self.assertAlmostEqual(19.941,vasprun_diel.dielectric[1][51][0])
        self.assertAlmostEqual(19.941,vasprun_diel.dielectric[1][51][1])
        self.assertAlmostEqual(19.941,vasprun_diel.dielectric[1][51][2])
        self.assertAlmostEqual(0.0,vasprun_diel.dielectric[1][51][3])
        self.assertAlmostEqual(34.186,vasprun_diel.dielectric[2][85][0])
        self.assertAlmostEqual(34.186,vasprun_diel.dielectric[2][85][1])
        self.assertAlmostEqual(34.186,vasprun_diel.dielectric[2][85][2])
        self.assertAlmostEqual(0.0,vasprun_diel.dielectric[2][85][3])

    def test_Xe(self):
        vr = Vasprun(os.path.join(test_dir, 'vasprun.xml.xe'), parse_potcar_file=False)
        self.assertEquals(vr.atomic_symbols, ['Xe'])

    def test_invalid_element(self):
        self.assertRaises(KeyError, Vasprun, os.path.join(test_dir, 'vasprun.xml.wrong_sp'))

    def test_as_dict(self):
        filepath = os.path.join(test_dir, 'vasprun.xml')
        vasprun = Vasprun(filepath,
                          parse_potcar_file=False)
        #Test that as_dict() is json-serializable
        self.assertIsNotNone(json.dumps(vasprun.as_dict()))
        self.assertEqual(
            vasprun.as_dict()["input"]["potcar_type"],
            ['PAW_PBE', 'PAW_PBE', 'PAW_PBE', 'PAW_PBE', 'PAW_PBE'])

    def test_get_band_structure(self):
        filepath = os.path.join(test_dir, 'vasprun_Si_bands.xml')
        vasprun = Vasprun(filepath, parse_potcar_file=False)
        bs = vasprun.get_band_structure(kpoints_filename=
                                        os.path.join(test_dir,
                                                     'KPOINTS_Si_bands'))
        cbm = bs.get_cbm()
        vbm = bs.get_vbm()
        self.assertEqual(cbm['kpoint_index'], [13], "wrong cbm kpoint index")
        self.assertAlmostEqual(cbm['energy'], 6.2301, "wrong cbm energy")
        self.assertEqual(cbm['band_index'], {Spin.up: [4], Spin.down: [4]},
                         "wrong cbm bands")
        self.assertEqual(vbm['kpoint_index'], [0, 63, 64],
                         "wrong vbm kpoint index")
        self.assertAlmostEqual(vbm['energy'], 5.6158, "wrong vbm energy")
        self.assertEqual(vbm['band_index'], {Spin.up: [1, 2, 3],
                                             Spin.down: [1, 2, 3]},
                         "wrong vbm bands")
        self.assertEqual(vbm['kpoint'].label, "\Gamma", "wrong vbm label")
        self.assertEqual(cbm['kpoint'].label, None, "wrong cbm label")

    def test_sc_step_overflow(self):
        filepath = os.path.join(test_dir, 'vasprun.xml.sc_overflow')
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            vasprun = Vasprun(filepath)
            self.assertEqual(len(w), 3)
        estep = vasprun.ionic_steps[0]['electronic_steps'][29]
        self.assertTrue(np.isnan(estep['e_wo_entrp']))

    def test_update_potcar(self):
        filepath = os.path.join(test_dir, 'vasprun.xml')
        potcar_path = os.path.join(test_dir, 'POTCAR.LiFePO4.gz')
        potcar_path2 = os.path.join(test_dir, 'POTCAR2.LiFePO4.gz')
        vasprun = Vasprun(filepath, parse_potcar_file=False)
        self.assertEqual(vasprun.potcar_spec, [{"titel": "PAW_PBE Li 17Jan2003", "hash": None},
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

        self.assertRaises(ValueError, Vasprun, filepath, parse_potcar_file=potcar_path2)

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
        #Ensure no potcar is found and nothing is updated
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            vasprun = Vasprun(filepath, parse_potcar_file='.')
            self.assertEqual(len(w), 1)
        self.assertEqual(vasprun.potcar_spec, [{"titel": "PAW_PBE Li 17Jan2003", "hash": None},
                                               {"titel": "PAW_PBE Fe 06Sep2000", "hash": None},
                                               {"titel": "PAW_PBE Fe 06Sep2000", "hash": None},
                                               {"titel": "PAW_PBE P 17Jan2003", "hash": None},
                                               {"titel": "PAW_PBE O 08Apr2002", "hash": None}])


class OutcarTest(unittest.TestCase):

    def test_init(self):
        for f in ['OUTCAR', 'OUTCAR.gz']:
            filepath = os.path.join(test_dir, f)
            outcar = Outcar(filepath)
            expected_mag = ({'d': 0.0, 'p': 0.003, 's': 0.002, 'tot': 0.005},
                             {'d': 0.798, 'p': 0.008, 's': 0.007, 'tot': 0.813},
                             {'d': 0.798, 'p': 0.008, 's': 0.007, 'tot': 0.813},
                             {'d': 0.0, 'p':-0.117, 's': 0.005, 'tot':-0.112},
                             {'d': 0.0, 'p':-0.165, 's': 0.004, 'tot':-0.162},
                             {'d': 0.0, 'p':-0.117, 's': 0.005, 'tot':-0.112},
                             {'d': 0.0, 'p':-0.165, 's': 0.004, 'tot':-0.162})
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

        filepath = os.path.join(test_dir, 'OUTCAR.stopped')
        outcar = Outcar(filepath)
        self.assertTrue(outcar.is_stopped)

        for f in ['OUTCAR.lepsilon', 'OUTCAR.lepsilon.gz']:
            filepath = os.path.join(test_dir, f)
            outcar = Outcar(filepath)

            outcar.read_lepsilon()
            outcar.read_lepsilon_ionic()

            self.assertAlmostEqual(outcar.dielectric_tensor[0][0], 3.716432)
            self.assertAlmostEqual(outcar.dielectric_tensor[0][1], -0.20464)
            self.assertAlmostEqual(outcar.dielectric_tensor[1][2], -0.20464)
            self.assertAlmostEqual(outcar.dielectric_ionic_tensor[0][0], 0.001419)
            self.assertAlmostEqual(outcar.dielectric_ionic_tensor[0][2], 0.001419)
            self.assertAlmostEqual(outcar.dielectric_ionic_tensor[2][2], 0.001419)
            self.assertAlmostEqual(outcar.piezo_tensor[0][0], 0.52799)
            self.assertAlmostEqual(outcar.piezo_tensor[1][3], 0.35998)
            self.assertAlmostEqual(outcar.piezo_tensor[2][5], 0.35997)
            self.assertAlmostEqual(outcar.piezo_ionic_tensor[0][0], 0.05868)
            self.assertAlmostEqual(outcar.piezo_ionic_tensor[1][3], 0.06241)
            self.assertAlmostEqual(outcar.piezo_ionic_tensor[2][5], 0.06242)
            self.assertAlmostEqual(outcar.born[0][1][2], -0.385)
            self.assertAlmostEqual(outcar.born[1][2][0], 0.36465)

    def test_elastic_tensor(self):
        filepath = os.path.join(test_dir, "OUTCAR.total_tensor.Li2O.gz")
        outcar = Outcar(filepath)

        elastic_tensor = outcar.elastic_tensor

        self.assertAlmostEqual(elastic_tensor[0][0], 1986.3391)
        self.assertAlmostEqual(elastic_tensor[0][1], 187.8324)
        self.assertAlmostEqual(elastic_tensor[3][3], 586.3034)

    def test_core_state_eigen(self):
        filepath = os.path.join(test_dir, "OUTCAR.CL")
        cl = Outcar(filepath).read_core_state_eigen()
        self.assertAlmostEqual(cl[6]["2s"][-1], -174.4779)

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


class ChgcarTest(unittest.TestCase):

    def test_init(self):
        filepath = os.path.join(test_dir, 'CHGCAR.nospin')
        chg = Chgcar.from_file(filepath)
        self.assertAlmostEqual(chg.get_integrated_diff(0, 2)[0, 1], 0)
        filepath = os.path.join(test_dir, 'CHGCAR.spin')
        chg = Chgcar.from_file(filepath)
        self.assertAlmostEqual(chg.get_integrated_diff(0, 1)[0, 1],
                               -0.0043896932237534022)
        #test sum
        chg += chg
        self.assertAlmostEqual(chg.get_integrated_diff(0, 1)[0, 1],
                               -0.0043896932237534022 * 2)

        filepath = os.path.join(test_dir, 'CHGCAR.Fe3O4')
        chg = Chgcar.from_file(filepath)
        ans = [1.93313368, 3.91201473, 4.11858277, 4.1240093, 4.10634989,
               3.38864822]
        myans = chg.get_integrated_diff(0, 3, 6)
        self.assertTrue(np.allclose(myans[:, 1], ans))


class ProcarTest(unittest.TestCase):

    def test_init(self):
        filepath = os.path.join(test_dir, 'PROCAR.simple')
        p = Procar(filepath)
        self.assertAlmostEqual(p.get_occupation(1, 'd'), 0)
        self.assertAlmostEqual(p.get_occupation(1, 's'), 0.3538125)
        self.assertAlmostEqual(p.get_occupation(1, 'p'), 1.19540625)
        self.assertRaises(ValueError, p.get_occupation, 1, 'm')
        self.assertEqual(p.nb_bands, 10)
        self.assertEqual(p.nb_kpoints, 10)
        lat = Lattice.cubic(3.)
        s = Structure(lat, ["Li", "Na", "K"], [[0., 0., 0.],
                                               [0.25, 0.25, 0.25],
                                               [0.75, 0.75, 0.75]])
        d = p.get_projection_on_elements(s)
        self.assertAlmostEqual(d[1][2][2], {'Na': 0.042, 'K': 0.646, 'Li': 0.042})
        filepath = os.path.join(test_dir, 'PROCAR')
        p = Procar(filepath)
        self.assertAlmostEqual(p.get_occupation(0, 'd'), 4.3698147704200059)
        self.assertAlmostEqual(p.get_occupation(0, 'dxy'), 0.85796295426000124)


class XdatcarTest(unittest.TestCase):

    def test_init(self):
        filepath = os.path.join(test_dir, 'XDATCAR_4')
        x = Xdatcar(filepath)
        structures = x.structures
        self.assertEqual(len(structures), 3)
        for s in structures:
            self.assertEqual(s.formula, "Li2 O1")

        filepath = os.path.join(test_dir, 'XDATCAR_5')
        x = Xdatcar(filepath)
        structures = x.structures
        self.assertEqual(len(structures), 3)
        for s in structures:
            self.assertEqual(s.formula, "Li2 O1")

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

if __name__ == "__main__":
    unittest.main()
