# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import os
import unittest
import numpy as np

from pymatgen.util.testing import PymatgenTest
from pymatgen.io.vasp import Vasprun
from pymatgen.analysis.defects.core import DefectEntry, Vacancy
from pymatgen.analysis.defects.corrections import FreysoldtCorrection,\
            BandFillingCorrection, BandEdgeShiftingCorrection

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", 'test_files')


class DefectsCorrectionsTest(PymatgenTest):
    def test_freysoldt(self):
        struc = PymatgenTest.get_structure("VO2")
        struc.make_supercell(3)
        struc = struc
        vac = Vacancy(struc, struc.sites[0], charge=-3)

        abc = struc.lattice.abc
        axisdata = [np.arange(0., lattval, 0.2) for lattval in abc]
        bldata = [np.array([1. for u in np.arange(0., lattval, 0.2)]) for lattval in abc]
        dldata = [
            np.array([(-1 - np.cos(2 * np.pi * u / lattval)) for u in np.arange(0., lattval, 0.2)]) for lattval in abc
        ]
        params = {'axis_grid': axisdata, 'bulk_planar_averages': bldata, 'defect_planar_averages': dldata}
        fc = FreysoldtCorrection(15)

        #test electrostatic correction
        es_corr = fc.perform_es_corr(struc.lattice, -3)
        self.assertAlmostEqual(es_corr, 0.975893)

        #test potential alignment method
        pot_corr = fc.perform_pot_corr(axisdata[0], bldata[0], dldata[0], struc.lattice, -3, vac.site.coords, 0)
        self.assertAlmostEqual(pot_corr, 2.836369987722345)

        #test entry full correction method
        de = DefectEntry(vac, 0., corrections={}, parameters=params, entry_id=None)
        val = fc.get_correction(de)
        self.assertAlmostEqual(val['freysoldt_electrostatic'], 0.975893)
        self.assertAlmostEqual(val['freysoldt_potential_alignment'], 4.4700574)

        #test the freysoldt plotter and that plot metadata exists
        pltsaver = []
        for ax in range(3):
            pltsaver.append(fc.plot(axis=ax))
        self.assertAlmostEqual(len(pltsaver), 3)

        #check that uncertainty metadata exists
        for ax in range(3):
            self.assertAlmostEqual(set(fc.metadata['pot_corr_uncertainty_md'][ax].keys()), set(['potcorr', 'stats']))

        #test a specified axis from entry
        fc = FreysoldtCorrection(15, axis=[1])
        val = fc.get_correction(de)
        self.assertAlmostEqual(val['freysoldt_potential_alignment'], 5.2869010593283132)

        #test a different charge
        #   for electrostatic correction
        es_corr = fc.perform_es_corr(struc.lattice, 2)
        self.assertAlmostEqual(es_corr, 0.43373)
        #   for potential alignment method
        pot_corr = fc.perform_pot_corr(axisdata[0], bldata[0], dldata[0], struc.lattice, 2, vac.site.coords, 0)
        self.assertAlmostEqual(pot_corr, -2.1375685936497768)

        #test an input anisotropic dielectric constant
        fc = FreysoldtCorrection([[1., 2., 3.], [0., 3., 5.], [4., 10., 8.]])
        self.assertAlmostEqual(fc.dielectric, 4.)
        val = fc.get_correction(de)
        self.assertAlmostEqual(val['freysoldt_electrostatic'], 3.659599)
        self.assertAlmostEqual(val['freysoldt_potential_alignment'], 3.3605255195745087)

        #test potalign being added to defect entry
        self.assertAlmostEqual(de.parameters['potalign'], 1.1201751731915028)

        #test that metadata entries exist in defect entry
        self.assertTrue('freysoldt_meta' in de.parameters.keys())
        self.assertAlmostEqual(
            set(de.parameters['freysoldt_meta'].keys()), set(['pot_plot_data', 'pot_corr_uncertainty_md']))

        #test a charge of zero
        vac = Vacancy(struc, struc.sites[0], charge=0)
        de = DefectEntry(vac, 0., corrections={}, parameters=params, entry_id=None)
        val = fc.get_correction(de)
        self.assertAlmostEqual(val['freysoldt_electrostatic'], 0.)
        self.assertAlmostEqual(val['freysoldt_potential_alignment'], 0.)

    def test_bandfilling(self):
        v = Vasprun(os.path.join(test_dir, 'vasprun.xml'))
        eigenvalues = v.eigenvalues.copy()
        kptweights = v.actual_kpoints_weights
        potalign = 0.
        vbm = v.eigenvalue_band_properties[2]
        cbm = v.eigenvalue_band_properties[1]
        params = {
            'eigenvalues': eigenvalues,
            'kpoint_weights': kptweights,
            'potalign': potalign,
            'vbm': vbm,
            'cbm': cbm
        }
        bfc = BandFillingCorrection()
        struc = PymatgenTest.get_structure("VO2")
        struc.make_supercell(3)
        vac = Vacancy(struc, struc.sites[0], charge=-3)

        #test trivial performing bandfilling correction
        bf_corr = bfc.perform_bandfill_corr(eigenvalues, kptweights, potalign, vbm, cbm)
        self.assertAlmostEqual(bf_corr, 0.)
        self.assertFalse(bfc.metadata['occupied_def_levels'])
        self.assertFalse(bfc.metadata['unoccupied_def_levels'])
        self.assertFalse(bfc.metadata['total_occupation_defect_levels'])
        self.assertFalse(bfc.metadata['num_elec_cbm'])
        self.assertFalse(bfc.metadata['num_hole_vbm'])
        self.assertFalse(bfc.metadata['potalign'])

        #test trivial full entry bandfill evaluation
        de = DefectEntry(vac, 0., corrections={}, parameters=params, entry_id=None)

        corr = bfc.get_correction(de)
        self.assertAlmostEqual(corr['bandfilling'], 0.)

        #modify the eigenvalue list to have free holes
        hole_eigenvalues = {}
        for spinkey, spinset in eigenvalues.items():
            hole_eigenvalues[spinkey] = []
            for kptset in spinset:
                hole_eigenvalues[spinkey].append([])
                for eig in kptset:
                    if (eig[0] < vbm) and (eig[0] > vbm - .8):
                        hole_eigenvalues[spinkey][-1].append([eig[0], 0.5])
                    else:
                        hole_eigenvalues[spinkey][-1].append(eig)

        hole_bf_corr = bfc.perform_bandfill_corr(hole_eigenvalues, kptweights, potalign, vbm, cbm)
        self.assertAlmostEqual(hole_bf_corr, -0.41138336)
        self.assertAlmostEqual(bfc.metadata['num_hole_vbm'], 0.8125000649)
        self.assertFalse(bfc.metadata['num_elec_cbm'])

        #modify the eigenvalue list to have free electrons
        elec_eigenvalues = {}
        for spinkey, spinset in eigenvalues.items():
            elec_eigenvalues[spinkey] = []
            for kptset in spinset:
                elec_eigenvalues[spinkey].append([])
                for eig in kptset:
                    if (eig[0] > cbm) and (eig[0] < cbm + .2):
                        elec_eigenvalues[spinkey][-1].append([eig[0], 0.5])
                    else:
                        elec_eigenvalues[spinkey][-1].append(eig)

        elec_bf_corr = bfc.perform_bandfill_corr(elec_eigenvalues, kptweights, potalign, vbm, cbm)
        self.assertAlmostEqual(elec_bf_corr, -0.0903187572254)
        self.assertAlmostEqual(bfc.metadata['num_elec_cbm'], 0.8541667349)
        self.assertFalse(bfc.metadata['num_hole_vbm'])

        #modify the potalignment and introduce new occupied defect levels from vbm states
        potalign = -0.1

        bf_corr = bfc.perform_bandfill_corr(eigenvalues, kptweights, potalign, vbm, cbm)
        self.assertAlmostEqual(bfc.metadata['num_hole_vbm'], 0.)
        self.assertAlmostEqual(bf_corr, 0.)
        occu = [[1.457,  0.0833333], [1.5204, 0.0833333], [1.53465, 0.0833333], [1.5498, 0.0416667]]
        self.assertArrayAlmostEqual(
            list(sorted(bfc.metadata['occupied_def_levels'], key=lambda x: x[0])), list(
                sorted(occu, key=lambda x: x[0])))
        self.assertAlmostEqual(bfc.metadata['total_occupation_defect_levels'], 0.29166669)
        self.assertFalse(bfc.metadata['unoccupied_def_levels'])

    def test_bandedgeshifting(self):
        struc = PymatgenTest.get_structure("VO2")
        struc.make_supercell(3)
        struc = struc
        vac = Vacancy(struc, struc.sites[0], charge=-3)

        besc = BandEdgeShiftingCorrection()
        params = {'hybrid_cbm': 1., 'hybrid_vbm': -1., 'vbm': -0.5, 'cbm': 0.6, 'num_hole_vbm': 0., 'num_elec_cbm': 0.}
        de = DefectEntry(vac, 0., corrections={}, parameters=params, entry_id=None)

        #test with no free carriers
        corr = besc.get_correction(de)
        self.assertEqual(corr['vbm_shift_correction'], 1.5)
        self.assertEqual(corr['elec_cbm_shift_correction'], 0.)
        self.assertEqual(corr['hole_vbm_shift_correction'], 0.)

        #test with free holes
        de.parameters.update({'num_hole_vbm': 1.})
        corr = besc.get_correction(de)
        self.assertEqual(corr['vbm_shift_correction'], 1.5)
        self.assertEqual(corr['elec_cbm_shift_correction'], 0.)
        self.assertEqual(corr['hole_vbm_shift_correction'], 0.5)

        #test with free electrons
        de.parameters.update({'num_hole_vbm': 0., 'num_elec_cbm': 1.})
        corr = besc.get_correction(de)
        self.assertEqual(corr['vbm_shift_correction'], 1.5)
        self.assertEqual(corr['elec_cbm_shift_correction'], 0.4)
        self.assertEqual(corr['hole_vbm_shift_correction'], 0.)


if __name__ == "__main__":
    unittest.main()
