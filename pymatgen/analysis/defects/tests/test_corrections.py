# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import os
import unittest
import numpy as np

from pymatgen.core import Lattice
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.vasp import Vasprun, Poscar, Outcar
from pymatgen.analysis.defects.core import DefectEntry, Vacancy
from pymatgen.analysis.defects.corrections import FreysoldtCorrection, \
    BandFillingCorrection, BandEdgeShiftingCorrection, KumagaiCorrection
from pymatgen.analysis.defects.utils import generate_R_and_G_vecs

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", 'test_files')


class DefectsCorrectionsTest(PymatgenTest):
    def test_freysoldt(self):
        struc = PymatgenTest.get_structure("VO2")
        struc.make_supercell(3)
        struc = struc
        vac = Vacancy(struc, struc.sites[0], charge=-3)
        ids = vac.generate_defect_structure(1)

        abc = struc.lattice.abc
        axisdata = [np.arange(0., lattval, 0.2) for lattval in abc]
        bldata = [np.array([1. for u in np.arange(0., lattval, 0.2)]) for lattval in abc]
        dldata = [
            np.array([(-1 - np.cos(2 * np.pi * u / lattval)) for u in np.arange(0., lattval, 0.2)]) for lattval in abc
        ]
        params = {'axis_grid': axisdata, 'bulk_planar_averages': bldata, 'defect_planar_averages': dldata,
                  'initial_defect_structure': ids, 'defect_frac_sc_coords': struc.sites[0].frac_coords}
        fc = FreysoldtCorrection(15)

        # test electrostatic correction
        es_corr = fc.perform_es_corr(struc.lattice, -3)
        self.assertAlmostEqual(es_corr, 0.975893)

        # test potential alignment method
        pot_corr = fc.perform_pot_corr(axisdata[0], bldata[0], dldata[0], struc.lattice, -3, vac.site.coords, 0)
        self.assertAlmostEqual(pot_corr, 2.836369987722345)

        # test entry full correction method
        de = DefectEntry(vac, 0., corrections={}, parameters=params, entry_id=None)
        val = fc.get_correction(de)
        self.assertAlmostEqual(val['freysoldt_electrostatic'], 0.975893)
        self.assertAlmostEqual(val['freysoldt_potential_alignment'], 4.4700574)

        # test the freysoldt plotter
        for ax in range(3):
            fcp = fc.plot(axis=ax)
            self.assertTrue(fcp)

        # check that uncertainty metadata exists
        for ax in range(3):
            self.assertAlmostEqual(set(fc.metadata['pot_corr_uncertainty_md'][ax].keys()), set(['potcorr', 'stats']))

        # test a specified axis from entry
        fc = FreysoldtCorrection(15, axis=[1])
        val = fc.get_correction(de)
        self.assertAlmostEqual(val['freysoldt_potential_alignment'], 5.2869010593283132)

        # test a different charge
        #   for electrostatic correction
        es_corr = fc.perform_es_corr(struc.lattice, 2)
        self.assertAlmostEqual(es_corr, 0.43373)
        #   for potential alignment method
        pot_corr = fc.perform_pot_corr(axisdata[0], bldata[0], dldata[0], struc.lattice, 2, vac.site.coords, 0)
        self.assertAlmostEqual(pot_corr, -2.1375685936497768)

        # test an input anisotropic dielectric constant
        fc = FreysoldtCorrection([[1., 2., 3.], [0., 3., 5.], [4., 10., 8.]])
        self.assertAlmostEqual(fc.dielectric, 4.)
        val = fc.get_correction(de)
        self.assertAlmostEqual(val['freysoldt_electrostatic'], 3.659599)
        self.assertAlmostEqual(val['freysoldt_potential_alignment'], 3.3605255195745087)

        # test potalign being added to defect entry
        self.assertAlmostEqual(de.parameters['potalign'], 1.1201751731915028)

        # test that metadata entries exist in defect entry
        self.assertTrue('freysoldt_meta' in de.parameters.keys())
        self.assertAlmostEqual(
            set(de.parameters['freysoldt_meta'].keys()), set(['pot_plot_data', 'pot_corr_uncertainty_md']))

        # test a charge of zero
        vac = Vacancy(struc, struc.sites[0], charge=0)
        de = DefectEntry(vac, 0., corrections={}, parameters=params, entry_id=None)
        val = fc.get_correction(de)
        self.assertAlmostEqual(val['freysoldt_electrostatic'], 0.)
        self.assertAlmostEqual(val['freysoldt_potential_alignment'], 0.)

    def test_kumagai(self):
        gamma = 0.19357221
        prec = 28
        lattice = Lattice([[4.692882, -8.12831, 0.],
                           [4.692882, 8.12831, 0.],
                           [0., 0., 10.03391]])

        # note that real/recip vector generation is not dependent on epsilon
        g_vecs, _, r_vecs, _ = generate_R_and_G_vecs(gamma, prec, lattice, 80. * np.identity(3))

        # test real space summation (bigger for large epsilon)
        kc_high_diel = KumagaiCorrection(80. * np.identity(3), gamma=gamma)
        real_sum = kc_high_diel.get_real_summation(gamma, r_vecs[0])
        self.assertAlmostEqual(real_sum, 0.00843104)

        # test recip space summation (bigger for small epsilon)
        kc_low_diel = KumagaiCorrection(0.1 * np.identity(3), gamma=gamma)
        recip_sum = kc_low_diel.get_recip_summation(gamma, g_vecs[0], lattice.volume)
        self.assertAlmostEqual(recip_sum, 0.31117099)

        # test self interaction
        si_corr = kc_low_diel.get_self_interaction(gamma)
        self.assertAlmostEqual(si_corr, -0.54965249)

        # test potenital shift interaction correction
        ps_corr = kc_low_diel.get_potential_shift(gamma, lattice.volume)
        self.assertAlmostEqual(ps_corr, -0.00871593)

        # """Test Defect Entry approach to correction """
        bulk_struc = Poscar.from_file(os.path.join(test_dir, 'defect', 'CONTCAR_bulk')).structure
        bulk_out = Outcar(os.path.join(test_dir, 'defect', 'OUTCAR_bulk.gz'))
        defect_out = Outcar(os.path.join(test_dir, 'defect', 'OUTCAR_vac_Ga_-3.gz'))
        epsilon = 18.118 * np.identity(3)
        vac = Vacancy(bulk_struc, bulk_struc.sites[0], charge=-3)
        defect_structure = vac.generate_defect_structure()
        defect_frac_coords = [0., 0., 0.]

        parameters = {'bulk_atomic_site_averages': bulk_out.electrostatic_potential,
                      'defect_atomic_site_averages': defect_out.electrostatic_potential,
                      'site_matching_indices': [[ind, ind - 1] for ind in range(len(bulk_struc))],
                      'initial_defect_structure': defect_structure,
                      'defect_frac_sc_coords': defect_frac_coords}
        dentry = DefectEntry(vac, 0., parameters=parameters)
        kc = KumagaiCorrection(epsilon)
        kcorr = kc.get_correction(dentry)
        self.assertAlmostEqual(kcorr['kumagai_electrostatic'], 0.88236299)
        self.assertAlmostEqual(kcorr['kumagai_potential_alignment'], 2.09704862)

        # test ES correction
        high_diel_es_corr = kc_high_diel.perform_es_corr(gamma, prec, lattice, -3.)
        self.assertAlmostEqual(high_diel_es_corr, 0.25176240)

        low_diel_es_corr = kc_low_diel.perform_es_corr(gamma, prec, lattice, -3.)
        self.assertAlmostEqual(low_diel_es_corr, 201.28810966)

        # test pot correction
        site_list = []
        for bs_ind, ds_ind in dentry.parameters['site_matching_indices']:
            Vqb = -(defect_out.electrostatic_potential[ds_ind] - bulk_out.electrostatic_potential[bs_ind])
            site_list.append([defect_structure[ds_ind], Vqb])

        sampling_radius = dentry.parameters["kumagai_meta"]["sampling_radius"]
        gamma = dentry.parameters["kumagai_meta"]["gamma"]
        q = -3
        g_vecs, _, r_vecs, _ = generate_R_and_G_vecs(gamma, 28, defect_structure.lattice, np.identity(3))
        high_diel_pot_corr = kc_high_diel.perform_pot_corr(defect_structure, defect_frac_coords,
                                                           site_list, sampling_radius, q,
                                                           r_vecs[0], g_vecs[0], gamma)
        self.assertAlmostEqual(high_diel_pot_corr, 2.35840716)
        low_diel_pot_corr = kc_low_diel.perform_pot_corr(defect_structure, defect_frac_coords,
                                                         site_list, sampling_radius, q,
                                                         r_vecs[0], g_vecs[0], gamma)
        self.assertAlmostEqual(low_diel_pot_corr, -58.83598095)

        # test the kumagai plotter
        kcp = kc.plot()
        self.assertTrue(kcp)

        # check that uncertainty metadata exists
        self.assertAlmostEqual(set(kc.metadata['pot_corr_uncertainty_md'].keys()), set(['number_sampled', 'stats']))

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

        # test trivial performing bandfilling correction
        bf_corr = bfc.perform_bandfill_corr(eigenvalues, kptweights, potalign, vbm, cbm)
        self.assertAlmostEqual(bf_corr, 0.)
        self.assertFalse(bfc.metadata['num_elec_cbm'])
        self.assertFalse(bfc.metadata['num_hole_vbm'])
        self.assertFalse(bfc.metadata['potalign'])

        # test trivial full entry bandfill evaluation
        de = DefectEntry(vac, 0., corrections={}, parameters=params, entry_id=None)

        corr = bfc.get_correction(de)
        self.assertAlmostEqual(corr['bandfilling_correction'], 0.)

        # modify the eigenvalue list to have free holes
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

        # test case with only one spin and eigen-occupations are 1.
        one_spin_eigen = hole_eigenvalues.copy()
        del one_spin_eigen[list(eigenvalues.keys())[0]]
        bf_corr = bfc.perform_bandfill_corr(one_spin_eigen, kptweights, potalign, vbm, cbm)
        self.assertAlmostEqual(bf_corr, -0.14487501159000005)

        # test case with only one spin and eigen-occupations are 2.
        one_spin_eigen_twooccu = one_spin_eigen.copy()
        for kptset in one_spin_eigen_twooccu.values():
            for bandset in kptset:
                for occuset in bandset:
                    if occuset[1] == 1.:
                        occuset[1] = 2.
                    elif occuset[1] == .5:
                        occuset[1] = 1.
        bf_corr = bfc.perform_bandfill_corr(one_spin_eigen_twooccu, kptweights, potalign, vbm, cbm)
        self.assertAlmostEqual(bf_corr, -0.14487501159000005)

    def test_bandedgeshifting(self):
        struc = PymatgenTest.get_structure("VO2")
        struc.make_supercell(3)
        struc = struc
        vac = Vacancy(struc, struc.sites[0], charge=-3)

        besc = BandEdgeShiftingCorrection()
        params = {'hybrid_cbm': 1., 'hybrid_vbm': -1., 'vbm': -0.5, 'cbm': 0.6}
        de = DefectEntry(vac, 0., corrections={}, parameters=params, entry_id=None)

        corr = besc.get_correction(de)
        self.assertEqual(corr['bandedgeshifting_correction'], 1.5)


if __name__ == "__main__":
    unittest.main()
