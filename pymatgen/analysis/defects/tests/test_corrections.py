# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import os
import unittest
import numpy as np

from pymatgen.util.testing import PymatgenTest
from pymatgen.io.vasp import Vasprun
from pymatgen.analysis.defects.core import DefectEntry, Vacancy
from pymatgen.analysis.defects.corrections import FreysoldtCorrection, freysoldt_plotter,\
            KumagaiCorrection, find_optimal_gamma, generate_g_sum, kumagai_plotter, \
            BandFillingCorrection, BandEdgeShiftingCorrection

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')


class DefectsCorrectionsTest(PymatgenTest):
    def test_freysoldt(self):
        struc = PymatgenTest.get_structure("VO2")
        struc.make_supercell(3)
        struc = struc
        vac = Vacancy(struc, struc.sites[0], charge=-3)

        abc = struc.lattice.abc
        axisdata = [np.arange(0., lattval, 0.2) for lattval in abc]
        bldata = [np.array([1. for u in np.arange(0., lattval, 0.2)])  for lattval in abc]
        dldata = [np.array([(-1-np.cos(2*np.pi*u/lattval)) for u in np.arange(0., lattval, 0.2)])  for lattval in abc]
        params = {'axis_grid':axisdata, 'bulk_planar_averages': bldata, 'defect_planar_averages':dldata}
        fc = FreysoldtCorrection(15)

        #test electrostatic correction
        es_corr = fc.perform_es_corr(struc.lattice, 15., -3)
        self.assertAlmostEqual(es_corr, 0.976412)

        #test potential alignment method
        pot_corr = fc.perform_pot_corr( axisdata[0], bldata[0], dldata[0], struc.lattice, 15., -3, vac.site.coords, 0)
        self.assertAlmostEqual(pot_corr, 2.836369987722345)

        #test entry full correction method
        de = DefectEntry( vac, 0., corrections={}, parameters=params, entry_id=None)
        val = fc.get_correction(de)
        self.assertAlmostEqual(val['freysoldt_electrostatic'], 0.976412)
        self.assertAlmostEqual(val['freysoldt_potential_alignment'], 4.4700574)

        #test the freysoldt plotter and that plot metadata exists
        pltsaver = []
        for ax in range(3):
            x = fc.metadata['pot_plot_data'][ax]['x']
            Vr = fc.metadata['pot_plot_data'][ax]['Vr']
            dft_diff = fc.metadata['pot_plot_data'][ax]['dft_diff']
            final_shift = fc.metadata['pot_plot_data'][ax]['final_shift']
            check = fc.metadata['pot_plot_data'][ax]['check']
            fp = freysoldt_plotter( x, Vr, dft_diff, final_shift, check, title='axis_'+str(ax+1), saved=False)
            if fp: #if plot exists then append it...
                pltsaver.append(fp)
        self.assertAlmostEqual(len(pltsaver), 3)

        #check that uncertainty metadata exists
        for ax in range(3):
            self.assertAlmostEqual(set(fc.metadata['pot_corr_uncertainty_md'][ax].keys()), set(['potcorr', 'stats']))

        #test a specified axis from entry
        fc = FreysoldtCorrection(15, axis = [1])
        val = fc.get_correction(de)
        self.assertAlmostEqual(val['freysoldt_potential_alignment'], 5.2869010593283132)

        #test a different charge
        #   for electrostatic correction
        es_corr = fc.perform_es_corr(struc.lattice, 15., 2)
        self.assertAlmostEqual(es_corr, 0.43396099999999999)
        #   for potential alignment method
        pot_corr = fc.perform_pot_corr( axisdata[0], bldata[0], dldata[0], struc.lattice, 15., 2, vac.site.coords, 0)
        self.assertAlmostEqual(pot_corr, -2.1375685936497768)

        #test an input anisotropic dielectric constant
        fc = FreysoldtCorrection([[1.,2.,3.],[0.,3.,5.],[4., 10., 8.]])
        self.assertAlmostEqual( fc.dielectric, 4.)
        val = fc.get_correction(de)
        self.assertAlmostEqual(val['freysoldt_electrostatic'], 3.6615440000000001)
        self.assertAlmostEqual(val['freysoldt_potential_alignment'], 3.3605255195745087)

        #test potalign being added to defect entry
        self.assertAlmostEqual( de.parameters['potalign'], 1.1201751731915028)

        #test that metadata entries exist in defect entry
        self.assertTrue( 'freysoldt_meta' in de.parameters.keys())
        self.assertAlmostEqual( set(de.parameters['freysoldt_meta'].keys()), set(['pot_plot_data', 'pot_corr_uncertainty_md']))

        #test a charge of zero
        vac = Vacancy(struc, struc.sites[0], charge=0)
        de = DefectEntry( vac, 0., corrections={}, parameters=params, entry_id=None)
        val = fc.get_correction(de)
        self.assertAlmostEqual(val['freysoldt_electrostatic'], 0.)
        self.assertAlmostEqual(val['freysoldt_potential_alignment'], 0.)


    def test_kumagai(self):
        struc = PymatgenTest.get_structure("VO2")
        struc.make_supercell(3)
        struc = struc
        vac = Vacancy(struc, struc.sites[0], charge=-3)

        dim = [100, 100, 100]
        bulk_atomic_site_averages = [0.4 for u in range(len(vac.bulk_structure))]
        defect_atomic_site_averages = [0.2  for u in range(len(vac.bulk_structure) - 1)]
        site_matching_indices = [[(u+1),u] for u in range(len(vac.bulk_structure) - 1)]
        params = {'dim': dim, 'bulk_atomic_site_averages': bulk_atomic_site_averages,
                  'defect_atomic_site_averages':defect_atomic_site_averages,
                  'site_matching_indices': site_matching_indices}
        epsilon = np.identity(3) * 15.
        kc = KumagaiCorrection(epsilon)

        #test find_optimal_gamma method
        gamma = find_optimal_gamma(vac.bulk_structure.lattice, epsilon)
        self.assertAlmostEqual(gamma, 2.006866136855667)

        #test generate_g_sum method
        g_sum = generate_g_sum( vac.bulk_structure.lattice, epsilon, dim, gamma)
        self.assertAlmostEqual( len(g_sum), 100)

        #test electrostatic correction
        es_corr = kc.perform_es_corr( vac.bulk_structure.lattice, -3, g_sum, gamma, kc.madelung_energy_tolerance)
        self.assertAlmostEqual(es_corr, 0.9764131640471401)

        #test potential alignment method
        sampling_radius = 4.5
        defect_structure = vac.bulk_structure.copy()
        del defect_structure[0]
        defect_position = vac.bulk_structure[0]
        site_list = [] #assemble site_list ...
        for bs_ind, ds_ind in site_matching_indices:
            Vqb = - (defect_atomic_site_averages[ds_ind] - bulk_atomic_site_averages[bs_ind])
            site_list.append([ vac.bulk_structure[bs_ind], defect_structure[ds_ind], Vqb])

        pot_corr = kc.perform_pot_corr( vac.bulk_structure, defect_structure, defect_position, site_list, sampling_radius, -3, g_sum, dim, gamma,  kc.madelung_energy_tolerance)
        self.assertAlmostEqual(pot_corr, 0.30603437757535096)

        #test entry full correction method
        de = DefectEntry( vac, 0., corrections={}, parameters=params, entry_id=None)
        kc = KumagaiCorrection( epsilon, gamma=gamma, g_sum=g_sum)
        val = kc.get_correction(de)
        self.assertAlmostEqual(val['kumagai_electrostatic'], 0.976413164047)
        self.assertAlmostEqual(val['kumagai_potential_alignment'], 0.2938097394999)

        #test wigner-seitz sampling radius method
        self.assertAlmostEqual(kc.metadata['sampling_radius'], 4.5531438299999)

        #test the kumagai plotter
        eltnames = []
        for bsind in kc.metadata['pot_plot_data'].keys():
            eltnames.append( vac.bulk_structure.sites[bsind].specie.symbol)

        eltnames = list(set(eltnames))
        eltkey = {sym: indnum for sym, indnum in zip(eltnames, range(len(eltnames)))}

        rset = [[] for tmp in range(len(eltnames))]
        Vqbset = [[] for tmp in range(len(eltnames))]
        Vpcset = [[] for tmp in range(len(eltnames))]

        for bsind, vals in kc.metadata['pot_plot_data'].items():
            elttype = vac.bulk_structure.sites[bsind].specie.symbol
            eltind = eltkey[elttype]
            rset[eltind].append(vals['dist_to_defect'])
            Vqbset[eltind].append(vals['Vqb'])
            Vpcset[eltind].append(vals['Vpc'])

        kp = kumagai_plotter(rset, Vqbset, Vpcset, eltnames, samplerad=4.5, title = 'test', saved=False)
        self.assertTrue( kp)

        #check that uncertainty metadata exists
        self.assertAlmostEqual(set(kc.metadata['pot_corr_uncertainty_md'].keys()), set(['number_sampled', 'stats']))
        self.assertAlmostEqual(kc.metadata['pot_corr_uncertainty_md']['number_sampled'], 125)

        #test a different sampling radius
        new_sampling_radius = 8.
        pot_corr = kc.perform_pot_corr( vac.bulk_structure, defect_structure, defect_position, site_list,
                                        new_sampling_radius, -3, g_sum, dim, gamma,  kc.madelung_energy_tolerance)
        self.assertAlmostEqual(pot_corr, 0.021267067525360898)
        self.assertAlmostEqual(kc.metadata['pot_corr_uncertainty_md']['number_sampled'], 21)
        #   also test sampling radius from entry
        kc = KumagaiCorrection( epsilon, sampling_radius=new_sampling_radius, gamma=gamma, g_sum=g_sum)
        val = kc.get_correction(de)
        self.assertAlmostEqual(val['kumagai_potential_alignment'], 0.021267067525360898)

        #test a different charge
        #   for electrostatic correction
        es_corr = kc.perform_es_corr( vac.bulk_structure.lattice, 2, g_sum, gamma, kc.madelung_energy_tolerance)
        self.assertAlmostEqual(es_corr, 0.43396140624317336)
        #   for potential alignment method
        pot_corr = kc.perform_pot_corr( vac.bulk_structure, defect_structure, defect_position, site_list,
                                        sampling_radius, 2, g_sum, dim, gamma,  kc.madelung_energy_tolerance)
        self.assertAlmostEqual(pot_corr, -0.53065138774428844)

        #test an input anisotropic dielectric constant
        aniso_dielectric = np.array([[15.,0,3.],[0,15.,0.],[0,0,10.]])
        aniso_gamma = 2.0068661368556668
        aniso_g_sum = generate_g_sum( vac.bulk_structure.lattice, aniso_dielectric, dim, aniso_gamma)
        kc = KumagaiCorrection(aniso_dielectric, gamma=aniso_gamma, g_sum=aniso_g_sum)
        for u,v in zip(kc.dielectric.flatten(), aniso_dielectric.flatten()):
            self.assertAlmostEqual(u, v)
        es_corr = kc.perform_es_corr( vac.bulk_structure.lattice, -3, aniso_g_sum, aniso_gamma, kc.madelung_energy_tolerance)
        self.assertAlmostEqual(es_corr, 1.027001049757768)
        pot_corr = kc.perform_pot_corr( vac.bulk_structure, defect_structure, defect_position,
                                        site_list, sampling_radius, -3, aniso_g_sum, dim, aniso_gamma,  kc.madelung_energy_tolerance)
        self.assertAlmostEqual(pot_corr, 0.2457398743896354)

        #test potalign being added to defect entry
        self.assertAlmostEqual( de.parameters['potalign'],0.0070890225084536329)

        #test that metadata entries exist in defect entry
        self.assertTrue( 'kumagai_meta' in de.parameters.keys())
        self.assertAlmostEqual( set(de.parameters['kumagai_meta'].keys()), set(['pot_plot_data', 'sampling_radius', 'gamma','pot_corr_uncertainty_md']))

        #test a charge of zero
        vac = Vacancy(struc, struc.sites[0], charge=0)
        de = DefectEntry( vac, 0., corrections={}, parameters=params, entry_id=None)
        val = kc.get_correction(de)
        self.assertAlmostEqual(val['kumagai_electrostatic'], 0.)
        self.assertAlmostEqual(val['kumagai_potential_alignment'], 0.)


    def test_bandfilling(self):
        v = Vasprun(os.path.join(test_dir, 'vasprun.xml'))
        eigenvalues = v.eigenvalues.copy()
        kptweights = v.actual_kpoints_weights
        potalign = 0.
        vbm = v.eigenvalue_band_properties[2]
        cbm = v.eigenvalue_band_properties[1]
        params = {'eigenvalues': eigenvalues, 'kpoint_weights': kptweights, 'potalign': potalign,
                  'vbm': vbm, 'cbm': cbm}
        bfc = BandFillingCorrection()
        struc = PymatgenTest.get_structure("VO2")
        struc.make_supercell(3)
        struc = struc
        vac = Vacancy(struc, struc.sites[0], charge=-3)

        #test trivial performing bandfilling correction
        bf_corr = bfc.perform_bandfill_corr( eigenvalues, kptweights, potalign, vbm, cbm)
        self.assertAlmostEqual(bf_corr, 0.)
        self.assertFalse(bfc.metadata['occupied_def_levels'])
        self.assertFalse(bfc.metadata['unoccupied_def_levels'])
        self.assertFalse(bfc.metadata['total_occupation_defect_levels'])
        self.assertFalse(bfc.metadata['num_elec_cbm'])
        self.assertFalse(bfc.metadata['num_hole_vbm'])
        self.assertFalse(bfc.metadata['potalign'])

        #test trivial full entry bandfill evaluation
        de = DefectEntry( vac, 0., corrections={}, parameters=params, entry_id=None)

        corr = bfc.get_correction( de)
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

        hole_bf_corr = bfc.perform_bandfill_corr( hole_eigenvalues, kptweights, potalign, vbm, cbm)
        self.assertAlmostEqual(hole_bf_corr, -0.82276673248)
        self.assertAlmostEqual(bfc.metadata['num_hole_vbm'], 1.6250001299)
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

        elec_bf_corr = bfc.perform_bandfill_corr( elec_eigenvalues, kptweights, potalign, vbm, cbm)
        self.assertAlmostEqual(elec_bf_corr, -0.18063751445099)
        self.assertAlmostEqual(bfc.metadata['num_elec_cbm'], 1.708333469999)
        self.assertFalse(bfc.metadata['num_hole_vbm'])

        #modify the potalignment and introduce new occupied defect levels from vbm states
        potalign = 0.1

        bf_corr = bfc.perform_bandfill_corr( eigenvalues, kptweights, potalign, vbm, cbm)
        self.assertAlmostEqual(bfc.metadata['num_hole_vbm'], 0.)
        self.assertAlmostEqual(bf_corr, 0.)
        occu = [[1.5569999999999999, 0.16666668000000001], [1.6346500000000002, 0.16666668000000001],
                [1.6498000000000002, 0.083333340000000006], [1.6204000000000001, 0.16666668000000001]]
        self.assertArrayAlmostEqual(list(sorted(bfc.metadata['occupied_def_levels'], key=lambda x: x[0])),
                                    list(sorted(occu, key=lambda x: x[0])))
        self.assertAlmostEqual(bfc.metadata['total_occupation_defect_levels'], 0.58333338)
        self.assertFalse(bfc.metadata['unoccupied_def_levels'])


    def test_bandedgeshifting(self):
        struc = PymatgenTest.get_structure("VO2")
        struc.make_supercell(3)
        struc = struc
        vac = Vacancy(struc, struc.sites[0], charge=-3)

        besc = BandEdgeShiftingCorrection()
        params = {'hybrid_cbm': 1., 'hybrid_vbm': -1.,
                  'vbm': -0.5, 'cbm': 0.6,
                  'num_hole_vbm': 0., 'num_elec_cbm': 0.}
        de = DefectEntry( vac, 0., corrections={}, parameters=params, entry_id=None)

        #test with no free carriers
        corr = besc.get_correction(de)
        self.assertEqual( corr['vbm_shift_correction'], 1.5)
        self.assertEqual( corr['elec_cbm_shift_correction'], 0.)
        self.assertEqual( corr['hole_vbm_shift_correction'], 0.)

        #test with free holes
        de.parameters.update( {'num_hole_vbm': 1.})
        corr = besc.get_correction(de)
        self.assertEqual( corr['vbm_shift_correction'], 1.5)
        self.assertEqual( corr['elec_cbm_shift_correction'], 0.)
        self.assertEqual( corr['hole_vbm_shift_correction'], 0.5)

        #test with free electrons
        de.parameters.update( {'num_hole_vbm': 0., 'num_elec_cbm': 1.})
        corr = besc.get_correction(de)
        self.assertEqual( corr['vbm_shift_correction'], 1.5)
        self.assertEqual( corr['elec_cbm_shift_correction'], 0.4)
        self.assertEqual( corr['hole_vbm_shift_correction'], 0.)




if __name__ == "__main__":
    unittest.main()
