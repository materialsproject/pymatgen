# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest
import numpy as np

from pymatgen.util.testing import PymatgenTest
from pymatgen.analysis.defects.core import DefectEntry, Vacancy
from pymatgen.analysis.defects.corrections import FreysoldtCorrection, freysoldt_plotter, KumagaiCorrection, find_optimal_gamma, generate_g_sum, kumagai_plotter


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
        self.assertEqual(es_corr, 0.976412)

        #test potential alignment method
        pot_corr = fc.perform_pot_corr( axisdata[0], bldata[0], dldata[0], struc.lattice, 15., -3, vac.site.coords, 0)
        self.assertEqual(pot_corr, 2.836369987722345)

        #test entry full correction method
        de = DefectEntry( vac, 0., corrections={}, parameters=params, entry_id=None)
        val = fc.get_correction(de)
        self.assertAlmostEqual(val['freysoldt_electrostatic'], 0.976412)
        self.assertAlmostEqual(val['freysoldt_potential_alignment'], 4.4700574)

        #test the freysoldt plotter
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
        self.assertEqual(len(pltsaver), 3)

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
        gamma = find_optimal_gamma(vac.bulk_structure, epsilon)
        self.assertEqual(gamma, 2.006866136855667)

        #test generate_g_sum method
        g_sum = generate_g_sum( vac.bulk_structure, epsilon, dim, gamma)
        self.assertEqual( len(g_sum), 100)

        #test electrostatic correction
        es_corr = kc.perform_es_corr( vac.bulk_structure, -3, g_sum, gamma, kc.madelung_energy_tolerance)
        self.assertEqual(es_corr, 0.9764131640471401)

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
        self.assertEqual(pot_corr, 0.30603437757535096)

        #test entry full correction method
        de = DefectEntry( vac, 0., corrections={}, parameters=params, entry_id=None)
        KC = KumagaiCorrection( epsilon, gamma=gamma, g_sum=g_sum)
        val = KC.get_correction(de)
        self.assertAlmostEqual(val['kumagai_electrostatic'], 0.976413164047)
        self.assertAlmostEqual(val['kumagai_potential_alignment'], 0.2938097394999)

        #test wigner-seitz sampling radius method
        self.assertAlmostEqual(KC.metadata['sampling_radius'], 4.5531438299999)


        #test the kumagai plotter
        eltnames = []
        for bsind in KC.metadata['pot_corr_uncertainty_md']['AllData'].keys():
            eltnames.append( vac.bulk_structure.sites[bsind].specie.symbol)

        eltnames = list(set(eltnames))
        eltkey = {sym: indnum for sym, indnum in zip(eltnames, range(len(eltnames)))}

        rset = [[] for tmp in range(len(eltnames))]
        Vqbset = [[] for tmp in range(len(eltnames))]
        Vpcset = [[] for tmp in range(len(eltnames))]

        for bsind, vals in KC.metadata['pot_corr_uncertainty_md']['AllData'].items():
            elttype = vac.bulk_structure.sites[bsind].specie.symbol
            eltind = eltkey[elttype]
            rset[eltind].append(vals['dist_to_defect'])
            Vqbset[eltind].append(vals['Vqb'])
            Vpcset[eltind].append(vals['Vpc'])
        
        kp = kumagai_plotter(rset, Vqbset, Vpcset, eltnames, samplerad=4.5, title = 'test', saved=False)
        self.assertTrue( kp)



if __name__ == "__main__":
    unittest.main()
