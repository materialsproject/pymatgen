# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest
import numpy as np

from pymatgen.util.testing import PymatgenTest
from pymatgen.analysis.defects.core import DefectEntry, Vacancy
from pymatgen.analysis.defects.corrections import FreysoldtCorrection


class DefectsCorrectionsTest(PymatgenTest):
    def test_freysoldt(self):
        #a trivial test for freysoldt correction
        struc = PymatgenTest.get_structure("VO2")
        struc.make_supercell(3)
        abc = struc.lattice.abc
        axisdata = [np.arange(0., lattval, 0.2) for lattval in abc]
        bldata = [np.array([1. for u in np.arange(0., lattval, 0.2)])  for lattval in abc]
        dldata = [np.array([(-1-np.cos(2*np.pi*u/lattval)) for u in np.arange(0., lattval, 0.2)])  for lattval in abc]
        params = {'axis_grid':axisdata, 'bulk_planar_averages': bldata, 'defect_planar_averages':dldata}
        vac = Vacancy(struc, struc.sites[0], charge=-3)
        fc = FreysoldtCorrection(15)

        #test electrostatic correction
        es_corr = fc.perform_es_corr(struc.lattice, 15., -3)
        self.assertEqual(es_corr, 0.976412)

        #test potential alignment method
        pot_corr = fc.perform_pot_corr( axisdata[0], bldata[0], dldata[0], struc.lattice, 15., -3, vac.site.coords, 0)
        self.assertEqual(pot_corr, 2.836369987722345)

        #test entry full correction method
        de = DefectEntry(vac, 0., corrections={}, parameters=params, entry_id=None)
        val = fc.get_correction(de)
        self.assertEqual(val['freysoldt_electrostatic'], 0.976412)
        self.assertEqual(val['freysoldt_potential_alignment'], 4.579907163820233)

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
        #a trivial test for kumagai correction

        #test find_optimal_gamma method
        # gam = find_optimal_gamma(structure, epsilon)

        #test generate_g_sum method
        # g_sum = generate_g_sum( structure, epsilon, dim, gamma)

        # KumagaiCorrection( dielectric_tensor)

        #test electrostatic correction

        #test potential alignment method

        #test entry full correction method (including

        #test the kumagai plotter
        # kumagai_plotter

        pass



if __name__ == "__main__":
    unittest.main()
