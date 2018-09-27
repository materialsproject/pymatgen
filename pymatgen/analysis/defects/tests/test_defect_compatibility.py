# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest

import os
import numpy as np

from pymatgen.core import PeriodicSite
from pymatgen.io.vasp import Vasprun
from pymatgen.analysis.defects.core import  Vacancy, Interstitial, DefectEntry
from pymatgen.analysis.defects.defect_compatibility import DefectCompatibility
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", 'test_files')


class DefectCompatibilityTest(PymatgenTest):

    def setUp(self):
        struc = PymatgenTest.get_structure("VO2")
        struc.make_supercell(3)
        struc = struc
        self.vac = Vacancy(struc, struc.sites[0], charge=-3)

        abc = self.vac.bulk_structure.lattice.abc
        axisdata = [np.arange(0., lattval, 0.2) for lattval in abc]
        bldata = [np.array([1. for u in np.arange(0., lattval, 0.2)]) for lattval in abc]
        dldata = [
            np.array([(-1 - np.cos(2 * np.pi * u / lattval)) for u in np.arange(0., lattval, 0.2)]) for lattval in abc
        ]
        self.frey_params = {'axis_grid': axisdata, 'bulk_planar_averages': bldata, 'defect_planar_averages': dldata, 'dielectric': 15}

        v = Vasprun(os.path.join(test_dir, 'vasprun.xml'))
        eigenvalues = v.eigenvalues.copy()
        kptweights = v.actual_kpoints_weights
        potalign = -0.1
        vbm = v.eigenvalue_band_properties[2]
        cbm = v.eigenvalue_band_properties[1]
        self.bandfill_params = { 'eigenvalues': eigenvalues,
                                 'kpoint_weights': kptweights,
                                 'potalign': potalign,
                                 'vbm': vbm, 'cbm': cbm }

        self.band_edge_params = {'hybrid_cbm': 1., 'hybrid_vbm': -1., 'vbm': -0.5,
                                 'cbm': 0.6, 'num_hole_vbm': 1., 'num_elec_cbm': 1.}

    def test_process_entry(self):

        # basic process with no corrections
        dentry = DefectEntry(self.vac, 0., corrections={}, parameters={'vbm': 0., 'cbm': 0.}, entry_id=None)
        dc = DefectCompatibility()
        dentry = dc.process_entry( dentry)
        self.assertIsNotNone( dentry)

        # process with corrections from parameters used in other unit tests
        params = self.frey_params.copy()
        params.update(self.bandfill_params)
        params.update({'hybrid_cbm': params['cbm'] + .2, 'hybrid_vbm': params['vbm'] - .4, })
        dentry = DefectEntry(self.vac, 0., corrections={}, parameters=params, entry_id=None)
        dc = DefectCompatibility()
        dentry = dc.process_entry( dentry)
        self.assertAlmostEqual( dentry.corrections['bandedgeshifting_correction'], 1.2)
        self.assertAlmostEqual( dentry.corrections['bandfilling_correction'], 0.0)
        self.assertAlmostEqual( dentry.corrections['charge_correction'], 5.44595036)

        #TODO: add correction analysis with Kumagai correction

        # test over delocalized free carriers which forces skipping charge correction
        # modify the eigenvalue list to have free holes
        hole_eigenvalues = {}
        for spinkey, spinset in params['eigenvalues'].items():
            hole_eigenvalues[spinkey] = []
            for kptset in spinset:
                hole_eigenvalues[spinkey].append([])
                for eig in kptset:
                    if (eig[0] < params['vbm']) and (eig[0] > params['vbm'] - .8):
                        hole_eigenvalues[spinkey][-1].append([eig[0], 0.5])
                    else:
                        hole_eigenvalues[spinkey][-1].append(eig)

        params.update( {'eigenvalues': hole_eigenvalues})
        dentry = DefectEntry(self.vac, 0., corrections={}, parameters=params, entry_id=None)
        dc = DefectCompatibility( free_chg_cutoff=0.8)
        dentry = dc.process_entry( dentry)
        self.assertAlmostEqual( dentry.corrections['bandedgeshifting_correction'], 1.52500002)
        self.assertAlmostEqual( dentry.corrections['bandfilling_correction'], -1.62202400)
        self.assertAlmostEqual( dentry.corrections['charge_correction'], 0.)

        # turn off band filling and band edge shifting
        dc = DefectCompatibility( free_chg_cutoff=0.8, use_bandfilling=False, use_bandedgeshift=False)
        dentry = dc.process_entry( dentry)
        self.assertAlmostEqual( dentry.corrections['bandedgeshifting_correction'], 0.)
        self.assertAlmostEqual( dentry.corrections['bandfilling_correction'], 0.)
        self.assertAlmostEqual( dentry.corrections['charge_correction'], 0.)


    def test_perform_all_corrections(self):

        #return entry even if insufficent values are provided
        # for freysoldt, kumagai, bandfilling, or band edge shifting
        de = DefectEntry(self.vac, 0., corrections={}, parameters={}, entry_id=None)
        dc = DefectCompatibility()
        dentry = dc.perform_all_corrections( de)
        self.assertIsNotNone( dentry)
        #all other correction applications are tested in unit tests below

    def test_perform_freysoldt(self):
        de = DefectEntry(self.vac, 0., corrections={}, parameters=self.frey_params, entry_id=None)
        dc = DefectCompatibility()
        dentry = dc.perform_freysoldt( de)

        val = dentry.parameters['freysoldt_meta']
        self.assertAlmostEqual(val['freysoldt_electrostatic'], 0.975893)
        self.assertAlmostEqual(val['freysoldt_potential_alignment_correction'], 4.4700574)
        self.assertAlmostEqual(val['freysoldt_potalign'], 1.4900191)
        self.assertTrue('pot_corr_uncertainty_md' in val.keys())
        self.assertTrue('pot_plot_data' in val.keys())

    def test_perform_kumagai(self):
        #TODO: add this once correction exists
        pass

    def test_run_bandfilling(self):
        de = DefectEntry(self.vac, 0., corrections={}, parameters=self.bandfill_params, entry_id=None)
        dc = DefectCompatibility()
        dentry = dc.perform_bandfilling( de)

        val = dentry.parameters['bandfilling_meta']
        self.assertAlmostEqual(val['num_hole_vbm'], 0.)
        self.assertAlmostEqual(val['num_elec_cbm'], 0.)
        self.assertAlmostEqual(val['bandfilling_correction'], 0.)
        occu = [[1.457, 0.0833333], [1.5204, 0.0833333], [1.53465, 0.0833333], [1.5498, 0.0416667]]
        self.assertArrayAlmostEqual(
            list(sorted(val['occupied_def_levels'], key=lambda x: x[0])), list(
                sorted(occu, key=lambda x: x[0])))
        self.assertAlmostEqual(val['total_occupation_defect_levels'], 0.29166669)
        self.assertFalse(val['unoccupied_def_levels'])

    def test_run_band_edge_shifting(self):
        de = DefectEntry(self.vac, 0., corrections={}, parameters=self.band_edge_params, entry_id=None)

        dc = DefectCompatibility()
        dentry = dc.perform_band_edge_shifting( de)
        val = dentry.parameters['bandshift_meta']
        self.assertEqual(val['vbmshift'], -0.5)
        self.assertEqual(val['cbmshift'], 0.4)
        self.assertEqual(val['vbm_shift_correction'], 1.5)
        self.assertEqual(val['hole_vbm_shift_correction'], 0.5)
        self.assertEqual(val['elec_cbm_shift_correction'], 0.4)

    def test_delocalization_analysis(self):
        #return entry even if insufficent values are provided
        # for delocalization analysis with freysoldt, kumagai,
        # bandfilling, or band edge shifting
        de = DefectEntry(self.vac, 0., corrections={}, parameters={}, entry_id=None)
        dc = DefectCompatibility()
        dentry = dc.delocalization_analysis( de)
        self.assertIsNotNone( dentry)
        #all other correction applications are tested in unit tests below

    def test_is_freysoldt_delocalized(self):
        de = DefectEntry(self.vac, 0., corrections={}, parameters=self.frey_params, entry_id=None)
        de.parameters.update( {'is_compatible': True}) #needs to be initialized with this here
        dc = DefectCompatibility( plnr_avg_var_tol=0.1, plnr_avg_minmax_tol=0.5)
        dentry = dc.perform_freysoldt( de)

        # check case which fits under compatibility constraints
        dentry = dc.is_freysoldt_delocalized( dentry)
        frey_delocal = dentry.parameters['delocalization_meta']['plnr_avg']
        self.assertTrue( frey_delocal['is_compatible'])
        ans_var = [0.00038993, 0.02119532, 0.02119532]
        ans_window = [0.048331509, 0.36797169, 0.36797169]
        for ax in range(3):
            ax_metadata = frey_delocal['metadata'][ax]
            self.assertTrue( ax_metadata['frey_variance_compatible'])
            self.assertAlmostEqual( ax_metadata['frey_variance'], ans_var[ax])
            self.assertTrue( ax_metadata['frey_minmax_compatible'])
            self.assertAlmostEqual( ax_metadata['frey_minmax_window'], ans_window[ax])

        self.assertTrue( dentry.parameters['is_compatible'])

        # check planar delocalization on 2nd and 3rd axes
        dc = DefectCompatibility( plnr_avg_var_tol=0.1, plnr_avg_minmax_tol=0.2)
        dentry.parameters.update( {'is_compatible': True})
        dentry = dc.is_freysoldt_delocalized( dentry)
        frey_delocal = dentry.parameters['delocalization_meta']['plnr_avg']
        self.assertFalse( frey_delocal['is_compatible'])
        ax_metadata = frey_delocal['metadata'][0]
        self.assertTrue( ax_metadata['frey_variance_compatible'])
        self.assertTrue( ax_metadata['frey_minmax_compatible'])
        for ax in [1,2]:
            ax_metadata = frey_delocal['metadata'][ax]
            self.assertTrue( ax_metadata['frey_variance_compatible'])
            self.assertFalse( ax_metadata['frey_minmax_compatible'])

        self.assertFalse( dentry.parameters['is_compatible'])

        # check variance based delocalization on 2nd and 3rd axes
        dc = DefectCompatibility( plnr_avg_var_tol=0.01, plnr_avg_minmax_tol=0.5)
        dentry.parameters.update( {'is_compatible': True})
        dentry = dc.is_freysoldt_delocalized( dentry)
        frey_delocal = dentry.parameters['delocalization_meta']['plnr_avg']
        self.assertFalse( frey_delocal['is_compatible'])
        ax_metadata = frey_delocal['metadata'][0]
        self.assertTrue( ax_metadata['frey_variance_compatible'])
        self.assertTrue( ax_metadata['frey_minmax_compatible'])
        for ax in [1,2]:
            ax_metadata = frey_delocal['metadata'][ax]
            self.assertFalse( ax_metadata['frey_variance_compatible'])
            self.assertTrue( ax_metadata['frey_minmax_compatible'])

        self.assertFalse( dentry.parameters['is_compatible'])

    def test_is_kumagai_delocalized(self):
        #TODO: add this once correction exists
        pass

    def test_is_final_relaxed_structure_delocalized(self):
        # test structure delocalization analysis
        # first test no movement in atoms
        initial_defect_structure = self.vac.generate_defect_structure()
        final_defect_structure = initial_defect_structure.copy()
        sampling_radius = 4.55

        params = {'initial_defect_structure': initial_defect_structure,
                  'final_defect_structure': final_defect_structure,
                  'sampling_radius': sampling_radius,
                  'is_compatible': True}
        dentry = DefectEntry(self.vac, 0., corrections={}, parameters=params, entry_id=None)

        dc = DefectCompatibility( tot_relax_tol=0.1, perc_relax_tol=0.1, defect_tot_relax_tol=0.1)
        dentry = dc.is_final_relaxed_structure_delocalized( dentry)

        struc_delocal = dentry.parameters['delocalization_meta']['structure_relax']
        self.assertTrue( dentry.parameters['is_compatible'])
        self.assertTrue( struc_delocal['is_compatible'])
        self.assertTrue( struc_delocal['metadata']['structure_tot_relax_compatible'])
        self.assertEqual( struc_delocal['metadata']['tot_relax_outside_rad'], 0.)
        self.assertTrue( struc_delocal['metadata']['structure_perc_relax_compatible'])
        self.assertEqual( struc_delocal['metadata']['perc_relax_outside_rad'], 0.)
        self.assertEqual( len(struc_delocal['metadata']['full_structure_relax_data']), len(initial_defect_structure))
        self.assertIsNone( struc_delocal['metadata']['defect_index'])

        defect_delocal = dentry.parameters['delocalization_meta']['defectsite_relax']
        self.assertTrue( defect_delocal['is_compatible'])
        self.assertIsNone( defect_delocal['metadata']['relax_amount'])


        # next test for when structure has delocalized outside of radius from defect
        pert_struct_fin_struct = initial_defect_structure.copy()
        pert_struct_fin_struct.perturb( 0.1)
        dentry.parameters.update( {'final_defect_structure': pert_struct_fin_struct})
        dentry = dc.is_final_relaxed_structure_delocalized( dentry)

        struc_delocal = dentry.parameters['delocalization_meta']['structure_relax']
        self.assertFalse( dentry.parameters['is_compatible'])
        self.assertFalse( struc_delocal['is_compatible'])
        self.assertFalse( struc_delocal['metadata']['structure_tot_relax_compatible'])
        self.assertAlmostEqual( struc_delocal['metadata']['tot_relax_outside_rad'], 12.5)
        self.assertFalse( struc_delocal['metadata']['structure_perc_relax_compatible'])
        self.assertAlmostEqual( struc_delocal['metadata']['perc_relax_outside_rad'], 77.63975155)


        # now test for when an interstitial defect has migrated too much
        inter_def_site = PeriodicSite('H',  [7.58857304, 11.70848069, 12.97817518],
                                self.vac.bulk_structure.lattice, to_unit_cell=True,
                                coords_are_cartesian=True)
        inter = Interstitial(self.vac.bulk_structure, inter_def_site, charge=0)

        initial_defect_structure = inter.generate_defect_structure()
        final_defect_structure = initial_defect_structure.copy()
        poss_deflist = sorted(
            final_defect_structure.get_sites_in_sphere(inter.site.coords,
                                                       2, include_index=True), key=lambda x: x[1])
        def_index = poss_deflist[0][2]
        final_defect_structure.translate_sites(indices=[def_index],
                                               vector=[0., 0., 0.008]) #fractional coords translation

        params = {'initial_defect_structure': initial_defect_structure,
                  'final_defect_structure': final_defect_structure,
                  'sampling_radius': sampling_radius,
                  'is_compatible': True}
        dentry = DefectEntry(inter, 0., corrections={}, parameters=params, entry_id=None)

        dentry = dc.is_final_relaxed_structure_delocalized( dentry)

        defect_delocal = dentry.parameters['delocalization_meta']['defectsite_relax']
        self.assertFalse( defect_delocal['is_compatible'])
        self.assertAlmostEqual( defect_delocal['metadata']['relax_amount'], 0.10836054)

if __name__ == "__main__":
    unittest.main()
