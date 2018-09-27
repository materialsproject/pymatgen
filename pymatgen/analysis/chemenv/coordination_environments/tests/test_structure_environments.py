#!/usr/bin/env python


__author__ = 'waroquiers'

import unittest
import os
import json
import numpy as np
import shutil
from pymatgen.util.testing import PymatgenTest

from pymatgen.analysis.chemenv.coordination_environments.structure_environments import StructureEnvironments
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
from pymatgen.core.periodic_table import Specie
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import MultiWeightsChemenvStrategy

se_files_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..",
                            'test_files', "chemenv", "structure_environments_files")


class StructureEnvironmentsTest(PymatgenTest):

    @classmethod
    def setUpClass(cls):
        os.makedirs('tmp_dir')

    def test_structure_environments(self):
        f = open("{}/{}".format(se_files_dir, 'se_mp-7000.json'), 'r')
        dd = json.load(f)
        f.close()

        se = StructureEnvironments.from_dict(dd)
        isite = 6
        csm_and_maps_fig, csm_and_maps_subplot = se.get_csm_and_maps(isite=isite)
        np.testing.assert_array_almost_equal(csm_and_maps_subplot.lines[0].get_xydata().flatten(), [0.0, 0.53499332])
        np.testing.assert_array_almost_equal(csm_and_maps_subplot.lines[1].get_xydata().flatten(), [1.0, 0.47026441])
        np.testing.assert_array_almost_equal(csm_and_maps_subplot.lines[2].get_xydata().flatten(), [2.0, 0.00988778])

        environments_figure, environments_subplot = se.get_environments_figure(isite=isite)
        np.testing.assert_array_almost_equal(np.array(environments_subplot.patches[0].get_xy()),
                                             [[1., 1.],
                                              [1., 0.99301365],
                                              [1.00179228, 0.99301365],
                                              [1.00179228, 1.],
                                              [1., 1.]])
        np.testing.assert_array_almost_equal(np.array(environments_subplot.patches[1].get_xy()),
                                             [[1., 0.99301365],
                                              [1., 0.],
                                              [1.00179228, 0.],
                                              [1.00179228, 0.99301365],
                                              [1., 0.99301365]])
        np.testing.assert_array_almost_equal(np.array(environments_subplot.patches[2].get_xy()),
                                             [[1.00179228, 1.],
                                              [1.00179228, 0.99301365],
                                              [2.25, 0.99301365],
                                              [2.25, 1.],
                                              [1.00179228, 1.]])
        np.testing.assert_array_almost_equal(np.array(environments_subplot.patches[3].get_xy()),
                                             [[1.00179228, 0.99301365],
                                              [1.00179228, 0.],
                                              [2.22376156, 0.],
                                              [2.22376156, 0.0060837],
                                              [2.25, 0.0060837],
                                              [2.25, 0.99301365],
                                              [1.00179228, 0.99301365]])
        np.testing.assert_array_almost_equal(np.array(environments_subplot.patches[4].get_xy()),
                                             [[2.22376156, 0.0060837],
                                              [2.22376156, 0.],
                                              [2.25, 0.],
                                              [2.25, 0.0060837],
                                              [2.22376156, 0.0060837]])

        imgpath = os.path.join('tmp_dir', 'image.png')
        se.save_environments_figure(isite=isite, imagename=imgpath)
        self.assertTrue(os.path.exists(imgpath))

        self.assertEqual(len(se.differences_wrt(se)), 0)

        self.assertFalse(se.__ne__(se))

        ce = se.ce_list[isite][4][0]

        self.assertTrue(ce.__len__(), 4)

        symbol, mingeom = ce.minimum_geometry(symmetry_measure_type='csm_wocs_ctwocc')
        self.assertEqual(symbol, 'T:4')
        self.assertAlmostEqual(mingeom['symmetry_measure'], 0.00988778424054)

        np.testing.assert_array_almost_equal(mingeom['other_symmetry_measures']['rotation_matrix_wcs_csc'],
                                             [[-0.8433079817973094, -0.19705747216466898, 0.5000000005010193],
                                              [0.4868840909509757, 0.11377118475194581, 0.8660254034951744],
                                              [-0.22754236927612112, 0.9737681809261427, 1.3979531202869064e-13]])
        self.assertEqual(mingeom['detailed_voronoi_index'], {'index': 0, 'cn': 4})
        self.assertAlmostEqual(mingeom['other_symmetry_measures']['scaling_factor_wocs_ctwocc'], 1.6270605877934026)

        ce_string = ce.__str__()

        self.assertTrue('csm1 (with central site) : 0.00988' in ce_string)
        self.assertTrue('csm2 (without central site) : 0.00981' in ce_string)
        self.assertTrue('csm1 (with central site) : 12.987' in ce_string)
        self.assertTrue('csm2 (without central site) : 11.827' in ce_string)
        self.assertTrue('csm1 (with central site) : 32.466' in ce_string)
        self.assertTrue('csm2 (without central site) : 32.466' in ce_string)
        self.assertTrue('csm1 (with central site) : 34.644' in ce_string)
        self.assertTrue('csm2 (without central site) : 32.466' in ce_string)

        mingeoms = ce.minimum_geometries(symmetry_measure_type='csm_wocs_ctwocc', max_csm=12.0)
        self.assertEqual(len(mingeoms), 2)
        mingeoms = ce.minimum_geometries(symmetry_measure_type='csm_wocs_ctwcc', max_csm=12.0)
        self.assertEqual(len(mingeoms), 1)
        mingeoms = ce.minimum_geometries(n=3)
        self.assertEqual(len(mingeoms), 3)

        ce2 = se.ce_list[7][4][0]

        self.assertTrue(ce.is_close_to(ce2, rtol=0.01, atol=1e-4))
        self.assertFalse(ce.is_close_to(ce2, rtol=0.0, atol=1e-8))

        self.assertFalse(ce.__eq__(ce2))
        self.assertTrue(ce.__ne__(ce2))

    def test_light_structure_environments(self):
        f = open("{}/{}".format(se_files_dir, 'se_mp-7000.json'), 'r')
        dd = json.load(f)
        f.close()

        se = StructureEnvironments.from_dict(dd)

        strategy = SimplestChemenvStrategy()
        lse = LightStructureEnvironments.from_structure_environments(structure_environments=se, strategy=strategy,
                                                                     valences='undefined')
        isite = 6
        nb_set = lse.neighbors_sets[isite][0]
        neighb_coords = [np.array([0.2443798, 1.80409653, -1.13218359]),
                         np.array([1.44020353, 1.11368738, 1.13218359]),
                         np.array([2.75513098, 2.54465207, -0.70467298]),
                         np.array([0.82616785, 3.65833945, 0.70467298])]
        neighb_indices = [0, 3, 5, 1]
        neighb_images = [[0, 0, -1], [0, 0, 0], [0, 0, -1], [0, 0, 0]]

        np.testing.assert_array_almost_equal(neighb_coords, nb_set.neighb_coords)
        np.testing.assert_array_almost_equal(neighb_coords, [s.coords for s in nb_set.neighb_sites])
        nb_sai = nb_set.neighb_sites_and_indices
        np.testing.assert_array_almost_equal(neighb_coords, [sai['site'].coords for sai in nb_sai])
        np.testing.assert_array_almost_equal(neighb_indices, [sai['index'] for sai in nb_sai])
        nb_iai = nb_set.neighb_indices_and_images
        np.testing.assert_array_almost_equal(neighb_indices, [iai['index'] for iai in nb_iai])
        np.testing.assert_array_equal(neighb_images, [iai['image_cell'] for iai in nb_iai])

        self.assertEqual(nb_set.__len__(), 4)
        self.assertEqual(nb_set.__hash__(), 4)

        self.assertFalse(nb_set.__ne__(nb_set))

        self.assertEqual(nb_set.__str__(), 'Neighbors Set for site #6 :\n'
                                           ' - Coordination number : 4\n'
                                           ' - Neighbors sites indices : 0, 1, 2, 3\n')

        stats = lse.get_statistics()

        neighbors = lse.strategy.get_site_neighbors(site=lse.structure[isite])
        self.assertArrayAlmostEqual(neighbors[0].coords, np.array([ 0.2443798, 1.80409653, -1.13218359]))
        self.assertArrayAlmostEqual(neighbors[1].coords, np.array([ 1.44020353, 1.11368738, 1.13218359]))
        self.assertArrayAlmostEqual(neighbors[2].coords, np.array([ 2.75513098, 2.54465207, -0.70467298]))
        self.assertArrayAlmostEqual(neighbors[3].coords, np.array([ 0.82616785, 3.65833945, 0.70467298]))

        equiv_site_index_and_transform = lse.strategy.equivalent_site_index_and_transform(neighbors[0])
        self.assertEqual(equiv_site_index_and_transform[0], 0)
        self.assertArrayAlmostEqual(equiv_site_index_and_transform[1], [0.0, 0.0, 0.0])
        self.assertArrayAlmostEqual(equiv_site_index_and_transform[2], [0.0, 0.0, -1.0])

        equiv_site_index_and_transform = lse.strategy.equivalent_site_index_and_transform(neighbors[1])
        self.assertEqual(equiv_site_index_and_transform[0], 3)
        self.assertArrayAlmostEqual(equiv_site_index_and_transform[1], [0.0, 0.0, 0.0])
        self.assertArrayAlmostEqual(equiv_site_index_and_transform[2], [0.0, 0.0, 0.0])

        self.assertEqual(stats['atom_coordination_environments_present'], {'Si': {'T:4': 3.0}})
        self.assertEqual(stats['coordination_environments_atom_present'], {'T:4': {'Si': 3.0}})
        self.assertEqual(stats['fraction_atom_coordination_environments_present'], {'Si': {'T:4': 1.0}})

        site_info_ce = lse.get_site_info_for_specie_ce(specie=Specie('Si', 4), ce_symbol='T:4')
        np.testing.assert_array_almost_equal(site_info_ce['fractions'], [1.0, 1.0, 1.0])
        np.testing.assert_array_almost_equal(site_info_ce['csms'],
                                             [0.009887784240541068, 0.009887786546730826, 0.009887787384385317])
        self.assertEqual(site_info_ce['isites'], [6, 7, 8])

        site_info_allces = lse.get_site_info_for_specie_allces(specie=Specie('Si', 4))

        self.assertEqual(site_info_allces['T:4'], site_info_ce)

        self.assertFalse(lse.contains_only_one_anion('I-'))
        self.assertFalse(lse.contains_only_one_anion_atom('I'))
        self.assertTrue(lse.site_contains_environment(isite=isite, ce_symbol='T:4'))
        self.assertFalse(lse.site_contains_environment(isite=isite, ce_symbol='S:4'))
        self.assertFalse(lse.structure_contains_atom_environment(atom_symbol='Si', ce_symbol='S:4'))
        self.assertTrue(lse.structure_contains_atom_environment(atom_symbol='Si', ce_symbol='T:4'))
        self.assertFalse(lse.structure_contains_atom_environment(atom_symbol='O', ce_symbol='T:4'))
        self.assertTrue(lse.uniquely_determines_coordination_environments)
        self.assertFalse(lse.__ne__(lse))

        envs = lse.strategy.get_site_coordination_environments(lse.structure[6])
        self.assertEqual(len(envs), 1)
        self.assertEqual(envs[0][0], 'T:4')

        multi_strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters()

        lse_multi = LightStructureEnvironments.from_structure_environments(strategy=multi_strategy,
                                                                           structure_environments=se,
                                                                           valences='undefined')
        self.assertAlmostEqual(lse_multi.coordination_environments[isite][0]['csm'], 0.009887784240541068)
        self.assertAlmostEqual(lse_multi.coordination_environments[isite][0]['ce_fraction'], 1.0)
        self.assertEqual(lse_multi.coordination_environments[isite][0]['ce_symbol'], 'T:4')

    @classmethod
    def tearDownClass(cls):
        # Remove the directory in which the temporary files have been created
        shutil.rmtree('tmp_dir')


if __name__ == "__main__":
    unittest.main()