#!/usr/bin/env python


__author__ = 'waroquiers'

import unittest
import os
import json
import shutil
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import StructureEnvironments
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import MultiWeightsChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import AngleNbSetWeight
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import CNBiasNbSetWeight
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import DeltaCSMNbSetWeight
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import DistanceAngleAreaNbSetWeight
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import NormalizedAngleDistanceNbSetWeight
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SelfCSMNbSetWeight
from pymatgen.analysis.chemenv.coordination_environments.voronoi import DetailedVoronoiContainer
from pymatgen.core.structure import Structure
import numpy as np

json_files_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..",
                              'test_files', "chemenv", "json_test_files")
se_files_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..",
                              'test_files', "chemenv", "structure_environments_files")


class ReadWriteChemenvTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.lgf = LocalGeometryFinder()
        cls.lgf.setup_parameters(centering_type='standard')
        os.makedirs('tmp_dir')

    def test_read_write_structure_environments(self):
        f = open("{}/{}".format(json_files_dir, 'test_T--4_FePO4_icsd_4266.json'), 'r')
        dd = json.load(f)
        f.close()

        atom_indices = dd['atom_indices']

        struct = Structure.from_dict(dd['structure'])
        self.lgf.setup_structure(struct)
        se = self.lgf.compute_structure_environments(only_indices=atom_indices,
                                                     maximum_distance_factor=2.25,
                                                     get_from_hints=True)

        f = open('tmp_dir/se.json', 'w')
        json.dump(se.as_dict(), f)
        f.close()

        f = open('tmp_dir/se.json', 'r')
        dd = json.load(f)
        f.close()

        se2 = StructureEnvironments.from_dict(dd)

        self.assertEqual(se, se2)

        strategy = SimplestChemenvStrategy()
        lse = LightStructureEnvironments.from_structure_environments(structure_environments=se, strategy=strategy,
                                                                     valences='undefined')

        f = open('tmp_dir/lse.json', 'w')
        json.dump(lse.as_dict(), f)
        f.close()

        f = open('tmp_dir/lse.json', 'r')
        dd = json.load(f)
        f.close()

        lse2 = LightStructureEnvironments.from_dict(dd)

        self.assertEqual(lse, lse2)

    def test_structure_environments_neighbors_sets(self):
        f = open("{}/{}".format(se_files_dir, 'se_mp-7000.json'), 'r')
        dd = json.load(f)
        f.close()

        se = StructureEnvironments.from_dict(dd)

        isite = 6
        nb_set = se.neighbors_sets[isite][4][0]

        nb_set_surface_points = np.array([[1.0017922780870239, 0.99301365328679292],
                                          [1.0017922780870239, 0.0],
                                          [2.2237615554448569, 0.0],
                                          [2.2237615554448569, 0.0060837],
                                          [2.25, 0.0060837],
                                          [2.25, 0.99301365328679292]])

        self.assertTrue(np.allclose(np.array(nb_set.voronoi_grid_surface_points()), nb_set_surface_points))

        neighb_sites = nb_set.neighb_sites
        coords = [np.array([0.2443798, 1.80409653, -1.13218359]), np.array([1.44020353, 1.11368738, 1.13218359]),
                  np.array([2.75513098, 2.54465207, -0.70467298]), np.array([0.82616785, 3.65833945, 0.70467298])]

        np.testing.assert_array_almost_equal(coords[0], neighb_sites[0].coords)
        np.testing.assert_array_almost_equal(coords[1], neighb_sites[1].coords)
        np.testing.assert_array_almost_equal(coords[2], neighb_sites[2].coords)
        np.testing.assert_array_almost_equal(coords[3], neighb_sites[3].coords)

        neighb_coords = nb_set.coords

        np.testing.assert_array_almost_equal(coords, neighb_coords[1:])
        np.testing.assert_array_almost_equal(nb_set.structure[nb_set.isite].coords, neighb_coords[0])

        normdist = nb_set.normalized_distances
        self.assertAlmostEqual(sorted(normdist),
                               sorted([1.0017922783963027, 1.0017922780870239, 1.000000000503177, 1.0]))
        normang = nb_set.normalized_angles
        self.assertAlmostEqual(sorted(normang),
                               sorted([0.9999999998419052, 1.0, 0.9930136530585189, 0.9930136532867929]))
        dist = nb_set.distances
        self.assertAlmostEqual(sorted(dist),
                               sorted([1.6284399814843944, 1.6284399809816534, 1.6255265861208676, 1.6255265853029401]))
        ang = nb_set.angles
        self.assertAlmostEqual(sorted(ang),
                               sorted([3.117389876236432, 3.117389876729275, 3.095610709498583, 3.0956107102102024]))

        nb_set_info = nb_set.info

        self.assertAlmostEqual(nb_set_info['normalized_angles_mean'], 0.996506826547)
        self.assertAlmostEqual(nb_set_info['normalized_distances_std'], 0.000896138995037)
        self.assertAlmostEqual(nb_set_info['angles_std'], 0.0108895833142)
        self.assertAlmostEqual(nb_set_info['distances_std'], 0.00145669776056)
        self.assertAlmostEqual(nb_set_info['distances_mean'], 1.62698328347)

        self.assertEqual(nb_set.__str__(), 'Neighbors Set for site #6 :\n'
                                           ' - Coordination number : 4\n'
                                           ' - Voronoi indices : 1, 4, 5, 6\n')

        self.assertFalse(nb_set.__ne__(nb_set))

        self.assertEqual(nb_set.__hash__(), 4)

    def test_strategies(self):
        simplest_strategy_1 = SimplestChemenvStrategy()
        simplest_strategy_2 = SimplestChemenvStrategy(distance_cutoff=1.5, angle_cutoff=0.5)
        self.assertFalse(simplest_strategy_1 == simplest_strategy_2)
        simplest_strategy_1_from_dict = SimplestChemenvStrategy.from_dict(simplest_strategy_1.as_dict())
        self.assertTrue(simplest_strategy_1, simplest_strategy_1_from_dict)

        effective_csm_estimator = {'function': 'power2_inverse_decreasing',
                                   'options': {'max_csm': 8.0}}
        self_csm_weight = SelfCSMNbSetWeight()
        surface_definition = {'type': 'standard_elliptic',
                              'distance_bounds': {'lower': 1.1, 'upper': 1.9},
                              'angle_bounds': {'lower': 0.1, 'upper': 0.9}}
        surface_definition_2 = {'type': 'standard_elliptic',
                              'distance_bounds': {'lower': 1.1, 'upper': 1.9},
                              'angle_bounds': {'lower': 0.1, 'upper': 0.95}}
        da_area_weight = DistanceAngleAreaNbSetWeight(weight_type='has_intersection',
                                                      surface_definition=surface_definition,
                                                      nb_sets_from_hints='fallback_to_source',
                                                      other_nb_sets='0_weight',
                                                      additional_condition=DistanceAngleAreaNbSetWeight.AC.ONLY_ACB)
        da_area_weight_2 = DistanceAngleAreaNbSetWeight(weight_type='has_intersection',
                                                      surface_definition=surface_definition_2,
                                                      nb_sets_from_hints='fallback_to_source',
                                                      other_nb_sets='0_weight',
                                                      additional_condition=DistanceAngleAreaNbSetWeight.AC.ONLY_ACB)
        weight_estimator = {'function': 'smootherstep',
                            'options': {'delta_csm_min': 0.5,
                                        'delta_csm_max': 3.0}}
        symmetry_measure_type = 'csm_wcs_ctwcc'
        delta_weight = DeltaCSMNbSetWeight(effective_csm_estimator=effective_csm_estimator,
                                           weight_estimator=weight_estimator,
                                           symmetry_measure_type=symmetry_measure_type)
        bias_weight = CNBiasNbSetWeight.linearly_equidistant(weight_cn1=1.0, weight_cn13=4.0)
        bias_weight_2 = CNBiasNbSetWeight.linearly_equidistant(weight_cn1=1.0, weight_cn13=5.0)
        angle_weight = AngleNbSetWeight()
        nad_weight = NormalizedAngleDistanceNbSetWeight(average_type='geometric', aa=1, bb=1)
        multi_weights_strategy_1 = MultiWeightsChemenvStrategy(dist_ang_area_weight=da_area_weight,
                                                               self_csm_weight=self_csm_weight,
                                                               delta_csm_weight=delta_weight,
                                                               cn_bias_weight=bias_weight,
                                                               angle_weight=angle_weight,
                                                               normalized_angle_distance_weight=nad_weight,
                                                               symmetry_measure_type=symmetry_measure_type)
        multi_weights_strategy_2 = MultiWeightsChemenvStrategy(dist_ang_area_weight=da_area_weight,
                                                               self_csm_weight=self_csm_weight,
                                                               delta_csm_weight=delta_weight,
                                                               cn_bias_weight=bias_weight_2,
                                                               angle_weight=angle_weight,
                                                               normalized_angle_distance_weight=nad_weight,
                                                               symmetry_measure_type=symmetry_measure_type)
        multi_weights_strategy_3 = MultiWeightsChemenvStrategy(dist_ang_area_weight=da_area_weight_2,
                                                               self_csm_weight=self_csm_weight,
                                                               delta_csm_weight=delta_weight,
                                                               cn_bias_weight=bias_weight,
                                                               angle_weight=angle_weight,
                                                               normalized_angle_distance_weight=nad_weight,
                                                               symmetry_measure_type=symmetry_measure_type)
        multi_weights_strategy_1_from_dict = MultiWeightsChemenvStrategy.from_dict(multi_weights_strategy_1.as_dict())

        self.assertTrue(multi_weights_strategy_1 == multi_weights_strategy_1_from_dict)
        self.assertFalse(simplest_strategy_1 == multi_weights_strategy_1)
        self.assertFalse(multi_weights_strategy_1 == multi_weights_strategy_2)
        self.assertFalse(multi_weights_strategy_1 == multi_weights_strategy_3)
        self.assertFalse(multi_weights_strategy_2 == multi_weights_strategy_3)

    def test_read_write_voronoi(self):
        f = open("{}/{}".format(json_files_dir, 'test_T--4_FePO4_icsd_4266.json'), 'r')
        dd = json.load(f)
        f.close()

        struct = Structure.from_dict(dd['structure'])

        valences = [site.specie.oxi_state for site in struct]

        detailed_voronoi_container = DetailedVoronoiContainer(structure=struct, valences=valences)

        f = open('tmp_dir/se.json', 'w')
        json.dump(detailed_voronoi_container.as_dict(), f)
        f.close()

        f = open('tmp_dir/se.json', 'r')
        dd = json.load(f)
        f.close()

        detailed_voronoi_container2 = DetailedVoronoiContainer.from_dict(dd)

        self.assertEqual(detailed_voronoi_container, detailed_voronoi_container2)

    @classmethod
    def tearDownClass(cls):
        #Remove the directory in which the temporary files have been created
        shutil.rmtree('tmp_dir')


if __name__ == "__main__":
    unittest.main()
