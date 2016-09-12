#!/usr/bin/env python


__author__ = 'waroquiers'

import unittest2
import os
import json
import numpy as np
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
# from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy
# from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimpleAbundanceChemenvStrategy
# from pymatgen.core.structure import Structure


json_files_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..",
                              'test_files', "chemenv", "json_test_files")


class CoordinationGeometryFinderTest(unittest2.TestCase):

    def setUp(self):
        self.lgf = LocalGeometryFinder()
        self.lgf.setup_parameters(centering_type='standard',
                                  structure_refinement=self.lgf.STRUCTURE_REFINEMENT_NONE)
    #     self.strategies = [SimplestChemenvStrategy(), SimpleAbundanceChemenvStrategy()]
    #
    # def _strategy_test(self, strategy):
    #     files = []
    #     for (dirpath, dirnames, filenames) in os.walk(json_files_dir):
    #         files.extend(filenames)
    #         break
    #
    #     for ifile, json_file in enumerate(files):
    #         with self.subTest(json_file=json_file):
    #             f = open("{}/{}".format(json_files_dir, json_file), 'r')
    #             dd = json.load(f)
    #             f.close()
    #
    #             atom_indices = dd['atom_indices']
    #             expected_geoms = dd['expected_geoms']
    #
    #             struct = Structure.from_dict(dd['structure'])
    #
    #             struct = self.lgf.setup_structure(struct)
    #             se = self.lgf.compute_structure_environments_detailed_voronoi(only_indices=atom_indices,
    #                                                                           maximum_distance_factor=1.5)
    #
    #             #All strategies should get the correct environment with their default parameters
    #             strategy.set_structure_environments(se)
    #             for ienv, isite in enumerate(atom_indices):
    #                 ce = strategy.get_site_coordination_environment(struct[isite])
    #                 try:
    #                     coord_env = ce[0]
    #                 except TypeError:
    #                     coord_env = ce
    #                 #Check that the environment found is the expected one
    #                 self.assertEqual(coord_env, expected_geoms[ienv])
    #
    # def test_simplest_chemenv_strategy(self):
    #     strategy = SimplestChemenvStrategy()
    #     self._strategy_test(strategy)
    #
    # def test_simple_abundance_chemenv_strategy(self):
    #     strategy = SimpleAbundanceChemenvStrategy()
    #     self._strategy_test(strategy)

    def test_perfect_environments(self):
        allcg = AllCoordinationGeometries()
        indices_CN = {1: [0],
                      2: [1, 0],
                      3: [1, 0, 2],
                      4: [2, 0, 3, 1],
                      5: [2, 3, 1, 0, 4],
                      6: [0, 2, 3, 1, 5, 4],
                      7: [2, 6, 0, 3, 4, 5, 1],
                      8: [1, 2, 6, 3, 7, 0, 4, 5],
                      9: [5, 2, 6, 0, 4, 7, 3, 8, 1],
                      10: [8, 5, 6, 3, 0, 7, 2, 4, 9, 1],
                      11: [7, 6, 4, 1, 2, 5, 0, 8, 9, 10, 3],
                      12: [5, 8, 9, 0, 3, 1, 4, 2, 6, 11, 10, 7],
                      13: [4, 11, 5, 12, 1, 2, 8, 3, 0, 6, 9, 7, 10],
                      }

        for coordination in range(1, 14):
            for mp_symbol in allcg.get_implemented_geometries(coordination=coordination,
                                                              returned='mp_symbol'):
                with self.subTest(msg=mp_symbol, mp_symbol=mp_symbol):
                    cg = allcg.get_geometry_from_mp_symbol(mp_symbol=mp_symbol)
                    self.lgf.allcg = AllCoordinationGeometries(only_symbols=[mp_symbol])
                    self.lgf.setup_test_perfect_environment(mp_symbol, randomness=False,
                                                            indices=indices_CN[coordination],
                                                            random_translation='NONE', random_rotation='NONE',
                                                            random_scale='NONE')
                    se = self.lgf.compute_structure_environments(only_indices=[0],
                                                                 maximum_distance_factor=1.01*cg.distfactor_max,
                                                                 min_cn=cg.coordination_number,
                                                                 max_cn=cg.coordination_number,
                                                                 only_symbols=[mp_symbol]
                                                                 )
                    self.assertAlmostEqual(se.get_csm(0, mp_symbol)['symmetry_measure'], 0.0, delta=1e-8,
                                           msg='Failed to get perfect environment with mp_symbol {}'.format(mp_symbol))


if __name__ == "__main__":
    unittest2.main()
