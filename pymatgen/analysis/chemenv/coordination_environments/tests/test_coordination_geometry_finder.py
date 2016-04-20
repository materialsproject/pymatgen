#!/usr/bin/env python


__author__ = 'waroquiers'

import unittest2
import os
import json
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimpleAbundanceChemenvStrategy
from pymatgen.core.structure import Structure


json_files_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..",
                              'test_files', "chemenv", "json_test_files")


class CoordinationGeometryFinderTest(unittest2.TestCase):

    def setUp(self):
        self.lgf = LocalGeometryFinder()
        self.lgf.setup_parameters(centering_type='standard')
        self.strategies = [SimplestChemenvStrategy(), SimpleAbundanceChemenvStrategy()]

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
        for coordination in range(1, 13):
            for mp_symbol in AllCoordinationGeometries().get_implemented_geometries(coordination=coordination,
                                                                                    returned='mp_symbol'):
                with self.subTest(msg=mp_symbol, mp_symbol=mp_symbol):
                    self.lgf.setup_test_perfect_environment(mp_symbol, randomness=False,
                                                            indices='RANDOM',
                                                            random_translation=True, random_rotation=True,
                                                            random_scale=True)
                    se = self.lgf.compute_structure_environments_detailed_voronoi(only_indices=[0],
                                                                                  maximum_distance_factor=1.5)
                    self.assertAlmostEqual(se.get_csm(0, mp_symbol)['symmetry_measure'], 0.0, 4)

if __name__ == "__main__":
    unittest2.main(verbosity=9)
