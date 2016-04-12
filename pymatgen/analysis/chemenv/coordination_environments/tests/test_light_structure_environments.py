#!/usr/bin/env python


__author__ = 'waroquiers'

import unittest2 as unittest
import os
import json
import shutil
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import StructureEnvironments
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimpleAbundanceChemenvStrategy
from pymatgen.core.structure import Structure


json_files_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..",
                              'test_files', "chemenv", "json_test_files")


class LightStructureEnvironmentsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        os.makedirs('tmp_dir')

    def setUp(self):
        self.lgf = LocalGeometryFinder()
        self.lgf.setup_parameters(centering_type='standard')
        self.strategies = [SimplestChemenvStrategy(), SimpleAbundanceChemenvStrategy()]

    def test_read_structure_environments(self):
        f = open("{}/{}".format(json_files_dir, 'test_T--4_FePO4_icsd_4266.json'), 'r')
        dd = json.load(f)
        f.close()

        atom_indices = dd['atom_indices']

        struct = Structure.from_dict(dd['structure'])
        self.lgf.setup_structure(struct)
        se = self.lgf.compute_structure_environments_detailed_voronoi(only_indices=atom_indices,
                                                                      maximum_distance_factor=2.25)

        f = open('tmp_dir/se.json', 'w')
        json.dump(se.as_dict(), f)
        f.close()

        f = open('tmp_dir/se.json', 'r')
        dd = json.load(f)
        f.close()

        se2 = StructureEnvironments.from_dict(dd)

        self.assertEqual(se, se2)

        _strategy = SimplestChemenvStrategy()
        light_se = LightStructureEnvironments(_strategy, se)

        f = open('tmp_dir/light_se.json', 'w')
        json.dump(light_se.as_dict(), f)
        f.close()

        f = open('tmp_dir/light_se.json', 'r')
        dd = json.load(f)
        f.close()

        light_se2 = LightStructureEnvironments.from_dict(dd)

        self.assertEqual(light_se, light_se2)

    @classmethod
    def tearDownClass(cls):
        #Remove the directory in which the temporary files have been created
        shutil.rmtree('tmp_dir')


if __name__ == "__main__":
    unittest.main()