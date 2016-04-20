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
from pymatgen.core.structure import Structure

json_files_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..",
                              'test_files', "chemenv", "json_test_files")


class ReadWriteChemenvTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.lgf = LocalGeometryFinder()
        cls.lgf.setup_parameters(centering_type='standard')
        os.makedirs('tmp_dir')

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

        self.assertEqual(light_se._strategy, light_se2._strategy)
        self.assertEqual(light_se._structure, light_se2._structure)
        self.assertEqual(light_se._bva_valences, light_se2._bva_valences)
        self.assertEqual(light_se._coordination_environments, light_se2._coordination_environments)
        self.assertEqual(light_se._neighbors, light_se2._neighbors)
        self.assertEqual(light_se, light_se2)

    @classmethod
    def tearDownClass(cls):
        #Remove the directory in which the temporary files have been created
        shutil.rmtree('tmp_dir')


if __name__ == "__main__":
    unittest.main()