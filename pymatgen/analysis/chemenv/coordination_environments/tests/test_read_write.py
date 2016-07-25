#!/usr/bin/env python


__author__ = 'waroquiers'

import unittest2 as unittest
import os
import json
import shutil
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import StructureEnvironments
from pymatgen.analysis.chemenv.coordination_environments.voronoi import DetailedVoronoiContainer
from pymatgen.core.structure import Structure

json_files_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..",
                              'test_files', "chemenv", "json_test_files")


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
                                                     maximum_distance_factor=2.25)

        f = open('tmp_dir/se.json', 'w')
        json.dump(se.as_dict(), f)
        f.close()

        f = open('tmp_dir/se.json', 'r')
        dd = json.load(f)
        f.close()

        se2 = StructureEnvironments.from_dict(dd)

        self.assertEqual(se, se2)

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