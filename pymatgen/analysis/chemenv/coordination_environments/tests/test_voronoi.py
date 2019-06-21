#!/usr/bin/env python


__author__ = 'waroquiers'

import unittest

from monty.tempfile import ScratchDir
from pymatgen.analysis.chemenv.coordination_environments.voronoi import DetailedVoronoiContainer
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
import os
import shutil
import random
import numpy as np
from pymatgen.util.testing import PymatgenTest

json_files_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..",
                              'test_files', "chemenv", "json_test_files")
img_files_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..",
                              'test_files', "chemenv", "images")



class VoronoiContainerTest(PymatgenTest):

    def test_voronoi(self):
        with ScratchDir("."):
            #Define a cubic lattice and a list of species (to be used for the fake structures)
            cubic_lattice = Lattice.cubic(10.0)
            species = ['Cu', 'O', 'O', 'O', 'O', 'O', 'O']
            valences = 'undefined'

            #First fake structure
            coords = [[5.0, 5.0, 5.0]]
            order_and_coords = [(1, [4.0, 5.0, 5.0]),
                                (2, [6.01, 5.0, 5.0]),
                                (3, [5.0, 3.98, 5.0]),
                                (4, [5.0, 6.03, 5.0]),
                                (5, [5.0, 5.0, 3.96]),
                                (6, [5.0, 5.0, 6.05])]
            random.shuffle(order_and_coords)
            sorted = np.argsort([oc[0] for oc in order_and_coords]) + 1
            coords.extend([oc[1] for oc in order_and_coords])
            fake_structure = Structure(cubic_lattice, species, coords, coords_are_cartesian=True)

            #First fake structure with a given normalized_distance_tolerance of 0.0100001
            detailed_voronoi_container = DetailedVoronoiContainer(structure=fake_structure, valences=valences,
                                                                  normalized_distance_tolerance=0.0100001, isites=[0])
            self.assertEqual(len(detailed_voronoi_container.voronoi_list2[0]), 6)
            neighbors = detailed_voronoi_container.neighbors(0, 1.0, 0.5, True)
            self.assertEqual(len(neighbors), 6)
            neighbors = detailed_voronoi_container.neighbors(0, 1.02, 0.5, True)
            self.assertEqual(len(neighbors), 6)
            neighbors = detailed_voronoi_container.neighbors(0, 1.026, 0.5, True)
            self.assertEqual(len(neighbors), 6)
            neighbors = detailed_voronoi_container.neighbors(0, 1.5, 0.5, True)
            self.assertEqual(len(neighbors), 6)

            #First fake structure with a given normalized_distance_tolerance of 0.001
            detailed_voronoi_container = DetailedVoronoiContainer(structure=fake_structure, valences=valences,
                                                                  normalized_distance_tolerance=0.001, isites=[0])
            self.assertEqual(len(detailed_voronoi_container.voronoi_list2[0]), 6)
            neighbors = detailed_voronoi_container.neighbors(0, 1.0, 0.5, True)
            self.assertEqual(len(neighbors), 1)
            self.assertEqual(neighbors[0]['site'], fake_structure[sorted[0]])
            neighbors = detailed_voronoi_container.neighbors(0, 1.02, 0.5, True)
            nbs = [nb['site'] for nb in neighbors]
            self.assertEqual(len(neighbors), 3)
            self.assertTrue(fake_structure[sorted[0]] in nbs)
            self.assertTrue(fake_structure[sorted[1]] in nbs)
            self.assertTrue(fake_structure[sorted[2]] in nbs)
            neighbors = detailed_voronoi_container.neighbors(0, 1.026, 0.5, True)
            nbs = [nb['site'] for nb in neighbors]
            self.assertEqual(len(neighbors), 3)
            self.assertTrue(fake_structure[sorted[0]] in nbs)
            self.assertTrue(fake_structure[sorted[1]] in nbs)
            self.assertTrue(fake_structure[sorted[2]] in nbs)
            neighbors = detailed_voronoi_container.neighbors(0, 1.5, 0.5, True)
            self.assertEqual(len(neighbors), 6)

            #Second fake structure
            coords2 = [[5.0, 5.0, 5.0]]
            order_and_coords = [(1, [4.0, 5.0, 5.0]),
                                (2, [6.01, 5.0, 5.0]),
                                (3, [5.0, 3.98, 5.0]),
                                (4, [5.0, 6.07, 5.0]),
                                (5, [5.0, 5.0, 3.92]),
                                (6, [5.0, 5.0, 6.09])]
            random.shuffle(order_and_coords)
            sorted = np.argsort([oc[0] for oc in order_and_coords]) + 1
            coords2.extend([oc[1] for oc in order_and_coords])
            fake_structure2 = Structure(cubic_lattice, species, coords2, coords_are_cartesian=True)

            #Second fake structure with a given normalized_distance_tolerance of 0.0100001
            detailed_voronoi_container = DetailedVoronoiContainer(structure=fake_structure2, valences=valences,
                                                                  normalized_distance_tolerance=0.0100001, isites=[0])
            self.assertEqual(len(detailed_voronoi_container.voronoi_list2[0]), 6)
            neighbors = detailed_voronoi_container.neighbors(0, 1.0, 0.5, True)
            nbs = [nb['site'] for nb in neighbors]
            self.assertEqual(len(neighbors), 3)
            self.assertTrue(fake_structure2[sorted[0]] in nbs)
            self.assertTrue(fake_structure2[sorted[1]] in nbs)
            self.assertTrue(fake_structure2[sorted[2]] in nbs)
            neighbors = detailed_voronoi_container.neighbors(0, 1.02, 0.5, True)
            nbs = [nb['site'] for nb in neighbors]
            self.assertEqual(len(neighbors), 3)
            self.assertTrue(fake_structure2[sorted[0]] in nbs)
            self.assertTrue(fake_structure2[sorted[1]] in nbs)
            self.assertTrue(fake_structure2[sorted[2]] in nbs)
            neighbors = detailed_voronoi_container.neighbors(0, 1.026, 0.5, True)
            nbs = [nb['site'] for nb in neighbors]
            self.assertEqual(len(neighbors), 3)
            self.assertTrue(fake_structure2[sorted[0]] in nbs)
            self.assertTrue(fake_structure2[sorted[1]] in nbs)
            self.assertTrue(fake_structure2[sorted[2]] in nbs)
            neighbors = detailed_voronoi_container.neighbors(0, 1.5, 0.5, True)
            self.assertEqual(len(neighbors), 6)

            species = ['Cu', 'Cu', 'O', 'O', 'O', 'Cu', 'O']
            valences = [2, 2, -2, -2, -2, 2, -2]

            #Third fake structure (test of the only_anion_cation_bonds)
            coords = [[5.0, 5.0, 5.0],
                      [6.01, 5.0, 5.0],
                      [5.0, 5.0, 3.96],
                      [4.0, 5.0, 5.0],
                      [5.0, 6.03, 5.0],
                      [5.0, 3.98, 5.0],
                      [5.0, 5.0, 6.05]]
            fake_structure3 = Structure(cubic_lattice, species, coords, coords_are_cartesian=True)
            detailed_voronoi_container = DetailedVoronoiContainer(structure=fake_structure3, valences=valences,
                                                                  normalized_distance_tolerance=0.0100001, isites=[0],
                                                                  additional_conditions=
                                                                  [DetailedVoronoiContainer.AC.ONLY_ACB])
            self.assertEqual(len(detailed_voronoi_container.voronoi_list2[0]), 6)
            neighbors = detailed_voronoi_container.neighbors(0, 1.01, 0.5, True)
            nbs = [nb['site'] for nb in neighbors]
            self.assertEqual(len(neighbors), 6)
            self.assertTrue(fake_structure3[1] in nbs)
            self.assertTrue(fake_structure3[2] in nbs)
            self.assertTrue(fake_structure3[3] in nbs)
            self.assertTrue(fake_structure3[4] in nbs)
            self.assertTrue(fake_structure3[5] in nbs)
            self.assertTrue(fake_structure3[6] in nbs)

            #Test of the as_dict() and from_dict() methods as well as __eq__ method
            other_detailed_voronoi_container = DetailedVoronoiContainer.from_dict(detailed_voronoi_container.as_dict())
            self.assertTrue(detailed_voronoi_container, other_detailed_voronoi_container)

    def test_get_vertices_dist_ang_indices(self):
        with ScratchDir("."):
            cubic_lattice = Lattice.cubic(10.0)
            species = ['Cu', 'O', 'O', 'O', 'O', 'O', 'O']
            valences = 'undefined'

            #First fake structure
            coords = [[5.0, 5.0, 5.0],
                      [6.01, 5.0, 5.0],
                      [5.0, 5.0, 3.96],
                      [4.0, 5.0, 5.0],
                      [5.0, 6.03, 5.0],
                      [5.0, 3.98, 5.0],
                      [5.0, 5.0, 6.05]]
            fake_structure = Structure(cubic_lattice, species, coords, coords_are_cartesian=True)

            #First fake structure with a given normalized_distance_tolerance of 0.0100001
            detailed_voronoi_container = DetailedVoronoiContainer(structure=fake_structure, valences=valences,
                                                                  normalized_distance_tolerance=0.0100001, isites=[0])
            fake_parameter_indices_list = []
            for ii in range(2, 5):
                for jj in range(7, 14):
                    fake_parameter_indices_list.append((ii, jj))
            for ii in range(5, 7):
                for jj in range(10, 14):
                    fake_parameter_indices_list.append((ii, jj))

            points = detailed_voronoi_container._get_vertices_dist_ang_indices(fake_parameter_indices_list)
            self.assertEqual(points[0], (2, 7))
            self.assertEqual(points[1], (4, 7))
            self.assertEqual(points[2], (4, 10))
            self.assertEqual(points[3], (6, 10))
            self.assertEqual(points[4], (6, 13))
            self.assertEqual(points[5], (2, 13))

if __name__ == "__main__":
    unittest.main()