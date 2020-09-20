# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import unittest
import os
import warnings

from pymatgen.analysis.ewald import EwaldSummation, EwaldMinimizer
from pymatgen.io.vasp.inputs import Poscar
import numpy as np

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class EwaldSummationTest(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore")
        filepath = os.path.join(test_dir, 'POSCAR')
        p = Poscar.from_file(filepath, check_for_POTCAR=False)
        self.original_s = p.structure
        self.s = self.original_s.copy()
        self.s.add_oxidation_state_by_element({"Li": 1, "Fe": 2,
                                               "P": 5, "O": -2})

    def tearDown(self):
        warnings.simplefilter("default")

    def test_init(self):
        ham = EwaldSummation(self.s, compute_forces=True)
        self.assertAlmostEqual(ham.real_space_energy, -502.23549897772602, 4)
        self.assertAlmostEqual(ham.reciprocal_space_energy, 6.1541071599534654, 4)
        self.assertAlmostEqual(ham.point_energy, -620.22598358035918, 4)
        self.assertAlmostEqual(ham.total_energy, -1123.00766, 1)
        self.assertAlmostEqual(ham.forces[0, 0], -1.98818620e-01, 4)
        self.assertAlmostEqual(sum(sum(abs(ham.forces))), 915.925354346, 4,
                               "Forces incorrect")
        self.assertAlmostEqual(sum(sum(ham.real_space_energy_matrix)),
                               ham.real_space_energy, 4)
        self.assertAlmostEqual(sum(sum(ham.reciprocal_space_energy_matrix)),
                               ham.reciprocal_space_energy, 4)
        self.assertAlmostEqual(sum(ham.point_energy_matrix), ham.point_energy,
                               4)
        self.assertAlmostEqual(sum(sum(ham.total_energy_matrix)) + ham._charged_cell_energy,
                               ham.total_energy, 2)

        self.assertRaises(ValueError, EwaldSummation, self.original_s)
        # try sites with charge.
        charges = []
        for site in self.original_s:
            if site.specie.symbol == "Li":
                charges.append(1)
            elif site.specie.symbol == "Fe":
                charges.append(2)
            elif site.specie.symbol == "P":
                charges.append(5)
            else:
                charges.append(-2)

        self.original_s.add_site_property('charge', charges)
        ham2 = EwaldSummation(self.original_s)
        self.assertAlmostEqual(ham2.real_space_energy, -502.23549897772602, 4)

    def test_from_dict(self):
        ham = EwaldSummation(self.s, compute_forces=True)
        ham2 = EwaldSummation.from_dict(ham.as_dict())
        self.assertIsNone(ham._real)
        self.assertFalse(ham._initialized)
        self.assertIsNone(ham2._real)
        self.assertFalse(ham2._initialized)
        self.assertTrue(np.array_equal(ham.total_energy_matrix,
                                       ham2.total_energy_matrix))
        # check lazy eval
        self.assertAlmostEqual(ham.total_energy, -1123.00766, 1)
        self.assertIsNotNone(ham._real)
        self.assertTrue(ham._initialized)
        ham2 = EwaldSummation.from_dict(ham.as_dict())
        self.assertIsNotNone(ham2._real)
        self.assertTrue(ham2._initialized)
        self.assertTrue(np.array_equal(ham.total_energy_matrix,
                                       ham2.total_energy_matrix))

    def test_as_dict(self):
        ham = EwaldSummation(self.s, compute_forces=True)
        d = ham.as_dict()
        self.assertTrue(d["compute_forces"])
        self.assertEqual(d["eta"], ham._eta)
        self.assertEqual(d["acc_factor"], ham._acc_factor)
        self.assertEqual(d["real_space_cut"], ham._rmax)
        self.assertEqual(d["recip_space_cut"], ham._gmax)
        self.assertEqual(ham.as_dict(),
                         EwaldSummation.from_dict(d).as_dict())


class EwaldMinimizerTest(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_init(self):
        matrix = np.array([[-3., 3., 4., -0., 3., 3., 1., 14., 9., -4.],
                           [1., -3., -3., 12., -4., -1., 5., 11., 1., 12.],
                           [14., 7., 13., 15., 13., 5., -5., 10., 14., -2.],
                           [9., 13., 4., 1., 3., -4., 7., 0., 6., -4.],
                           [4., -4., 6., 1., 12., -4., -2., 13., 0., 6.],
                           [13., 7., -4., 12., -2., 9., 8., -5., 3., 1.],
                           [8., 1., 10., -4., -2., 4., 13., 12., -3., 13.],
                           [2., 11., 8., 1., -1., 5., -3., 4., 5., 0.],
                           [-0., 14., 4., 3., -1., -5., 7., -1., -1., 3.],
                           [2., -2., 10., 1., 6., -5., -3., 12., 0., 13.]])

        m_list = [[.9, 4, [1, 2, 3, 4, 8], 'a'], [-1, 2, [5, 6, 7], 'b']]

        e_min = EwaldMinimizer(matrix, m_list, 50)

        self.assertEqual(len(e_min.output_lists), 15,
                         "Wrong number of permutations returned")
        self.assertAlmostEqual(e_min.minimized_sum, 111.63, 3,
                               "Returned wrong minimum value")
        self.assertEqual(len(e_min.best_m_list), 6,
                         "Returned wrong number of permutations")

    def test_site(self):
        """Test that uses an uncharged structure"""
        filepath = os.path.join(test_dir, 'POSCAR')
        p = Poscar.from_file(filepath, check_for_POTCAR=False)
        original_s = p.structure
        s = original_s.copy()
        s.add_oxidation_state_by_element({"Li": 1, "Fe": 3,
                                          "P": 5, "O": -2})

        # Comparison to LAMMPS result
        ham = EwaldSummation(s, compute_forces=True)
        self.assertAlmostEquals(-1226.3335, ham.total_energy, 3)
        self.assertAlmostEquals(-45.8338, ham.get_site_energy(0), 3)
        self.assertAlmostEquals(-27.2978, ham.get_site_energy(8), 3)


if __name__ == "__main__":
    unittest.main()
