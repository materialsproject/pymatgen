# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest
import os

from pymatgen.analysis.ewald import EwaldSummation, EwaldMinimizer
from pymatgen.io.vasp.inputs import Poscar
import numpy as np

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')

class EwaldSummationTest(unittest.TestCase):

    def test_init(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        p = Poscar.from_file(filepath)
        original_s = p.structure
        s = original_s.copy()
        s.add_oxidation_state_by_element({"Li": 1, "Fe": 2,
                                          "P": 5, "O": -2})
        ham = EwaldSummation(s)
        self.assertAlmostEqual(ham.real_space_energy, -354.91294268, 4,
                               "Real space energy incorrect!")
        self.assertAlmostEqual(ham.reciprocal_space_energy, 25.475754801, 4)
        self.assertAlmostEqual(ham.point_energy, -790.463835033, 4,
                               "Point space energy incorrect!")
        self.assertAlmostEqual(ham.total_energy, -1119.90102291, 2,
                               "Total space energy incorrect!")
        self.assertAlmostEqual(ham.forces[0,0], -1.98818620e-01, 4,
                               "Forces incorrect")
        self.assertAlmostEqual(sum(sum(abs(ham.forces))), 915.925354346, 4,
                               "Forces incorrect")
        self.assertAlmostEqual(sum(sum(ham.real_space_energy_matrix)),
                               - 354.91294268, 4,
                               "Real space energy matrix incorrect!")
        self.assertAlmostEqual(sum(sum(ham.reciprocal_space_energy_matrix)),
                               25.475754801, 4,
                               "Reciprocal space energy matrix incorrect!")
        self.assertAlmostEqual(sum(ham.point_energy_matrix), -790.463835033,
                               4, "Point space energy matrix incorrect!")
        self.assertAlmostEqual(sum(sum(ham.total_energy_matrix)),
                               - 1119.90102291, 2,
                               "Total space energy matrix incorrect!")
        #note that forces are not individually tested, but should work fine.

        self.assertRaises(ValueError, EwaldSummation, original_s)
        #try sites with charge.
        charges = []
        for site in original_s:
            if site.specie.symbol == "Li":
                charges.append(1)
            elif site.specie.symbol == "Fe":
                charges.append(2)
            elif site.specie.symbol == "P":
                charges.append(5)
            else:
                charges.append(-2)

        original_s.add_site_property('charge', charges)
        ham2 = EwaldSummation(original_s)
        self.assertAlmostEqual(ham2.real_space_energy, -354.91294268, 4,
                               "Real space energy incorrect!")


class EwaldMinimizerTest(unittest.TestCase):

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

if __name__ == "__main__":
    unittest.main()
