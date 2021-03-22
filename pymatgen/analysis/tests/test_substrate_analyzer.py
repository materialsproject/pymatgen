# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
TODO: Modify unittest doc.
"""

__author__ = "Shyam Dwaraknath"
__copyright__ = "Copyright 2016, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__date__ = "2/5/16"

import unittest

from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.substrate_analyzer import (
    SubstrateAnalyzer,
    ZSLGenerator,
    fast_norm,
    get_factors,
    reduce_vectors,
    vec_area,
)
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.testing import PymatgenTest


class ZSLGenTest(PymatgenTest):
    # Clean up test to be based on test structures

    def test_init(self):
        # Film VO2
        film = SpacegroupAnalyzer(self.get_structure("VO2"), symprec=0.1).get_conventional_standard_structure()

        # Substrate TiO2
        substrate = SpacegroupAnalyzer(self.get_structure("TiO2"), symprec=0.1).get_conventional_standard_structure()

        z = ZSLGenerator()

        self.assertAlmostEqual(fast_norm([3, 2, 1]), 3.7416573867739413)
        self.assertArrayEqual(reduce_vectors([1, 0, 0], [2, 2, 0]), [[1, 0, 0], [0, 2, 0]])
        self.assertEqual(vec_area([1, 0, 0], [0, 2, 0]), 2)
        self.assertArrayEqual(list(get_factors(18)), [1, 2, 3, 6, 9, 18])
        self.assertTrue(z.is_same_vectors([[1.01, 0, 0], [0, 2, 0]], [[1, 0, 0], [0, 2.01, 0]]))
        self.assertFalse(z.is_same_vectors([[1.01, 2, 0], [0, 2, 0]], [[1, 0, 0], [0, 2.01, 0]]))

        matches = list(z(film.lattice.matrix[:2], substrate.lattice.matrix[:2]))

        self.assertEqual(len(matches), 8)


class SubstrateAnalyzerTest(PymatgenTest):
    # Clean up test to be based on test structures

    def test_init(self):
        # Film VO2
        film = SpacegroupAnalyzer(self.get_structure("VO2"), symprec=0.1).get_conventional_standard_structure()

        # Substrate TiO2
        substrate = SpacegroupAnalyzer(self.get_structure("TiO2"), symprec=0.1).get_conventional_standard_structure()

        film_elac = ElasticTensor.from_voigt(
            [
                [324.32, 187.3, 170.92, 0.0, 0.0, 0.0],
                [187.3, 324.32, 170.92, 0.0, 0.0, 0.0],
                [170.92, 170.92, 408.41, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 150.73, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 150.73, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 238.74],
            ]
        )

        s = SubstrateAnalyzer()

        matches = list(s.calculate(film, substrate, film_elac))
        self.assertEqual(len(matches), 192)


if __name__ == "__main__":
    unittest.main()
