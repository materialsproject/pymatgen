# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest

from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.interfaces.substrate_analyzer import SubstrateAnalyzer
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.testing import PymatgenTest


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
        for match in matches:
            assert match is not None
            assert isinstance(match.match_area, float)


if __name__ == "__main__":
    unittest.main()
