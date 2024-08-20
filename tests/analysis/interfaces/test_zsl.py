from __future__ import annotations

import numpy as np
from numpy.testing import assert_array_equal
from pytest import approx

from pymatgen.analysis.interfaces.zsl import (
    ZSLGenerator,
    fast_norm,
    get_factors,
    is_same_vectors,
    reduce_vectors,
    vec_area,
)
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.testing import PymatgenTest

__author__ = "Shyam Dwaraknath"
__copyright__ = "Copyright 2016, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__date__ = "2/5/16"


class TestZSLGen(PymatgenTest):
    def setUp(self):
        # Film VO2
        self.film = SpacegroupAnalyzer(self.get_structure("VO2"), symprec=0.1).get_conventional_standard_structure()

        # Substrate TiO2
        self.substrate = SpacegroupAnalyzer(
            self.get_structure("TiO2"), symprec=0.1
        ).get_conventional_standard_structure()

    def test_zsl(self):
        zsl_gen = ZSLGenerator()

        assert fast_norm(np.array([3, 2, 1])) == approx(3.74165738)
        assert_array_equal(reduce_vectors(np.array([1, 0, 0]), np.array([2, 2, 0])), [[1, 0, 0], [0, 2, 0]])
        assert vec_area(np.array([1, 0, 0]), np.array([0, 2, 0])) == 2
        assert list(get_factors(18)) == [1, 2, 3, 6, 9, 18]
        assert is_same_vectors(
            np.array([[1.01, 0, 0], [0, 2, 0]], dtype=float), np.array([[1, 0, 0], [0, 2.01, 0]], dtype=float)
        )
        assert not is_same_vectors(
            np.array([[1.01, 2, 0], [0, 2, 0]], dtype=float), np.array([[1, 0, 0], [0, 2.01, 0]], dtype=float)
        )

        matches = list(zsl_gen(self.film.lattice.matrix[:2], self.substrate.lattice.matrix[:2]))
        assert len(matches) == 8

    def test_bidirectional(self):
        z = ZSLGenerator(max_area_ratio_tol=0.05, max_angle_tol=0.05, max_length_tol=0.05)

        matches = list(z(self.film.lattice.matrix[:2], self.substrate.lattice.matrix[:2]))
        assert len(matches) == 60

        matches = list(z(self.substrate.lattice.matrix[:2], self.film.lattice.matrix[:2]))
        assert len(matches) == 52

        z.bidirectional = True
        matches = list(z(self.substrate.lattice.matrix[:2], self.film.lattice.matrix[:2]))
        assert len(matches) == 60

        for match in matches:
            assert match is not None
            assert isinstance(match.match_area, float)
