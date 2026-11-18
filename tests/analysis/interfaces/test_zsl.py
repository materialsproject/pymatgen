from __future__ import annotations

import pickle
from itertools import product

import numpy as np
import pytest
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
from pymatgen.util.testing import MatSciTest

__author__ = "Shyam Dwaraknath"
__copyright__ = "Copyright 2016, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__date__ = "2/5/16"


class TestZSLGen(MatSciTest):
    def setup_method(self):
        # Film VO2
        self.film = SpacegroupAnalyzer(self.get_structure("VO2"), symprec=0.1).get_conventional_standard_structure()

        # Substrate TiO2
        self.substrate = SpacegroupAnalyzer(
            self.get_structure("TiO2"), symprec=0.1
        ).get_conventional_standard_structure()

    def test_zsl(self):
        zsl_gen = ZSLGenerator()

        assert fast_norm(np.array([3, 2, 1])) == approx(3.74165738)
        assert_array_equal(
            reduce_vectors(np.array([1, 0, 0]), np.array([2, 2, 0])),
            [[1, 0, 0], [0, 2, 0]],
        )
        assert vec_area(np.array([1, 0, 0]), np.array([0, 2, 0])) == 2
        assert list(get_factors(18)) == [1, 2, 3, 6, 9, 18]
        assert is_same_vectors(
            np.array([[1.01, 0, 0], [0, 2, 0]], dtype=float),
            np.array([[1, 0, 0], [0, 2.01, 0]], dtype=float),
        )
        assert not is_same_vectors(
            np.array([[1.01, 2, 0], [0, 2, 0]], dtype=float),
            np.array([[1, 0, 0], [0, 2.01, 0]], dtype=float),
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

    def test_get_equiv_transformations_matches_pairwise_check(self):
        film_vectors = self.film.lattice.matrix[:2]
        substrate_vectors = self.substrate.lattice.matrix[:2]
        for bidirectional in (False, True):
            z = ZSLGenerator(
                max_area_ratio_tol=0.05, max_angle_tol=0.05, max_length_tol=0.05, bidirectional=bidirectional
            )
            film_area, substrate_area = vec_area(*film_vectors), vec_area(*substrate_vectors)

            # Reference: reduce every super lattice and check each film/substrate pair one by one
            expected = []
            for film_trans, substrate_trans in z.generate_sl_transformation_sets(film_area, substrate_area):
                films = [reduce_vectors(*v) for v in np.dot(film_trans, film_vectors)]
                substrates = [reduce_vectors(*v) for v in np.dot(substrate_trans, substrate_vectors)]
                for (f_trans, s_trans), (f, s) in zip(
                    product(film_trans, substrate_trans), product(films, substrates), strict=True
                ):
                    if is_same_vectors(
                        np.array(f),
                        np.array(s),
                        bidirectional=bidirectional,
                        max_length_tol=z.max_length_tol,
                        max_angle_tol=z.max_angle_tol,
                    ):
                        expected.append([f, s, f_trans, s_trans])

            for _ in range(2):  # second pass reuses cached reductions
                transformation_sets = z.generate_sl_transformation_sets(film_area, substrate_area)
                matches = list(z.get_equiv_transformations(transformation_sets, film_vectors, substrate_vectors))
                assert len(matches) == len(expected) > 0
                for match, ref in zip(matches, expected, strict=True):
                    for arr, ref_arr in zip(match, ref, strict=True):
                        assert_array_equal(arr, ref_arr)
                # Yielded vectors must not share memory with cached reductions
                matches[0][0] *= 2

    def test_reduction_cache_is_local_to_generator(self):
        z = ZSLGenerator()
        other = ZSLGenerator()
        film_vectors = self.film.lattice.matrix[:2]
        substrate_vectors = self.substrate.lattice.matrix[:2]
        list(z(film_vectors, substrate_vectors))
        assert len(z._reduced_vector_sets) > 0
        assert len(other._reduced_vector_sets) == 0
        # the cache is not pickled, so generators still pickle (e.g. for multiprocessing)
        restored = pickle.loads(pickle.dumps(z))
        assert restored._reduced_vector_sets == {}
        assert restored.as_dict() == z.as_dict()

    def test_vector_shape(self):
        z = ZSLGenerator()
        vectors = self.film.lattice.matrix[:2]
        for film, substrate in ((vectors.reshape(3, 2), vectors), (vectors, vectors.reshape(3, 2))):
            with pytest.raises(ValueError, match=r"shape \(2, 3\)"):
                list(z.get_equiv_transformations([], film, substrate))
