from __future__ import annotations

import numpy as np
from numpy.testing import assert_allclose

from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.interfaces.substrate_analyzer import SubstrateAnalyzer
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.testing import MatSciTest

VO2 = MatSciTest.get_structure("VO2")
TiO2 = MatSciTest.get_structure("TiO2")

# Film VO2
film = SpacegroupAnalyzer(VO2, symprec=0.1).get_conventional_standard_structure()
# Substrate TiO2
substrate = SpacegroupAnalyzer(TiO2, symprec=0.1).get_conventional_standard_structure()


def test_substrate_analyzer_init():
    film_elastic_tensor = ElasticTensor.from_voigt(
        [
            [324.32, 187.3, 170.92, 0, 0, 0],
            [187.3, 324.32, 170.92, 0, 0, 0],
            [170.92, 170.92, 408.41, 0, 0, 0],
            [0, 0, 0, 150.73, 0, 0],
            [0, 0, 0, 0, 150.73, 0],
            [0, 0, 0, 0, 0, 238.74],
        ]
    )

    analyzer = SubstrateAnalyzer()

    matches = list(analyzer.calculate(film, substrate, film_elastic_tensor))
    assert len(matches) == 296
    for match in matches:
        assert isinstance(match.match_area, float)


def test_generate_surface_vectors():
    film_miller_indices = [(1, 0, 0)]
    substrate_miller_indices = [(1, 1, 1)]

    vector_sets = SubstrateAnalyzer().generate_surface_vectors(
        film, substrate, film_miller_indices, substrate_miller_indices
    )
    assert len(vector_sets) == 1
    film_vectors, substrate_vectors, film_millers, substrate_millers = vector_sets[0]

    assert [film_millers] == film_miller_indices
    assert [substrate_millers] == substrate_miller_indices
    assert_allclose(
        [np.linalg.norm(vector) for vector in film_vectors],
        [3.03542922, 4.51502263],
        atol=1e-6,
    )
    assert_allclose(
        [np.linalg.norm(vector) for vector in substrate_vectors],
        [7.613414525142487, 12.870736246973966],
        atol=1e-6,
    )
    assert_allclose(np.linalg.norm(np.cross(*film_vectors)), 13.705031620063249, atol=1e-6)
    assert_allclose(np.linalg.norm(np.cross(*substrate_vectors)), 97.52451926178126, atol=1e-6)
    assert_allclose(np.dot(*film_vectors), 0.0, atol=1e-6)
    assert_allclose(np.dot(*substrate_vectors), 9.542394617977376, atol=1e-6)
