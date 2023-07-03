from __future__ import annotations

from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.interfaces.substrate_analyzer import SubstrateAnalyzer
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.testing import test_structures


def test_substrate_analyzer_init():
    # Film VO2
    film = SpacegroupAnalyzer(test_structures.VO2, symprec=0.1).get_conventional_standard_structure()

    # Substrate TiO2
    substrate = SpacegroupAnalyzer(test_structures.TiO2, symprec=0.1).get_conventional_standard_structure()

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
