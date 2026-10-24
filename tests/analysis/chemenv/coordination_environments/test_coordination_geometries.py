from __future__ import annotations

import re

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import (
    EXPLICIT_PERMUTATIONS,
    SEPARATION_PLANE,
    AbstractChemenvAlgorithm,
    AllCoordinationGeometries,
    CoordinationGeometry,
    ExplicitPermutationsAlgorithm,
    SeparationPlane,
)
from pymatgen.util.testing import MatSciTest

__author__ = "waroquiers"

all_cg = AllCoordinationGeometries()


class FakeSite:
    def __init__(self, coords):
        self.coords = coords


class TestCoordinationGeometries(MatSciTest):
    def test_algorithms(self):
        expl_algo = ExplicitPermutationsAlgorithm(permutations=[[0, 1, 2], [1, 2, 3]])
        expl_algo2 = ExplicitPermutationsAlgorithm.from_dict(expl_algo.as_dict())
        assert expl_algo.permutations == expl_algo2.permutations

        sep_plane_algos_oct = all_cg["O:6"].algorithms
        assert len(sep_plane_algos_oct[0].safe_separation_permutations()) == 24
        assert len(sep_plane_algos_oct[1].safe_separation_permutations()) == 36

        sep_plane_algos_oct_0 = SeparationPlane.from_dict(sep_plane_algos_oct[0].as_dict())
        assert sep_plane_algos_oct[0].plane_points == sep_plane_algos_oct_0.plane_points
        assert sep_plane_algos_oct[0].mirror_plane == sep_plane_algos_oct_0.mirror_plane
        assert sep_plane_algos_oct[0].ordered_plane == sep_plane_algos_oct_0.ordered_plane
        assert sep_plane_algos_oct[0].point_groups == sep_plane_algos_oct_0.point_groups
        assert sep_plane_algos_oct[0].ordered_point_groups == sep_plane_algos_oct_0.ordered_point_groups
        assert all(
            np.array_equal(perm, sep_plane_algos_oct_0.explicit_optimized_permutations[idx])
            for idx, perm in enumerate(sep_plane_algos_oct[0].explicit_optimized_permutations)
        )

        expected_str = (
            "Separation plane algorithm with the following reference separation :\n[[4]] | [[0, 2, 1, 3]] | [[5]]"
        )
        assert str(sep_plane_algos_oct[0]) == expected_str

    def test_hints(self):
        hints = CoordinationGeometry.NeighborsSetsHints(hints_type="single_cap", options={"cap_index": 2, "csm_max": 8})
        csm_hints = hints.hints({"csm": 12.0})
        assert csm_hints == []

        hints2 = CoordinationGeometry.NeighborsSetsHints.from_dict(hints.as_dict())
        assert hints.hints_type == hints2.hints_type
        assert hints.options == hints2.options

    def test_coordination_geometry(self):
        cg_oct = all_cg["O:6"]
        cg_oct2 = CoordinationGeometry.from_dict(cg_oct.as_dict())

        assert cg_oct.central_site == approx(cg_oct2.central_site)
        for p1, p2 in zip(cg_oct.points, cg_oct2.points, strict=True):
            assert p1 == approx(p2)
        assert (
            str(cg_oct) == "Coordination geometry type : Octahedron (IUPAC: OC-6 || IUCr: [6o])\n"
            "\n"
            "  - coordination number : 6\n"
            "  - list of points :\n"
            "    - [0.0, 0.0, 1.0]\n"
            "    - [0.0, 0.0, -1.0]\n"
            "    - [1.0, 0.0, 0.0]\n"
            "    - [-1.0, 0.0, 0.0]\n"
            "    - [0.0, 1.0, 0.0]\n"
            "    - [0.0, -1.0, 0.0]\n"
            "------------------------------------------------------------\n"
        )

        assert len(cg_oct) == 6
        assert cg_oct.ce_symbol == cg_oct.mp_symbol
        assert cg_oct.is_implemented()
        assert cg_oct.get_name() == "Octahedron"
        assert cg_oct.IUPAC_symbol == "OC-6"
        assert cg_oct.IUPAC_symbol_str == "OC-6"
        assert cg_oct.IUCr_symbol == "[6o]"
        assert cg_oct.IUCr_symbol_str == "[6o]"

        cg_oct.permutations_safe_override = True
        assert cg_oct.number_of_permutations == 720
        assert cg_oct.ref_permutation([0, 3, 2, 4, 5, 1]) == (0, 3, 1, 5, 2, 4)

        sites = [FakeSite(coords=pp) for pp in cg_oct.points]
        faces = [
            [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
            [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
            [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0]],
            [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
            [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
            [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
            [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0]],
            [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
        ]
        assert_allclose(cg_oct.faces(sites=sites, permutation=[0, 3, 2, 4, 5, 1]), faces)

        faces = [
            [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
            [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
            [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
            [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
        ]
        assert_allclose(cg_oct.faces(sites=sites), faces)

        edges = [
            [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
            [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
            [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0]],
            [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]],
            [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
            [[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
            [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
            [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
            [[0.0, 1.0, 0.0], [0.0, -1.0, 0.0]],
            [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
        ]
        assert_allclose(cg_oct.edges(sites=sites, permutation=[0, 3, 2, 4, 5, 1]), edges)

        edges = [
            [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
            [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]],
            [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
            [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0]],
            [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0]],
            [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]],
            [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0]],
            [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0]],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
            [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
        ]
        assert_allclose(cg_oct.edges(sites=sites), edges)

        assert_allclose(
            cg_oct.solid_angles(),
            [2.0943951, 2.0943951, 2.0943951, 2.0943951, 2.0943951, 2.0943951],
        )

        pmeshes = cg_oct.get_pmeshes(sites=sites)
        assert (
            pmeshes[0]["pmesh_string"] == "14\n     0.00000000      0.00000000      1.00000000\n"
            "     0.00000000      0.00000000     -1.00000000\n"
            "     1.00000000      0.00000000      0.00000000\n"
            "    -1.00000000      0.00000000      0.00000000\n"
            "     0.00000000      1.00000000      0.00000000\n"
            "     0.00000000     -1.00000000      0.00000000\n"
            "     0.33333333      0.33333333      0.33333333\n"
            "     0.33333333     -0.33333333      0.33333333\n"
            "    -0.33333333      0.33333333      0.33333333\n"
            "    -0.33333333     -0.33333333      0.33333333\n"
            "     0.33333333      0.33333333     -0.33333333\n"
            "     0.33333333     -0.33333333     -0.33333333\n"
            "    -0.33333333      0.33333333     -0.33333333\n"
            "    -0.33333333     -0.33333333     -0.33333333\n"
            "8\n4\n0\n2\n4\n0\n4\n0\n2\n5\n0\n4\n0\n3\n4\n0\n"
            "4\n0\n3\n5\n0\n4\n1\n2\n4\n1\n4\n1\n2\n5\n1\n4\n"
            "1\n3\n4\n1\n4\n1\n3\n5\n1\n"
        )

        assert (
            "\n#=======================================================#\n"
            "# List of coordination geometries currently implemented #\n"
            "#=======================================================#\n"
            "\nCoordination geometry type : Single neighbor (IUCr: [1l])\n\n"
            "  - coordination number : 1\n"
            "  - list of points :\n"
            "    - [0.0, 0.0, 1.0]\n"
            "------------------------------------------------------------\n\n" in str(all_cg)
        )

        assert (
            "Coordination geometry type : Trigonal plane (IUPAC: TP-3 || IUCr: [3l])\n\n"
            "  - coordination number : 3\n"
            "  - list of points :\n" in str(all_cg)
        )

        all_symbols = [
            "S:1",
            "L:2",
            "A:2",
            "TL:3",
            "TY:3",
            "TS:3",
            "T:4",
            "S:4",
            "SY:4",
            "SS:4",
            "PP:5",
            "S:5",
            "T:5",
            "O:6",
            "T:6",
            "PP:6",
            "PB:7",
            "ST:7",
            "ET:7",
            "FO:7",  # codespell:ignore fo
            "C:8",
            "SA:8",
            "SBT:8",
            "TBT:8",
            "DD:8",
            "DDPN:8",
            "HB:8",
            "BO_1:8",
            "BO_2:8",
            "BO_3:8",
            "TC:9",
            "TT_1:9",
            "TT_2:9",
            "TT_3:9",
            "HD:9",
            "TI:9",
            "SMA:9",
            "SS:9",
            "TO_1:9",
            "TO_2:9",
            "TO_3:9",
            "PP:10",
            "PA:10",
            "SBSA:10",
            "MI:10",
            "S:10",
            "H:10",
            "BS_1:10",
            "BS_2:10",
            "TBSA:10",
            "PCPA:11",
            "H:11",
            "SH:11",
            "CO:11",
            "DI:11",
            "I:12",
            "PBP:12",
            "TT:12",
            "C:12",
            "AC:12",
            "SC:12",
            "S:12",
            "HP:12",
            "HA:12",
            "SH:13",
            "DD:20",
            "UNKNOWN",
            "UNCLEAR",
        ]

        assert len(all_cg.get_geometries()) == 68
        assert len(all_cg.get_geometries(coordination=3)) == 3
        assert sorted(all_cg.get_geometries(returned="mp_symbol")) == sorted(all_symbols)
        assert sorted(all_cg.get_geometries(returned="mp_symbol", coordination=3)) == [
            "TL:3",
            "TS:3",
            "TY:3",
        ]

        assert all_cg.get_symbol_name_mapping(coordination=3) == {
            "TY:3": "Triangular non-coplanar",
            "TL:3": "Trigonal plane",
            "TS:3": "T-shaped",
        }
        assert all_cg.get_symbol_cn_mapping(coordination=3) == {
            "TY:3": 3,
            "TL:3": 3,
            "TS:3": 3,
        }
        assert sorted(all_cg.get_implemented_geometries(coordination=4, returned="mp_symbol")) == [
            "S:4",
            "SS:4",
            "SY:4",
            "T:4",
        ]
        assert sorted(all_cg.get_not_implemented_geometries(returned="mp_symbol")) == [
            "CO:11",
            "H:10",
            "S:10",
            "S:12",
            "UNCLEAR",
            "UNKNOWN",
        ]

        assert all_cg.get_geometry_from_name("Octahedron").mp_symbol == cg_oct.mp_symbol
        with pytest.raises(LookupError, match=re.escape("No coordination geometry found with name 'Octahedran'")):
            all_cg.get_geometry_from_name("Octahedran")

        assert all_cg.get_geometry_from_IUPAC_symbol("OC-6").mp_symbol == cg_oct.mp_symbol
        with pytest.raises(LookupError, match=re.escape("No coordination geometry found with IUPAC symbol 'OC-7'")):
            all_cg.get_geometry_from_IUPAC_symbol("OC-7")

        assert all_cg.get_geometry_from_IUCr_symbol("[6o]").mp_symbol == cg_oct.mp_symbol
        with pytest.raises(LookupError, match=re.escape("No coordination geometry found with IUCr symbol '[6oct]'")):
            all_cg.get_geometry_from_IUCr_symbol("[6oct]")

        with pytest.raises(LookupError, match=re.escape("No coordination geometry found with mp_symbol 'O:7'")):
            all_cg.get_geometry_from_mp_symbol("O:7")

        assert (
            all_cg.pretty_print(maxcn=4)
            == "+-------------------------+\n| Coordination geometries |\n+-------------------------+\n"
            "\n==>> CN = 1 <<==\n - S:1 : Single neighbor\n\n"
            "==>> CN = 2 <<==\n"
            " - L:2 : Linear\n - A:2 : Angular\n\n"
            "==>> CN = 3 <<==\n"
            " - TL:3 : Trigonal plane\n - TY:3 : Triangular non-coplanar\n - TS:3 : T-shaped\n\n"
            "==>> CN = 4 <<==\n - T:4 : Tetrahedron\n - S:4 : Square plane\n"
            " - SY:4 : Square non-coplanar\n - SS:4 : See-saw\n\n"
        )
        assert (
            all_cg.pretty_print(maxcn=2, type="all_geometries_latex")
            == "\\subsection*{Coordination 1}\n\n\\begin{itemize}\n"
            "\\item S:1 $\\rightarrow$ Single neighbor (IUPAC : None - IUCr : $[$1l$]$)\n"
            "\\end{itemize}\n\n\\subsection*{Coordination 2}\n\n\\begin{itemize}\n"
            "\\item L:2 $\\rightarrow$ Linear (IUPAC : L-2 - IUCr : $[$2l$]$)\n"
            "\\item A:2 $\\rightarrow$ Angular (IUPAC : A-2 - IUCr : $[$2n$]$)\n"
            "\\end{itemize}\n\n"
        )
        assert (
            all_cg.pretty_print(maxcn=2, type="all_geometries_latex_images")
            == "\\section*{Coordination 1}\n\n\\subsubsection*{S:1 : Single neighbor}\n\n"
            "IUPAC : None\n\nIUCr : [1l]\n\n\\begin{center}\n"
            "\\includegraphics[scale=0.15]{images/S_1.png}\n"
            "\\end{center}\n\n\\section*{Coordination 2}\n\n"
            "\\subsubsection*{L:2 : Linear}\n\nIUPAC : L-2\n\n"
            "IUCr : [2l]\n\n\\begin{center}\n\\includegraphics[scale=0.15]{images/L_2.png}\n"
            "\\end{center}\n\n\\subsubsection*{A:2 : Angular}\n\nIUPAC : A-2\n\nIUCr : [2n]\n\n"
            "\\begin{center}\n\\includegraphics[scale=0.15]{images/A_2.png}\n\\end{center}\n\n"
        )
        assert all_cg.minpoints == {
            6: 2,
            7: 2,
            8: 2,
            9: 2,
            10: 2,
            11: 2,
            12: 2,
            13: 3,
            20: 2,
        }
        assert all_cg.maxpoints == {
            6: 5,
            7: 5,
            8: 6,
            9: 7,
            10: 6,
            11: 5,
            12: 8,
            13: 6,
            20: 10,
        }
        assert all_cg.maxpoints_inplane == {
            6: 5,
            7: 5,
            8: 6,
            9: 7,
            10: 6,
            11: 5,
            12: 8,
            13: 6,
            20: 10,
        }
        assert all_cg.separations_cg == {
            6: {
                (0, 3, 3): ["O:6", "T:6"],
                (1, 4, 1): ["O:6"],
                (0, 5, 1): ["PP:6"],
                (2, 2, 2): ["PP:6"],
                (0, 4, 2): ["T:6"],
            },
            7: {
                (1, 3, 3): ["ET:7", "FO:7"],  # codespell:ignore fo
                (2, 3, 2): ["PB:7", "ST:7", "ET:7"],
                (1, 4, 2): ["ST:7", "FO:7"],  # codespell:ignore fo
                (1, 5, 1): ["PB:7"],
            },
            8: {
                (1, 6, 1): ["HB:8"],
                (0, 4, 4): ["C:8", "SA:8", "SBT:8"],
                (1, 4, 3): ["SA:8", "SBT:8", "BO_2:8", "BO_3:8"],
                (2, 4, 2): [
                    "C:8",
                    "TBT:8",
                    "DD:8",
                    "DDPN:8",
                    "HB:8",
                    "BO_1:8",
                    "BO_1:8",
                    "BO_2:8",
                    "BO_2:8",
                    "BO_3:8",
                    "BO_3:8",
                ],
            },
            9: {
                (3, 3, 3): [
                    "TT_1:9",
                    "TT_1:9",
                    "TT_2:9",
                    "SMA:9",
                    "SMA:9",
                    "TO_1:9",
                    "TO_3:9",
                ],
                (0, 6, 3): ["TC:9"],
                (2, 4, 3): [
                    "TC:9",
                    "TT_2:9",
                    "TT_3:9",
                    "TI:9",
                    "SS:9",
                    "TO_1:9",
                    "TO_1:9",
                    "TO_2:9",
                    "TO_3:9",
                ],
                (1, 3, 5): ["TI:9"],
                (1, 4, 4): ["TT_1:9", "SMA:9", "SS:9"],
                (2, 3, 4): ["TC:9"],
                (2, 5, 2): ["TT_3:9", "SS:9", "TO_2:9"],
                (1, 7, 1): ["HD:9"],
            },
            10: {
                (0, 5, 5): ["PP:10", "PA:10"],
                (3, 4, 3): ["PA:10", "SBSA:10", "MI:10", "BS_2:10", "TBSA:10"],
                (2, 6, 2): ["BS_1:10"],
                (2, 4, 4): ["PP:10", "MI:10", "BS_2:10"],
                (3, 3, 4): ["SBSA:10"],
                (1, 4, 5): ["BS_2:10"],
                (0, 4, 6): ["BS_1:10", "TBSA:10"],
            },
            11: {
                (4, 3, 4): ["PCPA:11"],
                (3, 4, 4): ["DI:11"],
                (1, 5, 5): ["PCPA:11", "DI:11"],
                (3, 5, 3): ["H:11"],
            },
            12: {
                (3, 3, 6): ["TT:12"],
                (2, 4, 6): ["TT:12"],
                (0, 6, 6): ["HP:12", "HA:12"],
                (3, 6, 3): ["C:12", "AC:12"],
                (4, 4, 4): ["I:12", "PBP:12", "C:12", "HP:12"],
                (0, 8, 4): ["SC:12"],
            },
            13: {(0, 6, 7): ["SH:13"]},
            20: {(5, 10, 5): ["DD:20"]},
        }


class FakeNeighborsSet:
    """Minimal stand-in for StructureEnvironments.NeighborsSet used by the hints."""

    def __init__(self, site_voronoi_indices):
        self.site_voronoi_indices = site_voronoi_indices

    def get_neighb_voronoi_indices(self, permutation):
        return [self.site_voronoi_indices[ii] for ii in permutation]


class TestAlgorithms:
    def test_abstract_str(self):
        algo = ExplicitPermutationsAlgorithm(permutations=[[0, 1]])
        assert str(algo) == EXPLICIT_PERMUTATIONS
        assert AbstractChemenvAlgorithm.__str__(algo) == ""

    def test_separation_plane_invalid_point_groups(self):
        with pytest.raises(RuntimeError, match="The number of points in the first group should be"):
            SeparationPlane(plane_points=[0, 1, 2], point_groups=[[3, 4], [5]], minimum_number_of_points=2)

    def test_separation_plane_properties(self):
        sep_plane = SeparationPlane(
            plane_points=[0, 1, 2, 3], point_groups=[[4], [5, 6, 7]], ordered_plane=True, minimum_number_of_points=2
        )
        assert sep_plane.algorithm_type == SEPARATION_PLANE
        assert sep_plane.permutations is None
        assert sep_plane.ordered_point_groups == [False, False]
        assert sep_plane.separation == (1, 4, 3)
        assert sep_plane.maximum_number_of_points == 4
        assert sep_plane.ref_separation_perm == [4, 0, 1, 2, 3, 5, 6, 7]
        assert sep_plane.argsorted_ref_separation_perm == [1, 2, 3, 4, 0, 5, 6, 7]
        round_trip = SeparationPlane.from_dict(sep_plane.as_dict())
        assert round_trip.permutations is None
        assert round_trip.separation == (1, 4, 3)

    @pytest.mark.parametrize(
        ("ordered_plane", "ordered_point_groups", "add_opposite", "expected"),
        [
            (False, None, False, 1 * 24 * 6),
            # cyclic orderings of the 4 plane points, both directions: 8 instead of 4! = 24
            (True, None, False, 1 * 8 * 6),
            (True, [False, True], False, 1 * 8 * 6),
            (True, [True, True], True, 2 * 1 * 8 * 6),
        ],
    )
    def test_safe_separation_permutations(self, ordered_plane, ordered_point_groups, add_opposite, expected):
        sep_plane = SeparationPlane(
            plane_points=[0, 1, 2, 3],
            point_groups=[[4], [5, 6, 7]],
            ordered_plane=True,
            ordered_point_groups=[True, True],
            minimum_number_of_points=2,
        )
        perms = sep_plane.safe_separation_permutations(
            ordered_plane=ordered_plane, ordered_point_groups=ordered_point_groups, add_opposite=add_opposite
        )
        assert len(perms) == expected
        assert all(sorted(perm) == list(range(8)) for perm in perms)
        # results are cached for identical arguments
        assert (
            sep_plane.safe_separation_permutations(
                ordered_plane=ordered_plane, ordered_point_groups=ordered_point_groups, add_opposite=add_opposite
            )
            is perms
        )
        if ordered_plane:
            assert [0, 4, 3, 2, 1, 5, 6, 7] in perms
            assert [0, 2, 1, 3, 4, 5, 6, 7] not in perms
        if add_opposite:
            assert [5, 6, 7, 1, 2, 3, 4, 0] in perms


class TestNeighborsSetsHints:
    def test_invalid_type(self):
        with pytest.raises(ValueError, match="Type 'quadruple_cap' for NeighborsSetsHints is not allowed"):
            CoordinationGeometry.NeighborsSetsHints(hints_type="quadruple_cap", options={})

    @pytest.mark.parametrize(
        ("hints_type", "options", "expected"),
        [
            ("single_cap", {"cap_index": 1}, [[10, 11, 12, 14]]),
            (
                "double_cap",
                {"first_cap_index": 0, "second_cap_index": 2},
                [[10, 11, 12, 13], [10, 11, 13, 14], [10, 11, 13]],
            ),
            (
                "triple_cap",
                {"first_cap_index": 0, "second_cap_index": 2, "third_cap_index": 4},
                [
                    [10, 11, 12, 13],
                    [10, 11, 13, 14],
                    [11, 12, 13, 14],
                    [11, 13, 14],
                    [11, 12, 13],
                    [10, 11, 13],
                    [11, 13],
                ],
            ),
        ],
    )
    def test_hints(self, hints_type, options, expected):
        hints = CoordinationGeometry.NeighborsSetsHints(hints_type=hints_type, options={"csm_max": 8.0, **options})
        hints_info = {"csm": 2.0, "nb_set": FakeNeighborsSet([10, 11, 12, 13, 14]), "permutation": [4, 3, 2, 1, 0]}
        assert [list(hint) for hint in hints.hints(hints_info)] == expected
        # no hint beyond the maximum csm
        assert hints.hints({**hints_info, "csm": 8.5}) == []


class TestCoordinationGeometry:
    @pytest.mark.parametrize(
        ("mp_symbol", "expected"),
        [
            ("S:1", 0.0),
            ("L:2", 0.0),
            ("T:4", np.sqrt(1.5) - 1),
            ("O:6", np.sqrt(2) - 1),
            ("C:8", np.sqrt(3) - 1),
            # cuboctahedron: anion-anion and cation-anion distances are equal (close packing)
            ("C:12", 1.0),
        ],
    )
    def test_pauling_stability_ratio(self, mp_symbol, expected):
        cg = AllCoordinationGeometries(only_symbols=[mp_symbol])[mp_symbol]
        assert cg.pauling_stability_ratio == approx(expected)
        # cached value
        assert cg.pauling_stability_ratio == approx(expected)

    @pytest.mark.parametrize(
        ("IUPAC_symbol", "IUCr_symbol", "expected_symbol"),
        [
            ("OC-6", "[6o]", " (IUPAC: OC-6 || IUCr: [6o])"),
            ("OC-6", None, " (IUPAC: OC-6)"),
            (None, "[6o]", " (IUCr: [6o])"),
            (None, None, ""),
        ],
    )
    def test_str_repr(self, IUPAC_symbol, IUCr_symbol, expected_symbol):
        cg = CoordinationGeometry(
            "X:6", "Some geometry", IUPAC_symbol=IUPAC_symbol, IUCr_symbol=IUCr_symbol, coordination=6
        )
        assert not cg.is_implemented()
        header = f"Coordination geometry type : Some geometry{expected_symbol}\n\n  - coordination number : 6\n"
        assert repr(cg) == f"{header}------------------------------------------------------------\n"
        assert (
            str(cg)
            == f"{header}... not yet implemented\n------------------------------------------------------------\n"
        )

    def test_number_of_permutations(self):
        assert CoordinationGeometry("X:5", "Some geometry", coordination=5).number_of_permutations == 120
        assert AllCoordinationGeometries(only_symbols=["O:6"])["O:6"].number_of_permutations == 720
        cg = CoordinationGeometry("X:3", "Some geometry", coordination=3)
        cg.permutations = [[0, 1, 2], [2, 1, 0]]
        assert cg.number_of_permutations == 2

    def test_edges_and_solid_angles(self):
        cg = all_cg["T:4"]
        coords = [np.array(pt) for pt in cg.points]
        perm = [3, 2, 1, 0]
        edges_from_coords = cg.edges(coords, permutation=perm, input="coords")
        edges_from_sites = cg.edges([FakeSite(coords=pt) for pt in coords], permutation=perm)
        assert len(edges_from_coords) == 6
        assert_allclose(edges_from_coords, edges_from_sites)
        with pytest.raises(RuntimeError, match="Invalid input for edges"):
            cg.edges(coords, input="points")

        cg = all_cg["SS:4"]
        solid_angles = cg.solid_angles()
        assert cg.solid_angles(permutation=[3, 2, 1, 0]) == solid_angles[::-1]

    def test_get_pmeshes_polygonal_faces(self):
        # The pentagonal prism has two pentagonal faces, which are split in 5 quadrangles each
        cg = all_cg["PP:10"]
        sites = [FakeSite(coords=np.array(pt)) for pt in cg.points]
        lines = cg.get_pmeshes(sites=sites)[0]["pmesh_string"].splitlines()
        n_faces = len(cg._faces)
        assert n_faces == 7
        assert int(lines[0]) == 10 + n_faces
        assert int(lines[1 + 10 + n_faces]) == 5 + 5 + 5
        # the face centres are the centroids of the faces
        face_centres = np.array([[float(xx) for xx in line.split()] for line in lines[11 : 11 + n_faces]])
        assert_allclose(
            face_centres, [np.mean([cg.points[ii] for ii in face], axis=0) for face in cg._faces], atol=1e-8
        )


class TestAllCoordinationGeometries:
    def test_contains_and_mappings(self):
        assert "O:6" in all_cg
        assert "O:7" not in all_cg
        mapping = all_cg.get_symbol_name_mapping()
        assert len(mapping) == len(all_cg.cg_list) == 68
        assert mapping["O:6"] == "Octahedron"

    def test_only_symbols(self):
        allcg = AllCoordinationGeometries(only_symbols=["O:6", "T:4"])
        assert sorted(allcg.get_implemented_geometries(returned="mp_symbol")) == ["O:6", "T:4"]
        assert allcg.separations_cg == {6: {(0, 3, 3): ["O:6"], (1, 4, 1): ["O:6"]}}
        assert allcg.minpoints == {6: 2}
        assert allcg.maxpoints == {6: 4}
        assert allcg.maxpoints_inplane == {6: 4}

    def test_repr(self):
        string = repr(all_cg)
        assert string.startswith("\n#=================================#\n# List of coordination geometries #\n")
        assert string.count("Coordination geometry type :") == 68
        assert "Coordination geometry type : Octahedron (IUPAC: OC-6 || IUCr: [6o])\n" in string

    @pytest.mark.parametrize(
        ("kwargs", "expected"),
        [
            ({"mp_symbol": "O:6"}, True),
            ({"mp_symbol": "O:7"}, False),
            ({"mp_symbol": "O:6", "IUPAC_symbol": "OC-6", "IUCr_symbol": "[6o]", "cn": "6"}, True),
            ({"mp_symbol": "O:6", "IUPAC_symbol": "TPR-6"}, False),
            ({"mp_symbol": "O:6", "IUCr_symbol": "[6p]"}, False),
            ({"mp_symbol": "O:6", "cn": 5}, False),
            ({"IUPAC_symbol": "OC-6"}, True),
            ({"IUPAC_symbol": "OC-7"}, False),
            ({"IUPAC_symbol": "OC-6", "IUCr_symbol": "[6p]"}, False),
            ({"IUPAC_symbol": "OC-6", "cn": 6}, True),
            ({"IUPAC_symbol": "OC-6", "cn": 5}, False),
            ({"IUCr_symbol": "[6o]"}, True),
            ({"IUCr_symbol": "[6o]", "cn": 6}, True),
            ({"IUCr_symbol": "[6o]", "cn": 5}, False),
            # legacy behavior: an unknown IUCr symbol is considered valid
            ({"IUCr_symbol": "[99x]"}, True),
        ],
    )
    def test_is_a_valid_coordination_geometry(self, kwargs, expected):
        assert all_cg.is_a_valid_coordination_geometry(**kwargs) is expected

    def test_is_a_valid_coordination_geometry_errors(self):
        with pytest.raises(NotImplementedError, match="not implemented for the name"):
            all_cg.is_a_valid_coordination_geometry(name="Octahedron")
        with pytest.raises(SyntaxError, match="at least one of mp_symbol, IUPAC_symbol and IUCr_symbol"):
            all_cg.is_a_valid_coordination_geometry(cn=6)

    def test_pretty_print(self):
        # geometries with neighbors sets hints are marked with a star
        out = all_cg.pretty_print(type="implemented_geometries", maxcn=13, additional_info={"nb_hints": True})
        assert " - SH:13 : Square-face capped hexagonal prism *\n" in out
        assert " - O:6 : Octahedron\n" in out
        assert out.count(" *\n") == sum(cg.neighbors_sets_hints is not None for cg in all_cg.cg_list)
        out = all_cg.pretty_print(type="all_geometries", maxcn=10)
        # geometries that are not implemented are also listed
        assert " - S:10 : " in out
        assert " - H:10 : " in out
        assert sum(line.startswith(" - ") for line in out.splitlines()) == sum(
            cg.coordination_number is not None and cg.coordination_number <= 10 for cg in all_cg.cg_list
        )
        out = all_cg.pretty_print(type="all_geometries_latex", maxcn=10)
        assert "\\item S:10 $\\rightarrow$" in out
        assert all_cg.pretty_print(type="all_geometries_latex_images", maxcn=0) == "\n"
