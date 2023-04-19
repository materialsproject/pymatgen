from __future__ import annotations

import unittest

import numpy as np
import pytest

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import (
    AllCoordinationGeometries,
    CoordinationGeometry,
    ExplicitPermutationsAlgorithm,
    SeparationPlane,
)
from pymatgen.util.testing import PymatgenTest

__author__ = "waroquiers"

allcg = AllCoordinationGeometries()


class FakeSite:
    def __init__(self, coords):
        self.coords = coords


class CoordinationGeometriesTest(PymatgenTest):
    def test_algorithms(self):
        expl_algo = ExplicitPermutationsAlgorithm(permutations=[[0, 1, 2], [1, 2, 3]])
        expl_algo2 = ExplicitPermutationsAlgorithm.from_dict(expl_algo.as_dict)
        assert expl_algo.permutations == expl_algo2.permutations

        sepplane_algos_oct = allcg["O:6"].algorithms
        assert len(sepplane_algos_oct[0].safe_separation_permutations()) == 24
        assert len(sepplane_algos_oct[1].safe_separation_permutations()) == 36

        sepplane_algos_oct_0 = SeparationPlane.from_dict(sepplane_algos_oct[0].as_dict)
        assert sepplane_algos_oct[0].plane_points == sepplane_algos_oct_0.plane_points
        assert sepplane_algos_oct[0].mirror_plane == sepplane_algos_oct_0.mirror_plane
        assert sepplane_algos_oct[0].ordered_plane == sepplane_algos_oct_0.ordered_plane
        assert sepplane_algos_oct[0].point_groups == sepplane_algos_oct_0.point_groups
        assert sepplane_algos_oct[0].ordered_point_groups == sepplane_algos_oct_0.ordered_point_groups
        assert all(
            np.array_equal(
                perm,
                sepplane_algos_oct_0.explicit_optimized_permutations[iperm],
            )
            for iperm, perm in enumerate(sepplane_algos_oct[0].explicit_optimized_permutations)
        )

        assert (
            str(sepplane_algos_oct[0])
            == "Separation plane algorithm with the following reference separation :\n[[4]] | [[0, 2, 1, 3]] | [[5]]"
        )

    def test_hints(self):
        hints = CoordinationGeometry.NeighborsSetsHints(hints_type="single_cap", options={"cap_index": 2, "csm_max": 8})
        myhints = hints.hints({"csm": 12.0})
        assert myhints == []

        hints2 = CoordinationGeometry.NeighborsSetsHints.from_dict(hints.as_dict())
        assert hints.hints_type == hints2.hints_type
        assert hints.options == hints2.options

    def test_coordination_geometry(self):
        cg_oct = allcg["O:6"]
        cg_oct2 = CoordinationGeometry.from_dict(cg_oct.as_dict())

        self.assertArrayAlmostEqual(cg_oct.central_site, cg_oct2.central_site)
        self.assertArrayAlmostEqual(cg_oct.points, cg_oct2.points)
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
        assert cg_oct.number_of_permutations == 720.0
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
        self.assertArrayAlmostEqual(cg_oct.faces(sites=sites, permutation=[0, 3, 2, 4, 5, 1]), faces)

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
        self.assertArrayAlmostEqual(cg_oct.faces(sites=sites), faces)

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
        self.assertArrayAlmostEqual(cg_oct.edges(sites=sites, permutation=[0, 3, 2, 4, 5, 1]), edges)

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
        self.assertArrayAlmostEqual(cg_oct.edges(sites=sites), edges)

        self.assertArrayAlmostEqual(
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
            "------------------------------------------------------------\n\n" in str(allcg)
        )

        assert (
            "Coordination geometry type : Trigonal plane (IUPAC: TP-3 || IUCr: [3l])\n\n"
            "  - coordination number : 3\n"
            "  - list of points :\n" in str(allcg)
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
            "FO:7",
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

        assert len(allcg.get_geometries()) == 68
        assert len(allcg.get_geometries(coordination=3)) == 3
        assert sorted(allcg.get_geometries(returned="mp_symbol")) == sorted(all_symbols)
        assert sorted(allcg.get_geometries(returned="mp_symbol", coordination=3)) == ["TL:3", "TS:3", "TY:3"]

        assert allcg.get_symbol_name_mapping(coordination=3) == {
            "TY:3": "Triangular non-coplanar",
            "TL:3": "Trigonal plane",
            "TS:3": "T-shaped",
        }
        assert allcg.get_symbol_cn_mapping(coordination=3) == {"TY:3": 3, "TL:3": 3, "TS:3": 3}
        assert sorted(allcg.get_implemented_geometries(coordination=4, returned="mp_symbol")) == [
            "S:4",
            "SS:4",
            "SY:4",
            "T:4",
        ]
        assert sorted(allcg.get_not_implemented_geometries(returned="mp_symbol")) == [
            "CO:11",
            "H:10",
            "S:10",
            "S:12",
            "UNCLEAR",
            "UNKNOWN",
        ]

        assert allcg.get_geometry_from_name("Octahedron").mp_symbol == cg_oct.mp_symbol
        with pytest.raises(LookupError) as exc_info:
            allcg.get_geometry_from_name("Octahedran")
        assert str(exc_info.value) == "No coordination geometry found with name 'Octahedran'"

        assert allcg.get_geometry_from_IUPAC_symbol("OC-6").mp_symbol == cg_oct.mp_symbol
        with pytest.raises(LookupError) as exc_info:
            allcg.get_geometry_from_IUPAC_symbol("OC-7")
        assert str(exc_info.value) == "No coordination geometry found with IUPAC symbol 'OC-7'"

        assert allcg.get_geometry_from_IUCr_symbol("[6o]").mp_symbol == cg_oct.mp_symbol
        with pytest.raises(LookupError) as exc_info:
            allcg.get_geometry_from_IUCr_symbol("[6oct]")
        assert str(exc_info.value) == "No coordination geometry found with IUCr symbol '[6oct]'"

        with pytest.raises(LookupError) as exc_info:
            allcg.get_geometry_from_mp_symbol("O:7")
        assert str(exc_info.value) == "No coordination geometry found with mp_symbol 'O:7'"

        assert (
            allcg.pretty_print(maxcn=4)
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
            allcg.pretty_print(maxcn=2, type="all_geometries_latex")
            == "\\subsection*{Coordination 1}\n\n\\begin{itemize}\n"
            "\\item S:1 $\\rightarrow$ Single neighbor (IUPAC : None - IUCr : $[$1l$]$)\n"
            "\\end{itemize}\n\n\\subsection*{Coordination 2}\n\n\\begin{itemize}\n"
            "\\item L:2 $\\rightarrow$ Linear (IUPAC : L-2 - IUCr : $[$2l$]$)\n"
            "\\item A:2 $\\rightarrow$ Angular (IUPAC : A-2 - IUCr : $[$2n$]$)\n"
            "\\end{itemize}\n\n"
        )
        assert (
            allcg.pretty_print(maxcn=2, type="all_geometries_latex_images")
            == "\\section*{Coordination 1}\n\n\\subsubsection*{S:1 : Single neighbor}\n\n"
            "IUPAC : None\n\nIUCr : [1l]\n\n\\begin{center}\n"
            "\\includegraphics[scale=0.15]{images/S_1.png}\n"
            "\\end{center}\n\n\\section*{Coordination 2}\n\n"
            "\\subsubsection*{L:2 : Linear}\n\nIUPAC : L-2\n\n"
            "IUCr : [2l]\n\n\\begin{center}\n\\includegraphics[scale=0.15]{images/L_2.png}\n"
            "\\end{center}\n\n\\subsubsection*{A:2 : Angular}\n\nIUPAC : A-2\n\nIUCr : [2n]\n\n"
            "\\begin{center}\n\\includegraphics[scale=0.15]{images/A_2.png}\n\\end{center}\n\n"
        )
        assert allcg.minpoints == {6: 2, 7: 2, 8: 2, 9: 2, 10: 2, 11: 2, 12: 2, 13: 3, 20: 2}
        assert allcg.maxpoints == {6: 5, 7: 5, 8: 6, 9: 7, 10: 6, 11: 5, 12: 8, 13: 6, 20: 10}
        assert allcg.maxpoints_inplane == {6: 5, 7: 5, 8: 6, 9: 7, 10: 6, 11: 5, 12: 8, 13: 6, 20: 10}
        assert allcg.separations_cg == {
            6: {
                (0, 3, 3): ["O:6", "T:6"],
                (1, 4, 1): ["O:6"],
                (0, 5, 1): ["PP:6"],
                (2, 2, 2): ["PP:6"],
                (0, 4, 2): ["T:6"],
            },
            7: {
                (1, 3, 3): ["ET:7", "FO:7"],
                (2, 3, 2): ["PB:7", "ST:7", "ET:7"],
                (1, 4, 2): ["ST:7", "FO:7"],
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


if __name__ == "__main__":
    unittest.main()
