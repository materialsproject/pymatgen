#!/usr/bin/env python


__author__ = 'waroquiers'

import unittest
import numpy as np
from pymatgen.util.testing import PymatgenTest

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import ExplicitPermutationsAlgorithm
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import SeparationPlane
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import CoordinationGeometry


allcg = AllCoordinationGeometries()


class FakeSite(object):

    def __init__(self, coords):
        self.coords = coords


class CoordinationGeometriesTest(PymatgenTest):

    def test_algorithms(self):
        expl_algo = ExplicitPermutationsAlgorithm(permutations=[[0, 1, 2], [1, 2, 3]])
        expl_algo2 = ExplicitPermutationsAlgorithm.from_dict(expl_algo.as_dict)
        self.assertEqual(expl_algo.permutations, expl_algo2.permutations)

        sepplane_algos_oct = allcg['O:6'].algorithms
        self.assertEqual(len(sepplane_algos_oct[0].safe_separation_permutations()), 24)
        self.assertEqual(len(sepplane_algos_oct[1].safe_separation_permutations()), 36)

        sepplane_algos_oct_0 = SeparationPlane.from_dict(sepplane_algos_oct[0].as_dict)
        self.assertEqual(sepplane_algos_oct[0].plane_points, sepplane_algos_oct_0.plane_points)
        self.assertEqual(sepplane_algos_oct[0].mirror_plane, sepplane_algos_oct_0.mirror_plane)
        self.assertEqual(sepplane_algos_oct[0].ordered_plane, sepplane_algos_oct_0.ordered_plane)
        self.assertEqual(sepplane_algos_oct[0].point_groups, sepplane_algos_oct_0.point_groups)
        self.assertEqual(sepplane_algos_oct[0].ordered_point_groups, sepplane_algos_oct_0.ordered_point_groups)
        self.assertTrue(all([np.array_equal(perm, sepplane_algos_oct_0.explicit_optimized_permutations[iperm])
                             for iperm, perm in enumerate(sepplane_algos_oct[0].explicit_optimized_permutations)]))

        self.assertEqual(sepplane_algos_oct[0].__str__(),
                         'Separation plane algorithm with the following reference separation :\n'
                         '[[4]] | [[0, 2, 1, 3]] | [[5]]')

    def test_hints(self):
        hints = CoordinationGeometry.NeighborsSetsHints(hints_type='single_cap',
                                                        options={'cap_index': 2, 'csm_max': 8})
        myhints = hints.hints({'csm': 12.0})
        self.assertEqual(myhints, [])

        hints2 = CoordinationGeometry.NeighborsSetsHints.from_dict(hints.as_dict())
        self.assertEqual(hints.hints_type, hints2.hints_type)
        self.assertEqual(hints.options, hints2.options)

    def test_coordination_geometry(self):
        cg_oct = allcg['O:6']
        cg_oct2 = CoordinationGeometry.from_dict(cg_oct.as_dict())

        self.assertArrayAlmostEqual(cg_oct.central_site, cg_oct2.central_site)
        self.assertArrayAlmostEqual(cg_oct.points, cg_oct2.points)
        self.assertEqual(cg_oct.__str__(), 'Coordination geometry type : Octahedron (IUPAC: OC-6 || IUCr: [6o])\n'
                                           '\n'
                                           '  - coordination number : 6\n'
                                           '  - list of points :\n'
                                           '    - [0.0, 0.0, 1.0]\n'
                                           '    - [0.0, 0.0, -1.0]\n'
                                           '    - [1.0, 0.0, 0.0]\n'
                                           '    - [-1.0, 0.0, 0.0]\n'
                                           '    - [0.0, 1.0, 0.0]\n'
                                           '    - [0.0, -1.0, 0.0]\n'
                                           '------------------------------------------------------------\n')

        self.assertEqual(cg_oct.__len__(), 6)
        self.assertEqual(cg_oct.ce_symbol, cg_oct.mp_symbol)
        self.assertTrue(cg_oct.is_implemented())
        self.assertEqual(cg_oct.get_name(), 'Octahedron')
        self.assertEqual(cg_oct.IUPAC_symbol, 'OC-6')
        self.assertEqual(cg_oct.IUPAC_symbol_str, 'OC-6')
        self.assertEqual(cg_oct.IUCr_symbol, '[6o]')
        self.assertEqual(cg_oct.IUCr_symbol_str, '[6o]')

        cg_oct.permutations_safe_override = True
        self.assertEqual(cg_oct.number_of_permutations, 720.0)
        self.assertEqual(cg_oct.ref_permutation([0, 3, 2, 4, 5, 1]), (0, 3, 1, 5, 2, 4))

        sites = [FakeSite(coords=pp) for pp in cg_oct.points]
        faces = [[[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
                 [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                 [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0]],
                 [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                 [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
                 [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]],
                 [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0]],
                 [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]]
        self.assertArrayAlmostEqual(cg_oct.faces(sites=sites, permutation=[0, 3, 2, 4, 5, 1]), faces)

        faces = [[[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                 [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
                 [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                 [[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
                 [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                 [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
                 [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                 [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]]
        self.assertArrayAlmostEqual(cg_oct.faces(sites=sites), faces)

        edges = [[[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
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
                 [[0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]]
        self.assertArrayAlmostEqual(cg_oct.edges(sites=sites, permutation=[0, 3, 2, 4, 5, 1]), edges)

        edges = [[[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
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
                 [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]]
        self.assertArrayAlmostEqual(cg_oct.edges(sites=sites), edges)

        self.assertArrayAlmostEqual(cg_oct.solid_angles(),
                                    [2.0943951, 2.0943951, 2.0943951, 2.0943951, 2.0943951, 2.0943951])

        pmeshes = cg_oct.get_pmeshes(sites=sites)
        self.assertEqual(pmeshes[0]['pmesh_string'],
                         '14\n     0.00000000      0.00000000      1.00000000\n'
                         '     0.00000000      0.00000000     -1.00000000\n'
                         '     1.00000000      0.00000000      0.00000000\n'
                         '    -1.00000000      0.00000000      0.00000000\n'
                         '     0.00000000      1.00000000      0.00000000\n'
                         '     0.00000000     -1.00000000      0.00000000\n'
                         '     0.33333333      0.33333333      0.33333333\n'
                         '     0.33333333     -0.33333333      0.33333333\n'
                         '    -0.33333333      0.33333333      0.33333333\n'
                         '    -0.33333333     -0.33333333      0.33333333\n'
                         '     0.33333333      0.33333333     -0.33333333\n'
                         '     0.33333333     -0.33333333     -0.33333333\n'
                         '    -0.33333333      0.33333333     -0.33333333\n'
                         '    -0.33333333     -0.33333333     -0.33333333\n'
                         '8\n4\n0\n2\n4\n0\n4\n0\n2\n5\n0\n4\n0\n3\n4\n0\n'
                         '4\n0\n3\n5\n0\n4\n1\n2\n4\n1\n4\n1\n2\n5\n1\n4\n'
                         '1\n3\n4\n1\n4\n1\n3\n5\n1\n')

        allcg_str = allcg.__str__()
        self.assertTrue('\n#=======================================================#\n'
                        '# List of coordination geometries currently implemented #\n'
                        '#=======================================================#\n'
                        '\nCoordination geometry type : Single neighbor (IUCr: [1l])\n\n'
                        '  - coordination number : 1\n'
                        '  - list of points :\n'
                        '    - [0.0, 0.0, 1.0]\n'
                        '------------------------------------------------------------\n\n' in allcg_str)

        self.assertTrue('Coordination geometry type : Trigonal plane (IUPAC: TP-3 || IUCr: [3l])\n\n'
                        '  - coordination number : 3\n'
                        '  - list of points :\n' in allcg_str)

        all_symbols = [u'S:1', u'L:2', u'A:2', u'TL:3', u'TY:3', u'TS:3', u'T:4', u'S:4', u'SY:4', u'SS:4',
                       u'PP:5', u'S:5', u'T:5', u'O:6', u'T:6', u'PP:6', u'PB:7', u'ST:7', u'ET:7', u'FO:7',
                       u'C:8', u'SA:8', u'SBT:8', u'TBT:8', u'DD:8', u'DDPN:8', u'HB:8', u'BO_1:8', u'BO_2:8',
                       u'BO_3:8', u'TC:9', u'TT_1:9', u'TT_2:9', u'TT_3:9', u'HD:9', u'TI:9', u'SMA:9', u'SS:9',
                       u'TO_1:9', u'TO_2:9', u'TO_3:9', u'PP:10', u'PA:10', u'SBSA:10', u'MI:10', u'S:10',
                       u'H:10', u'BS_1:10', u'BS_2:10', u'TBSA:10', u'PCPA:11', u'H:11', u'SH:11', u'CO:11',
                       u'DI:11', u'I:12', u'PBP:12', u'TT:12', u'C:12', u'AC:12', u'SC:12', u'S:12', u'HP:12',
                       u'HA:12', u'SH:13', u'DD:20', u'UNKNOWN', u'UNCLEAR']

        self.assertEqual(len(allcg.get_geometries()), 68)
        self.assertEqual(len(allcg.get_geometries(coordination=3)), 3)
        self.assertEqual(sorted(allcg.get_geometries(returned='mp_symbol')), sorted(all_symbols))
        self.assertEqual(sorted(allcg.get_geometries(returned='mp_symbol', coordination=3)),
                         ['TL:3', 'TS:3', 'TY:3'])

        self.assertEqual(allcg.get_symbol_name_mapping(coordination=3),
                         {u'TY:3': u'Triangular non-coplanar', u'TL:3': u'Trigonal plane', u'TS:3': u'T-shaped'})
        self.assertEqual(allcg.get_symbol_cn_mapping(coordination=3),
                         {u'TY:3': 3, u'TL:3': 3, u'TS:3': 3})
        self.assertEqual(sorted(allcg.get_implemented_geometries(coordination=4, returned='mp_symbol')),
                         [u'S:4', u'SS:4', u'SY:4', u'T:4'])
        self.assertEqual(sorted(allcg.get_not_implemented_geometries(returned='mp_symbol')),
                         [u'CO:11', u'DD:20', u'H:10', u'S:10', u'S:12', u'UNCLEAR', u'UNKNOWN'])

        self.assertEqual(allcg.get_geometry_from_name('Octahedron').mp_symbol, cg_oct.mp_symbol)
        with self.assertRaises(LookupError) as cm:
            allcg.get_geometry_from_name('Octahedran')
        self.assertEqual(str(cm.exception), 'No coordination geometry found with name "Octahedran"')

        self.assertEqual(allcg.get_geometry_from_IUPAC_symbol('OC-6').mp_symbol, cg_oct.mp_symbol)
        with self.assertRaises(LookupError) as cm:
            allcg.get_geometry_from_IUPAC_symbol('OC-7')
        self.assertEqual(str(cm.exception), 'No coordination geometry found with IUPAC symbol "OC-7"')

        self.assertEqual(allcg.get_geometry_from_IUCr_symbol('[6o]').mp_symbol, cg_oct.mp_symbol)
        with self.assertRaises(LookupError) as cm:
            allcg.get_geometry_from_IUCr_symbol('[6oct]')
        self.assertEqual(str(cm.exception), 'No coordination geometry found with IUCr symbol "[6oct]"')

        with self.assertRaises(LookupError) as cm:
            allcg.get_geometry_from_mp_symbol('O:7')
        self.assertEqual(str(cm.exception), 'No coordination geometry found with mp_symbol "O:7"')

        self.assertEqual(allcg.pretty_print(maxcn=4),
                         '+-------------------------+\n| Coordination geometries |\n+-------------------------+\n'
                         '\n==>> CN = 1 <<==\n - S:1 : Single neighbor\n\n'
                         '==>> CN = 2 <<==\n'
                         ' - L:2 : Linear\n - A:2 : Angular\n\n'
                         '==>> CN = 3 <<==\n'
                         ' - TL:3 : Trigonal plane\n - TY:3 : Triangular non-coplanar\n - TS:3 : T-shaped\n\n'
                         '==>> CN = 4 <<==\n - T:4 : Tetrahedron\n - S:4 : Square plane\n'
                         ' - SY:4 : Square non-coplanar\n - SS:4 : See-saw\n\n'                        )
        self.assertEqual(allcg.pretty_print(maxcn=2, type='all_geometries_latex'),
                         '\\subsection*{Coordination 1}\n\n\\begin{itemize}\n'
                         '\\item S:1 $\\rightarrow$ Single neighbor (IUPAC : None - IUCr : $[$1l$]$)\n'
                         '\\end{itemize}\n\n\\subsection*{Coordination 2}\n\n\\begin{itemize}\n'
                         '\\item L:2 $\\rightarrow$ Linear (IUPAC : L-2 - IUCr : $[$2l$]$)\n'
                         '\\item A:2 $\\rightarrow$ Angular (IUPAC : A-2 - IUCr : $[$2n$]$)\n'
                         '\\end{itemize}\n\n')
        self.assertEqual(allcg.pretty_print(maxcn=2, type='all_geometries_latex_images'),
                         '\\section*{Coordination 1}\n\n\\subsubsection*{S:1 : Single neighbor}\n\n'
                         'IUPAC : None\n\nIUCr : [1l]\n\n\\begin{center}\n'
                         '\\includegraphics[scale=0.15]{images/S_1.png}\n'
                         '\\end{center}\n\n\\section*{Coordination 2}\n\n'
                         '\\subsubsection*{L:2 : Linear}\n\nIUPAC : L-2\n\n'
                         'IUCr : [2l]\n\n\\begin{center}\n\\includegraphics[scale=0.15]{images/L_2.png}\n'
                         '\\end{center}\n\n\\subsubsection*{A:2 : Angular}\n\nIUPAC : A-2\n\nIUCr : [2n]\n\n'
                         '\\begin{center}\n\\includegraphics[scale=0.15]{images/A_2.png}\n\\end{center}\n\n')
        self.assertDictEqual(allcg.minpoints, {6: 2, 7: 2, 8: 2, 9: 2, 10: 2, 11: 2, 12: 2, 13: 3})
        self.assertDictEqual(allcg.maxpoints, {6: 5, 7: 5, 8: 6, 9: 7, 10: 6, 11: 5, 12: 8, 13: 6})
        self.assertDictEqual(allcg.maxpoints_inplane, {6: 5, 7: 5, 8: 6, 9: 7, 10: 6, 11: 5, 12: 8, 13: 6})
        self.assertDictEqual(allcg.separations_cg, {6: {(0, 3, 3): [u'O:6', u'T:6'],
                                                        (1, 4, 1): [u'O:6'],
                                                        (0, 5, 1): [u'PP:6'],
                                                        (2, 2, 2): [u'PP:6'],
                                                        (0, 4, 2): [u'T:6']},
                                                    7: {(1, 3, 3): [u'ET:7', u'FO:7'],
                                                        (2, 3, 2): [u'PB:7', u'ST:7', u'ET:7'],
                                                        (1, 4, 2): [u'ST:7', u'FO:7'],
                                                        (1, 5, 1): [u'PB:7']},
                                                    8: {(1, 6, 1): [u'HB:8'],
                                                        (0, 4, 4):
                                                            [u'C:8', u'SA:8', u'SBT:8'],
                                                        (1, 4, 3): [u'SA:8', u'SBT:8', u'BO_2:8', u'BO_3:8'],
                                                        (2, 4, 2): [u'C:8', u'TBT:8', u'DD:8', u'DDPN:8', u'HB:8',
                                                                    u'BO_1:8', u'BO_1:8', u'BO_2:8', u'BO_2:8',
                                                                    u'BO_3:8', u'BO_3:8']},
                                                    9: {(3, 3, 3): [u'TT_1:9', u'TT_1:9', u'TT_2:9', u'SMA:9',
                                                                    u'SMA:9', u'TO_1:9', u'TO_3:9'],
                                                        (0, 6, 3): [u'TC:9'],
                                                        (2, 4, 3): [u'TC:9', u'TT_2:9', u'TT_3:9', u'TI:9',
                                                                    u'SS:9', u'TO_1:9', u'TO_1:9', u'TO_2:9',
                                                                    u'TO_3:9'],
                                                        (1, 3, 5): [u'TI:9'],
                                                        (1, 4, 4): [u'TT_1:9', u'SMA:9', u'SS:9'],
                                                        (2, 3, 4): [u'TC:9'],
                                                        (2, 5, 2): [u'TT_3:9', u'SS:9', u'TO_2:9'],
                                                        (1, 7, 1): [u'HD:9']},
                                                    10: {(0, 5, 5): [u'PP:10', u'PA:10'],
                                                         (3, 4, 3): [u'PA:10', u'SBSA:10', u'MI:10',
                                                                     u'BS_2:10', u'TBSA:10'],
                                                         (2, 6, 2): [u'BS_1:10'],
                                                         (2, 4, 4): [u'PP:10', u'MI:10', u'BS_2:10'],
                                                         (3, 3, 4): [u'SBSA:10'],
                                                         (1, 4, 5): [u'BS_2:10'],
                                                         (0, 4, 6): [u'BS_1:10', u'TBSA:10']},
                                                    11: {(4, 3, 4): [u'PCPA:11'],
                                                         (3, 4, 4): [u'DI:11'],
                                                         (1, 5, 5): [u'PCPA:11', u'DI:11'],
                                                         (3, 5, 3): [u'H:11']},
                                                    12: {(3, 3, 6): [u'TT:12'],
                                                         (2, 4, 6): [u'TT:12'],
                                                         (0, 6, 6): [u'HP:12', u'HA:12'],
                                                         (3, 6, 3): [u'C:12', u'AC:12'],
                                                         (4, 4, 4): [u'I:12', u'PBP:12', u'C:12', u'HP:12'],
                                                         (0, 8, 4): [u'SC:12']},
                                                    13: {(0, 6, 7): [u'SH:13']}})

if __name__ == "__main__":
    unittest.main()