# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Development script to get the multiplicity of the separation facets for some model coordination environments
"""

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries

if __name__ == '__main__':

    allcg = AllCoordinationGeometries()

    cg_symbol = 'I:12'
    all_plane_points = []
    cg = allcg[cg_symbol]

    # I:12
    if cg_symbol == 'I:12':
        opposite_points = {0: 3,
                           1: 2,
                           2: 1,
                           3: 0,
                           4: 7,
                           5: 6,
                           6: 5,
                           7: 4,
                           8: 11,
                           9: 10,
                           10: 9,
                           11: 8
                           }
        edges = cg._edges
        for edge in edges:
            opposite_edge = [opposite_points[edge[0]], opposite_points[edge[1]]]
            equiv_plane = list(edge)
            equiv_plane.extend(opposite_edge)
            equiv_plane.sort()
            equiv_plane = tuple(equiv_plane)
            all_plane_points.append(equiv_plane)
        all_plane_points = list(set(all_plane_points))
        all_plane_points = [list(equiv_plane) for equiv_plane in all_plane_points]

    print('All plane points ({:d}) for {} : '.format(len(all_plane_points), cg_symbol))
    print(all_plane_points)
