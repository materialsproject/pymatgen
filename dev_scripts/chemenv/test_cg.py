#!/usr/bin/env python

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import AbstractGeometry
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from pymatgen.analysis.chemenv.utils.scripts_utils import draw_cg
from pymatgen.vis.structure_vtk import StructureVis

import numpy as np
from random import shuffle


if __name__ == '__main__':

    # cg_symbol = raw_input('Enter symbol of the geometry you want to test : ')
    cg_symbol = 'PP:6'

    allcg = AllCoordinationGeometries()
    cg = allcg[cg_symbol]

    lgf = LocalGeometryFinder()
    lgf.setup_parameters(structure_refinement=lgf.STRUCTURE_REFINEMENT_NONE)


    myindices = range(cg.coordination_number)

    test_indices = range(cg.coordination_number)
    shuffle(test_indices)
    # test = raw_input('Enter "r" if you want to test one random permutation, "s" for the standard permutation '
    #                  '(i.e. "{}")\nor enter the permutation by hand '
    #                  '(e.g. "{}" or "{}") : '.format('-'.join([str(ii) for ii in myindices]),
    #                                                 ' '.join([str(ii) for ii in test_indices]),
    #                                                 '-'.join([str(ii) for ii in test_indices])))
    test = 's'

    if test == 'r':
        shuffle(myindices)
        permutation = myindices
    elif test == 's':
        permutation = myindices
    else:
        try:
            permutation = [int(part) for part in test.split(' ')]
        except:
            try:
                permutation = [int(part) for part in test.split('-')]
            except:
                raise ValueError('Could not turn {} into permutation ...'.format(test))

    print('Permutation : {}'.format('-'.join([str(ii) for ii in permutation])))

    lgf.setup_test_perfect_environment(cg_symbol, indices=permutation, randomness=True, max_random_dist=0.01)

    lgf.perfect_geometry = AbstractGeometry.from_cg(cg=cg)

    results = lgf.coordination_geometry_symmetry_measures_newpmg(coordination_geometry=cg,
                                                                 tested_permutations=False)
    csms, permutations, algos, local2perfect_maps, perfect2local_maps = results
    imin = np.argmin(csms)

    if not csms[imin] < 0.5:
        print('Following is not close 0.0 ...')
        raw_input(results)

    vis = StructureVis(show_polyhedron=False, show_unit_cell=True)
    vis.show_help = False
    vis.set_structure(lgf.structure)
    draw_cg(vis, lgf.structure[0], lgf.structure[1:], cg=cg,
            # perfect2local_map=ce[1]['perfect2local_map'])
            perm=permutations[imin])
    for ineighbor, neighbor in enumerate(lgf.structure[1:]):
        vis.add_text(neighbor.coords, '{}'.format(ineighbor), color=(0, 0, 0))
    print('Number of permutations tested : ', len(csms))
    vis.show()

