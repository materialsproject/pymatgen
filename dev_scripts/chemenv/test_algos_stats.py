#!/usr/bin/env python

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder, AbstractGeometry
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries, UNCLEAR_ENVIRONMENT_SYMBOL
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import StructureEnvironments
from pymatgen.analysis.chemenv.utils.scripts_utils import draw_cg
from pymatgen.analysis.chemenv.utils.scripts_utils import visualize
from pymatgen.analysis.chemenv.utils.chemenv_errors import NeighborsNotComputedChemenvError
from pymatgen.vis.structure_vtk import StructureVis
from pymatgen.io.cif import CifParser
from pymatgen.core.sites import PeriodicSite
from math import factorial

import numpy as np
import itertools
from random import shuffle


if __name__ == '__main__':

    cg_symbol = 'PP:6'

    allcg = AllCoordinationGeometries()
    cg = allcg[cg_symbol]

    for ialgo, algo in enumerate(cg.algorithms):
        algo._permutations = algo.explicit_permutations
        if algo.algorithm_type == 'EXPLICIT_PERMUTATIONS':
            raise ValueError('Do something for the explicit ones ... (these should anyway be by far ok!)')
        print('For ialgo {:d}, there are {:d} permutations'.format(ialgo, len(algo.permutations)))

    lgf = LocalGeometryFinder()
    lgf.setup_parameters(structure_refinement=lgf.STRUCTURE_REFINEMENT_NONE)


    myindices = range(cg.coordination_number)

    test = raw_input('Enter if you want to test all possible permutations ("all" or "a") or a given number of random permutations (i.e. "25")')

    if test == 'all' or test == 'a':
        perms_iterator = itertools.permutations(myindices)
        nperms = factorial(cg.coordination_number)
    else:
        try:
            nperms = int(test)
        except:
            raise ValueError('Could not turn {} into integer ...'.format(test))
        perms_iterator = []
        for ii in range(nperms):
            shuffle(myindices)
            perms_iterator.append(list(myindices))

    perms_used_algos = [dict() for algo in cg.algorithms]

    iperm = 1
    for indices_perm in perms_iterator:


        lgf.setup_test_perfect_environment(cg_symbol, indices=indices_perm, randomness=True, max_random_dist=0.1,
                                           random_translation=True, random_rotation=True, random_scale=True)

        lgf.perfect_geometry = AbstractGeometry.from_cg(cg=cg)

        print('Perm # {:d}/{:d} : '.format(iperm, nperms), indices_perm)

        algos_results = []
        for ialgo, algo in enumerate(cg.algorithms):
            if algo.algorithm_type == 'EXPLICIT_PERMUTATIONS':
                raise ValueError('Do something for the explicit ones ... (these should anyway be by far ok!)')

            results = lgf.coordination_geometry_symmetry_measures_separation_plane_newpmg(coordination_geometry=cg,
                                                                                          separation_plane_algo=algo,
                                                                                          testing=True)

            imin = np.argmin(results[0])
            tp = tuple(results[2][imin])
            if tp in perms_used_algos[ialgo]:
                perms_used_algos[ialgo][tp] += 1
            else:
                perms_used_algos[ialgo][tp] = 1

            print(results[0][imin])

            algos_results.append(min(results[0]))

            if not min(results[0]) < 1.0:
                print('Following is not close enough to 0.0 ...')
                raw_input(results)
        print('   => ', algos_results)
        iperm += 1
    algo_strings = []
    for ialgo, algo in enumerate(cg.algorithms):
        print('Permutations used for algorithm {:d} / {:d} ({:d} / {:d}) : '.format(ialgo, len(cg.algorithms),
                                                                                    len(perms_used_algos[ialgo]),
                                                                                    len(algo.permutations)))
        algo_string = '"explicit_optimized_permutations": ['
        pstrings = []
        for key, val in perms_used_algos[ialgo].items():
            print(key, val)
            pstring = '['
            pstring += ', '.join([str(ii) for ii in key])
            pstring += ']'
            pstrings.append(pstring)
        algo_string += ', '.join(pstrings)
        algo_string += ']'

        algo_strings.append(algo_string)

    print('Strings in algorithms for symbol {} :'.format(cg_symbol))
    for ialgo, algo in enumerate(cg.algorithms):
        print('Algorithm {:d} '.format(ialgo))
        print(algo_strings[ialgo])
