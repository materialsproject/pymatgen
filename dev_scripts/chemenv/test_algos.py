#!/usr/bin/env python

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder, AbstractGeometry
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries, UNCLEAR_ENVIRONMENT_SYMBOL
from math import factorial

import numpy as np
import itertools
from random import shuffle
import time


if __name__ == '__main__':

    cg_symbol = 'T:6'

    allcg = AllCoordinationGeometries()
    cg = allcg[cg_symbol]

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

    iperm = 1
    t1 = time.clock()
    for indices_perm in perms_iterator:


        lgf.setup_test_perfect_environment(cg_symbol, indices=indices_perm)

        lgf.perfect_geometry = AbstractGeometry.from_cg(cg=cg)

        print('Perm # {:d}/{:d} : '.format(iperm, nperms), indices_perm)

        algos_results = []
        for algo in cg.algorithms:
            print(algo)
            if algo.algorithm_type == 'EXPLICIT_PERMUTATIONS':
                raise ValueError('Do something for the explicit ones ... (these should anyway be by far ok!)')

            results = lgf.coordination_geometry_symmetry_measures_separation_plane_newpmg(coordination_geometry=cg,
                                                                                          separation_plane_algo=algo,
                                                                                          tested_permutations=False)
            print('Number of permutations tested : ', len(results[0]))
            algos_results.append(min(results[0]))

            if not np.isclose(min(results[0]), 0.0):
                print('Following is not 0.0 ...')
                raw_input(results)
        print('   => ', algos_results)
        iperm += 1
    t2 = time.clock()
    print('Time to test {:d} permutations for geometry "{}" (symbol "{}") : {:.2f} seconds'.format(nperms,
                                                                                                   cg.name,
                                                                                                   cg_symbol,
                                                                                                   t2-t1))