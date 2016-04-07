#!/usr/bin/env python

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import AbstractGeometry
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from pymatgen.analysis.chemenv.utils.scripts_utils import draw_cg
from math import factorial

import numpy as np
import itertools
from random import shuffle
import time


if __name__ == '__main__':

    cg_symbol = raw_input('Enter symbol of the geometry you want to test : ')

    allcg = AllCoordinationGeometries()
    cg = allcg[cg_symbol]

    lgf = LocalGeometryFinder()
    lgf.setup_parameters(structure_refinement=lgf.STRUCTURE_REFINEMENT_NONE)


    myindices = range(cg.coordination_number)

    test = raw_input('Enter if you want to test all possible permutations ("all" or "a") or a given number of random permutations (i.e. "25")')

    if test == 'all' or test == 'a':
        iterator_type = 'all'
        nperms = factorial(cg.coordination_number)
    else:
        iterator_type = 'random'
        try:
            nperms = int(test)
        except:
            raise ValueError('Could not turn {} into integer ...'.format(test))
        perms_iterator = []
        for ii in range(nperms):
            shuffle(myindices)
            perms_iterator.append(list(myindices))


    for tested_permutations in [True, False]:
        iperm = 1
        t1 = time.clock()
        if iterator_type == 'all':
            perms_iterator = itertools.permutations(myindices)
        for indices_perm in perms_iterator:


            lgf.setup_test_perfect_environment(cg_symbol, indices=indices_perm)

            lgf.perfect_geometry = AbstractGeometry.from_cg(cg=cg)

            print('Perm # {:d}/{:d} : '.format(iperm, nperms), indices_perm)

            results = lgf.coordination_geometry_symmetry_measures_newpmg(coordination_geometry=cg,
                                                                         tested_permutations=tested_permutations)
            if not np.isclose(min(results[0]), 0.0):
                print('Following is not 0.0 ...')
                raw_input(results)

            print('Number of permutations tested : ', len(results[0]))

            iperm += 1
        t2 = time.clock()
        print('tested_permutaitons = ', tested_permutations)
        print('Time to test {:d} permutations for geometry "{}" (symbol "{}") : {:.2f} seconds'.format(nperms,
                                                                                                       cg.name,
                                                                                                       cg_symbol,
                                                                                                       t2-t1))
        raw_input('WITH two algos')