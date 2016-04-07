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


    allcg = AllCoordinationGeometries()

    test = raw_input('Standard ("s", all permutations for cn <= 6, 500 random permutations for cn > 6) or on demand')
    if test == 's':
        standard = True
    else:
        standard = False

    for coordination in range(1, 9):
        print('IN COORDINATION {:d}'.format(coordination))
        symbol_name_mapping = allcg.get_symbol_name_mapping(coordination=coordination)

        if standard:
            if coordination > 6:
                test = '500'
            else:
                test = 'all'
        else:
            test = raw_input('Enter if you want to test all possible permutations ("all" or "a") or a given number of random permutations (i.e. "25")')
        myindices = range(coordination)

        if test == 'all' or test == 'a':
            perms_type = 'all'
            perms_iterator = itertools.permutations(myindices)
            nperms = factorial(coordination)
        else:
            perms_type = 'explicit'
            try:
                nperms = int(test)
            except:
                raise ValueError('Could not turn {} into integer ...'.format(test))
            perms_iterator = []
            for ii in range(nperms):
                shuffle(myindices)
                perms_iterator.append(list(myindices))

        for cg_symbol, cg_name in symbol_name_mapping.items():

            print('Testing {} ({})'.format(cg_symbol, cg_name))

            cg = allcg[cg_symbol]
            if cg.points is None:
                continue

            lgf = LocalGeometryFinder()
            lgf.setup_parameters(structure_refinement=lgf.STRUCTURE_REFINEMENT_NONE)



            # Reinitialize the itertools permutations
            if perms_type == 'all':
                perms_iterator = itertools.permutations(myindices)

            #Loop on the permutations
            iperm = 1
            for indices_perm in perms_iterator:


                lgf.setup_test_perfect_environment(cg_symbol, indices=indices_perm,
                                                   randomness=True, max_random_dist=0.1,
                                                   random_translation=True, random_rotation=True, random_scale=True)

                lgf.perfect_geometry = AbstractGeometry.from_cg(cg=cg)

                print('Perm # {:d}/{:d} : '.format(iperm, nperms), indices_perm)

                algos_results = []
                for algo in cg.algorithms:
                    if algo.algorithm_type == 'EXPLICIT_PERMUTATIONS':
                        results = lgf.coordination_geometry_symmetry_measures_newpmg(coordination_geometry=cg)
                        # raise ValueError('Do something for the explicit ones ... (these should anyway be by far ok!)')
                    else:
                        results = lgf.coordination_geometry_symmetry_measures_separation_plane_newpmg(coordination_geometry=cg,
                                                                                                      separation_plane_algo=algo)
                    algos_results.append(min(results[0]))

                    if not min(results[0]) < 1.5:
                        print('Following is not close to 0.0 ...')
                        raw_input(results)
                print('   => ', algos_results)
                iperm += 1