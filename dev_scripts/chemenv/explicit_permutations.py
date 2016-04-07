#!/usr/bin/env python

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder, AbstractGeometry
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries, CoordinationGeometry
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import ExplicitPermutationsAlgorithm


import numpy as np
import itertools
import json
import os


class Algo(object):
    pass


if __name__ == '__main__':

    cg_symbol = 'TY:3'

    allcg = AllCoordinationGeometries()
    cg = allcg[cg_symbol]

    algo = Algo()
    algo.permutations = []
    for perm in itertools.permutations(range(cg.coordination)):
        algo.permutations.append(perm)

    lgf = LocalGeometryFinder()
    lgf.setup_parameters(structure_refinement=lgf.STRUCTURE_REFINEMENT_NONE)
    lgf.setup_test_perfect_environment(cg_symbol, randomness=True, random_indices=False)

    lgf.perfect_geometry = AbstractGeometry.from_cg(cg=cg)

    (csms, perms, algos, local2perfect_maps, perfect2local_maps) = lgf.coordination_geometry_symmetry_measures_standard_newpmg(coordination_geometry=cg, algo=algo)


    csms_with_recorded_permutation = []
    explicit_permutations = []
    for icsm, csm in enumerate(csms):
        found = False
        for csm2 in csms_with_recorded_permutation:
            if np.isclose(csm, csm2):
                found = True
                break
        if not found:
            csms_with_recorded_permutation.append(csm)
            explicit_permutations.append(perms[icsm])

    print('Permutations found : ')
    print(explicit_permutations)

    print('Current algorithm(s) :')
    for algo in cg.algorithms:
        print(algo)
        if algo.algorithm_type == 'EXPLICIT_PERMUTATIONS':
            print(algo.permutations)
        else:
            raise ValueError('WRONG ALGORITHM !')

    test = raw_input('Save it ? ("y" to confirm)')
    if test == 'y':
        if len(cg.algorithms) != 1:
            raise ValueError('Multiple algorithms !')
        cg._algorithms = [ExplicitPermutationsAlgorithm(permutations=explicit_permutations)]
        cg_dict = cg.as_dict()
        f = open('../coordination_geometries_files_new/{}.json'.format(cg_symbol), 'w')
        json.dump(cg_dict, f)
        f.close()

        f = open('../coordination_geometries_files/{}.json'.format(cg_symbol), 'r')
        old_dict = json.load(f)
        f.close()

        print('COMPARISON :')
        keys_new = cg_dict.keys()
        keys_new.sort()
        keys_old = old_dict.keys()
        keys_old.sort()
        print('keys new : ', keys_new)
        print('keys old : ', keys_old)
        print('values :')
        for key, val in cg_dict.items():
            print(key)
            print(val)

            print(old_dict[key])
            print('\n\n')
    #
    # unique_csms = np.unique(csms)
    # print(csms)
    # print(unique_csms)
    #
    # prev = -1.0
    # ii = 0
    # for un in unique_csms:
    #     if not np.allclose(un, prev):
    #         ii+=1
    #         print(np.array(perms[np.argwhere(csms == un)[0]]), un)
    #     prev = un
    # print('ARRAY : ')
    # ii = 0
    # for un in unique_csms:
    #     if not np.allclose(un, prev):
    #         ii+=1
    #         perm = np.array(perms[np.argwhere(csms == un)[0]])
    #         print('np.array([{ii},'.format(ii=perm[0]))
    #         for ind in perm[1:-1]:
    #             print('{ii},'.format(ii=ind))
    #         print('{ii}], np.int),'.format(ii=perm[-1]))
    #     prev = un
    # print(ii)

#
#
#
#
# __author__ = 'waroquiers'
#
#
# import numpy as np
#
# from chemenv.coordination_environments.coordination_geometry_finder import coordination_geometries, LocalGeometryFinder
# from CoordinationGeometryUtils import Plane, OPlane, vectorsToMatrix, matrixTimesVector, rotateCoords, site_is_cation
# from CoordinationGeometryUtils import collinear, anticlockwise_sort, ChemicalEnvironments, ChemEnvAnalyzer, VoronoiContainer
#
#
# cg = coordination_geometries()
#
# mympsymbol = 'SS:4'
#
#
# lgf = LocalGeometryFinder()
# lgf.setup_parameters(centering_type='central_site')
#
# lgfsafe = LocalGeometryFinder(permutations_safe_override=True)
# lgfsafe.setup_parameters(centering_type='central_site')
#
#
# lgf.setup_test_perfect_structure(symbol=mympsymbol, randomness=True, randomness_percentage=10, symbol_type='MP', random_indices=True)
# lgf.setup_structure(lgf.get_structure())
# lgfsafe.setup_structure(lgf.get_structure())
# lgfsafe.setup_random_indices_local_geometry(cg.get_geometry_from_mp_symbol(mympsymbol).coordination_number)
# lgf.setup_voronoi(distfactor=1.4, angfactor=0.3, voronoi_cutoff=45)
# lgfsafe.setup_voronoi(distfactor=1.4, angfactor=0.3, voronoi_cutoff=45)
# #lgf.set_pmg_structure(lgfsafe.get_pmg_structure())
# #lgfsafe.set_pmg_structure(lgfsafe.get_pmg_structure())
# #lgf.setup_random_indices_local_geometry(cg.get_geometry_from_mp_symbol(mympsymbol).coordination_number)
#
# (csm, perms, algo) = lgf.coordination_geometry_symmetry_measures(cg.get_geometry_from_mp_symbol(mympsymbol))
# print(lgf.coordination_geometry_symmetry_measures(lgf.cg.get_geometry_from_mp_symbol(mympsymbol)))
# print(lgfsafe.coordination_geometry_symmetry_measures(lgfsafe.cg.get_geometry_from_mp_symbol(mympsymbol)))
# test = lgf.get_coordination_symmetry_measures(0)
# testsafe = lgfsafe.get_coordination_symmetry_measures(0)
#
# print(test)
# print(testsafe)
#
# input()
#
# print(csm)
# myunique = np.unique(csm)
# print(myunique)
# prev = -1.0
# ii = 0
# for un in myunique:
#     if not np.allclose(un, prev):
#         ii+=1
#         print(np.array(perms[np.argwhere(csm == un)[0]]), un)
#     prev = un
# print('ARRAY : ')
# ii = 0
# for un in myunique:
#     if not np.allclose(un, prev):
#         ii+=1
#         perm = np.array(perms[np.argwhere(csm == un)[0]])
#         print('np.array([{ii},'.format(ii=perm[0]), end=' ')
#         for ind in perm[1:-1]:
#             print('{ii},'.format(ii=ind), end=' ')
#         print('{ii}], np.int),'.format(ii=perm[-1]))
#     prev = un
# print(ii)