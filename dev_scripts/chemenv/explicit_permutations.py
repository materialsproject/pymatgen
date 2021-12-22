# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Development script of the ChemEnv utility to get the explicit permutations for coordination environments identified
with the explicit permutations algorithms (typically with coordination numbers <= 6)
"""

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import AbstractGeometry
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import ExplicitPermutationsAlgorithm

import numpy as np
import itertools
import json
import os


class Algo:
    pass


if __name__ == "__main__":

    # Choose the geometry
    allcg = AllCoordinationGeometries()
    while True:
        cg_symbol = input("Enter symbol of the geometry for which you want to get the explicit permutations : ")
        try:
            cg = allcg[cg_symbol]
            break
        except LookupError:
            print("Wrong geometry, try again ...")
            continue

    # Check if the algorithm currently defined for this geometry corresponds to the explicit permutation algorithm
    for algo in cg.algorithms:
        if algo.algorithm_type != "EXPLICIT_PERMUTATIONS":
            raise ValueError("WRONG ALGORITHM !")

    algo = Algo()
    algo.permutations = []
    for perm in itertools.permutations(range(cg.coordination)):
        algo.permutations.append(perm)

    lgf = LocalGeometryFinder()
    lgf.setup_parameters(structure_refinement=lgf.STRUCTURE_REFINEMENT_NONE)
    lgf.setup_test_perfect_environment(cg_symbol, randomness=True, indices="ORDERED")

    lgf.perfect_geometry = AbstractGeometry.from_cg(cg=cg)

    points_perfect = lgf.perfect_geometry.points_wocs_ctwocc()
    res = lgf.coordination_geometry_symmetry_measures_standard(
        coordination_geometry=cg, algo=algo, points_perfect=points_perfect
    )
    (csms, perms, algos, local2perfect_maps, perfect2local_maps) = res

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

    print("Permutations found : ")
    print(explicit_permutations)

    print("Current algorithm(s) :")
    for algo in cg.algorithms:
        print(algo)
        if algo.algorithm_type == "EXPLICIT_PERMUTATIONS":
            print(algo.permutations)
        else:
            raise ValueError("WRONG ALGORITHM !")

    test = input('Save it ? ("y" to confirm)')
    if test == "y":
        if len(cg.algorithms) != 1:
            raise ValueError("Multiple algorithms !")
        cg._algorithms = [ExplicitPermutationsAlgorithm(permutations=explicit_permutations)]
        newgeom_dir = "new_geometry_files"
        if not os.path.exists(newgeom_dir):
            os.makedirs(newgeom_dir)
        f = open(f"{newgeom_dir}/{cg_symbol}.json", "w")
        json.dump(cg.as_dict(), f)
        f.close()
