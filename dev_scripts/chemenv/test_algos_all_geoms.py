# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Development script to test the algorithms of all the model coordination environments
"""

from __future__ import annotations

import itertools
from math import factorial
from random import shuffle

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import (
    AllCoordinationGeometries,
)
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import (
    AbstractGeometry,
    LocalGeometryFinder,
)

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"

if __name__ == "__main__":
    allcg = AllCoordinationGeometries()

    test = input('Standard ("s", all permutations for cn <= 6, 500 random permutations for cn > 6) or on demand')
    if test == "s":
        perms_def = "standard"
    elif test == "o":
        perms_def = "on_demand"
    else:
        try:
            nperms = int(test)
            perms_def = "ndefined"
        except Exception:
            perms_def = "on_demand"

    for coordination in range(1, 13):
        print(f"IN COORDINATION {coordination:d}")
        symbol_name_mapping = allcg.get_symbol_name_mapping(coordination=coordination)

        if perms_def == "standard":
            test = "500" if coordination > 6 else "all"
        elif perms_def == "ndefined":
            test = nperms  # type: ignore[assignment]
        else:
            test = input(
                "Enter if you want to test all possible permutations ('all' or 'a') or "
                "a given number of random permutations (i.e. '25')"
            )
        myindices = range(coordination)

        if test == "all" or test == "a":
            perms_type = "all"
            perms_iterator = itertools.permutations(myindices)
            nperms = factorial(coordination)
        else:
            perms_type = "explicit"
            try:
                nperms = int(test)
            except Exception:
                raise ValueError(f"Could not turn {test} into integer ...")
            perms_iterator = []  # type: ignore[assignment]
            for _ in range(nperms):
                shuffle(myindices)  # type: ignore[arg-type]
                perms_iterator.append(list(myindices))  # type: ignore[attr-defined]

        for cg_symbol, cg_name in symbol_name_mapping.items():
            cg = allcg[cg_symbol]
            if cg.deactivate:
                continue

            print(f"Testing {cg_symbol} ({cg_name})")

            cg = allcg[cg_symbol]
            if cg.points is None:
                continue

            lgf = LocalGeometryFinder()
            lgf.setup_parameters(structure_refinement=lgf.STRUCTURE_REFINEMENT_NONE)

            # Reinitialize the itertools permutations
            if perms_type == "all":
                perms_iterator = itertools.permutations(myindices)

            # Loop on the permutations
            iperm = 1
            for indices_perm in perms_iterator:
                lgf.setup_test_perfect_environment(
                    cg_symbol,
                    indices=indices_perm,
                    randomness=True,
                    max_random_dist=0.1,
                    random_translation=True,
                    random_rotation=True,
                    random_scale=True,
                )

                lgf.perfect_geometry = AbstractGeometry.from_cg(cg=cg)
                points_perfect = lgf.perfect_geometry.points_wocs_ctwocc()

                print(f"Perm # {iperm:d}/{nperms:d} : ", indices_perm)

                algos_results = []
                for algo in cg.algorithms:
                    if algo.algorithm_type == "EXPLICIT_PERMUTATIONS":
                        results = lgf.coordination_geometry_symmetry_measures(
                            coordination_geometry=cg, points_perfect=points_perfect
                        )
                        # raise ValueError('Do something for the explicit ones ... (these should anyway be by far ok!)')
                    else:
                        results = lgf.coordination_geometry_symmetry_measures_separation_plane(
                            coordination_geometry=cg, separation_plane_algo=algo, points_perfect=points_perfect
                        )
                    algos_results.append(min(results[0]))

                    if not min(results[0]) < 1.5:
                        print("Following is not close to 0.0 ...")
                        input(results)
                print("   => ", algos_results)
                iperm += 1
