"""Development script to test the algorithms of a given model coordination environments."""

from __future__ import annotations

import itertools
import time
from math import factorial
from random import shuffle

import numpy as np

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
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
    all_cg = AllCoordinationGeometries()

    while True:
        cg_symbol = input("Enter symbol of the geometry for which you want to get the explicit permutations : ")
        try:
            cg = all_cg[cg_symbol]
            break
        except LookupError:
            print("Wrong geometry, try again ...")
            continue

    lgf = LocalGeometryFinder()
    lgf.setup_parameters(structure_refinement=lgf.STRUCTURE_REFINEMENT_NONE)

    my_indices = range(cg.coordination_number)

    test = input(
        'Enter if you want to test all possible permutations ("all" or "a") or a given number of '
        'random permutations (i.e. "25")'
    )

    if test in ("all", "a"):
        perms_iterator = itertools.permutations(my_indices)
        n_perms = factorial(cg.coordination_number)
    else:
        try:
            n_perms = int(test)
        except Exception:
            raise ValueError(f"Could not turn {test} into integer ...")
        perms_iterator = []  # type: ignore[assignment]
        for _ in range(n_perms):
            shuffle(my_indices)  # type: ignore[arg-type]
            perms_iterator.append(list(my_indices))  # type: ignore[attr-defined]

    idx_perm = 1
    t1 = time.perf_counter()
    for indices_perm in perms_iterator:
        lgf.setup_test_perfect_environment(cg_symbol, indices=indices_perm)

        lgf.perfect_geometry = AbstractGeometry.from_cg(cg=cg)
        points_perfect = lgf.perfect_geometry.points_wocs_ctwocc()

        print(f"Perm # {idx_perm}/{n_perms} : ", indices_perm)

        algos_results = []
        for algo in cg.algorithms:
            print(algo)
            if algo.algorithm_type == "EXPLICIT_PERMUTATIONS":
                raise ValueError("Do something for the explicit ones ... (these should anyway be by far ok!)")

            results = lgf.coordination_geometry_symmetry_measures_separation_plane(
                coordination_geometry=cg,
                separation_plane_algo=algo,
                tested_permutations=False,
                points_perfect=points_perfect,
            )
            print("Number of permutations tested : ", len(results[0]))
            algos_results.append(min(results[0]))

            if not np.isclose(min(results[0]), 0):
                print("Following is not 0 ...")
                input(results)
        print("   => ", algos_results)
        idx_perm += 1
    t2 = time.perf_counter()
    print(
        f"Time to test {n_perms} permutations for geometry {cg.name!r} (symbol {cg_symbol!r}) : {t2 - t1:.2f} seconds"
    )
