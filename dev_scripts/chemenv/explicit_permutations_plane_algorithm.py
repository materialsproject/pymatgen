# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Development script of the ChemEnv utility to get the explicit permutations for coordination environments identified
with the separation plane algorithms (typically with coordination numbers >= 6)
"""

from __future__ import annotations

import itertools
import json

import numpy as np

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import (
    AllCoordinationGeometries,
)
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import (
    AbstractGeometry,
    LocalGeometryFinder,
)
from pymatgen.analysis.chemenv.utils.coordination_geometry_utils import Plane, collinear

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
        if algo.algorithm_type != "SEPARATION_PLANE":
            raise ValueError("WRONG ALGORITHM !")

    new_algos = []

    ialgo = 1
    for sep_plane_algo in cg._algorithms:
        print(f"In {ialgo = :d}/{len(cg._algorithms):d}")
        ialgo += 1
        if sep_plane_algo.algorithm_type != "SEPARATION_PLANE":
            raise ValueError("Should all be separation plane")

        perms_on_file = f"Permutations on file in this algorithm ({len(sep_plane_algo._permutations):d}) "
        print(perms_on_file)
        print(sep_plane_algo._permutations)
        permutations = sep_plane_algo.safe_separation_permutations(
            ordered_plane=sep_plane_algo.ordered_plane, ordered_point_groups=sep_plane_algo.ordered_point_groups
        )

        sep_plane_algo._permutations = permutations

        print(f"Test permutations ({len(permutations):d}) :")
        print(permutations)

        lgf = LocalGeometryFinder()
        lgf.setup_parameters(structure_refinement=lgf.STRUCTURE_REFINEMENT_NONE)
        lgf.setup_test_perfect_environment(
            cg_symbol, randomness=True, indices=range(cg.coordination_number), max_random_dist=0.05
        )

        lgf.perfect_geometry = AbstractGeometry.from_cg(cg=cg)

        # Setting up the plane of separation
        local_plane = None
        found = False
        for n_points in range(
            sep_plane_algo.minimum_number_of_points, min(sep_plane_algo.maximum_number_of_points, 4) + 1
        ):
            if found:
                break
            for ipoints in itertools.combinations(sep_plane_algo.plane_points, n_points):
                points_combination = [lgf.local_geometry.coords[ipoint] for ipoint in ipoints]
                if n_points == 2:
                    if collinear(
                        points_combination[0], points_combination[1], lgf.local_geometry.central_site, tolerance=0.25
                    ):
                        continue
                    local_plane = Plane.from_3points(
                        points_combination[0], points_combination[1], lgf.local_geometry.central_site
                    )
                    found = True
                    break
                elif n_points == 3:
                    if collinear(points_combination[0], points_combination[1], points_combination[2], tolerance=0.25):
                        continue
                    local_plane = Plane.from_3points(
                        points_combination[0], points_combination[1], points_combination[2]
                    )
                    found = True
                    break
                elif n_points > 3:
                    local_plane = Plane.from_npoints(points_combination, best_fit="least_square_distance")
                    found = True
                    break
                else:
                    raise ValueError("Wrong number of points to initialize separation plane")

        points_perfect = lgf.perfect_geometry.points_wocs_ctwocc()
        # Actual test of the permutations
        cgsm = lgf._cg_csm_separation_plane(
            coordination_geometry=cg,
            sepplane=sep_plane_algo,
            local_plane=local_plane,
            plane_separations=[],
            dist_tolerances=[0.05, 0.1, 0.2, 0.3],
            testing=True,
            points_perfect=points_perfect,
        )

        print(cgsm)
        if cgsm[0] is None:
            print("IS NONE !")
            input()
            continue

        csms, perms, algos, sep_perms = cgsm[0], cgsm[1], cgsm[2], cgsm[3]

        print("Continuous symmetry measures")
        print(csms)
        csms_with_recorded_permutation: list[float] = []
        explicit_permutations = []
        for icsm, csm in enumerate(csms):
            found = False
            for csm2 in csms_with_recorded_permutation:
                if np.isclose(csm, csm2, rtol=0.0, atol=1.0e-6):
                    found = True
                    break
            if not found:
                print(perms[icsm], csm)
                csms_with_recorded_permutation.append(csm)
                explicit_permutations.append(sep_perms[icsm])

        print(perms_on_file)
        print(f"Permutations found ({len(explicit_permutations):d}) : ")
        print(explicit_permutations)
        sep_plane_algo.explicit_permutations = explicit_permutations
        new_algos.append(sep_plane_algo)

    # Write update geometry file ?
    test = input('Save it ? ("y" to confirm)')
    if test == "y":
        cg._algorithms = new_algos
        cg_dict = cg.as_dict()
        with open(f"../coordination_geometries_files_new/{cg_symbol}.json", "w") as f:
            json.dump(cg_dict, f)
