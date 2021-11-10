# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Development script of the ChemEnv utility to get the explicit permutations for coordination environments identified
with the separation plane algorithms (typically with coordination numbers >= 6)
"""

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import AbstractGeometry
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from pymatgen.analysis.chemenv.utils.coordination_geometry_utils import Plane, collinear

import numpy as np
import itertools
import json

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

    newalgos = []

    ialgo = 1
    for sepplanealgo in cg._algorithms:
        print(f"In ialgo = {ialgo:d}/{len(cg._algorithms):d}")
        ialgo += 1
        if sepplanealgo.algorithm_type != "SEPARATION_PLANE":
            raise ValueError("Should all be separation plane")

        permsonfile = f"Permutations on file in this algorithm ({len(sepplanealgo._permutations):d}) "
        print(permsonfile)
        print(sepplanealgo._permutations)
        permutations = sepplanealgo.safe_separation_permutations(
            ordered_plane=sepplanealgo.ordered_plane, ordered_point_groups=sepplanealgo.ordered_point_groups
        )

        sepplanealgo._permutations = permutations

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
        for npoints in range(sepplanealgo.minimum_number_of_points, min(sepplanealgo.maximum_number_of_points, 4) + 1):
            if found:
                break
            for ipoints in itertools.combinations(sepplanealgo.plane_points, npoints):
                points_combination = [lgf.local_geometry.coords[ipoint] for ipoint in ipoints]
                if npoints == 2:
                    if collinear(
                        points_combination[0], points_combination[1], lgf.local_geometry.central_site, tolerance=0.25
                    ):
                        continue
                    local_plane = Plane.from_3points(
                        points_combination[0], points_combination[1], lgf.local_geometry.central_site
                    )
                    found = True
                    break
                elif npoints == 3:
                    if collinear(points_combination[0], points_combination[1], points_combination[2], tolerance=0.25):
                        continue
                    local_plane = Plane.from_3points(
                        points_combination[0], points_combination[1], points_combination[2]
                    )
                    found = True
                    break
                elif npoints > 3:
                    local_plane = Plane.from_npoints(points_combination, best_fit="least_square_distance")
                    found = True
                    break
                else:
                    raise ValueError("Wrong number of points to initialize separation plane")

        points_perfect = lgf.perfect_geometry.points_wocs_ctwocc()
        # Actual test of the permutations
        cgsm = lgf._cg_csm_separation_plane(
            coordination_geometry=cg,
            sepplane=sepplanealgo,
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
        csms_with_recorded_permutation = []
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

        print(permsonfile)
        print(f"Permutations found ({len(explicit_permutations):d}) : ")
        print(explicit_permutations)
        sepplanealgo.explicit_permutations = explicit_permutations
        newalgos.append(sepplanealgo)

    # Write update geometry file ?
    test = input('Save it ? ("y" to confirm)')
    if test == "y":
        cg._algorithms = newalgos
        cg_dict = cg.as_dict()
        f = open(f"../coordination_geometries_files_new/{cg_symbol}.json", "w")
        json.dump(cg_dict, f)
        f.close()
