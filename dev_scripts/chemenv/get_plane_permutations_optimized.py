"""
Development script of the ChemEnv utility to get the optimized explicit permutations for coordination environments
identified with the separation plane algorithms (typically with coordination numbers >= 6).
"""

from __future__ import annotations

import itertools
import json
import os
import time
from math import factorial
from optparse import OptionParser
from random import shuffle

import numpy as np
import tabulate

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import (
    AbstractGeometry,
    LocalGeometryFinder,
)
from pymatgen.analysis.chemenv.utils.coordination_geometry_utils import Plane

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"


# Printing functions depending on the printing volume option
def prt1(string, printing_volume):
    if printing_volume >= 1:
        print(string)


def prt2(string, printing_volume):
    if printing_volume >= 2:
        print(string)


# Iterator function for the random permutations
def random_permutations_iterator(initial_permutation, n_permutations):
    """
    It takes a list and returns an iterator that yields random permutations of that list.

    Args:
        initial_permutation: the initial permutation of the data
        n_permutations: the number of permutations to generate
    """
    for _ in range(n_permutations):
        shuffle(initial_permutation)
        yield initial_permutation


if __name__ == "__main__":
    # Parse command line options
    option_parser = OptionParser()
    option_parser.add_option(
        "-v",
        "--printing_volume",
        action="store",
        type=int,
        default=1,
        help="Printing volume (0, 1 or 2)",
        dest="printing_volume",
    )
    option_parser.add_option(
        "-p",
        "--permutations_setup",
        action="store",
        type=str,
        default="y50",
        help="Setup of the permutations :\n"
        ' - "all" for all possible permutations,\n'
        ' - a number for a given number of random permutations (i.e. "25"),\n'
        ' - a number preceded by "x" for a number of random permutations equal to the '
        "maximum number of explicit permutations from all algorithms multiplied by this "
        "number or\n"
        ' - a number preceded by "y" for a number of random permutations equal to the '
        "maximum number of explicit permutations from all algorithms multiplied by this "
        "number except if this number exceeds the number of all possible permutations, in"
        "which case all possible permutations are used",
        dest="permutations_setup",
    )
    options, args = option_parser.parse_args()

    # Get the printing volume and permutations setup options
    printing_volume = options.printing_volume
    permutations_setup = options.permutations_setup
    n_perm_factor = None
    if permutations_setup == "all":
        permutations_setup_type = "all"
        n_permutations = None
    elif permutations_setup[0] == "x":
        permutations_setup_type = "x"
        n_permutations = None
        try:
            n_perm_factor = int(permutations_setup[1:])
        except Exception:
            raise ValueError("Wrong command line option for permutations_setup")
    elif permutations_setup[0] == "y":
        permutations_setup_type = "y"
        n_permutations = None
        try:
            n_perm_factor = int(permutations_setup[1:])
        except Exception:
            raise ValueError("Wrong command line option for permutations_setup")
    else:
        permutations_setup_type = "n"
        try:
            n_permutations = int(permutations_setup)
        except Exception:
            raise ValueError("Wrong command line option for permutations_setup")

    # Class containing all the coordination geometries
    all_cg = AllCoordinationGeometries()

    sep_plane_cgs = []
    for coordination in range(1, 21):
        symbol_name_mapping = all_cg.get_symbol_name_mapping(coordination=coordination)
        for symbol in symbol_name_mapping:
            cg = all_cg[symbol]
            if cg.points is None:
                continue
            if cg.algorithms[0].algorithm_type != "EXPLICIT_PERMUTATIONS":
                sep_plane_cgs.append(symbol)
                continue
    n_cols = 5
    n_lines = int(np.ceil(float(len(sep_plane_cgs)) / n_cols))
    sep_plane_cgs_grid = []
    for _ in range(n_lines):
        sep_plane_cgs_grid.append([""] * n_cols)
    for line_idx in range(n_lines):
        for col_idx in range(n_cols):
            ii = line_idx * n_cols + col_idx
            if ii >= len(sep_plane_cgs):
                break
            sep_plane_cgs_grid[line_idx][col_idx] = sep_plane_cgs[ii]

    while True:
        # Printing all symbols
        print("Coordination geometries using a separation plane algorithm :")
        print(tabulate.tabulate(sep_plane_cgs_grid, tablefmt="grid"))
        print()

        # Define the coordination geometry
        cg_symbol = input(
            'Enter symbol of the geometry for which you want to get the optimized permutations or "q" to quit : '
        )
        if cg_symbol == "q":
            break
        if cg_symbol not in sep_plane_cgs:
            print("Wrong geometry, try again ...")
            continue

        cg = all_cg[cg_symbol]

        print(f"Getting explicit permutations for geometry {cg.name!r} (symbol : {cg_symbol!r})\n")

        # Setup of the local geometry finder
        lgf = LocalGeometryFinder()
        lgf.setup_parameters(structure_refinement=lgf.STRUCTURE_REFINEMENT_NONE)

        # Setup the random environment
        lgf.setup_test_perfect_environment(
            cg_symbol, randomness=True, indices=range(cg.coordination_number), max_random_dist=0.05
        )
        lgf.perfect_geometry = AbstractGeometry.from_cg(cg=cg)
        points_perfect = lgf.perfect_geometry.points_wcs_ctwcc()

        # 1. Check the algorithms defined for this coordination geometry and get the explicit permutations
        original_n_explicit_perms = []
        original_n_explicit_optimized_perms = []
        for idx, algo in enumerate(cg.algorithms):
            algo._permutations = algo.explicit_permutations
            algo.minimum_number_of_points = 4
            if algo.algorithm_type == "EXPLICIT_PERMUTATIONS":
                raise ValueError("Do something for the explicit ones ... (these should anyway be by far ok!)")
            if algo.explicit_optimized_permutations is None:
                eop = "no"
            else:
                eop = str(len(algo.explicit_optimized_permutations))
            print(
                f"For {idx=}, plane_points are "
                f"[{', '.join(map(str, algo.plane_points))}], "
                f"side_0 is [{', '.join(map(str, algo.point_groups[0]))}] and "
                f"side_1 is [{', '.join(map(str, algo.point_groups[1]))}]."
            )
            original_n_explicit_perms.append(len(algo.explicit_permutations))
            original_n_explicit_optimized_perms.append(eop)
            print(
                f"  For this algorithm, there are {eop} optimized permutations and "
                f"{len(algo.explicit_permutations)} explicit permutations"
            )
            if algo.other_plane_points is None:
                input("Multiplicity and other plane points is not defined for this algorithm !")

            # Setup of safe permutations
            permutations = algo.safe_separation_permutations(
                ordered_plane=algo.ordered_plane, ordered_point_groups=algo.ordered_point_groups
            )
            algo._permutations = permutations
            print(f"Safe permutations found ({len(permutations)})")

            # Definition of the facets
            all_planes_point_indices = [algo.plane_points]
            if algo.other_plane_points is not None:
                all_planes_point_indices.extend(algo.other_plane_points)

            # Loop on the facets
            explicit_permutations_per_plane = []
            for idx_plane, plane_point_indices in enumerate(all_planes_point_indices):
                prt1(
                    string=f"In plane {idx_plane} ({'-'.join(str(pp) for pp in plane_point_indices)})",
                    printing_volume=printing_volume,
                )

                points_combination = [lgf.local_geometry._coords[ii] for ii in plane_point_indices]
                local_plane = Plane.from_npoints(points_combination, best_fit="least_square_distance")

                # Actual test of the permutations
                csms, perms, algos, sep_perms = lgf._cg_csm_separation_plane(
                    coordination_geometry=cg,
                    sep_plane=algo,
                    local_plane=local_plane,
                    plane_separations=[],
                    dist_tolerances=[0.05, 0.1, 0.2, 0.3, 0.5],
                    testing=True,
                    points_perfect=points_perfect,
                )

                sym_measures = [c["symmetry_measure"] for c in csms]
                prt1(string="Continuous symmetry measures", printing_volume=printing_volume)
                prt1(string=sym_measures, printing_volume=printing_volume)
                csms_with_recorded_permutation = []  # type: ignore
                explicit_permutations = []
                for icsm, csm in enumerate(csms):
                    found = False
                    for csm2 in csms_with_recorded_permutation:
                        if np.isclose(csm["symmetry_measure"], csm2["symmetry_measure"], rtol=0.0):
                            found = True
                            break
                    if not found:
                        prt1(
                            string=f" permutation {'-'.join(map(str, sep_perms[icsm]))} : {csm['symmetry_measure']}",
                            printing_volume=printing_volume,
                        )
                        csms_with_recorded_permutation.append(csm)
                        explicit_permutations.append(tuple(sep_perms[icsm]))

                prt1(string=explicit_permutations, printing_volume=printing_volume)
                explicit_permutations_per_plane.append(set(explicit_permutations))
                prt1(string="", printing_volume=printing_volume)
            # Check that the explicit permutations found are the same for each plane
            for ip1 in range(len(explicit_permutations_per_plane) - 1):
                ep_p1 = explicit_permutations_per_plane[ip1]
                for ip2 in range(1, len(explicit_permutations_per_plane)):
                    ep_p2 = explicit_permutations_per_plane[ip2]
                    if len(ep_p1 & ep_p2) != len(ep_p1):
                        print("Explicit permutations per plane :")
                        for eppp in explicit_permutations_per_plane:
                            print(eppp)
                        raise ValueError("Explicit permutations different from one plane to another !")
            algo.explicit_permutations = [list(perm) for perm in list(explicit_permutations_per_plane[0])]
            algo.explicit_permutations.sort()
            algo.explicit_permutations = np.array(algo.explicit_permutations)
            print(f"Explicit permutations found ({len(algo.explicit_permutations)})")
            print(algo.explicit_permutations)
            print()
            # Setup the permutations for the next optimization
            algo._permutations = algo.explicit_permutations

        while True:
            test = input(
                f"Get the explicit optimized permutations for geometry {cg.name!r} (symbol : "
                f'{cg_symbol!r}) ? ("y" to confirm, "q" to quit)\n'
            )
            if test not in ["y", "q"]:
                print("Wrong key, try again")
                continue
            if test == "y":
                break
            elif test == "q":
                raise SystemExit(0)
        # 2. Optimization of the permutations
        print(f"Getting explicit optimized permutations for geometry {cg.name!r} (symbol : {cg_symbol!r})\n")
        perms_used_algos: list[dict] = [{} for _ in cg.algorithms]

        # Loop on algorithms
        for idx, algo in enumerate(cg.algorithms):
            perms_used: dict[tuple, int] = {}
            print(
                f"In {idx=} (plane_points : "
                f"[{', '.join(map(str, algo.plane_points))}], "
                f"side_0 : [{', '.join(map(str, algo.point_groups[0]))}] and "
                f"side_1 : [{', '.join(map(str, algo.point_groups[1]))}])"
            )
            if algo.algorithm_type == "EXPLICIT_PERMUTATIONS":
                raise ValueError("Do something for the explicit ones ... (these should anyway be by far ok!)")

            # Definition of the facets
            all_planes_point_indices = [algo.plane_points]
            if algo.other_plane_points is not None:
                all_planes_point_indices.extend(algo.other_plane_points)

            # Setup of the permutations to be used for this algorithm

            indices = list(range(cg.coordination_number))
            if permutations_setup_type == "all":
                perms_iterator = itertools.permutations(indices)
                n_permutations = factorial(cg.coordination_number)
            elif permutations_setup_type == "n":
                if n_permutations and n_permutations >= factorial(cg.coordination_number):
                    perms_iterator = itertools.permutations(indices)
                    n_permutations = factorial(cg.coordination_number)
                else:
                    perms_iterator = random_permutations_iterator(
                        initial_permutation=indices, n_permutations=n_permutations
                    )
            elif permutations_setup_type in ["x", "y"] and n_perm_factor:
                n_permutations = n_perm_factor * len(algo.explicit_permutations)
                if permutations_setup_type == "y" and n_permutations >= factorial(cg.coordination_number):
                    perms_iterator = itertools.permutations(indices)
                    n_permutations = factorial(cg.coordination_number)
                else:
                    perms_iterator = random_permutations_iterator(
                        initial_permutation=indices, n_permutations=n_permutations
                    )
            else:
                raise ValueError("Permutation setup not allowed ...")

            # Loop on permutations
            idx_perm = 1
            t0 = time.process_time()
            time_left: float | str = "Unknown"
            for indices_perm in perms_iterator:
                prt1(
                    string=f"Perm # {idx_perm}/{n_permutations} : "
                    f"{'-'.join(map(str, indices_perm))} "
                    f"(est. rem. time : {time_left} sec)",
                    printing_volume=printing_volume,
                )
                # Setup of the local and perfect geometries
                lgf.setup_test_perfect_environment(
                    cg_symbol, indices=indices_perm, randomness=True, max_random_dist=0.02, random_rotation=True
                )
                lgf.perfect_geometry = AbstractGeometry.from_cg(cg=cg)
                points_perfect = lgf.perfect_geometry.points_wcs_ctwcc()

                # Loop on the facets
                separation_permutations = []
                for idx_plane, plane_point_indices in enumerate(all_planes_point_indices):
                    prt2(
                        string=f"In plane {idx_plane} ({'-'.join(str(pp) for pp in plane_point_indices)})",
                        printing_volume=printing_volume,
                    )

                    # Setup of separation plane
                    perm_plane_points_indices = [indices_perm.index(plane_point) for plane_point in plane_point_indices]
                    shuffle(perm_plane_points_indices)
                    points_combination = [lgf.local_geometry._coords[ii] for ii in perm_plane_points_indices]
                    local_plane = Plane.from_npoints(points_combination, best_fit="least_square_distance")

                    # Get the results for this algorithm and plane
                    csms, perms, algos, sep_perms = lgf._cg_csm_separation_plane(
                        coordination_geometry=cg,
                        sep_plane=algo,
                        local_plane=local_plane,
                        plane_separations=[],
                        dist_tolerances=[0.05, 0.1, 0.2, 0.3, 0.5],
                        testing=True,
                        points_perfect=points_perfect,
                    )

                    sym_measures = [c["symmetry_measure"] for c in csms]
                    imin = np.argmin(sym_measures)
                    min_csm = min(sym_measures)
                    if not min_csm < 1.0:
                        print("Following is not close enough to 0 ...")
                        input(sym_measures)
                    mincsm_indices = []
                    for icsm, csm in enumerate(sym_measures):
                        if np.isclose(min_csm, csm, rtol=0.0):
                            mincsm_indices.append(icsm)
                    this_plane_sep_perm = tuple(sep_perms[imin])
                    prt2(
                        string=f"  permutation {'-'.join(map(str, this_plane_sep_perm))} "
                        f"gives csm={sym_measures[imin]:.6f}",
                        printing_volume=printing_volume,
                    )

                    separation_permutations.append(this_plane_sep_perm)
                for some_perm in separation_permutations:
                    if some_perm in perms_used:
                        perms_used[some_perm] += 1
                    else:
                        perms_used[some_perm] = 1
                tcurrent = time.process_time()
                time_left = (n_permutations - idx_perm) * (tcurrent - t0) / idx_perm  # type: ignore
                time_left = f"{time_left:.1f}"
                idx_perm += 1
            print(
                f"Optimized permutations {len(perms_used)}/{len(algo.permutations)}"
                f"(old : {original_n_explicit_optimized_perms[idx]}/{original_n_explicit_perms[idx]}) : "
            )
            for perm, number in perms_used.items():
                print(f" - permutation {'-'.join(map(str, perm))} : {number}")
            print(
                f"For {idx=} (plane_points : [{', '.join(map(str, algo.plane_points))}], "
                f"side_0 : [{', '.join(map(str, algo.point_groups[0]))}] and "
                f"side_1 : [{', '.join(map(str, algo.point_groups[1]))}]),\n"
                f"Optimized perturbations {len(perms_used)}/{len(algo.permutations)} (old : "
                f"{original_n_explicit_optimized_perms[idx]}/{original_n_explicit_perms[idx]}) are :"
            )
            # print(f"Optimized permutations ({len(perms_used)}/{len(algo.permutations)}) : ")
            explicit_optimized_permutations = [list(perm) for perm in perms_used]
            explicit_optimized_permutations.sort()
            print(explicit_optimized_permutations)
            print()
            test = input(f'Set optimized permutations for algorithm {idx} ? ("y" to confirm)')
            if test == "y":
                algo.explicit_optimized_permutations = np.array(explicit_optimized_permutations)

        test = input(
            f"Save coordination geometry {cg.name!r} (symbol {cg_symbol!r}) and new explicit and optimized "
            'permutations ? ("y" to confirm)'
        )
        if test == "y":
            new_geom_dir = "new_geometry_files"
            os.makedirs(new_geom_dir, exist_ok=True)
            with open(f"{new_geom_dir}/{cg_symbol}.json", "w") as file:
                json.dump(cg.as_dict(), file)
