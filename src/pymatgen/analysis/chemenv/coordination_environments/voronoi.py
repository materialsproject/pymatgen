"""This module contains the object used to describe the possible bonded atoms based on a Voronoi analysis."""

from __future__ import annotations

import itertools
import logging
import time
from typing import TYPE_CHECKING, Any

import matplotlib.pyplot as plt
import numpy as np
from monty.json import MSONable
from scipy.spatial import Voronoi

from pymatgen.analysis.chemenv.utils.coordination_geometry_utils import (
    get_lower_and_upper_f,
    rectangle_surface_intersection,
    solid_angle,
)
from pymatgen.analysis.chemenv.utils.defs_utils import AdditionalConditions
from pymatgen.analysis.chemenv.utils.math_utils import normal_cdf_step
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure

if TYPE_CHECKING:
    from typing_extensions import Self

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"

logger = logging.getLogger(__name__)


def from_bson_voronoi_list2(bson_nb_voro_list2: list[PeriodicSite], structure: Structure):
    """Get the voronoi_list needed for the VoronoiContainer object from a BSON-encoded voronoi_list."""
    # Pre-fetch per-site data once to avoid repeated attribute lookups
    lattice = structure._lattice
    species_list = [s._species for s in structure]
    frac_coords_list = [s.frac_coords for s in structure]  # numpy arrays
    props_list = [s.properties for s in structure]

    voronoi_list: list[list[dict] | None] = [None] * len(bson_nb_voro_list2)

    PSite = PeriodicSite  # local bind for speed

    for i, voro in enumerate(bson_nb_voro_list2):
        if voro is None or voro == "None":
            continue

        # Build the site list in one pass; avoid touching structure inside the loop
        site_entries: list[dict[str, Any]] = []
        append = site_entries.append  # local bind
        for psd, dct in voro:
            # psd == [index_in_structure, image_shift_frac]
            sidx = dct["index"]  # index of the symmetry-equivalent site in the original structure
            # Compose fractional coords in the target image cell
            # (frac_coords_list[sidx] is ndarray; psd[1] is a short list of floats)
            fcoords = frac_coords_list[sidx] + psd[1]

            # Avoid mutating the original dict object if it might be reused elsewhere
            nd = dct.copy()
            nd["site"] = PSite(
                species_list[sidx],
                fcoords,
                lattice,
                properties=props_list[sidx],
            )
            append(nd)

        voronoi_list[i] = site_entries

    return voronoi_list


class DetailedVoronoiContainer(MSONable):
    """Store the full Voronoi of a given structure."""

    AC = AdditionalConditions()
    default_voronoi_cutoff = 10.0
    default_normalized_distance_tolerance = 1e-5
    default_normalized_angle_tolerance = 1e-3

    def __init__(
        self,
        structure=None,
        voronoi_list2=None,
        voronoi_cutoff=default_voronoi_cutoff,
        isites=None,
        normalized_distance_tolerance=default_normalized_distance_tolerance,
        normalized_angle_tolerance=default_normalized_angle_tolerance,
        additional_conditions=None,
        valences=None,
        maximum_distance_factor=None,
        minimum_angle_factor=None,
    ):
        """
        Constructor for the VoronoiContainer object. Either a structure is given, in which case the Voronoi is
        computed, or the different components of the VoronoiContainer are given (used in the from_dict method).

        Args:
            structure: Structure for which the Voronoi is computed.
            voronoi_list2: List of voronoi polyhedrons for each site.
            voronoi_cutoff: cutoff used for the voronoi.
            isites: indices of sites for which the Voronoi has to be computed.
            normalized_distance_tolerance: Tolerance for two normalized distances to be considered equal.
            normalized_angle_tolerance:Tolerance for two normalized angles to be considered equal.
            additional_conditions: Additional conditions to be used.
            valences: Valences of all the sites in the structure (used when additional conditions require it).
            maximum_distance_factor: The maximum distance factor to be considered.
            minimum_angle_factor: The minimum angle factor to be considered.

        Raises:
            RuntimeError if the Voronoi cannot be constructed.
        """
        self.normalized_distance_tolerance = normalized_distance_tolerance
        self.normalized_angle_tolerance = normalized_angle_tolerance
        if additional_conditions is None:
            self.additional_conditions = [self.AC.NONE, self.AC.ONLY_ACB]
        else:
            self.additional_conditions = additional_conditions
        self.valences = valences
        self.maximum_distance_factor = maximum_distance_factor
        self.minimum_angle_factor = minimum_angle_factor
        indices = list(range(len(structure))) if isites is None else isites
        self.structure = structure
        logger.debug("Setting Voronoi list")
        if voronoi_list2 is not None:
            self.voronoi_list2 = voronoi_list2
        else:
            self.setup_voronoi_list(indices=indices, voronoi_cutoff=voronoi_cutoff)
        logger.debug("Setting neighbors distances and angles")
        t1 = time.process_time()
        self.setup_neighbors_distances_and_angles(indices=indices)
        t2 = time.process_time()
        logger.debug(f"Neighbors distances and angles set up in {t2 - t1:.2f} seconds")

    def setup_voronoi_list(self, indices, voronoi_cutoff):
        """Set up of the voronoi list of neighbors by calling qhull.

        Args:
            indices: indices of the sites for which the Voronoi is needed.
            voronoi_cutoff: Voronoi cutoff for the search of neighbors.

        Raises:
            RuntimeError: If an infinite vertex is found in the voronoi construction.
        """
        self.voronoi_list2 = [None] * len(self.structure)
        self.voronoi_list_coords = [None] * len(self.structure)
        logger.debug("Getting all neighbors in structure")
        struct_neighbors = self.structure.get_all_neighbors(voronoi_cutoff, include_index=True)
        size_neighbors = [(not len(neigh) > 3) for neigh in struct_neighbors]
        if np.any(size_neighbors):
            logger.debug("Please consider increasing voronoi_distance_cutoff")
        t1 = time.process_time()
        logger.debug("Setting up Voronoi list :")
        for jj, isite in enumerate(indices, start=1):
            logger.debug(f"  - Voronoi analysis for site #{isite} ({jj}/{len(indices)})")
            site = self.structure[isite]
            neighbors1 = [(site, 0.0, isite)]
            neighbors1.extend(struct_neighbors[isite])
            distances = [i[1] for i in sorted(neighbors1, key=lambda s: s[1])]
            neighbors = [i[0] for i in sorted(neighbors1, key=lambda s: s[1])]
            qvoronoi_input = [s.coords for s in neighbors]
            voro = Voronoi(points=qvoronoi_input, qhull_options="o Fv")
            all_vertices = voro.vertices

            results2 = []
            max_angle = 0.0
            min_dist = 10000.0
            for idx, ridge_points in enumerate(voro.ridge_points):
                if 0 in ridge_points:
                    ridge_vertices_indices = voro.ridge_vertices[idx]
                    if -1 in ridge_vertices_indices:
                        raise RuntimeError(
                            "This structure is pathological, infinite vertex in the voronoi construction"
                        )

                    ridge_point2 = max(ridge_points)
                    facets = [all_vertices[i] for i in ridge_vertices_indices]
                    sa = solid_angle(site.coords, facets)
                    max_angle = max([sa, max_angle])

                    min_dist = min([min_dist, distances[ridge_point2]])
                    for iii, sss in enumerate(self.structure):
                        if neighbors[ridge_point2].is_periodic_image(sss, tolerance=1e-6):
                            idx = iii
                            break
                    results2.append(
                        {
                            "site": neighbors[ridge_point2],
                            "angle": sa,
                            "distance": distances[ridge_point2],
                            "index": idx,
                        }
                    )
            for dd in results2:
                dd["normalized_angle"] = dd["angle"] / max_angle
                dd["normalized_distance"] = dd["distance"] / min_dist
            self.voronoi_list2[isite] = results2
            self.voronoi_list_coords[isite] = np.array([dd["site"].coords for dd in results2])
        t2 = time.process_time()
        logger.debug(f"Voronoi list set up in {t2 - t1:.2f} seconds")

    def setup_neighbors_distances_and_angles(self, indices):
        """Initialize the angle and distance separations.

        Args:
            indices: Indices of the sites for which the Voronoi is needed.
        """
        nsites = len(self.structure)
        self.neighbors_distances = [None] * nsites
        self.neighbors_normalized_distances = [None] * nsites
        self.neighbors_angles = [None] * nsites
        self.neighbors_normalized_angles = [None] * nsites

        nd_tol = self.normalized_distance_tolerance
        na_tol = self.normalized_angle_tolerance
        max_dfact = self.maximum_distance_factor
        min_afact = self.minimum_angle_factor
        default_cutoff = self.default_voronoi_cutoff

        for site_idx in indices:
            results = self.voronoi_list2[site_idx]
            if results is None:
                continue

            # Prepare arrays once
            ndists = np.array([nb["normalized_distance"] for nb in results], dtype=float)
            dists = np.array([nb["distance"] for nb in results], dtype=float)
            order_d = np.argsort(ndists)

            # Initialize group containers
            nn_dist_groups = []  # normalized distance groups (dicts)
            dist_groups = []  # real distance groups (dicts)

            # Seed first group
            first_idx = int(order_d[0])
            current_min_nd = ndists[first_idx]
            current_max_nd = ndists[first_idx]
            current_min_d = dists[first_idx]
            current_max_d = dists[first_idx]
            nb_indices_set = {first_idx}  # all indices up to current group
            dnb_set = {first_idx}  # indices equal (within tol) to current group's edge

            # Scan sorted normalized distances
            for id_np in order_d:  # already ascending
                idist = int(id_np)  # ensure Python int (BSON-safe)
                wd = ndists[idist]

                # early stop if a max distance factor is set
                if (max_dfact is not None) and (wd > max_dfact):
                    # close current group
                    nn_dist_groups.append(
                        {
                            "min": current_min_nd,
                            "max": current_max_nd,
                            "nb_indices": list(nb_indices_set),
                            "dnb_indices": list(dnb_set),
                        }
                    )
                    dist_groups.append(
                        {
                            "min": current_min_d,
                            "max": current_max_d,
                            "nb_indices": list(nb_indices_set),
                            "dnb_indices": list(dnb_set),
                        }
                    )
                    break

                if np.isclose(wd, current_max_nd, rtol=0.0, atol=nd_tol):
                    # still in the current plateau => extend max and dnb_set
                    current_max_nd = wd
                    current_max_d = dists[idist]
                    dnb_set.add(idist)
                else:
                    # finalize previous group
                    nn_dist_groups.append(
                        {
                            "min": current_min_nd,
                            "max": current_max_nd,
                            "nb_indices": list(nb_indices_set),
                            "dnb_indices": list(dnb_set),
                        }
                    )
                    dist_groups.append(
                        {
                            "min": current_min_d,
                            "max": current_max_d,
                            "nb_indices": list(nb_indices_set),
                            "dnb_indices": list(dnb_set),
                        }
                    )
                    # start new group
                    current_min_nd = current_max_nd = wd
                    current_min_d = current_max_d = dists[idist]
                    dnb_set = {idist}

                nb_indices_set.add(idist)
            else:
                # loop ended normally: close last group
                nn_dist_groups.append(
                    {
                        "min": current_min_nd,
                        "max": current_max_nd,
                        "nb_indices": list(nb_indices_set),
                        "dnb_indices": list(dnb_set),
                    }
                )
                dist_groups.append(
                    {
                        "min": current_min_d,
                        "max": current_max_d,
                        "nb_indices": list(nb_indices_set),
                        "dnb_indices": list(dnb_set),
                    }
                )

            # Fill "next" for distance groups
            for i in range(len(dist_groups) - 1):
                dist_groups[i]["next"] = dist_groups[i + 1]["min"]
                nn_dist_groups[i]["next"] = nn_dist_groups[i + 1]["min"]

            # Tail "next" based on cutoff logic
            if max_dfact is not None:
                dfact = max_dfact
            else:
                # avoid repeated indexing
                first_min_dist = dist_groups[0]["min"]
                dfact = default_cutoff / first_min_dist
            nn_dist_groups[-1]["next"] = dfact
            dist_groups[-1]["next"] = dfact * dist_groups[0]["min"]

            self.neighbors_normalized_distances[site_idx] = nn_dist_groups
            self.neighbors_distances[site_idx] = dist_groups

            nangs = np.array([nb["normalized_angle"] for nb in results], dtype=float)
            angs = np.array([nb["angle"] for nb in results], dtype=float)
            order_a = np.argsort(nangs)[::-1]  # descending

            nn_ang_groups = []
            ang_groups = []

            first_a_idx = int(order_a[0])
            current_max_na = nangs[first_a_idx]
            current_min_na = nangs[first_a_idx]
            current_max_a = angs[first_a_idx]
            current_min_a = angs[first_a_idx]
            nb_indices_set = {first_a_idx}
            dnb_set = {first_a_idx}

            for ia_np in order_a:
                iang = int(ia_np)
                wa = nangs[iang]

                # early stop if a minimum angle factor is set
                if (min_afact is not None) and (wa < min_afact):
                    nn_ang_groups.append(
                        {
                            "max": current_max_na,
                            "min": current_min_na,
                            "nb_indices": list(nb_indices_set),
                            "dnb_indices": list(dnb_set),
                        }
                    )
                    ang_groups.append(
                        {
                            "max": current_max_a,
                            "min": current_min_a,
                            "nb_indices": list(nb_indices_set),
                            "dnb_indices": list(dnb_set),
                        }
                    )
                    break

                if np.isclose(wa, current_min_na, rtol=0.0, atol=na_tol):
                    # staying on the current lower edge (since we traverse from high to low)
                    current_min_na = wa
                    current_min_a = angs[iang]
                    dnb_set.add(iang)
                else:
                    # finalize current group
                    nn_ang_groups.append(
                        {
                            "max": current_max_na,
                            "min": current_min_na,
                            "nb_indices": list(nb_indices_set),
                            "dnb_indices": list(dnb_set),
                        }
                    )
                    ang_groups.append(
                        {
                            "max": current_max_a,
                            "min": current_min_a,
                            "nb_indices": list(nb_indices_set),
                            "dnb_indices": list(dnb_set),
                        }
                    )
                    # start new group anchored at this lower value
                    current_max_na = current_min_na = wa
                    current_max_a = current_min_a = angs[iang]
                    dnb_set = {iang}

                nb_indices_set.add(iang)
            else:
                # loop ended normally: close last group
                nn_ang_groups.append(
                    {
                        "max": current_max_na,
                        "min": current_min_na,
                        "nb_indices": list(nb_indices_set),
                        "dnb_indices": list(dnb_set),
                    }
                )
                ang_groups.append(
                    {
                        "max": current_max_a,
                        "min": current_min_a,
                        "nb_indices": list(nb_indices_set),
                        "dnb_indices": list(dnb_set),
                    }
                )

            # Fill "next" for angle groups
            for i in range(len(ang_groups) - 1):
                ang_groups[i]["next"] = ang_groups[i + 1]["max"]
                nn_ang_groups[i]["next"] = nn_ang_groups[i + 1]["max"]

            afact_tail = min_afact if (min_afact is not None) else 0.0
            nn_ang_groups[-1]["next"] = afact_tail
            # follow existing convention for real angle "next"
            ang_groups[-1]["next"] = afact_tail * ang_groups[0]["max"]

            self.neighbors_normalized_angles[site_idx] = nn_ang_groups
            self.neighbors_angles[site_idx] = ang_groups

    def _precompute_additional_conditions(self, ivoronoi, voronoi, valences):
        additional_conditions = {}
        check = self.AC.check_condition
        structure = self.structure
        nb_indices = [vals["index"] for _, vals in voronoi]

        for ac in self.additional_conditions:
            additional_conditions[ac] = [
                check(
                    condition=ac,
                    structure=structure,
                    parameters={
                        "valences": valences,
                        "neighbor_index": nb_idx,
                        "site_index": ivoronoi,
                    },
                )
                for nb_idx in nb_indices
            ]
        return additional_conditions

    def _precompute_angle_conditions(self, ivoronoi, voronoi):
        angle_conditions = []
        # Vectorize the list of voronoi normalized angles once
        v_angles = np.fromiter((vals["normalized_angle"] for _, vals in voronoi), dtype=float)
        atol = self.normalized_angle_tolerance / 2.0

        for ap_dict in self.neighbors_normalized_angles[ivoronoi]:
            ap = ap_dict["max"]
            # (v_angles >= ap) OR close to ap within tolerance
            cond = (v_angles >= ap) | np.isclose(v_angles, ap, rtol=0.0, atol=atol)
            angle_conditions.append(cond.tolist())  # ensure plain Python bools

        return angle_conditions

    def neighbors_surfaces(self, isite, surface_calculation_type=None, max_dist=2.0):
        """Get the different surfaces corresponding to the different distance-angle cutoffs for a given site.

        Args:
            isite: Index of the site
            surface_calculation_type: How to compute the surface.
            max_dist: The maximum distance factor to be considered.

        Returns:
            Surfaces for each distance-angle cutoff.
        """
        if self.voronoi_list2[isite] is None:
            return None
        bounds_and_limits = self.voronoi_parameters_bounds_and_limits(isite, surface_calculation_type, max_dist)
        distance_bounds = bounds_and_limits["distance_bounds"]
        angle_bounds = bounds_and_limits["angle_bounds"]
        surfaces = np.zeros((len(distance_bounds), len(angle_bounds)), float)
        for idp in range(len(distance_bounds) - 1):
            this_dist_plateau = distance_bounds[idp + 1] - distance_bounds[idp]
            for iap in range(len(angle_bounds) - 1):
                this_ang_plateau = angle_bounds[iap + 1] - angle_bounds[iap]
                surfaces[idp][iap] = np.absolute(this_dist_plateau * this_ang_plateau)
        return surfaces

    def neighbors_surfaces_bounded(self, isite, surface_calculation_options=None):
        """Get the different surfaces (using boundaries) corresponding to the different distance-angle cutoffs
        for a given site.
        """
        if self.voronoi_list2[isite] is None:
            return None

        # Defaults
        if surface_calculation_options is None:
            surface_calculation_options = {
                "type": "standard_elliptic",
                "distance_bounds": {"lower": 1.2, "upper": 1.8},
                "angle_bounds": {"lower": 0.1, "upper": 0.8},
            }

        s_type = surface_calculation_options["type"]
        if s_type in ("standard_elliptic", "standard_diamond", "standard_spline"):
            plot_type = {
                "distance_parameter": ("initial_normalized", None),
                "angle_parameter": ("initial_normalized", None),
            }
        else:
            raise ValueError(f"Type {s_type!r} for the surface calculation in DetailedVoronoiContainer is invalid")

        # Hoist & cache lookups
        d_lower = surface_calculation_options["distance_bounds"]["lower"]
        d_upper = surface_calculation_options["distance_bounds"]["upper"]
        a_lower = surface_calculation_options["angle_bounds"]["lower"]
        a_upper = surface_calculation_options["angle_bounds"]["upper"]

        max_dist = d_upper + 0.1
        bounds_and_limits = self.voronoi_parameters_bounds_and_limits(
            isite=isite, plot_type=plot_type, max_dist=max_dist
        )

        distance_bounds = bounds_and_limits["distance_bounds"]
        angle_bounds = bounds_and_limits["angle_bounds"]

        # Functions and constant bounds for intersection
        lu = get_lower_and_upper_f(surface_calculation_options=surface_calculation_options)
        f_lower = lu["lower"]
        f_upper = lu["upper"]
        rect_bounds_lower = (d_lower, d_upper)
        rect_bounds_upper = (d_lower, d_upper)

        # Preallocate result
        n_d = len(distance_bounds)
        n_a = len(angle_bounds)
        surfaces = np.zeros((n_d, n_a), float)

        # Iterate over consecutive bound pairs without repeated indexing
        for idp, (dp1, dp2) in enumerate(itertools.pairwise(distance_bounds)):
            # Quick reject by distance window
            if dp2 < d_lower or dp1 > d_upper:
                continue
            # Clamp to window
            d1 = max(d_lower, dp1)
            d2 = min(d_upper, dp2)

            for iap, (ap1_raw, ap2_raw) in enumerate(itertools.pairwise(angle_bounds)):
                # Maintain original safety swap in case of non-monotonic inputs
                if ap1_raw <= ap2_raw:
                    ap1, ap2 = ap1_raw, ap2_raw
                else:
                    ap1, ap2 = ap2_raw, ap1_raw

                # Quick reject by angle window
                if ap2 < a_lower or ap1 > a_upper:
                    continue
                # Clamp to window
                a1 = max(a_lower, ap1)
                a2 = min(a_upper, ap2)

                # Compute intersection area
                inter, _ = rectangle_surface_intersection(
                    rectangle=((d1, d2), (a1, a2)),
                    f_lower=f_lower,
                    f_upper=f_upper,
                    bounds_lower=rect_bounds_lower,
                    bounds_upper=rect_bounds_upper,
                    check=False,
                )
                surfaces[idp, iap] = inter

        return surfaces

    @staticmethod
    def _get_vertices_dist_ang_indices(parameter_indices_list):
        pp0 = [pp[0] for pp in parameter_indices_list]
        pp1 = [pp[1] for pp in parameter_indices_list]
        min_idist = min(pp0)
        min_iang = min(pp1)
        max_idist = max(pp0)
        max_iang = max(pp1)
        i_min_angs = np.argwhere(np.array(pp1) == min_iang)
        i_max_dists = np.argwhere(np.array(pp0) == max_idist)
        pp0_at_min_iang = [pp0[ii[0]] for ii in i_min_angs]
        pp1_at_max_idist = [pp1[ii[0]] for ii in i_max_dists]
        max_idist_at_min_iang = max(pp0_at_min_iang)
        min_iang_at_max_idist = min(pp1_at_max_idist)

        p1 = (min_idist, min_iang)
        p2 = (max_idist_at_min_iang, min_iang)
        p3 = (max_idist_at_min_iang, min_iang_at_max_idist)
        p4 = (max_idist, min_iang_at_max_idist)
        p5 = (max_idist, max_iang)
        p6 = (min_idist, max_iang)

        return [p1, p2, p3, p4, p5, p6]

    def maps_and_surfaces(
        self,
        isite,
        surface_calculation_type=None,
        max_dist=2.0,
        additional_conditions=None,
    ):
        """Get the different surfaces and their cn_map corresponding to the different distance-angle cutoffs
        for a given site.

        Args:
            isite: Index of the site
            surface_calculation_type: How to compute the surface.
            max_dist: The maximum distance factor to be considered.
            additional_conditions: If additional conditions have to be considered.

        Returns:
            Surfaces and cn_map's for each distance-angle cutoff.
        """
        if self.voronoi_list2[isite] is None:
            return None
        if additional_conditions is None:
            additional_conditions = [self.AC.ONLY_ACB]
        surfaces = self.neighbors_surfaces(
            isite=isite,
            surface_calculation_type=surface_calculation_type,
            max_dist=max_dist,
        )
        maps_and_surfaces = []
        for cn, value in self._unique_coordinated_neighbors_parameters_indices[isite].items():
            for imap, list_parameters_indices in enumerate(value):
                this_surf = 0.0
                for idp, iap, iacb in list_parameters_indices:
                    if iacb in additional_conditions:
                        this_surf += surfaces[idp, iap]
                maps_and_surfaces.append(
                    {
                        "map": (cn, imap),
                        "surface": this_surf,
                        "parameters_indices": list_parameters_indices,
                    }
                )
        return maps_and_surfaces

    def maps_and_surfaces_bounded(self, isite, surface_calculation_options=None, additional_conditions=None):
        """Get the different surfaces (using boundaries) and their cn_map corresponding to the different
        distance-angle cutoffs for a given site.

        Args:
            isite: Index of the site
            surface_calculation_options: Options for the boundaries.
            additional_conditions: If additional conditions have to be considered.

        Returns:
            Surfaces and cn_map's for each distance-angle cutoff.
        """
        if self.voronoi_list2[isite] is None:
            return None
        if additional_conditions is None:
            additional_conditions = [self.AC.ONLY_ACB]
        surfaces = self.neighbors_surfaces_bounded(isite=isite, surface_calculation_options=surface_calculation_options)
        maps_and_surfaces = []
        for cn, value in self._unique_coordinated_neighbors_parameters_indices[isite].items():
            for imap, list_parameters_indices in enumerate(value):
                this_surf = 0.0
                for idp, iap, iacb in list_parameters_indices:
                    if iacb in additional_conditions:
                        this_surf += surfaces[idp, iap]
                maps_and_surfaces.append(
                    {
                        "map": (cn, imap),
                        "surface": this_surf,
                        "parameters_indices": list_parameters_indices,
                    }
                )
        return maps_and_surfaces

    def neighbors(self, isite, distfactor, angfactor, additional_condition=None):
        """Get the neighbors of a given site corresponding to a given distance and angle factor.

        Args:
            isite: Index of the site.
            distfactor: Distance factor.
            angfactor: Angle factor.
            additional_condition: Additional condition to be used (currently not implemented).

        Returns:
            List of neighbors of the given site for the given distance and angle factors.
        """
        idist = dfact = None
        for iwd, wd in enumerate(self.neighbors_normalized_distances[isite]):
            if distfactor >= wd["min"]:
                idist = iwd
                dfact = wd["max"]
            else:
                break
        iang = afact = None
        for iwa, wa in enumerate(self.neighbors_normalized_angles[isite]):
            if angfactor <= wa["max"]:
                iang = iwa
                afact = wa["min"]
            else:
                break
        if idist is None or iang is None:
            raise ValueError("Distance or angle parameter not found ...")

        return [
            nb
            for nb in self.voronoi_list2[isite]
            if nb["normalized_distance"] <= dfact and nb["normalized_angle"] >= afact
        ]

    def voronoi_parameters_bounds_and_limits(self, isite, plot_type, max_dist):
        """Get the different boundaries and limits of the distance and angle factors for the given site.

        Args:
            isite: Index of the site.
            plot_type: Types of distance/angle parameters to get.
            max_dist: Maximum distance factor.

        Returns:
            Distance and angle bounds and limits.
        """
        # Initializes the distance and angle parameters
        if self.voronoi_list2[isite] is None:
            return None
        if plot_type is None:
            plot_type = {
                "distance_parameter": ("initial_inverse_opposite", None),
                "angle_parameter": ("initial_opposite", None),
            }
        dd = [dist["min"] for dist in self.neighbors_normalized_distances[isite]]
        dd[0] = 1.0
        if plot_type["distance_parameter"][0] == "initial_normalized":
            dd.append(max_dist)
            distance_bounds = np.array(dd)
            dist_limits = [1.0, max_dist]
        elif plot_type["distance_parameter"][0] == "initial_inverse_opposite":
            ddinv = [1.0 / dist for dist in dd]
            ddinv.append(0.0)
            distance_bounds = np.array([1.0 - invdist for invdist in ddinv])
            dist_limits = [0.0, 1.0]
        elif plot_type["distance_parameter"][0] == "initial_inverse3_opposite":
            ddinv = [1.0 / dist**3.0 for dist in dd]
            ddinv.append(0.0)
            distance_bounds = np.array([1.0 - invdist for invdist in ddinv])
            dist_limits = [0.0, 1.0]
        else:
            raise NotImplementedError(
                f"Plotting type {plot_type['distance_parameter']!r} for the distance is not implemented"
            )
        if plot_type["angle_parameter"][0] == "initial_normalized":
            aa = [0.0]
            aa.extend([ang["max"] for ang in self.neighbors_normalized_angles[isite]])
            angle_bounds = np.array(aa)
        elif plot_type["angle_parameter"][0] == "initial_opposite":
            aa = [0.0]
            aa.extend([ang["max"] for ang in self.neighbors_normalized_angles[isite]])
            aa = [1.0 - ang for ang in aa]
            angle_bounds = np.array(aa)
        else:
            raise NotImplementedError(
                f"Plotting type {plot_type['angle_parameter']!r} for the angle is not implemented"
            )
        ang_limits = [0.0, 1.0]
        return {
            "distance_bounds": distance_bounds,
            "distance_limits": dist_limits,
            "angle_bounds": angle_bounds,
            "angle_limits": ang_limits,
        }

    def is_close_to(self, other, rtol=0.0, atol=1e-8) -> bool:
        if not (
            np.isclose(self.normalized_angle_tolerance, other.normalized_angle_tolerance, rtol=rtol, atol=atol)
            and np.isclose(
                self.normalized_distance_tolerance, other.normalized_distance_tolerance, rtol=rtol, atol=atol
            )
            and self.additional_conditions == other.additional_conditions
            and self.valences == other.valences
        ):
            return False

        # Cheap shape check
        if len(self.voronoi_list2) != len(other.voronoi_list2):
            return False

        for site_idx, voronoi_site in enumerate(self.voronoi_list2):
            other_site = other.voronoi_list2[site_idx]

            # Handle None symmetry
            if voronoi_site is None or other_site is None:
                if voronoi_site is not other_site:
                    return False
                continue

            # Quick length check
            if len(voronoi_site) != len(other_site):
                return False

            # Build an index -> neighbor dict for O(1) lookup on the "other" side.
            # Index should be unique per site; it is compared below anyway.
            other_by_index = {nb2["index"]: nb2 for nb2 in other_site}

            for nb in voronoi_site:
                if nb is None:
                    # Mirror None handling (shouldn't normally happen once list is non-None)
                    if any(x is not None for x in other_site):
                        return False
                    continue

                idx = nb["index"]
                nb_other = other_by_index.get(idx)
                if nb_other is None:
                    return False

                # Compare sites and numeric fields with tolerances
                if nb["site"] != nb_other["site"]:
                    return False

                if not np.isclose(nb["distance"], nb_other["distance"], rtol=rtol, atol=atol):
                    return False
                if not np.isclose(nb["angle"], nb_other["angle"], rtol=rtol, atol=atol):
                    return False
                if not np.isclose(nb["normalized_distance"], nb_other["normalized_distance"], rtol=rtol, atol=atol):
                    return False
                if not np.isclose(nb["normalized_angle"], nb_other["normalized_angle"], rtol=rtol, atol=atol):
                    return False

        return True

    def get_rdf_figure(self, isite, normalized=True, figsize=None, step_function=None):
        """Get the Radial Distribution Figure for a given site.

        Args:
            isite: Index of the site.
            normalized: Whether to normalize distances.
            figsize: Size of the figure.
            step_function: Type of step function to be used for the RDF.

        Returns:
            plt.figure: Matplotlib figure.
        """

        def dp_func(dp):
            return 1.0 - 1.0 / np.power(dp, 3.0)

        if step_function is None:
            step_function = {"type": "normal_cdf", "scale": 0.0001}

        # Initializes the figure
        fig = plt.figure() if figsize is None else plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        dists = self.neighbors_normalized_distances[isite] if normalized else self.neighbors_distances[isite]

        if step_function["type"] == "step_function":
            isorted = np.argsort([dd["min"] for dd in dists])
            sorted_dists = [dists[ii]["min"] for ii in isorted]
            dnb_dists = [len(dists[ii]["dnb_indices"]) for ii in isorted]
            xx = [0.0]
            yy = [0.0]
            for idist, dist in enumerate(sorted_dists):
                xx.extend((dist, dist))
                yy.extend((yy[-1], yy[-1] + dnb_dists[idist]))
            xx.append(1.1 * xx[-1])
            yy.append(yy[-1])
        elif step_function["type"] == "normal_cdf":
            scale = step_function["scale"]
            _dists = [dp_func(dd["min"]) for dd in dists]
            _dcns = [len(dd["dnb_indices"]) for dd in dists]
            xx = np.linspace(0.0, 1.1 * max(_dists), num=500)
            yy = np.zeros_like(xx)
            for idist, dist in enumerate(_dists):
                yy += _dcns[idist] * normal_cdf_step(xx, mean=dist, scale=scale)
        else:
            raise ValueError(f"Step function of type {step_function['type']!r} is not allowed")
        ax.plot(xx, yy)

        return fig

    def get_sadf_figure(self, isite, normalized=True, figsize=None, step_function=None):
        """Get the Solid Angle Distribution Figure for a given site.

        Args:
            isite: Index of the site.
            normalized: Whether to normalize angles.
            figsize: Size of the figure.
            step_function: Type of step function to be used for the SADF.

        Returns:
            plt.figure: matplotlib figure.
        """

        def ap_func(ap):
            return np.power(ap, -0.1)

        if step_function is None:
            step_function = {"type": "step_function", "scale": 0.0001}

        # Initializes the figure
        fig = plt.figure() if figsize is None else plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        angs = self.neighbors_normalized_angles[isite] if normalized else self.neighbors_angles[isite]

        if step_function["type"] == "step_function":
            isorted = np.argsort([ap_func(aa["min"]) for aa in angs])
            sorted_angs = [ap_func(angs[ii]["min"]) for ii in isorted]
            dnb_angs = [len(angs[ii]["dnb_indices"]) for ii in isorted]
            xx = [0.0]
            yy = [0.0]
            for iang, ang in enumerate(sorted_angs):
                xx.extend((ang, ang))
                yy.extend((yy[-1], yy[-1] + dnb_angs[iang]))
            xx.append(1.1 * xx[-1])
            yy.append(yy[-1])
        elif step_function["type"] == "normal_cdf":
            scale = step_function["scale"]
            _angles = [ap_func(aa["min"]) for aa in angs]
            _dcns = [len(dd["dnb_indices"]) for dd in angs]
            xx = np.linspace(0.0, 1.1 * max(_angles), num=500)
            yy = np.zeros_like(xx)
            for iang, ang in enumerate(_angles):
                yy += _dcns[iang] * normal_cdf_step(xx, mean=ang, scale=scale)
        else:
            raise ValueError(f"Step function of type {step_function['type']!r} is not allowed")
        ax.plot(xx, yy)

        return fig

    def __eq__(self, other: object) -> bool:
        needed_attrs = (
            "normalized_angle_tolerance",
            "normalized_distance_tolerance",
            "additional_conditions",
            "valences",
            "voronoi_list2",
            "structure",
        )
        if not all(hasattr(other, attr) for attr in needed_attrs):
            return NotImplemented
        return all(getattr(self, attr) == getattr(other, attr) for attr in needed_attrs)

    def to_bson_voronoi_list2(self):
        """
        Transforms the voronoi_list into a vlist + bson_nb_voro_list, that are BSON-encodable.

        Returns:
            [vlist, bson_nb_voro_list], to be used in the as_dict method.
        """
        bson_nb_voro_list2 = [None] * len(self.voronoi_list2)
        for ivoro, voro in enumerate(self.voronoi_list2):
            if voro is None or voro == "None":
                continue
            site_voro = []
            # {'site': neighbors[nn[1]],
            #  'angle': sa,
            #  'distance': distances[nn[1]],
            #  'index': myindex}
            for nb_dict in voro:
                site = nb_dict["site"]
                site_dict = {key: val for key, val in nb_dict.items() if key != "site"}
                # site_voro.append([ps.as_dict(), dd]) [float(c) for c in self.frac_coords]
                diff = site.frac_coords - self.structure[nb_dict["index"]].frac_coords
                site_voro.append([[nb_dict["index"], [float(c) for c in diff]], site_dict])
            bson_nb_voro_list2[ivoro] = site_voro
        return bson_nb_voro_list2

    def as_dict(self):
        """
        Bson-serializable dict representation of the VoronoiContainer.

        Returns:
            dictionary that is BSON-encodable.
        """
        bson_nb_voro_list2 = self.to_bson_voronoi_list2()
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "bson_nb_voro_list2": bson_nb_voro_list2,
            # "neighbors_lists": self.neighbors_lists,
            "structure": self.structure.as_dict(),
            "normalized_angle_tolerance": self.normalized_angle_tolerance,
            "normalized_distance_tolerance": self.normalized_distance_tolerance,
            "additional_conditions": self.additional_conditions,
            "valences": self.valences,
            "maximum_distance_factor": self.maximum_distance_factor,
            "minimum_angle_factor": self.minimum_angle_factor,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Reconstructs the VoronoiContainer object from a dict representation of the VoronoiContainer created using
        the as_dict method.

        Args:
            dct: dict representation of the VoronoiContainer object.

        Returns:
            VoronoiContainer object.
        """
        structure = Structure.from_dict(dct["structure"])
        voronoi_list2 = from_bson_voronoi_list2(dct["bson_nb_voro_list2"], structure)
        maximum_distance_factor = dct.get("maximum_distance_factor")
        minimum_angle_factor = dct.get("minimum_angle_factor")

        return cls(
            structure=structure,
            voronoi_list2=voronoi_list2,
            # neighbors_lists=neighbors_lists,
            normalized_angle_tolerance=dct["normalized_angle_tolerance"],
            normalized_distance_tolerance=dct["normalized_distance_tolerance"],
            additional_conditions=dct["additional_conditions"],
            valences=dct["valences"],
            maximum_distance_factor=maximum_distance_factor,
            minimum_angle_factor=minimum_angle_factor,
        )
