"""This module contains some script utils that are used in the chemenv package."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import (
    SimpleAbundanceChemenvStrategy,
    SimplestChemenvStrategy,
    TargetedPenaltiedAbundanceChemenvStrategy,
)
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import (
    UNCLEAR_ENVIRONMENT_SYMBOL,
    AllCoordinationGeometries,
)
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import (
    AbstractGeometry,
    LocalGeometryFinder,
)
from pymatgen.analysis.chemenv.utils.chemenv_errors import NeighborsNotComputedChemenvError
from pymatgen.analysis.chemenv.utils.coordination_geometry_utils import rotateCoords
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Molecule
from pymatgen.io.cif import CifParser

try:
    from pymatgen.vis.structure_vtk import StructureVis

except ImportError:
    StructureVis = None  # type: ignore[misc]

if TYPE_CHECKING:
    from typing import Any

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"

strategies_class_lookup: dict[str, Any] = {
    "SimplestChemenvStrategy": SimplestChemenvStrategy,
    "SimpleAbundanceChemenvStrategy": SimpleAbundanceChemenvStrategy,
    "TargetedPenaltiedAbundanceChemenvStrategy": TargetedPenaltiedAbundanceChemenvStrategy,
}


def draw_cg(
    vis,
    site,
    neighbors,
    cg=None,
    perm=None,
    perfect2local_map=None,
    show_perfect=False,
    csm_info=None,
    symmetry_measure_type="csm_wcs_ctwcc",
    perfect_radius=0.1,
    show_distorted=True,
    faces_color_override=None,
):
    """
    Draw cg.

    Args:
        site:
        vis:
        neighbors:
        cg:
        perm:
        perfect2local_map:
        show_perfect:
        csm_info:
        symmetry_measure_type:
        perfect_radius:
        show_distorted:
        faces_color_override:
    """
    csm_suffix = ""
    perf_radius = 0
    if show_perfect:
        if csm_info is None:
            raise ValueError("Not possible to show perfect environment without csm_info")
        csm_suffix = symmetry_measure_type[4:]
        perf_radius = (perfect_radius - 0.2) / 0.002
    if perm is not None and perfect2local_map is not None:
        raise ValueError('Only "perm" or "perfect2local_map" should be provided in draw_cg, not both')
    if show_distorted:
        vis.add_bonds(neighbors, site)
        for n in neighbors:
            vis.add_site(n)
    if len(neighbors) < 3:
        if show_distorted:
            vis.add_bonds(neighbors, site, color=[0.0, 1.0, 0.0], opacity=0.4, radius=0.175)
        if show_perfect and len(neighbors) == 2:
            perfect_geometry = AbstractGeometry.from_cg(cg)
            trans = csm_info["other_symmetry_measures"][f"translation_vector_{csm_suffix}"]
            rot = csm_info["other_symmetry_measures"][f"rotation_matrix_{csm_suffix}"]
            scale = csm_info["other_symmetry_measures"][f"scaling_factor_{csm_suffix}"]
            points = perfect_geometry.points_wcs_ctwcc()
            rotated_points = rotateCoords(points, rot)
            points = [scale * pp + trans for pp in rotated_points]
            ef_points = points[1:] if "wcs" in csm_suffix else points
            edges = cg.edges(ef_points, input="coords")
            vis.add_edges(edges, color=[1.0, 0.0, 0.0])
            for point in points:
                vis.add_partial_sphere(
                    coords=point,
                    radius=perf_radius,
                    color=[0.0, 0.0, 0.0],
                    start=0,
                    end=360,
                    opacity=1,
                )
    else:
        if show_distorted:
            if perm is not None:
                faces = cg.faces(neighbors, permutation=perm)
                edges = cg.edges(neighbors, permutation=perm)
            elif perfect2local_map is not None:
                faces = cg.faces(neighbors, perfect2local_map=perfect2local_map)
                edges = cg.edges(neighbors, perfect2local_map=perfect2local_map)
            else:
                faces = cg.faces(neighbors)
                edges = cg.edges(neighbors)
            symbol = next(iter(site.species)).symbol
            color = faces_color_override or [float(i) / 255 for i in vis.el_color_mapping[symbol]]
            vis.add_faces(faces, color, opacity=0.4)
            vis.add_edges(edges)
        if show_perfect:
            perfect_geometry = AbstractGeometry.from_cg(cg)
            trans = csm_info["other_symmetry_measures"][f"translation_vector_{csm_suffix}"]
            rot = csm_info["other_symmetry_measures"][f"rotation_matrix_{csm_suffix}"]
            scale = csm_info["other_symmetry_measures"][f"scaling_factor_{csm_suffix}"]
            points = perfect_geometry.points_wcs_ctwcc()
            rotated_points = rotateCoords(points, rot)
            points = [scale * pp + trans for pp in rotated_points]
            ef_points = points[1:] if "wcs" in csm_suffix else points
            edges = cg.edges(ef_points, input="coords")
            vis.add_edges(edges, color=[1.0, 0.0, 0.0])
            for point in points:
                vis.add_partial_sphere(
                    coords=point,
                    radius=perf_radius,
                    color=[0.0, 0.0, 0.0],
                    start=0,
                    end=360,
                    opacity=1,
                )


def visualize(cg, zoom=None, vis=None, factor=1.0, view_index=True, faces_color_override=None):
    """
    Visualizing a coordination geometry
    Args:
        cg:
        zoom:
        vis:
        factor:
        view_index:
        faces_color_override:
    """
    if vis is None and StructureVis is not None:
        vis = StructureVis(show_polyhedron=False, show_unit_cell=False)
    species = ["O"] * (cg.coordination_number + 1)
    species[0] = "Cu"
    coords = [np.zeros(3, float) + cg.central_site]

    for pp in cg.points:
        coords.append(np.array(pp) + cg.central_site)
    coords = [cc * factor for cc in coords]
    structure = Molecule(species=species, coords=coords)
    vis.set_structure(structure=structure, reset_camera=True)
    draw_cg(
        vis,
        site=structure[0],
        neighbors=structure[1:],
        cg=cg,
        faces_color_override=faces_color_override,
    )
    if view_index:
        for nbr_idx, neighbor in enumerate(structure[1:]):
            vis.add_text(neighbor.coords, f"{nbr_idx}", color=(0, 0, 0))
    if zoom is not None:
        vis.zoom(zoom)
    return vis


def compute_environments(chemenv_configuration):
    """
    Compute the environments.

    Args:
        chemenv_configuration:
    """
    string_sources = {
        "cif": {"string": "a Cif file", "regexp": r".*\.cif$"},
        "mp": {"string": "the Materials Project database", "regexp": r"mp-[0-9]+$"},
    }
    questions = {"c": "cif"}
    questions["m"] = "mp"
    lgf = LocalGeometryFinder()
    lgf.setup_parameters()
    all_cg = AllCoordinationGeometries()
    strategy_class = strategies_class_lookup[chemenv_configuration.package_options["default_strategy"]["strategy"]]
    # TODO: Add the possibility to change the parameters and save them in the chemenv_configuration
    default_strategy = strategy_class()
    default_strategy.setup_options(chemenv_configuration.package_options["default_strategy"]["strategy_options"])
    max_dist_factor = chemenv_configuration.package_options["default_max_distance_factor"]
    first_time = True
    test = None
    while True:
        if len(questions) > 1:
            found = False
            print("Enter the source from which the structure is coming or <q> to quit :")
            for key_character, qq in questions.items():
                print(f" - <{key_character}> for a structure from {string_sources[qq]['string']}")
            test = input(" ... ")
            source_type = ""
            if test == "q":
                break
            if test not in list(questions):
                for qq in questions.values():
                    if re.match(string_sources[qq]["regexp"], str(test)) is not None:
                        found = True
                        source_type = qq
                if not found:
                    print("Wrong key, try again ...")
                    continue
            else:
                source_type = questions[test]

        else:
            found = False
            source_type = next(iter(questions.values()))

        input_source = ""
        if found and len(questions) > 1:
            input_source = test

        structure = None
        if source_type == "cif":
            if not found:
                input_source = input("Enter path to CIF file : ")
            parser = CifParser(input_source)
            structure = parser.parse_structures(primitive=True)[0]

        elif source_type == "mp":
            if not found:
                input_source = input('Enter materials project id (e.g. "mp-1902") : ')
            from pymatgen.ext.matproj import MPRester

            with MPRester() as mpr:
                structure = mpr.get_structure_by_material_id(input_source)

        lgf.setup_structure(structure)
        print(f"Computing environments for {structure.reduced_formula} ... ")
        se = lgf.compute_structure_environments(maximum_distance_factor=max_dist_factor)
        print("Computing environments finished")
        while True:
            test = input(
                "See list of environments determined for each (inequivalent) site ? "
                '("y" or "n", "d" with details, "g" to see the grid) : '
            )
            strategy = default_strategy
            if test in {"y", "d", "g"}:
                strategy.set_structure_environments(se)
                for equiv_list in se.equivalent_sites:
                    site = equiv_list[0]
                    site_idx = se.structure.index(site)
                    try:
                        if strategy.uniquely_determines_coordination_environments:
                            ces = strategy.get_site_coordination_environments(site)
                        else:
                            ces = strategy.get_site_coordination_environments_fractions(site)
                    except NeighborsNotComputedChemenvError:
                        continue
                    if ces is None:
                        continue
                    if len(ces) == 0:
                        continue
                    comp = site.species
                    # ce = strategy.get_site_coordination_environment(site)
                    reduced_formula = comp.get_reduced_formula_and_factor()[0]
                    the_cg = None
                    if strategy.uniquely_determines_coordination_environments:
                        ce = ces[0]
                        if ce is None:
                            continue
                        the_cg = all_cg.get_geometry_from_mp_symbol(ce[0])
                        msg = f"Environment for site #{site_idx} {reduced_formula} ({comp}) : {the_cg.name} ({ce[0]})\n"
                    else:
                        msg = f"Environments for site #{site_idx} {reduced_formula} ({comp}) : \n"
                        for ce in ces:
                            cg = all_cg.get_geometry_from_mp_symbol(ce[0])
                            csm = ce[1]["other_symmetry_measures"]["csm_wcs_ctwcc"]
                            msg += f" - {cg.name} ({cg.mp_symbol}): {ce[2]:.2%} (csm : {csm:2f})\n"
                    if (
                        test in ["d", "g"]
                        and strategy.uniquely_determines_coordination_environments
                        and the_cg.mp_symbol != UNCLEAR_ENVIRONMENT_SYMBOL
                    ):
                        msg += "  <Continuous symmetry measures>  "
                        min_geoms = se.ce_list[site_idx][the_cg.coordination_number][0].minimum_geometries()
                        for min_geom in min_geoms:
                            csm = min_geom[1]["other_symmetry_measures"]["csm_wcs_ctwcc"]
                            msg += f"{min_geom[0]} : {csm:.2f}       "
                    print(msg)
            if test == "g":
                while True:
                    test = input(
                        "Enter index of site(s) (e.g. 0 1 2, separated by spaces) for which you want to see the grid "
                        "of parameters : "
                    )
                    try:
                        indices = [int(x) for x in test.split()]
                        print(str(indices))
                        for site_idx in indices:
                            if site_idx < 0:
                                raise IndexError
                            se.plot_environments(site_idx)
                        break
                    except ValueError:
                        print("This is not a valid site")
                    except IndexError:
                        print("This site is out of the site range")

            if StructureVis is None:
                test = input('Go to next structure ? ("y" to do so)')
                if test == "y":
                    break
                continue
            test = input('View structure with environments ? ("y" for the unit cell or "m" for a supercell or "n") : ')
            if test in ["y", "m"]:
                if test == "m":
                    deltas = []
                    while True:
                        try:
                            test = input("Enter multiplicity (e.g. 3 2 2) : ")
                            nns = test.split()
                            for i0 in range(int(nns[0])):
                                for i1 in range(int(nns[1])):
                                    for i2 in range(int(nns[2])):
                                        deltas.append(np.array([1.0 * i0, 1.0 * i1, 1.0 * i2], float))
                            break

                        except (ValueError, IndexError):
                            print("Not a valid multiplicity")
                else:
                    deltas = [np.zeros(3, float)]
                if first_time and StructureVis is not None:
                    vis = StructureVis(show_polyhedron=False, show_unit_cell=True)
                    vis.show_help = first_time = False
                else:
                    vis = None  # TODO: following code logic seems buggy

                vis.set_structure(se.structure)
                strategy.set_structure_environments(se)
                for site in se.structure:
                    try:
                        ces = strategy.get_site_coordination_environments(site)
                    except NeighborsNotComputedChemenvError:
                        continue
                    if len(ces) == 0:
                        continue
                    ce = strategy.get_site_coordination_environment(site)
                    if ce is not None and ce[0] != UNCLEAR_ENVIRONMENT_SYMBOL:
                        for delta in deltas:
                            p_site = PeriodicSite(
                                site.species,
                                site.frac_coords + delta,
                                site.lattice,
                                properties=site.properties,
                            )
                            vis.add_site(p_site)
                            neighbors = strategy.get_site_neighbors(p_site)
                            draw_cg(
                                vis,
                                p_site,
                                neighbors,
                                cg=lgf.allcg.get_geometry_from_mp_symbol(ce[0]),
                                perm=ce[1]["permutation"],
                            )
                vis.show()
            test = input('Go to next structure ? ("y" to do so) : ')
            if test == "y":
                break
        print()
