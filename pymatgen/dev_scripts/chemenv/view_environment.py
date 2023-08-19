"""Script to visualize the model coordination environments."""

from __future__ import annotations

import numpy as np

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import (
    SEPARATION_PLANE,
    AllCoordinationGeometries,
)
from pymatgen.analysis.chemenv.utils.coordination_geometry_utils import Plane
from pymatgen.analysis.chemenv.utils.scripts_utils import visualize

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"

if __name__ == "__main__":
    print(
        "+-------------------------------------------------------+\n"
        "| Development script of the ChemEnv utility of pymatgen |\n"
        "| Visualization of the model coordination environments  |\n"
        "+-------------------------------------------------------+\n"
    )
    allcg = AllCoordinationGeometries()
    vis = None
    while True:
        cg_symbol = input(
            'Enter symbol of the geometry you want to see, "l" to see the list '
            'of existing geometries or "q" to quit : '
        )
        if cg_symbol == "q":
            break
        if cg_symbol == "l":
            print(allcg.pretty_print(maxcn=20, additional_info={"nb_hints": True}))
            continue
        try:
            cg = allcg[cg_symbol]
        except LookupError:
            print("Wrong geometry, try again ...")
            continue
        print(cg.name)
        for idx, point in enumerate(cg.points):
            print(f"Point #{idx} : {point[0]!r} {point[1]!r} {point[2]!r}")
        print("Algorithms used :")
        for idx, algo in enumerate(cg.algorithms):
            print(f"Algorithm #{idx} :")
            print(algo)
            print()
        # Visualize the separation plane of a given algorithm
        sep_plane = False
        if any(algo.algorithm_type == SEPARATION_PLANE for algo in cg.algorithms):
            test = input("Enter index of the algorithm for which you want to visualize the plane : ")
            if test != "":
                try:
                    idx = int(test)
                    algo = cg.algorithms[idx]
                    sep_plane = True
                except Exception:
                    print(
                        "Unable to determine the algorithm/separation_plane you want "
                        "to visualize for this geometry. Continues without ..."
                    )
        my_factor = 3
        if vis is None:
            vis = visualize(cg=cg, zoom=1, myfactor=my_factor)
        else:
            vis = visualize(cg=cg, vis=vis, myfactor=my_factor)
        cg_points = [my_factor * np.array(pp) for pp in cg.points]
        cg_central_site = my_factor * np.array(cg.central_site)
        if sep_plane:
            pts = [cg_points[ii] for ii in algo.plane_points]
            if algo.minimum_number_of_points == 2:
                pts.append(cg_central_site)
                centre = cg_central_site
            else:
                centre = np.sum(pts, axis=0) / len(pts)

            factor = 1.5
            target_dist = max(np.dot(pp - centre, pp - centre) for pp in cg_points)
            current_dist = np.dot(pts[0] - centre, pts[0] - centre)
            factor = factor * target_dist / current_dist
            plane = Plane.from_npoints(points=pts)
            p1 = centre + factor * (pts[0] - centre)
            perp = factor * np.cross(pts[0] - centre, plane.normal_vector)
            p2 = centre + perp
            p3 = centre - factor * (pts[0] - centre)
            p4 = centre - perp

            vis.add_faces([[p1, p2, p3, p4]], [1, 0, 0], opacity=0.5)

            target_radius = 0.25
            radius = 1.5 * target_radius

            if algo.minimum_number_of_points == 2:
                vis.add_partial_sphere(
                    coords=cg_central_site, radius=radius, color=[1, 0, 0], start=0, end=360, opacity=0.5
                )
            for pp in pts:
                vis.add_partial_sphere(coords=pp, radius=radius, color=[1, 0, 0], start=0, end=360, opacity=0.5)

            ps1 = [cg_points[ii] for ii in algo.point_groups[0]]
            ps2 = [cg_points[ii] for ii in algo.point_groups[1]]

            for pp in ps1:
                vis.add_partial_sphere(coords=pp, radius=radius, color=[0, 1, 0], start=0, end=360, opacity=0.5)
            for pp in ps2:
                vis.add_partial_sphere(coords=pp, radius=radius, color=[0, 0, 1], start=0, end=360, opacity=0.5)
        vis.show()
