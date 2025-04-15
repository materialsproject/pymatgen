"""
This module define a WulffShape class to generate the Wulff shape from
a lattice, a list of indices and their corresponding surface energies,
and the total area and volume of the Wulff shape, the weighted surface energy,
the anisotropy and shape_factor can also be calculated.
In support of plotting from a given view in terms of miller index.

The lattice is from the conventional unit cell, and (hkil) for hexagonal
lattices.

If you use this code extensively, consider citing the following:

Tran, R.; Xu, Z.; Radhakrishnan, B.; Winston, D.; Persson, K. A.; Ong, S. P.
(2016). Surface energies of elemental crystals. Scientific Data.
"""

from __future__ import annotations

import itertools
import logging
import warnings
from typing import TYPE_CHECKING

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from scipy.spatial import ConvexHull

from pymatgen.core.structure import Structure
from pymatgen.util.coord import get_angle
from pymatgen.util.string import unicodeify_spacegroup

if TYPE_CHECKING:
    from pymatgen.core.lattice import Lattice

__author__ = "Zihan Xu, Richard Tran, Shyue Ping Ong"
__copyright__ = "Copyright 2013, The Materials Virtual Lab"
__version__ = "0.1"
__maintainer__ = "Zihan Xu"
__email__ = "zix009@eng.ucsd.edu"
__date__ = "May 5 2016"

logger = logging.getLogger(__name__)


def hkl_tuple_to_str(hkl):
    """
    Prepare for display on plots "(hkl)" for surfaces.

    Args:
        hkl: in the form of [h, k, l] or (h, k, l).
    """
    out = "".join(f"\\overline{{{-x}}}" if x < 0 else str(x) for x in hkl)
    return f"(${out}$)"


def get_tri_area(pts):
    """
    Given a list of coords for 3 points,
    Compute the area of this triangle.

    Args:
        pts: [a, b, c] three points
    """
    a, b, c = pts
    v1 = np.array(b) - np.array(a)
    v2 = np.array(c) - np.array(a)
    return abs(np.linalg.norm(np.cross(v1, v2)) / 2)


class WulffFacet:
    """Helper container for each Wulff plane."""

    def __init__(self, normal, e_surf, normal_pt, dual_pt, index, m_ind_orig, miller):
        """
        Args:
            normal:
            e_surf:
            normal_pt:
            dual_pt:
            index:
            m_ind_orig:
            miller:
        """
        self.normal = normal
        self.e_surf = e_surf
        self.normal_pt = normal_pt
        self.dual_pt = dual_pt
        self.index = index
        self.m_ind_orig = m_ind_orig
        self.miller = miller
        self.points: list = []
        self.outer_lines: list = []


class WulffShape:
    """
    Generate Wulff Shape from list of miller index and surface energies,
    with given conventional unit cell.
    surface energy (Jm^2) is the length of normal.

    Wulff shape is the convex hull.
    Based on:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html

    Process:
        1. get Wulff simplices
        2. label with color
        3. get wulff_area and other properties

    Attributes:
        debug (bool): Whether to print debug information.
        alpha (float): Transparency of the Wulff shape.
        color_set (list): colors to use for facets.
        grid_off (bool): Whether to turn off the grid.
        axis_off (bool): Whether to turn off the axis.
        show_area (bool): Whether to show the area of each facet.
        off_color (str): Color of facets not on the Wulff shape.
        structure (Structure): Input conventional unit cell (with H) from lattice.
        miller_list (list): input Miller indices, for hcp in the form of hkil.
        hkl_list (list): Modified Miller indices in the same order as input_miller.
        e_surf_list (list): input surface energies in the same order as input_miller.
        lattice (Lattice): Input lattice for the conventional unit cell.
        facets (list): WulffFacet objects considering symmetry.
        dual_cv_simp (list): Simplices from the dual convex hull (dual_pt).
        wulff_pt_list (list): Wulff points.
        wulff_cv_simp (list): Simplices from the convex hull of wulff_pt_list.
        on_wulff (list): List for all input_miller, True if on the Wulff shape.
        color_area (list): List for all input_miller, total area on the Wulff shape, off_wulff = 0.
        miller_area (dict): Dictionary of Miller indices and their corresponding areas.
    """

    def __init__(self, lattice: Lattice, miller_list, e_surf_list, symprec=1e-5):
        """
        Args:
            lattice: Lattice object of the conventional unit cell
            miller_list ([(hkl), ...]: list of hkl or hkil for hcp
            e_surf_list ([float]): list of corresponding surface energies
            symprec (float): for reciprocal lattice operation, default is 1e-5.
        """
        if any(se < 0 for se in e_surf_list):
            warnings.warn("Unphysical (negative) surface energy detected.", stacklevel=2)

        self.color_ind = list(range(len(miller_list)))

        self.input_miller_fig = [hkl_tuple_to_str(x) for x in miller_list]
        # store input data
        self.structure = Structure(lattice, ["H"], [[0, 0, 0]])
        self.miller_list = tuple(tuple(x) for x in miller_list)
        self.hkl_list = tuple((x[0], x[1], x[-1]) for x in miller_list)
        self.e_surf_list = tuple(e_surf_list)
        self.lattice = lattice
        self.symprec = symprec

        # 2. get all the data for wulff construction
        # get all the surface normal from get_all_miller_e()
        self.facets = self._get_all_miller_e()
        logger.debug(len(self.facets))

        # 3. consider the dual condition
        dual_pts = [facet.dual_pt for facet in self.facets]
        dual_convex = ConvexHull(dual_pts)
        dual_cv_simp = dual_convex.simplices
        # simplices	(ndarray of ints, shape (n_facet, n_dim))
        # list of [i, j, k] , n_dim = 3
        # i, j, k: ind for normal_e_m
        # recalculate the dual of dual, get the Wulff shape.
        # corner <-> surface
        # get cross point from the simplices of the dual convex hull
        wulff_pt_list = [self._get_cross_pt_dual_simp(dual_simp) for dual_simp in dual_cv_simp]

        wulff_convex = ConvexHull(wulff_pt_list)
        wulff_cv_simp = wulff_convex.simplices
        logger.debug(", ".join(str(len(x)) for x in wulff_cv_simp))

        # store simplices and convex
        self.dual_cv_simp = dual_cv_simp
        self.wulff_pt_list = wulff_pt_list
        self.wulff_cv_simp = wulff_cv_simp
        self.wulff_convex = wulff_convex

        self.on_wulff, self.color_area = self._get_simpx_plane()

        miller_area = []
        for miller, in_mill_fig in enumerate(self.input_miller_fig):
            miller_area.append(f"{in_mill_fig} : {round(self.color_area[miller], 4)}")
        self.miller_area = miller_area

    def _get_all_miller_e(self):
        """From self: get miller_list(unique_miller), e_surf_list and symmetry operations(symm_ops)
        according to lattice apply symm_ops to get all the miller index, then get normal, get
        all the facets functions for Wulff shape calculation: |normal| = 1, e_surf is plane's
        distance to (0, 0, 0), normal[0]x + normal[1]y + normal[2]z = e_surf.

        Returns:
            [WulffFacet]
        """
        all_hkl = []
        color_ind = self.color_ind
        planes = []
        recp = self.structure.lattice.reciprocal_lattice_crystallographic
        recp_symm_ops = self.lattice.get_recp_symmetry_operation(self.symprec)

        for idx, (hkl, energy) in enumerate(zip(self.hkl_list, self.e_surf_list, strict=True)):
            for op in recp_symm_ops:
                miller = tuple(int(x) for x in op.operate(hkl))
                if miller not in all_hkl:
                    all_hkl.append(miller)
                    normal = recp.get_cartesian_coords(miller)
                    normal /= np.linalg.norm(normal)
                    normal_pt = [x * energy for x in normal]
                    dual_pt = [x / energy for x in normal]
                    color_plane = color_ind[divmod(idx, len(color_ind))[1]]
                    planes.append(WulffFacet(normal, energy, normal_pt, dual_pt, color_plane, idx, hkl))

        # sort by e_surf
        planes.sort(key=lambda x: x.e_surf)
        return planes

    def _get_cross_pt_dual_simp(self, dual_simp):
        """
        |normal| = 1, e_surf is plane's distance to (0, 0, 0),
        plane function:
            normal[0]x + normal[1]y + normal[2]z = e_surf.

        from self:
            normal_e_m to get the plane functions
            dual_simp: (i, j, k) simplices from the dual convex hull
                i, j, k: plane index(same order in normal_e_m)
        """
        matrix_surfs = [self.facets[dual_simp[i]].normal for i in range(3)]
        matrix_e = [self.facets[dual_simp[i]].e_surf for i in range(3)]
        return np.dot(np.linalg.inv(matrix_surfs), matrix_e)

    def _get_simpx_plane(self):
        """
        Locate the plane for simpx of on wulff_cv, by comparing the center of
        the simpx triangle with the plane functions.
        """
        on_wulff = [False] * len(self.miller_list)
        surface_area = [0.0] * len(self.miller_list)
        for simpx in self.wulff_cv_simp:
            pts = [self.wulff_pt_list[simpx[i]] for i in range(3)]
            center = np.sum(pts, 0) / 3.0
            # check whether the center of the simplices is on one plane
            for plane in self.facets:
                abs_diff = abs(np.dot(plane.normal, center) - plane.e_surf)
                if abs_diff < 1e-5:
                    on_wulff[plane.index] = True
                    surface_area[plane.index] += get_tri_area(pts)

                    plane.points.append(pts)
                    plane.outer_lines.append([simpx[0], simpx[1]])
                    plane.outer_lines.append([simpx[1], simpx[2]])
                    plane.outer_lines.append([simpx[0], simpx[2]])
                    # already find the plane, move to the next simplices
                    break
        for plane in self.facets:
            plane.outer_lines.sort()
            plane.outer_lines = [line for line in plane.outer_lines if plane.outer_lines.count(line) != 2]
        return on_wulff, surface_area

    def _get_colors(self, color_set, alpha, off_color, custom_colors=None):
        """
        Assign colors according to the surface energies of on_wulff facets.

        Returns:
            tuple: color_list, color_proxy, color_proxy_on_wulff, miller_on_wulff, e_surf_on_wulff_list
        """
        color_list = [off_color] * len(self.hkl_list)
        color_proxy_on_wulff = []
        miller_on_wulff = []
        e_surf_on_wulff = [(idx, e_surf) for idx, e_surf in enumerate(self.e_surf_list) if self.on_wulff[idx]]

        c_map = plt.get_cmap(color_set)
        e_surf_on_wulff.sort(key=lambda x: x[1], reverse=False)
        e_surf_on_wulff_list = [x[1] for x in e_surf_on_wulff]
        if len(e_surf_on_wulff) > 1:
            color_norm = mpl.colors.Normalize(vmin=min(e_surf_on_wulff_list), vmax=max(e_surf_on_wulff_list))
        else:
            # if there is only one hkl on wulff, choose the color of the median
            color_norm = mpl.colors.Normalize(
                vmin=min(e_surf_on_wulff_list) - 0.1,
                vmax=max(e_surf_on_wulff_list) + 0.1,
            )
        scalar_map = mpl.cm.ScalarMappable(norm=color_norm, cmap=c_map)

        for idx, e_surf in e_surf_on_wulff:
            color_list[idx] = scalar_map.to_rgba(e_surf, alpha=alpha)
            if tuple(self.miller_list[idx]) in custom_colors:
                color_list[idx] = custom_colors[tuple(self.miller_list[idx])]
            color_proxy_on_wulff.append(plt.Rectangle((2, 2), 1, 1, fc=color_list[idx], alpha=alpha))
            miller_on_wulff.append(self.input_miller_fig[idx])
        scalar_map.set_array([x[1] for x in e_surf_on_wulff])
        color_proxy = [plt.Rectangle((2, 2), 1, 1, fc=x, alpha=alpha) for x in color_list]

        return (
            color_list,
            color_proxy,
            color_proxy_on_wulff,
            miller_on_wulff,
            e_surf_on_wulff_list,
        )

    def show(self, *args, **kwargs):
        """Show the Wulff plot.

        Args:
            *args: Passed to get_plot.
            **kwargs: Passed to get_plot.
        """
        self.get_plot(*args, **kwargs).get_figure().show()

    def get_line_in_facet(self, facet):
        """Get the sorted pts in a facet used to draw a line."""
        lines = list(facet.outer_lines)
        pt = []
        prev = line = None
        while len(lines) > 0:
            if prev is None:
                line = lines.pop(0)
            else:
                for idx, line in enumerate(lines):
                    if prev in line:
                        line = lines.pop(idx)
                        if line[1] == prev:
                            line.reverse()
                        break
            # make sure the lines are connected one by one.
            # find the way covering all pts and facets
            pt.extend(
                (
                    self.wulff_pt_list[line[0]].tolist(),
                    self.wulff_pt_list[line[1]].tolist(),
                )
            )
            prev = line[1]

        return pt

    def get_plot(
        self,
        color_set="PuBu",
        grid_off=True,
        axis_off=True,
        show_area=False,
        alpha=1,
        off_color="red",
        direction=None,
        bar_pos=(0.75, 0.15, 0.05, 0.65),
        bar_on=False,
        units_in_JPERM2=True,
        legend_on=True,
        aspect_ratio=(8, 8),
        custom_colors=None,
    ):
        """Get the Wulff shape plot.

        Args:
            color_set: default is 'PuBu'
            grid_off (bool): default is True
            axis_off (bool): default is True
            show_area (bool): default is False
            alpha (float): chosen from 0 to 1 (float), default is 1
            off_color: Default color for facets not present on the Wulff shape.
            direction: default is (1, 1, 1)
            bar_pos: default is [0.75, 0.15, 0.05, 0.65]
            bar_on (bool): default is False
            legend_on (bool): default is True
            aspect_ratio: default is (8, 8)
            custom_colors ({(h,k,l}: [r,g,b,alpha]}): Customize color of each
                facet with a dictionary. The key is the corresponding Miller
                index and value is the color. Undefined facets will use default
                color site. Note: If you decide to set your own colors, it
                probably won't make any sense to have the color bar on.
            units_in_JPERM2 (bool): Units of surface energy, defaults to
                Joules per square meter (True)

        Returns:
            mpl_toolkits.mplot3d.Axes3D: 3D plot of the Wulff shape.
        """
        from mpl_toolkits.mplot3d import art3d

        colors = self._get_colors(color_set, alpha, off_color, custom_colors=custom_colors or {})
        (
            color_list,
            color_proxy,
            color_proxy_on_wulff,
            miller_on_wulff,
            e_surf_on_wulff,
        ) = colors

        if not direction:
            # If direction is not specified, use the miller indices of
            # maximum area.
            direction = max(self.area_fraction_dict.items(), key=lambda x: x[1])[0]

        fig = plt.figure()
        fig.set_size_inches(aspect_ratio[0], aspect_ratio[1])
        azim, elev = self._get_azimuth_elev([direction[0], direction[1], direction[-1]])

        wulff_pt_list = self.wulff_pt_list

        ax_3d = fig.add_subplot(projection="3d")
        ax_3d.view_init(azim=azim, elev=elev)
        fig.add_axes(ax_3d)

        for plane in self.facets:
            # check whether [pts] is empty
            if len(plane.points) < 1:
                # empty, plane is not on_wulff.
                continue
            # assign the color for on_wulff facets according to its
            # index and the color_list for on_wulff
            plane_color = color_list[plane.index]
            pt = self.get_line_in_facet(plane)
            # plot from the sorted pts from [simpx]
            tri = art3d.Poly3DCollection([pt])
            tri.set_color(plane_color)
            tri.set_edgecolor("#808080")
            ax_3d.add_collection3d(tri)

        # set ranges of x, y, z
        # find the largest distance between on_wulff pts and the origin,
        # to ensure complete and consistent display for all directions
        r_range = max(np.linalg.norm(x) for x in wulff_pt_list)
        ax_3d.set_xlim([-r_range * 1.1, r_range * 1.1])
        ax_3d.set_ylim([-r_range * 1.1, r_range * 1.1])
        ax_3d.set_zlim([-r_range * 1.1, r_range * 1.1])
        # add legend
        if legend_on:
            if show_area:
                ax_3d.legend(
                    color_proxy,
                    self.miller_area,
                    loc="upper left",
                    bbox_to_anchor=(0, 1),
                    fancybox=True,
                    shadow=False,
                )
            else:
                ax_3d.legend(
                    color_proxy_on_wulff,
                    miller_on_wulff,
                    loc="upper center",
                    bbox_to_anchor=(0.5, 1),
                    ncol=3,
                    fancybox=True,
                    shadow=False,
                )
        ax_3d.set(xlabel="x", ylabel="y", zlabel="z")

        # Add color bar
        if bar_on:
            cmap = plt.get_cmap(color_set)
            cmap.set_over("0.25")
            cmap.set_under("0.75")
            bounds = [round(ene, 2) for ene in e_surf_on_wulff]
            bounds.append(1.2 * bounds[-1])
            norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
            # display surface energies
            ax1 = fig.add_axes(bar_pos)
            cbar = mpl.colorbar.ColorbarBase(
                ax1,
                cmap=cmap,
                norm=norm,
                boundaries=[0, *bounds, 10],
                extend="both",
                ticks=bounds[:-1],
                spacing="proportional",
                orientation="vertical",
            )
            units = "$J/m^2$" if units_in_JPERM2 else r"$eV/\AA^2$"
            cbar.set_label(f"Surface Energies ({units})", fontsize=25)

        if grid_off:
            ax_3d.grid("off")
        if axis_off:
            ax_3d.axis("off")
        return ax_3d

    def get_plotly(
        self,
        color_set="PuBu",
        off_color="red",
        alpha=1,
        custom_colors=None,
        units_in_JPERM2=True,
    ):
        """Get the Wulff shape as a plotly Figure object.

        Args:
            color_set: default is 'PuBu'
            alpha (float): chosen from 0 to 1 (float), default is 1
            off_color: Default color for facets not present on the Wulff shape.
            custom_colors ({(h,k,l}: [r,g,b,alpha}): Customize color of each
                facet with a dictionary. The key is the corresponding Miller
                index and value is the color. Undefined facets will use default
                color site. Note: If you decide to set your own colors, it
                probably won't make any sense to have the color bar on.
            units_in_JPERM2 (bool): Units of surface energy, defaults to
                Joules per square meter (True)

        Returns:
            (plotly.graph_objects.Figure)
        """
        units = "Jm⁻²" if units_in_JPERM2 else "eVÅ⁻²"
        (
            color_list,
            _color_proxy,
            _color_proxy_on_wulff,
            _miller_on_wulff,
            e_surf_on_wulff,
        ) = self._get_colors(color_set, alpha, off_color, custom_colors=custom_colors or {})

        planes_data, color_scale, ticktext, tickvals = [], [], [], []
        for plane in self.facets:
            if len(plane.points) < 1:
                # empty, plane is not on_wulff.
                continue

            plane_color = color_list[plane.index]
            plane_color = (1, 0, 0, 1) if plane_color == off_color else plane_color  # set to red for now

            pt = self.get_line_in_facet(plane)
            x_pts, y_pts, z_pts = [], [], []
            for p in pt:
                x_pts.append(p[0])
                y_pts.append(p[1])
                z_pts.append(p[2])

            # remove duplicate x y z pts to save time
            all_xyz = []

            [all_xyz.append(list(coord)) for coord in np.array([x_pts, y_pts, z_pts]).T if list(coord) not in all_xyz]
            all_xyz = np.array(all_xyz).T
            x_pts, y_pts, z_pts = all_xyz[0], all_xyz[1], all_xyz[2]
            index_list = [int(i) for i in np.linspace(0, len(x_pts) - 1, len(x_pts))]

            tri_indices = np.array(list(itertools.combinations(index_list, 3))).T
            hkl = self.miller_list[plane.index]
            hkl = unicodeify_spacegroup(f"({'%s' * len(hkl) % hkl})")
            cs = tuple(np.array(plane_color) * 255)
            color = f"rgba({cs[0]:.5f}, {cs[1]:.5f}, {cs[2]:.5f}, {cs[3]:.5f})"

            # note hoverinfo is incompatible with latex, need unicode instead
            planes_data.append(
                go.Mesh3d(
                    x=x_pts,
                    y=y_pts,
                    z=z_pts,
                    i=tri_indices[0],
                    j=tri_indices[1],
                    k=tri_indices[2],
                    hovertemplate=f"<br>%{{text}}<br>y={plane.e_surf:.3f} {units}<br>",
                    color=color,
                    text=[f"Miller index: {hkl}"] * len(x_pts),
                    hoverinfo="name",
                    name="",
                )
            )

            # normalize surface energy from a scale of 0 to 1 for colorbar
            if max(e_surf_on_wulff) == min(e_surf_on_wulff):
                norm_e = 1
            else:
                norm_e = (plane.e_surf - min(e_surf_on_wulff)) / (max(e_surf_on_wulff) - min(e_surf_on_wulff))
            c = [norm_e, color]
            if c not in color_scale:
                color_scale.append(c)
                ticktext.append(f"{plane.e_surf:.3f}")
                tickvals.append(norm_e)

        # Add colorbar
        color_scale = sorted(color_scale, key=lambda c: c[0])
        colorbar = go.Mesh3d(
            x=[0],
            y=[0],
            z=[0],
            colorbar=go.mesh3d.ColorBar(
                title={
                    "text": f"Surface energy {units}",
                    "side": "right",
                    "font": {"size": 25},
                },
                ticktext=ticktext,
                tickvals=tickvals,
            ),
            colorscale=[[0, "rgb(255,255,255, 255)"], *color_scale],  # fix the scale
            intensity=[0, 0.33, 0.66, 1],
            i=[0],
            j=[0],
            k=[0],
            name="y",
            showscale=True,
        )
        planes_data.append(colorbar)

        # Format aesthetics: background, axis, etc.
        axis_dict = {
            "title": "",
            "autorange": True,
            "showgrid": False,
            "zeroline": False,
            "ticks": "",
            "showline": False,
            "showticklabels": False,
            "showbackground": False,
        }
        fig = go.Figure(data=planes_data)
        fig.layout.update(
            showlegend=True,
            scene={"xaxis": axis_dict, "yaxis": axis_dict, "zaxis": axis_dict},
        )

        return fig

    def _get_azimuth_elev(self, miller_index):
        """
        Args:
            miller_index: viewing direction.

        Returns:
            azim, elev for plotting
        """
        if miller_index in [(0, 0, 1), (0, 0, 0, 1)]:
            return 0, 90
        cart = self.lattice.get_cartesian_coords(miller_index)
        azim = get_angle([cart[0], cart[1], 0], (1, 0, 0))
        v = [cart[0], cart[1], 0]
        elev = get_angle(cart, v)
        return azim, elev

    @property
    def volume(self) -> float:
        """Volume of the Wulff shape."""
        return self.wulff_convex.volume

    @property
    def miller_area_dict(self) -> dict[tuple, float]:
        """{hkl: area_hkl on wulff}."""
        return dict(zip(self.miller_list, self.color_area, strict=True))

    @property
    def miller_energy_dict(self) -> dict[tuple, float]:
        """{hkl: surface energy_hkl}."""
        return dict(zip(self.miller_list, self.e_surf_list, strict=True))

    @property
    def surface_area(self) -> float:
        """Total surface area of Wulff shape."""
        return sum(self.miller_area_dict.values())

    @property
    def weighted_surface_energy(self) -> float:
        """
        Returns:
            sum(surface_energy_hkl * area_hkl)/ sum(area_hkl).
        """
        return self.total_surface_energy / self.surface_area

    @property
    def area_fraction_dict(self) -> dict[tuple, float]:
        """
        Returns:
            dict: {hkl: area_hkl/total area on wulff}.
        """
        return {hkl: area / self.surface_area for hkl, area in self.miller_area_dict.items()}

    @property
    def anisotropy(self) -> float:
        """
        Returns:
            float: Coefficient of Variation from weighted surface energy. The ideal sphere is 0.
        """
        square_diff_energy = 0.0
        weighted_energy = self.weighted_surface_energy
        area_frac_dict = self.area_fraction_dict
        miller_energy_dict = self.miller_energy_dict

        for hkl, energy in miller_energy_dict.items():
            square_diff_energy += (energy - weighted_energy) ** 2 * area_frac_dict[hkl]
        return np.sqrt(square_diff_energy) / weighted_energy

    @property
    def shape_factor(self) -> float:
        """Determine the critical nucleus size.
        A large shape factor indicates great anisotropy.
        See Ballufi, R. W., Allen, S. M. & Carter, W. C. Kinetics
            of Materials. (John Wiley & Sons, 2005), p.461.

        Returns:
            float: Shape factor.
        """
        return self.surface_area / (self.volume ** (2 / 3))

    @property
    def effective_radius(self) -> float:
        """
        Radius of the WulffShape (in Angstroms) when the WulffShape is approximated as a sphere.

        Returns:
            float: radius R_eff
        """
        return ((3 / 4) * (self.volume / np.pi)) ** (1 / 3)

    @property
    def total_surface_energy(self) -> float:
        """
        Total surface energy of the Wulff shape.

        Returns:
            float: sum(surface_energy_hkl * area_hkl)
        """
        tot_surface_energy = 0.0
        for hkl, energy in self.miller_energy_dict.items():
            tot_surface_energy += energy * self.miller_area_dict[hkl]
        return tot_surface_energy

    @property
    def tot_corner_sites(self):
        """The number of vertices in the convex hull.
        Useful for identifying catalytically active sites.
        """
        return len(self.wulff_convex.vertices)

    @property
    def tot_edges(self):
        """The number of edges in the convex hull.
        Useful for identifying catalytically active sites.
        """
        all_edges = []
        for facet in self.facets:
            edges = []
            pt = self.get_line_in_facet(facet)

            lines = []
            for idx in range(len(pt)):
                if idx == len(pt) / 2:
                    break
                lines.append(tuple(sorted((tuple(pt[idx * 2]), tuple(pt[idx * 2 + 1])))))

            for p in lines:
                if p not in all_edges:
                    edges.append(p)

            all_edges.extend(edges)

        return len(all_edges)
