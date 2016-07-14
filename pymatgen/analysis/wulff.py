# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module define a WulffShape class to generate the Wulff shape from
a lattice, a list of indices and their corresponding surface energies,
and the total area and volume of the wulff shape,the weighted surface energy,
the anisotropy and shape_factor can also be calculated.
In support of plotting from a given view in terms of miller index.

The lattice is from the conventional unit cell, and (hkil) for hexagonal lattices.
"""

from __future__ import division, unicode_literals
from pymatgen.core.structure import Structure
from pymatgen.core.surface import get_recp_symmetry_operation
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord_utils import get_angle
import numpy as np
import scipy as scp
from scipy.spatial import ConvexHull
import copy
import logging

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
import matplotlib.cm as cm
import mpl_toolkits.mplot3d as a3

__author__ = 'Zihan Xu, Richard Tran'
__copyright__ = 'Copyright 2013, The Materials Virtual Lab'
__version__ = '0.1'
__maintainer__ = 'Zihan Xu'
__email__ = 'zix009@eng.ucsd.edu'
__date__ = 'May 5 2016'

logger = logging.getLogger(__name__)


def hkl_tuple_to_str(hkl):
    """
    Prepare for display on plots
    "(hkl)" for surfaces
    Agrs:
        hkl: in the form of [h, k, l] or (h, k, l)
    """
    str_format = '($'
    for x in hkl:
        if x < 0:
            str_format += '\overline{' + str(-x) + '}'
        else:
            str_format += str(x)
    str_format += '$)'
    return str_format


def get_tri_area(pts):
    """
    Given a list of coords for 3 points,
    Compute the area of this triangle.
    Args:
        pts: [a, b, c] three points
    """
    a, b, c = pts[0], pts[1], pts[2]
    v1 = np.array(b) - np.array(a)
    v2 = np.array(c) - np.array(a)
    area_tri = abs(scp.linalg.norm(scp.cross(v1, v2)) / 2)
    return area_tri


class WulffShape(object):
    """
    Generate Wulff Shape from list of miller index and surface energies,
    with given conventional unit cell.
    surface energy (Jm^2) is the length of normal.

    Wulff shape is the convex hull.
    Based on:
    http://scipy.github.io/devdocs/generated/scipy.spatial.ConvexHull.html

    Process:
        1. get wulff simplices
        2. label with color
        3. get wulff_area and other properties

    .. attribute:: debug (bool)

    .. attribute:: alpha
        transparency

    .. attribute:: color_set

    .. attribute:: grid_off (bool)

    .. attribute:: axis_off (bool)

    .. attribute:: show_area

    .. attribute:: off_color
        color of planes off wulff

    .. attribute:: input_miller_fig
        ($input_miller$) for legend on figures

    .. attribute:: structure
        Structure object, input conventional unit cell (with H ) from lattice

    .. attribute:: input_miller
        list of input miller index, for hcp in the form of hkil

    .. attribute:: input_hkl
        modify hkill to hkl, in the same order with input_miller

    .. attribute:: input_e_surf
        list of input surface energies, in the same order with input_miller

    .. attribute:: latt
        Lattice object, the input lattice for the conventional unit cell

    .. attribute:: recp
        reciprocal_lattice_crystallographic of the input structure

    .. attribute:: recp_symmops
        list of symmetric operations for recp

    .. attribute:: cart_symmops
        list of symmetric operations for input structure (cartesian)

    .. attribute:: normal_e_m: item:
        [normal, e_surf, normal_pt, dual_pt, color_plane, m_ind_orig, miller]
        for all planes considering symm

    .. attribute:: dual_cv_simp
        simplices from the dual convex hull (dual_pt)

    .. attribute:: wulff_pt_list

    .. attribute:: wulff_cv_simp
        simplices from the convex hull of wulff_pt_list

    .. attribute:: simpx_info

    .. attribute:: plane_wulff_info: item:
        [normal, e_surf, pts, simpx, color_plane, m_ind_orig, miller]

    .. attribute:: on_wulff
        list for all input_miller, True is on wulff.

    .. attribute:: color_area
        list for all input_miller, total area on wulff, off_wulff = 0.

    .. attribute:: color_list
        color all input_miller, on_wulff ones according to e_surf

    .. attribute:: color_proxy
        color_proxy for all input_miller

    .. attribute:: color_proxy_on_wulff
        color_proxy for on wulff input_miller only

    .. attribute:: miller_on_wulff
        input_miller_fig for on wulff input_miller only

    .. attribute:: e_surf_on_wulff
        e_surf for on wulff input_miller only

    .. attribute:: miller_area
        ($hkl$): area for all input_miller

    """

    def __init__(self, lattice, miller_list, e_surf_list, color_set='PuBu',
                 grid_off=True, axis_off=True, show_area=False, alpha=1,
                 off_color='red', symprec=0.00001):
        """
        Args:
            lattice: Lattice object of the conventional unit cell
            miller_list: list of hkl or hkil for hcp
            e_surf_list: list of corresponding surface energies
            color_set: default is 'PuBu'
            grid_off(bool): default is True
            axis_off(bool): default is Ture
            show_area(bool): default is False
            alpha: chosen from 0 to 1 (float), default is 1
            off_color: color_legend for off_wulff planes on show_area legend
            symprec: for recp_operation, default is 0.01
        """
        latt = lattice.scale(1)
        structure = Structure(latt, ["H"], [[0, 0, 0]])
        # 1. store input args:
        # store plot settings:
        self.alpha = alpha
        self.color_set = color_set
        self.color_ind = list(range(len(miller_list)))
        self.grid_off = grid_off
        self.axis_off = axis_off
        self.show_area = show_area
        self.off_color = off_color
        self.input_miller_fig = [hkl_tuple_to_str(x) for x in miller_list]
        # store input data
        self.structure = structure
        self.input_miller = [list(x) for x in miller_list]
        self.input_hkl = copy.copy([[x[0], x[1], x[-1]] for x in miller_list])
        self.input_e_surf = copy.copy(e_surf_list)
        self.latt = latt
        self.recp = structure.lattice.reciprocal_lattice_crystallographic
        self.operation = get_recp_symmetry_operation
        self.recp_symmops = self.operation(structure, symprec)
        self.cart_symmops = self.symmop_cartesian(symprec)

        # 2. get all the data for wulff construction
        # get all the surface normal from get_all_miller_e()
        normal_e_m = self.get_all_miller_e()
        # [normal, e_surf, normal_pt, dual_pt, color_plane, m_ind_orig, miller]
        logger.debug(len(normal_e_m))
        self.normal_e_m = normal_e_m

        # 3. consider the dual condition
        dual_pts = [x[3] for x in normal_e_m]
        dual_convex = ConvexHull(dual_pts)
        dual_cv_simp = dual_convex.simplices
        # simplices	(ndarray of ints, shape (nfacet, ndim))
        # list of [i, j, k] , ndim = 3
        # i, j, k: ind for normal_e_m
        # recalculate the dual of dual, get the wulff shape.
        # conner <-> surface
        # get cross point from the simplices of the dual convex hull
        wulff_pt_list = [self.get_cross_pt_dual_simp(dual_simp)
                         for dual_simp in dual_cv_simp]

        wulff_convex = ConvexHull(wulff_pt_list)
        wulff_cv_simp = wulff_convex.simplices
        logger.debug(", ".join([str(len(x)) for x in wulff_cv_simp]))

        # store simplices and convex
        self.dual_cv_simp = dual_cv_simp
        self.wulff_pt_list = wulff_pt_list
        self.wulff_cv_simp = wulff_cv_simp

        # 4. get wulff info
        # return (simpx_info, plane_wulff_info, on_wulff, surface_area)
        wulff_info = self.get_simpx_plane()
        self.simpx_info = wulff_info[0]
        # need to update the color
        # plan_wulff_info: [normal, e_surf, pts, simpx,
        #   color_plane, m_ind_orig, miller]
        self.plane_wulff_info = wulff_info[1]
        self.on_wulff = wulff_info[2]
        self.color_area = wulff_info[3]

        # 5. assign color for on_wulff plane
        # return (color_list, color_proxy, color_proxy_on_wulff,
        # miller_on_wulff, e_surf_on_wulff_list)
        color_info = self.get_colors()
        self.color_list = color_info[0]
        self.color_proxy = color_info[1]
        self.color_proxy_on_wulff = color_info[2]
        self.miller_on_wulff = color_info[3]
        self.e_surf_on_wulff = color_info[4]

        miller_area = []
        for m, in_mill_fig in enumerate(self.input_miller_fig):
            miller_area.append(
                in_mill_fig + ' : ' + str(round(self.color_area[m], 4)))
        self.miller_area = miller_area

    def symmop_cartesian(self, symmprec):
        structure = self.structure
        space_group_analyzer = SpacegroupAnalyzer(structure, symmprec)
        symm_ops = space_group_analyzer.get_point_group_operations(
            cartesian=True)
        return symm_ops

    def get_all_miller_e(self):
        """
        from self:
            get miller_list(unique_miller), e_surf_list and symmetry
            operations(symmops) according to lattice
        apply symmops to get all the miller index, then get normal,
        get all the planes functions for wulff shape calculation:
            |normal| = 1, e_surf is plane's distance to (0, 0, 0),
            normal[0]x + normal[1]y + normal[2]z = e_surf

        return:
            normal_e_m, item: [normal, e_surf, normal_pt, dual_pt,
            color_plane, m_ind_orig, miller]
        """
        all_hkl = copy.copy(self.input_hkl)
        all_hkl_ind = list(enumerate(all_hkl))
        e_surf_list = copy.copy(self.input_e_surf)
        symmops = self.recp_symmops
        recp = self.recp
        color_ind = self.color_ind
        normal_e_m = []
        color = copy.copy(color_ind)
        miller_ind_orig = [x[0] for x in all_hkl_ind]

        for i, hkl in enumerate(all_hkl):
            for op in symmops:
                miller = list(op.operate(hkl))
                miller = [int(x) for x in miller]
                if miller in all_hkl:
                    continue
                else:
                    all_hkl.append(miller)
                    e_surf_list.append(e_surf_list[i])
                    miller_ind_orig.append(i)
                    color.append(color_ind[divmod(i, len(color_ind))[1]])

        for i, hkl in enumerate(all_hkl):
            # get normal (length=1)
            normal = recp.get_cartesian_coords(hkl)
            normal /= scp.linalg.norm(normal)
            e_surf = e_surf_list[i]
            normal_pt = [x * e_surf for x in normal]
            dual_pt = [x / e_surf for x in normal]
            # the index for color and plane
            color_plane = color[i]
            m_ind_orig = miller_ind_orig[i]
            normal_e_m.append([normal, e_surf, normal_pt, dual_pt,
                               color_plane, m_ind_orig, hkl])

        # sorted by e_surf
        normal_e_m.sort(key=lambda x: x[1])
        return normal_e_m

    def get_cross_pt_dual_simp(self, dual_simp):
        """
        |normal| = 1, e_surf is plane's distance to (0, 0, 0),
        plane function:
            normal[0]x + normal[1]y + normal[2]z = e_surf
        normal_e_m, item: [normal, e_surf, normal_pt, dual_pt,
            color_plane, m_ind_orig, miller]

        from self:
            normal_e_m to get the plane functions
            dual_simp: (i, j, k) simplices from the dual convex hull
                i, j, k: plane index(same order in normal_e_m)
        """
        i, j, k = dual_simp[0], dual_simp[1], dual_simp[2]
        normal_e_m = self.normal_e_m
        matrix_surfs = [normal_e_m[i][0], normal_e_m[j][0], normal_e_m[k][0]]
        matrix_e = [normal_e_m[i][1], normal_e_m[j][1], normal_e_m[k][1]]
        cross_pt = scp.dot(scp.linalg.inv(matrix_surfs), matrix_e)
        return cross_pt

    def get_simpx_plane(self):
        """
        local the plane for simpx of on wulff_cv,
        by comparing the center of the simpx triangle
        with the plane functions.

        based on: wulff_cv_simp, normal_e_m, wulff_pt_list
        """
        wulff_simpx = self.wulff_cv_simp
        # normal_e_m, item: [normal, e_surf, normal_pt, dual_pt,
        #     color_plane, m_ind_orig, miller]
        normal_e_m = self.normal_e_m
        wulff_pt_list = self.wulff_pt_list
        simpx_info = []
        on_wulff = [False] * len(self.input_miller)
        surface_area = [0.0] * len(self.input_miller)
        plane_wulff_info = [[x[0], x[1], [], [], x[4], x[5], x[6]] for x in
                            normal_e_m]
        # each simpx (i,j,k) from self.wulff_cv_simp
        #  forms a triangle on the wulff shape.
        # check which surface it belongs to
        for simpx in wulff_simpx:
            pts = [wulff_pt_list[simpx[0]], wulff_pt_list[simpx[1]],
                   wulff_pt_list[simpx[2]]]
            center = np.sum(pts, 0) / 3.0
            # check whether the center of the simplices is on one plane
            for i, plane in enumerate(normal_e_m):
                normal, e_surf = plane[0], plane[1]
                abs_diff = abs(np.dot(normal, center) - e_surf)
                if abs_diff < 10e-6:
                    # assign values only if the simplices is on the plane.
                    plane_ind = plane[4]  # the ind of input_miller
                    on_wulff[plane_ind] = True
                    surface_area[plane_ind] += get_tri_area(pts)
                    simp_inf = [plane[0], plane[1], pts, simpx, plane[4],
                                plane[5], plane[6]]
                    # build a list in the same order of simpx
                    # with related plane information
                    simpx_info.append(simp_inf)

                    plane_wulff_info[i][2].append(pts)
                    # outer_lines
                    plane_wulff_info[i][3].append([simpx[0], simpx[1]])
                    plane_wulff_info[i][3].append([simpx[1], simpx[2]])
                    plane_wulff_info[i][3].append([simpx[0], simpx[2]])
                    # already find the plane, move to the next simplices
                    break
        for plane in plane_wulff_info:
            if len(plane[3]):
                plane[3].sort()
                outer_lines = []
                for line in plane[3]:
                    if plane[3].count(line) == 2:
                        continue
                    outer_lines.append(line)
                plane[3] = outer_lines
        # plan_wulff_info: [normal, e_surf, pts, simpx,
        #     color_plane, m_ind_orig, miller]
        return simpx_info, plane_wulff_info, on_wulff, surface_area

    def get_colors(self):
        """
        assign colors according to the surface energies of on_wulff planes.

        return:
            (color_list, color_proxy, color_proxy_on_wulff, miller_on_wulff,
            e_surf_on_wulff_list)
        """
        on_wulff = self.on_wulff
        input_hkl = self.input_hkl
        input_e_surf = copy.copy(self.input_e_surf)
        input_miller_fig = self.input_miller_fig
        color_set = self.color_set
        alpha = self.alpha
        color_list = [self.off_color] * len(input_hkl)
        color_proxy_on_wulff = []
        miller_on_wulff = []
        e_surf_on_wulff = []
        # get list of on_wulff surface energies (with ind)
        for i, e_surf in enumerate(input_e_surf):
            if on_wulff[i]:
                e_surf_on_wulff.append((i, e_surf))

        c_map = plt.get_cmap(color_set)
        e_surf_on_wulff.sort(key=lambda x: x[1], reverse=False)
        e_surf_on_wulff_list = [x[1] for x in e_surf_on_wulff]
        if len(e_surf_on_wulff) > 1:
            cnorm = colors.Normalize(vmin=min(e_surf_on_wulff_list),
                                     vmax=max(e_surf_on_wulff_list))
        else:
            # if there is only one hkl on wulff, choose the color of the median
            cnorm = colors.Normalize(vmin=min(e_surf_on_wulff_list) - 0.1,
                                     vmax=max(e_surf_on_wulff_list) + 0.1)
        scalar_map = cm.ScalarMappable(norm=cnorm, cmap=c_map)
        logger.debug(e_surf_on_wulff)
        # prepare for
        for i, e_surf in e_surf_on_wulff:
            plane_color = scalar_map.to_rgba(e_surf, alpha=alpha)
            color_list[i] = plane_color
            color_proxy_on_wulff.append(
                plt.Rectangle((2, 2), 1, 1, fc=plane_color, alpha=alpha))
            miller_on_wulff.append(input_miller_fig[i])
        scalar_map.set_array([x[1] for x in e_surf_on_wulff])
        color_proxy = [plt.Rectangle((2, 2), 1, 1, fc=x, alpha=alpha)
                       for x in color_list]

        return color_list, color_proxy, color_proxy_on_wulff, miller_on_wulff, \
            e_surf_on_wulff_list

    def plot_wf_simpx(self, direction=None, bar_pos=(0.75, 0.15, 0.05, 0.65),
                      bar_on=False, legend_on=True, aspect_ratio=(8, 8)):
        """
        plot the wulff shape from self.wulff_pt_list, self.plane_wulff_info

        Args:
            direction: default is (1, 1, 1)
            bar_pos: default is [0.75, 0.15, 0.05, 0.65]
            bar_on (bool): default is False
            legend_on (bool): default is True
            aspect_ratio: default is (8, 8)

        """

        if not direction:
            area_dict = self.area_fraction_dict
            miller_indices, areas = [], []
            for hkl in area_dict.keys():
                miller_indices.append(hkl)
                areas.append(area_dict[hkl])
            direction = miller_indices[areas.index(max(areas))]
        fig = plt.figure()
        fig.set_size_inches(aspect_ratio[0],
                            aspect_ratio[1])
        azim, elev = self.get_azimuth_elev(direction)

        wulff_pt_list = self.wulff_pt_list
        plane_wulff_info = self.plane_wulff_info
        # [normal, e_surf, [pts], [simpx],
        #     color_plane, m_ind_orig, miller]
        ax = a3.Axes3D(fig, azim=azim, elev=elev)

        for plane in plane_wulff_info:
            # check whether [pts] is empty
            if len(plane[2]) < 1:
                # empty, plane is not on_wulff.
                continue
            # assign the color for on_wulff planes according to its
            # color_plane index and the color_list for on_wulff
            plane_color = self.color_list[plane[4]]
            # plane[3]: [simpx]
            lines = list(plane[3])
            pt = []
            prev = None
            while len(lines) > 0:
                if prev is None:
                    l = lines.pop(0)
                else:
                    for i, l in enumerate(lines):
                        if prev in l:
                            l = lines.pop(i)
                            if l[1] == prev:
                                l.reverse()
                            break
                # make sure the lines are connected one by one.
                # find the way covering all pts and facets
                pt.append(self.wulff_pt_list[l[0]].tolist())
                pt.append(self.wulff_pt_list[l[1]].tolist())
                prev = l[1]
            # plot from the sorted pts from [simpx]
            tri = a3.art3d.Poly3DCollection([pt])
            tri.set_color(plane_color)
            # "#808080" is the default edge color.
            tri.set_edgecolor("#808080")
            ax.add_collection3d(tri)
        # set ranges of x, y, z
        # find the largest distance between on_wulff pts and the origin,
        # to ensure complete and consistent display for all directions
        r_range = max([np.linalg.norm(x) for x in wulff_pt_list])
        ax.set_xlim([-r_range * 1.1, r_range * 1.1])
        ax.set_ylim([-r_range * 1.1, r_range * 1.1])
        ax.set_zlim([-r_range * 1.1, r_range * 1.1])
        # add legend
        if legend_on:
            color_proxy = self.color_proxy
            if self.show_area:
                ax.legend(color_proxy, self.miller_area, loc='upper left',
                          bbox_to_anchor=(0, 1), fancybox=True, shadow=False)
            else:
                ax.legend(self.color_proxy_on_wulff, self.miller_on_wulff,
                          loc='upper center',
                          bbox_to_anchor=(0.5, 1), ncol=3, fancybox=True,
                          shadow=False)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        # add colorbar
        # cbar = fig.colorbar(self.scalarcm, alpha=self.alpha)
        # Add an axes at position rect [left, bottom, width, height]
        cmap = plt.get_cmap(self.color_set)
        cmap.set_over('0.25')
        cmap.set_under('0.75')
        bounds = [round(e, 2) for e in self.e_surf_on_wulff]
        bounds.append(1.2 * bounds[-1])
        norm = colors.BoundaryNorm(bounds, cmap.N)
        if bar_on:
            # display surface energies
            ax1 = fig.add_axes(bar_pos)
            cbar = colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm,
                                         boundaries=[0] + bounds + [10],
                                         extend='both',
                                         ticks=bounds[:-1],  # optional
                                         spacing='proportional',
                                         orientation='vertical'
                                         )
            cbar.set_label('Surface Energies ($J/m^2$)', fontsize=100)
        # [normal, e_surf, normal_pt, dual_pt, color_plane, m_ind_orig, miller]
        if self.grid_off:
            ax.grid('off')
        if self.axis_off:
            ax.axis('off')
        plt.draw()
        return plt

    def get_azimuth_elev(self, miller_index):
        """
        :param miller_index: viewing direction
        :return: azim, elev for plotting
        """

        cart = self.latt.get_cartesian_coords(miller_index)
        azim = get_angle([cart[0], cart[1], 0], (1, 0, 0))
        v = [cart[0], cart[1], 0]
        elev = get_angle(cart, v)
        if miller_index == (0, 0, 1) or miller_index == (0, 0, 0, 1):
            return 0, 90
        else:
            return azim, elev

    @property
    def wulff_volume(self):
        """
        :return:
            the volume of the wulff shape
        """
        return ConvexHull(self.wulff_pt_list).volume

    @property
    def miller_area_dict(self):
        """
        :return:
            (dict): {hkl: area_hkl on wulff}
        """
        miller_area_dict = {}
        for i, hkl in enumerate(self.input_miller):
            miller_area_dict[tuple(hkl)] = self.color_area[i]
        return miller_area_dict

    @property
    def miller_energy_dict(self):
        """
        :return:
            (dict): {hkl: surface energy_hkl}
        """
        miller_energy_dict = {}
        for i, hkl in enumerate(self.input_miller):
            miller_energy_dict[tuple(hkl)] = self.input_e_surf[i]
        return miller_energy_dict

    @property
    def total_surface_area(self):
        """
        :return:
            total area on wulff
        """
        tot_area = 0
        for hkl in self.miller_area_dict.keys():
            tot_area += self.miller_area_dict[hkl]
        return tot_area

    @property
    def weighted_surface_energy(self):
        """
        :return:
            sum(surface_energy_hkl * area_hkl)/ sum(area_hkl)
        """
        tot_area_energy = 0
        for hkl in self.miller_energy_dict.keys():
            tot_area_energy += self.miller_energy_dict[hkl] * \
                               self.miller_area_dict[hkl]
        return tot_area_energy / self.total_surface_area

    @property
    def area_fraction_dict(self):
        """
        :return:
            (dict): {hkl: area_hkl/total area on wulff}
        """
        area_fraction_dict = {}
        for hkl in self.miller_area_dict.keys():
            area_fraction_dict[hkl] = self.miller_area_dict[hkl] / \
                                      self.total_surface_area
        return area_fraction_dict

    @property
    def anisotropy(self):
        """
        :return:
            variation from weighted surface energy
            The ideal sphere is 0.
        """
        square_diff_energy = 0
        weighted_energy = self.weighted_surface_energy
        area_frac_dict = self.area_fraction_dict
        miller_energy_dict = self.miller_energy_dict

        for hkl in miller_energy_dict.keys():
            square_diff_energy += ((miller_energy_dict[
                                        hkl] - weighted_energy) ** 2) * \
                                  area_frac_dict[hkl]
        return np.sqrt(square_diff_energy) / weighted_energy

    @property
    def shape_factor(self):
        """
        This is useful for determining the critical nucleus size.
        A large shape factor indicates great anisotropy.
        See Ballufi, R. W., Allen, S. M. & Carter, W. C. Kinetics
            of Materials. (John Wiley & Sons, 2005), p.461

        :return:
            variation from weighted surface energy
            The ideal sphere is 0.
        """
        return self.total_surface_area / (self.wulff_volume ** (2 / 3))
