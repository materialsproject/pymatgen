# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module provides classes for plotting PhaseDiagram objects.
"""

from six.moves import zip

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Jun 15, 2012"

import math
import numpy as np
import itertools

from pymatgen.phasediagram.analyzer import PDAnalyzer
from pymatgen.util.string_utils import latexify
from pymatgen.util.plotting_utils import get_publication_quality_plot
from pymatgen.util.coord_utils import in_coord_list


class PDPlotter(object):
    """
    A plotter class for phase diagrams.

    Args:
        phasediagram: PhaseDiagram object.
        show_unstable (float): Whether unstable phases will be plotted as
            well as red crosses. If a number > 0 is entered, all phases with
            ehull < show_unstable will be shown.
    """

    def __init__(self, phasediagram, show_unstable=0):
        self._pd = phasediagram
        self._dim = len(self._pd.elements)
        if self._dim > 4:
            raise ValueError("Only 1-4 components supported!")
        self.lines = uniquelines(self._pd.facets) if self._dim > 1 else \
            [[self._pd.facets[0][0], self._pd.facets[0][0]]]
        self.show_unstable = show_unstable

    @property
    def pd_plot_data(self):
        """
        Plot data for phase diagram.
        2-comp - Full hull with energies
        3/4-comp - Projection into 2D or 3D Gibbs triangle.

        Returns:
            (lines, stable_entries, unstable_entries):
            - lines is a list of list of coordinates for lines in the PD.
            - stable_entries is a {coordinate : entry} for each stable node
            in the phase diagram. (Each coordinate can only have one
            stable phase)
            - unstable_entries is a {entry: coordinates} for all unstable
            nodes in the phase diagram.
        """
        pd = self._pd
        entries = pd.qhull_entries
        data = np.array(pd.qhull_data)
        lines = []
        stable_entries = {}
        for line in self.lines:
            entry1 = entries[line[0]]
            entry2 = entries[line[1]]
            if self._dim < 3:
                x = [data[line[0]][0], data[line[1]][0]]
                y = [pd.get_form_energy_per_atom(entry1),
                     pd.get_form_energy_per_atom(entry2)]
                coord = [x, y]
            elif self._dim == 3:
                coord = triangular_coord(data[line, 0:2])
            else:
                coord = tet_coord(data[line, 0:3])
            lines.append(coord)
            labelcoord = list(zip(*coord))
            stable_entries[labelcoord[0]] = entry1
            stable_entries[labelcoord[1]] = entry2

        all_entries = pd.all_entries
        all_data = np.array(pd.all_entries_hulldata)
        unstable_entries = dict()
        stable = pd.stable_entries
        for i in range(0, len(all_entries)):
            entry = all_entries[i]
            if entry not in stable:
                if self._dim < 3:
                    x = [all_data[i][0], all_data[i][0]]
                    y = [pd.get_form_energy_per_atom(entry),
                         pd.get_form_energy_per_atom(entry)]
                    coord = [x, y]
                elif self._dim == 3:
                    coord = triangular_coord([all_data[i, 0:2],
                                              all_data[i, 0:2]])
                else:
                    coord = tet_coord([all_data[i, 0:3], all_data[i, 0:3],
                                       all_data[i, 0:3]])
                labelcoord = list(zip(*coord))
                unstable_entries[entry] = labelcoord[0]

        return lines, stable_entries, unstable_entries

    def get_plot(self, label_stable=True, label_unstable=True, ordering=None,
                 energy_colormap=None, process_attributes=False):
        if self._dim < 4:
            plt = self._get_2d_plot(label_stable, label_unstable, ordering,
                                    energy_colormap,
                                    process_attributes=process_attributes)
        elif self._dim == 4:
            plt = self._get_3d_plot(label_stable)

        return plt

    def show(self, *args, **kwargs):
        """
        Draws the phase diagram using Matplotlib and show it.

        Args:
            \*args: Passed to get_plot.
            \*\*kwargs: Passed to get_plot.
        """
        self.get_plot(*args, **kwargs).show()

    def _get_2d_plot(self, label_stable=True, label_unstable=True,
                     ordering=None, energy_colormap=None, vmin_mev=-60.0,
                     vmax_mev=60.0, show_colorbar=True,
                     process_attributes=False):
        """
        Shows the plot using pylab.  Usually I won't do imports in methods,
        but since plotting is a fairly expensive library to load and not all
        machines have matplotlib installed, I have done it this way.
        """

        plt = get_publication_quality_plot(8, 6)
        from matplotlib.font_manager import FontProperties
        if ordering is None:
            (lines, labels, unstable) = self.pd_plot_data
        else:
            (_lines, _labels, _unstable) = self.pd_plot_data
            (lines, labels, unstable) = order_phase_diagram(
                _lines, _labels, _unstable, ordering)
        if energy_colormap is None:
            if process_attributes:
                for x, y in lines:
                    plt.plot(x, y, "k-", linewidth=3, markeredgecolor="k")
                # One should think about a clever way to have "complex"
                # attributes with complex processing options but with a clear
                #  logic. At this moment, I just use the attributes to know
                # whether an entry is a new compound or an existing (from the
                #  ICSD or from the MP) one.
                for x, y in labels.keys():
                    if labels[(x, y)].attribute is None or \
                            labels[(x, y)].attribute == "existing":
                        plt.plot(x, y, "ko", linewidth=3, markeredgecolor="k",
                                 markerfacecolor="b", markersize=12)
                    else:
                        plt.plot(x, y, "k*", linewidth=3, markeredgecolor="k",
                                 markerfacecolor="g", markersize=18)
            else:
                for x, y in lines:
                    plt.plot(x, y, "ko-", linewidth=3, markeredgecolor="k",
                             markerfacecolor="b", markersize=15)
        else:
            from matplotlib.colors import Normalize, LinearSegmentedColormap
            from matplotlib.cm import ScalarMappable
            pda = PDAnalyzer(self._pd)
            for x, y in lines:
                plt.plot(x, y, "k-", linewidth=3, markeredgecolor="k")
            vmin = vmin_mev / 1000.0
            vmax = vmax_mev / 1000.0
            if energy_colormap == 'default':
                mid = - vmin / (vmax - vmin)
                cmap = LinearSegmentedColormap.from_list(
                    'my_colormap', [(0.0, '#005500'), (mid, '#55FF55'),
                                    (mid, '#FFAAAA'), (1.0, '#FF0000')])
            else:
                cmap = energy_colormap
            norm = Normalize(vmin=vmin, vmax=vmax)
            _map = ScalarMappable(norm=norm, cmap=cmap)
            _energies = [pda.get_equilibrium_reaction_energy(entry)
                         for coord, entry in labels.items()]
            energies = [en if en < 0.0 else -0.00000001 for en in _energies]
            vals_stable = _map.to_rgba(energies)
            ii = 0
            if process_attributes:
                for x, y in labels.keys():
                    if labels[(x, y)].attribute is None or \
                            labels[(x, y)].attribute == "existing":
                        plt.plot(x, y, "o", markerfacecolor=vals_stable[ii],
                                 markersize=12)
                    else:
                        plt.plot(x, y, "*", markerfacecolor=vals_stable[ii],
                                 markersize=18)
                    ii += 1
            else:
                for x, y in labels.keys():
                    plt.plot(x, y, "o", markerfacecolor=vals_stable[ii],
                             markersize=15)
                    ii += 1

        font = FontProperties()
        font.set_weight("bold")
        font.set_size(24)

        # Sets a nice layout depending on the type of PD. Also defines a
        # "center" for the PD, which then allows the annotations to be spread
        # out in a nice manner.
        if len(self._pd.elements) == 3:
            plt.axis("equal")
            plt.xlim((-0.1, 1.2))
            plt.ylim((-0.1, 1.0))
            plt.axis("off")
            center = (0.5, math.sqrt(3) / 6)
        else:
            all_coords = labels.keys()
            miny = min([c[1] for c in all_coords])
            ybuffer = max(abs(miny) * 0.1, 0.1)
            plt.xlim((-0.1, 1.1))
            plt.ylim((miny - ybuffer, ybuffer))
            center = (0.5, miny / 2)
            plt.xlabel("Fraction", fontsize=28, fontweight='bold')
            plt.ylabel("Formation energy (eV/fu)", fontsize=28,
                       fontweight='bold')

        for coords in sorted(labels.keys(), key=lambda x: -x[1]):
            entry = labels[coords]
            label = entry.name

            # The follow defines an offset for the annotation text emanating
            # from the center of the PD. Results in fairly nice layouts for the
            # most part.
            vec = (np.array(coords) - center)
            vec = vec / np.linalg.norm(vec) * 10 if np.linalg.norm(vec) != 0 \
                else vec
            valign = "bottom" if vec[1] > 0 else "top"
            if vec[0] < -0.01:
                halign = "right"
            elif vec[0] > 0.01:
                halign = "left"
            else:
                halign = "center"
            if label_stable:
                if process_attributes and entry.attribute == 'new':
                    plt.annotate(latexify(label), coords, xytext=vec,
                                 textcoords="offset points",
                                 horizontalalignment=halign,
                                 verticalalignment=valign,
                                 fontproperties=font,
                                 color='g')
                else:
                    plt.annotate(latexify(label), coords, xytext=vec,
                                 textcoords="offset points",
                                 horizontalalignment=halign,
                                 verticalalignment=valign,
                                 fontproperties=font)

        if self.show_unstable:
            font = FontProperties()
            font.set_size(16)
            pda = PDAnalyzer(self._pd)
            energies_unstable = [pda.get_e_above_hull(entry)
                                 for entry, coord in unstable.items()]
            if energy_colormap is not None:
                energies.extend(energies_unstable)
                vals_unstable = _map.to_rgba(energies_unstable)
            ii = 0
            for entry, coords in unstable.items():
                ehull = pda.get_e_above_hull(entry)
                if ehull < self.show_unstable:
                    vec = (np.array(coords) - center)
                    vec = vec / np.linalg.norm(vec) * 10 \
                        if np.linalg.norm(vec) != 0 else vec
                    label = entry.name
                    if energy_colormap is None:
                        plt.plot(coords[0], coords[1], "ks", linewidth=3,
                                 markeredgecolor="k", markerfacecolor="r",
                                 markersize=8)
                    else:
                        plt.plot(coords[0], coords[1], "s", linewidth=3,
                                 markeredgecolor="k",
                                 markerfacecolor=vals_unstable[ii],
                                 markersize=8)
                    if label_unstable:
                        plt.annotate(latexify(label), coords, xytext=vec,
                                     textcoords="offset points",
                                     horizontalalignment=halign, color="b",
                                     verticalalignment=valign,
                                     fontproperties=font)
                    ii += 1
        if energy_colormap is not None and show_colorbar:
            _map.set_array(energies)
            cbar = plt.colorbar(_map)
            cbar.set_label(
                'Energy [meV/at] above hull (in red)\nInverse energy ['
                'meV/at] above hull (in green)',
                rotation=-90, ha='left', va='center')
            ticks = cbar.ax.get_yticklabels()
            # cbar.ax.set_yticklabels(['${v}$'.format(
            #     v=float(t.get_text().strip('$'))*1000.0) for t in ticks])
        f = plt.gcf()
        f.set_size_inches((8, 6))
        plt.subplots_adjust(left=0.09, right=0.98, top=0.98, bottom=0.07)
        return plt

    def _get_3d_plot(self, label_stable=True):
        """
        Shows the plot using pylab.  Usually I won"t do imports in methods,
        but since plotting is a fairly expensive library to load and not all
        machines have matplotlib installed, I have done it this way.
        """
        import matplotlib.pyplot as plt
        import mpl_toolkits.mplot3d.axes3d as p3
        from matplotlib.font_manager import FontProperties
        fig = plt.figure()
        ax = p3.Axes3D(fig)
        font = FontProperties()
        font.set_weight("bold")
        font.set_size(20)
        (lines, labels, unstable) = self.pd_plot_data
        count = 1
        newlabels = list()
        for x, y, z in lines:
            ax.plot(x, y, z, "bo-", linewidth=3, markeredgecolor="b",
                    markerfacecolor="r", markersize=10)
        for coords in sorted(labels.keys()):
            entry = labels[coords]
            label = entry.name
            if label_stable:
                if len(entry.composition.elements) == 1:
                    ax.text(coords[0], coords[1], coords[2], label)
                else:
                    ax.text(coords[0], coords[1], coords[2], str(count))
                    newlabels.append("{} : {}".format(count, latexify(label)))
                    count += 1
        plt.figtext(0.01, 0.01, "\n".join(newlabels))
        ax.axis("off")
        return plt

    def write_image(self, stream, image_format="svg", label_stable=True,
                    label_unstable=True, ordering=None,
                    energy_colormap=None, process_attributes=False):
        """
        Writes the phase diagram to an image in a stream.

        Args:
            stream:
                stream to write to. Can be a file stream or a StringIO stream.
            image_format
                format for image. Can be any of matplotlib supported formats.
                Defaults to svg for best results for vector graphics.
        """
        plt = self.get_plot(
            label_stable=label_stable, label_unstable=label_unstable,
            ordering=ordering, energy_colormap=energy_colormap,
            process_attributes=process_attributes)

        f = plt.gcf()
        f.set_size_inches((12, 10))

        plt.savefig(stream, format=image_format)

    def plot_chempot_range_map(self, elements, referenced=True):
        """
        Plot the chemical potential range _map. Currently works only for
        3-component PDs.

        Args:
            elements: Sequence of elements to be considered as independent
                variables. E.g., if you want to show the stability ranges of
                all Li-Co-O phases wrt to uLi and uO, you will supply
                [Element("Li"), Element("O")]
            referenced: if True, gives the results with a reference being the
                        energy of the elemental phase. If False, gives absolute values.
        """
        self.get_chempot_range_map_plot(elements, referenced=referenced).show()

    def get_chempot_range_map_plot(self, elements,referenced=True):
        """
        Returns a plot of the chemical potential range _map. Currently works
        only for 3-component PDs.

        Args:
            elements: Sequence of elements to be considered as independent
                variables. E.g., if you want to show the stability ranges of
                all Li-Co-O phases wrt to uLi and uO, you will supply
                [Element("Li"), Element("O")]
            referenced: if True, gives the results with a reference being the
                        energy of the elemental phase. If False, gives absolute values.

        Returns:
            A matplotlib plot object.
        """

        plt = get_publication_quality_plot(12, 8)
        analyzer = PDAnalyzer(self._pd)
        chempot_ranges = analyzer.get_chempot_range_map(
            elements, referenced=referenced)
        missing_lines = {}
        excluded_region = []
        for entry, lines in chempot_ranges.items():
            comp = entry.composition
            center_x = 0
            center_y = 0
            coords = []
            contain_zero = any([comp.get_atomic_fraction(el) == 0
                                for el in elements])
            is_boundary = (not contain_zero) and \
                sum([comp.get_atomic_fraction(el) for el in elements]) == 1
            for line in lines:
                (x, y) = line.coords.transpose()
                plt.plot(x, y, "k-")

                for coord in line.coords:
                    if not in_coord_list(coords, coord):
                        coords.append(coord.tolist())
                        center_x += coord[0]
                        center_y += coord[1]
                if is_boundary:
                    excluded_region.extend(line.coords)

            if coords and contain_zero:
                missing_lines[entry] = coords
            else:
                xy = (center_x / len(coords), center_y / len(coords))
                plt.annotate(latexify(entry.name), xy, fontsize=22)

        ax = plt.gca()
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        #Shade the forbidden chemical potential regions.
        excluded_region.append([xlim[1], ylim[1]])
        excluded_region = sorted(excluded_region, key=lambda c: c[0])
        (x, y) = np.transpose(excluded_region)
        plt.fill(x, y, "0.80")

        #The hull does not generate the missing horizontal and vertical lines.
        #The following code fixes this.
        el0 = elements[0]
        el1 = elements[1]
        for entry, coords in missing_lines.items():
            center_x = sum([c[0] for c in coords])
            center_y = sum([c[1] for c in coords])
            comp = entry.composition
            is_x = comp.get_atomic_fraction(el0) < 0.01
            is_y = comp.get_atomic_fraction(el1) < 0.01
            n = len(coords)
            if not (is_x and is_y):
                if is_x:
                    coords = sorted(coords, key=lambda c: c[1])
                    for i in [0, -1]:
                        x = [min(xlim), coords[i][0]]
                        y = [coords[i][1], coords[i][1]]
                        plt.plot(x, y, "k")
                        center_x += min(xlim)
                        center_y += coords[i][1]
                elif is_y:
                    coords = sorted(coords, key=lambda c: c[0])
                    for i in [0, -1]:
                        x = [coords[i][0], coords[i][0]]
                        y = [coords[i][1], min(ylim)]
                        plt.plot(x, y, "k")
                        center_x += coords[i][0]
                        center_y += min(ylim)
                xy = (center_x / (n + 2), center_y / (n + 2))
            else:
                center_x = sum(coord[0] for coord in coords) + xlim[0]
                center_y = sum(coord[1] for coord in coords) + ylim[0]
                xy = (center_x / (n + 1), center_y / (n + 1))

            plt.annotate(latexify(entry.name), xy,
                         horizontalalignment="center",
                         verticalalignment="center", fontsize=22)

        plt.xlabel("$\mu_{{{0}}} - \mu_{{{0}}}^0$ (eV)"
                   .format(el0.symbol))
        plt.ylabel("$\mu_{{{0}}} - \mu_{{{0}}}^0$ (eV)"
                   .format(el1.symbol))
        plt.tight_layout()
        return plt

    def get_contour_pd_plot(self):
        """
        Plot a contour phase diagram plot, where phase triangles are colored
        according to degree of instability by interpolation. Currently only
        works for 3-component phase diagrams.

        Returns:
            A matplotlib plot object.
        """
        from scipy import interpolate
        from matplotlib import cm

        pd = self._pd
        entries = pd.qhull_entries
        data = np.array(pd.qhull_data)

        plt = self._get_2d_plot()
        analyzer = PDAnalyzer(pd)
        data[:, 0:2] = triangular_coord(data[:, 0:2]).transpose()
        for i, e in enumerate(entries):
            data[i, 2] = analyzer.get_e_above_hull(e)

        gridsize = 0.005
        xnew = np.arange(0, 1., gridsize)
        ynew = np.arange(0, 1, gridsize)

        f = interpolate.LinearNDInterpolator(data[:, 0:2], data[:, 2])
        znew = np.zeros((len(ynew), len(xnew)))
        for (i, xval) in enumerate(xnew):
            for (j, yval) in enumerate(ynew):
                znew[j, i] = f(xval, yval)

        plt.contourf(xnew, ynew, znew, 1000, cmap=cm.autumn_r)

        plt.colorbar()
        return plt


def uniquelines(q):
    """
    Given all the facets, convert it into a set of unique lines.  Specifically
    used for converting convex hull facets into line pairs of coordinates.

    Args:
        q: A 2-dim sequence, where each row represents a facet. E.g.,
            [[1,2,3],[3,6,7],...]

    Returns:
        setoflines:
            A set of tuple of lines.  E.g., ((1,2), (1,3), (2,3), ....)
    """
    setoflines = set()
    for facets in q:
        for line in itertools.combinations(facets, 2):
            setoflines.add(tuple(sorted(line)))
    return setoflines


def triangular_coord(coord):
    """
    Convert a 2D coordinate into a triangle-based coordinate system for a
    prettier phase diagram.

    Args:
        coordinate: coordinate used in the convex hull computation.

    Returns:
        coordinates in a triangular-based coordinate system.
    """
    unitvec = np.array([[1, 0], [0.5, math.sqrt(3) / 2]])
    result = np.dot(np.array(coord), unitvec)
    return result.transpose()


def tet_coord(coord):
    """
    Convert a 3D coordinate into a tetrahedron based coordinate system for a
    prettier phase diagram.

    Args:
        coordinate: coordinate used in the convex hull computation.

    Returns:
        coordinates in a tetrahedron-based coordinate system.
    """
    unitvec = np.array([[1, 0, 0], [0.5, math.sqrt(3) / 2, 0],
                        [0.5, 1.0 / 3.0 * math.sqrt(3) / 2, math.sqrt(6) / 3]])
    result = np.dot(np.array(coord), unitvec)
    return result.transpose()


def order_phase_diagram(lines, stable_entries, unstable_entries, ordering):
    """
    Orders the entries (their coordinates) in a phase diagram plot according
    to the user specified ordering.
    Ordering should be given as ['Up', 'Left', 'Right'], where Up,
    Left and Right are the names of the entries in the upper, left and right
    corners of the triangle respectively.

    Args:
        lines: list of list of coordinates for lines in the PD.
        stable_entries: {coordinate : entry} for each stable node in the
            phase diagram. (Each coordinate can only have one stable phase)
        unstable_entries: {entry: coordinates} for all unstable nodes in the
            phase diagram.
        ordering: Ordering of the phase diagram, given as a list ['Up',
            'Left','Right']

    Returns:
        (newlines, newstable_entries, newunstable_entries):
        - newlines is a list of list of coordinates for lines in the PD.
        - newstable_entries is a {coordinate : entry} for each stable node
        in the phase diagram. (Each coordinate can only have one
        stable phase)
        - newunstable_entries is a {entry: coordinates} for all unstable
        nodes in the phase diagram.
    """
    yup = -1000.0
    xleft = 1000.0
    xright = -1000.0

    for coord in stable_entries:
        if coord[0] > xright:
            xright = coord[0]
            nameright = stable_entries[coord].name
        if coord[0] < xleft:
            xleft = coord[0]
            nameleft = stable_entries[coord].name
        if coord[1] > yup:
            yup = coord[1]
            nameup = stable_entries[coord].name

    if (not nameup in ordering) or (not nameright in ordering) or \
            (not nameleft in ordering):
        raise ValueError(
            'Error in ordering_phase_diagram : \n"{up}", "{left}" and "{'
            'right}"'
            ' should be in ordering : {ord}'.format(up=nameup, left=nameleft,
                                                    right=nameright,
                                                    ord=ordering))

    cc = np.array([0.5, np.sqrt(3.0) / 6.0], np.float)

    if nameup == ordering[0]:
        if nameleft == ordering[1]:
            # The coordinates were already in the user ordering
            return lines, stable_entries, unstable_entries
        else:
            newlines = [[np.array(1.0 - x), y] for x, y in lines]
            newstable_entries = {(1.0 - c[0], c[1]): entry
                                 for c, entry in stable_entries.items()}
            newunstable_entries = {entry: (1.0 - c[0], c[1])
                                   for entry, c in
                                   unstable_entries.items()}
            return newlines, newstable_entries, newunstable_entries
    elif nameup == ordering[1]:
        if nameleft == ordering[2]:
            c120 = np.cos(2.0 * np.pi / 3.0)
            s120 = np.sin(2.0 * np.pi / 3.0)
            newlines = []
            for x, y in lines:
                newx = np.zeros_like(x)
                newy = np.zeros_like(y)
                for ii, xx in enumerate(x):
                    newx[ii] = c120 * (xx - cc[0]) - s120 * (y[ii] - cc[1]) + \
                               cc[0]
                    newy[ii] = s120 * (xx - cc[0]) + c120 * (y[ii] - cc[1]) + \
                               cc[1]
                newlines.append([newx, newy])
            newstable_entries = {
                (c120 * (c[0] - cc[0]) - s120 * (c[1] - cc[1]) + cc[0],
                 s120 * (c[0] - cc[0]) + c120 * (c[1] - cc[1]) + cc[1]): entry
                for c, entry in stable_entries.items()}
            newunstable_entries = {
                entry: (c120 * (c[0] - cc[0]) - s120 * (c[1] - cc[1]) + cc[0],
                        s120 * (c[0] - cc[0]) + c120 * (c[1] - cc[1]) + cc[1])
                for entry, c in unstable_entries.items()}
            return newlines, newstable_entries, newunstable_entries
        else:
            c120 = np.cos(2.0 * np.pi / 3.0)
            s120 = np.sin(2.0 * np.pi / 3.0)
            newlines = []
            for x, y in lines:
                newx = np.zeros_like(x)
                newy = np.zeros_like(y)
                for ii, xx in enumerate(x):
                    newx[ii] = -c120 * (xx - 1.0) - s120 * y[ii] + 1.0
                    newy[ii] = -s120 * (xx - 1.0) + c120 * y[ii]
                newlines.append([newx, newy])
            newstable_entries = {(-c120 * (c[0] - 1.0) - s120 * c[1] + 1.0,
                                  -s120 * (c[0] - 1.0) + c120 * c[1]): entry
                                 for c, entry in stable_entries.items()}
            newunstable_entries = {
            entry: (-c120 * (c[0] - 1.0) - s120 * c[1] + 1.0,
                    -s120 * (c[0] - 1.0) + c120 * c[1])
            for entry, c in unstable_entries.items()}
            return newlines, newstable_entries, newunstable_entries
    elif nameup == ordering[2]:
        if nameleft == ordering[0]:
            c240 = np.cos(4.0 * np.pi / 3.0)
            s240 = np.sin(4.0 * np.pi / 3.0)
            newlines = []
            for x, y in lines:
                newx = np.zeros_like(x)
                newy = np.zeros_like(y)
                for ii, xx in enumerate(x):
                    newx[ii] = c240 * (xx - cc[0]) - s240 * (y[ii] - cc[1]) + \
                               cc[0]
                    newy[ii] = s240 * (xx - cc[0]) + c240 * (y[ii] - cc[1]) + \
                               cc[1]
                newlines.append([newx, newy])
            newstable_entries = {
                (c240 * (c[0] - cc[0]) - s240 * (c[1] - cc[1]) + cc[0],
                 s240 * (c[0] - cc[0]) + c240 * (c[1] - cc[1]) + cc[1]): entry
                for c, entry in stable_entries.items()}
            newunstable_entries = {
                entry: (c240 * (c[0] - cc[0]) - s240 * (c[1] - cc[1]) + cc[0],
                        s240 * (c[0] - cc[0]) + c240 * (c[1] - cc[1]) + cc[1])
                for entry, c in unstable_entries.items()}
            return newlines, newstable_entries, newunstable_entries
        else:
            c240 = np.cos(4.0 * np.pi / 3.0)
            s240 = np.sin(4.0 * np.pi / 3.0)
            newlines = []
            for x, y in lines:
                newx = np.zeros_like(x)
                newy = np.zeros_like(y)
                for ii, xx in enumerate(x):
                    newx[ii] = -c240 * xx - s240 * y[ii]
                    newy[ii] = -s240 * xx + c240 * y[ii]
                newlines.append([newx, newy])
            newstable_entries = {(-c240 * c[0] - s240 * c[1],
                                  -s240 * c[0] + c240 * c[1]): entry
                                 for c, entry in stable_entries.items()}
            newunstable_entries = {entry: (-c240 * c[0] - s240 * c[1],
                                           -s240 * c[0] + c240 * c[1])
                                   for entry, c in unstable_entries.items()}
            return newlines, newstable_entries, newunstable_entries
