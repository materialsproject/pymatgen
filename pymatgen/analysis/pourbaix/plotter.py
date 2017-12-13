# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import platform

from matplotlib.font_manager import FontProperties
from matplotlib.ticker import FormatStrFormatter

"""
This module provides classes for plotting Pourbaix objects.
"""

import six
from six.moves import map
from six.moves import zip
import warnings

__author__ = "Sai Jayaraman"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Sai Jayaraman"
__email__ = "sjayaram@mit.edu"
__status__ = "Production"
__date__ = "Jan 26, 2012"

import numpy as np
import re
import collections


from pymatgen.analysis.pourbaix.analyzer import PourbaixAnalyzer
from pymatgen.analysis.pourbaix.maker import PREFAC
from pymatgen.analysis.pourbaix.entry import MultiEntry

from pymatgen.analysis.phase_diagram import uniquelines
from pymatgen.util.string import latexify
from pymatgen.util.plotting import pretty_plot
from pymatgen.util.coord import in_coord_list


class PourbaixPlotter(object):
    """
    A plotter class for phase diagrams.

    Args:
        phasediagram: A PhaseDiagram object.
        show_unstable: Whether unstable phases will be plotted as well as
            red crosses. Defaults to False.
    """

    def __init__(self, pourbaixdiagram, show_unstable=False):
        self._pd = pourbaixdiagram
        self.lines = uniquelines(self._pd.facets)
        self.show_unstable = show_unstable

    @property
    def pourbaix_hull_plot_data(self):
        """
        Pourbaix diagram convex hull data.

        Returns:
            (lines, stable_entries, unstable_entries)
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
        facetlines = self.lines
        lines = list()
        stable_entries = dict()
        for line in facetlines:
            entry1 = entries[line[0]]
            entry2 = entries[line[1]]
            x = [data[line[0]][0], data[line[1]][0]]
            y = [data[line[0]][1], data[line[1]][1]]
            z = [data[line[0]][2], data[line[1]][2]]
            coord = [x, y, z]
            lines.append(coord)
            labelcoord = list(zip(*coord))
            stable_entries[labelcoord[0]] = entry1
            stable_entries[labelcoord[1]] = entry2

        allentries = pd.all_entries
        alldata = np.array(pd.qhull_data)
        unstable_entries = dict()
        stable = pd.stable_entries
        for i in range(len(allentries)):
            entry = allentries[i]
            if entry not in stable:
                x = [alldata[i][0], alldata[i][0]]
                y = [alldata[i][1], alldata[i][1]]
                z = [alldata[i][2], alldata[i][2]]
                coord = [x, y, z]
                labelcoord = list(zip(*coord))
                unstable_entries[entry] = labelcoord[0]

        return lines, stable_entries, unstable_entries

    def show(self, label_stable=True, label_unstable=False, filename=""):
        """
        Draws the convex hull diagram using Matplotlib and show it.
        """
        plt = self._get_plot(label_stable=label_stable,
                             label_unstable=label_unstable)
        if filename == "":
            plt.show()
        else:
            plt.savefig(filename, bbox_inches=0)

    def _get_plot(self, label_stable=True, label_unstable=False):
        """
        Plot convex hull of Pourbaix Diagram entries
        """
        import matplotlib.pyplot as plt
        import mpl_toolkits.mplot3d.axes3d as p3
        from matplotlib.font_manager import FontProperties
        fig = plt.figure()
        ax = p3.Axes3D(fig)
        font = FontProperties()
        font.set_weight("bold")
        font.set_size(14)
        (lines, labels, unstable) = self.pourbaix_hull_plot_data
        count = 1
        newlabels = list()
        for x, y, z in lines:
            ax.plot(x, y, z, "bo-", linewidth=3, markeredgecolor="b",
                    markerfacecolor="r", markersize=10)
        for coords in sorted(labels.keys()):
            entry = labels[coords]
            label = self.print_name(entry)
            if label_stable:
                ax.text(coords[0], coords[1], coords[2], str(count))
                newlabels.append("{} : {}".format(
                    count, latexify_ion(latexify(label))))
                count += 1

        if label_unstable:
            for entry in unstable.keys():
                label = self.print_name(entry)
                coords = unstable[entry]
                ax.plot([coords[0], coords[0]], [coords[1], coords[1]],
                        [coords[2], coords[2]], "bo", markerfacecolor="g",
                        markersize=10)
                ax.text(coords[0], coords[1], coords[2], str(count))
                newlabels.append("{} : {}".format(
                    count, latexify_ion(latexify(label))))
                count += 1

        plt.figtext(0.01, 0.01, "\n".join(newlabels))
        plt.xlabel("pH")
        plt.ylabel("V")
        return plt

    def plot_planes(self):
        """
        Plot the free energy facets as a function of pH and V
        """
        if self.show_unstable:
            entries = self._pd._all_entries
        else:
            entries = self._pd.stable_entries
        num_plots = len(entries)
        import matplotlib.pyplot as plt
        colormap = plt.cm.gist_ncar
        fig = plt.figure().gca(projection='3d')
        color_array = [colormap(i) for i in np.linspace(0, 0.9, num_plots)]
        labels = []
        color_index = -1
        for entry in entries:
            normal = np.array([-PREFAC * entry.npH, -entry.nPhi, +1])
            d = entry.g0
            color_index += 1
            pH, V = np.meshgrid(np.linspace(-10, 28, 100),
                                np.linspace(-3, 3, 100))
            g = (-normal[0] * pH - normal[1] * V + d) / normal[2]
            lbl = latexify_ion(
                latexify(entry._entry.composition.reduced_formula))
            labels.append(lbl)
            fig.plot_surface(pH, V, g, color=color_array[color_index],
                             label=lbl)
        plt.legend(labels)
        plt.xlabel("pH")
        plt.ylabel("E (V)")
        plt.show()

    def plot_chempot_range_map(self, limits=None, title="", filename=""):
        self.plot_pourbaix(limits, title, filename)

    def plot_pourbaix(self, limits=None, title="", filename="", label_domains=True):
        plt = self.get_pourbaix_plot(limits=limits, title=title, label_domains=label_domains)
        if filename == "":
            plt.show()
        else:
            f = plt.gcf()
            f.set_size_inches((11.5, 9))
            plt.tight_layout(pad=1.09)

    def plot_entry_stability(self, entry, resolution=100, e_hull_max=1,
                             **kwargs):
        # plot the Pourbaix diagram
        plt = self.get_pourbaix_plot(**kwargs)
        ax = plt.gca()
        ph_min, ph_max, v_min, v_max = ax.get_xlim() + ax.get_ylim()

        # Get stability map
        phs, vs = np.mgrid[ph_min:ph_max:resolution*1j, 
                           v_min:v_max:resolution*1j]
        if self._pd._multielement and not isinstance(entry, MultiEntry):
            _, _, m_entries = self._analyzer.\
                    get_all_decomp_and_e_above_hull(entry)
            warnings.warn("{} is not a multi-entry, calculating stability of "
                          "representative multientries".format(entry.name))
            stability = [self._analyzer.get_entry_stability(m_entries[0], ph, v)
                         for ph, v in zip(phs.ravel(), vs.ravel())]
            for m_entry in m_entries[1:]:
                new_stab = [self._analyzer.get_entry_stability(m_entry, ph, v)
                            for ph, v in zip(phs.ravel(), vs.ravel())]
                stability = np.min(np.array([stability, new_stab]), axis=0)
        else:
            stability = [self._analyzer.get_entry_stability(entry, ph, v)
                         for ph, v in zip(phs.ravel(), vs.ravel())]
        stability = np.array(stability).reshape(phs.shape)

        # Plot stability map
        plt.pcolor(phs, vs, stability, cmap='RdYlBu_r',
                   vmin=0, vmax=e_hull_max)
        cbar = plt.colorbar()
        cbar.set_label("Stability of {} (eV)".format(self.print_name(entry)))

        # Set ticklabels
        ticklabels = [t.get_text() for t in cbar.ax.get_yticklabels()]
        ticklabels[-1] = '>={}'.format(ticklabels[-1])
        cbar.ax.set_yticklabels(ticklabels)

        return plt

    def pourbaix_plot_data(self, limits=None):
        """
        Get data required to plot Pourbaix diagram.

        Args:
            limits: 2D list containing limits of the Pourbaix diagram
                of the form [[xlo, xhi], [ylo, yhi]]

        Returns:
            stable_entries, unstable_entries
            stable_entries: dict of lines. The keys are Pourbaix Entries, and
            lines are in the form of a list
            unstable_entries: list of unstable entries
        """

        analyzer = PourbaixAnalyzer(self._pd)
        self._analyzer = analyzer
        if limits:
            analyzer.chempot_limits = limits
        chempot_ranges = analyzer.get_chempot_range_map(limits)
        self.chempot_ranges = chempot_ranges
        stable_entries_list = collections.defaultdict(list)

        for entry in chempot_ranges:
            for line in chempot_ranges[entry]:
                x = [line.coords[0][0], line.coords[1][0]]
                y = [line.coords[0][1], line.coords[1][1]]
                coords = [x, y]
                stable_entries_list[entry].append(coords)

        unstable_entries_list = [entry for entry in self._pd.all_entries
                                 if entry not in self._pd.stable_entries]

        return stable_entries_list, unstable_entries_list

    def get_center(self, lines):
        """
        Returns coordinates of center of a domain. Useful
        for labeling a Pourbaix plot.

        Args:
            lines:
                Lines corresponding to a domain
            limits:
                Limits of Pourbaix diagram

        Returns:
            center_x, center_y:
                x,y coordinate of center of domain. If domain lies
                outside limits, center will lie on the boundary.
        """
        center_x = 0.0
        center_y = 0.0
        coords = []
        count_center = 0.0
        for line in lines:
            for coord in np.array(line).T:
                if not in_coord_list(coords, coord):
                    coords.append(coord.tolist())
                    cx = coord[0]
                    cy = coord[1]
                    center_x += cx
                    center_y += cy
                    count_center += 1.0
        if count_center == 0.0:
            count_center = 1.0
        center_x /= count_center
        center_y /= count_center
        return center_x, center_y

    def get_distribution_corrected_center(self, lines, h2o_h_line=None, h2o_o_line=None, radius=None):
        """
        Returns coordinates of distribution corrected center of a domain. Similar to get_center(), but
        considers the distance to the surronding lines that mostly affects the feeling of "center".
        This function will also try avoid overalapping the text babel with H2O stability line if H2O
        stability line is provided. Useful for labeling a Pourbaix plot.

        Args:
            lines:
                Lines corresponding to a domain
            limits:
                Limits of Pourbaix diagram
            h2o_h_line: Hydrogen line of H2O stability
            h2o_o_line: Oxygen line of H2O stablity
            radius: Half height of the text label.

        Returns:
            center_x, center_y:
                x,y coordinate of center of domain. If domain lies
                outside limits, center will lie on the boundary.
        """
        coords = []
        pts_x = []
        pts_y = []
        for line in lines:
            for coord in np.array(line).T:
                if not in_coord_list(coords, coord):
                    coords.append(coord.tolist())
                    cx = coord[0]
                    cy = coord[1]
                    pts_x.append(cx)
                    pts_y.append(cy)
        if len(pts_x) < 1:
            return 0.0, 0.0
        cx_1 = (max(pts_x) + min(pts_x)) / 2.0
        cy_1 = (max(pts_y) + min(pts_y)) / 2.0
        mid_x_list = []
        mid_y_list = []
        # move the center to the center of surrounding lines
        for line in lines:
            (x1, y1), (x2, y2) = np.array(line).T
            if (x1 - cx_1) * (x2 - cx_1) <= 0.0:
                # horizontal line
                mid_y = ((y2 - y1) / (x2 - x1)) * (cx_1 - x1) + y1
                assert (y2 - mid_y) * (y1 - mid_y) <= 0.0
                mid_y_list.append(mid_y)
            if (y1 - cy_1) * (y2 - cy_1) <= 0.0:
                # vertical line
                mid_x = ((x2 - x1) / (y2 - y1)) * (cy_1 - y1) + x1
                assert (x2 - mid_x) * (x1 - mid_x) <= 0.0
                mid_x_list.append(mid_x)
        upper_y = sorted([y for y in mid_y_list if y >= cy_1] + [cy_1])[0]
        lower_y = sorted([y for y in mid_y_list if y < cy_1] + [cy_1])[-1]
        left_x = sorted([x for x in mid_x_list if x <= cx_1] + [cx_1])[-1]
        right_x = sorted([x for x in mid_x_list if x > cx_1] + [cx_1])[0]
        center_x = (left_x + right_x) / 2.0
        center_y = (upper_y + lower_y) / 2.0
        if h2o_h_line is not None:
            (h2o_h_x1, h2o_h_y1), (h2o_h_x2, h2o_h_y2) = h2o_h_line.T
            h_slope = (h2o_h_y2 - h2o_h_y1) / (h2o_h_x2 - h2o_h_x1)
            (h2o_o_x1, h2o_o_y1), (h2o_o_x2, h2o_o_y2) = h2o_o_line.T
            o_slope = (h2o_o_y2 - h2o_o_y1) / (h2o_o_x2 - h2o_o_x1)
            h_y = h_slope * (cx_1 - h2o_h_x1) + h2o_h_y1
            o_y = o_slope * (cx_1 - h2o_o_x1) + h2o_o_y1
            h2o_y = None
            if abs(center_y - h_y) < radius:
                h2o_y = h_y
            elif  abs(center_y - o_y) < radius:
                h2o_y = o_y
            if h2o_y is not None:
                if (upper_y - lower_y) / 2.0 > radius * 2.0:
                    # The space can hold the whole text (radius * 2.0)
                    if h2o_y > center_y:
                        center_y = h2o_y - radius
                    else:
                        center_y = h2o_y + radius
        return center_x, center_y

    def get_pourbaix_plot(self, limits=[[-2, 14], [-3, 3]],
                          title="", label_domains=True):
        """
        Plot Pourbaix diagram.

        Args:
            limits: 2D list containing limits of the Pourbaix diagram
                of the form [[xlo, xhi], [ylo, yhi]]

        Returns:
            plt:
                matplotlib plot object
        """
        plt = pretty_plot(16)
        (stable, unstable) = self.pourbaix_plot_data(limits)
        if limits:
            xlim = limits[0]
            ylim = limits[1]
        else:
            xlim = self._analyzer.chempot_limits[0]
            ylim = self._analyzer.chempot_limits[1]

        h_line = np.transpose([[xlim[0], -xlim[0] * PREFAC],
                               [xlim[1], -xlim[1] * PREFAC]])
        o_line = np.transpose([[xlim[0], -xlim[0] * PREFAC + 1.23],
                               [xlim[1], -xlim[1] * PREFAC + 1.23]])
        neutral_line = np.transpose([[7, ylim[0]], [7, ylim[1]]])
        V0_line = np.transpose([[xlim[0], 0], [xlim[1], 0]])

        ax = plt.gca()
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        lw = 3
        plt.plot(h_line[0], h_line[1], "r--", linewidth=lw)
        plt.plot(o_line[0], o_line[1], "r--", linewidth=lw)
        plt.plot(neutral_line[0], neutral_line[1], "k-.", linewidth=lw)
        plt.plot(V0_line[0], V0_line[1], "k-.", linewidth=lw)

        for entry, lines in stable.items():
            center_x = 0.0
            center_y = 0.0
            coords = []
            count_center = 0.0
            for line in lines:
                (x, y) = line
                plt.plot(x, y, "k-", linewidth=lw)
                for coord in np.array(line).T:
                    if not in_coord_list(coords, coord):
                        coords.append(coord.tolist())
                        cx = coord[0]
                        cy = coord[1]
                        center_x += cx
                        center_y += cy
                        count_center += 1.0
            if count_center == 0.0:
                count_center = 1.0
            center_x /= count_center
            center_y /= count_center
            if ((center_x <= xlim[0]) | (center_x >= xlim[1]) |
                    (center_y <= ylim[0]) | (center_y >= ylim[1])):
                continue
            xy = (center_x, center_y)
            if label_domains:
                plt.annotate(self.print_name(entry), xy, ha='center',
                             va='center', fontsize=20, color="b")

        plt.xlabel("pH")
        plt.ylabel("E (V)")
        plt.title(title, fontsize=20, fontweight='bold')
        return plt

    def print_name(self, entry):
        """
        Print entry name if single, else print multientry
        """
        str_name = ""
        if isinstance(entry, MultiEntry):
            if len(entry.entrylist) > 2:
                return str(self._pd.qhull_entries.index(entry))
            for e in entry.entrylist:
                str_name += latexify_ion(latexify(e.name)) + " + "
            str_name = str_name[:-3]
            return str_name
        else:
            return latexify_ion(latexify(entry.name))

    def legend(self, label_unstable=False, legend_file=""):
        if self._pd._multielement:
            unprocessed_entries = self._pd.unprocessed_entries
            set_of_entries = set()
            list_of_entries = {}
            for entry in self._pd.stable_entries:
                index_ent = self._pd.qhull_entries.index(entry)
                str_ename = ""
                for e in entry.entrylist:
                    str_ename += e.name + " + "
                    for ent in unprocessed_entries:
                        if ent.name == e.name:
                            indx = unprocessed_entries.index(ent)
                            set_of_entries.add(indx)
                            continue
                str_ename = str_ename[:-3]
                list_of_entries[index_ent] = str_ename
            if label_unstable:
                for entry in [entry for entry in self._pd.all_entries
                              if entry not in self._pd.stable_entries]:
                    for e in entry.entrylist:
                        indx = unprocessed_entries.index(e)
                        set_of_entries.add(indx)
            str_labels = " Species: \n"
            if legend_file:
                f = open(legend_file, 'w')
                for i in list_of_entries.keys():
                    str_labels += str(i) + " : " + list_of_entries[i] + "\n"
                f.write(str_labels)
                f.close()
            return str_labels

    def write_image(self, plt, stream, image_format="svg"):
        """
        Writes the phase diagram to an image in a stream.

        Args:
            plt:
                matplotlib plot
            stream:
                stream to write to. Can be a file stream or a StringIO stream.
            image_format
                format for image. Can be any of matplotlib supported formats.
                Defaults to svg for best results for vector graphics.
        """
        f = plt.gcf()
        f.set_size_inches((12, 10))
        plt.tight_layout(pad=1.09)
        plt.savefig(stream, format=image_format)

    def domain_vertices(self, entry):
        """
        Returns the vertices of the Pourbaix domain.

        Args:
            entry: Entry for which domain vertices are desired

        Returns:
            list of vertices
        """
        if entry not in self._analyzer.pourbaix_domain_vertices.keys():
            return []
        return self._analyzer.pourbaix_domain_vertices[entry]

    def get_pourbaix_plot_colorfill_by_element(self, limits=None, title="",
                                                label_domains=True, element=None):
        """
        Color domains by element
        """
        from matplotlib.patches import Polygon

        entry_dict_of_multientries = collections.defaultdict(list)
        plt = pretty_plot(16)
        optim_colors = ['#0000FF', '#FF0000', '#00FF00', '#FFFF00', '#FF00FF',
                         '#FF8080', '#DCDCDC', '#800000', '#FF8000']
        optim_font_color = ['#FFFFA0', '#00FFFF', '#FF00FF', '#0000FF', '#00FF00',
                            '#007F7F', '#232323', '#7FFFFF', '#007FFF']
        hatch = ['/', '\\', '|', '-', '+', 'o', '*']
        (stable, unstable) = self.pourbaix_plot_data(limits)
        num_of_overlaps = {key: 0 for key in stable.keys()}
        for entry in stable:
            if isinstance(entry, MultiEntry):
                for e in entry.entrylist:
                    if element in e.composition.elements:
                        entry_dict_of_multientries[e.name].append(entry)
                        num_of_overlaps[entry] += 1
            else:
                entry_dict_of_multientries[entry.name].append(entry)
        if limits:
            xlim = limits[0]
            ylim = limits[1]
        else:
            xlim = self._analyzer.chempot_limits[0]
            ylim = self._analyzer.chempot_limits[1]

        h_line = np.transpose([[xlim[0], -xlim[0] * PREFAC],
                               [xlim[1], -xlim[1] * PREFAC]])
        o_line = np.transpose([[xlim[0], -xlim[0] * PREFAC + 1.23],
                               [xlim[1], -xlim[1] * PREFAC + 1.23]])
        neutral_line = np.transpose([[7, ylim[0]], [7, ylim[1]]])
        V0_line = np.transpose([[xlim[0], 0], [xlim[1], 0]])

        ax = plt.gca()
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        from pymatgen import Composition, Element
        from pymatgen.core.ion import Ion

        def len_elts(entry):
            if "(s)" in entry:
                comp = Composition(entry[:-3])
            else:
                comp = Ion.from_formula(entry)
            return len([el for el in comp.elements if el not in
                        [Element("H"), Element("O")]])

        sorted_entry = entry_dict_of_multientries.keys()
        sorted_entry.sort(key=len_elts)
        i = -1
        label_chr = map(chr, list(range(65, 91)))
        for entry in sorted_entry:
            color_indx = 0
            x_coord = 0.0
            y_coord = 0.0
            npts = 0
            i += 1
            for e in entry_dict_of_multientries[entry]:
                hc = 0
                fc = 0
                bc = 0
                xy = self.domain_vertices(e)
                c = self.get_center(stable[e])
                x_coord += c[0]
                y_coord += c[1]
                npts += 1
                color_indx = i
                if "(s)" in entry:
                    comp = Composition(entry[:-3])
                else:
                    comp = Ion.from_formula(entry)
                if len([el for el in comp.elements if el not in
                         [Element("H"), Element("O")]]) == 1:
                    if color_indx >= len(optim_colors):
                        color_indx = color_indx -\
                         int(color_indx / len(optim_colors)) * len(optim_colors)
                    patch = Polygon(xy, facecolor=optim_colors[color_indx],
                                     closed=True, lw=3.0, fill=True)
                    bc = optim_colors[color_indx]
                else:
                    if color_indx >= len(hatch):
                        color_indx = color_indx - int(color_indx / len(hatch)) * len(hatch)
                    patch = Polygon(xy, hatch=hatch[color_indx], closed=True, lw=3.0, fill=False)
                    hc = hatch[color_indx]
                ax.add_patch(patch)

            xy_center = (x_coord / npts, y_coord / npts)
            if label_domains:
                if color_indx >= len(optim_colors):
                    color_indx = color_indx -\
                        int(color_indx / len(optim_colors)) * len(optim_colors)
                fc = optim_font_color[color_indx]
                if bc and not hc:
                    bbox = dict(boxstyle="round", fc=fc)
                if hc and not bc:
                    bc = 'k'
                    fc = 'w'
                    bbox = dict(boxstyle="round", hatch=hc, fill=False)
                if bc and hc:
                    bbox = dict(boxstyle="round", hatch=hc, fc=fc)
#                 bbox.set_path_effects([PathEffects.withSimplePatchShadow()])
                plt.annotate(latexify_ion(latexify(entry)), xy_center,
                            color=bc, fontsize=30, bbox=bbox)
#                 plt.annotate(label_chr[i], xy_center,
#                               color=bc, fontsize=30, bbox=bbox)

        lw = 3
        plt.plot(h_line[0], h_line[1], "r--", linewidth=lw)
        plt.plot(o_line[0], o_line[1], "r--", linewidth=lw)
        plt.plot(neutral_line[0], neutral_line[1], "k-.", linewidth=lw)
        plt.plot(V0_line[0], V0_line[1], "k-.", linewidth=lw)

        plt.xlabel("pH")
        plt.ylabel("E (V)")
        plt.title(title, fontsize=20, fontweight='bold')
        return plt

    def get_pourbaix_plot_colorfill_by_domain_name(self, limits=None, title="",
            label_domains=True, label_color='k', domain_color=None, domain_fontsize=None,
            domain_edge_lw=0.5, bold_domains=None, cluster_domains=(),
            add_h2o_stablity_line=True, add_center_line=False, h2o_lw=0.5,
            fill_domain=True, width=8, height=None, font_family='Times New Roman'):
        """
        Color domains by the colors specific by the domain_color dict

        Args:
            limits: 2D list containing limits of the Pourbaix diagram
                of the form [[xlo, xhi], [ylo, yhi]]
            lable_domains (Bool): whether add the text lable for domains
            label_color (str): color of domain lables, defaults to be black
            domain_color (dict): colors of each domain e.g {"Al(s)": "#FF1100"}. If set
                to None default color set will be used.
            domain_fontsize (int): Font size used in domain text labels.
            domain_edge_lw (int): line width for the boundaries between domains.
            bold_domains (list): List of domain names to use bold text style for domain
                lables. If set to False, no domain will be bold.
            cluster_domains (list): List of domain names in cluster phase
            add_h2o_stablity_line (Bool): whether plot H2O stability line
            add_center_line (Bool): whether plot lines shows the center coordinate
            h2o_lw (int): line width for H2O stability line and center lines
            fill_domain (bool): a version without color will be product if set 
                to False.
            width (float): Width of plot in inches. Defaults to 8in.
                height (float): Height of plot in inches. Defaults to width * golden
                ratio.
            font_family (str): Font family of the labels
        """

        def special_lines(xlim, ylim):
            h_line = np.transpose([[xlim[0], -xlim[0] * PREFAC],
                                   [xlim[1], -xlim[1] * PREFAC]])
            o_line = np.transpose([[xlim[0], -xlim[0] * PREFAC + 1.23],
                                   [xlim[1], -xlim[1] * PREFAC + 1.23]])
            neutral_line = np.transpose([[7, ylim[0]], [7, ylim[1]]])
            V0_line = np.transpose([[xlim[0], 0], [xlim[1], 0]])
            return h_line, o_line, neutral_line, V0_line

        from matplotlib.patches import Polygon
        import copy

        default_domain_font_size = 12
        default_solid_phase_color = '#b8f9e7'    # this slighly darker than the MP scheme, to
        default_cluster_phase_color = '#d0fbef'  # avoid making the cluster phase too light

        plt = pretty_plot(width=width, height=height, dpi=300)

        (stable, unstable) = self.pourbaix_plot_data(limits)

        xlim, ylim = limits[:2] if limits else self._analyzer.chempot_limits[:2]
        h_line, o_line, neutral_line, V0_line = special_lines(xlim, ylim)
        ax = plt.gca()
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.tick_params(direction='out')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        sorted_entry = list(stable.keys())
        sorted_entry.sort(key=lambda en: en.energy, reverse=True)

        if domain_fontsize is None:
            domain_fontsize = {en.name: default_domain_font_size for en in sorted_entry}
        elif not isinstance(domain_fontsize, dict):
            domain_fontsize = {en.name: domain_fontsize for en in sorted_entry}
        if domain_color is None:
            domain_color = {en.name: default_solid_phase_color if '(s)' in en.name else
                            (default_cluster_phase_color if en.name in cluster_domains else 'w')
                            for i, en in enumerate(sorted_entry)}
        else:
            domain_color = {en.name: domain_color[en.name] if en.name in domain_color else "w"
                            for i, en in enumerate(sorted_entry)}
        if bold_domains is None:
            bold_domains = [en.name for en in sorted_entry if '(s)' not in en.name]
        if bold_domains == False:
            bold_domains = []

        for entry in sorted_entry:
            xy = self.domain_vertices(entry)
            if add_h2o_stablity_line:
                c = self.get_distribution_corrected_center(stable[entry], h_line, o_line, 0.3)
            else:
                c = self.get_distribution_corrected_center(stable[entry])
            patch = Polygon(xy, facecolor=domain_color[entry.name], edgecolor="black",
                            closed=True, lw=domain_edge_lw, fill=fill_domain, antialiased=True)
            ax.add_patch(patch)
            if label_domains:
                if platform.system() == 'Darwin' and font_family == "Times New Roman":
                    # Have to hack to the hard coded font path to get current font On Mac OS X
                    if entry.name in bold_domains:
                        font = FontProperties(fname='/Library/Fonts/Times New Roman Bold.ttf',
                                              size=domain_fontsize[entry.name])
                    else:
                        font = FontProperties(fname='/Library/Fonts/Times New Roman.ttf',
                                              size=domain_fontsize[entry.name])
                else:
                    if entry.name in bold_domains:
                        font = FontProperties(family=font_family,
                                              weight='bold',
                                              size=domain_fontsize[entry.name])
                    else:
                        font = FontProperties(family=font_family,
                                              weight='regular',
                                              size=domain_fontsize[entry.name])
                plt.text(*c, s=self.print_name(entry),
                         fontproperties=font,
                         horizontalalignment="center", verticalalignment="center",
                         multialignment="center", color=label_color)

        if add_h2o_stablity_line:
            dashes = (3, 1.5)
            line, = plt.plot(h_line[0], h_line[1], "k--", linewidth=h2o_lw, antialiased=True)
            line.set_dashes(dashes)
            line, = plt.plot(o_line[0], o_line[1], "k--", linewidth=h2o_lw, antialiased=True)
            line.set_dashes(dashes)
        if add_center_line:
            plt.plot(neutral_line[0], neutral_line[1], "k-.", linewidth=h2o_lw, antialiased=False)
            plt.plot(V0_line[0], V0_line[1], "k-.", linewidth=h2o_lw, antialiased=False)

        plt.xlabel("pH", fontname=font_family, fontsize=18)
        plt.ylabel("E (V vs SHE)", fontname=font_family, fontsize=18)
        plt.xticks(fontname=font_family, fontsize=16)
        plt.yticks(fontname=font_family, fontsize=16)
        plt.title(title, fontsize=20, fontweight='bold', fontname=font_family)
        return plt

    def get_pourbaix_mark_passive(self, limits=None, title="", label_domains=True, passive_entry=None):
        """
        Color domains by element
        """
        from matplotlib.patches import Polygon
        from pymatgen import Element
        from itertools import chain
        import operator

        plt = pretty_plot(16)
        optim_colors = ['#0000FF', '#FF0000', '#00FF00', '#FFFF00', '#FF00FF',
                        '#FF8080', '#DCDCDC', '#800000', '#FF8000']
        optim_font_colors = ['#FFC000', '#00FFFF', '#FF00FF', '#0000FF', '#00FF00',
                            '#007F7F', '#232323', '#7FFFFF', '#007FFF']
        (stable, unstable) = self.pourbaix_plot_data(limits)
        mark_passive = {key: 0 for key in stable.keys()}

        if self._pd._elt_comp:
            maxval = max(six.iteritems(self._pd._elt_comp), key=operator.itemgetter(1))[1]
            key = [k for k, v in self._pd._elt_comp.items() if v == maxval]
        passive_entry = key[0]

        def list_elts(entry):
            elts_list = set()
            if isinstance(entry, MultiEntry):
                for el in chain.from_iterable([[el for el in e.composition.elements]
                                                for e in entry.entrylist]):
                    elts_list.add(el)
            else:
                elts_list = entry.composition.elements
            return elts_list

        for entry in stable:
            if passive_entry + str("(s)") in entry.name:
                mark_passive[entry] = 2
                continue
            if "(s)" not in entry.name:
                continue
            elif len(set([Element("O"), Element("H")]).intersection(set(list_elts(entry)))) > 0:
                mark_passive[entry] = 1

        if limits:
            xlim = limits[0]
            ylim = limits[1]
        else:
            xlim = self._analyzer.chempot_limits[0]
            ylim = self._analyzer.chempot_limits[1]

        h_line = np.transpose([[xlim[0], -xlim[0] * PREFAC],
                               [xlim[1], -xlim[1] * PREFAC]])
        o_line = np.transpose([[xlim[0], -xlim[0] * PREFAC + 1.23],
                               [xlim[1], -xlim[1] * PREFAC + 1.23]])
        neutral_line = np.transpose([[7, ylim[0]], [7, ylim[1]]])
        V0_line = np.transpose([[xlim[0], 0], [xlim[1], 0]])

        ax = plt.gca()
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        for e in stable.keys():
            xy = self.domain_vertices(e)
            c = self.get_center(stable[e])
            if mark_passive[e] == 1:
                color = optim_colors[0]
                fontcolor = optim_font_colors[0]
                colorfill = True
            elif mark_passive[e] == 2:
                color = optim_colors[1]
                fontcolor = optim_font_colors[1]
                colorfill = True
            else:
                color = "w"
                colorfill = False
                fontcolor = "k"
            patch = Polygon(xy, facecolor=color, closed=True, lw=3.0, fill=colorfill)
            ax.add_patch(patch)
            if label_domains:
                plt.annotate(self.print_name(e), c, color=fontcolor, fontsize=20)

        lw = 3
        plt.plot(h_line[0], h_line[1], "r--", linewidth=lw)
        plt.plot(o_line[0], o_line[1], "r--", linewidth=lw)
        plt.plot(neutral_line[0], neutral_line[1], "k-.", linewidth=lw)
        plt.plot(V0_line[0], V0_line[1], "k-.", linewidth=lw)

        plt.xlabel("pH")
        plt.ylabel("E (V)")
        plt.title(title, fontsize=20, fontweight='bold')
        return plt


def latexify_ion(formula):
    return re.sub(r"()\[([^)]*)\]", r"\1$^{\2}$", formula)
