#!/usr/bin/env python

"""
This module provides classes for plotting Pourbaix objects.
"""

from __future__ import division

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

from pymatgen.phasediagram.plotter import uniquelines
from pymatgen.util.string_utils import latexify
from pymatgen.util.plotting_utils import get_publication_quality_plot
from pymatgen.util.coord_utils import in_coord_list


class PourbaixPlotter(object):
    """
    A plotter class for phase diagrams.
    """

    def __init__(self, pourbaixdiagram, show_unstable=False):
        """
        Args:
            phasediagram:
                A PhaseDiagram object.
            show_unstable:
                Whether unstable phases will be plotted as well as red crosses.
                Defaults to False.
        """
        self._pd = pourbaixdiagram
        self.lines = uniquelines(self._pd.facets)
        self.show_unstable = show_unstable

    @property
    def pourbaix_hull_plot_data(self):
        """
        Pourbaix diagram convex hull data.

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
        for i in xrange(len(allentries)):
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
        plt = self._get_plot(label_stable, label_unstable)
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

        if self.show_unstable:
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
#        plt.tight_layout()
        return plt

    def plot_planes(self):
        """
        Plot the free energy planes as a function of pH and V
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

    def plot_pourbaix(self, limits=None, title="", filename=""):
        plt = self.get_pourbaix_plot(limits=limits, title=title)
        if filename == "":
            plt.show()
        else:
            plt.savefig(filename, bbox_inches=0)

    def pourbaix_plot_data(self, limits=None):
        """
        Get data required to plot Pourbaix diagram.

        Args:
            limits:
                2D list containing limits of the Pourbaix diagram
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

    def get_pourbaix_plot(self, limits=None, title=""):
        """
        Plot Pourbaix diagram.

        Args:
            limits:
                2D list containing limits of the Pourbaix diagram
                of the form [[xlo, xhi], [ylo, yhi]]

        Returns:
            plt:
                matplotlib plot object
        """
#        plt = get_publication_quality_plot(24, 14.4)
        plt = get_publication_quality_plot(16)
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
            plt.annotate(self.print_name(entry), xy, fontsize=20, color="b")

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
            entry:
                Entry for which domain vertices are desired

        Returns:
            list of vertices
        """
        if entry not in self._analyzer.pourbaix_domain_vertices.keys():
            return []
        return self._analyzer.pourbaix_domain_vertices[entry]


def latexify_ion(formula):
    return re.sub(r"()\[([^)]*)\]", r"\1$^{\2}$", formula)
