#!/usr/bin/env python

"""
This module provides classes for plotting PhaseDiagram objects.
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
        for i in xrange(0, len(allentries)):
            entry = allentries[i]
            if entry not in stable:
                x = [alldata[i][0], alldata[i][0]]
                y = [alldata[i][1], alldata[i][1]]
                z = [alldata[i][2], alldata[i][2]]
                coord = [x, y, z]
                labelcoord = list(zip(*coord))
                unstable_entries[entry] = labelcoord[0]

        return (lines, stable_entries, unstable_entries)

    def show(self, label_stable=True, label_unstable=False):
        """
        Draws the convex diagram using Matplotlib and show it.
        """
        plt = self._get_plot(label_stable, label_unstable)
        plt.show()

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
        font.set_size(20)
        (lines, labels, unstable) = self.pd_plot_data
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
                newlabels.append("{} : {}".format(count, latexify_ion(latexify(label))))
                count += 1

        if self.show_unstable:
            for entry in unstable.keys():
                label = self.print_name(entry)
                coords = unstable[entry]
                ax.plot([coords[0], coords[0]], [coords[1], coords[1]],\
                         [coords[2], coords[2]], "bo", markerfacecolor="g",\
                          markersize = 10)
                ax.text(coords[0], coords[1], coords[2], str(count))
                newlabels.append("{} : {}".format(count, latexify_ion(latexify(label))))
                count += 1

        plt.figtext(0.01, 0.01, "\n".join(newlabels))
        plt.xlabel("npH")
        plt.ylabel("nel")
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
            pH, V = np.meshgrid(np.linspace(-10, 28, 100), np.linspace(-3, 3, 100))
            g = (-normal[0] * pH - normal[1] * V + d) / normal[2]
            lbl = latexify_ion(latexify(entry._entry.composition.reduced_formula))
            labels.append(lbl)
            fig.plot_surface(pH, V, g, color=color_array[color_index], label=lbl)
        plt.legend(labels)
        plt.xlabel("pH")
        plt.ylabel("E (V)")
        plt.show()

    def plot_chempot_range_map(self, limits=None):
        """
        Plot pourbaix diagram
        """
        elements = ["H+", 'V']
        plt = get_publication_quality_plot(12, 8)
        analyzer = PourbaixAnalyzer(self._pd)
        chempot_ranges = analyzer.get_chempot_range_map()
        if (limits):
            xlim = limits[0]
            ylim = limits[1]
        else:
            xlim = analyzer.chempot_limits[0]
            ylim = analyzer.chempot_limits[1]

        h_line = np.transpose([[xlim[0], -xlim[0] * PREFAC],\
                                [xlim[1], -xlim[1] * PREFAC]])
        o_line = np.transpose([[xlim[0], -xlim[0] * PREFAC + 1.23],\
                                [xlim[1], -xlim[1] * PREFAC + 1.23]])
        neutral_line = np.transpose([[7, ylim[0]], [7, ylim[1]]])
        V0_line = np.transpose([[xlim[0], 0], [xlim[1], 0]])

        ax = plt.gca()
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        plt.plot(h_line[0], h_line[1], "r--")
        plt.plot(o_line[0], o_line[1], "r--")
        plt.plot(neutral_line[0], neutral_line[1], "k-.")
        plt.plot(V0_line[0], V0_line[1], "k-.")
        for entry, lines in chempot_ranges.items():
            region = []
            center_x = 0.0
            center_y = 0.0
            coords = []
            count_center = 0.0
            for line in lines:
                (x, y) = line.coords.transpose()
                plt.plot(x, y, "k-")
                for coord in line.coords:
                    if not in_coord_list(coords, coord):
                        coords.append(coord.tolist())
                        cx = coord[0]
                        cy = coord[1]
                        if cx < xlim[0]:
                            cx = xlim[0]
                        if cx > xlim[1]:
                            cx = xlim[1]
                        if cy < ylim[0]:
                            cy = ylim[0]
                        if cy > ylim[1]:
                            cy = ylim[1]
                        center_x += cx
                        center_y += cy
                        count_center += 1.0
                region.append(line.coords)
            if count_center == 0.0:
                count_center = 1.0
            center_x /= count_center
            center_y /= count_center
            if (((center_x == xlim[0]) | (center_x == xlim[1])) |\
                 ((center_y == ylim[0]) | (center_y == ylim[1]))):
                continue
            xy = (center_x, center_y)
            plt.annotate(self.print_name(entry), xy, fontsize=22)

        plt.xlabel("pH")
        plt.ylabel("E (V)")
        plt.tight_layout()
        plt.show()

    def print_name(self, entry):
        """
        Print entry name if single, else print multientry
        """
        str_name = ""
        if isinstance(entry, MultiEntry):
            for e in entry.entrylist:
                indx = self._pd.unprocessed_entries.index(e)
                str_name += str(indx) + " + "
            str_name = str_name[:-3]
            return str_name
        else:
            return latexify_ion(latexify(entry.name))

    def legend(self, label_unstable = False):
        import matplotlib.pyplot as plt
        if self._pd._multielement:
            fig = plt.figure(facecolor='white')
            fig.suptitle('Legend', fontsize=14, fontweight='bold')
            ax1 = plt.axes(frameon=False)
            ax1.get_xaxis().set_visible(False)
            ax1.axes.get_yaxis().set_visible(False)
            unprocessed_entries = self._pd.unprocessed_entries
            set_of_entries = set()
            for entry in self._pd.qhull_entries:
                for e in entry.entrylist:
                    for ent in unprocessed_entries: 
                        if ent.name == e.name:
                            indx = unprocessed_entries.index(ent)
                            set_of_entries.add(indx)
                            continue
            if (label_unstable):
                for entry in [entry for entry in self._pd.all_entries if entry not in self._pd.stable_entries]:
                    for e in entry.entrylist:
                        indx = unprocessed_entries.index(e)
                        set_of_entries.add(indx)
            str_labels = " Species: \n"
            f = open("Legends_file", 'w')
            for i in set_of_entries:
                str_labels += str(i) + " : " + str(unprocessed_entries[i].name) + "\n"
            f.write(str_labels)
            f.close()
            plt.text(0, 0, str_labels, fontsize=15)
            plt.show()


def latexify_ion(formula):
    return re.sub(r"()\[([^)]*)\]", r"\1$^{\2}$", formula)
