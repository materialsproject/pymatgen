# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function
import logging
from collections import OrderedDict

import numpy as np

from monty.json import jsanitize
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.util.plotting import pretty_plot
from pymatgen.electronic_structure.plotter import plot_brillouin_zone

"""
This module implements plotter for DOS and band structure.
"""

logger = logging.getLogger(__name__)


class PhononDosPlotter(object):
    """
    Class for plotting phonon DOSs. Note that the interface is extremely flexible
    given that there are many different ways in which people want to view
    DOS. The typical usage is::

        # Initializes plotter with some optional args. Defaults are usually
        # fine,
        plotter = PhononDosPlotter()

        # Adds a DOS with a label.
        plotter.add_dos("Total DOS", dos)

        # Alternatively, you can add a dict of DOSs. This is the typical
        # form returned by CompletePhononDos.get_element_dos().

    Args:
        stack: Whether to plot the DOS as a stacked area graph
        key_sort_func: function used to sort the dos_dict keys.
        sigma: A float specifying a standard deviation for Gaussian smearing
            the DOS for nicer looking plots. Defaults to None for no
            smearing.
    """

    def __init__(self, stack=False, sigma=None):
        self.stack = stack
        self.sigma = sigma
        self._doses = OrderedDict()

    def add_dos(self, label, dos):
        """
        Adds a dos for plotting.

        Args:
            label:
                label for the DOS. Must be unique.
            dos:
                PhononDos object
        """

        densities = dos.get_smeared_densities(self.sigma) if self.sigma \
            else dos.densities
        self._doses[label] = {'frequencies': dos.frequencies, 'densities': densities}

    def add_dos_dict(self, dos_dict, key_sort_func=None):
        """
        Add a dictionary of doses, with an optional sorting function for the
        keys.

        Args:
            dos_dict: dict of {label: Dos}
            key_sort_func: function used to sort the dos_dict keys.
        """
        if key_sort_func:
            keys = sorted(dos_dict.keys(), key=key_sort_func)
        else:
            keys = dos_dict.keys()
        for label in keys:
            self.add_dos(label, dos_dict[label])

    def get_dos_dict(self):
        """
        Returns the added doses as a json-serializable dict. Note that if you
        have specified smearing for the DOS plot, the densities returned will
        be the smeared densities, not the original densities.

        Returns:
            Dict of dos data. Generally of the form, {label: {'frequencies':..,
            'densities': ...}}
        """
        return jsanitize(self._doses)

    def get_plot(self, xlim=None, ylim=None):
        """
        Get a matplotlib plot showing the DOS.

        Args:
            xlim: Specifies the x-axis limits. Set to None for automatic
                determination.
            ylim: Specifies the y-axis limits.
        """
        import prettyplotlib as ppl
        from prettyplotlib import brewer2mpl

        ncolors = max(3, len(self._doses))
        ncolors = min(9, ncolors)
        colors = brewer2mpl.get_map('Set1', 'qualitative', ncolors).mpl_colors

        y = None
        alldensities = []
        allfrequencies = []
        plt = pretty_plot(12, 8)

        # Note that this complicated processing of frequencies is to allow for
        # stacked plots in matplotlib.
        for key, dos in self._doses.items():
            frequencies = dos['frequencies']
            densities = dos['densities']
            if y is None:
                y = np.zeros(frequencies.shape)
            if self.stack:
                y += densities
                newdens = y.copy()
            else:
                newdens = densities
            allfrequencies.append(frequencies)
            alldensities.append(newdens)

        keys = list(self._doses.keys())
        keys.reverse()
        alldensities.reverse()
        allfrequencies.reverse()
        allpts = []
        for i, (key, frequencies, densities) in enumerate(zip(keys, allfrequencies, alldensities)):
            allpts.extend(list(zip(frequencies, densities)))
            if self.stack:
                plt.fill(frequencies, densities, color=colors[i % ncolors],
                         label=str(key))
            else:
                ppl.plot(frequencies, densities, color=colors[i % ncolors],
                         label=str(key), linewidth=3)

        if xlim:
            plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)
        else:
            xlim = plt.xlim()
            relevanty = [p[1] for p in allpts
                         if xlim[0] < p[0] < xlim[1]]
            plt.ylim((min(relevanty), max(relevanty)))

        ylim = plt.ylim()
        plt.plot([0, 0], ylim, 'k--', linewidth=2)

        plt.xlabel('Frequencies (THz)')
        plt.ylabel('Density of states')

        plt.legend()
        leg = plt.gca().get_legend()
        ltext = leg.get_texts()  # all the text.Text instance in the legend
        plt.setp(ltext, fontsize=30)
        plt.tight_layout()
        return plt

    def save_plot(self, filename, img_format="eps", xlim=None, ylim=None):
        """
        Save matplotlib plot to a file.

        Args:
            filename: Filename to write to.
            img_format: Image format to use. Defaults to EPS.
            xlim: Specifies the x-axis limits. Set to None for automatic
                determination.
            ylim: Specifies the y-axis limits.
        """
        plt = self.get_plot(xlim, ylim)
        plt.savefig(filename, format=img_format)

    def show(self, xlim=None, ylim=None):
        """
        Show the plot using matplotlib.

        Args:
            xlim: Specifies the x-axis limits. Set to None for automatic
                determination.
            ylim: Specifies the y-axis limits.
        """
        plt = self.get_plot(xlim, ylim)
        plt.show()


class PhononBSPlotter(object):
    """
    Class to plot or get data to facilitate the plot of band structure objects.

    Args:
        bs: A BandStructureSymmLine object.
    """

    def __init__(self, bs):
        if not isinstance(bs, PhononBandStructureSymmLine):
            raise ValueError(
                "PhononBSPlotter only works with PhononBandStructureSymmLine objects. "
                "A PhononBandStructure object (on a uniform grid for instance and "
                "not along symmetry lines won't work)")
        self._bs = bs
        self._nb_bands = self._bs.nb_bands

    def _maketicks(self, plt):
        """
        utility private method to add ticks to a band structure
        """
        ticks = self.get_ticks()
        # Sanitize only plot the uniq values
        uniq_d = []
        uniq_l = []
        temp_ticks = list(zip(ticks['distance'], ticks['label']))
        for i in range(len(temp_ticks)):
            if i == 0:
                uniq_d.append(temp_ticks[i][0])
                uniq_l.append(temp_ticks[i][1])
                logger.debug("Adding label {l} at {d}".format(
                    l=temp_ticks[i][0], d=temp_ticks[i][1]))
            else:
                if temp_ticks[i][1] == temp_ticks[i - 1][1]:
                    logger.debug("Skipping label {i}".format(
                        i=temp_ticks[i][1]))
                else:
                    logger.debug("Adding label {l} at {d}".format(
                        l=temp_ticks[i][0], d=temp_ticks[i][1]))
                    uniq_d.append(temp_ticks[i][0])
                    uniq_l.append(temp_ticks[i][1])

        logger.debug("Unique labels are %s" % list(zip(uniq_d, uniq_l)))
        plt.gca().set_xticks(uniq_d)
        plt.gca().set_xticklabels(uniq_l)

        for i in range(len(ticks['label'])):
            if ticks['label'][i] is not None:
                # don't print the same label twice
                if i != 0:
                    if ticks['label'][i] == ticks['label'][i - 1]:
                        logger.debug("already print label... "
                                     "skipping label {i}".format(
                            i=ticks['label'][i]))
                    else:
                        logger.debug("Adding a line at {d}"
                                     " for label {l}".format(
                            d=ticks['distance'][i], l=ticks['label'][i]))
                        plt.axvline(ticks['distance'][i], color='k')
                else:
                    logger.debug("Adding a line at {d} for label {l}".format(
                        d=ticks['distance'][i], l=ticks['label'][i]))
                    plt.axvline(ticks['distance'][i], color='k')
        return plt

    def bs_plot_data(self):

        """
        Get the data nicely formatted for a plot

        Returns:
            A dict of the following format:
            ticks: A dict with the 'distances' at which there is a qpoint (the
            x axis) and the labels (None if no label)
            frequencies: A list (one element for each branch) of frequencies for
            each qpoint: [branch][qpoint][mode]. The data is
            stored by branch to facilitate the plotting
            lattice: The reciprocal lattice.
        """
        distance = []
        frequency = []

        ticks = self.get_ticks()

        for b in self._bs.branches:

            frequency.append([])
            distance.append([self._bs.distance[j]
                             for j in range(b['start_index'],
                                            b['end_index'] + 1)])

            for i in range(self._nb_bands):
                frequency[-1].append(
                    [self._bs.bands[i][j]
                     for j in range(b['start_index'], b['end_index'] + 1)])

        return {'ticks': ticks, 'distances': distance, 'frequency': frequency,
                'lattice': self._bs.lattice_rec.as_dict()}

    def get_plot(self, ylim=None):
        """
        Get a matplotlib object for the bandstructure plot.

        Args:
            ylim: Specify the y-axis (frequency) limits; by default None let
                the code choose.
        """
        plt = pretty_plot(12, 8)
        from matplotlib import rc
        import scipy.interpolate as scint
        try:
            rc('text', usetex=True)
        except:
            # Fall back on non Tex if errored.
            rc('text', usetex=False)

        band_linewidth = 1

        data = self.bs_plot_data()
        for d in range(len(data['distances'])):
            for i in range(self._nb_bands):
                plt.plot(data['distances'][d],
                         [data['frequency'][d][i][j]
                          for j in range(len(data['distances'][d]))], 'b-',
                         linewidth=band_linewidth)


        self._maketicks(plt)

        # plot y=0 line
        plt.axhline(0, linewidth=1, color='k')

        # Main X and Y Labels
        plt.xlabel(r'$\mathrm{Wave\ Vector}$', fontsize=30)
        ylabel = r'$\mathrm{Frequency\ (THz)}$'
        plt.ylabel(ylabel, fontsize=30)

        # X range (K)
        # last distance point
        x_max = data['distances'][-1][-1]
        plt.xlim(0, x_max)

        if ylim is not None:
            plt.ylim(ylim)

        plt.tight_layout()

        return plt

    def show(self, ylim=None):
        """
        Show the plot using matplotlib.

        Args:
            ylim: Specify the y-axis (frequency) limits; by default None let
                the code choose.
        """
        plt = self.get_plot(ylim)
        plt.show()

    def save_plot(self, filename, img_format="eps", ylim=None):
        """
        Save matplotlib plot to a file.

        Args:
            filename: Filename to write to.
            img_format: Image format to use. Defaults to EPS.
            ylim: Specifies the y-axis limits.
        """
        plt = self.get_plot(ylim=ylim)
        plt.savefig(filename, format=img_format)
        plt.close()

    def get_ticks(self):
        """
        Get all ticks and labels for a band structure plot.

        Returns:
            A dict with 'distance': a list of distance at which ticks should
            be set and 'label': a list of label for each of those ticks.
        """
        tick_distance = []
        tick_labels = []
        previous_label = self._bs.qpoints[0].label
        previous_branch = self._bs.branches[0]['name']
        for i, c in enumerate(self._bs.qpoints):
            if c.label is not None:
                tick_distance.append(self._bs.distance[i])
                this_branch = None
                for b in self._bs.branches:
                    if b['start_index'] <= i <= b['end_index']:
                        this_branch = b['name']
                        break
                if c.label != previous_label \
                        and previous_branch != this_branch:
                    label1 = c.label
                    if label1.startswith("\\") or label1.find("_") != -1:
                        label1 = "$" + label1 + "$"
                    label0 = previous_label
                    if label0.startswith("\\") or label0.find("_") != -1:
                        label0 = "$" + label0 + "$"
                    tick_labels.pop()
                    tick_distance.pop()
                    tick_labels.append(label0 + "$\mid$" + label1)
                else:
                    if c.label.startswith("\\") or c.label.find("_") != -1:
                        tick_labels.append("$" + c.label + "$")
                    else:
                        tick_labels.append(c.label)
                previous_label = c.label
                previous_branch = this_branch
        return {'distance': tick_distance, 'label': tick_labels}

    def plot_compare(self, other_plotter):
        """
        plot two band structure for comparison. One is in red the other in blue.
        The two band structures need to be defined on the same symmetry lines!
        and the distance between symmetry lines is
        the one of the band structure used to build the PhononBSPlotter

        Args:
            another PhononBSPlotter object defined along the same symmetry lines

        Returns:
            a matplotlib object with both band structures

        """

        data_orig = self.bs_plot_data()
        data = other_plotter.bs_plot_data()

        if len(data_orig['distances']) != len(data['distances']):
            raise ValueError('The two objects are not compatible.')

        plt = self.get_plot()
        band_linewidth = 1
        for i in range(other_plotter._nb_bands):
            for d in range(len(data_orig['distances'])):
                plt.plot(data_orig['distances'][d],
                         [e[i] for e in data['frequency']][d],
                         'r-', linewidth=band_linewidth)

        return plt

    def plot_brillouin(self):
        """
        plot the Brillouin zone
        """

        # get labels and lines
        labels = {}
        for q in self._bs.qpoints:
            if q.label:
                labels[q.label] = q.frac_coords

        lines = []
        for b in self._bs.branches:
            lines.append([self._bs.qpoints[b['start_index']].frac_coords,
                          self._bs.qpoints[b['end_index']].frac_coords])

        plot_brillouin_zone(self._bs.lattice_rec, lines=lines, labels=labels)

