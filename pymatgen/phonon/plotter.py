# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function
import logging
from collections import OrderedDict, namedtuple

import numpy as np
import scipy.constants as const

from monty.json import jsanitize
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.util.plotting import pretty_plot, add_fig_kwargs, get_ax_fig_plt
from pymatgen.electronic_structure.plotter import plot_brillouin_zone

"""
This module implements plotter for DOS and band structure.
"""

logger = logging.getLogger(__name__)

FreqUnits = namedtuple("FreqUnits", ["factor", "label"])


def freq_units(units):
    """
    Returns conversion factor from THz to the requred units and the label in the form of a namedtuple
    Accepted values: thz, ev, mev, ha, cm-1, cm^-1
    """

    d = {"thz": FreqUnits(1, "THz"),
         "ev": FreqUnits(const.value("hertz-electron volt relationship") * const.tera, "eV"),
         "mev": FreqUnits(const.value("hertz-electron volt relationship") * const.tera / const.milli, "meV"),
         "ha": FreqUnits(const.value("hertz-hartree relationship") * const.tera, "Ha"),
         "cm-1": FreqUnits(const.value("hertz-inverse meter relationship") * const.tera * const.centi, "cm$^{-1}$"),
         'cm^-1': FreqUnits(const.value("hertz-inverse meter relationship") * const.tera * const.centi, "cm$^{-1}$")
         }
    try:
        return d[units.lower().strip()]
    except KeyError:
        raise KeyError('Value for units `{}` unknown\nPossible values are:\n {}'.format(units, list(d.keys())))


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

    def get_plot(self, xlim=None, ylim=None, units="thz"):
        """
        Get a matplotlib plot showing the DOS.

        Args:
            xlim: Specifies the x-axis limits. Set to None for automatic
                determination.
            ylim: Specifies the y-axis limits.
            units: units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.
        """

        u = freq_units(units)

        ncolors = max(3, len(self._doses))
        ncolors = min(9, ncolors)

        import palettable

        colors = palettable.colorbrewer.qualitative.Set1_9.mpl_colors

        y = None
        alldensities = []
        allfrequencies = []
        plt = pretty_plot(12, 8)

        # Note that this complicated processing of frequencies is to allow for
        # stacked plots in matplotlib.
        for key, dos in self._doses.items():
            frequencies = dos['frequencies'] * u.factor
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

        plt.xlabel('Frequencies ({})'.format(u.label))
        plt.ylabel('Density of states')

        plt.legend()
        leg = plt.gca().get_legend()
        ltext = leg.get_texts()  # all the text.Text instance in the legend
        plt.setp(ltext, fontsize=30)
        plt.tight_layout()
        return plt

    def save_plot(self, filename, img_format="eps", xlim=None, ylim=None, units="thz"):
        """
        Save matplotlib plot to a file.

        Args:
            filename: Filename to write to.
            img_format: Image format to use. Defaults to EPS.
            xlim: Specifies the x-axis limits. Set to None for automatic
                determination.
            ylim: Specifies the y-axis limits.
            units: units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1
        """
        plt = self.get_plot(xlim, ylim, units=units)
        plt.savefig(filename, format=img_format)

    def show(self, xlim=None, ylim=None, units="thz"):
        """
        Show the plot using matplotlib.

        Args:
            xlim: Specifies the x-axis limits. Set to None for automatic
                determination.
            ylim: Specifies the y-axis limits.
            units: units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.
        """
        plt = self.get_plot(xlim, ylim, units=units)
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

    def get_plot(self, ylim=None, units="thz"):
        """
        Get a matplotlib object for the bandstructure plot.

        Args:
            ylim: Specify the y-axis (frequency) limits; by default None let
                the code choose.
            units: units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.
        """

        u = freq_units(units)

        plt = pretty_plot(12, 8)

        band_linewidth = 1

        data = self.bs_plot_data()
        for d in range(len(data['distances'])):
            for i in range(self._nb_bands):
                plt.plot(data['distances'][d],
                         [data['frequency'][d][i][j] * u.factor
                          for j in range(len(data['distances'][d]))], 'b-',
                         linewidth=band_linewidth)

        self._maketicks(plt)

        # plot y=0 line
        plt.axhline(0, linewidth=1, color='k')

        # Main X and Y Labels
        plt.xlabel(r'$\mathrm{Wave\ Vector}$', fontsize=30)
        ylabel = r'Frequency ({})'.format(u.label)
        plt.ylabel(ylabel, fontsize=30)

        # X range (K)
        # last distance point
        x_max = data['distances'][-1][-1]
        plt.xlim(0, x_max)

        if ylim is not None:
            plt.ylim(ylim)

        plt.tight_layout()

        return plt

    def show(self, ylim=None, units="thz"):
        """
        Show the plot using matplotlib.

        Args:
            ylim: Specify the y-axis (frequency) limits; by default None let
                the code choose.
            units: units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.
        """
        plt = self.get_plot(ylim, units=units)
        plt.show()

    def save_plot(self, filename, img_format="eps", ylim=None, units="thz"):
        """
        Save matplotlib plot to a file.

        Args:
            filename: Filename to write to.
            img_format: Image format to use. Defaults to EPS.
            ylim: Specifies the y-axis limits.
            units: units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.
        """
        plt = self.get_plot(ylim=ylim, units=units)
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
                    tick_labels.append(label0 + "$\\mid$" + label1)
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


class ThermoPlotter(object):
    """
    Plotter for thermodynamic properties obtained from phonon DOS.
    If the structure corresponding to the DOS, it will be used to extract the forumla unit and provide
    the plots in units of mol instead of mole-cell
    """

    def __init__(self, dos, structure=None):
        """
        Args:
            dos: A PhononDos object.
            structure: A Structure object corresponding to the structure used for the calculation.
        """
        self.dos = dos
        self.structure = structure

    def _plot_thermo(self, func, temperatures, factor=1, ax=None, ylabel=None, label=None, ylim=None, **kwargs):
        """
        Plots a thermodynamic property for a generic function from a PhononDos instance.

        Args:
            func: the thermodynamic function to be used to calculate the property
            temperatures: a list of temperatures
            factor: a multiplicative factor applied to the thermodynamic property calculated. Used to change
                the units.
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            ylabel: label for the y axis
            label: label of the plot
            ylim: tuple specifying the y-axis limits.
            kwargs: kwargs passed to the matplotlib function 'plot'.
        Returns:
            matplotlib figure
        """

        ax, fig, plt = get_ax_fig_plt(ax)

        values = []

        for t in temperatures:
            values.append(func(t, structure=self.structure) * factor)

        ax.plot(temperatures, values, label=label, **kwargs)

        if ylim:
            ax.set_ylim(ylim)

        ax.set_xlim((np.min(temperatures), np.max(temperatures)))
        ylim = plt.ylim()
        if ylim[0] < 0 < ylim[1]:
            plt.plot(plt.xlim(), [0, 0], 'k-', linewidth=1)

        ax.set_xlabel(r"$T$ (K)")
        if ylabel:
            ax.set_ylabel(ylabel)

        return fig

    @add_fig_kwargs
    def plot_cv(self, tmin, tmax, ntemp, ylim=None, **kwargs):
        """
        Plots the constant volume specific heat C_v in a temperature range.

        Args:
            tmin: minimum temperature
            tmax: maximum temperature
            ntemp: number of steps
            ylim: tuple specifying the y-axis limits.
            kwargs: kwargs passed to the matplotlib function 'plot'.
        Returns:
            matplotlib figure
        """
        temperatures = np.linspace(tmin, tmax, ntemp)

        if self.structure:
            ylabel = r"$C_v$ (J/K/mol)"
        else:
            ylabel = r"$C_v$ (J/K/mol-c)"

        fig = self._plot_thermo(self.dos.cv, temperatures, ylabel=ylabel, ylim=ylim, **kwargs)

        return fig

    @add_fig_kwargs
    def plot_entropy(self, tmin, tmax, ntemp, ylim=None, **kwargs):
        """
        Plots the vibrational entrpy in a temperature range.

        Args:
            tmin: minimum temperature
            tmax: maximum temperature
            ntemp: number of steps
            ylim: tuple specifying the y-axis limits.
            kwargs: kwargs passed to the matplotlib function 'plot'.
        Returns:
            matplotlib figure
        """
        temperatures = np.linspace(tmin, tmax, ntemp)

        if self.structure:
            ylabel = r"$S$ (J/K/mol)"
        else:
            ylabel = r"$S$ (J/K/mol-c)"

        fig = self._plot_thermo(self.dos.entropy, temperatures, ylabel=ylabel, ylim=ylim, **kwargs)

        return fig

    @add_fig_kwargs
    def plot_internal_energy(self, tmin, tmax, ntemp, ylim=None, **kwargs):
        """
        Plots the vibrational internal energy in a temperature range.

        Args:
            tmin: minimum temperature
            tmax: maximum temperature
            ntemp: number of steps
            ylim: tuple specifying the y-axis limits.
            kwargs: kwargs passed to the matplotlib function 'plot'.
        Returns:
            matplotlib figure
        """
        temperatures = np.linspace(tmin, tmax, ntemp)

        if self.structure:
            ylabel = r"$\Delta E$ (kJ/K/mol)"
        else:
            ylabel = r"$\Delta E$ (kJ/K/mol-c)"

        fig = self._plot_thermo(self.dos.internal_energy, temperatures, ylabel=ylabel, ylim=ylim,
                                factor=1e-3, **kwargs)

        return fig

    @add_fig_kwargs
    def plot_helmholtz_free_energy(self, tmin, tmax, ntemp, ylim=None, **kwargs):
        """
        Plots the vibrational contribution to the Helmoltz free energy in a temperature range.

        Args:
            tmin: minimum temperature
            tmax: maximum temperature
            ntemp: number of steps
            ylim: tuple specifying the y-axis limits.
            kwargs: kwargs passed to the matplotlib function 'plot'.
        Returns:
            matplotlib figure
        """
        temperatures = np.linspace(tmin, tmax, ntemp)

        if self.structure:
            ylabel = r"$\Delta F$ (kJ/K/mol)"
        else:
            ylabel = r"$\Delta F$ (kJ/K/mol-c)"

        fig = self._plot_thermo(self.dos.helmholtz_free_energy, temperatures, ylabel=ylabel, ylim=ylim,
                                factor=1e-3, **kwargs)

        return fig

    @add_fig_kwargs
    def plot_thermodynamic_properties(self, tmin, tmax, ntemp, ylim=None, **kwargs):
        """
        Plots all the thermodynamic properties in a temperature range.

        Args:
            tmin: minimum temperature
            tmax: maximum temperature
            ntemp: number of steps
            ylim: tuple specifying the y-axis limits.
            kwargs: kwargs passed to the matplotlib function 'plot'.
        Returns:
            matplotlib figure
        """
        temperatures = np.linspace(tmin, tmax, ntemp)

        mol = "" if self.structure else "-c"

        fig = self._plot_thermo(self.dos.cv, temperatures, ylabel="Thermodynamic properties", ylim=ylim,
                                label=r"$C_v$ (J/K/mol{})".format(mol), **kwargs)
        self._plot_thermo(self.dos.entropy, temperatures, ylim=ylim, ax=fig.axes[0],
                          label=r"$S$ (J/K/mol{})".format(mol), **kwargs)
        self._plot_thermo(self.dos.internal_energy, temperatures, ylim=ylim, ax=fig.axes[0], factor=1e-3,
                          label=r"$\Delta E$ (kJ/K/mol{})".format(mol), **kwargs)
        self._plot_thermo(self.dos.helmholtz_free_energy, temperatures, ylim=ylim, ax=fig.axes[0], factor=1e-3,
                          label=r"$\Delta F$ (kJ/K/mol{})".format(mol), **kwargs)

        fig.axes[0].legend(loc="best")

        return fig
