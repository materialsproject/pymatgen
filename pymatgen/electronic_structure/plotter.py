# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import logging
import math
import itertools
import warnings
from collections import OrderedDict

import numpy as np

from monty.json import jsanitize

from pymatgen.core.periodic_table import Element
from pymatgen.electronic_structure.core import Spin, Orbital, OrbitalType
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.util.plotting import pretty_plot, \
    add_fig_kwargs, get_ax3d_fig_plt
from collections import Counter
import copy

from pymatgen.electronic_structure.boltztrap import BoltztrapError
from pymatgen.symmetry.bandstructure import HighSymmKpath

"""
This module implements plotter for DOS and band structure.
"""

__author__ = "Shyue Ping Ong, Geoffroy Hautier, Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "May 1, 2012"

logger = logging.getLogger(__name__)


class DosPlotter:
    """
    Class for plotting DOSs. Note that the interface is extremely flexible
    given that there are many different ways in which people want to view
    DOS. The typical usage is::

        # Initializes plotter with some optional args. Defaults are usually
        # fine,
        plotter = DosPlotter()

        # Adds a DOS with a label.
        plotter.add_dos("Total DOS", dos)

        # Alternatively, you can add a dict of DOSs. This is the typical
        # form returned by CompleteDos.get_spd/element/others_dos().
        plotter.add_dos_dict({"dos1": dos1, "dos2": dos2})
        plotter.add_dos_dict(complete_dos.get_spd_dos())

    Args:
        zero_at_efermi: Whether to shift all Dos to have zero energy at the
            fermi energy. Defaults to True.
        stack: Whether to plot the DOS as a stacked area graph
        key_sort_func: function used to sort the dos_dict keys.
        sigma: A float specifying a standard deviation for Gaussian smearing
            the DOS for nicer looking plots. Defaults to None for no
            smearing.
    """

    def __init__(self, zero_at_efermi=True, stack=False, sigma=None):
        self.zero_at_efermi = zero_at_efermi
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
                Dos object
        """
        energies = dos.energies - dos.efermi if self.zero_at_efermi \
            else dos.energies
        densities = dos.get_smeared_densities(self.sigma) if self.sigma \
            else dos.densities
        efermi = dos.efermi
        self._doses[label] = {'energies': energies, 'densities': densities,
                              'efermi': efermi}

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
            dict: Dict of dos data. Generally of the form
            {label: {'energies':..., 'densities': {'up':...}, 'efermi':efermi}}
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

        ncolors = max(3, len(self._doses))
        ncolors = min(9, ncolors)

        import palettable

        colors = palettable.colorbrewer.qualitative.Set1_9.mpl_colors

        y = None
        alldensities = []
        allenergies = []
        plt = pretty_plot(12, 8)

        # Note that this complicated processing of energies is to allow for
        # stacked plots in matplotlib.
        for key, dos in self._doses.items():
            energies = dos['energies']
            densities = dos['densities']
            if not y:
                y = {Spin.up: np.zeros(energies.shape),
                     Spin.down: np.zeros(energies.shape)}
            newdens = {}
            for spin in [Spin.up, Spin.down]:
                if spin in densities:
                    if self.stack:
                        y[spin] += densities[spin]
                        newdens[spin] = y[spin].copy()
                    else:
                        newdens[spin] = densities[spin]
            allenergies.append(energies)
            alldensities.append(newdens)

        keys = list(self._doses.keys())
        keys.reverse()
        alldensities.reverse()
        allenergies.reverse()
        allpts = []
        for i, key in enumerate(keys):
            x = []
            y = []
            for spin in [Spin.up, Spin.down]:
                if spin in alldensities[i]:
                    densities = list(int(spin) * alldensities[i][spin])
                    energies = list(allenergies[i])
                    if spin == Spin.down:
                        energies.reverse()
                        densities.reverse()
                    x.extend(energies)
                    y.extend(densities)
            allpts.extend(list(zip(x, y)))
            if self.stack:
                plt.fill(x, y, color=colors[i % ncolors],
                         label=str(key))
            else:
                plt.plot(x, y, color=colors[i % ncolors],
                         label=str(key), linewidth=3)
            if not self.zero_at_efermi:
                ylim = plt.ylim()
                plt.plot([self._doses[key]['efermi'],
                          self._doses[key]['efermi']], ylim,
                         color=colors[i % ncolors],
                         linestyle='--', linewidth=2)

        if xlim:
            plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)
        else:
            xlim = plt.xlim()
            relevanty = [p[1] for p in allpts
                         if xlim[0] < p[0] < xlim[1]]
            plt.ylim((min(relevanty), max(relevanty)))

        if self.zero_at_efermi:
            ylim = plt.ylim()
            plt.plot([0, 0], ylim, 'k--', linewidth=2)

        plt.xlabel('Energies (eV)')
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


class BSPlotter:
    """
    Class to plot or get data to facilitate the plot of band structure objects.

    Args:
        bs: A BandStructureSymmLine object.
    """

    def __init__(self, bs):
        if not isinstance(bs, BandStructureSymmLine):
            raise ValueError(
                "BSPlotter only works with BandStructureSymmLine objects. "
                "A BandStructure object (on a uniform grid for instance and "
                "not along symmetry lines won't work)")
        self._bs = bs
        # TODO: come with an intelligent way to cut the highest unconverged
        # bands
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

    def bs_plot_data(self, zero_to_efermi=True):

        """
        Get the data nicely formatted for a plot

        Args:
            zero_to_efermi: Automatically subtract off the Fermi energy from the
                eigenvalues and plot.

        Returns:
            dict: A dictionary of the following format:
            ticks: A dict with the 'distances' at which there is a kpoint (the
            x axis) and the labels (None if no label).
            energy: A dict storing bands for spin up and spin down data
            [{Spin:[band_index][k_point_index]}] as a list (one element
            for each branch) of energy for each kpoint. The data is
            stored by branch to facilitate the plotting.
            vbm: A list of tuples (distance,energy) marking the vbms. The
            energies are shifted with respect to the fermi level is the
            option has been selected.
            cbm: A list of tuples (distance,energy) marking the cbms. The
            energies are shifted with respect to the fermi level is the
            option has been selected.
            lattice: The reciprocal lattice.
            zero_energy: This is the energy used as zero for the plot.
            band_gap:A string indicating the band gap and its nature (empty if
            it's a metal).
            is_metal: True if the band structure is metallic (i.e., there is at
            least one band crossing the fermi level).
        """
        distance = []
        energy = []
        if self._bs.is_metal():
            zero_energy = self._bs.efermi
        else:
            zero_energy = self._bs.get_vbm()['energy']

        if not zero_to_efermi:
            zero_energy = 0.0

        for b in self._bs.branches:

            if self._bs.is_spin_polarized:
                energy.append({str(Spin.up): [], str(Spin.down): []})
            else:
                energy.append({str(Spin.up): []})
            distance.append([self._bs.distance[j]
                             for j in range(b['start_index'],
                                            b['end_index'] + 1)])
            ticks = self.get_ticks()

            for i in range(self._nb_bands):
                energy[-1][str(Spin.up)].append(
                    [self._bs.bands[Spin.up][i][j] - zero_energy
                     for j in range(b['start_index'], b['end_index'] + 1)])
            if self._bs.is_spin_polarized:
                for i in range(self._nb_bands):
                    energy[-1][str(Spin.down)].append(
                        [self._bs.bands[Spin.down][i][j] - zero_energy
                         for j in range(b['start_index'], b['end_index'] + 1)])

        vbm = self._bs.get_vbm()
        cbm = self._bs.get_cbm()

        vbm_plot = []
        cbm_plot = []

        for index in cbm['kpoint_index']:
            cbm_plot.append((self._bs.distance[index],
                             cbm['energy'] - zero_energy if zero_to_efermi
                             else cbm['energy']))

        for index in vbm['kpoint_index']:
            vbm_plot.append((self._bs.distance[index],
                             vbm['energy'] - zero_energy if zero_to_efermi
                             else vbm['energy']))

        bg = self._bs.get_band_gap()
        direct = "Indirect"
        if bg['direct']:
            direct = "Direct"

        return {'ticks': ticks, 'distances': distance, 'energy': energy,
                'vbm': vbm_plot, 'cbm': cbm_plot,
                'lattice': self._bs.lattice_rec.as_dict(),
                'zero_energy': zero_energy, 'is_metal': self._bs.is_metal(),
                'band_gap': "{} {} bandgap = {}".format(direct,
                                                        bg['transition'],
                                                        bg['energy'])
                if not self._bs.is_metal() else ""}

    def get_plot(self, zero_to_efermi=True, ylim=None, smooth=False,
                 vbm_cbm_marker=False, smooth_tol=None):
        """
        Get a matplotlib object for the bandstructure plot.
        Blue lines are up spin, red lines are down
        spin.

        Args:
            zero_to_efermi: Automatically subtract off the Fermi energy from
                the eigenvalues and plot (E-Ef).
            ylim: Specify the y-axis (energy) limits; by default None let
                the code choose. It is vbm-4 and cbm+4 if insulator
                efermi-10 and efermi+10 if metal
            smooth: interpolates the bands by a spline cubic
            smooth_tol (float) : tolerance for fitting spline to band data.
                Default is None such that no tolerance will be used.
        """
        plt = pretty_plot(12, 8)
        from matplotlib import rc
        import scipy.interpolate as scint

        # main internal config options
        e_min = -4
        e_max = 4
        if self._bs.is_metal():
            e_min = -10
            e_max = 10
        # band_linewidth = 3
        band_linewidth = 1

        data = self.bs_plot_data(zero_to_efermi)
        if not smooth:
            for d in range(len(data['distances'])):
                for i in range(self._nb_bands):
                    plt.plot(data['distances'][d],
                             [data['energy'][d][str(Spin.up)][i][j]
                              for j in range(len(data['distances'][d]))], 'b-',
                             linewidth=band_linewidth)
                    if self._bs.is_spin_polarized:
                        plt.plot(data['distances'][d],
                                 [data['energy'][d][str(Spin.down)][i][j]
                                  for j in range(len(data['distances'][d]))],
                                 'r--', linewidth=band_linewidth)
        else:
            # Interpolation failure can be caused by trying to fit an entire
            # band with one spline rather than fitting with piecewise splines
            # (splines are ill-suited to fit discontinuities).
            #
            # The number of splines used to fit a band is determined by the 
            # number of branches (high symmetry lines) defined in the 
            # BandStructureSymmLine object (see BandStructureSymmLine._branches). 

            warning = "WARNING! Distance / branch {d}, band {i} cannot be " + \
                      "interpolated.\n" + \
                      "See full warning in source.\n" + \
                      "If this is not a mistake, try increasing " + \
                      "smooth_tol.\nCurrent smooth_tol is {s}."

            for d in range(len(data['distances'])):
                for i in range(self._nb_bands):
                    tck = scint.splrep(
                        data['distances'][d],
                        [data['energy'][d][str(Spin.up)][i][j]
                         for j in range(len(data['distances'][d]))],
                        s=smooth_tol)
                    step = (data['distances'][d][-1]
                            - data['distances'][d][0]) / 1000

                    xs = [x * step + data['distances'][d][0]
                          for x in range(1000)]

                    ys = [scint.splev(x * step + data['distances'][d][0],
                                      tck, der=0)
                          for x in range(1000)]

                    for y in ys:
                        if np.isnan(y):
                            print(warning.format(d=str(d), i=str(i),
                                                 s=str(smooth_tol)))
                            break

                    plt.plot(xs, ys, 'b-', linewidth=band_linewidth)

                    if self._bs.is_spin_polarized:
                        tck = scint.splrep(
                            data['distances'][d],
                            [data['energy'][d][str(Spin.down)][i][j]
                             for j in range(len(data['distances'][d]))],
                            s=smooth_tol)
                        step = (data['distances'][d][-1]
                                - data['distances'][d][0]) / 1000

                        xs = [x * step + data['distances'][d][0]
                              for x in range(1000)]

                        ys = [scint.splev(
                            x * step + data['distances'][d][0],
                            tck, der=0)
                            for x in range(1000)]

                        for y in ys:
                            if np.isnan(y):
                                print(warning.format(d=str(d), i=str(i),
                                                     s=str(smooth_tol)))
                                break

                        plt.plot(xs, ys, 'r--', linewidth=band_linewidth)

        self._maketicks(plt)

        # Main X and Y Labels
        plt.xlabel(r'$\mathrm{Wave\ Vector}$', fontsize=30)
        ylabel = r'$\mathrm{E\ -\ E_f\ (eV)}$' if zero_to_efermi \
            else r'$\mathrm{Energy\ (eV)}$'
        plt.ylabel(ylabel, fontsize=30)

        # Draw Fermi energy, only if not the zero
        if not zero_to_efermi:
            ef = self._bs.efermi
            plt.axhline(ef, linewidth=2, color='k')

        # X range (K)
        # last distance point
        x_max = data['distances'][-1][-1]
        plt.xlim(0, x_max)

        if ylim is None:
            if self._bs.is_metal():
                # Plot A Metal
                if zero_to_efermi:
                    plt.ylim(e_min, e_max)
                else:
                    plt.ylim(self._bs.efermi + e_min, self._bs.efermi + e_max)
            else:
                if vbm_cbm_marker:
                    for cbm in data['cbm']:
                        plt.scatter(cbm[0], cbm[1], color='r', marker='o',
                                    s=100)
                    for vbm in data['vbm']:
                        plt.scatter(vbm[0], vbm[1], color='g', marker='o',
                                    s=100)
                plt.ylim(data['vbm'][0][1] + e_min,
                         data['cbm'][0][1] + e_max)
        else:
            plt.ylim(ylim)
            if not self._bs.is_metal() and vbm_cbm_marker:
                for cbm in data['cbm']:
                        plt.scatter(cbm[0], cbm[1], color='r', marker='o',
                                    s=100)
                for vbm in data['vbm']:
                        plt.scatter(vbm[0], vbm[1], color='g', marker='o',
                                    s=100)

        plt.tight_layout()

        return plt

    def show(self, zero_to_efermi=True, ylim=None, smooth=False,
             smooth_tol=None):
        """
        Show the plot using matplotlib.

        Args:
            zero_to_efermi: Automatically subtract off the Fermi energy from
                the eigenvalues and plot (E-Ef).
            ylim: Specify the y-axis (energy) limits; by default None let
                the code choose. It is vbm-4 and cbm+4 if insulator
                efermi-10 and efermi+10 if metal
            smooth: interpolates the bands by a spline cubic
            smooth_tol (float) : tolerance for fitting spline to band data.
                Default is None such that no tolerance will be used.
        """
        plt = self.get_plot(zero_to_efermi, ylim, smooth)
        plt.show()

    def save_plot(self, filename, img_format="eps", ylim=None,
                  zero_to_efermi=True, smooth=False):
        """
        Save matplotlib plot to a file.

        Args:
            filename: Filename to write to.
            img_format: Image format to use. Defaults to EPS.
            ylim: Specifies the y-axis limits.
        """
        plt = self.get_plot(ylim=ylim, zero_to_efermi=zero_to_efermi,
                            smooth=smooth)
        plt.savefig(filename, format=img_format)
        plt.close()

    def get_ticks(self):
        """
        Get all ticks and labels for a band structure plot.

        Returns:
            dict: A dictionary with 'distance': a list of distance at which
            ticks should be set and 'label': a list of label for each of those
            ticks.
        """
        tick_distance = []
        tick_labels = []
        previous_label = self._bs.kpoints[0].label
        previous_branch = self._bs.branches[0]['name']
        for i, c in enumerate(self._bs.kpoints):
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

    def plot_compare(self, other_plotter, legend=True):
        """
        plot two band structure for comparison. One is in red the other in blue
        (no difference in spins). The two band structures need to be defined
        on the same symmetry lines! and the distance between symmetry lines is
        the one of the band structure used to build the BSPlotter

        Args:
            another band structure object defined along the same symmetry lines

        Returns:
            a matplotlib object with both band structures

        """
        # TODO: add exception if the band structures are not compatible
        import matplotlib.lines as mlines
        plt = self.get_plot()
        data_orig = self.bs_plot_data()
        data = other_plotter.bs_plot_data()
        band_linewidth = 1
        for i in range(other_plotter._nb_bands):
            for d in range(len(data_orig['distances'])):
                plt.plot(data_orig['distances'][d],
                         [e[str(Spin.up)][i] for e in data['energy']][d],
                         'c-', linewidth=band_linewidth)
                if other_plotter._bs.is_spin_polarized:
                    plt.plot(data_orig['distances'][d],
                             [e[str(Spin.down)][i] for e in data['energy']][d],
                             'm--', linewidth=band_linewidth)
        if legend:
            handles = [mlines.Line2D([], [], linewidth=2,
                                     color='b', label='bs 1 up'),
                       mlines.Line2D([], [], linewidth=2,
                                     color='r', label='bs 1 down',
                                     linestyle="--"),
                       mlines.Line2D([], [], linewidth=2,
                                     color='c', label='bs 2 up'),
                       mlines.Line2D([], [], linewidth=2,
                                     color='m', linestyle="--",
                                     label='bs 2 down')]

            plt.legend(handles=handles)
        return plt

    def plot_brillouin(self):
        """
        plot the Brillouin zone
        """

        # get labels and lines
        labels = {}
        for k in self._bs.kpoints:
            if k.label:
                labels[k.label] = k.frac_coords

        lines = []
        for b in self._bs.branches:
            lines.append([self._bs.kpoints[b['start_index']].frac_coords,
                          self._bs.kpoints[b['end_index']].frac_coords])

        plot_brillouin_zone(self._bs.lattice_rec, lines=lines, labels=labels)


class BSPlotterProjected(BSPlotter):
    """
    Class to plot or get data to facilitate the plot of band structure objects
    projected along orbitals, elements or sites.

    Args:
        bs: A BandStructureSymmLine object with projections.
    """

    def __init__(self, bs):
        if len(bs.projections) == 0:
            raise ValueError("try to plot projections"
                             " on a band structure without any")
        super(BSPlotterProjected, self).__init__(bs)

    def _get_projections_by_branches(self, dictio):
        proj = self._bs.get_projections_on_elements_and_orbitals(dictio)
        proj_br = []
        for b in self._bs.branches:
            if self._bs.is_spin_polarized:
                proj_br.append(
                    {str(Spin.up): [[] for l in range(self._nb_bands)],
                     str(Spin.down): [[] for l in range(self._nb_bands)]})
            else:
                proj_br.append(
                    {str(Spin.up): [[] for l in range(self._nb_bands)]})

            for i in range(self._nb_bands):
                for j in range(b['start_index'], b['end_index'] + 1):
                    proj_br[-1][str(Spin.up)][i].append(
                        {e: {o: proj[Spin.up][i][j][e][o]
                             for o in proj[Spin.up][i][j][e]}
                         for e in proj[Spin.up][i][j]})
            if self._bs.is_spin_polarized:
                for b in self._bs.branches:
                    for i in range(self._nb_bands):
                        for j in range(b['start_index'], b['end_index'] + 1):
                            proj_br[-1][str(Spin.down)][i].append(
                                {e: {o: proj[Spin.down][i][j][e][o]
                                     for o in proj[Spin.down][i][j][e]}
                                 for e in proj[Spin.down][i][j]})
        return proj_br

    def get_projected_plots_dots(self, dictio, zero_to_efermi=True, ylim=None,
                                 vbm_cbm_marker=False):
        """
        Method returning a plot composed of subplots along different elements
        and orbitals.

        Args:
            dictio: The element and orbitals you want a projection on. The
                format is {Element:[Orbitals]} for instance
                {'Cu':['d','s'],'O':['p']} will give projections for Cu on
                d and s orbitals and on oxygen p.

        Returns:
            a pylab object with different subfigures for each projection
            The blue and red colors are for spin up and spin down.
            The bigger the red or blue dot in the band structure the higher
            character for the corresponding element and orbital.
        """
        band_linewidth = 1.0
        fig_number = sum([len(v) for v in dictio.values()])
        proj = self._get_projections_by_branches(dictio)
        data = self.bs_plot_data(zero_to_efermi)
        plt = pretty_plot(12, 8)
        e_min = -4
        e_max = 4
        if self._bs.is_metal():
            e_min = -10
            e_max = 10
        count = 1

        for el in dictio:
            for o in dictio[el]:
                plt.subplot(100 * math.ceil(fig_number / 2) + 20 + count)
                self._maketicks(plt)
                for b in range(len(data['distances'])):
                    for i in range(self._nb_bands):
                        plt.plot(data['distances'][b],
                                 [data['energy'][b][str(Spin.up)][i][j]
                                  for j in range(len(data['distances'][b]))],
                                 'b-',
                                 linewidth=band_linewidth)
                        if self._bs.is_spin_polarized:
                            plt.plot(data['distances'][b],
                                     [data['energy'][b][str(Spin.down)][i][j]
                                      for j in
                                      range(len(data['distances'][b]))],
                                     'r--', linewidth=band_linewidth)
                            for j in range(
                                    len(data['energy'][b][str(Spin.up)][i])):
                                plt.plot(data['distances'][b][j],
                                         data['energy'][b][str(Spin.down)][i][
                                             j], 'ro',
                                         markersize=
                                         proj[b][str(Spin.down)][i][j][str(el)][
                                             o] * 15.0)
                        for j in range(len(data['energy'][b][str(Spin.up)][i])):
                            plt.plot(data['distances'][b][j],
                                     data['energy'][b][str(Spin.up)][i][j],
                                     'bo',
                                     markersize=
                                     proj[b][str(Spin.up)][i][j][str(el)][
                                         o] * 15.0)
                if ylim is None:
                    if self._bs.is_metal():
                        if zero_to_efermi:
                            plt.ylim(e_min, e_max)
                        else:
                            plt.ylim(self._bs.efermi + e_min, self._bs.efermi
                                     + e_max)
                    else:
                        if vbm_cbm_marker:
                            for cbm in data['cbm']:
                                plt.scatter(cbm[0], cbm[1], color='r',
                                            marker='o',
                                            s=100)

                            for vbm in data['vbm']:
                                plt.scatter(vbm[0], vbm[1], color='g',
                                            marker='o',
                                            s=100)

                        plt.ylim(data['vbm'][0][1] + e_min, data['cbm'][0][1]
                                 + e_max)
                else:
                    plt.ylim(ylim)
                plt.title(str(el) + " " + str(o))
                count += 1
        return plt

    def get_elt_projected_plots(self, zero_to_efermi=True, ylim=None,
                                vbm_cbm_marker=False):
        """
        Method returning a plot composed of subplots along different elements

        Returns:
            a pylab object with different subfigures for each projection
            The blue and red colors are for spin up and spin down
            The bigger the red or blue dot in the band structure the higher
            character for the corresponding element and orbital
        """
        band_linewidth = 1.0
        proj = self._get_projections_by_branches({e.symbol: ['s', 'p', 'd']
                                                  for e in
                                                  self._bs.structure.composition.elements})
        data = self.bs_plot_data(zero_to_efermi)
        plt = pretty_plot(12, 8)
        e_min = -4
        e_max = 4
        if self._bs.is_metal():
            e_min = -10
            e_max = 10
        count = 1
        for el in self._bs.structure.composition.elements:
            plt.subplot(220 + count)
            self._maketicks(plt)
            for b in range(len(data['distances'])):
                for i in range(self._nb_bands):
                    plt.plot(data['distances'][b],
                             [data['energy'][b][str(Spin.up)][i][j]
                              for j in range(len(data['distances'][b]))], '-',
                             color=[192 / 255, 192 / 255, 192 / 255],
                             linewidth=band_linewidth)
                    if self._bs.is_spin_polarized:
                        plt.plot(data['distances'][b],
                                 [data['energy'][b][str(Spin.down)][i][j]
                                  for j in range(len(data['distances'][b]))],
                                 '--', color=[128 / 255, 128 / 255, 128 / 255],
                                 linewidth=band_linewidth)
                        for j in range(len(data['energy'][b][str(Spin.up)][i])):
                            markerscale = sum([proj[b][str(Spin.down)][i][
                                                   j][str(el)][o] for o in
                                               proj[b]
                                               [str(Spin.down)][i][j][
                                                   str(el)]])
                            plt.plot(data['distances'][b][j],
                                     data['energy'][b][str(Spin.down)][i][j],
                                     'bo',
                                     markersize=markerscale * 15.0,
                                     color=[markerscale, 0.3 * markerscale,
                                            0.4 * markerscale])
                    for j in range(len(data['energy'][b][str(Spin.up)][i])):
                        markerscale = sum(
                            [proj[b][str(Spin.up)][i][j][str(el)][o]
                             for o in proj[b]
                             [str(Spin.up)][i][j][str(el)]])
                        plt.plot(data['distances'][b][j],
                                 data['energy'][b][str(Spin.up)][i][j], 'o',
                                 markersize=markerscale * 15.0,
                                 color=[markerscale, 0.3 * markerscale,
                                        0.4 * markerscale])
            if ylim is None:
                if self._bs.is_metal():
                    if zero_to_efermi:
                        plt.ylim(e_min, e_max)
                    else:
                        plt.ylim(self._bs.efermi + e_min, self._bs.efermi
                                 + e_max)
                else:
                    if vbm_cbm_marker:
                        for cbm in data['cbm']:
                            plt.scatter(cbm[0], cbm[1], color='r', marker='o',
                                        s=100)

                        for vbm in data['vbm']:
                            plt.scatter(vbm[0], vbm[1], color='g', marker='o',
                                        s=100)

                    plt.ylim(data['vbm'][0][1] + e_min, data['cbm'][0][1]
                             + e_max)
            else:
                plt.ylim(ylim)
            plt.title(str(el))
            count += 1

        return plt

    def get_elt_projected_plots_color(self, zero_to_efermi=True,
                                      elt_ordered=None):
        """
        returns a pylab plot object with one plot where the band structure
        line color depends on the character of the band (along different
        elements). Each element is associated with red, green or blue
        and the corresponding rgb color depending on the character of the band
        is used. The method can only deal with binary and ternary compounds

        spin up and spin down are differientiated by a '-' and a '--' line

        Args:
            elt_ordered: A list of Element ordered. The first one is red,
                second green, last blue

        Returns:
            a pylab object

        """
        band_linewidth = 3.0
        if len(self._bs.structure.composition.elements) > 3:
            raise ValueError
        if elt_ordered is None:
            elt_ordered = self._bs.structure.composition.elements
        proj = self._get_projections_by_branches(
            {e.symbol: ['s', 'p', 'd']
             for e in self._bs.structure.composition.elements})
        data = self.bs_plot_data(zero_to_efermi)
        plt = pretty_plot(12, 8)

        spins = [Spin.up]
        if self._bs.is_spin_polarized:
            spins = [Spin.up, Spin.down]
        self._maketicks(plt)
        for s in spins:
            for b in range(len(data['distances'])):
                for i in range(self._nb_bands):
                    for j in range(len(data['energy'][b][str(s)][i]) - 1):
                        sum_e = 0.0
                        for el in elt_ordered:
                            sum_e = sum_e + \
                                    sum([proj[b][str(s)][i][j][str(el)][o]
                                         for o
                                         in proj[b][str(s)][i][j][str(el)]])
                        if sum_e == 0.0:
                            color = [0.0] * len(elt_ordered)
                        else:
                            color = [sum([proj[b][str(s)][i][j][str(el)][o]
                                          for o
                                          in proj[b][str(s)][i][j][str(el)]])
                                     / sum_e
                                     for el in elt_ordered]
                        if len(color) == 2:
                            color.append(0.0)
                            color[2] = color[1]
                            color[1] = 0.0
                        sign = '-'
                        if s == Spin.down:
                            sign = '--'
                        plt.plot([data['distances'][b][j],
                                  data['distances'][b][j + 1]],
                                 [data['energy'][b][str(s)][i][j],
                                  data['energy'][b][str(s)][i][j + 1]], sign,
                                 color=color, linewidth=band_linewidth)

        if self._bs.is_metal():
            if zero_to_efermi:
                e_min = -10
                e_max = 10
                plt.ylim(e_min, e_max)
                plt.ylim(self._bs.efermi + e_min, self._bs.efermi + e_max)
        else:
            plt.ylim(data['vbm'][0][1] - 4.0, data['cbm'][0][1] + 2.0)
        return plt

    def _get_projections_by_branches_patom_pmorb(self, dictio, dictpa,
                                                 sum_atoms, sum_morbs,
                                                 selected_branches):
        import copy
        setos = {'s': 0, 'py': 1, 'pz': 2, 'px': 3, 'dxy': 4, 'dyz': 5,
                 'dz2': 6, 'dxz': 7,
                 'dx2': 8, 'f_3': 9, 'f_2': 10, 'f_1': 11, 'f0': 12, 'f1': 13,
                 'f2': 14, 'f3': 15}

        num_branches = len(self._bs.branches)
        if selected_branches is not None:
            indices = []
            if not isinstance(selected_branches, list):
                raise TypeError(
                    "You do not give a correct type of 'selected_branches'. It should be 'list' type.")
            elif len(selected_branches) == 0:
                raise ValueError(
                    "The 'selected_branches' is empty. We cannot do anything.")
            else:
                for index in selected_branches:
                    if not isinstance(index, int):
                        raise ValueError(
                            "You do not give a correct type of index of symmetry lines. It should be "
                            "'int' type")
                    elif index > num_branches or index < 1:
                        raise ValueError(
                            "You give a incorrect index of symmetry lines: %s. The index should be in "
                            "range of [1, %s]." % (
                            str(index), str(num_branches)))
                    else:
                        indices.append(index - 1)
        else:
            indices = range(0, num_branches)

        proj = self._bs.projections
        proj_br = []
        for index in indices:
            b = self._bs.branches[index]
            print(b)
            if self._bs.is_spin_polarized:
                proj_br.append(
                    {str(Spin.up): [[] for l in range(self._nb_bands)],
                     str(Spin.down): [[] for l in range(self._nb_bands)]})
            else:
                proj_br.append(
                    {str(Spin.up): [[] for l in range(self._nb_bands)]})

            for i in range(self._nb_bands):
                for j in range(b['start_index'], b['end_index'] + 1):
                    edict = {}
                    for elt in dictpa:
                        for anum in dictpa[elt]:
                            edict[elt + str(anum)] = {}
                            for morb in dictio[elt]:
                                edict[elt + str(anum)][morb] = \
                                proj[Spin.up][i][j][setos[morb]][anum - 1]
                    proj_br[-1][str(Spin.up)][i].append(edict)

            if self._bs.is_spin_polarized:
                for i in range(self._nb_bands):
                    for j in range(b['start_index'], b['end_index'] + 1):
                        edict = {}
                        for elt in dictpa:
                            for anum in dictpa[elt]:
                                edict[elt + str(anum)] = {}
                                for morb in dictio[elt]:
                                    edict[elt + str(anum)][morb] = \
                                    proj[Spin.up][i][j][setos[morb]][anum - 1]
                        proj_br[-1][str(Spin.down)][i].append(edict)

        # Adjusting  projections for plot
        dictio_d, dictpa_d = self._summarize_keys_for_plot(dictio, dictpa,
                                                           sum_atoms, sum_morbs)
        print('dictio_d: %s' % str(dictio_d))
        print('dictpa_d: %s' % str(dictpa_d))

        if (sum_atoms is None) and (sum_morbs is None):
            proj_br_d = copy.deepcopy(proj_br)
        else:
            proj_br_d = []
            branch = -1
            for index in indices:
                branch += 1
                br = self._bs.branches[index]
                if self._bs.is_spin_polarized:
                    proj_br_d.append(
                        {str(Spin.up): [[] for l in range(self._nb_bands)],
                         str(Spin.down): [[] for l in range(self._nb_bands)]})
                else:
                    proj_br_d.append(
                        {str(Spin.up): [[] for l in range(self._nb_bands)]})

                if (sum_atoms is not None) and (sum_morbs is None):
                    for i in range(self._nb_bands):
                        for j in range(br['end_index'] - br['start_index'] + 1):
                            atoms_morbs = copy.deepcopy(
                                proj_br[branch][str(Spin.up)][i][j])
                            edict = {}
                            for elt in dictpa:
                                if elt in sum_atoms:
                                    for anum in dictpa_d[elt][:-1]:
                                        edict[elt + anum] = copy.deepcopy(
                                            atoms_morbs[elt + anum])
                                    edict[elt + dictpa_d[elt][-1]] = {}
                                    for morb in dictio[elt]:
                                        sprojection = 0.0
                                        for anum in sum_atoms[elt]:
                                            sprojection += \
                                            atoms_morbs[elt + str(anum)][morb]
                                        edict[elt + dictpa_d[elt][-1]][
                                            morb] = sprojection
                                else:
                                    for anum in dictpa_d[elt]:
                                        edict[elt + anum] = copy.deepcopy(
                                            atoms_morbs[elt + anum])
                            proj_br_d[-1][str(Spin.up)][i].append(edict)
                    if self._bs.is_spin_polarized:
                        for i in range(self._nb_bands):
                            for j in range(br['end_index'] - br[
                                'start_index'] + 1):
                                atoms_morbs = copy.deepcopy(
                                    proj_br[branch][str(Spin.down)][i][j])
                                edict = {}
                                for elt in dictpa:
                                    if elt in sum_atoms:
                                        for anum in dictpa_d[elt][:-1]:
                                            edict[elt + anum] = copy.deepcopy(
                                                atoms_morbs[elt + anum])
                                        edict[elt + dictpa_d[elt][-1]] = {}
                                        for morb in dictio[elt]:
                                            sprojection = 0.0
                                            for anum in sum_atoms[elt]:
                                                sprojection += \
                                                atoms_morbs[elt + str(anum)][
                                                    morb]
                                            edict[elt + dictpa_d[elt][-1]][
                                                morb] = sprojection
                                    else:
                                        for anum in dictpa_d[elt]:
                                            edict[elt + anum] = copy.deepcopy(
                                                atoms_morbs[elt + anum])
                                proj_br_d[-1][str(Spin.down)][i].append(edict)

                elif (sum_atoms is None) and (sum_morbs is not None):
                    for i in range(self._nb_bands):
                        for j in range(br['end_index'] - br['start_index'] + 1):
                            atoms_morbs = copy.deepcopy(
                                proj_br[branch][str(Spin.up)][i][j])
                            edict = {}
                            for elt in dictpa:
                                if elt in sum_morbs:
                                    for anum in dictpa_d[elt]:
                                        edict[elt + anum] = {}
                                        for morb in dictio_d[elt][:-1]:
                                            edict[elt + anum][morb] = \
                                            atoms_morbs[elt + anum][morb]
                                        sprojection = 0.0
                                        for morb in sum_morbs[elt]:
                                            sprojection += \
                                            atoms_morbs[elt + anum][morb]
                                        edict[elt + anum][
                                            dictio_d[elt][-1]] = sprojection
                                else:
                                    for anum in dictpa_d[elt]:
                                        edict[elt + anum] = copy.deepcopy(
                                            atoms_morbs[elt + anum])
                            proj_br_d[-1][str(Spin.up)][i].append(edict)
                    if self._bs.is_spin_polarized:
                        for i in range(self._nb_bands):
                            for j in range(br['end_index'] - br[
                                'start_index'] + 1):
                                atoms_morbs = copy.deepcopy(
                                    proj_br[branch][str(Spin.down)][i][j])
                                edict = {}
                                for elt in dictpa:
                                    if elt in sum_morbs:
                                        for anum in dictpa_d[elt]:
                                            edict[elt + anum] = {}
                                            for morb in dictio_d[elt][:-1]:
                                                edict[elt + anum][morb] = \
                                                atoms_morbs[elt + anum][morb]
                                            sprojection = 0.0
                                            for morb in sum_morbs[elt]:
                                                sprojection += \
                                                atoms_morbs[elt + anum][morb]
                                            edict[elt + anum][
                                                dictio_d[elt][-1]] = sprojection
                                    else:
                                        for anum in dictpa_d[elt]:
                                            edict[elt + anum] = copy.deepcopy(
                                                atoms_morbs[elt + anum])
                                proj_br_d[-1][str(Spin.down)][i].append(edict)

                else:
                    for i in range(self._nb_bands):
                        for j in range(br['end_index'] - br['start_index'] + 1):
                            atoms_morbs = copy.deepcopy(
                                proj_br[branch][str(Spin.up)][i][j])
                            edict = {}
                            for elt in dictpa:
                                if (elt in sum_atoms) and (elt in sum_morbs):
                                    for anum in dictpa_d[elt][:-1]:
                                        edict[elt + anum] = {}
                                        for morb in dictio_d[elt][:-1]:
                                            edict[elt + anum][morb] = \
                                            atoms_morbs[elt + anum][morb]
                                        sprojection = 0.0
                                        for morb in sum_morbs[elt]:
                                            sprojection += \
                                            atoms_morbs[elt + anum][morb]
                                        edict[elt + anum][
                                            dictio_d[elt][-1]] = sprojection

                                    edict[elt + dictpa_d[elt][-1]] = {}
                                    for morb in dictio_d[elt][:-1]:
                                        sprojection = 0.0
                                        for anum in sum_atoms[elt]:
                                            sprojection += \
                                            atoms_morbs[elt + str(anum)][morb]
                                        edict[elt + dictpa_d[elt][-1]][
                                            morb] = sprojection

                                    sprojection = 0.0
                                    for anum in sum_atoms[elt]:
                                        for morb in sum_morbs[elt]:
                                            sprojection += \
                                            atoms_morbs[elt + str(anum)][morb]
                                    edict[elt + dictpa_d[elt][-1]][
                                        dictio_d[elt][-1]] = sprojection

                                elif (elt in sum_atoms) and (
                                    elt not in sum_morbs):
                                    for anum in dictpa_d[elt][:-1]:
                                        edict[elt + anum] = copy.deepcopy(
                                            atoms_morbs[elt + anum])
                                    edict[elt + dictpa_d[elt][-1]] = {}
                                    for morb in dictio[elt]:
                                        sprojection = 0.0
                                        for anum in sum_atoms[elt]:
                                            sprojection += \
                                            atoms_morbs[elt + str(anum)][morb]
                                        edict[elt + dictpa_d[elt][-1]][
                                            morb] = sprojection

                                elif (elt not in sum_atoms) and (
                                    elt in sum_morbs):
                                    for anum in dictpa_d[elt]:
                                        edict[elt + anum] = {}
                                        for morb in dictio_d[elt][:-1]:
                                            edict[elt + anum][morb] = \
                                            atoms_morbs[elt + anum][morb]
                                        sprojection = 0.0
                                        for morb in sum_morbs[elt]:
                                            sprojection += \
                                            atoms_morbs[elt + anum][morb]
                                        edict[elt + anum][
                                            dictio_d[elt][-1]] = sprojection

                                else:
                                    for anum in dictpa_d[elt]:
                                        edict[elt + anum] = {}
                                        for morb in dictio_d[elt]:
                                            edict[elt + anum][morb] = \
                                            atoms_morbs[elt + anum][morb]
                            proj_br_d[-1][str(Spin.up)][i].append(edict)

                    if self._bs.is_spin_polarized:
                        for i in range(self._nb_bands):
                            for j in range(br['end_index'] - br[
                                'start_index'] + 1):
                                atoms_morbs = copy.deepcopy(
                                    proj_br[branch][str(Spin.down)][i][j])
                                edict = {}
                                for elt in dictpa:
                                    if (elt in sum_atoms) and (
                                        elt in sum_morbs):
                                        for anum in dictpa_d[elt][:-1]:
                                            edict[elt + anum] = {}
                                            for morb in dictio_d[elt][:-1]:
                                                edict[elt + anum][morb] = \
                                                atoms_morbs[elt + anum][morb]
                                            sprojection = 0.0
                                            for morb in sum_morbs[elt]:
                                                sprojection += \
                                                atoms_morbs[elt + anum][morb]
                                            edict[elt + anum][
                                                dictio_d[elt][-1]] = sprojection

                                        edict[elt + dictpa_d[elt][-1]] = {}
                                        for morb in dictio_d[elt][:-1]:
                                            sprojection = 0.0
                                            for anum in sum_atoms[elt]:
                                                sprojection += \
                                                atoms_morbs[elt + str(anum)][
                                                    morb]
                                            edict[elt + dictpa_d[elt][-1]][
                                                morb] = sprojection

                                        sprojection = 0.0
                                        for anum in sum_atoms[elt]:
                                            for morb in sum_morbs[elt]:
                                                sprojection += \
                                                atoms_morbs[elt + str(anum)][
                                                    morb]
                                        edict[elt + dictpa_d[elt][-1]][
                                            dictio_d[elt][-1]] = sprojection

                                    elif (elt in sum_atoms) and (
                                        elt not in sum_morbs):
                                        for anum in dictpa_d[elt][:-1]:
                                            edict[elt + anum] = copy.deepcopy(
                                                atoms_morbs[elt + anum])
                                        edict[elt + dictpa_d[elt][-1]] = {}
                                        for morb in dictio[elt]:
                                            sprojection = 0.0
                                            for anum in sum_atoms[elt]:
                                                sprojection += \
                                                atoms_morbs[elt + str(anum)][
                                                    morb]
                                            edict[elt + dictpa_d[elt][-1]][
                                                morb] = sprojection

                                    elif (elt not in sum_atoms) and (
                                        elt in sum_morbs):
                                        for anum in dictpa_d[elt]:
                                            edict[elt + anum] = {}
                                            for morb in dictio_d[elt][:-1]:
                                                edict[elt + anum][morb] = \
                                                atoms_morbs[elt + anum][morb]
                                            sprojection = 0.0
                                            for morb in sum_morbs[elt]:
                                                sprojection += \
                                                atoms_morbs[elt + anum][morb]
                                            edict[elt + anum][
                                                dictio_d[elt][-1]] = sprojection

                                    else:
                                        for anum in dictpa_d[elt]:
                                            edict[elt + anum] = {}
                                            for morb in dictio_d[elt]:
                                                edict[elt + anum][morb] = \
                                                atoms_morbs[elt + anum][morb]
                                proj_br_d[-1][str(Spin.down)][i].append(edict)

        return proj_br_d, dictio_d, dictpa_d, indices

    def get_projected_plots_dots_patom_pmorb(self, dictio, dictpa,
                                             sum_atoms=None, sum_morbs=None,
                                             zero_to_efermi=True, ylim=None,
                                             vbm_cbm_marker=False,
                                             selected_branches=None,
                                             w_h_size=(12, 8), num_column=None):
        """
        Method returns a plot composed of subplots for different atoms and
        orbitals (subshell orbitals such as 's', 'p', 'd' and 'f' defined by
        azimuthal quantum numbers l = 0, 1, 2 and 3, respectively or
        individual orbitals like 'px', 'py' and 'pz' defined by magnetic
        quantum numbers m = -1, 1 and 0, respectively).
        This is an extension of "get_projected_plots_dots" method.

        Args:
            dictio: The elements and the orbitals you need to project on. The
                format is {Element:[Orbitals]}, for instance:
                {'Cu':['dxy','s','px'],'O':['px','py','pz']} will give
                projections for Cu on orbitals dxy, s, px and
                for O on orbitals px, py, pz. If you want to sum over all
                individual orbitals of subshell orbitals,
                for example, 'px', 'py' and 'pz' of O, just simply set
                {'Cu':['dxy','s','px'],'O':['p']} and set sum_morbs (see
                explanations below) as {'O':[p],...}.
                Otherwise, you will get an error.
            dictpa: The elements and their sites (defined by site numbers) you
                need to project on. The format is
                {Element: [Site numbers]}, for instance: {'Cu':[1,5],'O':[3,4]}
                will give projections for Cu on site-1
                and on site-5, O on site-3 and on site-4 in the cell.
                Attention:
                The correct site numbers of atoms are consistent with
                themselves in the structure computed. Normally,
                the structure should be totally similar with POSCAR file,
                however, sometimes VASP can rotate or
                translate the cell. Thus, it would be safe if using Vasprun
                class to get the final_structure and as a
                result, correct index numbers of atoms.
            sum_atoms: Sum projection of the similar atoms together (e.g.: Cu
                on site-1 and Cu on site-5). The format is
                {Element: [Site numbers]}, for instance:
                 {'Cu': [1,5], 'O': [3,4]} means summing projections over Cu on
                 site-1 and Cu on site-5 and O on site-3
                 and on site-4. If you do not want to use this functional, just
                 turn it off by setting sum_atoms = None.
            sum_morbs: Sum projections of individual orbitals of similar atoms
                together (e.g.: 'dxy' and 'dxz'). The
                format is {Element: [individual orbitals]}, for instance:
                {'Cu': ['dxy', 'dxz'], 'O': ['px', 'py']} means summing
                projections over 'dxy' and 'dxz' of Cu and 'px'
                and 'py' of O. If you do not want to use this functional, just
                turn it off by setting sum_morbs = None.
            selected_branches: The index of symmetry lines you chose for
                plotting. This can be useful when the number of
                symmetry lines (in KPOINTS file) are manny while you only want
                to show for certain ones. The format is
                [index of line], for instance:
                [1, 3, 4] means you just need to do projection along lines
                number 1, 3 and 4 while neglecting lines
                number 2 and so on. By default, this is None type and all
                symmetry lines will be plotted.
            w_h_size: This variable help you to control the width and height
                of figure. By default, width = 12 and
                height = 8 (inches). The width/height ratio is kept the same
                for subfigures and the size of each depends
                on how many number of subfigures are plotted.
            num_column: This variable help you to manage how the subfigures are
                arranged in the figure by setting
                up the number of columns of subfigures. The value should be an
                int number. For example, num_column = 3
                means you want to plot subfigures in 3 columns. By default,
                num_column = None and subfigures are
                aligned in 2 columns.

        Returns:
            A pylab object with different subfigures for different projections.
            The blue and red colors lines are bands
            for spin up and spin down. The green and cyan dots are projections
            for spin up and spin down. The bigger
            the green or cyan dots in the projected band structures, the higher
            character for the corresponding elements
            and orbitals. List of individual orbitals and their numbers (set up
            by VASP and no special meaning):
            s = 0; py = 1 pz = 2 px = 3; dxy = 4 dyz = 5 dz2 = 6 dxz = 7 dx2 = 8;
            f_3 = 9 f_2 = 10 f_1 = 11 f0 = 12 f1 = 13 f2 = 14 f3 = 15
        """
        dictio, sum_morbs = self._Orbitals_SumOrbitals(dictio, sum_morbs)
        dictpa, sum_atoms, number_figs = self._number_of_subfigures(dictio,
                                                                    dictpa,
                                                                    sum_atoms,
                                                                    sum_morbs)
        print('Number of subfigures: %s' % str(number_figs))
        if number_figs > 9:
            print(
                "The number of sub-figures %s might be too manny and the implementation might take a long time.\n"
                "A smaller number or a plot with selected symmetry lines (selected_branches) might be better.\n"
                % str(number_figs))
        import math
        from pymatgen.util.plotting import pretty_plot
        band_linewidth = 0.5
        plt = pretty_plot(w_h_size[0], w_h_size[1])
        proj_br_d, dictio_d, dictpa_d, branches = self._get_projections_by_branches_patom_pmorb(
            dictio, dictpa,
            sum_atoms, sum_morbs, selected_branches)
        data = self.bs_plot_data(zero_to_efermi)
        e_min = -4
        e_max = 4
        if self._bs.is_metal():
            e_min = -10
            e_max = 10

        count = 0
        for elt in dictpa_d:
            for numa in dictpa_d[elt]:
                for o in dictio_d[elt]:

                    count += 1
                    if num_column is None:
                        if number_figs == 1:
                            plt.subplot(1, 1, 1)
                        else:
                            row = number_figs / 2
                            if number_figs % 2 == 0:
                                plt.subplot(row, 2, count)
                            else:
                                plt.subplot(row + 1, 2, count)
                    elif isinstance(num_column, int):
                        row = number_figs / num_column
                        if number_figs % num_column == 0:
                            plt.subplot(row, num_column, count)
                        else:
                            plt.subplot(row + 1, num_column, count)
                    else:
                        raise ValueError(
                            "The invalid 'num_column' is assigned. It should be an integer.")

                    plt, shift = self._maketicks_selected(plt, branches)
                    br = -1
                    for b in branches:
                        br += 1
                        for i in range(self._nb_bands):
                            plt.plot(list(map(lambda x: x - shift[br], data['distances'][b])),
                                     [data['energy'][b][str(Spin.up)][i][j]
                                      for j in
                                      range(len(data['distances'][b]))],
                                     'b-', linewidth=band_linewidth)

                            if self._bs.is_spin_polarized:
                                plt.plot(list(map(lambda x: x - shift[br], data['distances'][b])),
                                         [data['energy'][b][str(Spin.down)][i][
                                              j]
                                          for j in
                                          range(len(data['distances'][b]))],
                                         'r--', linewidth=band_linewidth)
                                for j in range(len(
                                        data['energy'][b][str(Spin.up)][i])):
                                    plt.plot(
                                        data['distances'][b][j] - shift[br],
                                        data['energy'][b][str(Spin.down)][i][j],
                                        'co', markersize= \
                                            proj_br_d[br][str(Spin.down)][i][j][
                                                elt + numa][o] * 15.0)

                            for j in range(
                                    len(data['energy'][b][str(Spin.up)][i])):
                                plt.plot(data['distances'][b][j] - shift[br],
                                         data['energy'][b][str(Spin.up)][i][j],
                                         'go', markersize= \
                                             proj_br_d[br][str(Spin.up)][i][j][
                                                 elt + numa][o] * 15.0)

                    if ylim is None:
                        if self._bs.is_metal():
                            if zero_to_efermi:
                                plt.ylim(e_min, e_max)
                            else:
                                plt.ylim(self._bs.efermi + e_min,
                                         self._bs._efermi
                                         + e_max)
                        else:
                            if vbm_cbm_marker:
                                for cbm in data['cbm']:
                                    plt.scatter(cbm[0], cbm[1], color='r',
                                                marker='o',
                                                s=100)

                                for vbm in data['vbm']:
                                    plt.scatter(vbm[0], vbm[1], color='g',
                                                marker='o',
                                                s=100)

                            plt.ylim(data['vbm'][0][1] + e_min,
                                     data['cbm'][0][1]
                                     + e_max)
                    else:
                        plt.ylim(ylim)
                    plt.title(elt + " " + numa + " " + str(o))

        return plt

    def _Orbitals_SumOrbitals(self, dictio, sum_morbs):
        all_orbitals = ['s', 'p', 'd', 'f', 'px', 'py', 'pz', 'dxy', 'dyz',
                        'dxz', 'dx2', 'dz2',
                        'f_3', 'f_2', 'f_1', 'f0', 'f1', 'f2', 'f3']
        individual_orbs = {'p': ['px', 'py', 'pz'],
                           'd': ['dxy', 'dyz', 'dxz', 'dx2', 'dz2'],
                           'f': ['f_3', 'f_2', 'f_1', 'f0', 'f1', 'f2', 'f3']}

        if not isinstance(dictio, dict):
            raise TypeError(
                "The invalid type of 'dictio' was bound. It should be dict type.")
        elif len(dictio.keys()) == 0:
            raise KeyError("The 'dictio' is empty. We cannot do anything.")
        else:
            for elt in dictio:
                if Element.is_valid_symbol(elt):
                    if isinstance(dictio[elt], list):
                        if len(dictio[elt]) == 0:
                            raise ValueError(
                                "The dictio[%s] is empty. We cannot do anything" % elt)
                        for orb in dictio[elt]:
                            if not isinstance(orb, str):
                                raise ValueError(
                                    "The invalid format of orbitals is in 'dictio[%s]': %s. "
                                    "They should be string." % (elt, str(orb)))
                            elif orb not in all_orbitals:
                                raise ValueError(
                                    "The invalid name of orbital is given in 'dictio[%s]'." % elt)
                            else:
                                if orb in individual_orbs.keys():
                                    if len(set(dictio[elt]).intersection(
                                            individual_orbs[orb])) != 0:
                                        raise ValueError(
                                            "The 'dictio[%s]' contains orbitals repeated." % elt)
                        nelems = Counter(dictio[elt]).values()
                        if sum(nelems) > len(nelems):
                            raise ValueError(
                                "You put in at least two similar orbitals in dictio[%s]." % elt)
                    else:
                        raise TypeError(
                            "The invalid type of value was put into 'dictio[%s]'. It should be list "
                            "type." % elt)
                else:
                    raise KeyError(
                        "The invalid element was put into 'dictio' as a key: %s" % elt)

        if sum_morbs is None:
            print("You do not want to sum projection over orbitals.")
        elif not isinstance(sum_morbs, dict):
            raise TypeError(
                "The invalid type of 'sum_orbs' was bound. It should be dict or 'None' type.")
        elif len(sum_morbs.keys()) == 0:
            raise KeyError("The 'sum_morbs' is empty. We cannot do anything")
        else:
            for elt in sum_morbs:
                if Element.is_valid_symbol(elt):
                    if isinstance(sum_morbs[elt], list):
                        for orb in sum_morbs[elt]:
                            if not isinstance(orb, str):
                                raise TypeError(
                                    "The invalid format of orbitals is in 'sum_morbs[%s]': %s. "
                                    "They should be string." % (elt, str(orb)))
                            elif orb not in all_orbitals:
                                raise ValueError(
                                    "The invalid name of orbital in 'sum_morbs[%s]' is given." % elt)
                            else:
                                if orb in individual_orbs.keys():
                                    if len(set(sum_morbs[elt]).intersection(
                                            individual_orbs[orb])) != 0:
                                        raise ValueError(
                                            "The 'sum_morbs[%s]' contains orbitals repeated." % elt)
                        nelems = Counter(sum_morbs[elt]).values()
                        if sum(nelems) > len(nelems):
                            raise ValueError(
                                "You put in at least two similar orbitals in sum_morbs[%s]." % elt)
                    else:
                        raise TypeError(
                            "The invalid type of value was put into 'sum_morbs[%s]'. It should be list "
                            "type." % elt)
                    if elt not in dictio.keys():
                        raise ValueError(
                            "You cannot sum projection over orbitals of atoms '%s' because they are not "
                            "mentioned in 'dictio'." % elt)
                else:
                    raise KeyError(
                        "The invalid element was put into 'sum_morbs' as a key: %s" % elt)

        for elt in dictio:
            if len(dictio[elt]) == 1:
                if len(dictio[elt][0]) > 1:
                    if elt in sum_morbs.keys():
                        raise ValueError(
                            "You cannot sum projection over one individual orbital '%s' of '%s'." %
                            (dictio[elt][0], elt))
                else:
                    if sum_morbs is None:
                        pass
                    elif elt not in sum_morbs.keys():
                        print(
                            "You do not want to sum projection over orbitals of element: %s" % elt)
                    else:
                        if len(sum_morbs[elt]) == 0:
                            raise ValueError(
                                "The empty list is an invalid value for sum_morbs[%s]." % elt)
                        elif len(sum_morbs[elt]) > 1:
                            for orb in sum_morbs[elt]:
                                if dictio[elt][0] not in orb:
                                    raise ValueError(
                                        "The invalid orbital '%s' was put into 'sum_morbs[%s]'." %
                                        (orb, elt))
                        else:
                            if orb == 's' or len(orb) > 1:
                                raise ValueError(
                                    "The invalid orbital '%s' was put into sum_orbs['%s']." % (
                                    orb, elt))
                            else:
                                sum_morbs[elt] = individual_orbs[dictio[elt][0]]
                                dictio[elt] = individual_orbs[dictio[elt][0]]
            else:
                duplicate = copy.deepcopy(dictio[elt])
                for orb in dictio[elt]:
                    if orb in individual_orbs.keys():
                        duplicate.remove(orb)
                        for o in individual_orbs[orb]:
                            duplicate.append(o)
                dictio[elt] = copy.deepcopy(duplicate)

                if sum_morbs is None:
                    pass
                elif elt not in sum_morbs.keys():
                    print(
                        "You do not want to sum projection over orbitals of element: %s" % elt)
                else:
                    if len(sum_morbs[elt]) == 0:
                        raise ValueError(
                            "The empty list is an invalid value for sum_morbs[%s]." % elt)
                    elif len(sum_morbs[elt]) == 1:
                        orb = sum_morbs[elt][0]
                        if orb == 's':
                            raise ValueError(
                                "We do not sum projection over only 's' orbital of the same "
                                "type of element.")
                        elif orb in individual_orbs.keys():
                            sum_morbs[elt].pop(0)
                            for o in individual_orbs[orb]:
                                sum_morbs[elt].append(o)
                        else:
                            raise ValueError(
                                "You never sum projection over one orbital in sum_morbs[%s]" % elt)
                    else:
                        duplicate = copy.deepcopy(sum_morbs[elt])
                        for orb in sum_morbs[elt]:
                            if orb in individual_orbs.keys():
                                duplicate.remove(orb)
                                for o in individual_orbs[orb]:
                                    duplicate.append(o)
                        sum_morbs[elt] = copy.deepcopy(duplicate)

                    for orb in sum_morbs[elt]:
                        if orb not in dictio[elt]:
                            raise ValueError(
                                "The orbitals of sum_morbs[%s] conflict with those of dictio[%s]." %
                                (elt, elt))

        return dictio, sum_morbs

    def _number_of_subfigures(self, dictio, dictpa, sum_atoms, sum_morbs):
        from pymatgen.core.periodic_table import Element
        from collections import Counter

        if (not isinstance(dictpa, dict)):
            raise TypeError(
                "The invalid type of 'dictpa' was bound. It should be dict type.")
        elif len(dictpa.keys()) == 0:
            raise KeyError("The 'dictpa' is empty. We cannot do anything.")
        else:
            for elt in dictpa:
                if Element.is_valid_symbol(elt):
                    if isinstance(dictpa[elt], list):
                        if len(dictpa[elt]) == 0:
                            raise ValueError(
                                "The dictpa[%s] is empty. We cannot do anything" % elt)
                        _sites = self._bs.structure.sites
                        indices = []
                        for i in range(0, len(_sites)):
                            if list(_sites[i]._species.keys())[0].__eq__(
                                    Element(elt)):
                                indices.append(i + 1)
                        for number in dictpa[elt]:
                            if isinstance(number, str):
                                if 'all' == number.lower():
                                    dictpa[elt] = indices
                                    print(
                                        "You want to consider all '%s' atoms." % elt)
                                    break
                                else:
                                    raise ValueError(
                                        "You put wrong site numbers in 'dictpa[%s]': %s." %
                                        (elt, str(number)))
                            elif isinstance(number, int):
                                if number not in indices:
                                    raise ValueError(
                                        "You put wrong site numbers in 'dictpa[%s]': %s." %
                                        (elt, str(number)))
                            else:
                                raise ValueError(
                                    "You put wrong site numbers in 'dictpa[%s]': %s." % (
                                    elt, str(number)))
                        nelems = Counter(dictpa[elt]).values()
                        if sum(nelems) > len(nelems):
                            raise ValueError(
                                "You put at least two similar site numbers into 'dictpa[%s]'." % elt)
                    else:
                        raise TypeError(
                            "The invalid type of value was put into 'dictpa[%s]'. It should be list "
                            "type." % elt)
                else:
                    raise KeyError(
                        "The invalid element was put into 'dictpa' as a key: %s" % elt)

        if len(list(dictio.keys())) != len(list(dictpa.keys())):
            raise KeyError(
                "The number of keys in 'dictio' and 'dictpa' are not the same.")
        else:
            for elt in dictio.keys():
                if elt not in dictpa.keys(): raise KeyError(
                    "The element '%s' is not in both dictpa and dictio." % elt)
            for elt in dictpa.keys():
                if elt not in dictio.keys(): raise KeyError(
                    "The element '%s' in not in both dictpa and dictio." % elt)

        if sum_atoms is None:
            print("You do not want to sum projection over atoms.")
        elif (not isinstance(sum_atoms, dict)):
            raise TypeError(
                "The invalid type of 'sum_atoms' was bound. It should be dict type.")
        elif len(sum_atoms.keys()) == 0:
            raise KeyError("The 'sum_atoms' is empty. We cannot do anything.")
        else:
            for elt in sum_atoms:
                if Element.is_valid_symbol(elt):
                    if isinstance(sum_atoms[elt], list):
                        if len(sum_atoms[elt]) == 0:
                            raise ValueError(
                                "The sum_atoms[%s] is empty. We cannot do anything" % elt)
                        _sites = self._bs.structure.sites
                        indices = []
                        for i in range(0, len(_sites)):
                            if list(_sites[i]._species.keys())[0].__eq__(
                                    Element(elt)):
                                indices.append(i + 1)
                        for number in sum_atoms[elt]:
                            if isinstance(number, str):
                                if 'all' == number.lower():
                                    sum_atoms[elt] = indices
                                    print(
                                        "You want to sum projection over all '%s' atoms." % elt)
                                    break
                                else:
                                    raise ValueError(
                                        "You put wrong site numbers in 'sum_atoms[%s]'." % elt)
                            elif isinstance(number, int):
                                if number not in indices:
                                    raise ValueError(
                                        "You put wrong site numbers in 'sum_atoms[%s]'." % elt)
                                elif number not in dictpa[elt]:
                                    raise ValueError(
                                        "You cannot sum projection with atom number '%s' because it is not "
                                        "metioned in dicpta[%s]" % (
                                        str(number), elt))
                            else:
                                raise ValueError(
                                    "You put wrong site numbers in 'sum_atoms[%s]'." % elt)
                        nelems = Counter(sum_atoms[elt]).values()
                        if sum(nelems) > len(nelems):
                            raise ValueError(
                                "You put at least two similar site numbers into 'sum_atoms[%s]'." % elt)
                    else:
                        raise TypeError(
                            "The invalid type of value was put into 'sum_atoms[%s]'. It should be list "
                            "type." % elt)
                    if elt not in dictpa.keys():
                        raise ValueError(
                            "You cannot sum projection over atoms '%s' because it is not "
                            "mentioned in 'dictio'." % elt)
                else:
                    raise KeyError(
                        "The invalid element was put into 'sum_atoms' as a key: %s" % elt)
                if len(sum_atoms[elt]) == 1:
                    raise ValueError(
                        "We do not sum projection over only one atom: %s" % elt)

        max_number_figs = 0
        decrease = 0
        for elt in dictio:
            max_number_figs += len(dictio[elt]) * len(dictpa[elt])

        if (sum_atoms is None) and (sum_morbs is None):
            number_figs = max_number_figs
        elif (sum_atoms is not None) and (sum_morbs is None):
            for elt in sum_atoms:
                decrease += (len(sum_atoms[elt]) - 1) * len(dictio[elt])
            number_figs = max_number_figs - decrease
        elif (sum_atoms is None) and (sum_morbs is not None):
            for elt in sum_morbs:
                decrease += (len(sum_morbs[elt]) - 1) * len(dictpa[elt])
            number_figs = max_number_figs - decrease
        elif (sum_atoms is not None) and (sum_morbs is not None):
            for elt in sum_atoms:
                decrease += (len(sum_atoms[elt]) - 1) * len(dictio[elt])
            for elt in sum_morbs:
                if elt in sum_atoms:
                    decrease += (len(sum_morbs[elt]) - 1) * (
                    len(dictpa[elt]) - len(sum_atoms[elt]) + 1)
                else:
                    decrease += (len(sum_morbs[elt]) - 1) * len(dictpa[elt])
            number_figs = max_number_figs - decrease
        else:
            raise ValueError("Invalid format of 'sum_atoms' and 'sum_morbs'.")

        return dictpa, sum_atoms, number_figs

    def _summarize_keys_for_plot(self, dictio, dictpa, sum_atoms, sum_morbs):
        from pymatgen.core.periodic_table import Element
        individual_orbs = {'p': ['px', 'py', 'pz'],
                           'd': ['dxy', 'dyz', 'dxz', 'dx2', 'dz2'],
                           'f': ['f_3', 'f_2', 'f_1', 'f0', 'f1', 'f2', 'f3']}

        def number_label(list_numbers):
            list_numbers = sorted(list_numbers)
            divide = [[]]
            divide[0].append(list_numbers[0])
            group = 0
            for i in range(1, len(list_numbers)):
                if list_numbers[i] == list_numbers[i - 1] + 1:
                    divide[group].append(list_numbers[i])
                else:
                    group += 1
                    divide.append([list_numbers[i]])
            label = ""
            for elem in divide:
                if len(elem) > 1:
                    label += str(elem[0]) + "-" + str(elem[-1]) + ","
                else:
                    label += str(elem[0]) + ","
            return label[:-1]

        def orbital_label(list_orbitals):
            divide = {}
            for orb in list_orbitals:
                if orb[0] in divide:
                    divide[orb[0]].append(orb)
                else:
                    divide[orb[0]] = []
                    divide[orb[0]].append(orb)
            label = ""
            for elem in divide:
                if elem == 's':
                    label += "s" + ","
                else:
                    if len(divide[elem]) == len(individual_orbs[elem]):
                        label += elem + ","
                    else:
                        l = [o[1:] for o in divide[elem]]
                        label += elem + str(l).replace("['", "").replace("']",
                                                                         "").replace(
                            "', '", "-") + ","
            return label[:-1]

        if (sum_atoms is None) and (sum_morbs is None):
            dictio_d = dictio
            dictpa_d = {elt: [str(anum) for anum in dictpa[elt]] for elt in
                        dictpa}

        elif (sum_atoms is not None) and (sum_morbs is None):
            dictio_d = dictio
            dictpa_d = {}
            for elt in dictpa:
                dictpa_d[elt] = []
                if elt in sum_atoms:
                    _sites = self._bs.structure.sites
                    indices = []
                    for i in range(0, len(_sites)):
                        if list(_sites[i]._species.keys())[0].__eq__(Element(elt)):
                            indices.append(i + 1)
                    flag_1 = len(set(dictpa[elt]).intersection(indices))
                    flag_2 = len(set(sum_atoms[elt]).intersection(indices))
                    if flag_1 == len(indices) and flag_2 == len(indices):
                        dictpa_d[elt].append('all')
                    else:
                        for anum in dictpa[elt]:
                            if anum not in sum_atoms[elt]:
                                dictpa_d[elt].append(str(anum))
                        label = number_label(sum_atoms[elt])
                        dictpa_d[elt].append(label)
                else:
                    for anum in dictpa[elt]:
                        dictpa_d[elt].append(str(anum))

        elif (sum_atoms is None) and (sum_morbs is not None):
            dictio_d = {}
            for elt in dictio:
                dictio_d[elt] = []
                if elt in sum_morbs:
                    for morb in dictio[elt]:
                        if morb not in sum_morbs[elt]:
                            dictio_d[elt].append(morb)
                    label = orbital_label(sum_morbs[elt])
                    dictio_d[elt].append(label)
                else:
                    dictio_d[elt] = dictio[elt]
            dictpa_d = {elt: [str(anum) for anum in dictpa[elt]] for elt in
                        dictpa}

        else:
            dictio_d = {}
            for elt in dictio:
                dictio_d[elt] = []
                if elt in sum_morbs:
                    for morb in dictio[elt]:
                        if morb not in sum_morbs[elt]:
                            dictio_d[elt].append(morb)
                    label = orbital_label(sum_morbs[elt])
                    dictio_d[elt].append(label)
                else:
                    dictio_d[elt] = dictio[elt]
            dictpa_d = {}
            for elt in dictpa:
                dictpa_d[elt] = []
                if elt in sum_atoms:
                    _sites = self._bs.structure.sites
                    indices = []
                    for i in range(0, len(_sites)):
                        if list(_sites[i]._species.keys())[0].__eq__(Element(elt)):
                            indices.append(i + 1)
                    flag_1 = len(set(dictpa[elt]).intersection(indices))
                    flag_2 = len(set(sum_atoms[elt]).intersection(indices))
                    if flag_1 == len(indices) and flag_2 == len(indices):
                        dictpa_d[elt].append('all')
                    else:
                        for anum in dictpa[elt]:
                            if anum not in sum_atoms[elt]:
                                dictpa_d[elt].append(str(anum))
                        label = number_label(sum_atoms[elt])
                        dictpa_d[elt].append(label)
                else:
                    for anum in dictpa[elt]:
                        dictpa_d[elt].append(str(anum))

        return dictio_d, dictpa_d

    def _maketicks_selected(self, plt, branches):
        """
        utility private method to add ticks to a band structure with selected branches
        """
        ticks = self.get_ticks()
        distance = []
        label = []
        rm_elems = []
        for i in range(1, len(ticks['distance'])):
            if ticks['label'][i] == ticks['label'][i - 1]:
                rm_elems.append(i)
        for i in range(len(ticks['distance'])):
            if i not in rm_elems:
                distance.append(ticks['distance'][i])
                label.append(ticks['label'][i])
        l_branches = [distance[i] - distance[i - 1] for i in
                      range(1, len(distance))]
        n_distance = []
        n_label = []
        for branch in branches:
            n_distance.append(l_branches[branch])
            if ("$\\mid$" not in label[branch]) and (
                "$\\mid$" not in label[branch + 1]):
                n_label.append([label[branch], label[branch + 1]])
            elif ("$\\mid$" in label[branch]) and (
                "$\\mid$" not in label[branch + 1]):
                n_label.append(
                    [label[branch].split("$")[-1], label[branch + 1]])
            elif ("$\\mid$" not in label[branch]) and (
                "$\\mid$" in label[branch + 1]):
                n_label.append([label[branch], label[branch + 1].split("$")[0]])
            else:
                n_label.append([label[branch].split("$")[-1],
                                label[branch + 1].split("$")[0]])

        f_distance = []
        rf_distance = []
        f_label = []
        f_label.append(n_label[0][0])
        f_label.append(n_label[0][1])
        f_distance.append(0.0)
        f_distance.append(n_distance[0])
        rf_distance.append(0.0)
        rf_distance.append(n_distance[0])
        length = n_distance[0]
        for i in range(1, len(n_distance)):
            if n_label[i][0] == n_label[i - 1][1]:
                f_distance.append(length)
                f_distance.append(length + n_distance[i])
                f_label.append(n_label[i][0])
                f_label.append(n_label[i][1])
            else:
                f_distance.append(length + n_distance[i])
                f_label[-1] = n_label[i - 1][1] + "$\\mid$" + n_label[i][0]
                f_label.append(n_label[i][1])
            rf_distance.append(length + n_distance[i])
            length += n_distance[i]

        n_ticks = {'distance': f_distance, 'label': f_label}
        uniq_d = []
        uniq_l = []
        temp_ticks = list(zip(n_ticks['distance'], n_ticks['label']))
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

        for i in range(len(n_ticks['label'])):
            if n_ticks['label'][i] is not None:
                # don't print the same label twice
                if i != 0:
                    if n_ticks['label'][i] == n_ticks['label'][i - 1]:
                        logger.debug("already print label... "
                                     "skipping label {i}".format(
                            i=n_ticks['label'][i]))
                    else:
                        logger.debug("Adding a line at {d}"
                                     " for label {l}".format(
                            d=n_ticks['distance'][i], l=n_ticks['label'][i]))
                        plt.axvline(n_ticks['distance'][i], color='k')
                else:
                    logger.debug("Adding a line at {d} for label {l}".format(
                        d=n_ticks['distance'][i], l=n_ticks['label'][i]))
                    plt.axvline(n_ticks['distance'][i], color='k')

        shift = []
        br = -1
        for branch in branches:
            br += 1
            shift.append(distance[branch] - rf_distance[br])

        return plt, shift


class BSDOSPlotter:
    """
    A joint, aligned band structure and density of states plot. Contributions 
    from Jan Pohls as well as the online example from Germain Salvato-Vallverdu:
    http://gvallver.perso.univ-pau.fr/?p=587
    """

    def __init__(self, bs_projection="elements", dos_projection="elements",
                 vb_energy_range=4, cb_energy_range=4, fixed_cb_energy=False,
                 egrid_interval=1, font="Times New Roman", axis_fontsize=20,
                 tick_fontsize=15, legend_fontsize=14, bs_legend="best",
                 dos_legend="best", rgb_legend=True, fig_size=(11, 8.5)):

        """
        Instantiate plotter settings.

        Args:
            bs_projection (str): "elements" or None
            dos_projection (str): "elements", "orbitals", or None
            vb_energy_range (float): energy in eV to show of valence bands
            cb_energy_range (float): energy in eV to show of conduction bands
            fixed_cb_energy (bool): If true, the cb_energy_range will be interpreted
                as constant (i.e., no gap correction for cb energy)
            egrid_interval (float): interval for grid marks
            font (str): font family
            axis_fontsize (float): font size for axis
            tick_fontsize (float): font size for axis tick labels
            legend_fontsize (float): font size for legends
            bs_legend (str): matplotlib string location for legend or None
            dos_legend (str): matplotlib string location for legend or None
            rgb_legend (bool): (T/F) whether to draw RGB triangle/bar for element proj.
            fig_size(tuple): dimensions of figure size (width, height)
        """
        self.bs_projection = bs_projection
        self.dos_projection = dos_projection
        self.vb_energy_range = vb_energy_range
        self.cb_energy_range = cb_energy_range
        self.fixed_cb_energy = fixed_cb_energy
        self.egrid_interval = egrid_interval
        self.font = font
        self.axis_fontsize = axis_fontsize
        self.tick_fontsize = tick_fontsize
        self.legend_fontsize = legend_fontsize
        self.bs_legend = bs_legend
        self.dos_legend = dos_legend
        self.rgb_legend = rgb_legend
        self.fig_size = fig_size

    def get_plot(self, bs, dos=None):
        """
        Get a matplotlib plot object.
        Args:
            bs (BandStructureSymmLine): the bandstructure to plot. Projection 
                data must exist for projected plots.
            dos (Dos): the Dos to plot. Projection data must exist (i.e., 
                CompleteDos) for projected plots.

        Returns:
            matplotlib.pyplot object on which you can call commands like show() 
            and savefig()
        """
        import matplotlib.lines as mlines
        from matplotlib.gridspec import GridSpec
        import matplotlib.pyplot as mplt

        # make sure the user-specified band structure projection is valid
        bs_projection = self.bs_projection
        if dos:
            elements = [e.symbol for e in dos.structure.composition.elements]
        elif bs_projection and bs.structure:
            elements = [e.symbol for e in bs.structure.composition.elements]
        else:
            elements = []

        rgb_legend = self.rgb_legend and bs_projection and \
                     bs_projection.lower() == "elements" and \
                     len(elements) in [2, 3]

        if bs_projection and bs_projection.lower() == "elements" and \
                (len(elements) not in [2, 3] or
                     not bs.get_projection_on_elements()):
            warnings.warn(
                "Cannot get element projected data; either the projection data "
                "doesn't exist, or you don't have a compound with exactly 2 "
                "or 3 unique elements.")
            bs_projection = None

        # specify energy range of plot
        emin = -self.vb_energy_range
        emax = self.cb_energy_range if self.fixed_cb_energy else \
            self.cb_energy_range + bs.get_band_gap()["energy"]

        # initialize all the k-point labels and k-point x-distances for bs plot
        xlabels = []  # all symmetry point labels on x-axis
        xlabel_distances = []  # positions of symmetry point x-labels
        x_distances = []  # x positions of kpoint data
        prev_right_klabel = None  # used to determine which branches require a midline separator

        for idx, l in enumerate(bs.branches):
            # get left and right kpoint labels of this branch
            left_k, right_k = l["name"].split("-")

            # add $ notation for LaTeX kpoint labels
            if left_k[0] == "\\" or "_" in left_k:
                left_k = "$" + left_k + "$"
            if right_k[0] == "\\" or "_" in right_k:
                right_k = "$" + right_k + "$"

            # add left k label to list of labels
            if prev_right_klabel is None:
                xlabels.append(left_k)
                xlabel_distances.append(0)
            elif prev_right_klabel != left_k:  # used for pipe separator
                xlabels[-1] = xlabels[-1] + "$\\mid$ " + left_k

            # add right k label to list of labels
            xlabels.append(right_k)
            prev_right_klabel = right_k

            # add x-coordinates for labels
            left_kpoint = bs.kpoints[l["start_index"]].cart_coords
            right_kpoint = bs.kpoints[l["end_index"]].cart_coords
            distance = np.linalg.norm(right_kpoint - left_kpoint)
            xlabel_distances.append(xlabel_distances[-1] + distance)

            # add x-coordinates for kpoint data
            npts = l["end_index"] - l["start_index"]
            distance_interval = distance / npts
            x_distances.append(xlabel_distances[-2])
            for i in range(npts):
                x_distances.append(x_distances[-1] + distance_interval)

        # set up bs and dos plot
        gs = GridSpec(1, 2, width_ratios=[2, 1]) if dos else GridSpec(1, 1)

        fig = mplt.figure(figsize=self.fig_size)
        fig.patch.set_facecolor('white')
        bs_ax = mplt.subplot(gs[0])
        if dos:
            dos_ax = mplt.subplot(gs[1])

        # set basic axes limits for the plot
        bs_ax.set_xlim(0, x_distances[-1])
        bs_ax.set_ylim(emin, emax)
        if dos:
            dos_ax.set_ylim(emin, emax)

        # add BS xticks, labels, etc.
        bs_ax.set_xticks(xlabel_distances)
        bs_ax.set_xticklabels(xlabels, size=self.tick_fontsize)
        bs_ax.set_xlabel('Wavevector $k$', fontsize=self.axis_fontsize,
                         family=self.font)
        bs_ax.set_ylabel('$E-E_F$ / eV', fontsize=self.axis_fontsize,
                         family=self.font)

        # add BS fermi level line at E=0 and gridlines
        bs_ax.hlines(y=0, xmin=0, xmax=x_distances[-1], color="k", lw=2)
        bs_ax.set_yticks(np.arange(emin, emax + 1E-5, self.egrid_interval))
        bs_ax.set_yticklabels(np.arange(emin, emax + 1E-5, self.egrid_interval),
                              size=self.tick_fontsize)
        bs_ax.set_axisbelow(True)
        bs_ax.grid(color=[0.5, 0.5, 0.5], linestyle='dotted', linewidth=1)
        if dos:
            dos_ax.set_yticks(np.arange(emin, emax + 1E-5, self.egrid_interval))
            dos_ax.set_yticklabels([])
            dos_ax.grid(color=[0.5, 0.5, 0.5], linestyle='dotted', linewidth=1)

        # renormalize the band energy to the Fermi level
        band_energies = {}
        for spin in (Spin.up, Spin.down):
            if spin in bs.bands:
                band_energies[spin] = []
                for band in bs.bands[spin]:
                    band_energies[spin].append([e - bs.efermi for e in band])

        # renormalize the DOS energies to Fermi level
        if dos:
            dos_energies = [e - dos.efermi for e in dos.energies]

        # get the projection data to set colors for the band structure
        colordata = self._get_colordata(bs, elements, bs_projection)

        # plot the colored band structure lines
        for spin in (Spin.up, Spin.down):
            if spin in band_energies:
                linestyles = "solid" if spin == Spin.up else "dotted"
                for band_idx, band in enumerate(band_energies[spin]):
                    self._rgbline(bs_ax, x_distances, band,
                                  colordata[spin][band_idx, :, 0],
                                  colordata[spin][band_idx, :, 1],
                                  colordata[spin][band_idx, :, 2],
                                  linestyles=linestyles)

        if dos:
            # Plot the DOS and projected DOS
            for spin in (Spin.up, Spin.down):
                if spin in dos.densities:
                    # plot the total DOS
                    dos_densities = dos.densities[spin] * int(spin)
                    label = "total" if spin == Spin.up else None
                    dos_ax.plot(dos_densities, dos_energies,
                                color=(0.6, 0.6, 0.6), label=label)
                    dos_ax.fill_betweenx(dos_energies, 0,dos_densities,
                                        color=(0.7, 0.7, 0.7),
                                        facecolor=(0.7, 0.7, 0.7))

                    if self.dos_projection is None:
                        pass

                    elif self.dos_projection.lower() == "elements":
                        # plot the atom-projected DOS
                        colors = ['b', 'r', 'g', 'm', 'y', 'c', 'k', 'w']
                        el_dos = dos.get_element_dos()
                        for idx, el in enumerate(elements):
                            dos_densities = el_dos[Element(el)].densities[
                                                spin] * int(spin)
                            label = el if spin == Spin.up else None
                            dos_ax.plot(dos_densities, dos_energies,
                                        color=colors[idx], label=label)

                    elif self.dos_projection.lower() == "orbitals":
                        # plot each of the atomic projected DOS
                        colors = ['b', 'r', 'g', 'm']
                        spd_dos = dos.get_spd_dos()
                        for idx, orb in enumerate([OrbitalType.s,
                                                   OrbitalType.p,
                                                   OrbitalType.d,
                                                   OrbitalType.f]):
                            if orb in spd_dos:
                                dos_densities = spd_dos[orb].densities[spin] * \
                                                int(spin)
                                label = orb if spin == Spin.up else None
                                dos_ax.plot(dos_densities, dos_energies,
                                            color=colors[idx], label=label)

            # get index of lowest and highest energy being plotted, used to help auto-scale DOS x-axis
            emin_idx = next(x[0] for x in enumerate(dos_energies) if
                            x[1] >= emin)
            emax_idx = len(dos_energies) - \
                       next(x[0] for x in enumerate(reversed(dos_energies))
                            if x[1] <= emax)

            # determine DOS x-axis range
            dos_xmin = 0 if Spin.down not in dos.densities else -max(
                dos.densities[Spin.down][emin_idx:emax_idx + 1] * 1.05)
            dos_xmax = max([max(dos.densities[Spin.up][emin_idx:emax_idx]) *
                            1.05, abs(dos_xmin)])

            # set up the DOS x-axis and add Fermi level line
            dos_ax.set_xlim(dos_xmin, dos_xmax)
            dos_ax.set_xticklabels([])
            dos_ax.hlines(y=0, xmin=dos_xmin, xmax=dos_xmax, color="k", lw=2)
            dos_ax.set_xlabel('DOS', fontsize=self.axis_fontsize,
                              family=self.font)

        # add legend for band structure
        if self.bs_legend and not rgb_legend:
            handles = []

            if bs_projection is None:
                handles = [mlines.Line2D([], [], linewidth=2,
                                         color='k', label='spin up'),
                           mlines.Line2D([], [], linewidth=2,
                                         color='b', linestyle="dotted",
                                         label='spin down')]

            elif bs_projection.lower() == "elements":
                colors = ['b', 'r', 'g']
                for idx, el in enumerate(elements):
                    handles.append(mlines.Line2D([], [],
                                                 linewidth=2,
                                                 color=colors[idx], label=el))

            bs_ax.legend(handles=handles, fancybox=True,
                         prop={'size': self.legend_fontsize,
                               'family': self.font}, loc=self.bs_legend)

        elif self.bs_legend and rgb_legend:
            if len(elements) == 2:
                self._rb_line(bs_ax, elements[1], elements[0],
                              loc=self.bs_legend)
            elif len(elements) == 3:
                self._rgb_triangle(bs_ax, elements[1], elements[2], elements[0],
                                   loc=self.bs_legend)

        # add legend for DOS
        if dos and self.dos_legend:
            dos_ax.legend(fancybox=True, prop={'size': self.legend_fontsize,
                                               'family': self.font},
                          loc=self.dos_legend)

        mplt.subplots_adjust(wspace=0.1)
        return mplt

    @staticmethod
    def _rgbline(ax, k, e, red, green, blue, alpha=1, linestyles="solid"):
        """
        An RGB colored line for plotting.
        creation of segments based on:
        http://nbviewer.ipython.org/urls/raw.github.com/dpsanders/matplotlib-examples/master/colorline.ipynb
        Args:
            ax: matplotlib axis
            k: x-axis data (k-points)
            e: y-axis data (energies)
            red: red data
            green: green data
            blue: blue data
            alpha: alpha values data
            linestyles: linestyle for plot (e.g., "solid" or "dotted")
        """
        from matplotlib.collections import LineCollection

        pts = np.array([k, e]).T.reshape(-1, 1, 2)
        seg = np.concatenate([pts[:-1], pts[1:]], axis=1)

        nseg = len(k) - 1
        r = [0.5 * (red[i] + red[i + 1]) for i in range(nseg)]
        g = [0.5 * (green[i] + green[i + 1]) for i in range(nseg)]
        b = [0.5 * (blue[i] + blue[i + 1]) for i in range(nseg)]
        a = np.ones(nseg, np.float) * alpha
        lc = LineCollection(seg, colors=list(zip(r, g, b, a)),
                            linewidth=2, linestyles=linestyles)
        ax.add_collection(lc)

    @staticmethod
    def _get_colordata(bs, elements, bs_projection):
        """
        Get color data, including projected band structures
        Args:
            bs: Bandstructure object
            elements: elements (in desired order) for setting to blue, red, green
            bs_projection: None for no projection, "elements" for element projection

        Returns:

        """
        contribs = {}
        if bs_projection and bs_projection.lower() == "elements":
            projections = bs.get_projection_on_elements()

        for spin in (Spin.up, Spin.down):
            if spin in bs.bands:
                contribs[spin] = []
                for band_idx in range(bs.nb_bands):
                    colors = []
                    for k_idx in range(len(bs.kpoints)):
                        if bs_projection and bs_projection.lower() == "elements":
                            c = [0, 0, 0]
                            projs = projections[spin][band_idx][k_idx]
                            # note: squared color interpolations are smoother
                            # see: https://youtu.be/LKnqECcg6Gw
                            projs = dict(
                                [(k, v ** 2) for k, v in projs.items()])
                            total = sum(projs.values())
                            if total > 0:
                                for idx, e in enumerate(elements):
                                    c[idx] = math.sqrt(projs[
                                                           e] / total)  # min is to handle round errors

                            c = [c[1], c[2],
                                 c[0]]  # prefer blue, then red, then green

                        else:
                            c = [0, 0, 0] if spin == Spin.up \
                                else [0, 0,
                                      1]  # black for spin up, blue for spin down

                        colors.append(c)

                    contribs[spin].append(colors)
                contribs[spin] = np.array(contribs[spin])

        return contribs

    @staticmethod
    def _rgb_triangle(ax, r_label, g_label, b_label, loc):
        """
        Draw an RGB triangle legend on the desired axis
        """
        if not loc in range(1, 11):
            loc = 2

        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        inset_ax = inset_axes(ax, width=1, height=1, loc=loc)
        mesh = 35
        x = []
        y = []
        color = []
        for r in range(0, mesh):
            for g in range(0, mesh):
                for b in range(0, mesh):
                    if not (r == 0 and b == 0 and g == 0):
                        r1 = r / (r + g + b)
                        g1 = g / (r + g + b)
                        b1 = b / (r + g + b)
                        x.append(0.33 * (2. * g1 + r1) / (r1 + b1 + g1))
                        y.append(0.33 * np.sqrt(3) * r1 / (r1 + b1 + g1))
                        rc = math.sqrt(r ** 2 / (r ** 2 + g ** 2 + b ** 2))
                        gc = math.sqrt(g ** 2 / (r ** 2 + g ** 2 + b ** 2))
                        bc = math.sqrt(b ** 2 / (r ** 2 + g ** 2 + b ** 2))
                        color.append([rc, gc, bc])

        # x = [n + 0.25 for n in x]  # nudge x coordinates
        # y = [n + (max_y - 1) for n in y]  # shift y coordinates to top
        # plot the triangle
        inset_ax.scatter(x, y, s=7, marker='.', edgecolor=color)
        inset_ax.set_xlim([-0.35, 1.00])
        inset_ax.set_ylim([-0.35, 1.00])

        # add the labels
        inset_ax.text(0.70, -0.2, g_label, fontsize=13,
                      family='Times New Roman', color=(0, 0, 0),
                      horizontalalignment='left')
        inset_ax.text(0.325, 0.70, r_label, fontsize=13,
                      family='Times New Roman', color=(0, 0, 0),
                      horizontalalignment='center')
        inset_ax.text(-0.05, -0.2, b_label, fontsize=13,
                      family='Times New Roman', color=(0, 0, 0),
                      horizontalalignment='right')

        inset_ax.get_xaxis().set_visible(False)
        inset_ax.get_yaxis().set_visible(False)

    @staticmethod
    def _rb_line(ax, r_label, b_label, loc):
        # Draw an rb bar legend on the desired axis

        if not loc in range(1, 11):
            loc = 2
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        inset_ax = inset_axes(ax, width=1.2, height=0.4, loc=loc)

        x = []
        y = []
        color = []
        for i in range(0, 1000):
            x.append(i / 1800. + 0.55)
            y.append(0)
            color.append([math.sqrt(c) for c in
                          [1 - (i / 1000) ** 2, 0, (i / 1000) ** 2]])

        # plot the bar
        inset_ax.scatter(x, y, s=250., marker='s', edgecolor=color)
        inset_ax.set_xlim([-0.1, 1.7])
        inset_ax.text(1.35, 0, b_label, fontsize=13,
                      family='Times New Roman', color=(0, 0, 0),
                      horizontalalignment="left", verticalalignment="center")
        inset_ax.text(0.30, 0, r_label, fontsize=13,
                      family='Times New Roman', color=(0, 0, 0),
                      horizontalalignment="right", verticalalignment="center")

        inset_ax.get_xaxis().set_visible(False)
        inset_ax.get_yaxis().set_visible(False)


class BoltztrapPlotter:
    # TODO: We need a unittest for this. Come on folks.
    """
    class containing methods to plot the data from Boltztrap.

    Args:
        bz: a BoltztrapAnalyzer object
    """

    def __init__(self, bz):
        self._bz = bz

    def _plot_doping(self, temp):
        import matplotlib.pyplot as plt
        if len(self._bz.doping) != 0:
            limit = 2.21e15
            plt.axvline(self._bz.mu_doping['n'][temp][0], linewidth=3.0,
                        linestyle="--")
            plt.text(self._bz.mu_doping['n'][temp][0] + 0.01,
                     limit,
                     "$n$=10$^{" + str(
                         math.log10(self._bz.doping['n'][0])) + "}$",
                     color='b')
            plt.axvline(self._bz.mu_doping['n'][temp][-1], linewidth=3.0,
                        linestyle="--")
            plt.text(self._bz.mu_doping['n'][temp][-1] + 0.01,
                     limit,
                     "$n$=10$^{" + str(math.log10(self._bz.doping['n'][-1]))
                     + "}$", color='b')
            plt.axvline(self._bz.mu_doping['p'][temp][0], linewidth=3.0,
                        linestyle="--")
            plt.text(self._bz.mu_doping['p'][temp][0] + 0.01,
                     limit,
                     "$p$=10$^{" + str(
                         math.log10(self._bz.doping['p'][0])) + "}$",
                     color='b')
            plt.axvline(self._bz.mu_doping['p'][temp][-1], linewidth=3.0,
                        linestyle="--")
            plt.text(self._bz.mu_doping['p'][temp][-1] + 0.01,
                     limit, "$p$=10$^{" +
                     str(math.log10(self._bz.doping['p'][-1])) + "}$",
                     color='b')

    def _plot_bg_limits(self):
        import matplotlib.pyplot as plt
        plt.axvline(0.0, color='k', linewidth=3.0)
        plt.axvline(self._bz.gap, color='k', linewidth=3.0)

    def plot_seebeck_eff_mass_mu(self, temps=[300], output='average',
                                 Lambda=0.5):
        """
        Plot respect to the chemical potential of the Seebeck effective mass
        calculated as explained in Ref.
        Gibbs, Z. M. et al., Effective mass and fermi surface complexity factor 
        from ab initio band structure calculations. 
        npj Computational Materials 3, 8 (2017).

        Args:
            output: 'average' returns the seebeck effective mass calculated
                using the average of the three diagonal components of the
                seebeck tensor. 'tensor' returns the seebeck effective mass
                respect to the three diagonal components of the seebeck tensor.
            temps:  list of temperatures of calculated seebeck.
            Lambda: fitting parameter used to model the scattering (0.5 means
                constant relaxation time).
        Returns:
            a matplotlib object
        """

        import matplotlib.pyplot as plt
        plt.figure(figsize=(9, 7))
        for T in temps:
            sbk_mass = self._bz.get_seebeck_eff_mass(output=output, temp=T,
                                                     Lambda=0.5)
            # remove noise inside the gap
            start = self._bz.mu_doping['p'][T][0]
            stop = self._bz.mu_doping['n'][T][0]
            mu_steps_1 = []
            mu_steps_2 = []
            sbk_mass_1 = []
            sbk_mass_2 = []
            for i, mu in enumerate(self._bz.mu_steps):
                if mu <= start:
                    mu_steps_1.append(mu)
                    sbk_mass_1.append(sbk_mass[i])
                elif mu >= stop:
                    mu_steps_2.append(mu)
                    sbk_mass_2.append(sbk_mass[i])

            plt.plot(mu_steps_1, sbk_mass_1, label=str(T) + 'K', linewidth=3.0)
            plt.plot(mu_steps_2, sbk_mass_2, linewidth=3.0)
            if output == 'average':
                plt.gca().get_lines()[1].set_c(plt.gca().get_lines()[0].get_c())
            elif output == 'tensor':
                plt.gca().get_lines()[3].set_c(plt.gca().get_lines()[0].get_c())
                plt.gca().get_lines()[4].set_c(plt.gca().get_lines()[1].get_c())
                plt.gca().get_lines()[5].set_c(plt.gca().get_lines()[2].get_c())

        plt.xlabel("E-E$_f$ (eV)", fontsize=30)
        plt.ylabel("Seebeck effective mass", fontsize=30)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        if output == 'tensor':
            plt.legend([str(i) + '_' + str(T) + 'K' for T in temps for i in
                        ('x', 'y', 'z')], fontsize=20)
        elif output == 'average':
            plt.legend(fontsize=20)
        plt.tight_layout()
        return plt

    def plot_complexity_factor_mu(self, temps=[300], output='average',
                                  Lambda=0.5):
        """
        Plot respect to the chemical potential of the Fermi surface complexity 
        factor calculated as explained in Ref.
        Gibbs, Z. M. et al., Effective mass and fermi surface complexity factor 
        from ab initio band structure calculations. 
        npj Computational Materials 3, 8 (2017).
        
        Args:
            output: 'average' returns the complexity factor calculated using the average
                    of the three diagonal components of the seebeck and conductivity tensors.
                    'tensor' returns the complexity factor respect to the three
                    diagonal components of seebeck and conductivity tensors.
            temps:  list of temperatures of calculated seebeck and conductivity.
            Lambda: fitting parameter used to model the scattering (0.5 means constant 
                    relaxation time).
        Returns:
            a matplotlib object
        """
        import matplotlib.pyplot as plt
        plt.figure(figsize=(9, 7))
        for T in temps:
            cmplx_fact = self._bz.get_complexity_factor(output=output, temp=T,
                                                        Lambda=Lambda)
            start = self._bz.mu_doping['p'][T][0]
            stop = self._bz.mu_doping['n'][T][0]
            mu_steps_1 = []
            mu_steps_2 = []
            cmplx_fact_1 = []
            cmplx_fact_2 = []
            for i, mu in enumerate(self._bz.mu_steps):
                if mu <= start:
                    mu_steps_1.append(mu)
                    cmplx_fact_1.append(cmplx_fact[i])
                elif mu >= stop:
                    mu_steps_2.append(mu)
                    cmplx_fact_2.append(cmplx_fact[i])

            plt.plot(mu_steps_1, cmplx_fact_1, label=str(T) + 'K',
                     linewidth=3.0)
            plt.plot(mu_steps_2, cmplx_fact_2, linewidth=3.0)
            if output == 'average':
                plt.gca().get_lines()[1].set_c(plt.gca().get_lines()[0].get_c())
            elif output == 'tensor':
                plt.gca().get_lines()[3].set_c(plt.gca().get_lines()[0].get_c())
                plt.gca().get_lines()[4].set_c(plt.gca().get_lines()[1].get_c())
                plt.gca().get_lines()[5].set_c(plt.gca().get_lines()[2].get_c())

        plt.xlabel("E-E$_f$ (eV)", fontsize=30)
        plt.ylabel("Complexity Factor", fontsize=30)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        if output == 'tensor':
            plt.legend([str(i) + '_' + str(T) + 'K' for T in temps for i in
                        ('x', 'y', 'z')], fontsize=20)
        elif output == 'average':
            plt.legend(fontsize=20)
        plt.tight_layout()
        return plt

    def plot_seebeck_mu(self, temp=600, output='eig', xlim=None):
        """
        Plot the seebeck coefficient in function of Fermi level

        Args:
            temp:
                the temperature
            xlim:
                a list of min and max fermi energy by default (0, and band gap)
        Returns:
            a matplotlib object
        """
        import matplotlib.pyplot as plt
        plt.figure(figsize=(9, 7))
        seebeck = self._bz.get_seebeck(output=output, doping_levels=False)[
            temp]
        plt.plot(self._bz.mu_steps, seebeck,
                 linewidth=3.0)

        self._plot_bg_limits()
        self._plot_doping(temp)
        if output == 'eig':
            plt.legend(['S$_1$', 'S$_2$', 'S$_3$'])
        if xlim is None:
            plt.xlim(-0.5, self._bz.gap + 0.5)
        else:
            plt.xlim(xlim[0], xlim[1])
        plt.ylabel("Seebeck \n coefficient  ($\\mu$V/K)", fontsize=30.0)
        plt.xlabel("E-E$_f$ (eV)", fontsize=30)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.tight_layout()
        return plt

    def plot_conductivity_mu(self, temp=600, output='eig',
                             relaxation_time=1e-14, xlim=None):
        """
        Plot the conductivity in function of Fermi level. Semi-log plot

        Args:
            temp: the temperature
            xlim: a list of min and max fermi energy by default (0, and band
                gap)
            tau: A relaxation time in s. By default none and the plot is by
               units of relaxation time

        Returns:
            a matplotlib object
        """
        import matplotlib.pyplot as plt
        cond = self._bz.get_conductivity(relaxation_time=relaxation_time,
                                         output=output, doping_levels=False)[
            temp]
        plt.figure(figsize=(9, 7))
        plt.semilogy(self._bz.mu_steps, cond, linewidth=3.0)
        self._plot_bg_limits()
        self._plot_doping(temp)
        if output == 'eig':
            plt.legend(['$\\Sigma_1$', '$\\Sigma_2$', '$\\Sigma_3$'])
        if xlim is None:
            plt.xlim(-0.5, self._bz.gap + 0.5)
        else:
            plt.xlim(xlim)
        plt.ylim([1e13 * relaxation_time, 1e20 * relaxation_time])
        plt.ylabel("conductivity,\n $\\Sigma$ (1/($\\Omega$ m))", fontsize=30.0)
        plt.xlabel("E-E$_f$ (eV)", fontsize=30.0)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.tight_layout()
        return plt

    def plot_power_factor_mu(self, temp=600, output='eig',
                             relaxation_time=1e-14, xlim=None):
        """
        Plot the power factor in function of Fermi level. Semi-log plot

        Args:
            temp: the temperature
            xlim: a list of min and max fermi energy by default (0, and band
                gap)
            tau: A relaxation time in s. By default none and the plot is by
               units of relaxation time

        Returns:
            a matplotlib object
        """
        import matplotlib.pyplot as plt
        plt.figure(figsize=(9, 7))
        pf = self._bz.get_power_factor(relaxation_time=relaxation_time,
                                       output=output, doping_levels=False)[
            temp]
        plt.semilogy(self._bz.mu_steps, pf, linewidth=3.0)
        self._plot_bg_limits()
        self._plot_doping(temp)
        if output == 'eig':
            plt.legend(['PF$_1$', 'PF$_2$', 'PF$_3$'])
        if xlim is None:
            plt.xlim(-0.5, self._bz.gap + 0.5)
        else:
            plt.xlim(xlim)
        plt.ylabel("Power factor, ($\\mu$W/(mK$^2$))", fontsize=30.0)
        plt.xlabel("E-E$_f$ (eV)", fontsize=30.0)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.tight_layout()
        return plt

    def plot_zt_mu(self, temp=600, output='eig', relaxation_time=1e-14,
                   xlim=None):
        """
        Plot the ZT in function of Fermi level.

        Args:
            temp: the temperature
            xlim: a list of min and max fermi energy by default (0, and band
                gap)
            tau: A relaxation time in s. By default none and the plot is by
               units of relaxation time

        Returns:
            a matplotlib object
        """
        import matplotlib.pyplot as plt
        plt.figure(figsize=(9, 7))
        zt = self._bz.get_zt(relaxation_time=relaxation_time, output=output,
                             doping_levels=False)[temp]
        plt.plot(self._bz.mu_steps, zt, linewidth=3.0)
        self._plot_bg_limits()
        self._plot_doping(temp)
        if output == 'eig':
            plt.legend(['ZT$_1$', 'ZT$_2$', 'ZT$_3$'])
        if xlim is None:
            plt.xlim(-0.5, self._bz.gap + 0.5)
        else:
            plt.xlim(xlim)
        plt.ylabel("ZT", fontsize=30.0)
        plt.xlabel("E-E$_f$ (eV)", fontsize=30.0)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.tight_layout()
        return plt

    def plot_seebeck_temp(self, doping='all', output='average'):
        """
        Plot the Seebeck coefficient in function of temperature for different 
        doping levels.

        Args:
            dopings: the default 'all' plots all the doping levels in the analyzer.
                     Specify a list of doping levels if you want to plot only some.
            output: with 'average' you get an average of the three directions
                    with 'eigs' you get all the three directions.
        Returns:
            a matplotlib object
        """

        import matplotlib.pyplot as plt
        if output == 'average':
            sbk = self._bz.get_seebeck(output='average')
        elif output == 'eigs':
            sbk = self._bz.get_seebeck(output='eigs')

        plt.figure(figsize=(22, 14))
        tlist = sorted(sbk['n'].keys())
        doping = self._bz.doping['n'] if doping == 'all' else doping
        for i, dt in enumerate(['n', 'p']):
            plt.subplot(121 + i)
            for dop in doping:
                d = self._bz.doping[dt].index(dop)
                sbk_temp = []
                for temp in tlist:
                    sbk_temp.append(sbk[dt][temp][d])
                if output == 'average':
                    plt.plot(tlist, sbk_temp, marker='s',
                             label=str(dop) + ' $cm^{-3}$')
                elif output == 'eigs':
                    for xyz in range(3):
                        plt.plot(tlist, zip(*sbk_temp)[xyz], marker='s',
                                 label=str(xyz) + ' ' + str(dop) + ' $cm^{-3}$')
            plt.title(dt + '-type', fontsize=20)
            if i == 0:
                plt.ylabel("Seebeck \n coefficient  ($\\mu$V/K)", fontsize=30.0)
            plt.xlabel('Temperature (K)', fontsize=30.0)

            p = 'lower right' if i == 0 else ''
            plt.legend(loc=p, fontsize=15)
            plt.grid()
            plt.xticks(fontsize=25)
            plt.yticks(fontsize=25)

        plt.tight_layout()

        return plt

    def plot_conductivity_temp(self, doping='all', output='average',
                               relaxation_time=1e-14):
        """
        Plot the conductivity in function of temperature for different doping levels.

        Args:
            dopings: the default 'all' plots all the doping levels in the analyzer.
                     Specify a list of doping levels if you want to plot only some.
            output: with 'average' you get an average of the three directions
                    with 'eigs' you get all the three directions.
            relaxation_time: specify a constant relaxation time value
        
        Returns:
            a matplotlib object
        """
        import matplotlib.pyplot as plt
        import matplotlib.ticker as mtick

        if output == 'average':
            cond = self._bz.get_conductivity(relaxation_time=relaxation_time,
                                             output='average')
        elif output == 'eigs':
            cond = self._bz.get_conductivity(relaxation_time=relaxation_time,
                                             output='eigs')

        plt.figure(figsize=(22, 14))
        tlist = sorted(cond['n'].keys())
        doping = self._bz.doping['n'] if doping == 'all' else doping
        for i, dt in enumerate(['n', 'p']):
            plt.subplot(121 + i)
            for dop in doping:
                d = self._bz.doping[dt].index(dop)
                cond_temp = []
                for temp in tlist:
                    cond_temp.append(cond[dt][temp][d])
                if output == 'average':
                    plt.plot(tlist, cond_temp, marker='s',
                             label=str(dop) + ' $cm^{-3}$')
                elif output == 'eigs':
                    for xyz in range(3):
                        plt.plot(tlist, zip(*cond_temp)[xyz], marker='s',
                                 label=str(xyz) + ' ' + str(dop) + ' $cm^{-3}$')
            plt.title(dt + '-type', fontsize=20)
            if i == 0:
                plt.ylabel("conductivity $\\sigma$ (1/($\\Omega$ m))",
                           fontsize=30.0)
            plt.xlabel('Temperature (K)', fontsize=30.0)

            p = ''  # 'lower right' if i == 0 else ''
            plt.legend(loc=p, fontsize=15)
            plt.grid()
            plt.xticks(fontsize=25)
            plt.yticks(fontsize=25)
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

        plt.tight_layout()

        return plt

    def plot_power_factor_temp(self, doping='all', output='average',
                               relaxation_time=1e-14):
        """
        Plot the Power Factor in function of temperature for different doping levels.

        Args:
            dopings: the default 'all' plots all the doping levels in the analyzer.
                     Specify a list of doping levels if you want to plot only some.
            output: with 'average' you get an average of the three directions
                    with 'eigs' you get all the three directions.
            relaxation_time: specify a constant relaxation time value
        
        Returns:
            a matplotlib object
        """

        import matplotlib.pyplot as plt
        if output == 'average':
            pf = self._bz.get_power_factor(relaxation_time=relaxation_time,
                                           output='average')
        elif output == 'eigs':
            pf = self._bz.get_power_factor(relaxation_time=relaxation_time,
                                           output='eigs')

        plt.figure(figsize=(22, 14))
        tlist = sorted(pf['n'].keys())
        doping = self._bz.doping['n'] if doping == 'all' else doping
        for i, dt in enumerate(['n', 'p']):
            plt.subplot(121 + i)
            for dop in doping:
                d = self._bz.doping[dt].index(dop)
                pf_temp = []
                for temp in tlist:
                    pf_temp.append(pf[dt][temp][d])
                if output == 'average':
                    plt.plot(tlist, pf_temp, marker='s',
                             label=str(dop) + ' $cm^{-3}$')
                elif output == 'eigs':
                    for xyz in range(3):
                        plt.plot(tlist, zip(*pf_temp)[xyz], marker='s',
                                 label=str(xyz) + ' ' + str(dop) + ' $cm^{-3}$')
            plt.title(dt + '-type', fontsize=20)
            if i == 0:
                plt.ylabel("Power Factor ($\\mu$W/(mK$^2$))", fontsize=30.0)
            plt.xlabel('Temperature (K)', fontsize=30.0)

            p = ''  # 'lower right' if i == 0 else ''
            plt.legend(loc=p, fontsize=15)
            plt.grid()
            plt.xticks(fontsize=25)
            plt.yticks(fontsize=25)
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

        plt.tight_layout()
        return plt

    def plot_zt_temp(self, doping='all', output='average',
                     relaxation_time=1e-14):
        """
        Plot the figure of merit zT in function of temperature for different doping levels.

        Args:
            dopings: the default 'all' plots all the doping levels in the analyzer.
                     Specify a list of doping levels if you want to plot only some.
            output: with 'average' you get an average of the three directions
                    with 'eigs' you get all the three directions.
            relaxation_time: specify a constant relaxation time value
        
        Returns:
            a matplotlib object
        """

        import matplotlib.pyplot as plt
        if output == 'average':
            zt = self._bz.get_zt(relaxation_time=relaxation_time,
                                 output='average')
        elif output == 'eigs':
            zt = self._bz.get_zt(relaxation_time=relaxation_time, output='eigs')

        plt.figure(figsize=(22, 14))
        tlist = sorted(zt['n'].keys())
        doping = self._bz.doping['n'] if doping == 'all' else doping
        for i, dt in enumerate(['n', 'p']):
            plt.subplot(121 + i)
            for dop in doping:
                d = self._bz.doping[dt].index(dop)
                zt_temp = []
                for temp in tlist:
                    zt_temp.append(zt[dt][temp][d])
                if output == 'average':
                    plt.plot(tlist, zt_temp, marker='s',
                             label=str(dop) + ' $cm^{-3}$')
                elif output == 'eigs':
                    for xyz in range(3):
                        plt.plot(tlist, zip(*zt_temp)[xyz], marker='s',
                                 label=str(xyz) + ' ' + str(dop) + ' $cm^{-3}$')
            plt.title(dt + '-type', fontsize=20)
            if i == 0:
                plt.ylabel("zT", fontsize=30.0)
            plt.xlabel('Temperature (K)', fontsize=30.0)

            p = ''  # 'lower right' if i == 0 else ''
            plt.legend(loc=p, fontsize=15)
            plt.grid()
            plt.xticks(fontsize=25)
            plt.yticks(fontsize=25)

        plt.tight_layout()
        return plt

    def plot_eff_mass_temp(self, doping='all', output='average'):
        """
        Plot the average effective mass in function of temperature 
        for different doping levels.

        Args:
            dopings: the default 'all' plots all the doping levels in the analyzer.
                     Specify a list of doping levels if you want to plot only some.
            output: with 'average' you get an average of the three directions
                    with 'eigs' you get all the three directions.

        Returns:
            a matplotlib object
        """

        import matplotlib.pyplot as plt
        if output == 'average':
            em = self._bz.get_average_eff_mass(output='average')
        elif output == 'eigs':
            em = self._bz.get_average_eff_mass(output='eigs')

        plt.figure(figsize=(22, 14))
        tlist = sorted(em['n'].keys())
        doping = self._bz.doping['n'] if doping == 'all' else doping
        for i, dt in enumerate(['n', 'p']):
            plt.subplot(121 + i)
            for dop in doping:
                d = self._bz.doping[dt].index(dop)
                em_temp = []
                for temp in tlist:
                    em_temp.append(em[dt][temp][d])
                if output == 'average':
                    plt.plot(tlist, em_temp, marker='s',
                             label=str(dop) + ' $cm^{-3}$')
                elif output == 'eigs':
                    for xyz in range(3):
                        plt.plot(tlist, zip(*em_temp)[xyz], marker='s',
                                 label=str(xyz) + ' ' + str(dop) + ' $cm^{-3}$')
            plt.title(dt + '-type', fontsize=20)
            if i == 0:
                plt.ylabel("Effective mass (m$_e$)", fontsize=30.0)
            plt.xlabel('Temperature (K)', fontsize=30.0)

            p = ''  # 'lower right' if i == 0 else ''
            plt.legend(loc=p, fontsize=15)
            plt.grid()
            plt.xticks(fontsize=25)
            plt.yticks(fontsize=25)

        plt.tight_layout()
        return plt

    def plot_seebeck_dop(self, temps='all', output='average'):
        """
        Plot the Seebeck in function of doping levels for different temperatures.

        Args:
            temps: the default 'all' plots all the temperatures in the analyzer.
                   Specify a list of temperatures if you want to plot only some.
            output: with 'average' you get an average of the three directions
                    with 'eigs' you get all the three directions.

        Returns:
            a matplotlib object
        """

        import matplotlib.pyplot as plt
        if output == 'average':
            sbk = self._bz.get_seebeck(output='average')
        elif output == 'eigs':
            sbk = self._bz.get_seebeck(output='eigs')

        tlist = sorted(sbk['n'].keys()) if temps == 'all' else temps
        plt.figure(figsize=(22, 14))
        for i, dt in enumerate(['n', 'p']):
            plt.subplot(121 + i)
            for temp in tlist:
                if output == 'eigs':
                    for xyz in range(3):
                        plt.semilogx(self._bz.doping[dt],
                                     zip(*sbk[dt][temp])[xyz],
                                     marker='s',
                                     label=str(xyz) + ' ' + str(temp) + ' K')
                elif output == 'average':
                    plt.semilogx(self._bz.doping[dt], sbk[dt][temp],
                                 marker='s', label=str(temp) + ' K')
            plt.title(dt + '-type', fontsize=20)
            if i == 0:
                plt.ylabel("Seebeck coefficient ($\\mu$V/K)", fontsize=30.0)
            plt.xlabel('Doping concentration (cm$^{-3}$)', fontsize=30.0)

            p = 'lower right' if i == 0 else ''
            plt.legend(loc=p, fontsize=15)
            plt.grid()
            plt.xticks(fontsize=25)
            plt.yticks(fontsize=25)

        plt.tight_layout()

        return plt

    def plot_conductivity_dop(self, temps='all', output='average',
                              relaxation_time=1e-14):
        """
        Plot the conductivity in function of doping levels for different
        temperatures.

        Args:
            temps: the default 'all' plots all the temperatures in the analyzer.
                   Specify a list of temperatures if you want to plot only some.
            output: with 'average' you get an average of the three directions
                    with 'eigs' you get all the three directions.
            relaxation_time: specify a constant relaxation time value
        
        Returns:
            a matplotlib object
        """
        import matplotlib.pyplot as plt
        if output == 'average':
            cond = self._bz.get_conductivity(relaxation_time=relaxation_time,
                                             output='average')
        elif output == 'eigs':
            cond = self._bz.get_conductivity(relaxation_time=relaxation_time,
                                             output='eigs')

        tlist = sorted(cond['n'].keys()) if temps == 'all' else temps
        plt.figure(figsize=(22, 14))
        for i, dt in enumerate(['n', 'p']):
            plt.subplot(121 + i)
            for temp in tlist:
                if output == 'eigs':
                    for xyz in range(3):
                        plt.semilogx(self._bz.doping[dt],
                                     zip(*cond[dt][temp])[xyz],
                                     marker='s',
                                     label=str(xyz) + ' ' + str(temp) + ' K')
                elif output == 'average':
                    plt.semilogx(self._bz.doping[dt], cond[dt][temp],
                                 marker='s', label=str(temp) + ' K')
            plt.title(dt + '-type', fontsize=20)
            if i == 0:
                plt.ylabel("conductivity $\\sigma$ (1/($\\Omega$ m))",
                           fontsize=30.0)
            plt.xlabel('Doping concentration ($cm^{-3}$)', fontsize=30.0)
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            plt.legend(fontsize=15)
            plt.grid()
            plt.xticks(fontsize=25)
            plt.yticks(fontsize=25)

        plt.tight_layout()

        return plt

    def plot_power_factor_dop(self, temps='all', output='average',
                              relaxation_time=1e-14):
        """
        Plot the Power Factor in function of doping levels for different temperatures.

        Args:
            temps: the default 'all' plots all the temperatures in the analyzer.
                   Specify a list of temperatures if you want to plot only some.
            output: with 'average' you get an average of the three directions
                    with 'eigs' you get all the three directions.
            relaxation_time: specify a constant relaxation time value
        
        Returns:
            a matplotlib object
        """
        import matplotlib.pyplot as plt
        if output == 'average':
            pf = self._bz.get_power_factor(relaxation_time=relaxation_time,
                                           output='average')
        elif output == 'eigs':
            pf = self._bz.get_power_factor(relaxation_time=relaxation_time,
                                           output='eigs')

        tlist = sorted(pf['n'].keys()) if temps == 'all' else temps
        plt.figure(figsize=(22, 14))
        for i, dt in enumerate(['n', 'p']):
            plt.subplot(121 + i)
            for temp in tlist:
                if output == 'eigs':
                    for xyz in range(3):
                        plt.semilogx(self._bz.doping[dt],
                                     zip(*pf[dt][temp])[xyz],
                                     marker='s',
                                     label=str(xyz) + ' ' + str(temp) + ' K')
                elif output == 'average':
                    plt.semilogx(self._bz.doping[dt], pf[dt][temp],
                                 marker='s', label=str(temp) + ' K')
            plt.title(dt + '-type', fontsize=20)
            if i == 0:
                plt.ylabel("Power Factor  ($\\mu$W/(mK$^2$))", fontsize=30.0)
            plt.xlabel('Doping concentration ($cm^{-3}$)', fontsize=30.0)
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            p = ''  # 'lower right' if i == 0 else ''
            plt.legend(loc=p, fontsize=15)
            plt.grid()
            plt.xticks(fontsize=25)
            plt.yticks(fontsize=25)

        plt.tight_layout()

        return plt

    def plot_zt_dop(self, temps='all', output='average', relaxation_time=1e-14):
        """
        Plot the figure of merit zT in function of doping levels for different
        temperatures.

        Args:
            temps: the default 'all' plots all the temperatures in the analyzer.
                   Specify a list of temperatures if you want to plot only some.
            output: with 'average' you get an average of the three directions
                    with 'eigs' you get all the three directions.
            relaxation_time: specify a constant relaxation time value
        
        Returns:
            a matplotlib object
        """
        import matplotlib.pyplot as plt
        if output == 'average':
            zt = self._bz.get_zt(relaxation_time=relaxation_time,
                                 output='average')
        elif output == 'eigs':
            zt = self._bz.get_zt(relaxation_time=relaxation_time, output='eigs')

        tlist = sorted(zt['n'].keys()) if temps == 'all' else temps
        plt.figure(figsize=(22, 14))
        for i, dt in enumerate(['n', 'p']):
            plt.subplot(121 + i)
            for temp in tlist:
                if output == 'eigs':
                    for xyz in range(3):
                        plt.semilogx(self._bz.doping[dt],
                                     zip(*zt[dt][temp])[xyz],
                                     marker='s',
                                     label=str(xyz) + ' ' + str(temp) + ' K')
                elif output == 'average':
                    plt.semilogx(self._bz.doping[dt], zt[dt][temp],
                                 marker='s', label=str(temp) + ' K')
            plt.title(dt + '-type', fontsize=20)
            if i == 0:
                plt.ylabel("zT", fontsize=30.0)
            plt.xlabel('Doping concentration ($cm^{-3}$)', fontsize=30.0)

            p = 'lower right' if i == 0 else ''
            plt.legend(loc=p, fontsize=15)
            plt.grid()
            plt.xticks(fontsize=25)
            plt.yticks(fontsize=25)

        plt.tight_layout()

        return plt

    def plot_eff_mass_dop(self, temps='all', output='average'):
        """
        Plot the average effective mass in function of doping levels 
        for different temperatures.

        Args:
            temps: the default 'all' plots all the temperatures in the analyzer.
                   Specify a list of temperatures if you want to plot only some.
            output: with 'average' you get an average of the three directions
                    with 'eigs' you get all the three directions.
            relaxation_time: specify a constant relaxation time value
        
        Returns:
            a matplotlib object
        """

        import matplotlib.pyplot as plt
        if output == 'average':
            em = self._bz.get_average_eff_mass(output='average')
        elif output == 'eigs':
            em = self._bz.get_average_eff_mass(output='eigs')

        tlist = sorted(em['n'].keys()) if temps == 'all' else temps
        plt.figure(figsize=(22, 14))
        for i, dt in enumerate(['n', 'p']):
            plt.subplot(121 + i)
            for temp in tlist:
                if output == 'eigs':
                    for xyz in range(3):
                        plt.semilogx(self._bz.doping[dt],
                                     zip(*em[dt][temp])[xyz],
                                     marker='s',
                                     label=str(xyz) + ' ' + str(temp) + ' K')
                elif output == 'average':
                    plt.semilogx(self._bz.doping[dt], em[dt][temp],
                                 marker='s', label=str(temp) + ' K')
            plt.title(dt + '-type', fontsize=20)
            if i == 0:
                plt.ylabel("Effective mass (m$_e$)", fontsize=30.0)
            plt.xlabel('Doping concentration ($cm^{-3}$)', fontsize=30.0)

            p = 'lower right' if i == 0 else ''
            plt.legend(loc=p, fontsize=15)
            plt.grid()
            plt.xticks(fontsize=25)
            plt.yticks(fontsize=25)

        plt.tight_layout()

        return plt

    def plot_dos(self, sigma=0.05):
        """
        plot dos

        Args:
            sigma: a smearing

        Returns:
            a matplotlib object
        """
        plotter = DosPlotter(sigma=sigma)
        plotter.add_dos("t", self._bz.dos)
        return plotter.get_plot()

    def plot_carriers(self, temp=300):
        """
        Plot the carrier concentration in function of Fermi level

        Args:
            temp: the temperature

        Returns:
            a matplotlib object
        """
        import matplotlib.pyplot as plt
        plt.semilogy(self._bz.mu_steps,
                     abs(self._bz._carrier_conc[temp] / (self._bz.vol * 1e-24)),
                     linewidth=3.0, color='r')
        self._plot_bg_limits()
        self._plot_doping(temp)
        plt.xlim(-0.5, self._bz.gap + 0.5)
        plt.ylim(1e14, 1e22)
        plt.ylabel("carrier concentration (cm-3)", fontsize=30.0)
        plt.xlabel("E-E$_f$ (eV)", fontsize=30)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        return plt

    def plot_hall_carriers(self, temp=300):
        """
        Plot the Hall carrier concentration in function of Fermi level

        Args:
            temp: the temperature

        Returns:
            a matplotlib object
        """
        import matplotlib.pyplot as plt
        hall_carriers = [abs(i) for i in
                         self._bz.get_hall_carrier_concentration()[temp]]
        plt.semilogy(self._bz.mu_steps,
                     hall_carriers,
                     linewidth=3.0, color='r')
        self._plot_bg_limits()
        self._plot_doping(temp)
        plt.xlim(-0.5, self._bz.gap + 0.5)
        plt.ylim(1e14, 1e22)
        plt.ylabel("Hall carrier concentration (cm-3)", fontsize=30.0)
        plt.xlabel("E-E$_f$ (eV)", fontsize=30)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        return plt


class CohpPlotter:
    """
    Class for plotting crystal orbital Hamilton populations (COHPs) or
    crystal orbital overlap populations (COOPs). It is modeled after the
    DosPlotter object.

    Args/attributes:
        zero_at_efermi: Whether to shift all populations to have zero
            energy at the Fermi level. Defaults to True.

        are_coops: Switch to indicate that these are COOPs, not COHPs.
            Defaults to False for COHPs.

    """
    def __init__(self, zero_at_efermi=True, are_coops=False):
        self.zero_at_efermi = zero_at_efermi
        self.are_coops = are_coops
        self._cohps = OrderedDict()

    def add_cohp(self, label, cohp):
        """
        Adds a COHP for plotting.

        Args:
            label: Label for the COHP. Must be unique.

            cohp: COHP object.
        """
        energies = cohp.energies - cohp.efermi if self.zero_at_efermi \
            else cohp.energies
        populations = cohp.get_cohp()
        int_populations = cohp.get_icohp()
        self._cohps[label] = {"energies": energies, "COHP": populations,
                              "ICOHP": int_populations, "efermi": cohp.efermi}

    def add_cohp_dict(self, cohp_dict, key_sort_func=None):
        """
        Adds a dictionary of COHPs with an optional sorting function
        for the keys.

        Args:
            cohp_dict: dict of the form {label: Cohp}

            key_sort_func: function used to sort the cohp_dict keys.
        """
        if key_sort_func:
            keys = sorted(cohp_dict.keys(), key=key_sort_func)
        else:
            keys = cohp_dict.keys()
        for label in keys:
            self.add_cohp(label, cohp_dict[label])

    def get_cohp_dict(self):
        """
        Returns the added COHPs as a json-serializable dict. Note that if you
        have specified smearing for the COHP plot, the populations returned
        will be the smeared and not the original populations.

        Returns:
            dict: Dict of COHP data of the form {label: {"efermi": efermi,
            "energies": ..., "COHP": {Spin.up: ...}, "ICOHP": ...}}.
        """
        return jsanitize(self._cohps)

    def get_plot(self, xlim=None, ylim=None, plot_negative=None,
                 integrated=False, invert_axes=True):
        """
        Get a matplotlib plot showing the COHP.

        Args:
            xlim: Specifies the x-axis limits. Defaults to None for
                automatic determination.

            ylim: Specifies the y-axis limits. Defaults to None for
                automatic determination.

            plot_negative: It is common to plot -COHP(E) so that the
                sign means the same for COOPs and COHPs. Defaults to None
                for automatic determination: If are_coops is True, this
                will be set to False, else it will be set to True.

            integrated: Switch to plot ICOHPs. Defaults to False.

            invert_axes: Put the energies onto the y-axis, which is
                common in chemistry.

        Returns:
            A matplotlib object.
        """
        if self.are_coops:
            cohp_label = "COOP"
        else:
            cohp_label = "COHP"

        if plot_negative is None:
            plot_negative = True if not self.are_coops else False

        if integrated:
            cohp_label = "I" + cohp_label + " (eV)"

        if plot_negative:
            cohp_label = "-" + cohp_label

        if self.zero_at_efermi:
            energy_label = "$E - E_f$ (eV)"
        else:
            energy_label = "$E$ (eV)"

        ncolors = max(3, len(self._cohps))
        ncolors = min(9, ncolors)

        import palettable

        colors = palettable.colorbrewer.qualitative.Set1_9.mpl_colors

        plt = pretty_plot(12, 8)

        allpts = []
        keys = self._cohps.keys()
        for i, key in enumerate(keys):
            energies = self._cohps[key]["energies"]
            if not integrated:
                populations = self._cohps[key]["COHP"]
            else:
                populations = self._cohps[key]["ICOHP"]
            for spin in [Spin.up, Spin.down]:
                if spin in populations:
                    if invert_axes:
                        x = -populations[spin] if plot_negative \
                            else populations[spin]
                        y = energies
                    else:
                        x = energies
                        y = -populations[spin] if plot_negative \
                            else populations[spin]
                    allpts.extend(list(zip(x, y)))
                    if spin == Spin.up:
                        plt.plot(x, y, color=colors[i % ncolors],
                                 linestyle='-', label=str(key), linewidth=3)
                    else:
                        plt.plot(x, y, color=colors[i % ncolors],
                                 linestyle='--', linewidth=3)

        if xlim:
            plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)
        else:
            xlim = plt.xlim()
            relevanty = [p[1] for p in allpts if xlim[0] < p[0] < xlim[1]]
            plt.ylim((min(relevanty), max(relevanty)))

        xlim = plt.xlim()
        ylim = plt.ylim()
        if not invert_axes:
            plt.plot(xlim, [0, 0], "k-", linewidth=2)
            if self.zero_at_efermi:
                plt.plot([0, 0], ylim, "k--", linewidth=2)
            else:
                plt.plot([self._cohps[key]['efermi'],
                         self._cohps[key]['efermi']], ylim,
                         color=colors[i % ncolors],
                         linestyle='--', linewidth=2)
        else:
            plt.plot([0, 0], ylim, "k-", linewidth=2)
            if self.zero_at_efermi:
                plt.plot(xlim, [0, 0], "k--", linewidth=2)
            else:
                plt.plot(xlim, [self._cohps[key]['efermi'],
                         self._cohps[key]['efermi']],
                         color=colors[i % ncolors],
                         linestyle='--', linewidth=2)

        if invert_axes:
            plt.xlabel(cohp_label)
            plt.ylabel(energy_label)
        else:
            plt.xlabel(energy_label)
            plt.ylabel(cohp_label)

        plt.legend()
        leg = plt.gca().get_legend()
        ltext = leg.get_texts()
        plt.setp(ltext, fontsize=30)
        plt.tight_layout()
        return plt

    def save_plot(self, filename, img_format="eps", xlim=None, ylim=None):
        """
        Save matplotlib plot to a file.

        Args:
            filename: File name to write to.
            img_format: Image format to use. Defaults to EPS.
            xlim: Specifies the x-axis limits. Defaults to None for
                automatic determination.
            ylim: Specifies the y-axis limits. Defaults to None for
                automatic determination.
        """
        plt = self.get_plot(xlim, ylim)
        plt.savefig(filename, format=img_format)

    def show(self, xlim=None, ylim=None):
        """
        Show the plot using matplotlib.

        Args:
            xlim: Specifies the x-axis limits. Defaults to None for
                automatic determination.
            ylim: Specifies the y-axis limits. Defaults to None for
                automatic determination.
        """
        plt = self.get_plot(xlim, ylim)
        plt.show()


def plot_fermi_surface(data, structure, cbm, energy_levels=[],
                       multiple_figure=True,
                       mlab_figure=None, kpoints_dict={}, color=(0, 0, 1),
                       transparency_factor=[], labels_scale_factor=0.05,
                       points_scale_factor=0.02, interative=True):
    """
    Plot the Fermi surface at specific energy value.

    Args:
        data: energy values in a 3D grid from a CUBE file 
              via read_cube_file function, or from a 
              BoltztrapAnalyzer.fermi_surface_data
        structure: structure object of the material
        energy_levels: list of energy value of the fermi surface. 
                       By default 0 eV correspond to the VBM, as in 
                       the plot of band structure along symmetry line.
                       Default: max energy value + 0.01 eV
        cbm: Boolean value to specify if the considered band is 
             a conduction band or not
        multiple_figure: if True a figure for each energy level will be shown.
                         If False all the surfaces will be shown in the same figure.
                         In this las case, tune the transparency factor.
        mlab_figure: provide a previous figure to plot a new surface on it.
        kpoints_dict: dictionary of kpoints to show in the plot.
                      example: {"K":[0.5,0.0,0.5]}, 
                      where the coords are fractional.
        color: tuple (r,g,b) of integers to define the color of the surface.
        transparency_factor: list of values in the range [0,1] to tune
                             the opacity of the surfaces.
        labels_scale_factor: factor to tune the size of the kpoint labels
        points_scale_factor: factor to tune the size of the kpoint points
        interative: if True an interactive figure will be shown.
                    If False a non interactive figure will be shown, but
                    it is possible to plot other surfaces on the same figure.
                    To make it interactive, run mlab.show().
        
    Returns:
        a Mayavi figure and a mlab module to control the plot.

    Note: Experimental. 
          Please, double check the surface shown by using some 
          other software and report issues.
    """

    try:
        from mayavi import mlab
    except ImportError:
        raise BoltztrapError(
            "Mayavi package should be installed to use this function")

    bz = structure.lattice.reciprocal_lattice.get_wigner_seitz_cell()
    cell = structure.lattice.reciprocal_lattice.matrix

    fact = 1 if cbm == False else -1
    en_min = np.min(fact * data.ravel())
    en_max = np.max(fact * data.ravel())

    if energy_levels == []:
        energy_levels = [en_min + 0.01] if cbm == True else \
            [en_max - 0.01]
        print("Energy level set to: " + str(energy_levels[0]) + " eV")

    else:
        for e in energy_levels:
            if e > en_max or e < en_min:
                raise BoltztrapError("energy level " + str(e) +
                                     " not in the range of possible energies: [" +
                                     str(en_min) + ", " + str(en_max) + "]")

    if transparency_factor == []:
        transparency_factor = [1] * len(energy_levels)

    if mlab_figure:
        fig = mlab_figure

    if mlab_figure == None and not multiple_figure:
        fig = mlab.figure(size=(1024, 768), bgcolor=(1, 1, 1))
        for iface in range(len(bz)):
            for line in itertools.combinations(bz[iface], 2):
                for jface in range(len(bz)):
                    if iface < jface and any(np.all(line[0] == x)
                                             for x in bz[jface]) and \
                            any(np.all(line[1] == x)
                                for x in bz[jface]):
                        mlab.plot3d(*zip(line[0], line[1]), color=(0, 0, 0),
                                    tube_radius=None, figure=fig)
        for label, coords in kpoints_dict.iteritems():
            label_coords = structure.lattice.reciprocal_lattice \
                .get_cartesian_coords(coords)
            mlab.points3d(*label_coords, scale_factor=points_scale_factor,
                          color=(0, 0, 0), figure=fig)
            mlab.text3d(*label_coords, text=label, scale=labels_scale_factor,
                        color=(0, 0, 0), figure=fig)

    for isolevel, alpha in zip(energy_levels, transparency_factor):
        if multiple_figure:
            fig = mlab.figure(size=(1024, 768), bgcolor=(1, 1, 1))

            for iface in range(len(bz)):
                for line in itertools.combinations(bz[iface], 2):
                    for jface in range(len(bz)):
                        if iface < jface and any(np.all(line[0] == x)
                                                 for x in bz[jface]) and \
                                any(np.all(line[1] == x)
                                    for x in bz[jface]):
                            mlab.plot3d(*zip(line[0], line[1]), color=(0, 0, 0),
                                        tube_radius=None, figure=fig)

            for label, coords in kpoints_dict.iteritems():
                label_coords = structure.lattice.reciprocal_lattice \
                    .get_cartesian_coords(coords)
                mlab.points3d(*label_coords, scale_factor=points_scale_factor,
                              color=(0, 0, 0), figure=fig)
                mlab.text3d(*label_coords, text=label,
                            scale=labels_scale_factor, color=(0, 0, 0),
                            figure=fig)

        cp = mlab.contour3d(fact * data, contours=[isolevel], transparent=True,
                            colormap='hot', color=color, opacity=alpha,
                            figure=fig)

        polydata = cp.actor.actors[0].mapper.input
        pts = np.array(polydata.points)  # - 1
        polydata.points = np.dot(pts,
                                 cell / np.array(data.shape)[:, np.newaxis])

        cx, cy, cz = [np.mean(np.array(polydata.points)[:, i])
                      for i in range(3)]

        polydata.points = (np.array(polydata.points) - [cx, cy, cz]) * 2

        mlab.view(distance='auto')
    if interative == True:
        mlab.show()

    return fig, mlab


def plot_wigner_seitz(lattice, ax=None, **kwargs):
    """
    Adds the skeleton of the Wigner-Seitz cell of the lattice to a matplotlib Axes

    Args:
        lattice: Lattice object
        ax: matplotlib :class:`Axes` or None if a new figure should be created.
        kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to black
            and linewidth to 1.

    Returns:
        matplotlib figure and matplotlib ax
    """
    ax, fig, plt = get_ax3d_fig_plt(ax)

    if "color" not in kwargs:
        kwargs["color"] = "k"
    if "linewidth" not in kwargs:
        kwargs["linewidth"] = 1

    bz = lattice.get_wigner_seitz_cell()
    ax, fig, plt = get_ax3d_fig_plt(ax)
    for iface in range(len(bz)):
        for line in itertools.combinations(bz[iface], 2):
            for jface in range(len(bz)):
                if iface < jface and any(
                        np.all(line[0] == x) for x in bz[jface]) \
                        and any(np.all(line[1] == x) for x in bz[jface]):
                    ax.plot(*zip(line[0], line[1]), **kwargs)

    return fig, ax


def plot_lattice_vectors(lattice, ax=None, **kwargs):
    """
    Adds the basis vectors of the lattice provided to a matplotlib Axes

    Args:
        lattice: Lattice object
        ax: matplotlib :class:`Axes` or None if a new figure should be created.
        kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to green
            and linewidth to 3.

    Returns:
        matplotlib figure and matplotlib ax
    """
    ax, fig, plt = get_ax3d_fig_plt(ax)

    if "color" not in kwargs:
        kwargs["color"] = "g"
    if "linewidth" not in kwargs:
        kwargs["linewidth"] = 3

    vertex1 = lattice.get_cartesian_coords([0.0, 0.0, 0.0])
    vertex2 = lattice.get_cartesian_coords([1.0, 0.0, 0.0])
    ax.plot(*zip(vertex1, vertex2), **kwargs)
    vertex2 = lattice.get_cartesian_coords([0.0, 1.0, 0.0])
    ax.plot(*zip(vertex1, vertex2), **kwargs)
    vertex2 = lattice.get_cartesian_coords([0.0, 0.0, 1.0])
    ax.plot(*zip(vertex1, vertex2), **kwargs)

    return fig, ax


def plot_path(line, lattice=None, coords_are_cartesian=False, ax=None,
              **kwargs):
    """
    Adds a line passing through the coordinates listed in 'line' to a matplotlib Axes

    Args:
        line: list of coordinates.
        lattice: Lattice object used to convert from reciprocal to cartesian coordinates
        coords_are_cartesian: Set to True if you are providing
            coordinates in cartesian coordinates. Defaults to False.
            Requires lattice if False.
        ax: matplotlib :class:`Axes` or None if a new figure should be created.
        kwargs: kwargs passed to the matplotlib function 'plot'. Color defaults to red
            and linewidth to 3.

    Returns:
        matplotlib figure and matplotlib ax
    """

    ax, fig, plt = get_ax3d_fig_plt(ax)

    if "color" not in kwargs:
        kwargs["color"] = "r"
    if "linewidth" not in kwargs:
        kwargs["linewidth"] = 3

    for k in range(1, len(line)):
        vertex1 = line[k - 1]
        vertex2 = line[k]
        if not coords_are_cartesian:
            if lattice is None:
                raise ValueError(
                    "coords_are_cartesian False requires the lattice")
            vertex1 = lattice.get_cartesian_coords(vertex1)
            vertex2 = lattice.get_cartesian_coords(vertex2)
        ax.plot(*zip(vertex1, vertex2), **kwargs)

    return fig, ax


def plot_labels(labels, lattice=None, coords_are_cartesian=False, ax=None,
                **kwargs):
    """
    Adds labels to a matplotlib Axes

    Args:
        labels: dict containing the label as a key and the coordinates as value.
        lattice: Lattice object used to convert from reciprocal to cartesian coordinates
        coords_are_cartesian: Set to True if you are providing.
            coordinates in cartesian coordinates. Defaults to False.
            Requires lattice if False.
        ax: matplotlib :class:`Axes` or None if a new figure should be created.
        kwargs: kwargs passed to the matplotlib function 'text'. Color defaults to blue
            and size to 25.

    Returns:
        matplotlib figure and matplotlib ax
    """
    ax, fig, plt = get_ax3d_fig_plt(ax)

    if "color" not in kwargs:
        kwargs["color"] = "b"
    if "size" not in kwargs:
        kwargs["size"] = 25

    for k, coords in labels.items():
        label = k
        if k.startswith("\\") or k.find("_") != -1:
            label = "$" + k + "$"
        off = 0.01
        if coords_are_cartesian:
            coords = np.array(coords)
        else:
            if lattice is None:
                raise ValueError(
                    "coords_are_cartesian False requires the lattice")
            coords = lattice.get_cartesian_coords(coords)
        ax.text(*(coords + off), s=label, **kwargs)

    return fig, ax


def fold_point(p, lattice, coords_are_cartesian=False):
    """
    Folds a point with coordinates p inside the first Brillouin zone of the lattice.

    Args:
        p: coordinates of one point
        lattice: Lattice object used to convert from reciprocal to cartesian coordinates
        coords_are_cartesian: Set to True if you are providing
            coordinates in cartesian coordinates. Defaults to False.

    Returns:
        The cartesian coordinates folded inside the first Brillouin zone
    """

    if coords_are_cartesian:
        p = lattice.get_fractional_coords(p)
    else:
        p = np.array(p)

    p = np.mod(p + 0.5 - 1e-10, 1) - 0.5 + 1e-10
    p = lattice.get_cartesian_coords(p)

    closest_lattice_point = None
    smallest_distance = 10000
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                lattice_point = np.dot((i, j, k), lattice.matrix)
                dist = np.linalg.norm(p - lattice_point)
                if closest_lattice_point is None or dist < smallest_distance:
                    closest_lattice_point = lattice_point
                    smallest_distance = dist

    if not np.allclose(closest_lattice_point, (0, 0, 0)):
        p = p - closest_lattice_point

    return p


def plot_points(points, lattice=None, coords_are_cartesian=False, fold=False,
                ax=None, **kwargs):
    """
    Adds Points to a matplotlib Axes

    Args:
        points: list of coordinates
        lattice: Lattice object used to convert from reciprocal to cartesian coordinates
        coords_are_cartesian: Set to True if you are providing
            coordinates in cartesian coordinates. Defaults to False.
            Requires lattice if False.
        fold: whether the points should be folded inside the first Brillouin Zone.
            Defaults to False. Requires lattice if True.
        ax: matplotlib :class:`Axes` or None if a new figure should be created.
        kwargs: kwargs passed to the matplotlib function 'scatter'. Color defaults to blue

    Returns:
        matplotlib figure and matplotlib ax
    """
    ax, fig, plt = get_ax3d_fig_plt(ax)

    if "color" not in kwargs:
        kwargs["color"] = "b"

    if (not coords_are_cartesian or fold) and lattice is None:
        raise ValueError(
            "coords_are_cartesian False or fold True require the lattice")

    for p in points:

        if fold:
            p = fold_point(p, lattice,
                           coords_are_cartesian=coords_are_cartesian)

        elif not coords_are_cartesian:
            p = lattice.get_cartesian_coords(p)

        ax.scatter(*p, **kwargs)

    return fig, ax


@add_fig_kwargs
def plot_brillouin_zone_from_kpath(kpath, ax=None, **kwargs):
    """
    Gives the plot (as a matplotlib object) of the symmetry line path in
        the Brillouin Zone.

    Args:
        kpath (HighSymmKpath): a HighSymmKPath object
        ax: matplotlib :class:`Axes` or None if a new figure should be created.
        **kwargs: provided by add_fig_kwargs decorator

    Returns:
        matplotlib figure

    """
    lines = [[kpath.kpath['kpoints'][k] for k in p]
             for p in kpath.kpath['path']]
    return plot_brillouin_zone(bz_lattice=kpath.prim_rec, lines=lines, ax=ax,
                               labels=kpath.kpath['kpoints'], **kwargs)


@add_fig_kwargs
def plot_brillouin_zone(bz_lattice, lines=None, labels=None, kpoints=None,
                        fold=False, coords_are_cartesian=False,
                        ax=None, **kwargs):
    """
    Plots a 3D representation of the Brillouin zone of the structure.
    Can add to the plot paths, labels and kpoints

    Args:
        bz_lattice: Lattice object of the Brillouin zone
        lines: list of lists of coordinates. Each list represent a different path
        labels: dict containing the label as a key and the coordinates as value.
        kpoints: list of coordinates
        fold: whether the points should be folded inside the first Brillouin Zone.
            Defaults to False. Requires lattice if True.
        coords_are_cartesian: Set to True if you are providing
            coordinates in cartesian coordinates. Defaults to False.
        ax: matplotlib :class:`Axes` or None if a new figure should be created.
        kwargs: provided by add_fig_kwargs decorator

    Returns:
        matplotlib figure
    """

    fig, ax = plot_lattice_vectors(bz_lattice, ax=ax)
    plot_wigner_seitz(bz_lattice, ax=ax)
    if lines is not None:
        for line in lines:
            plot_path(line, bz_lattice,
                      coords_are_cartesian=coords_are_cartesian, ax=ax)

    if labels is not None:
        plot_labels(labels, bz_lattice,
                    coords_are_cartesian=coords_are_cartesian, ax=ax)
        plot_points(labels.values(), bz_lattice,
                    coords_are_cartesian=coords_are_cartesian,
                    fold=False, ax=ax)

    if kpoints is not None:
        plot_points(kpoints, bz_lattice,
                    coords_are_cartesian=coords_are_cartesian,
                    ax=ax, fold=fold)

    ax.set_xlim3d(-1, 1)
    ax.set_ylim3d(-1, 1)
    ax.set_zlim3d(-1, 1)

    ax.set_aspect('equal')
    ax.axis("off")

    return fig


def plot_ellipsoid(hessian, center, lattice=None, rescale=1.0, ax=None,
                   coords_are_cartesian=False, arrows=False, **kwargs):
    """
    Plots a 3D ellipsoid rappresenting the Hessian matrix in input.
    Useful to get a graphical visualization of the effective mass
    of a band in a single k-point.
    
    Args:
        hessian: the Hessian matrix
        center: the center of the ellipsoid in reciprocal coords (Default)
        lattice: Lattice object of the Brillouin zone
        rescale: factor for size scaling of the ellipsoid
        ax: matplotlib :class:`Axes` or None if a new figure should be created.
        coords_are_cartesian: Set to True if you are providing a center in
                              cartesian coordinates. Defaults to False.
        kwargs: kwargs passed to the matplotlib function 'plot_wireframe'. 
                Color defaults to blue, rstride and cstride
                default to 4, alpha defaults to 0.2.
    Returns:
        matplotlib figure and matplotlib ax
    Example of use:
        fig,ax=plot_wigner_seitz(struct.reciprocal_lattice)
        plot_ellipsoid(hessian,[0.0,0.0,0.0], struct.reciprocal_lattice,ax=ax)
    """

    if (not coords_are_cartesian) and lattice is None:
        raise ValueError(
            "coords_are_cartesian False or fold True require the lattice")

    if not coords_are_cartesian:
        center = lattice.get_cartesian_coords(center)

    if "color" not in kwargs:
        kwargs["color"] = "b"
    if "rstride" not in kwargs:
        kwargs["rstride"] = 4
    if "cstride" not in kwargs:
        kwargs["cstride"] = 4
    if "alpha" not in kwargs:
        kwargs["alpha"] = 0.2

    # calculate the ellipsoid
    # find the rotation matrix and radii of the axes
    U, s, rotation = np.linalg.svd(hessian)
    radii = 1.0 / np.sqrt(s)

    # from polar coordinates
    u = np.linspace(0.0, 2.0 * np.pi, 100)
    v = np.linspace(0.0, np.pi, 100)
    x = radii[0] * np.outer(np.cos(u), np.sin(v))
    y = radii[1] * np.outer(np.sin(u), np.sin(v))
    z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
    for i in range(len(x)):
        for j in range(len(x)):
            [x[i, j], y[i, j], z[i, j]] = np.dot([x[i, j], y[i, j], z[i, j]],
                                                 rotation) * rescale + center

    # add the ellipsoid to the current axes
    ax, fig, plt = get_ax3d_fig_plt(ax)
    ax.plot_wireframe(x, y, z, **kwargs)

    if arrows:
        color = ('b', 'g', 'r')
        em = np.zeros((3, 3))
        for i in range(3):
            em[i, :] = rotation[i, :] / np.linalg.norm(rotation[i, :])
        for i in range(3):
            ax.quiver3D(center[0], center[1], center[2], em[i, 0], em[i, 1],
                        em[i, 2], pivot='tail',
                        arrow_length_ratio=0.2, length=radii[i] * rescale,
                        color=color[i])

    return fig, ax
