# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

"""
This module implements plotter for DOS and band structure.
"""

__author__ = "Shyue Ping Ong, Geoffroy Hautier"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "May 1, 2012"

import logging
import math
import itertools
from collections import OrderedDict

import numpy as np

from monty.json import jsanitize
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine

logger = logging.getLogger('BSPlotter')


class DosPlotter(object):
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
            Dict of dos data. Generally of the form, {label: {'energies':..,
            'densities': {'up':...}, 'efermi':efermi}}
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
        from pymatgen.util.plotting_utils import get_publication_quality_plot
        ncolors = max(3, len(self._doses))
        ncolors = min(9, ncolors)
        colors = brewer2mpl.get_map('Set1', 'qualitative', ncolors).mpl_colors

        y = None
        alldensities = []
        allenergies = []
        plt = get_publication_quality_plot(12, 8)

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
                ppl.plot(x, y, color=colors[i % ncolors],
                         label=str(key), linewidth=3)
            if not self.zero_at_efermi:
                ylim = plt.ylim()
                ppl.plot([self._doses[key]['efermi'],
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


class BSPlotter(object):
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
        self._nb_bands = self._bs._nb_bands

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
            A dict of the following format:
            ticks: A dict with the 'distances' at which there is a kpoint (the
            x axis) and the labels (None if no label)
            energy: A dict storing bands for spin up and spin down data
            [{Spin:[band_index][k_point_index]}] as a list (one element
            for each branch) of energy for each kpoint. The data is
            stored by branch to facilitate the plotting
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

        for b in self._bs._branches:

            if self._bs.is_spin_polarized:
                energy.append({str(Spin.up): [], str(Spin.down): []})
            else:
                energy.append({str(Spin.up): []})
            distance.append([self._bs._distance[j]
                             for j in range(b['start_index'],
                                            b['end_index'] + 1)])
            ticks = self.get_ticks()

            for i in range(self._nb_bands):
                energy[-1][str(Spin.up)].append(
                    [self._bs._bands[Spin.up][i][j] - zero_energy
                     for j in range(b['start_index'], b['end_index'] + 1)])
            if self._bs.is_spin_polarized:
                for i in range(self._nb_bands):
                    energy[-1][str(Spin.down)].append(
                        [self._bs._bands[Spin.down][i][j] - zero_energy
                         for j in range(b['start_index'], b['end_index'] + 1)])

        vbm = self._bs.get_vbm()
        cbm = self._bs.get_cbm()

        vbm_plot = []
        cbm_plot = []

        for index in cbm['kpoint_index']:
            cbm_plot.append((self._bs._distance[index],
                             cbm['energy'] - zero_energy if zero_to_efermi
                             else cbm['energy']))

        for index in vbm['kpoint_index']:
            vbm_plot.append((self._bs._distance[index],
                             vbm['energy'] - zero_energy if zero_to_efermi
                             else vbm['energy']))

        bg = self._bs.get_band_gap()
        direct = "Indirect"
        if bg['direct']:
            direct = "Direct"

        return {'ticks': ticks, 'distances': distance, 'energy': energy,
                'vbm': vbm_plot, 'cbm': cbm_plot,
                'lattice': self._bs._lattice_rec.as_dict(),
                'zero_energy': zero_energy, 'is_metal': self._bs.is_metal(),
                'band_gap': "{} {} bandgap = {}".format(direct,
                                                        bg['transition'],
                                                        bg['energy'])
                if not self._bs.is_metal() else ""}

    def get_plot(self, zero_to_efermi=True, ylim=None, smooth=False,
                 vbm_cbm_marker=False):
        """
        get a matplotlib object for the bandstructure plot.
        Blue lines are up spin, red lines are down
        spin.

        Args:
            zero_to_efermi: Automatically subtract off the Fermi energy from
                the eigenvalues and plot (E-Ef).
            ylim: Specify the y-axis (energy) limits; by default None let
                the code choose. It is vbm-4 and cbm+4 if insulator
                efermi-10 and efermi+10 if metal
            smooth: interpolates the bands by a spline cubic
        """
        from pymatgen.util.plotting_utils import get_publication_quality_plot
        plt = get_publication_quality_plot(12, 8)
        from matplotlib import rc
        import scipy.interpolate as scint
        rc('text', usetex=True)

        # main internal config options
        e_min = -4
        e_max = 4
        if self._bs.is_metal():
            e_min = -10
            e_max = 10
        #band_linewidth = 3
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
            for d in range(len(data['distances'])):
                for i in range(self._nb_bands):
                    tck = scint.splrep(
                        data['distances'][d],
                        [data['energy'][d][str(Spin.up)][i][j]
                         for j in range(len(data['distances'][d]))])
                    step = (data['distances'][d][-1]
                            - data['distances'][d][0]) / 1000

                    plt.plot([x * step + data['distances'][d][0]
                              for x in range(1000)],
                             [scint.splev(x * step + data['distances'][d][0],
                                          tck, der=0)
                              for x in range(1000)], 'b-',
                             linewidth=band_linewidth)

                    if self._bs.is_spin_polarized:
                        tck = scint.splrep(
                            data['distances'][d],
                            [data['energy'][d][str(Spin.down)][i][j]
                             for j in range(len(data['distances'][d]))])
                        step = (data['distances'][d][-1]
                                - data['distances'][d][0]) / 1000

                        plt.plot([x * step + data['distances'][d][0]
                                  for x in range(1000)],
                                 [scint.splev(
                                     x * step + data['distances'][d][0],
                                     tck, der=0)
                                  for x in range(1000)], 'r--',
                                 linewidth=band_linewidth)
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
                    plt.ylim(self._bs.efermi + e_min, self._bs._efermi + e_max)
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
           
        plt.tight_layout()

        return plt

    def show(self, zero_to_efermi=True, ylim=None, smooth=False):
        """
        Show the plot using matplotlib.

        Args:
            zero_to_efermi: Automatically subtract off the Fermi energy from
                the eigenvalues and plot (E-Ef).
            ylim: Specify the y-axis (energy) limits; by default None let
                the code choose. It is vbm-4 and cbm+4 if insulator
                efermi-10 and efermi+10 if metal
            smooth: interpolates the bands by a spline cubic
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
            A dict with 'distance': a list of distance at which ticks should
            be set and 'label': a list of label for each of those ticks.
        """
        tick_distance = []
        tick_labels = []
        previous_label = self._bs._kpoints[0].label
        previous_branch = self._bs._branches[0]['name']
        for i, c in enumerate(self._bs._kpoints):
            if c.label is not None:
                tick_distance.append(self._bs._distance[i])
                this_branch = None
                for b in self._bs._branches:
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
        plt = self.get_plot()
        data_orig = self.bs_plot_data()
        data = other_plotter.bs_plot_data()
        band_linewidth = 3
        for i in range(other_plotter._nb_bands):
            plt.plot(data_orig['distances'],
                     [e for e in data['energy'][str(Spin.up)][i]],
                     'r-', linewidth=band_linewidth)
            if other_plotter._bs.is_spin_polarized:
                plt.plot(data_orig['distances'],
                         [e for e in data['energy'][str(Spin.down)][i]],
                         'r-', linewidth=band_linewidth)
        return plt

    def plot_brillouin(self):
        """
        plot the Brillouin zone
        """
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        mpl.rcParams['legend.fontsize'] = 10

        fig = plt.figure()
        ax = Axes3D(fig)
        vec1 = self._bs.lattice.matrix[0]
        vec2 = self._bs.lattice.matrix[1]
        vec3 = self._bs.lattice.matrix[2]

        # make the grid
        max_x = -1000
        max_y = -1000
        max_z = -1000
        min_x = 1000
        min_y = 1000
        min_z = 1000
        list_k_points = []
        for i in [-1, 0, 1]:
            for j in [-1, 0, 1]:
                for k in [-1, 0, 1]:
                    list_k_points.append(i * vec1 + j * vec2 + k * vec3)
                    if list_k_points[-1][0] > max_x:
                        max_x = list_k_points[-1][0]
                    if list_k_points[-1][1] > max_y:
                        max_y = list_k_points[-1][1]
                    if list_k_points[-1][2] > max_z:
                        max_z = list_k_points[-1][0]

                    if list_k_points[-1][0] < min_x:
                        min_x = list_k_points[-1][0]
                    if list_k_points[-1][1] < min_y:
                        min_y = list_k_points[-1][1]
                    if list_k_points[-1][2] < min_z:
                        min_z = list_k_points[-1][0]

        vertex = _qvertex_target(list_k_points, 13)
        lines = get_lines_voronoi(vertex)

        for i in range(len(lines)):
            vertex1 = lines[i]['start']
            vertex2 = lines[i]['end']
            ax.plot([vertex1[0], vertex2[0]], [vertex1[1], vertex2[1]],
                    [vertex1[2], vertex2[2]], color='k')

        for b in self._bs._branches:
            vertex1 = self._bs.kpoints[b['start_index']].cart_coords
            vertex2 = self._bs.kpoints[b['end_index']].cart_coords
            ax.plot([vertex1[0], vertex2[0]], [vertex1[1], vertex2[1]],
                    [vertex1[2], vertex2[2]], color='r', linewidth=3)

        for k in self._bs.kpoints:
            if k.label:
                label = k.label
                if k.label.startswith("\\") or k.label.find("_") != -1:
                    label = "$" + k.label + "$"
                off = 0.01
                ax.text(k.cart_coords[0] + off, k.cart_coords[1] + off,
                        k.cart_coords[2] + off, label, color='b', size='25')
                ax.scatter([k.cart_coords[0]], [k.cart_coords[1]],
                           [k.cart_coords[2]], color='b')

        # make ticklabels and ticklines invisible
        for a in ax.w_xaxis.get_ticklines() + ax.w_xaxis.get_ticklabels():
            a.set_visible(False)
        for a in ax.w_yaxis.get_ticklines() + ax.w_yaxis.get_ticklabels():
            a.set_visible(False)
        for a in ax.w_zaxis.get_ticklines() + ax.w_zaxis.get_ticklabels():
            a.set_visible(False)

        ax.grid(False)

        plt.show()
        ax.axis("off")


class BSPlotterProjected(BSPlotter):
    """
    Class to plot or get data to facilitate the plot of band structure objects
    projected along orbitals, elements or sites.

    Args:
        bs: A BandStructureSymmLine object with projections.
    """

    def __init__(self, bs):
        if len(bs._projections) == 0:
            raise ValueError("try to plot projections"
                             " on a band structure without any")
        super(BSPlotterProjected, self).__init__(bs)

    def _get_projections_by_branches(self, dictio):
        proj = self._bs.get_projections_on_elts_and_orbitals(dictio)
        proj_br = []
        print(len(proj[Spin.up]))
        print(len(proj[Spin.up][0]))
        for c in proj[Spin.up][0]:
            print(c)
        for b in self._bs._branches:
            print(b)
            if self._bs.is_spin_polarized:
                proj_br.append(
                    {str(Spin.up): [[] for l in range(self._nb_bands)],
                     str(Spin.down): [[] for l in range(self._nb_bands)]})
            else:
                proj_br.append(
                    {str(Spin.up): [[] for l in range(self._nb_bands)]})
            print((len(proj_br[-1][str(Spin.up)]), self._nb_bands))

            for i in range(self._nb_bands):
                for j in range(b['start_index'], b['end_index'] + 1):
                    proj_br[-1][str(Spin.up)][i].append(
                        {e: {o: proj[Spin.up][i][j][e][o]
                             for o in proj[Spin.up][i][j][e]}
                         for e in proj[Spin.up][i][j]})
            if self._bs.is_spin_polarized:
                for b in self._bs._branches:
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
        from pymatgen.util.plotting_utils import get_publication_quality_plot
        band_linewidth = 1.0
        fig_number = sum([len(v) for v in dictio.values()])
        proj = self._get_projections_by_branches(dictio)
        data = self.bs_plot_data(zero_to_efermi)
        plt = get_publication_quality_plot(12, 8)
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
                            plt.ylim(self._bs.efermi + e_min, self._bs._efermi
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
                                                  self._bs._structure.composition.elements})
        data = self.bs_plot_data(zero_to_efermi)
        from pymatgen.util.plotting_utils import get_publication_quality_plot
        plt = get_publication_quality_plot(12, 8)
        e_min = -4
        e_max = 4
        if self._bs.is_metal():
            e_min = -10
            e_max = 10
        count = 1
        for el in self._bs._structure.composition.elements:
            plt.subplot(220 + count)
            self._maketicks(plt)
            for b in range(len(data['distances'])):
                for i in range(self._nb_bands):
                    plt.plot(data['distances'][b],
                             [data['energy'][b][str(Spin.up)][i][j]
                              for j in range(len(data['distances'][b]))], 'b-',
                             linewidth=band_linewidth)
                    if self._bs.is_spin_polarized:
                        plt.plot(data['distances'][b],
                                 [data['energy'][b][str(Spin.down)][i][j]
                                  for j in range(len(data['distances'][b]))],
                                 'r--', linewidth=band_linewidth)
                        for j in range(len(data['energy'][b][str(Spin.up)][i])):
                            plt.plot(data['distances'][b][j],
                                     data['energy'][b][str(Spin.down)][i][j],
                                     'ro',
                                     markersize=sum([proj[b][str(Spin.down)][i][
                                                         j][str(el)][o] for o in
                                                     proj[b]
                                                     [str(Spin.down)][i][j][
                                                         str(el)]]) * 15.0)
                    for j in range(len(data['energy'][b][str(Spin.up)][i])):
                        plt.plot(data['distances'][b][j],
                                 data['energy'][b][str(Spin.up)][i][j], 'bo',
                                 markersize=sum(
                                     [proj[b][str(Spin.up)][i][j][str(el)][o]
                                      for o in proj[b]
                                      [str(Spin.up)][i][j][str(el)]]) * 15.0)
            if ylim is None:
                if self._bs.is_metal():
                    if zero_to_efermi:
                        plt.ylim(e_min, e_max)
                    else:
                        plt.ylim(self._bs.efermi + e_min, self._bs._efermi
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
        if len(self._bs._structure.composition.elements) > 3:
            raise ValueError
        if elt_ordered is None:
            elt_ordered = self._bs._structure.composition.elements
        proj = self._get_projections_by_branches(
            {e.symbol: ['s', 'p', 'd']
             for e in self._bs._structure.composition.elements})
        data = self.bs_plot_data(zero_to_efermi)
        from pymatgen.util.plotting_utils import get_publication_quality_plot
        plt = get_publication_quality_plot(12, 8)

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

        plt.ylim(data['vbm'][0][1] - 4.0, data['cbm'][0][1] + 2.0)
        return plt


def _qvertex_target(data, index):
    """
    Input data should be in the form of a list of a list of floats.
    index is the index of the targeted point
    Returns the vertices of the voronoi construction around this target point.
    """
    from pyhull import qvoronoi
    output = qvoronoi("p QV" + str(index), data)
    output.pop(0)
    output.pop(0)
    return [[float(i) for i in row.split()] for row in output]


def get_lines_voronoi(data):
    from pyhull import qconvex
    output = qconvex("o", data)

    nb_points = int(output[1].split(" ")[0])
    list_lines = []
    list_points = []
    for i in range(2, 2 + nb_points):
        list_points.append([float(c) for c in output[i].strip().split()])
    facets = []
    for i in range(2 + nb_points, len(output)):
        if output[i] != '':
            tmp = output[i].strip().split(" ")
            facets.append([int(tmp[j]) for j in range(1, len(tmp))])

    for i in range(len(facets)):
        for line in itertools.combinations(facets[i], 2):
            for j in range(len(facets)):
                if i != j and line[0] in facets[j] and line[1] in facets[j]:
                    # check if the two facets i and j are not coplanar
                    vector1 = np.array(list_points[facets[j][0]]) \
                              - np.array(list_points[facets[j][1]])
                    vector2 = np.array(list_points[facets[j][0]]) \
                              - np.array(list_points[facets[j][2]])
                    n1 = np.cross(vector1, vector2)
                    vector1 = np.array(list_points[facets[i][0]]) \
                              - np.array(list_points[facets[i][1]])
                    vector2 = np.array(list_points[facets[i][0]]) \
                              - np.array(list_points[facets[i][2]])
                    n2 = np.cross(vector1, vector2)

                    dot = math.fabs(np.dot(n1, n2) / (np.linalg.norm(n1)
                                                      * np.linalg.norm(n2)))
                    if 1.05 > dot > 0.95:
                        continue
                    list_lines.append({'start': list_points[line[0]],
                                       'end': list_points[line[1]]})
                    break
    return list_lines
