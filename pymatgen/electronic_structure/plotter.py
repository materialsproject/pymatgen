#!/usr/bin/env python

'''
This module implements plotter for DOS and band structure.
'''

from __future__ import division

__author__ = "Shyue Ping Ong, Geoffroy Hautier"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "May 1, 2012"


from collections import OrderedDict

import numpy as np
import logging
import math

from pymatgen.electronic_structure.core import Spin
from pymatgen.util.io_utils import clean_json
from pymatgen.command_line.qhull_caller import qvertex_target, \
    get_lines_voronoi

logger = logging.getLogger('BSPlotter')


class DosPlotter(object):
    """
    Class for plotting DOSes.
    """

    def __init__(self, zero_at_efermi=True, stack=False, sigma=None):
        """
        Args:
            zero_at_efermi:
                Whether to shift all Dos to have zero energy at the fermi
                energy. Defaults to True.
            stack:
                Whether to plot the DOS as a stacked area graph
            key_sort_func:
                function used to sort the dos_dict keys.
            sigma:
                A float specifying a standard deviation for Gaussian smearing
                the DOS for nicer looking plots. Defaults to None for no
                smearing.
        """
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
            dos_dict:
                dict of {label: Dos}
            key_sort_func:
                function used to sort the dos_dict keys.
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
        return clean_json(self._doses)

    def get_plot(self, xlim=None, ylim=None):
        """
        Get a matplotlib plot showing the DOS.

        Args:
            xlim:
                Specifies the x-axis limits. Set to None for automatic
                determination.
            ylim:
                Specifies the y-axis limits.
        """
        from pymatgen.util.plotting_utils import get_publication_quality_plot
        plt = get_publication_quality_plot(12, 8)
        color_order = ['r', 'b', 'g', 'c', 'm', 'k']

        y = None
        alldensities = []
        allenergies = []
        """
        Note that this complicated processing of energies is to allow for
        stacked plots in matplotlib.
        """
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
            allpts.extend(zip(x, y))
            if self.stack:
                plt.fill(x, y, color=color_order[i % len(color_order)],
                         label=str(key))
            else:
                plt.plot(x, y, color=color_order[i % len(color_order)],
                         label=str(key))
            if not self.zero_at_efermi:
                ylim = plt.ylim()
                plt.plot([self._doses[key]['efermi'],
                          self._doses[key]['efermi']], ylim,
                         color_order[i % 4] + '--', linewidth=2)

        plt.xlabel('Energies (eV)')
        plt.ylabel('Density of states')
        if xlim:
            plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)
        else:
            xlim = plt.xlim()
            relevanty = [p[1] for p in allpts
                         if p[0] > xlim[0] and p[0] < xlim[1]]
            plt.ylim((min(relevanty), max(relevanty)))

        if self.zero_at_efermi:
            ylim = plt.ylim()
            plt.plot([0, 0], ylim, 'k--', linewidth=2)

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
            filename:
                Filename to write to.
            img_format:
                Image format to use. Defaults to EPS.
            xlim:
                Specifies the x-axis limits. Set to None for automatic
                determination.
            ylim:
                Specifies the y-axis limits.
        """
        plt = self.get_plot(xlim, ylim)
        plt.savefig(filename, format=img_format)

    def show(self, xlim=None, ylim=None):
        """
        Show the plot using matplotlib.

        Args:
            xlim:
                Specifies the x-axis limits. Set to None for automatic
                determination.
            ylim:
                Specifies the y-axis limits.
        """
        plt = self.get_plot(xlim, ylim)
        plt.show()


class BSPlotter(object):
    """
    Class to plot or get data to facilitate the plot of band structure objects.
    """

    def __init__(self, bs):
        """
        Args:
            bs:
                A BandStructureSymmLine object.
        """
        self._bs = bs
        #Many ab initio codes do not give good results for the highest
        #occupied bands, we therefore only give 90% of the bands for plotting
        self._nb_bands = int(math.floor(self._bs._nb_bands * 0.9))

    def bs_plot_data(self, zero_to_efermi=True):

        """
        Get the data nicely formatted for a plot

        Args:
            zero_to_efermi:
                Automatically subtract off the Fermi energy from the
                eigenvalues and plot.

        Returns:
            A dict of the following format:
                ticks:
                    A dict with the 'distances' at which there is a kpoint (the
                    x axis) and the labels (None if no label)
                energy:
                    A dict storing bands for spin up and spin down data
                    {Spin:[band_index][k_point_index]} as a list (one element
                    for each band) of energy for each kpoint.
                vbm:
                    A list of tuples (distance,energy) marking the vbms. The
                    energies are shifted with respect to the fermi level is the
                    option has been selected.
                cbm:
                    A list of tuples (distance,energy) marking the cbms. The
                    energies are shifted with respect to the fermi level is the
                    option has been selected.
                lattice:
                    The reciprocal lattice.
                zero_energy:
                    This is the energy used as zero for the plot.
                band_gap:
                    A string indicating the band gap and its nature (empty if
                    it's a metal).
                is_metal:
                    True if the band structure is metallic (i.e., there is at
                    least one band crossing the fermi level).
        """
        zero_energy = None

        if self._bs.is_metal():
            zero_energy = self._bs.efermi
        else:
            zero_energy = self._bs.get_vbm()['energy']

        if not zero_to_efermi:
            zero_energy = 0.0

        energy = {str(Spin.up): []}
        kpoints = self._bs._kpoints
        if self._bs.is_spin_polarized:
            energy = {str(Spin.up): [], str(Spin.down): []}
        distance = [self._bs._distance[j]
                    for j in range(len(kpoints))]
        ticks = self.get_ticks()
        for i in range(self._nb_bands):
            energy[str(Spin.up)].append([self._bs._bands[Spin.up][i][j]
                                         - zero_energy
                                         for j in range(len(kpoints))])
        if self._bs.is_spin_polarized:
            for i in range(self._nb_bands):
                energy[str(Spin.down)].append([self._bs._bands[Spin.down][i][j] - zero_energy for j in range(len(self._bs._kpoints))])


        vbm = self._bs.get_vbm()
        cbm = self._bs.get_cbm()

        vbm_plot = []
        cbm_plot = []

        for index in cbm['kpoint_index']:
            cbm_plot.append((self._bs._distance[index],
                             cbm['energy'] - zero_energy if zero_to_efermi \
                             else cbm['energy']))

        for index in vbm['kpoint_index']:
            vbm_plot.append((self._bs._distance[index],
                             vbm['energy'] - zero_energy if zero_to_efermi \
                             else vbm['energy']))

        bg = self._bs.get_band_gap()
        direct = "Indirect"
        if bg['direct']:
            direct = "Direct"

        return {'ticks': ticks, 'distances': distance, 'energy': energy,
                'vbm': vbm_plot, 'cbm': cbm_plot,
                'lattice': self._bs._lattice_rec.to_dict,
                'zero_energy': zero_energy, 'is_metal': self._bs.is_metal(),
                'band_gap': "{} {} bandgap = {}".format(direct,
                                                        bg['transition'],
                                                        bg['energy'])
                if not self._bs.is_metal() else ""}

    def get_plot(self, zero_to_efermi=True, ylim=None):
        """
        get a matplotlib object for the bandstructure plot. 
        Blue lines are up spin, red lines are down
        spin.

        Args:
            zero_to_efermi:
                Automatically subtract off the Fermi energy from the eigenvalues
                and plot (E-Ef).
            ylim
                specify the y-axis (energy) limits; by default None let the code choose.
                It is vbm-4 and cbm+4 if insulator efermi-10 and efermi+10 if metal
        """
        from pymatgen.util.plotting_utils import get_publication_quality_plot
        plt = get_publication_quality_plot(12, 8)
        from matplotlib import rc

        rc('text', usetex=True)

        #main internal config options
        e_min = -4
        e_max = 4
        if self._bs.is_metal():
            e_min = -10
            e_max = 10
        band_linewidth = 3

        #pylab.figure
        data = self.bs_plot_data(zero_to_efermi)
        for i in range(self._nb_bands):
                plt.plot(data['distances'],
                         [e for e in data['energy'][str(Spin.up)][i]],
                         'b-', linewidth=band_linewidth)
                if self._bs.is_spin_polarized:
                    plt.plot(data['distances'],
                             [e for e in data['energy'][str(Spin.down)][i]],
                             'r-', linewidth=band_linewidth)

        ticks = self.get_ticks()
        # ticks is dict wit keys: distances (array floats), labels (array str)
        logger.debug("ticks {t}".format(t=ticks))
        logger.debug("ticks has {n} distances and {m} labels"
                     .format(n=len(ticks['distance']), m=len(ticks['label'])))
        # Draw lines for BZ boundries
        for i in range(len(ticks['label'])):
            if ticks['label'][i] is not None:
                # don't print the same label twice
                if i != 0:
                    if (ticks['label'][i] == ticks['label'][i - 1]):
                        logger.debug("already printed... skipping label {i}"
                                     .format(i=ticks['label'][i]))
                    else:
                        logger.debug("Adding a line at {d} for label {l}"
                                     .format(d=ticks['distance'][i],
                                             l=ticks['label'][i]))
                        plt.axvline(ticks['distance'][i], color='k')
                else:
                    logger.debug("Adding a line at {d} for label {l}"
                                 .format(d=ticks['distance'][i],
                                         l=ticks['label'][i]))
                    plt.axvline(ticks['distance'][i], color='k')

        #Sanitize only plot the uniq values
        uniq_d = []
        uniq_l = []
        temp_ticks = zip(ticks['distance'], ticks['label'])
        for i in xrange(len(temp_ticks)):
            if i == 0:
                uniq_d.append(temp_ticks[i][0])
                uniq_l.append(temp_ticks[i][1])
                logger.debug("Adding label {l} at {d}"
                             .format(l=temp_ticks[i][0], d=temp_ticks[i][1]))
            else:
                if temp_ticks[i][1] == temp_ticks[i - 1][1]:
                    logger.debug("Skipping label {i}"
                                 .format(i=temp_ticks[i][1]))
                else:
                    logger.debug("Adding label {l} at {d}"
                                 .format(l=temp_ticks[i][0],
                                         d=temp_ticks[i][1]))
                    uniq_d.append(temp_ticks[i][0])
                    uniq_l.append(temp_ticks[i][1])

        logger.debug("Unique labels are {i}".format(i=zip(uniq_d, uniq_l)))
        #pylab.gca().set_xticks(ticks['distance'])
        #pylab.gca().set_xticklabels(ticks['label'])
        plt.gca().set_xticks(uniq_d)
        plt.gca().set_xticklabels(uniq_l)

        #Main X and Y Labels
        plt.xlabel(r'$\mathrm{Wave\ Vector}$', fontsize=30)
        ylabel = r'$\mathrm{E\ -\ E_f\ (eV)}$' if zero_to_efermi \
            else r'$\mathrm{Energy\ (eV)}$'
        plt.ylabel(ylabel, fontsize=30)

        # Draw Fermi energy, only if not the zero
        if not zero_to_efermi:
            ef = self._bs.efermi
            plt.axhline(ef, linewidth=2, color='k')

        # X range (K)
        #last distance point
        x_max = data['distances'][-1]
        plt.xlim(0, x_max)

        if ylim == None:
            if self._bs.is_metal():
                # Plot A Metal
                if zero_to_efermi:
                    plt.ylim(e_min, e_max)
                else:
                    plt.ylim(self._bs.efermi + e_min, self._bs._efermi + e_max)
            else:

                for cbm in data['cbm']:
                    plt.scatter(cbm[0], cbm[1], color='r', marker='o', s=100)

                for vbm in data['vbm']:
                    plt.scatter(vbm[0], vbm[1], color='g', marker='o', s=100)

                plt.ylim(data['vbm'][0][1] + e_min, data['cbm'][0][1] + e_max)
        else:
            plt.ylim(ylim)

        plt.tight_layout()

        return plt

    def show(self, zero_to_efermi=True, ylim=None):
        """
        Show the plot using matplotlib.
        
        Args:
            zero_to_efermi:
                Automatically subtract off the Fermi energy from the eigenvalues
                and plot (E-Ef).
            ylim
                specify the y-axis (energy) limits; by default None let the code choose.
                It is vbm-4 and cbm+4 if insulator efermi-10 and efermi+10 if metal
        """
        plt = self.get_plot(zero_to_efermi, ylim)
        plt.show()

    def save_plot(self, filename, img_format="eps", ylim=None, zero_to_efermi=True):
        """
        Save matplotlib plot to a file.
        
        Args:
            filename:
                Filename to write to.
            img_format:
                Image format to use. Defaults to EPS.
            ylim:
                Specifies the y-axis limits. 
        """
        plt = self.get_plot(ylim=ylim, zero_to_efermi=zero_to_efermi)
        plt.savefig(filename, format=img_format)
        plt.close()

    def get_ticks(self):
        """
        Get all ticks and labels for a band structure plot.

        Returns:
            A dict with
                'distance': a list of distance at which ticks should be set.
                'label': a list of label for each of those ticks.
        """
        tick_distance = []
        tick_labels = []
        previous_label = self._bs._kpoints[0].label
        previous_branch = self._bs._branches[0]['name']
        for i, c in enumerate(self._bs._kpoints):
            if c.label != None:
                tick_distance.append(self._bs._distance[i])
                this_branch = None
                for b in self._bs._branches:
                    if i >= b['start_index'] and i <= b['end_index']:
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
                    tick_labels.append(label0 + "$|$" + label1)
                    #print label0+","+label1
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
        (no difference in spins)
        TODO: still a lot of work to do that nicely!
        """
        import pylab
        data = self.bs_plot_data
        data_other = other_plotter.bs_plot_data
        for spin in data:
            for i in range(self._nb_bands):
                pylab.plot(data['distances'], data[spin][i], 'b-', linewidth=3)
        for spin in data_other:
            for i in range(self._nb_bands):
                pylab.plot(data['distances'], data_other[spin][i], 'r--',
                           linewidth=3)

        ticks = self.get_ticks()

        pylab.gca().set_xticks(ticks['distance'])
        pylab.gca().set_xticklabels(ticks['label'])
        pylab.xlabel('Kpoints', fontsize='large')
        pylab.ylabel('Energy(eV)', fontsize='large')
        for i in range(len(ticks['label'])):
            if ticks['label'][i]:
                pylab.axvline(ticks['distance'][i], color='k')
        pylab.show()
        pylab.legend()

    def plot_brillouin(self):
        """
            plot the Brillouin zone
        """
        import pymatgen.command_line.qhull_caller
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        mpl.rcParams['legend.fontsize'] = 10

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        vec1 = self._bs._lattice_rec.matrix[0]
        vec2 = self._bs._lattice_rec.matrix[1]
        vec3 = self._bs._lattice_rec.matrix[2]

        #make the grid
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

        vertex = qvertex_target(list_k_points, 13)
        lines = get_lines_voronoi(vertex)

        for i in range(len(lines)):
            vertex1 = lines[i]['start']
            vertex2 = lines[i]['end']
            ax.plot([vertex1[0], vertex2[0]], [vertex1[1], vertex2[1]],
                    [vertex1[2], vertex2[2]], color='k')

        for b in self._bs._branches:
            vertex1 = self._bs._kpoints[b['start_index']].cart_coords
            vertex2 = self._bs._kpoints[b['end_index']].cart_coords
            ax.plot([vertex1[0], vertex2[0]], [vertex1[1], vertex2[1]],
                    [vertex1[2], vertex2[2]], color='r', linewidth=3)

        for k in self._bs._kpoints:
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
    Class to plot or get data to facilitate the plot of band structure objects projected along orbitals, elements or sites
    """

    def __init__(self, bs):
        """
        Args:
            bs:
                A BandStructureSymmLine object with projections.
        """
        BSPlotter.__init__(self, bs)

    def _maketicks(self, plt):
        """
        utility private method to add ticks to a band structure
        """
        ticks = self.get_ticks()
        #Sanitize only plot the uniq values
        uniq_d = []
        uniq_l = []
        temp_ticks = zip(ticks['distance'], ticks['label'])
        for i in xrange(len(temp_ticks)):
            if i == 0:
                uniq_d.append(temp_ticks[i][0])
                uniq_l.append(temp_ticks[i][1])
                logger.debug("Adding label {l} at {d}".format(l=temp_ticks[i][0], d=temp_ticks[i][1]))
            else:
                if temp_ticks[i][1] == temp_ticks[i - 1][1]:
                    logger.debug("Skipping label {i}".format(i=temp_ticks[i][1]))
                else:
                    logger.debug("Adding label {l} at {d}".format(l=temp_ticks[i][0], d=temp_ticks[i][1]))
                    uniq_d.append(temp_ticks[i][0])
                    uniq_l.append(temp_ticks[i][1])

        logger.debug("Unique labels are {i}".format(i=zip(uniq_d, uniq_l)))
        plt.gca().set_xticks(uniq_d)
        plt.gca().set_xticklabels(uniq_l)

        for i in range(len(ticks['label'])):
            if ticks['label'][i] is not None:
                # don't print the same label twice
                if i != 0:
                    if (ticks['label'][i] == ticks['label'][i - 1]):
                        logger.debug("already print label... skipping label {i}".format(i=ticks['label'][i]))
                    else:
                        logger.debug("Adding a line at {d} for label {l}".format(d=ticks['distance'][i], l=ticks['label'][i]))
                        plt.axvline(ticks['distance'][i], color='k')
                else:
                    logger.debug("Adding a line at {d} for label {l}".format(d=ticks['distance'][i], l=ticks['label'][i]))
                    plt.axvline(ticks['distance'][i], color='k')
        return plt

    def get_projected_plots_dots(self, dictio, zero_to_efermi=True):
        """
        Method returning a plot composed of subplots along different elements
        and orbitals.

        Args:
            dictio:
                the element and orbitals you want a projection on. The format
                is {Element:[Orbitals]} for instance {'Cu':['d','s'],'O':['p']}
                will give projections for Cu on d and s orbitals and on oxygen
                p.

        Returns:
            a pylab object with different subfigures for each projection
            The blue and red colors are for spin up and spin down.
            The bigger the red or blue dot in the band structure the higher
            character for the corresponding element and orbital.
        """
        from pymatgen.util.plotting_utils import get_publication_quality_plot
        fig_number = 0
        for e in dictio:
            for o in dictio[e]:
                fig_number = fig_number + 1
        proj = self._bs.get_projections_on_elts_and_orbitals(dictio)
        data = self.bs_plot_data(zero_to_efermi)
        plt = get_publication_quality_plot(12, 8)
        count = 1
        for el in dictio:
            for o in dictio[el]:
                plt.subplot(100 * math.ceil(fig_number / 2) + 20 + count)
                self._maketicks(plt)
                for i in range(self._nb_bands):
                    plt.plot(data['distances'],
                             [e for e in data['energy'][str(Spin.up)][i]],
                             'b-')
                    if self._bs.is_spin_polarized:
                        plt.plot(data['distances'],
                                 [e for e in data['energy'][str(Spin.down)][i]],
                                 'r-')
                    for j in range(len(data['energy'][str(Spin.up)][i])):
                        plt.plot(data['distances'][j],
                                     data['energy'][str(Spin.up)][i][j], 'bo',
                                     markersize=proj[Spin.up][i][j][str(el)][o] * 15.0)
                    if self._bs.is_spin_polarized:
                        for j in range(len(data['energy'][str(Spin.down)][i])):
                            plt.plot(data['distances'][j],
                                         data['energy'][str(Spin.down)][i][j], 'ro',
                                         markersize=proj[Spin.down][i][j][str(el)][o] * 15.0)
                plt.ylim(data['vbm'][0][1] - 4.0, data['cbm'][0][1] + 4.0)
                plt.title(str(el) + " " + str(o))
                count = count + 1
        return plt

    def get_elt_projected_plots(self, zero_to_efermi=True):
        """
        Method returning a plot composed of subplots along different elements
        
        Returns:
            a pylab object with different subfigures for each projection
            The blue and red colors are for spin up and spin down
            The bigger the red or blue dot in the band structure the higher 
            character for the corresponding element and orbital
        """
        proj = self._bs.get_projection_on_elements()
        data = self.bs_plot_data(zero_to_efermi)
        from pymatgen.util.plotting_utils import get_publication_quality_plot
        plt = get_publication_quality_plot(12, 8)
        ticks = self.get_ticks()
        #Sanitize only plot the uniq values
        uniq_d = []
        uniq_l = []
        temp_ticks = zip(ticks['distance'], ticks['label'])
        for i in xrange(len(temp_ticks)):
            if i == 0:
                uniq_d.append(temp_ticks[i][0])
                uniq_l.append(temp_ticks[i][1])
                logger.debug("Adding label {l} at {d}".format(l=temp_ticks[i][0], d=temp_ticks[i][1]))
            else:
                if temp_ticks[i][1] == temp_ticks[i - 1][1]:
                    logger.debug("Skipping label {i}".format(i=temp_ticks[i][1]))
                else:
                    logger.debug("Adding label {l} at {d}".format(l=temp_ticks[i][0], d=temp_ticks[i][1]))
                    uniq_d.append(temp_ticks[i][0])
                    uniq_l.append(temp_ticks[i][1])

        logger.debug("Unique labels are {i}".format(i=zip(uniq_d, uniq_l)))
        plt.gca().set_xticks(uniq_d)
        plt.gca().set_xticklabels(uniq_l)

        #Main X and Y Labels
        plt.xlabel(r'Wave vector', fontsize=30)
        #ylabel = r'$\mathrm{E\ -\ E_f\ (eV)}$' if zero_to_efermi else r'$\mathrm{Energy\ (eV)}$'
        ylabel = 'Energy (eV)'
        plt.ylabel(ylabel, fontsize=30)
        for i in range(len(ticks['label'])):
            if ticks['label'][i] is not None:
                # don't print the same label twice
                if i != 0:
                    if (ticks['label'][i] == ticks['label'][i - 1]):
                        logger.debug("already print label... skipping label {i}".format(i=ticks['label'][i]))
                    else:
                        logger.debug("Adding a line at {d} for label {l}".format(d=ticks['distance'][i], l=ticks['label'][i]))
                        plt.axvline(ticks['distance'][i], color='k')
                else:
                    logger.debug("Adding a line at {d} for label {l}".format(d=ticks['distance'][i], l=ticks['label'][i]))
                    plt.axvline(ticks['distance'][i], color='k')
        count = 1
        for el in self._bs._structure.composition.elements:
            plt.subplot(220 + count)
            print count
            self._maketicks(plt)
            for i in range(self._nb_bands):
                plt.plot(data['distances'],
                             [e for e in data['energy'][str(Spin.up)][i]],
                             'b-')
                if self._bs.is_spin_polarized:
                    plt.plot(data['distances'],
                             [e for e in data['energy'][str(Spin.down)][i]],
                             'r-')
                for j in range(len(data['energy'][str(Spin.up)][i])):
                    plt.plot(data['distances'][j],
                                 data['energy'][str(Spin.up)][i][j], 'bo',
                                 markersize=proj[Spin.up][i][j][str(el)] * 15.0)
                if self._bs.is_spin_polarized:
                    for j in range(len(data['energy'][str(Spin.down)][i])):
                        plt.plot(data['distances'][j],
                                     data['energy'][str(Spin.down)][i][j], 'ro',
                                     markersize=proj[Spin.down][i][j][str(el)] * 15.0)
            plt.ylim(data['vbm'][0][1] - 4.0, data['cbm'][0][1] + 4.0)
            plt.title(str(el))
            count = count + 1
        return plt

    def get_elt_projected_plots_color(self, zero_to_efermi=True, elt_ordered=None):
        """
        returns a pylab plot object with one plot where the band structure line color depends
        on the character of the band (along different elements). Each element is associated with red, green or blue
        and the corresponding rgb color depending on the character of the band is used. The method can only deal with binary and ternary compounds
        
        the method does not make a difference for now between spin up and spin down
        
        Args:
            elt_ordered:
                a list of Element ordered. The first one is red, second green, last blue
        
        Returns:
            a pylab object
        
        """
        if len(self._bs._structure.composition.elements) > 3:
            raise ValueError
        if elt_ordered == None:
            elt_ordered = self._bs._structure.composition.elements
        proj = self._bs.get_projection_on_elements()
        data = self.bs_plot_data(zero_to_efermi)
        from pymatgen.util.plotting_utils import get_publication_quality_plot
        plt = get_publication_quality_plot(12, 8)
        count = 1
        ticks = self.get_ticks()
        #Sanitize only plot the uniq values
        uniq_d = []
        uniq_l = []
        temp_ticks = zip(ticks['distance'], ticks['label'])
        for i in xrange(len(temp_ticks)):
            if i == 0:
                uniq_d.append(temp_ticks[i][0])
                uniq_l.append(temp_ticks[i][1])
                logger.debug("Adding label {l} at {d}".format(l=temp_ticks[i][0], d=temp_ticks[i][1]))
            else:
                if temp_ticks[i][1] == temp_ticks[i - 1][1]:
                    logger.debug("Skipping label {i}".format(i=temp_ticks[i][1]))
                else:
                    logger.debug("Adding label {l} at {d}".format(l=temp_ticks[i][0], d=temp_ticks[i][1]))
                    uniq_d.append(temp_ticks[i][0])
                    uniq_l.append(temp_ticks[i][1])

        logger.debug("Unique labels are {i}".format(i=zip(uniq_d, uniq_l)))
        plt.gca().set_xticks(uniq_d)
        plt.gca().set_xticklabels(uniq_l)

        #Main X and Y Labels
        plt.xlabel(r'Wave vector', fontsize=30)
        #ylabel = r'$\mathrm{E\ -\ E_f\ (eV)}$' if zero_to_efermi else r'$\mathrm{Energy\ (eV)}$'
        ylabel = 'Energy (eV)'
        plt.ylabel(ylabel, fontsize=30)
        for i in range(len(ticks['label'])):
            if ticks['label'][i] is not None:
                # don't print the same label twice
                if i != 0:
                    if (ticks['label'][i] == ticks['label'][i - 1]):
                        logger.debug("already print label... skipping label {i}".format(i=ticks['label'][i]))
                    else:
                        logger.debug("Adding a line at {d} for label {l}".format(d=ticks['distance'][i], l=ticks['label'][i]))
                        plt.axvline(ticks['distance'][i], color='k')
                else:
                    logger.debug("Adding a line at {d} for label {l}".format(d=ticks['distance'][i], l=ticks['label'][i]))
                    plt.axvline(ticks['distance'][i], color='k')
        spins = [Spin.up]
        if self._bs.is_spin_polarized:
            spins = [Spin.up, Spin.down]
        for s in spins:
            for i in range(self._nb_bands):
                for j in range(len(data['energy'][str(s)][i]) - 1):
                    #if j%2!=0:
                    #    continue
                    #print j
                    sum = 0.0
                    for el in elt_ordered:
                        sum = sum + proj[s][i][j][str(el)]
                    if sum == 0.0:
                        color = [0.0 for e in elt_ordered]
                    else:
                        color = [proj[s][i][j][str(el)] / sum for el in elt_ordered]
                    if len(color) == 2:
                        color.append(0.0)
                        color[2] = color[1]
                        color[1] = 0.0
                    plt.plot([data['distances'][j], data['distances'][j + 1]],
                                 [data['energy'][str(s)][i][j], data['energy'][str(s)][i][j + 1]],
                                 color=color
                                 , linewidth=3)
        plt.ylim(data['vbm'][0][1] - 4.0, data['cbm'][0][1] + 2.0)
        count = count + 1
        return plt


