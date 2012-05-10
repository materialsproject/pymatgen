#!/usr/bin/env python
"""
This module provides classes and utilities to plot band structures
"""

__author__ = "Geoffroy Hautier, Shyue Ping Ong, Michael Kocher"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "March 14, 2012"

import logging
import math
from pymatgen.electronic_structure.core import Spin

log = logging.getLogger('BSPlotter')


class BSPlotter(object):

    """
    class used to plot or get data to facilitate the plot of band structure line objects
    """

    def __init__(self, bs):
        """
        Arguments:
            bs:
                A Bandstructure_line object.
        """
        self._bs = bs
        #Many ab initio codes do not give good results for the highest occupied bands, we therefore only
        #give 90% of the bands for plotting
        self._nb_bands = int(math.floor(self._bs._nb_bands * 0.9))

    def bs_plot_data(self, zero_to_efermi=True):

        """
        Get the data nicely formatted for a plot
        
        Args:
            zero_to_efermi:
                Automatically subtract off the Fermi energy from the eigenvalues and plot (E-Ef)
        
        Returns:
            A dict of the following format:
                'ticks': a dictionary with the 'distances' at which there is a 
                kpoint (the x axis) and the labels (None if no label)
                'energy': a dictionnary storing bands for spin up and spin down data {Spin:[band_index][k_point_index]} as a list (one element for each band) of energy for 
                each kpoint
                'vbm': a list of tuples (distance,energy) marking the vbms. The energies are shifted with respect to the fermi level is the option has been selected.
                'cbm': a list of tuples (distance,energy) marking the cbms. The energies are shifted with respect to the fermi level is the option has been selected.
                'lattice': the reciprocal lattice
                'zero_energy': this is the energy used as zero for the plot
                'band_gap': a string indicating the band gap and it's nature (empty if it's a metal)
        """
        zero_energy = None
        
        if self._bs.is_metal():
            zero_energy = self._bs.efermi
        else:
            zero_energy = self._bs.get_vbm()['energy']

        if not zero_to_efermi:
            zero_energy = 0.0
        
        energy = {str(Spin.up): []}
        if self._bs.is_spin_polarized:
            energy = {str(Spin.up): [], str(Spin.down): []}
        distance = [self._bs._distance[j] for j in range(len(self._bs._kpoints))]
        ticks = self.get_ticks()
        for i in range(self._nb_bands):
            energy[str(Spin.up)].append([self._bs._bands[Spin.up][i][j] - zero_energy for j in range(len(self._bs._kpoints))])
        if self._bs.is_spin_polarized:
            for i in range(self._nb_bands):
                energy[str(Spin.down)].append([self._bs._bands[Spin.down][i][j] - zero_energy for j in range(len(self._bs._kpoints))])
        
        vbm = self._bs.get_vbm()
        cbm = self._bs.get_cbm()
        
        vbm_plot=[]
        cbm_plot=[]
        
        for index in cbm['kpoint_index']:
            cbm_plot.append((self._bs._distance[index],cbm['energy'] - zero_energy if zero_to_efermi else cbm['energy']))
        
        for index in vbm['kpoint_index']:
            vbm_plot.append((self._bs._distance[index],vbm['energy'] - zero_energy if zero_to_efermi else vbm['energy']))
        
        bg=self._bs.get_band_gap()
        direct="Indirect"
        if bg['direct']:
            direct="Direct"
            
        return {'ticks': ticks, 'distances': distance, 'energy': energy, 'vbm':vbm_plot, 'cbm':cbm_plot, 
                'lattice':self._bs._lattice_rec.to_dict, 'zero_energy':zero_energy, 'band_gap':direct+" "+bg['transition']+" band gap="+str(bg['energy']) if self._bs.is_metal()==False else ""}

    def show(self, file_name=None, zero_to_efermi=True):
        """
        Show the bandstrucure plot. Blue lines are up spin, red lines are down spin
        
        Args:
            file_name:
                File name to write image to (e.g., plot.eps). If None no image is created
            zero_to_efermi:
                Automatically subtract off the Fermi energy from the eigenvalues and plot (E-Ef)
        """
        import pylab
        from matplotlib import rc

        #rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica'], 'size': 20})
        rc('text', usetex=True)

        #main internal config options
        e_min = -4
        e_max = 4
        band_linewidth = 3

        pylab.figure
        data = self.bs_plot_data(zero_to_efermi)
        for i in range(self._nb_bands):
                pylab.plot(data['distances'], [e for e in data['energy'][Spin.up][i]], 'b-', linewidth=band_linewidth)
                if self._bs.is_spin_polarized:
                    pylab.plot(data['distances'], [e for e in data['energy'][Spin.down][i]], 'r-', linewidth=band_linewidth)

        ticks = self.get_ticks()
        # ticks is dict wit keys: distances (array floats), labels (array str)
        log.debug("ticks {t}".format(t=ticks))
        log.debug("ticks has {n} distances and {m} labels".format(n=len(ticks['distance']), m=len(ticks['label'])))
        # Draw lines for BZ boundries
        for i in range(len(ticks['label'])):
            if ticks['label'][i] is not None:
                # don't print the same label twice
                if i != 0:
                    if (ticks['label'][i] == ticks['label'][i - 1]):
                        log.debug("already print label... skipping label {i}".format(i=ticks['label'][i]))
                    else:
                        log.debug("Adding a line at {d} for label {l}".format(d=ticks['distance'][i], l=ticks['label'][i]))
                        pylab.axvline(ticks['distance'][i], color='k')
                else:
                    log.debug("Adding a line at {d} for label {l}".format(d=ticks['distance'][i], l=ticks['label'][i]))
                    pylab.axvline(ticks['distance'][i], color='k')

        #Sanitize only plot the uniq values
        uniq_d = []
        uniq_l = []
        temp_ticks = zip(ticks['distance'], ticks['label'])
        for i in xrange(len(temp_ticks)):
            if i == 0:
                uniq_d.append(temp_ticks[i][0])
                uniq_l.append(temp_ticks[i][1])
                log.debug("Adding label {l} at {d}".format(l=temp_ticks[i][0], d=temp_ticks[i][1]))
            else:
                if temp_ticks[i][1] == temp_ticks[i - 1][1]:
                    log.debug("Skipping label {i}".format(i=temp_ticks[i][1]))
                else:
                    log.debug("Adding label {l} at {d}".format(l=temp_ticks[i][0], d=temp_ticks[i][1]))
                    uniq_d.append(temp_ticks[i][0])
                    uniq_l.append(temp_ticks[i][1])

        log.debug("Unique labels are {i}".format(i=zip(uniq_d, uniq_l)))
        #pylab.gca().set_xticks(ticks['distance'])
        #pylab.gca().set_xticklabels(ticks['label'])
        pylab.gca().set_xticks(uniq_d)
        pylab.gca().set_xticklabels(uniq_l)

        #Main X and Y Labels
        pylab.xlabel(r'$\mathrm{Wave\ Vector}$', fontsize='large')
        ylabel = r'$\mathrm{E\ -\ E_f\ (eV)}$' if zero_to_efermi else r'$\mathrm{Energy\ (eV)}$'
        pylab.ylabel(ylabel, fontsize='large')

        # Draw Fermi energy, only if not the zero
        if not zero_to_efermi:
            ef = self._bs.efermi
            pylab.axhline(ef, linewidth=2, color='k')

        # X range (K)
        #last distance point
        x_max = data['distances'][-1]
        pylab.xlim(0, x_max)

        if self._bs.is_metal():
            # Plot A Metal
            if zero_to_efermi:
                pylab.ylim(e_min, e_max)
            else:
                pylab.ylim(self._bs.efermi + e_min, self._bs._efermi + e_max)
        else:
            # Semiconductor, or Insulator
            # cbm, vbm are dict with keys: kpoint, energy, is_direct
            vbm = self._bs.get_vbm()
            cbm = self._bs.get_cbm()

            e_cbm = cbm['energy'] - self._bs.efermi if zero_to_efermi else cbm['energy']
            e_vbm = vbm['energy'] - self._bs.efermi if zero_to_efermi else vbm['energy']

            for index in cbm['kpoint_index']:
                pylab.scatter(self._bs._distance[index], e_cbm, color='r', marker='o', s=100)

            for index in vbm['kpoint_index']:
                pylab.scatter(self._bs._distance[index], e_vbm, color='g', marker='o', s=100)

            pylab.ylim(e_vbm + e_min, e_cbm + e_max)

        pylab.legend()
        if file_name is not None:
            pylab.plot()
            pylab.savefig(file_name)
            pylab.close()
        else:
            pylab.show()

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
                if c.label != previous_label and previous_branch != this_branch:
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
        plot two band structure for comparison. One is in red the other in blue (no difference in spins)
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
                pylab.plot(data['distances'], data_other[spin][i], 'r--', linewidth=3)


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

        vertex = pymatgen.command_line.qhull_caller.qvertex_target(list_k_points, 13)
        lines = pymatgen.command_line.qhull_caller.get_lines_voronoi(vertex)

        for i in range(len(lines)):
            vertex1 = lines[i]['start']
            vertex2 = lines[i]['end']
            ax.plot([vertex1[0], vertex2[0]], [vertex1[1], vertex2[1]], [vertex1[2], vertex2[2]], color='k')

        for b in self._bs._branches:
            vertex1 = self._bs._kpoints[b['start_index']].cart_coords
            vertex2 = self._bs._kpoints[b['end_index']].cart_coords
            ax.plot([vertex1[0], vertex2[0]], [vertex1[1], vertex2[1]], [vertex1[2], vertex2[2]], color='r', linewidth=3)

        for k in self._bs._kpoints:
            if k.label:
                label = k.label
                if k.label.startswith("\\") or k.label.find("_") != -1:
                    label = "$" + k.label + "$"
                off = 0.01
                ax.text(k.cart_coords[0] + off, k.cart_coords[1] + off, k.cart_coords[2] + off, label, color='b', size='25')
                ax.scatter([k.cart_coords[0]], [k.cart_coords[1]], [k.cart_coords[2]], color='b')

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

