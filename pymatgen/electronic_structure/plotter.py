#!/usr/bin/env python

'''
This module implements plotter for DOS.
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "May 1, 2012"


from collections import OrderedDict

import numpy as np

from pymatgen.electronic_structure.core import Spin
from pymatgen.util.io_utils import clean_json


class DosPlotter(object):
    """
    Class for plotting DOSes.
    """

    def __init__(self, zero_at_efermi=True, stack=False,
                 sigma=None):
        """
        Args:
            zero_at_efermi:
                Whether to shift all Dos to have zero energy at the fermi energy.
                Defaults to True.
            stack:
                Whether to plot the DOS as a stacked area graph
            key_sort_func:
                function used to sort the dos_dict keys.
            sigma:
                A float specifying a standard deviation for Gaussian smearing the 
                DOS for nicer looking plots. Defaults to None for no smearing.
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
        energies = dos.energies - dos.efermi if self.zero_at_efermi else dos.energies
        densities = dos.get_smeared_densities(self.sigma) if self.sigma else dos.densities
        efermi = dos.efermi
        self._doses[label] = {'energies':energies, 'densities': densities, 'efermi':efermi}

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
        color_order = ['r', 'b', 'g', 'c']

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
                y = {Spin.up: np.zeros(energies.shape), Spin.down: np.zeros(energies.shape)}
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
                plt.fill(x, y, color=color_order[i % 4], label=str(key))
            else:
                plt.plot(x, y, color=color_order[i % 4], label=str(key))
            if not self.zero_at_efermi:
                ylim = plt.ylim()
                plt.plot([self._doses[key]['efermi'], self._doses[key]['efermi']], ylim, color_order[i % 4] + '--', linewidth=2)

        plt.xlabel('Energies (eV)')
        plt.ylabel('Density of states')
        if xlim:
            plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)
        else:
            xlim = plt.xlim()
            relevanty = [p[1] for p in allpts if p[0] > xlim[0] and p[0] < xlim[1]]
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

