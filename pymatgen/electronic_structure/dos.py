#!/usr/bin/env python

"""
This module defines classes to represent the density of states, etc.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 20, 2012"

import numpy as np

from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.core.structure import Structure

class Dos(object):

    def __init__(self, efermi, energies, densities):
        self._efermi = efermi
        self._energies = np.array(energies)
        self._dos = {k:np.array(d) for k, d in densities.items()} # should be a dict having as keys the Spin.up and Spin.down objects

    def get_densities(self, spin = None):
        if self._dos == None:
            result = None
        elif spin == None:
            if Spin.down in self._dos:
                result = self._dos[Spin.up] + self._dos[Spin.down]
            else:
                result = self._dos[Spin.up]
        else:
            result = self._dos[spin]
        return result

    @property
    def densities(self):
        return self._dos

    @property
    def efermi(self):
        return self._efermi

    @property
    def energies(self):
        return self._energies

    def __add__(self, other):
        """
        Adds two DOS together. Checks that energy scales are the same. Otherwise, a ValueError is thrown.
        """
        if not (self.energies == other.energies).all():
            raise ValueError("Energies of both DOS are not compatible!")
        densities = {spin: self._dos[spin] + other._dos[spin] for spin in self._dos.keys()}
        return Dos(self.efermi, self.energies, densities)

    def get_interpolated_value(self, energy):
        f = {}
        import scipy.interpolate as spint
        for spin in self._dos.keys():
            f[spin] = spint.interp1d(self._energies, self._dos[spin])(energy)
        return f

    def get_interpolated_gap(self, tol = 0.001, abs_tol = False, spin = None):
        """
        Expects a DOS object and finds the gap
        
        Args:
            tol:
                tolerance in occupations for determining the gap
            abs_tol:
                tolerance an absolute tolerance (True) and a relative one (False)
            spin:
                Possible values are:
                    None - finds the ap in the summed densities
                    Up - finds the gap in the up spin channel
                    Down - finds the gap in teh down spin channel
        
        Returns:
            (gap, cbm, vbm) :  tuple of floats in eV corresponding to the gap, cbm and vbm
        """

        tdos = self.get_densities(spin)
        if abs_tol == False:
            tol = tol * tdos.sum() / tdos.shape[0]
        energies = self._energies
        below_fermi = [i for i in xrange(len(energies)) if energies[i] < self._efermi and tdos[i] > tol]
        above_fermi = [i for i in xrange(len(energies)) if energies[i] > self._efermi and tdos[i] > tol]
        vbm_start = max(below_fermi)
        cbm_start = min(above_fermi)
        import scipy.interpolate as spint
        if vbm_start == cbm_start:
            return 0.0, self._efermi, self._efermi
        else:
            # Interpolate between adjacent values
            f = spint.interp1d(tdos[vbm_start:vbm_start + 2][::-1], energies[vbm_start:vbm_start + 2][::-1])
            start = f(tol)
            f = spint.interp1d(tdos[cbm_start - 1:cbm_start + 1], energies[cbm_start - 1:cbm_start + 1])
            end = f(tol)
            return end - start, end, start

    def get_cbm_vbm(self, tol = 0.001, abs_tol = False, spin = None):
        """
        Expects a DOS object and finds the cbm and vbm.
        
        Args:
            tol: 
                tolerance in occupations for determining the gap
            abs_tol: 
                an absolute tolerance (True) and a relative one (False)
            spin:
                Possible values are:
                    None - finds the gap in the summed densities
                    Up - finds the gap in the up spin channel
                    Down - finds the gap in teh down spin channel
        
        Returns:
            (cbm, vbm): float in eV corresponding to the gap
        """
        #determine tolerance
        tdos = self.get_densities(spin)
        if abs_tol == False:
            tol = tol * tdos.sum() / tdos.shape[0]

        # find index of fermi energy
        i_fermi = 0
        while (self._energies[i_fermi] <= self._efermi):
            i_fermi += 1

        # work backwards until tolerance is reached
        i_gap_start = i_fermi
        while i_gap_start - 1 >= 0 and tdos[i_gap_start - 1] <= tol :
            i_gap_start -= 1

        # work forwards until tolerance is reached
        i_gap_end = i_gap_start
        while i_gap_end < tdos.shape[0] and tdos[i_gap_end] <= tol :
            i_gap_end += 1
        i_gap_end -= 1
        return (self._energies[i_gap_end], self._energies[i_gap_start])

    def get_gap(self, tol = 0.001, abs_tol = False, spin = None):
        """
        Expects a DOS object and finds the gap.
        
        Args:
            tol: 
                tolerance in occupations for determining the gap
            abs_tol: 
                an absolute tolerance (True) and a relative one (False)
            spin:
                Possible values are:
                    None - finds the gap in the summed densities
                    Up - finds the gap in the up spin channel
                    Down - finds the gap in teh down spin channel
        
        Returns:
            gap in eV
        """
        (cbm, vbm) = self.get_cbm_vbm(tol, abs_tol, spin)
        return max(cbm - vbm, 0.0)

    def __str__(self):
        """
        Returns a string which can be easily plotted
        """
        if Spin.down in self._dos:
            stringarray = ["#%30s %30s %30s" % ('Energy', 'DensityUp', 'DensityDown')]
            stringarray.extend(["%.5f %.5f %.5f" % (self._energies[i], self._dos[Spin.up][i], self._dos[Spin.down][i]) for i in range(len(self._energies))])
        else:
            stringarray = ["#%30s %30s" % ('Energy', 'DensityUp')]
            stringarray.extend(["%.5f %.5f" % (self._energies[i], self._dos[Spin.up][i]) for i in range(len(self._energies))])
        return "\n".join(stringarray)

    @staticmethod
    def from_dict(d):
        return Dos(d['efermi'], d['energies'], { Spin.from_int(int(k)):v for k, v in d['densities'].items()})

    @property
    def to_dict(self):
        d = {}
        d['efermi'] = self._efermi
        d['energies'] = list(self._energies)
        d['densities'] = { str(int(spin)) : list(dens) for spin , dens in self._dos.items() }
        return d

class PDos(Dos):
    """
    Projected DOS for a specific orbital. Extends the Dos object.
    """
    def __init__(self, efermi, energies, densities, orbital):
        Dos.__init__(self, efermi, energies, densities)
        self.orbital = orbital

    def __str__(self):
        return "#" + str(self.orbital) + "\n" + super(PDos, self).__str__()


class CompleteDos(Dos):
    """
    This wrapper class defines a total dos, and also provides a list of PDos.
    Mainly used by pymatgen.io.vaspio.Vasprun to create a complete Dos from
    a vasprun.xml file. 
    """

    def __init__(self, structure, total_dos, pdoss):
        """
        Args:
            structure:
                Structure associated with this particular DOS.
            total_dos:
                total Dos for structure
            pdoss:
                a list of array of Pdos.  pdoss corresponds to site order in structure 
        """
        self._efermi = total_dos.efermi
        self._energies = total_dos.energies
        self._dos = total_dos.densities
        self._pdos = pdoss
        self._structure = structure

    @property
    def structure(self):
        return self._structure

    def get_site_orbital_dos(self, site, orbital):
        return self._pdos[site][orbital]

    def get_site_dos(self, site):
        site_dos = None
        for pdos in self._pdos[site].values():
            if site_dos == None:
                site_dos = Dos(pdos.efermi, pdos.energies, pdos.densities)
            else:
                site_dos += pdos
        return site_dos

    def get_spd_dos(self):
        """
        Get orbital projected Dos.
        
        Returns:
            dict of {orbital: Dos}, e.g. {'s': Dos object, ...}
        """
        spd_dos = dict()
        for atom_dos in self._pdos.values():
            for pdos in atom_dos.values():
                orbital_type = pdos.orbital.orbital_type
                if orbital_type not in spd_dos:
                    spd_dos[orbital_type] = Dos(pdos.efermi, pdos.energies, pdos.densities)
                else:
                    spd_dos[orbital_type] += pdos
        return spd_dos

    def get_element_dos(self):
        """
        Get element projected Dos.
        
        Returns:
            dict of {Element: Dos}
        """

        el_dos = dict()
        for site, atom_dos in self._pdos.items():
            el = site.specie
            for pdos in atom_dos.values():
                if el not in el_dos:
                    el_dos[el] = Dos(pdos.efermi, pdos.energies, pdos.densities)
                else:
                    el_dos[el] += pdos
        return el_dos

    @staticmethod
    def from_dict(d):
        tdos = Dos.from_dict(d)
        struct = Structure.from_dict(d['structure'])
        pdoss = {}
        for i in xrange(len(d['pdos'])):
            at = struct[i]
            orb_dos = {}
            for orb_str, odos in d['pdos'][i].items():
                orb = Orbital.from_string(orb_str)
                orb_dos[orb] = PDos(odos['efermi'], odos['energies'], { Spin.from_int(int(k)):v for k, v in odos['densities'].items()}, orb)
            pdoss[at] = orb_dos
        return CompleteDos(struct, tdos, pdoss)

    @property
    def to_dict(self):
        d = {}
        d['efermi'] = self._efermi
        d['structure'] = self._structure.to_dict
        d['energies'] = list(self._energies)
        d['densities'] = { str(int(spin)) : list(dens) for spin , dens in self._dos.items() }
        d['pdos'] = []
        print self._pdos.keys()
        if len(self._pdos) > 0:
            for at in self._structure:
                dd = dict()
                for pdos in self._pdos[at].values():
                    dd[str(pdos.orbital)] = {'efermi' : pdos.efermi, 'energies': list(pdos.energies), 'densities' : { str(int(spin)) : list(dens) for spin , dens in pdos.densities.items() }}
                d['pdos'].append(dd)
            d['atom_dos'] = {str(at) : dos.to_dict for at, dos in self.get_element_dos().items()}
            d['spd_dos'] = {str(orb) : dos.to_dict for orb, dos in self.get_spd_dos().items()}
        return d

    def __str__(self):
        return "Complete DOS for " + str(self._structure)


def plot_dos(dos_dict, zero_at_efermi = True, stack = False, key_sort_func = None, xlim = None, ylim = None):
    """
    Plots a series of Dos using matplotlib.
    
    Args:
        dos_dict:
            dict of {label: Dos}
        zero_at_efermi:
            Whether to shift all Dos to have zero energy at the fermi energy.
            Defaults to True.
        stack:
            Whether to plot the DOS as a stacked area graph
        key_sort_func:
            function used to sort the dos_dict keys.
        xlim:
            Specifies the x-axis limits.
        ylim:
            Specifies the y-axis limits.
    """
    from pymatgen.util.plotting_utils import get_publication_quality_plot
    plt = get_publication_quality_plot(12, 8)
    color_order = ['r', 'b', 'g', 'c']

    y = None
    alldensities = []
    allenergies = []
    keys = list(dos_dict.keys())
    if key_sort_func:
        keys = sorted(keys, key = key_sort_func)
    for key in keys:
        dos = dos_dict[key]
        energies = dos.energies - dos.efermi if zero_at_efermi else dos.energies
        densities = dos.densities
        if not y:
            y = {Spin.up: np.zeros(energies.shape), Spin.down: np.zeros(energies.shape)}
        newdens = {}
        for spin in [Spin.up, Spin.down]:
            if spin in dos.densities:
                if stack:
                    y[spin] += densities[spin]
                    newdens[spin] = y[spin].copy()
                else:
                    newdens[spin] = densities[spin]
        allenergies.append(energies)
        alldensities.append(newdens)

    keys = list(keys)
    keys.reverse()
    alldensities.reverse()

    allpts = []
    for i, key in enumerate(keys):
        x = []
        y = []
        for spin in [Spin.up, Spin.down]:
            if spin in alldensities[i]:
                x.extend(allenergies[i])
                y.extend(int(spin) * alldensities[i][spin])
        allpts.extend(zip(x, y))
        if stack:
            plt.fill(x, y, color = color_order[i % 4], label = str(key))
        else:
            plt.plot(x, y, color = color_order[i % 4], label = str(key))

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

    if zero_at_efermi:
        ylim = plt.ylim()
        plt.plot([0, 0], ylim, 'k--', linewidth = 2)

    plt.legend()
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize = 30)
    plt.tight_layout()
    plt.show()

