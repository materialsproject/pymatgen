#!/usr/bin/env python

"""
This module provides classes to define electronic structure, such as the density of states, etc.
"""

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ ="$Sep 23, 2011M$"

import numpy as np
import scipy.interpolate as spint

class _SpinImpl(object):
    """
    Internal representation of a Spin. 
    Do not use directly.
    """
       
    def __init__(self,name):
        self._name = name
            
    def __int__(self):
        return 1 if self._name == "up" else -1
        
    def __repr__(self):
        return self._name
        
    def __eq__(self, other):
        if other == None:
            return False
        return self._name == other._name
        
    def __hash__(self):
        return self.__int__()
           
    def __str__(self):
        return self._name

class Spin(object):
    """
    Enum type for Spin.  Only up and down.  Design follows somewhat the familiar Java syntax.
    """
    
    up = _SpinImpl("up")
    down = _SpinImpl("down")
    all_spins = (up, down)
    
    @staticmethod
    def from_int(i):
        if i == 1:
            return Spin.up
        elif i == -1:
            return Spin.down
        else:
            raise ValueError("Spin integers must be 1 or -1")
        
class _OrbitalImpl(object):
    """
    Internal representation of an orbital.  Do not use directly. 
    Use the Orbital class enum types.
    """
    
    def __init__(self,name, vasp_index):
        self._name = name
        self._vasp_index = vasp_index
            
    def __int__(self):
        return self._vasp_index
        
    def __repr__(self):
        return self._name
        
    def __eq__(self, other):
        if other == None:
            return False
        return self._name == other._name
        
    def __hash__(self):
        return self.__int__()
        
    @property
    def orbital_type(self):
        return self._name[0].upper()
        
    def __str__(self):
        return self._name

class Orbital(object):
    """
    Enum type for OrbitalType. Indices are basically the azimutal quantum number, l.
    Design follows somewhat the familiar Java syntax.
    """
    
    s = _OrbitalImpl("s",0)
    py = _OrbitalImpl("py",1)
    pz = _OrbitalImpl("pz",2)
    px = _OrbitalImpl("px",3)
    dxy = _OrbitalImpl("dxy",4)
    dyz = _OrbitalImpl("dyz",5)
    dz2 = _OrbitalImpl("dz2",6)
    dxz = _OrbitalImpl("dxz",7)
    dx2 = _OrbitalImpl("dx2",8)
    all_orbitals = (s, py, pz, px, dxy, dyz, dz2, dxz, dx2)
            
    @staticmethod
    def from_vasp_index(i):
        for orb in Orbital.all_orbitals:
            if int(orb) == i:
                return orb
    
    @staticmethod
    def from_string(orb_str):
        for orb in Orbital.all_orbitals:
            if str(orb) == orb_str:
                return orb
        raise ValueError("Illegal orbital definition!")

class Dos(object):

    def __init__(self, efermi, energies, densities):
        self._efermi = efermi
        self._energies = np.array(energies)
        self._dos = {k:np.array(d) for k,d in densities.items()} # should be a dict having as keys the Spin.up and Spin.down objects

    def get_densities(self, spin=None):
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

    def __add__(self,other):
        """
        Adds two DOS together. Checks that energy scales are the same. Otherwise, a ValueError is thrown.
        """
        if not (self.energies == other.energies).all():
            raise ValueError("Energies of both DOS are not compatible!")
        densities = {spin: self._dos[spin] + other._dos[spin] for spin in self.densities.keys()}
        return Dos(self.efermi, self.energies, densities)

    def get_interpolated_value(self,energy):
        f = {}
        for spin in self._dos.keys():
            f[spin] = spint.interp1d(self._energies, self._dos[spin])(energy)
        return f
    
    def get_interpolated_gap(self,tol=0.001,abs_tol=False,spin=None):
        """
        Expects a DOS object and finds the gap
        
        Arguments:
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
        
        if vbm_start == cbm_start:
            return 0.0,  self._efermi,self._efermi
        else:
            # Interpolate between adjacent values
            f = spint.interp1d(tdos[vbm_start:vbm_start+2][::-1], energies[vbm_start:vbm_start+2][::-1])
            start = f(tol)
            f = spint.interp1d(tdos[cbm_start-1:cbm_start+1], energies[cbm_start-1:cbm_start+1])
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
        i_fermi=0
        while (self._energies[i_fermi] <= self._efermi):
            i_fermi+=1
    
        # work backwards until tolerance is reached
        i_gap_start=i_fermi
        while i_gap_start-1 >=0 and tdos[i_gap_start-1] <= tol :
            i_gap_start -= 1
    
        # work forwards until tolerance is reached
        i_gap_end=i_gap_start
        while i_gap_end < tdos.shape[0] and tdos[i_gap_end] <= tol :
            i_gap_end += 1
        i_gap_end -= 1
        return (self._energies[i_gap_end],self._energies[i_gap_start])

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
        (cbm,vbm) = self.get_cbm_vbm(tol,abs_tol,spin)
        return max(cbm - vbm,0.0)

    def __str__(self):
        """
        Returns a string which can be easily plotted
        """
        if Spin.down in self._dos:
            stringarray = ["#%30s %30s %30s" %('Energy','DensityUp','DensityDown')]
            stringarray.extend(["%.5f %.5f %.5f" % (self._energies[i], self._dos[Spin.up][i], self._dos[Spin.down][i]) for i in range(len(self._energies))])
        else:
            stringarray = ["#%30s %30s" %('Energy','DensityUp')]
            stringarray.extend(["%.5f %.5f" % (self._energies[i], self._dos[Spin.up][i]) for i in range(len(self._energies))])
        return "\n".join(stringarray)

class PDos(Dos):
    """
    Projected DOS for a specific orbital. Extends the Dos object.
    """
    def __init__(self,efermi,energies,densities,orbital):
        Dos.__init__(self,efermi,energies,densities)
        self.orbital = orbital
        
    def __str__(self):
        return "#"+str(self.orbital) + "\n" + super(PDos,self).__str__()

class CompleteDos(Dos):
    """
    This wrapper class defines a total dos, and also provides a list of Pdos.
    Mainly used by pymatgen.io.vaspio.Vasprun to create a complete Dos from
    a vasprun.xml file. 
    """
    
    def __init__(self, structure, total_dos, pdoss):
        """
        Arguments:
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
        self._pdos = {structure[i]:{Orbital.from_vasp_index(j) : pdoss[i][j] for j in range(len(pdoss[i]))} for i in range(structure.num_sites)}
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
                site_dos = Dos(pdos.efermi, pdos.energies,pdos.densities)
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
                    spd_dos[orbital_type] = Dos(pdos.efermi, pdos.energies,pdos.densities)
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
                    el_dos[el] = Dos(pdos.efermi, pdos.energies,pdos.densities)
                else:
                    el_dos[el] += pdos
        return el_dos
    
    def __str__(self):
        return "Complete DOS for "+str(self._structure)
    
def plot_dos(dos_dict, zero_at_efermi = True):
    """
    Plots a series of Dos using matplotlib.
    
    Args:
        dos_dict:
            dict of {label: Dos}
        zero_at_efermi:
            Whether to shift all Dos to have zero energy at the fermi energy.
            Defaults to True.
    """
    import pylab
    color_order = ['r', 'b', 'g', 'y']
    count = 0
    for key,dos in dos_dict.items():
        energies = dos.energies - dos.efermi if zero_at_efermi else dos.energies
        densities = dos.densities
        if Spin.up in densities:
            pylab.plot(energies, densities[Spin.up], color_order[count % 4], label=str(key) + ' up')
        if Spin.down in densities:
            pylab.plot(energies, - densities[Spin.down], color_order[count % 4], label=str(key) + ' down')
        count += 1
    
    pylab.xlabel('Energies (eV)', fontsize = 'large')
    pylab.ylabel('Density of states', fontsize = 'large')
    pylab.legend()
    pylab.show()

    
    