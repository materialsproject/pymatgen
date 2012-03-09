#!/usr/bin/env python

'''
This module provides conversion between the Atomic Simulation Environment
Atoms object and pymatgen Structure objects.  It also includes an adaptor
for spglib for spacegroup determination
'''

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 8, 2012"

import re

from pymatgen.core.structure import Structure

try:
    from ase import Atoms
    ase_loaded = True
except ImportError:
    ase_loaded = False

try:
    from pyspglib import spglib
    spglib_loaded = True
except:
    spglib_loaded = False
    
class AseAtomsAdaptor(object):
    '''
    classdocs
    '''
    
    @staticmethod
    def get_atoms(structure):
        """
        Returns ASE Atoms object from pymatgen structure
        """
        if not structure.is_ordered:
            raise ValueError('ASE Atoms only supports ordered structures')
        symbols = [str(site.specie.symbol) for site in structure]
        positions = [site.coords for site in structure]
        cell = structure.lattice.matrix
        return Atoms(symbols=symbols, positions=positions, pbc = True, cell=cell)

    @staticmethod
    def get_structure(atoms):
        """
        Returns pymatgen structure from ASE Atoms
        """
        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()
        lattice = atoms.get_cell()
        return Structure(lattice, symbols, positions, coords_are_cartesian = True)

class SpglibAdaptor(object):
    
    def __init__(self, structure, symprec = 1e-5):
        self._symprec = symprec
        self._atoms = AseAtomsAdaptor.get_atoms(structure)
        self._spacegroup = spglib.get_spacegroup(self._atoms, symprec = self._symprec)
        
    def get_spacegroup(self):
        return self._spacegroup

    def get_spacegroup_symbol(self):
        return re.split("\s+", self._spacegroup)[0]

    def get_spacegroup_number(self):
        sgnum = re.split("\s+", self._spacegroup)[1]
        sgnum = int(re.sub("\D", "", sgnum))
        return sgnum