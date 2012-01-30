#!/usr/bin/env python

"""
Wrapper classes for Cif input and output from pymatgen.core.structure.Structures.
"""

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ ="Sep 23, 2011"

import re
import StringIO
import math
import warnings
from collections import OrderedDict

import CifFile
import numpy as np

from pymatgen.core.periodic_table import Element, Specie
from pymatgen.util.io_utils import file_open_zip_aware
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure, Composition

class CifParser:
    '''
    A wrapper class around PyCifRW to read Cif and convert into a pymatgen Structure object.
    '''
    
    def __init__(self, filename):
        """
        Args:
            filename - cif file name.  bzipped or gzipped cifs are fine too.
        """
        if isinstance(filename, basestring):
            with file_open_zip_aware(filename, "r") as f:
                self._cif = CifFile.ReadCif(f)
        else:
            self._cif = CifFile.ReadCif(filename)
                    

    @staticmethod
    def from_string(cif_string):
        output = StringIO.StringIO()
        output.write(cif_string)
        output.seek(0)
        return CifParser(output)
        
    def _unique_coords(self,coord_in,sympos, primitive, lattice, primlattice):
        """
        Generate unique coordinates using coord and symmetry positions.
        """
        coords = list()
        (x, y, z) = coord_in
        for gen in sympos:
            coord = eval(gen)
            if primitive:
                cart = lattice.get_cartesian_coords(np.array(coord))
                coord = primlattice.get_fractional_coords(cart)
            coord = np.array([i - math.floor(i) for i in coord])
            if not coord_in_list(coord,coords,1e-3):
                coords.append(coord)
        return coords
    
    def _get_structure(self, data, primitive):
        """
        Generate structure from part of the cif.
        """
        spacegroup = data['_symmetry_space_group_name_H-M']
        
        if len(spacegroup) == 0:
            latt_type = "P"
        else:
            latt_type = spacegroup[0]
        lengths = [float_from_string(data['_cell_length_'+i]) for i in ['a','b','c']]
        angles = [float_from_string(data['_cell_angle_'+i]) for i in ['alpha','beta','gamma']]
        lattice = Lattice.from_lengths_and_angles(lengths,angles)
        primlattice = lattice.get_primitive_lattice(latt_type)
        try:
            sympos =  data['_symmetry_equiv_pos_as_xyz']
        except:
            try:
                sympos =  data['_symmetry_equiv_pos_as_xyz_']
            except:
                warnings.warn("No _symmetry_equiv_pos_as_xyz type key found. Defaulting to P1.")
                sympos
        def parse_symbol(sym):
            m = re.search("([A-Z][a-z]*)",sym)
            if m:
                return m.group(1)
            return ''
        
        #oxi_states = None
        try:
            oxi_states = dict()
            for i in xrange(len(data['_atom_type_symbol'])):
                oxi_states[data['_atom_type_symbol'][i]] = float_from_string(data['_atom_type_oxidation_number'][i])
        except:
            oxi_states = None
            
        coord_to_species = OrderedDict()
        
        for i in xrange(len(data['_atom_site_type_symbol'])):
            symbol = parse_symbol(data['_atom_site_type_symbol'][i])
            if oxi_states != None:
                el = Specie(symbol, oxi_states[data['_atom_site_type_symbol'][i]])
            else:
                el = Element(symbol)
            x = float_from_string(data['_atom_site_fract_x'][i])
            y = float_from_string(data['_atom_site_fract_y'][i])
            z = float_from_string(data['_atom_site_fract_z'][i])
            try:
                occu = float_from_string(data['_atom_site_occupancy'][i])
            except:
                occu = 1
            coord = (x,y,z)
            if coord not in coord_to_species:
                coord_to_species[coord] = {el:occu}
            else:
                coord_to_species[coord][el] = occu
        
        allspecies = list()
        allcoords = list()
            
        for coord, species in coord_to_species.items():
            coords = self._unique_coords(coord,sympos, primitive, lattice, primlattice)
            allcoords.extend(coords)
            allspecies.extend(len(coords) * [species])
                   
        if primitive:
            return Structure(primlattice,allspecies,allcoords)
        else:
            return Structure(lattice,allspecies,allcoords)
    
    def get_structures(self, primitive=True):
        '''
        Return list of structures in CIF file. primitive boolean sets whether a
        conventional cell structure or primitive cell structure is returned.
        
        Arguments:
            primitive:
                Set to False to return conventional unit cells.  Defaults to True.
        
        Returns:
            List of Structures.
        '''
        return [self._get_structure(v, primitive) for k, v in self._cif.items()]
    
class CifWriter:
    '''
    A wrapper around PyCifRW to write CIF files from pymatgen structures.
    '''
    
    def __init__(self, struct):
        """
        Arguments:
            struct:
                A pymatgen.core.structure.Structure object.
        """
        block = CifFile.CifBlock()
        latt = struct.lattice
        comp = struct.composition
        no_oxi_comp = Composition.from_formula(comp.formula)
        block['_symmetry_space_group_name_H-M'] = 'P 1'
        block['_cell_length_a'] = str(latt.a)
        block['_cell_length_b'] = str(latt.b)
        block['_cell_length_c'] = str(latt.c)
        block['_cell_angle_alpha'] = str(latt.alpha)
        block['_cell_angle_beta'] = str(latt.beta)
        block['_cell_angle_gamma'] = str(latt.gamma)
        block['_chemical_name_systematic'] = "Generated by pymatgen"
        block['_symmetry_Int_Tables_number'] = 1
        block['_chemical_formula_structural'] = str(no_oxi_comp.reduced_formula)
        block['_chemical_formula_sum'] = str(no_oxi_comp.formula)
        block['_cell_volume'] = str(latt.volume)
        
        reduced_comp = Composition.from_formula(no_oxi_comp.reduced_formula)
        el = no_oxi_comp.elements[0]
        amt = comp[el]
        fu = int(amt/reduced_comp[Element(el.symbol)])
        
        block['_cell_formula_units_Z'] = str(fu)
        block.AddCifItem(([['_symmetry_equiv_pos_site_id','_symmetry_equiv_pos_as_xyz']], [[['1'],['x, y, z']]]))
        
        contains_oxidation = True
        symbol_to_oxinum = dict()
        try:
            symbol_to_oxinum = {str(el):el.oxi_state for el in comp.elements}
        except:
            symbol_to_oxinum = {el.symbol:0 for el in comp.elements}
            contains_oxidation = False
        if contains_oxidation:
            block.AddCifItem(([['_atom_type_symbol','_atom_type_oxidation_number']], [[symbol_to_oxinum.keys(),symbol_to_oxinum.values()]]))
        
        atom_site_type_symbol = []
        atom_site_symmetry_multiplicity = []
        atom_site_fract_x = []
        atom_site_fract_y = []
        atom_site_fract_z = []
        atom_site_attached_hydrogens = []
        atom_site_B_iso_or_equiv = []
        atom_site_label = []
        atom_site_occupancy = []
        count = 1
        for site in struct:
            for sp, occu in site.species_and_occu.items():
                atom_site_type_symbol.append(str(sp))
                atom_site_symmetry_multiplicity.append('1')
                atom_site_fract_x.append(str(site.a))
                atom_site_fract_y.append(str(site.b))
                atom_site_fract_z.append(str(site.c))
                atom_site_attached_hydrogens.append('0')
                atom_site_B_iso_or_equiv.append('.')
                atom_site_label.append(str(sp.symbol)+str(count))
                atom_site_occupancy.append(str(occu))
                count += 1
                
        block['_atom_site_type_symbol'] = atom_site_type_symbol
        block.AddToLoop('_atom_site_type_symbol', {'_atom_site_symmetry_multiplicity':atom_site_symmetry_multiplicity})
        block.AddToLoop('_atom_site_type_symbol', {'_atom_site_fract_x':atom_site_fract_x})
        block.AddToLoop('_atom_site_type_symbol', {'_atom_site_fract_y':atom_site_fract_y})
        block.AddToLoop('_atom_site_type_symbol', {'_atom_site_fract_z':atom_site_fract_z})
        block.AddToLoop('_atom_site_type_symbol', {'_atom_site_attached_hydrogens':atom_site_attached_hydrogens})
        block.AddToLoop('_atom_site_type_symbol', {'_atom_site_B_iso_or_equiv':atom_site_B_iso_or_equiv})
        block.AddToLoop('_atom_site_type_symbol', {'_atom_site_label':atom_site_label})
        block.AddToLoop('_atom_site_type_symbol', {'_atom_site_occupancy':atom_site_occupancy})
        
        self._cf = CifFile.CifFile()
        self._cf[comp.reduced_formula] = block
    
    def __str__(self):
        '''
        Returns the cif as a string.
        '''
        return str(self._cf)
    
    def write_file(self, filename):
        '''
        Write the cif file.
        '''
        with open(filename,'w') as f:
            f.write(self.__str__())

def around_diff_num(a,b):
    """
    Used to compare differences in fractional coordinates, taking into account PBC. 
    """
    diff_num = abs(a-b)
    return diff_num if diff_num < 0.5 else abs(1 - diff_num)

def coord_in_list(coord, coord_list, tol):
    """
    Helper method to check if coord is already in a list of coords, subject to a tolerance.
    """ 
    for c in coord_list:
        diff = np.array([around_diff_num(c[i], coord[i]) for i in xrange(3)])
        if (diff < tol).all():
            return True
    return False

def float_from_string(text):
    '''
    Remove uncertainty brackets from strings and return the float.
    '''
    return float(re.sub('\(\d+\)','',text))
