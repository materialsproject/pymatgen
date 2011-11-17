#!/usr/bin/env python

'''
This module defines the VaspParameterSet abstract base class and
a concrete implementation for the Materials Project.  The basic
concept behind a parameter set is to specify a scheme to generate
a consistent set of Vasp inputs from a structure without further
user intervention.
'''

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Nov 16, 2011"

import os
import abc
import ConfigParser
import json

from pymatgen.io.vaspio import Incar, Poscar, Potcar, Kpoints

class VaspParameterSet(object):
    """
    Abstract base class representing a set of Vasp input parameters.
    The idea is that using a VaspParameterSet, a complete set of input files (INPUT, KPOINTS, POSCAR and POTCAR)
    can be generated in an automated fashion for any structure.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def get_poscar(self, structure):
        '''
        Returns Poscar from a structure.
        '''
        return Poscar(structure)
    
    @abc.abstractmethod 
    def get_kpoints(self, structure):
        '''
        Returns Kpoints from a structure.
        
        Arguments:
            structure:
                Structure object
        '''
        return
    
    @abc.abstractmethod 
    def get_incar(self, structure):
        '''
        Returns Incar from a structure.
        
        Arguments:
            structure:
                Structure object
        '''
        return
    
    @abc.abstractmethod 
    def get_potcar(self, structure):
        '''
        Returns Potcar from a structure.
        
        Arguments:
            structure:
                Structure object
        '''
        return
    
    def get_all_vasp_input(self, structure):
        '''
        Returns all input files
        '''
        return {'INCAR':self.get_incar(structure), 'KPOINTS':self.get_kpoints(structure), 
                'POSCAR': self.get_poscar(structure), 'POTCAR': self.get_potcar(structure)}
    
    def write_input(self, structure, output_dir, make_dir_if_not_present = True):
        if make_dir_if_not_present and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        for k,v in self.get_all_vasp_input(structure).items():
            v.write_file(os.path.join(output_dir, k))
    

class MaterialsProjectVaspParameterSet(VaspParameterSet):
    """
    Class representing a set of Vasp input parameters.
    The idea is that using a VaspParameterSet, a complete set of input files (INPUT, KPOINTS, POSCAR and POTCAR)
    can be generated in an automated fashion for any structure.
    """
    def __init__(self):
        """
        Args:
            name: 
                Name of parameters set to use.  Defaults to MaterialsProject parameter set.
        """
        self.name = "MaterialsProject"
        module_dir = os.path.dirname(os.path.abspath(__file__))
        self._config = ConfigParser.SafeConfigParser()
        self._config.optionxform = str
        self._config.readfp(open(os.path.join(module_dir, "VaspParameterSets.cfg")))
        self.potcar_settings = dict(self._config.items(self.name + 'POTCAR'))
        self.kpoints_settings = dict(self._config.items(self.name + 'KPOINTS'))
        self.incar_settings = dict(self._config.items(self.name+'INCAR'))
        for key in ['MAGMOM', 'LDAUU', 'LDAUJ', 'LDAUL']:
            self.incar_settings[key] = json.loads(self.incar_settings[key])

    def get_incar(self, structure):
        incar = Incar()
        poscar = Poscar(structure)
        for key, setting in self.incar_settings.items():
            if key == "MAGMOM":
                incar[key] =  [setting.get(site.specie.symbol, 0.6) for site in structure]
            elif key in ['LDAUU', 'LDAUJ', 'LDAUL']:
                incar[key] =  [setting.get(sym, 0) for sym in poscar.site_symbols]
            else:
                incar[key] = setting
                
        has_u = sum(incar['LDAUU']) > 0
        if not has_u:
            for key in incar.keys():
                if key.startswith('LDAU'):
                    del incar[key]
        
        return incar
    
    def get_poscar(self, structure):
        return super(MaterialsProjectVaspParameterSet, self).get_poscar(structure)

    def get_potcar(self, structure):               
        p = self.get_poscar(structure)
        elements = p.site_symbols
        potcar_symbols = []
        for el in elements:
            potcar_symbols.append(self.potcar_settings[el] if el in self.potcar_settings else el)
        return Potcar(potcar_symbols)
        

    def get_kpoints(self, structure):
        '''
        Writes out a KPOINTS file using the fully automated grid method
        The parameter @param tetra controls whether the KPOINTS file instructs
        VASP to use the tetrahedron method or not, @param kppra is the desired
        number of kpoints * the number of atoms. 
        
        NOTE: uses Gamma centered meshes for hexagonal cells and Monkhorst-Pack grids otherwise.
        
        Algorithm: 
            Uses a simple approach scaling the number of divisions along each 
            reciprocal lattice vector proportional to its length. So of N_{grid}
            grid points are desired, then 
            N_{grid} = n_{a}*n_{b}*n_{c}
            n_{b} = l_{b}/l_{a} * n_{a}
            n_{c} = l_{c}/l_{a} * n_{a}
            
            and
            
            n_{a} = round( ( (l_{a}**2 / (l_{b}*l_{c}) )*N_{grid} )**(1/3) )
            n_{a} = (n_{a} == 0 ? 1 : n_{a}) 
        '''
        
        latt = structure.lattice
        lengths = latt.abc
        ngrid = int(self.kpoints_settings['grid_density']) / structure.num_sites
        
        mult = (ngrid*lengths[0]*lengths[1] * lengths[2]) ** (1/3)
        
        num_div = [int(round(1/lengths[i] * mult)) for i in xrange(3)]
        #ensure that numDiv[i] > 0

        num_div = [i if i > 0 else 1 for i in num_div]
        
        angles = latt.angles
        right_angles = [i for i in xrange(3) if abs(angles[i] - 90) < 5]
        hex_angles = [i for i in xrange(3) if abs(angles[i] - 60) < 5 or abs(angles[i] - 120) < 5]
        
        is_hexagonal = (len(right_angles) == 2 and len(hex_angles) == 1 and lengths[right_angles[0]] == lengths[right_angles[1]])
        
        style = 'Gamma'
        if not is_hexagonal:
            num_div = [i + i % 2 for i in num_div]
            style = 'Monk'
        comment = "pymatgen generated Materials Project kpoints with grid density = " + self.kpoints_settings['grid_density'] + ' per atom.'
        num_kpts = 0
        return Kpoints(comment, num_kpts, style, [num_div], [0,0,0])        
