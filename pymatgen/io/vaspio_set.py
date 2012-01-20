#!/usr/bin/env python

'''
This module defines the VaspInputSet abstract base class and
a concrete implementation for the Materials Project.  The basic
concept behind an input set is to specify a scheme to generate
a consistent set of Vasp inputs from a structure without further
user intervention. This ensures comparability across runs.
'''

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Nov 16, 2011"

import os
import abc
import ConfigParser
import json

from pymatgen.io.vaspio import Incar, Poscar, Potcar, Kpoints

class AbstractVaspInputSet(object):
    """
    Abstract base class representing a set of Vasp input parameters.
    The idea is that using a VaspInputSet, a complete set of input files (INPUT, KPOINTS, POSCAR and POTCAR)
    can be generated in an automated fashion for any structure.
    """
    __metaclass__ = abc.ABCMeta

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
    
    @abc.abstractmethod 
    def get_potcar_symbols(self, structure):
        '''
        Returns Potcar from a structure.
        
        Arguments:
            structure:
                Structure object
        '''
        return
    
    def get_all_vasp_input(self, structure):
        '''
        Returns all input files as a dict of {filename: file_as_string}
        
        Arguments:
            structure:
                Structure object
                
        Returns:
            dict of {filename: file_as_string}, e.g., {'INCAR':'EDIFF=1e-4...'}
        '''
        return {'INCAR':self.get_incar(structure), 'KPOINTS':self.get_kpoints(structure), 
                'POSCAR': self.get_poscar(structure), 'POTCAR': self.get_potcar(structure)}
    
    def write_input(self, structure, output_dir, make_dir_if_not_present = True):
        """
        Writes a set of VASP input to a directory.
        
        Arguments:
            structure: 
                Structure object
            output_dir:
                Directory to output the VASP input files
            make_dir_if_not_present:
                Set to True if you want the directory (and the whole path) to be created if it is not present.
        """
        if make_dir_if_not_present and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        for k,v in self.get_all_vasp_input(structure).items():
            v.write_file(os.path.join(output_dir, k))  

class MITVaspInputSet(AbstractVaspInputSet):
    """
    Standard implementation of VaspInputSet utilizing parameters in the MIT High-throughput project.
    The parameters are chosen specifically for a high-throughput project, which means in general smaller
    pseudopotentials were chosen.
    
    Please refer to A Jain, G. Hautier, C. Moore, S. P. Ong, C. Fischer, T. Mueller, K. A. Persson, G. Ceder (2011). 
    A high-throughput infrastructure for density functional theory calculations. Computational Materials Science, 50(8), 
    2295-2310. doi:10.1016/j.commatsci.2011.02.023 for more information.
    """
    
    def __init__(self):
        self.name = "MITMatgen"
        module_dir = os.path.dirname(os.path.abspath(__file__))
        self._config = ConfigParser.SafeConfigParser()
        self._config.optionxform = str
        self._config.readfp(open(os.path.join(module_dir, "VaspInputSets.cfg")))
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
            elif key == "EDIFF":
                incar[key] =  float(setting) * structure.num_sites
            else:
                incar[key] = setting
                
        has_u = sum(incar['LDAUU']) > 0
        if not has_u:
            for key in incar.keys():
                if key.startswith('LDAU'):
                    del incar[key]
        
        return incar
    
    #get_poscar method inherited from AbstractVaspInputSet

    def get_potcar(self, structure):
        return Potcar(self.get_potcar_symbols(structure))
    
    def get_potcar_symbols(self, structure):               
        p = self.get_poscar(structure)
        elements = p.site_symbols
        potcar_symbols = []
        for el in elements:
            potcar_symbols.append(self.potcar_settings[el] if el in self.potcar_settings else el)
        return potcar_symbols
        
    def get_kpoints(self, structure):
        '''
        Writes out a KPOINTS file using the fully automated grid method. Uses Gamma centered meshes 
        for hexagonal cells and Monkhorst-Pack grids otherwise.
        
        Algorithm: 
            Uses a simple approach scaling the number of divisions along each 
            reciprocal lattice vector proportional to its length. 
        '''
        
        latt = structure.lattice
        lengths = latt.abc
        ngrid = int(self.kpoints_settings['grid_density']) / structure.num_sites
        
        mult = (ngrid*lengths[0]*lengths[1] * lengths[2]) ** (1/3)
        
        num_div = [int(round(1/lengths[i] * mult)) for i in xrange(3)]
        #ensure that numDiv[i] > 0

        num_div = [i if i > 0 else 1 for i in num_div]
        
        angles = latt.angles
        hex_angle_tol = 5 #in degrees
        hex_length_tol = 0.01 #in angstroms
        right_angles = [i for i in xrange(3) if abs(angles[i] - 90) < hex_angle_tol]
        hex_angles = [i for i in xrange(3) if abs(angles[i] - 60) < hex_angle_tol or abs(angles[i] - 120) < hex_angle_tol]
        
        is_hexagonal = (len(right_angles) == 2 and len(hex_angles) == 1 and abs(lengths[right_angles[0]] == lengths[right_angles[1]]) < hex_length_tol)
        
        style = 'Gamma'
        if not is_hexagonal:
            num_div = [i + i % 2 for i in num_div]
            style = 'Monk'
        comment = "pymatgen generated Materials Project kpoints with grid density = " + self.kpoints_settings['grid_density'] + ' per atom.'
        num_kpts = 0
        return Kpoints(comment, num_kpts, style, [num_div], [0,0,0])    
    
    def __str__(self):
        output = [self.name]
        output.append("")
        section_names = ['INCAR settings', 'KPOINTS settings', 'POTCAR settings']
        count = 0
        for d in [self.incar_settings, self.kpoints_settings, self.potcar_settings]:
            output.append(section_names[count])
            for k, v in d.items():
                output.append("%s = %s" % (k, str(v)))
            output.append("")
            count += 1
                
        return "\n".join(output)

class MaterialsProjectVaspInputSet(MITVaspInputSet):
    """
    Implementation of VaspInputSet utilizing parameters in the public Materials Project.
    Typically, the psuedopotentials chosen contain more electrons than the MIT parameters,
    and the k-point grid is ~50% more dense.  The LDAUU parameters are also different 
    due to the different psps used, which result in different fitted values (even though
    the methodology of fitting is exactly the same as the MIT scheme).
    """
    def __init__(self):
        self.name = "MaterialsProject"
        module_dir = os.path.dirname(os.path.abspath(__file__))
        self._config = ConfigParser.SafeConfigParser()
        self._config.optionxform = str
        self._config.readfp(open(os.path.join(module_dir, "VaspInputSets.cfg")))
        self.potcar_settings = dict(self._config.items(self.name + 'POTCAR'))
        self.kpoints_settings = dict(self._config.items(self.name + 'KPOINTS'))
        self.incar_settings = dict(self._config.items(self.name+'INCAR'))
        for key in ['MAGMOM', 'LDAUU', 'LDAUJ', 'LDAUL']:
            self.incar_settings[key] = json.loads(self.incar_settings[key])
