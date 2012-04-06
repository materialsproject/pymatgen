#!/usr/bin/env python

'''
This module defines the VaspInputSet abstract base class and
a concrete implementation for the Materials Project.  The basic
concept behind an input set is to specify a scheme to generate
a consistent set of Vasp inputs from a structure without further
user intervention. This ensures comparability across runs.
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Nov 16, 2011"

import os
import abc
import ConfigParser
import json

from pymatgen.core.periodic_table import Element
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

    def get_all_vasp_input(self, structure, generate_potcar=True):
        '''
        Returns all input files as a dict of {filename: vaspio object}
        
        Arguments:
            structure:
                Structure object
            generate_potcar:
                Set to False to generate a POTCAR.spec file instead of a POTCAR,
                which contains the POTCAR labels but not the actual POTCAR. Defaults
                to True.
                
        Returns:
            dict of {filename: file_as_string}, e.g., {'INCAR':'EDIFF=1e-4...'}
        '''
        if generate_potcar:
            return {'INCAR':self.get_incar(structure), 'KPOINTS':self.get_kpoints(structure),
                    'POSCAR': self.get_poscar(structure), 'POTCAR': self.get_potcar(structure)}
        else:
            return {'INCAR':self.get_incar(structure), 'KPOINTS':self.get_kpoints(structure),
                    'POSCAR': self.get_poscar(structure), 'POTCAR.spec': "\n".join(self.get_potcar_symbols(structure))}

    def write_input(self, structure, output_dir, make_dir_if_not_present=True):
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
        for k, v in self.get_all_vasp_input(structure).items():
            v.write_file(os.path.join(output_dir, k))


class VaspInputSet(AbstractVaspInputSet):
    """
    Standard implementation of VaspInputSet, which can be extended by specific
    implementations.
    """

    def __init__(self, name):
        self.name = name
        module_dir = os.path.dirname(os.path.abspath(__file__))
        self._config = ConfigParser.SafeConfigParser()
        self._config.optionxform = str
        self._config.readfp(open(os.path.join(module_dir, "VaspInputSets.cfg")))
        self.potcar_settings = dict(self._config.items(self.name + 'POTCAR'))
        self.kpoints_settings = dict(self._config.items(self.name + 'KPOINTS'))
        self.incar_settings = dict(self._config.items(self.name + 'INCAR'))
        for key in ['MAGMOM', 'LDAUU', 'LDAUJ', 'LDAUL']:
            if key in self.incar_settings:
                self.incar_settings[key] = json.loads(self.incar_settings[key])

    def get_incar(self, structure):
        incar = Incar()
        comp = structure.composition
        poscar = Poscar(structure)
        for key, setting in self.incar_settings.items():
            if key == "MAGMOM":
                incar[key] = [setting.get(site.specie.symbol, 0.6) for site in structure]
            elif key in ['LDAUU', 'LDAUJ', 'LDAUL']:
                if comp[Element("O")] > 0 or comp[Element("F")] > 0:
                    incar[key] = [setting.get(sym, 0) for sym in poscar.site_symbols]
                else:
                    incar[key] = [0 for sym in poscar.site_symbols]
            elif key == "EDIFF":
                incar[key] = float(setting) * structure.num_sites
            else:
                incar[key] = setting

        has_u = ("LDAUU" in incar and sum(incar['LDAUU']) > 0)
        if has_u:
            # modify LMAXMIX if LSDA+U and you have d or f electrons
            # note that if the user explicitly sets LMAXMIX in settings it will override this logic
            if 'LMAXMIX' not in self.incar_settings.keys():
                if any([el.Z > 56 for el in structure.composition]):  # contains f-electrons
                    incar['LMAXMIX'] = 6
                elif any([el.Z > 20 for el in structure.composition]):  # contains d-electrons
                    incar['LMAXMIX'] = 4
        else:
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
        return Kpoints.automatic_density(structure, int(self.kpoints_settings['grid_density']))

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


class MITVaspInputSet(VaspInputSet):
    """
    Standard implementation of VaspInputSet utilizing parameters in the MIT 
    High-throughput project.
    The parameters are chosen specifically for a high-throughput project, 
    which means in general smaller pseudopotentials were chosen.
    
    Please refer to A Jain, G. Hautier, C. Moore, S. P. Ong, C. Fischer, 
    T. Mueller, K. A. Persson, G. Ceder (2011). 
    A high-throughput infrastructure for density functional theory calculations. 
    Computational Materials Science, 50(8), 2295-2310. 
    doi:10.1016/j.commatsci.2011.02.023 for more information.
    """
    def __init__(self):
        super(MITVaspInputSet, self).__init__("MITMatgen")


class MITHSEVaspInputSet(VaspInputSet):
    """
    Typical implementation of input set for a HSE run.
    """
    def __init__(self):
        super(MITHSEVaspInputSet, self).__init__("MITHSE")


class MaterialsProjectVaspInputSet(VaspInputSet):
    """
    Implementation of VaspInputSet utilizing parameters in the public 
    Materials Project. Typically, the psuedopotentials chosen contain more 
    electrons than the MIT parameters, and the k-point grid is ~50% more dense.  
    The LDAUU parameters are also different due to the different psps used, 
    which result in different fitted values (even though the methodology of 
    fitting is exactly the same as the MIT scheme).
    """
    def __init__(self):
        super(MaterialsProjectVaspInputSet, self).__init__("MaterialsProject")
