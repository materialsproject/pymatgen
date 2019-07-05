# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License


import f90nml
import numpy as np
import re

"""
This module defines IO ShengBTE for reading, updating, and writing the
CONTROL input file
"""

__author__ = "Rees Chang"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.1"
__email__ = "rc564@cornell.edu"
__date__ = "June 27, 2019"


class Control:

    """
    Class for reading, updating, and writing ShengBTE CONTROL files.
    Currently only supports ShengBTE options relevant to CSLD.
    """

    def __init__(self,
                 alloc_dict={},
                 crystal_dict={},
                 params_dict={},
                 flags_dict={}):
        """
        Args:
            alloc_dict (dict): ShengBTE 'allocations' parameters
            crystal_dict (dict): ShengBTE 'crystal' parameters
            params_dict (dict): ShengBTE 'parameters' parameters
            flags_dict (dict): ShengBTE 'flags' parameters
        """

        self.alloc_dict = {
            'nelements': None,
            'natoms': None,
            'ngrid': None,
            'norientations': 0,
        }
        self.alloc_dict.update(alloc_dict)

        self.crystal_dict = {
            'lfactor': 0.1,
            'lattvec': None, #required
            'types': None, #required
            'elements': None, #required
            'positions': None, #required
            'masses': None,
            'gfactors': None,
            'epsilon': None,
            'born': None,
            'scell': None, #required
            'orientations': None
        }
        self.crystal_dict.update(crystal_dict)

        self.params_dict = {
            'T': 300, #required
            'T_min': None,
            'T_max': None,
            'T_step': None,
            'omega_max': None,
            'scalebroad': 0.5, #required
            'rmin': None,
            'rmax': None,
            'dr': None,
            'maxiter': None,
            'nticks': None,
            'eps': None
        }
        self.params_dict.update(params_dict)

        self.flags_dict = {
            'nonanalytic': None,
            'convergence': None,
            'isotopes': None,
            'autoisotopes': None,
            'nanowires': None,
            'onlyharmonic': None,
            'espresso': None
        }
        self.flags_dict.update(flags_dict)

        def check_required_params():
            """
            Raise error if any required parameters are missing
            """
            required_params = {'crystal':
                                   ['lattvec',
                                    'types',
                                    'elements',
                                    'positions',
                                    'scell'],
                               'parameters':
                                   ['T',
                                    'scalebroad'],
                               }
            required_namelists = list(required_params.keys())
            for required_namelist in required_namelists:
                required_namelist_params = required_params[required_namelist]
                for required_namelist_param in required_namelist_params:
                    if required_namelist == 'allocations' and \
                            self.alloc_dict[required_namelist_param] is None:
                        raise AttributeError('Missing argument: {}>{}'.format(required_namelist,
                                                                              required_namelist_param))
                    elif required_namelist == 'crystal' and \
                            self.crystal_dict[required_namelist_param] is None:
                        raise AttributeError('Missing argument: {}>{}'.format(required_namelist,
                                                                              required_namelist_param))
                    elif required_namelist == 'parameters' and \
                            self.params_dict[required_namelist_param] is None:
                        raise AttributeError('Missing argument: {}>{}'.format(required_namelist,
                                                                              required_namelist_param))
                    elif required_namelist == 'flags' and \
                            self.flags_dict[required_namelist_param] is None:
                        raise AttributeError('Missing argument: {}>{}'.format(required_namelist,
                                                                              required_namelist_param))
        check_required_params()

    @classmethod
    def from_file(cls, filepath):
        """
        Read a CONTROL namelist file and output a 'Control' object

        Args:
            filepath (String): Path of the CONTROL file.

        Returns:
            'Control' object with parameters instantiated.
        """
        nml = f90nml.read(filepath)
        sdict = nml.todict()
        if 't' in sdict['parameters']:
            sdict['parameters']['T'] = sdict['parameters']['t']
            del sdict['parameters']['t']
        if 't_min' in sdict['parameters']:
            sdict['parameters']['T_min'] = sdict['parameters']['t_min']
            del sdict['parameters']['t_min']
        if 't_max' in sdict['parameters']:
            sdict['parameters']['T_max'] = sdict['parameters']['t_max']
            del sdict['parameters']['t_max']
        if 't_step' in sdict['parameters']:
            sdict['parameters']['T_step'] = sdict['parameters']['t_step']
            del sdict['parameters']['t_step']

        alloc_dict = sdict['allocations']
        crystal_dict = sdict['crystal']
        params_dict = sdict['parameters']
        flags_dict = sdict['flags']

        return cls(alloc_dict, crystal_dict, params_dict, flags_dict)

    @classmethod
    def from_dict(cls, sdict):
        """
        Write a CONTROL file from a Python dictionary.
        Description and default parameters can be found at
        https://bitbucket.org/sousaw/shengbte/src/master/.
        Note some parameters are mandatory. Optional parameters
        default here to None and will not be written to file.

        Args:
            dict: A Python dictionary of ShengBTE input parameters.
            filename: Filename to save the CONTROL file
        """

        try:
            alloc_dict = sdict['allocations']
        except:
            alloc_dict = {}
        try:
            crystal_dict = sdict['crystal']
        except:
            crystal_dict = {}
        try:
            params_dict = sdict['parameters']
        except:
            params_dict = {}
        try:
            flags_dict = sdict['flags']
        except:
            flags_dict = {}

        return cls(alloc_dict, crystal_dict, params_dict, flags_dict)

    def to_file(self, filename):
        """
        Writes ShengBTE CONTROL file from 'Control' object
        """
        positions = np.asarray(self.crystal_dict['positions'])
        num_sites, _ = positions.shape

        nelements = str(self.alloc_dict['nelements'])
        natoms = str(self.alloc_dict['natoms'])
        ngrid = self.alloc_dict['ngrid']
        norientations = str(self.alloc_dict['norientations'])

        lfactor = str(self.crystal_dict['lfactor'])
        lattvec1 = self.crystal_dict['lattvec'][0]
        lattvec2 = self.crystal_dict['lattvec'][1]
        lattvec3 = self.crystal_dict['lattvec'][2]
        elements = self.crystal_dict['elements']
        types = self.crystal_dict['types']
        scell = self.crystal_dict['scell']
        # new from here
        if self.crystal_dict['epsilon'] is not None:
            epsilon1 = self.crystal_dict['epsilon'][0]
            epsilon2 = self.crystal_dict['epsilon'][1]
            epsilon3 = self.crystal_dict['epsilon'][2]
        else:
            epsilon1 = np.full(3, None)
            epsilon2 = np.full(3, None)
            epsilon3 = np.full(3, None)
        if self.crystal_dict['born'] is not None:
            born = np.asarray(self.crystal_dict['born'])
        else:
            born = np.full((num_sites, 3, 3), None)
        orientations = np.asarray(self.crystal_dict['orientations'])

        temperature = str(int(self.params_dict['T']))
        scalebroad = str(self.params_dict['scalebroad'])
        t_min = str(self.params_dict['T_min'])
        t_max = str(self.params_dict['T_max'])
        t_step = str(self.params_dict['T_step'])
        omega_max = str(self.params_dict['omega_max'])
        rmin = str(self.params_dict['rmin'])
        rmax = str(self.params_dict['rmax'])
        dr = str(self.params_dict['dr'])
        maxiter = str(self.params_dict['maxiter'])
        nticks = str(self.params_dict['nticks'])
        eps = str(self.params_dict['eps'])

        onlyharmonic = self.flags_dict['onlyharmonic']
        isotopes = self.flags_dict['isotopes']
        nonanalytic = self.flags_dict['nonanalytic']
        nanowires = self.flags_dict['nanowires']
        convergence = self.flags_dict['convergence']
        autoisotopes = self.flags_dict['autoisotopes']
        espresso = self.flags_dict['espresso']

        def boolean_to_string(boolean):
            if boolean is not None:
                if boolean is True:
                    return '.TRUE.'
                else:
                    return '.FALSE.'
            else:
                return 'None'

        #Write strings for types, positions, and born
        indent = '        '
        types_string = 'types='
        positions_string = ''
        born_string = ''
        for line in range(num_sites):
            if line != num_sites-1:
                types_string += str(types[line])+' '
            else:
                types_string += str(types[line])+',\n'

            positions_string += indent+'positions(:,' + str(line+1) + ')=' + str(positions[line,0]) + '  ' \
                                + str(positions[line,1]) + '  ' + str(positions[line,2]) + ',\n'

            for i in range(3):
                born_string += indent+'born(:,'+str(i+1)+','+str(line+1)+')='+str(born[line][i][0])+' '\
                               +str(born[line][i][1])+' '+str(born[line][i][2])+',\n'

        #Write string for orientations
        num_orientations = self.alloc_dict['norientations']
        orientations_string = ''
        for o in range(num_orientations):
            if o != num_orientations-1:
                orientations_string += indent+'orientations(:,'+str(o+1)+')='+str(orientations[o][0])+' '+\
                                       str(orientations[o][1])+' '+str(orientations[o][2])+',\n'
            else:
                orientations_string += indent + 'orientations(:,' + str(o + 1) + ')=' + str(orientations[o][0]) + ' ' + \
                                       str(orientations[o][1]) + ' ' + str(orientations[o][2]) + '\n'

        #masses, gfactors


        full_string = '&allocations\n'+indent+'nelements='+nelements+',\n'
        full_string += indent+'natoms='+natoms+',\n'
        full_string += indent+'ngrid(:)='+str(ngrid[0])+' '+str(ngrid[1])+' '+str(ngrid[2])+'\n'
        full_string += indent+'norientations='+norientations+'\n'

        full_string += '&end\n&crystal\n'
        full_string += indent+'lfactor='+lfactor+',\n'
        full_string += indent+'lattvec(:,1)='+str(lattvec1[0])+'  '+str(lattvec1[1])+'  '+str(lattvec1[2])+',\n'
        full_string += indent+'lattvec(:,2)='+str(lattvec2[0])+'  '+str(lattvec2[1])+'  '+str(lattvec2[2])+',\n'
        full_string += indent+'lattvec(:,3)='+str(lattvec3[0])+'  '+str(lattvec3[1])+'  '+str(lattvec3[2])+',\n'
        full_string += indent+'elements='
        if isinstance(elements, list):
            for i in range(len(elements)):
                full_string += '\"'+elements[i]+str('\"')
                if i != (len(elements)-1):
                    full_string += ' '
                else:
                    full_string += '\n'
        else:
            full_string += '\"'+elements+str('\"\n')
        full_string += indent+types_string
        full_string += positions_string
        full_string += indent+'epsilon(:,1)='+str(epsilon1[0])+' '+str(epsilon1[1])+' '+str(epsilon1[2])+',\n'
        full_string += indent+'epsilon(:,2)='+str(epsilon2[0])+' '+str(epsilon2[1])+' '+str(epsilon2[2])+',\n'
        full_string += indent+'epsilon(:,3)='+str(epsilon3[0])+' '+str(epsilon3[1])+' '+str(epsilon3[2])+',\n'
        full_string += born_string
        full_string += indent+'scell(:)='+str(scell[0])+' '+str(scell[1])+' '+str(scell[2])+'\n'
        full_string += orientations_string

        full_string += '&end\n&parameters\n'
        full_string += indent+'T='+temperature+'\n'
        full_string += indent+'scalebroad='+scalebroad+'\n'
        full_string += indent+'T_min='+t_min+'\n'
        full_string += indent+'T_max='+t_max+'\n'
        full_string += indent+'T_step='+t_step+'\n'
        full_string += indent+'omega_max='+omega_max+'\n'
        full_string += indent+'rmin='+rmin+'\n'
        full_string += indent+'rmax='+rmax+'\n'
        full_string += indent+'dr='+dr+'\n'
        full_string += indent+'maxiter='+maxiter+'\n'
        full_string += indent+'nticks='+nticks+'\n'
        full_string += indent+'eps='+eps+'\n'

        full_string += '&end\n&flags\n'
        full_string += indent+'isotopes='+boolean_to_string(isotopes)+'\n'
        full_string += indent+'onlyharmonic='+boolean_to_string(onlyharmonic)+'\n'
        full_string += indent+'nonanalytic='+boolean_to_string(nonanalytic)+'\n'
        full_string += indent+'nanowires='+boolean_to_string(nanowires)+'\n'
        full_string += indent + 'convergence=' + boolean_to_string(convergence) + '\n'
        full_string += indent + 'autoisotopes=' + boolean_to_string(autoisotopes) + '\n'
        full_string += indent + 'espresso=' + boolean_to_string(espresso) + '\n'
        full_string += '&end'

        def remove_substring(substring, string):
            #Removes lines from 'string' containing 'substring'
            return re.sub('.*'+substring+'.*\n?', '', string)

        full_string = remove_substring('None', full_string)
        file = open(filename, 'w+')
        file.write(full_string)
        file.close()