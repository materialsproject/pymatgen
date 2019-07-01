# pymatgen>pymatgen>io

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


class ShengBTE_CONTROL_IO:

    """
    Class for reading, updating, and writing ShengBTE CONTROL files.
    Currently only supports ShengBTE options relevant to CSLD.
    """

    def read_CONTROL(self, filename):
        """
        Read a CONTROL namelist file and output a namelist object
        :param filename: Name of the CONTROL file if in current directory.
            If not, use full path.
        :return: Dictionary of CONTROL parameters.
        """
        nml = f90nml.read(filename)
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
        return sdict

    def file_writer_helper_func(self, dict, filename):
        """
        Helper function, not meant to be called directly
        """
        nelements = str(dict['allocations']['nelements'])
        natoms = str(dict['allocations']['natoms'])
        ngrid = dict['allocations']['ngrid']
        norientations = str(dict['allocations']['norientations'])

        lfactor = str(dict['crystal']['lfactor'])
        lattvec1 = dict['crystal']['lattvec'][0]
        lattvec2 = dict['crystal']['lattvec'][1]
        lattvec3 = dict['crystal']['lattvec'][2]
        elements = dict['crystal']['elements']
        types = dict['crystal']['types']
        positions = np.asarray(dict['crystal']['positions'])
        scell = dict['crystal']['scell']
        # new from here
        epsilon1 = dict['crystal']['epsilon'][0]
        epsilon2 = dict['crystal']['epsilon'][1]
        epsilon3 = dict['crystal']['epsilon'][2]
        born = np.asarray(dict['crystal']['born'])
        orientations = np.asarray(dict['crystal']['orientations'])

        temperature = str(int(dict['parameters']['T']))
        scalebroad = str(dict['parameters']['scalebroad'])
        t_min = str(dict['parameters']['T_min'])
        t_max = str(dict['parameters']['T_max'])
        t_step = str(dict['parameters']['T_step'])
        omega_max = str(dict['parameters']['omega_max'])
        rmin = str(dict['parameters']['rmin'])
        rmax = str(dict['parameters']['rmax'])
        dr = str(dict['parameters']['dr'])
        maxiter = str(dict['parameters']['maxiter'])
        nticks = str(dict['parameters']['nticks'])
        eps = str(dict['parameters']['eps'])

        onlyharmonic = dict['flags']['onlyharmonic']
        isotopes = dict['flags']['isotopes']
        nonanalytic = dict['flags']['nonanalytic']
        nanowires = dict['flags']['nanowires']
        convergence = dict['flags']['convergence']
        autoisotopes = dict['flags']['autoisotopes']
        espresso = dict['flags']['espresso']

        def boolean_to_string(boolean):
            if boolean is not None:
                if boolean is True:
                    return '.TRUE.'
                else:
                    return '.FALSE.'
            else:
                return 'None'

        #Write strings for types, positions, and born
        num_sites, _ = positions.shape
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
        num_orientations = dict['allocations']['norientations']
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

    def write_CONTROL_from_dict(self, dict, filename='CONTROL'):
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
        new_dict = {'allocations':
                            {'nelements': dict.get('allocations', None).get('nelements', None),
                             'natoms': dict.get('allocations', None).get('natoms', None),
                             'ngrid': dict.get('allocations', None).get('ngrid', None),
                             'norientations': dict.get('allocations', 0).get('norientations', 0)},
                    'crystal':
                            {'lfactor': dict.get('crystal', 1.0).get('lfactor', 1.0),
                             'lattvec': dict.get('crystal', 'mandatory').get('lattvec', 'mandatory'),
                             'types': dict.get('crystal', 'mandatory').get('types', 'mandatory'),
                             'elements': dict.get('crystal', 'mandatory').get('elements', 'mandatory'),
                             'positions': dict.get('crystal', 'mandatory').get('positions', 'mandatory'),
                             'masses': dict.get('crystal', None).get('masses', None), #new
                             'gfactors': dict.get('crystal', None).get('gfactors', None), #new
                             'epsilon': dict.get('crystal', None).get('epsilon', None),
                             'born': dict.get('crystal', None).get('born', None),
                             'scell': dict.get('crystal', 'mandatory').get('scell', 'mandatory'),
                             'orientations': dict.get('crystal', None).get('orientations', None)},
                    'parameters':
                            {'T': dict.get('parameters', 'mandatory').get('T', 'mandatory'),
                             'T_min': dict.get('parameters', None).get('T_min', None), #new
                             'T_max': dict.get('parameters', None).get('T_max', None), #new
                             'T_step': dict.get('parameters', None).get('T_step', None), #new
                             'omega_max': dict.get('parameters', None).get('omega_max', None), #new
                             'scalebroad': dict.get('parameters', 1.0).get('scalebroad', 1.0),
                             'rmin': dict.get('parameters', None).get('rmin', None), #new
                             'rmax': dict.get('parameters', None).get('rmax', None), #new
                             'dr': dict.get('parameters', None).get('dr', None), #new
                             'maxiter': dict.get('parameters', None).get('maxiter', None), #new
                             'nticks': dict.get('parameters', None).get('nticks', None), #new
                             'eps': dict.get('parameters', None).get('eps', None)}, #new
                    'flags':
                            {'nonanalytic': dict.get('flags', None).get('nonanalytic', None),
                             'convergence': dict.get('flags', None).get('convergence', None),
                             'isotopes': dict.get('flags', None).get('isotopes', None),
                             'autoisotopes': dict.get('flags', None).get('autoisotopes', None),
                             'nanowires': dict.get('flags', None).get('nanowires', None),
                             'onlyharmonic': dict.get('flags', None).get('onlyharmonic', None),
                             'espresso': dict.get('flags', None).get('espresso', None)}
                    }
        self.file_writer_helper_func(new_dict, filename=filename)


def main():
    io = ShengBTE_CONTROL_IO()
    sdict = io.read_CONTROL('CONTROL')
    print(sdict)
    print(isinstance(sdict, dict))
    print(sdict['crystal']['born'])
    print(sdict['crystal']['born'][1][0][0])

    print(isinstance(sdict['crystal']['types'],list))
    print(sdict['parameters']['T'])
    print(sdict['crystal']['types'])
    print(type(sdict['crystal']['types'][0]))
    print(type(sdict['crystal']['types'][1]))

    io.write_CONTROL_from_dict(sdict, filename='CONTROL_test1')


if __name__ == '__main__':
    main()