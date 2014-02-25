#!/usr/bin/env python

"""
Script to test writing GW Input for VASP.
Reads the POSCAR_name in the the current folder and outputs GW input to
subfolders name
"""

from __future__ import division

__author__ = "Michiel van Setten"
__copyright__ = " "
__version__ = "0.9"
__maintainer__ = "Michiel van Setten"
__email__ = "mjvansetten@gmail.com"
__date__ = "Oct 23, 2013"

import os

from pymatgen.io.abinitio.netcdf import NetcdfReader
from pymatgen.io.vaspio.vasp_output import Vasprun
from pymatgen.core.units import Ha_to_eV

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


class GWConvergenceData():
    def __init__(self, structure, spec):
        self.structure = structure
        self.spec = spec
        self.data = {}
        self.name = structure.composition.reduced_formula

    def read(self):
        if self.spec['code'] == 'ABINIT':
            read_next = True
            n = 3
            while read_next:
                output = os.path.join(self.name,  'work_0', 'task_' + str(n), 'outdata', 'out_SIGRES.nc')
                if os.path.isfile(output):
                    n += 1
                    data = NetcdfReader(output)
                    data.print_tree()
                    self.data.update({n: {'ecuteps': Ha_to_eV * data.read_value('ecuteps'), 'nbands': data.read_value('sigma_nband'), 'gwgap': data.read_value('egwgap')[0][0]}})
                    data.close()
                else:
                    read_next = False
        elif self.spec['code'] == 'VASP':
            tree = os.walk('.')
            n = 0
            for dirs in tree:
                if "/" + self.name + "." in dirs[0] and ('G0W0' in dirs[0] or 'GW0' in dirs[0] or 'scGW0' in dirs[0]):
                    run = os.path.join(dirs[0], 'vasprun.xml')
                    print run
                    kpoints = os.path.join(dirs[0], 'IBZKPT')
                    if os.path.isfile(run):
                        try:
                            data = Vasprun(run, ionic_step_skip=1)
                            parameters = data.__getattribute__('incar').to_dict
                            bandstructure = data.get_band_structure(kpoints)
                            self.data.update({n: {'ecuteps': parameters['ENCUTGW'], 'nbands': parameters['NBANDS'],
                                              'nomega': parameters['NOMEGA'], 'gwgap': bandstructure.get_band_gap()['energy']}})
                            n += 1
                        except BaseException:
                            pass

    def print_plot_data(self):
        data_list = []
        for k in self.data:
            if self.spec['code'] == 'VASP':
                data_list.append([self.data[k]['nbands'], self.data[k]['ecuteps'], self.data[k]['nomega'], self.data[k]['gwgap']])
            else:
                data_list.append([self.data[k]['nbands'], self.data[k]['ecuteps'], self.data[k]['gwgap']])

        data_file = self.name + '.data'
        f = open(data_file, mode='w')
        tmp = sorted(data_list)[0][0]
        for data in sorted(data_list):
            if tmp != data[0]:
                f.write('\n')
            tmp = data[0]
            if self.spec['code'] == 'VASP':
                f.write(str(data[0]) + ' ' + str(data[1]) + ' ' + str(data[2]) + ' ' + str(data[3]) + '\n')
            else:
                f.write(str(data[0]) + ' ' + str(data[1]) + ' ' + str(data[2]) + '\n')
        f.close()

        '''
        [u'space_group', u'primitive_vectors', u'reduced_symmetry_matrices', u'reduced_symmetry_translations',
         u'atom_species', u'reduced_atom_positions', u'valence_charges', u'atomic_numbers', u'atom_species_names',
         u'chemical_symbols', u'pseudopotential_types', u'symafm', u'kpoint_grid_shift', u'kpoint_grid_vectors',
         u'monkhorst_pack_folding', u'reduced_coordinates_of_kpoints', u'kpoint_weights', u'number_of_electrons',
         u'exchange_functional', u'correlation_functional', u'fermi_energy', u'smearing_scheme', u'smearing_width',
          u'number_of_states', u'eigenvalues', u'occupations', u'tphysel', u'occopt', u'istwfk', u'shiftk',
         u'kptrlatt', u'ngkpt_shiftk', u'ecutwfn', u'ecuteps', u'ecutsigx', u'sigma_nband', u'gwcalctyp',
         u'omegasrdmax', u'deltae', u'omegasrmax', u'scissor_ene', u'usepawu', u'kptgw', u'minbnd', u'maxbnd',
         u'degwgap', u'egwgap', u'en_qp_diago', u'e0', u'e0gap', u'sigxme', u'vxcme', u'vUme', u'degw', u'dsigmee0',
         u'egw', u'eigvec_qp', u'hhartree', u'sigmee', u'sigcmee0', u'sigcmesi', u'sigcme4sd', u'sigxcme4sd',
         u'ze0', u'omega4sd'] for later use
        '''