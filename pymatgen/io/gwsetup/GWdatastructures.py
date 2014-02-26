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
#from numpy import linspace

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


def get_derivatives(xs, ys):
    try:
        from scipy.interpolate import UnivariateSpline
        spline = UnivariateSpline(xs, ys)
        d = spline.derivative(1)
        raise NotImplementedError
    except ImportError:
        d = []
        m, left, right = 0, 0, 0
        for n in range(0, len(xs), 1):
            try:
                left = (ys[n] - ys[n-1]) / (xs[n] - xs[n-1])
                m += 1
            except IndexError:
                pass
            try:
                right = (ys[n+1] - ys[n]) / (xs[n+1] - xs[n])
                m += 1
            except IndexError:
                pass
            d.append(left + right / m)
    return d


def test_conv(xs, ys, tol):
    conv = False
    x_value = float('inf')
    y_value = None
    n_value = None
    ds = get_derivatives(xs, ys)
    for n in range(0, len(ds), 1):
        if abs(ds[n]) < tol:
            conv = True
            if xs[n] < x_value:
                x_value = xs[n]
                y_value = ys[n]
                n_value = n
        else:
            conv = False
            x_value = float('inf')
    return [conv, x_value, y_value, n_value]


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

    def get_var_range(self, var):
        var_range = []
        for data_point in self.data.values():
            if data_point[var] not in var_range:
                var_range.append(data_point[var])
        return sorted(var_range)

    def get_conv_pars(self, tol=0.0001):
        ecuteps_l = False
        nbands_l = False
        ecuteps_c = 0
        nbands_c = 0
        y_conv = []
        z_conv = []
        xs = self.get_var_range('nbands')
        ys = self.get_var_range('ecuteps')
        zd = self.get_data_array()
        for x in xs:
            zs = []
            for y in ys:
                zs.append(zd[x][y])
            conv_data = test_conv(ys, zs, tol)
            if conv_data[0]:
                y_conv.append(conv_data[1])
                z_conv.append(conv_data[2])
                ecuteps_l = conv_data[0]
            else:
                y_conv.append(None)
                z_conv.append(None)
        if ecuteps_l:
            conv_data = test_conv(xs, z_conv, tol)
            print conv_data
            if conv_data[0]:
                nbands_l = conv_data[0]
                nbands_c = conv_data[1]
                ecuteps_c = y_conv[conv_data[3]]
                print "splot '"+self.name+".data' u 1:2:4 w pm3d, '< echo "+'"', nbands_c, ecuteps_c, conv_data[2], '"'+"' w p"

        return {'control': {'ecuteps': ecuteps_l, 'nbands': nbands_l}, 'values': {'ecuteps': ecuteps_c, 'nbands': nbands_c}}

    def get_sorted_data_list(self):
        data_list = []
        for k in self.data:
            if self.spec['code'] == 'VASP':
                data_list.append([self.data[k]['nbands'], self.data[k]['ecuteps'], self.data[k]['nomega'], self.data[k]['gwgap']])
            else:
                data_list.append([self.data[k]['nbands'], self.data[k]['ecuteps'], self.data[k]['gwgap']])
        return sorted(data_list)

    def get_data_array(self):
        data_array = {}
        for k in self.data:
            try:
                data_array[self.data[k]['nbands']].update({self.data[k]['ecuteps']: self.data[k]['gwgap']})
            except KeyError:
                data_array.update({self.data[k]['nbands']: {self.data[k]['ecuteps']: self.data[k]['gwgap']}})
        return data_array

    def print_plot_data(self):
        data_file = self.name + '.data'
        f = open(data_file, mode='w')
        try:
            tmp = self.get_sorted_data_list()[0][0]
        except IndexError:
            tmp = 0
            pass
        for data in self.get_sorted_data_list():
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