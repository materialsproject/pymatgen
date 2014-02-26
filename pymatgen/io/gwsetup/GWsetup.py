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
import stat
import os.path
import ast
import pymatgen as pmg

from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen.matproj.rest import MPRester
from pymatgen.serializers.json_coders import MSONable
from pymatgen.io.gwsetup.GWvaspinputsets import SingleVaspGWWork
from pymatgen.io.gwsetup.GWdatastructures import GWConvergenceData
from pymatgen.io.gwsetup.GWworkflows import SingleAbinitGWWorkFlow, VaspGWFWWorkFlow
from pymatgen.io.gwsetup.GWvaspinputsets import MPGWscDFTPrepVaspInputSet, MPGWDFTDiagVaspInputSet, MPGWG0W0VaspInputSet

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

"""
MPGWVaspInputSet.joson contains the standards for GW calculations. This set contains all
parameters for the first sc dft calculation. The modifications for the subsequent
sub calculations are made below.
For many settings the number of cores on which the calculations will be run is needed, this
number is assumed to be on the environment variable NPARGWCALC.
"""


class GWSpecs(MSONable):
    """
    main program
    Class for GW specifications
    """
    def __init__(self):
        self.data = {'mode': 'ceci', 'jobs': ['prep', 'G0W0'], 'test': False, 'source': 'mp-vasp', 'code': 'VASP',
                     'functional': 'PBE', 'kp_grid_dens': 500, 'prec': 'm', 'converge': False}
        self.warnings = []
        self.errors = []
        self.fw_specs = []

    def __getitem__(self, item):
        return self.data[item]

    def to_dict(self):
        return self.data

    def from_dict(self, data):
        self.data = data
        self.test()

    @staticmethod
    def get_npar(self, structure):
        return MPGWG0W0VaspInputSet(structure, self).get_npar(structure)

    def write_to_file(self, filename):
        f = open(filename, mode='w')
        f.write(str(self.to_dict()))
        f.close()

    def read_from_file(self, filename):
        try:
            f = open(filename, mode='r')
            self.data = ast.literal_eval(f.read())
        except OSError:
            print 'Inputfile ', filename, ' not found exiting.'
            exit()

    def reset_job_collection(self):
        if 'ceci' in self.data['mode']:
            if os.path.isfile('job_collection'):
                os.remove('job_collection')
            if 'ABINIT' in self.data['code']:
                job_file = open('job_collection', mode='w')
                job_file.write('module load abinit \n')
                job_file.close()

    def update_interactive(self):
        """
        method to make changes to the GW input setting interactively
        """
        key = 'tmp'
        while len(key) != 0:
            print self.data
            key = raw_input('enter key to change: ')
            if key in self.data.keys():
                value = raw_input('enter new value: ')
                if key == 'jobs':
                    if len(value) == 0:
                        print 'removed', self.data['jobs'].pop(-1)
                    else:
                        self.data['jobs'].append(value)
                elif key == 'test':
                    if value.lower() in ['true', 't']:
                        self.data['test'] = True
                    elif value.lower() in ['false', 'f']:
                        self.data['test'] = False
                    else:
                        print 'undefined value, test should be True or False'
                elif key == 'converge':
                    if value.lower() in ['true', 't']:
                        self.data['converge'] = True
                    elif value.lower() in ['false', 'f']:
                        self.data['converge'] = False
                    else:
                        print 'undefined value, test should be True or False'
                elif key in 'kp_grid_dens':
                    self.data[key] = int(value)
                else:
                    self.data[key] = value
            elif key in ['help', 'h']:
                print "source:       poscar, mp-vasp, any other will be interpreted as a filename to read mp-id's from"
                print "              poscar will read files starting with POSCAR_ in the working folder"
                print 'mode:         input, ceci, fw'
                print 'functional:   PBE, LDA'
                print 'jobs:         prep, G0W0, GW0, scGW0'
                print 'code:         VASP, ABINIT'
                print 'kp_grid_dens: usually 500 - 1000'
                print 'prec:         l, m, h NOT IMPLEMENTED YET'
            elif len(key) == 0:
                print 'setup finished'
            else:
                print 'undefined key'
        self.data['functional'] = self.data['functional'].upper()

    def get_code(self):
        return self['code']

    def test(self):
        if self.data['mode'].lower() not in ['input', 'ceci', 'fw']:
            self.errors.append('unspecified mode')
        if self.data['code'] == 'VASP':
            if self.data['functional'] not in ['PBE', 'LDA']:
                self.errors.append(str(self.data['functional'] + 'not defined for VASP yet'))
        elif self.data['code'] == 'ABINIT':
            if self.data['converge'] and self.data['code'] == 'ABINIT':
                self.warnings.append('converge defined for abinit')
            if self.data['functional'] not in ['PBE']:
                self.errors.append(str(self.data['functional'] + 'not defined for ABINIT yet'))
        else:
            self.errors.append('unknown code')
        if self.data["source"] not in ['poscar', 'mp-vasp']:
            if not os.path.isfile(self.data['source']):
                self.warnings.append('no structures defined')
        if self.data["test"] and self.data["converge"]:
            self.errors.append('both "test" and "converge" are specified, only one can be done at the time')
        if self.data["converge"] and not self.data['mode'] == 'fw':
            self.warnings.append('only real converging in fw mode, for other modes ALL convergence steps are created')
        if len(self.errors) > 0:
            print str(len(self.errors)) + ' error(s) found:'
            print self.errors
            exit()
        if len(self.warnings) > 0:
            print str(len(self.warnings)) + ' warning(s) found:'
            print self.warnings
        self.reset_job_collection()

    def excecute_flow(self, structure):
        """
        excecute spec prepare input/jobfiles or submit to fw for a given structure
        for vasp the different jobs are created into a flow
        for abinit a flow is created using abinitio
        """
        if self.get_code() == 'VASP':
            if self.data['mode'] == 'fw':
                fw_work_flow = VaspGWFWWorkFlow(self.fw_specs)
            else:
                fw_work_flow = []
            if self.data['test'] or self.data['converge']:
                if self.data['test']:
                    tests_prep = MPGWscDFTPrepVaspInputSet(structure, self).tests
                    tests_prep.update(MPGWDFTDiagVaspInputSet(structure, self).tests)
                else:
                    tests_prep = MPGWscDFTPrepVaspInputSet(structure, self).convs
                    tests_prep.update(MPGWDFTDiagVaspInputSet(structure, self).convs)
                for test_prep in tests_prep:
                    print 'setting up test for: ' + test_prep
                    for value_prep in tests_prep[test_prep]['test_range']:
                        print "**" + str(value_prep) + "**"
                        option = {'test_prep': test_prep, 'value_prep': value_prep}
                        self.create_job(structure, 'prep', fw_work_flow, option)
                        for job in self.data['jobs'][1:]:
                            if job == 'G0W0':
                                if self.data['test']:
                                    tests = MPGWG0W0VaspInputSet(structure, self).tests
                                else:
                                    tests = MPGWG0W0VaspInputSet(structure, self).convs
                            if job in ['GW0', 'scGW0']:
                                input_set = MPGWG0W0VaspInputSet(structure, self)
                                input_set.gw0_on()
                                if self.data['test']:
                                    tests = input_set.tests
                                else:
                                    tests = input_set.tests
                            for test in tests:
                                print '    setting up test for: ' + test
                                for value in tests[test]['test_range']:
                                    print "    **" + str(value) + "**"
                                    option.update({'test': test, 'value': value})
                                    self.create_job(structure, job, fw_work_flow, option)
            else:
                for job in self['jobs']:
                    self.create_job(structure, job, fw_work_flow)
            if self.data['mode'] == 'fw':
                fw_work_flow.create()
                fw_work_flow.add_to_db()
        elif self.get_code() == 'ABINIT':
            work_flow = SingleAbinitGWWorkFlow(structure, self)
            flow = work_flow.create()
            flow.build_and_pickle_dump()
            work_flow.create_job_file()
        else:
            print 'unspecified code, actually this should have been catched earlier .. '
            exit()

    def create_job(self, structure, job, fw_work_flow, option=None):
        work = SingleVaspGWWork(structure, job, self, option)
        if 'input' in self.data['mode'] or 'ceci' in self.data['mode']:
            work.create_input()
            if 'ceci' in self.data['mode']:
                work.create_job_script()
        if 'fw' in self.data['mode']:
            structure_dict = structure.to_dict
            band_structure_dict = {'vbm_l': structure.vbm_l, 'cbm_l': structure.cbm_l, 'vbm_a': structure.vbm[0],
                                   'vbm_b': structure.vbm[1], 'vbm_c': structure.vbm[2], 'cbm_a': structure.cbm[0],
                                   'cbm_b': structure.cbm[1], 'cbm_c': structure.cbm[2]}
            parameters = {'structure': structure_dict, 'band_structure': band_structure_dict, 'job': job,
                          'spec': self.to_dict(), 'option': option}
            fw_work_flow.add_work(parameters)

    def process_data(self, structure):
        data = GWConvergenceData(spec=self, structure=structure)
        data.read()
        print data.get_conv_pars()
        data.print_plot_data()

    def loop_structures(self, mode='i'):
        """
        reading the structures specified in spec, add special points, and excecute the specs
        """

        mp_key = os.environ['MP_KEY']

        mp_list_vasp = ['mp-149', 'mp-2534', 'mp-8062', 'mp-2469', 'mp-1550', 'mp-830', 'mp-510626', 'mp-10695', 'mp-66',
                        'mp-1639', 'mp-1265', 'mp-1138', 'mp-23155', 'mp-111']

        if self.data['source'] == 'mp-vasp':
            items_list = mp_list_vasp
        elif self.data['source'] == 'poscar':
            files = os.listdir('.')
            items_list = files
        else:
            items_list = [line.strip() for line in open(self.data['source'])]

        for item in items_list:
            if item.startswith('POSCAR_'):
                structure = pmg.read_structure(item)
                comment = Poscar.from_file(item).comment
                print comment
                if comment.startswith("gap"):
                    structure.vbm_l = comment.split(" ")[1]
                    structure.vbm = (comment.split(" ")[2], comment.split(" ")[3], comment.split(" ")[4])
                    structure.cbm_l = comment.split(" ")[5]
                    structure.cbm = (comment.split(" ")[6], comment.split(" ")[7], comment.split(" ")[8])
                else:
                    print "no bandstructure information available, adding GG as 'gap'"
                    structure.vbm_l = "G"
                    structure.cbm_l = "G"
                    structure.cbm = (0.0, 0.0, 0.0)
                    structure.vbm = (0.0, 0.0, 0.0)
            elif item.startswith('mp-'):
                print item
                with MPRester(mp_key) as mp_database:
                    structure = mp_database.get_structure_by_material_id(item, final=True)
                    bandstructure = mp_database.get_bandstructure_by_material_id(item)
                    structure.vbm_l = bandstructure.kpoints[bandstructure.get_vbm()['kpoint_index'][0]].label
                    structure.cbm_l = bandstructure.kpoints[bandstructure.get_cbm()['kpoint_index'][0]].label
                    structure.cbm = tuple(bandstructure.kpoints[bandstructure.get_cbm()['kpoint_index'][0]].frac_coords)
                    structure.vbm = tuple(bandstructure.kpoints[bandstructure.get_vbm()['kpoint_index'][0]].frac_coords)
            else:
                next(item)
            print structure.composition.reduced_formula
            if mode == 'i':
                self.excecute_flow(structure)
            elif mode == 'o':
                self.process_data(structure)

        if 'ceci' in self.data['mode'] and mode == 'i':
            os.chmod("job_collection", stat.S_IRWXU)

if __name__ == "__main__":
    spec_in = GWSpecs()
    spec_in.update_interactive()
    spec_in.test()
    spec_in.write_to_file('spec.in')
    spec_in.loop_structures('i')