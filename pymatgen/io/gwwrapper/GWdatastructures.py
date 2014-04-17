"""
Classes for writing GW Input and analizing GW data. The undelying classe can handle the use of VASP and ABINIT
Reads the POSCAR_name in the the current folder and outputs GW input to subfolders name or lists of structures
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
import math
import pymatgen as pmg

from abc import abstractproperty, abstractmethod, ABCMeta
from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen.matproj.rest import MPRester
from pymatgen.serializers.json_coders import MSONable
from pymatgen.io.gwwrapper.GWvaspinputsets import SingleVaspGWWork
from pymatgen.io.gwwrapper.GWworkflows import SingleAbinitGWWorkFlow, VaspGWFWWorkFlow
from pymatgen.io.gwwrapper.GWvaspinputsets import MPGWscDFTPrepVaspInputSet, MPGWDFTDiagVaspInputSet, MPGWG0W0VaspInputSet
from pymatgen.io.abinitio.netcdf import NetcdfReader
from pymatgen.io.vaspio.vasp_output import Vasprun
from pymatgen.core.units import Ha_to_eV, eV_to_Ha
from pymatgen.io.gwwrapper.GWhelpers import test_conv, print_gnuplot_header, s_name

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

"""
MPGWVaspInputSet.joson contains the standards for GW calculations. This set contains all
parameters for the first sc dft calculation. The modifications for the subsequent
sub calculations are made below.
For many settings the number of cores on which the calculations will be run is needed, this
number is assumed to be on the environment variable NPARGWCALC.
"""


class AbstractAbinitioSprec():
    """

    """

    def __init__(self):
        self.data = {}

    def __getitem__(self, item):
        return self.data[item]

    @abstractproperty
    def help(self):
        """
        property containing a help string for the interactive update
        """
        return str

    @abstractmethod
    def test(self):
        """
        method to test the input for consistency. needs to be implemented by the concrete class
        """

    @abstractmethod
    def excecute_flow(self, structure):
        """
        method called in loopstructures in 'i' input mode, this method schould generate input, job script files etc
         or create fire_work workflows and put them in a database
        """

    @abstractmethod
    def process_data(self, structure):
        """
        method called in loopstructures in 'o' input mode, this method should take the results from the calcultations
        and perform analysis'
        """

    def get_code(self):
        return self['code']

    def to_dict(self):
        return self.data

    def from_dict(self, data):
        self.data = data
        self.test()

    def write_to_file(self, filename):
        f = open(filename, mode='w')
        f.write(str(self.to_dict()))
        f.close()

    def read_from_file(self, filename):
        try:
            f = open(filename, mode='r')
            self.data = ast.literal_eval(f.read())
            f.close()
        except OSError:
            print 'Inputfile ', filename, ' not found exiting.'
            exit()

    def update_interactive(self):
        """
        method to make changes to the GW input setting interactively
        """
        key = 'tmp'
        while len(key) != 0:
            print self
            key = raw_input('enter key to change (h for help): ')
            if key in self.data.keys():
                value = raw_input('enter new value: ')
                if isinstance(self.data[key], list):                        # list
                    if len(value) == 0:
                        print 'removed', self.data[key].pop(-1)
                    else:
                        self.data[key].append(value)
                elif isinstance(self.data[key], bool):                      # logical
                    if value.lower() in ['true', 't']:
                        self.data[key] = True
                    elif value.lower() in ['false', 'f']:
                        self.data[key] = False
                    else:
                        print 'undefined value, should be True or False'
                elif isinstance(self.data[key], int):                       # integer
                    self.data[key] = int(value)
                elif isinstance(self.data[key], float):                     # float
                    self.data[key] = float(value)
                elif isinstance(self.data[key], str):                       # string
                    self.data[key] = value
            elif key in ['help', 'h']:
                print self.help
            elif len(key) == 0:
                print 'setup finished'
            else:
                print 'undefined key'
        self.data['functional'] = self.data['functional'].upper()

    def loop_structures(self, mode='i'):
        """
        reading the structures specified in spec, add special points, and excecute the specs
        mode:
        i: loop structures for input generation
        o: loop structures for output parsing
        """

        mp_key = os.environ['MP_KEY']

        mp_list_vasp = ['mp-149', 'mp-2534', 'mp-8062', 'mp-2469', 'mp-1550', 'mp-830', 'mp-1986', 'mp-10695', 'mp-66',
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
                # print comment
                if comment.startswith("gap"):
                    structure.vbm_l = comment.split(" ")[1]
                    structure.vbm = (comment.split(" ")[2], comment.split(" ")[3], comment.split(" ")[4])
                    structure.cbm_l = comment.split(" ")[5]
                    structure.cbm = (comment.split(" ")[6], comment.split(" ")[7], comment.split(" ")[8])
                else:
                    # print "no bandstructure information available, adding GG as 'gap'"
                    structure.vbm_l = "G"
                    structure.cbm_l = "G"
                    structure.cbm = (0.0, 0.0, 0.0)
                    structure.vbm = (0.0, 0.0, 0.0)
            elif item.startswith('mp-'):
                with MPRester(mp_key) as mp_database:
                    structure = mp_database.get_structure_by_material_id(item, final=True)
                    bandstructure = mp_database.get_bandstructure_by_material_id(item)
                    structure.vbm_l = bandstructure.kpoints[bandstructure.get_vbm()['kpoint_index'][0]].label
                    structure.cbm_l = bandstructure.kpoints[bandstructure.get_cbm()['kpoint_index'][0]].label
                    structure.cbm = tuple(bandstructure.kpoints[bandstructure.get_cbm()['kpoint_index'][0]].frac_coords)
                    structure.vbm = tuple(bandstructure.kpoints[bandstructure.get_vbm()['kpoint_index'][0]].frac_coords)
            else:
                continue
            print "\n", item, s_name(structure)
            if mode == 'i':
                self.excecute_flow(structure)
            elif mode == 'o':
                if os.path.isdir(s_name(structure)) or os.path.isdir(s_name(structure)+'.conv'):
                    self.process_data(structure)

        if 'ceci' in self.data['mode'] and mode == 'i':
            os.chmod("job_collection", stat.S_IRWXU)

    def reset_job_collection(self):
        if 'ceci' in self.data['mode']:
            if os.path.isfile('job_collection'):
                os.remove('job_collection')


class GWSpecs(MSONable, AbstractAbinitioSprec):
    """
    Class for GW specifications.
    """
    def __init__(self):
        self.data = {'mode': 'ceci', 'jobs': ['prep', 'G0W0'], 'test': False, 'source': 'mp-vasp', 'code': 'VASP',
                     'functional': 'PBE', 'kp_grid_dens': 500, 'prec': 'm', 'converge': False, 'tol': 0.0001}
        self.warnings = []
        self.errors = []
        self.fw_specs = []
        self._message = str

    def __str__(self):
        self._message = '%s  %s\n' \
                        '  code         : %s \n' \
                        '  source       : %s \n' \
                        '  jobs         : %s \n' \
                        '  mode         : %s \n' \
                        '  functional   : %s \n' \
                        '  kp_grid_dens : %s \n' \
                        '  prec         : %s \n' \
                        '  tol          : %s \n' \
                        '  test         : %s \n' \
                        '  converge     : %s' % (self.__class__.__name__, self.__doc__, self.get_code(), self.data['source'],
                                                 self['jobs'], self['mode'], self['functional'], self['kp_grid_dens'],
                                                 self['prec'], self['tol'], self['test'], self['converge'])
        return self._message

    @staticmethod
    def get_npar(self, structure):
        return MPGWG0W0VaspInputSet(structure, self).get_npar(structure)

    @property
    def help(self):
        return "source:       poscar, mp-vasp, any other will be interpreted as a filename to read mp-id's from \n" \
               '              poscar will read files starting with POSCAR_ in the working folder \n' \
               'mode:         input, ceci, fw \n' \
               'functional:   PBE, LDA \n' \
               'jobs:         prep, G0W0, GW0, scGW0, no option, just enter, remove the last option \n' \
               'code:         VASP, ABINIT \n' \
               'kp_grid_dens: usually 500 - 1000, 1 gives gamma only, 2 gives a 2x2x2 mesh \n' \
               'tol:          tolerance for determining convergence d(gap)/d(encuteps) and d(gap)/d(nbands) < tol \n' \
               '              for negative values the data is extra polated agains 1/n and covergence wrt the \n' \
               '              complete limit is tested \n' \
               'test:         run all settings defined in test the calculation specific tests \n' \
               'converge:     run all settings defined in convs at low kp mesh, determine converged values, \n' \
               '              rerun all test defined in tests relative to the converged valuues with provided \n' \
               '              kp_grid dens' \
               'prec:         m, h'

    def test(self):
        """
        test if there are inherent errors in the input
        """
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
            self.errors.append('both "test" and "converge" are specified, only one can be done at the same time')
        if self.data["converge"] and not self.data['mode'] == 'fw':
            self.warnings.append('only real converging in fw mode, for other modes ALL convergence steps are created')
        if len(self.errors) > 0:
            print str(len(self.errors)) + ' error(s) found:'
            for error in self.errors:
                print ' > ' + error
            exit()
        if len(self.warnings) > 0:
            print str(len(self.warnings)) + ' warning(s) found:'
            for warning in self.warnings:
                print ' > ' + warning
        self.reset_job_collection()

    def is_converged(self, structure, return_values=False):
        filename = s_name(structure) + ".conv_res"
        try:
            f = open(filename, mode='r')
            conv_res = ast.literal_eval(f.read())
            f.close()
            converged = conv_res['control']['nbands']
        except (IOError, OSError):
            if return_values:
                print 'Inputfile ', filename, ' not found, the convergence calculation did not finish properly or was not' \
                                              ' parsed ...'
            converged = False
        if return_values and converged:
            if self.get_code() == 'ABINIT':
                conv_res['values']['ecuteps'] = 4 * math.ceil(conv_res['values']['ecuteps'] * eV_to_Ha / 4)
            return conv_res['values']
        else:
            return converged

    def get_conv_res_test(self, structure):
        """
        return test sets for the tests in test relative to the convergence results
        """
        tests_conv = {}
        tests_prep_conv = {}
        tests_prep = MPGWscDFTPrepVaspInputSet(structure, self).tests
        tests_prep.update(MPGWDFTDiagVaspInputSet(structure, self).tests)
        tests = MPGWG0W0VaspInputSet(structure, self).tests
        conv_res = self.is_converged(structure, return_values=True)
        for test in conv_res.keys():
            if test in tests_prep.keys():
                rel = tests_prep[test]['test_range'][1] - tests_prep[test]['test_range'][0]
                value = conv_res[test]
                tests_prep_conv.update({test: tests_prep[test]})
                tests_prep_conv[test].update({'test_range': (value, value + rel)})
            elif test in tests.keys():
                rel = tests[test]['test_range'][1] - tests[test]['test_range'][0]
                value = conv_res[test]
                tests_conv.update({test: tests[test]})
                tests_conv[test].update({'test_range': (value, value + rel)})
        return {'tests': tests_conv, 'tests_prep': tests_prep_conv}

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
                elif self.data['converge'] and self.is_converged(structure):
                    tests_prep = self.get_conv_res_test(structure)['tests_prep']
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
                                elif self.data['converge'] and self.is_converged(structure):
                                    tests = self.get_conv_res_test(structure)['tests']
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
            if self['converge'] and self.is_converged(structure):
                option = self.is_converged(structure, return_values=True)
            else:
                option = None
            work_flow = SingleAbinitGWWorkFlow(structure, self, option)
            flow = work_flow.create()
            if flow is not None:
                flow.build_and_pickle_dump()
                work_flow.create_job_file()
        else:
            print 'unspecified code, actually this should have been catched earlier .. '
            exit()

    def create_job(self, structure, job, fw_work_flow, option=None):
        work = SingleVaspGWWork(structure, job, self.data, option=option, converged=self.is_converged(structure))
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
        """
        Process the data of a set of GW calculations:
        for 'single' and 'test' calculations the data is read and outputted
        for the parameter scanning part of a convergence calculation the data is read and parameters that provide
        converged results are determined
        for the 'full' part of a convergence calculation the data is read and it is tested if the slopes are in
        agreement with the scanning part
        """
        data = GWConvergenceData(spec=self, structure=structure)
        if self.data['converge']:
            done = False
            try:
                data.read_full_res_from_file()
                if data.full_res['all_done']:
                    done = True
                    print '| no action needed al is done already'
            except (IOError, OSError, SyntaxError):
                pass

            data.set_type()
            while not done:
                if data.type['parm_scr']:
                    data.read()
                    # determine the parameters that give converged results
                    extrapolated = data.find_conv_pars(self['tol'])
                    # if converged ok, if not increase the grid parameter of the next set of calculations
                    if data.conv_res['control']['nbands']:
                        print '| parm_scr type calculation, converged values found, extrapolated value: ', extrapolated[4]
                    else:
                        print '| parm_scr type calculation, no converged values found, increasing grid'
                        data.full_res['grid'] += 1
                    data.print_full_res()
                    data.print_conv_res()
                    # plot data:
                    print_gnuplot_header('plots', s_name(structure)+' tol = '+str(self['tol']))
                    data.print_gnuplot_line('plots')
                    data.print_plot_data()
                    done = True
                elif data.type['full']:
                    if not data.read_conv_res_from_file(s_name(structure)+'.conv_res'):
                        print '| Full type calculation but the conv_res file is not available, trying to reconstruct'
                        data.read()
                        data.find_conv_pars(self['tol'])
                        data.print_conv_res()
                    data.read(subset='.conv')
                    if len(data.data) == 0:
                        print '| Full type calculation but no data found.'
                        break
                    if data.test_full_kp_results(tol_rel=1, tol_abs=0.001):
                        print '| Full type calculation and the full results agree with the parm_scr. All_done for this compound.'
                        data.full_res.update({'all_done': True})
                        data.print_full_res()
                        done = True
                        data.print_plot_data()
                    else:
                        print '| Full type calculation but the full results do not agree with the parm_scr.'
                        print '|   Increase the tol to find beter converged parameters and test the full grid again.'
                        print '|   TODO'
                        # read the system specific tol for Sytem.conv_res
                        # if it's not there create it from the global tol
                        # reduce tol
                        # set data.type to convergence
                        # loop
                        done = True
                        pass
        elif self.data['test']:
            data.read()
            data.set_type()
            data.print_plot_data()
        else:
            data.read()
            data.set_type()
            data.print_plot_data()


class GWConvergenceData():
    """
    Class for GW data, reading, plotting and testing for convergence
    """
    def __init__(self, structure, spec):
        self.structure = structure
        self.spec = spec
        self.data = {}
        self.conv_res = {'control': {}}
        self.full_res = {'all_done': False, 'grid': 0}
        self.name = s_name(structure)
        self.type = {'parm_scr': False, 'full': False, 'single': False, 'test': False}

    def read_conv_res_from_file(self, filename):
        """
        Read the results of a previous paramert screening set of calculations from file
        """
        success = False
        try:
            f = open(filename, mode='r')
            self.conv_res = ast.literal_eval(f.read())
            f.close()
            success = True
        except SyntaxError:
            print 'Problems reading ', filename
        except (OSError, IOError):
            print 'Inputfile ', filename, ' not found exiting.'
        return success

    def read_full_res_from_file(self):
        """
        Read the results of a full set of calculations from file
        """
        filename = self.name+'.full_res'
        success = False
        try:
            f = open(filename, mode='r')
            self.full_res = ast.literal_eval(f.read())
            f.close()
            success = True
        except SyntaxError:
            print 'Problems reading ', filename
        except (OSError, IOError):
            print 'Inputfile ', filename, ' not found exiting.'
        return success

    def read(self, subset=''):
        """
        Read the G-G gaps from a GW calculation ABINIT / VASP.
        subset: read from System.subset (the . 'dot' should be included in subset)
        Data is placed in self.data, combined with the values of key parameters
        """
        self.data = {}
        if self.spec['code'] == 'ABINIT':
            read_next = True
            n = 3
            while read_next:
                output = os.path.join(self.name + subset,  'work_0', 'task_' + str(n), 'outdata', 'out_SIGRES.nc')
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
                if "/" + self.name + subset + "." in dirs[0] and ('G0W0' in dirs[0] or 'GW0' in dirs[0] or 'scGW0' in dirs[0]):
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

    def set_type(self):
        """
        determine the type of data we are looking at.
        single   : Single calculation with standard parameters
        test     : Test calculation, results of a set of calculations specified in TESTS
        parm_scr : Convergence calculation first part, screening the parameters nbands and ecuteps at a low kp density
        full     : Convergence calculation second part, testing the obtained parameters values at the provided kp density
        """
        name = self.name
        if self.spec['converge']:
            if os.path.isdir(name) and not os.path.isdir(name+'.conv'):
                # convergence setting in spec, but only the low grid dirs exist
                self.type['parm_scr'] = True
            if not os.path.isdir(name) and os.path.isdir(name+'.conv'):
                # test case, separate folder was made for the .conv caluclations
                self.type['full'] = True
            elif os.path.isdir(name) and os.path.isdir(name+'.conv'):
                # both convergence and full dirs exists
                a = max(os.path.getatime(name), os.path.getctime(name), os.path.getmtime(name))
                b = max(os.path.getatime(name+'.conv'), os.path.getctime(name+'.conv'), os.path.getmtime(name+'.conv'))
                if b > a:
                    # full is newer
                    self.type['full'] = True
                elif a > b:
                    # convergence on low grid is newer
                    self.type['parm_scr'] = True
        elif self.spec['test']:
            self.type['test'] = True
        else:
            self.type['single'] = True
        print self.type

    def find_conv_pars(self, tol=0.0001):
        """
        find the pair of smallest values of ecuteps and nbands that give a gamma - gamma gap converged within tol
        positive tol ensures the dirivative is smaller than tol
        negative tol ensures the value is closer than -tol to the assymtotic limmit of a 'A + B / x ^ N' fit
        """
        ecuteps_l = False
        nbands_l = False
        ecuteps_c = 0
        nbands_c = 0
        ecuteps_d = 0
        nbands_d = 0
        gap = None
        y_conv = []
        z_conv = []
        y_conv_der = []
        extrapolated = []
        xs = self.get_var_range('nbands')
        ys = self.get_var_range('ecuteps')
        zd = self.get_data_array()
        print 'plot "'+self.name+'condat'+'"'
        for x in xs:
            zs = []
            for y in ys:
                try:
                    zs.append(zd[x][y])
                except KeyError:
                    pass
            conv_data = test_conv(ys, zs, tol, file_name=self.name+'condat')
            extrapolated.append(conv_data[4])
            if conv_data[0]:
                y_conv.append(conv_data[1])
                y_conv_der.append(conv_data[5])
                z_conv.append(conv_data[2])
                ecuteps_l = conv_data[0]
            else:
                y_conv.append(None)
                z_conv.append(None)
        if ecuteps_l:
            conv_data = test_conv(xs, z_conv, tol, file_name=self.name+'condat')
            if conv_data[0]:
                nbands_l = conv_data[0]
                nbands_c = conv_data[1]
                gap = conv_data[2]
                ecuteps_c = y_conv[conv_data[3]]
                nbands_d = conv_data[5]
                ecuteps_d = y_conv_der[conv_data[3]]
        self.conv_res['control'].update({'ecuteps': ecuteps_l, 'nbands': nbands_l})
        self.conv_res.update({'values': {'ecuteps': ecuteps_c, 'nbands': nbands_c, 'gap': gap},
                              'derivatives': {'ecuteps': ecuteps_d, 'nbands': nbands_d}})
        return test_conv(xs, extrapolated, -1, file_name=self.name+'condat')

    def test_full_kp_results(self, tol_rel=0.5, tol_abs=0.001):
        """
        test if the slopes of the gap data at the full kp mesh are 'comparable' to those of the low kp mesh
        """
        print 'test full kp results'
        self.read_conv_res_from_file(self.name+'.conv_res')
        nbs = self.get_var_range('nbands')
        ecs = self.get_var_range('ecuteps')
        zd = self.get_data_array()
        nb_slope = (zd[nbs[-1]][ecs[-1]] - zd[nbs[0]][ecs[-1]]) / (nbs[-1] - nbs[0])
        ec_slope = (zd[nbs[-1]][ecs[-1]] - zd[nbs[-1]][ecs[0]]) / (ecs[-1] - ecs[0])
        print '        parm_scan   full'
        lnb = abs(nb_slope) < (1 + tol_rel) * abs(self.conv_res['derivatives']['nbands']) and abs(nb_slope) < tol_abs
        print 'nbands  %0.5f     %0.5f %r' % (abs(self.conv_res['derivatives']['nbands']), abs(nb_slope), lnb)
        lec = abs(ec_slope) < (1 + tol_rel) * abs(self.conv_res['derivatives']['ecuteps']) and abs(ec_slope) < tol_abs
        print 'ecuteps %0.5f     %0.5f %r' % (abs(self.conv_res['derivatives']['ecuteps']), abs(ec_slope), lec)
        if lnb and lec:
            return True
        else:
            return False

    def get_sorted_data_list(self):
        data_list = []
        for k in self.data:
            if self.spec['code'] == 'VASP':
                data_list.append([self.data[k]['nbands'], self.data[k]['ecuteps'], self.data[k]['gwgap'], self.data[k]['nomega']])
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

    def get_var_range(self, var):
        var_range = []
        for data_point in self.data.values():
            if data_point[var] not in var_range:
                var_range.append(data_point[var])
        return sorted(var_range)

    def print_gnuplot_line(self, filename):
        """
        print plotting instructions for plotting the G-G gap data
        """
        string1 = "%s%s%s" % ("set output '", self.name, ".jpeg'\n")
        if self.conv_res['control']['nbands']:
            string2 = "%s%s%s%s%s%s%s%s%s%s%s" % ("splot '", self.name, ".data' u 1:2:3 w pm3d, '< echo ", '" ',
                       str(self.conv_res['values']['nbands']), ' ', str(self.conv_res['values']['ecuteps']), ' ',
                       str(self.conv_res['values']['gap']), ' "', "' w p\n")
        else:
            string2 = "%s%s%s" % ("splot '", self.name, ".data' u 1:2:3 w pm3d\n")
        f = open(filename, mode='a')
        f.write(string1)
        f.write(string2)
        f.close()

    def print_plot_data(self):
        """
        print the gap data in a way for 3d plotting using gnuplot to file
        """
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
                f.write('%10.5f %10.5f %10.5f %10.5f \n' % (data[0], data[1], data[2], data[3]))
            else:
                f.write('%10.5f %10.5f %10.5f \n' % (data[0], data[1], data[2]))
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

    def print_conv_res(self):
        """
        print the results of a paramter screening calculation to file
        this file is later used in a subsequent 'full' calculation to perform the calculations at a higher kp-mesh
        """
        # print self.conv_res['control']['nbands']
        if self.conv_res['control']['nbands'] or True:
            filename = self.name + '.conv_res'
            f = open(filename, mode='w')
            string = "{'control': "+str(self.conv_res['control'])+", 'values': "
            if self.spec['code'] == 'VASP':
                string += str({'NBANDS': self.conv_res['values']['nbands'], 'ENCUTGW': self.conv_res['values']['ecuteps']})
                pass
            elif self.spec['code'] == 'ABINIT':
                string += str({'nscf_nbands': self.conv_res['values']['nbands'], 'ecuteps': self.conv_res['values']['ecuteps']})
                pass
            else:
                string = 'undefined code'
            string += ", 'derivatives': "+str(self.conv_res['derivatives'])
            string += '}'
            f.write(string)
            f.close()

    def print_full_res(self):
        """
        print the results of full calculation to fiule
        """
        print self.full_res
        filename = self.name + '.full_res'
        f = open(filename, mode='w')
        string = str(self.full_res)
        f.write(string)
        f.close()