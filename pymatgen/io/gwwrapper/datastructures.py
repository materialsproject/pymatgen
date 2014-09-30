# coding: utf-8

from __future__ import unicode_literals, division, print_function

"""
Classes for writing GW Input and analyzing GW data. The underlying classes can handle the use of VASP and ABINIT via the
code interfaces provided in codeinterfaces.
Reads the POSCAR_name in the the current folder and outputs GW input to subfolders name or lists of structures
test 3
"""

__author__ = "Michiel van Setten"
__copyright__ = " "
__version__ = "0.9"
__maintainer__ = "Michiel van Setten"
__email__ = "mjvansetten@gmail.com"
__date__ = "May 2014"

import os
import stat
import os.path
import ast
import pymatgen as pmg
import pymongo
import copy
import gridfs
import six
import numpy as np

from abc import abstractproperty, abstractmethod, ABCMeta
from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen.matproj.rest import MPRester, MPRestError
from pymatgen.serializers.json_coders import MSONable
from pymatgen.io.gwwrapper.convergence import test_conv
from pymatgen.io.gwwrapper.helpers import print_gnuplot_header, s_name, add_gg_gap, refine_structure, now
from pymatgen.io.gwwrapper.codeinterfaces import get_code_interface
from pymatgen.core.structure import Structure
from pymatgen.core.units import eV_to_Ha

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


@six.add_metaclass(ABCMeta)
class AbstractAbinitioSpec(MSONable):
    """
    Contains all non GW specific methods
    todo for some reason I can not make this class have both a metaclass and subcalss from msonable ...
    """

    def __init__(self):
        self.data = {'code': 'ABINIT',
                     'source': 'mp-vasp',
                     'mode': 'ceci',
                     'test': False,
                     'converge': False,
                     'functional': 'PBE',
                     'prec': 'm',
                     'kp_grid_dens': 500,
                     'tol': 0.0001}
        self.warnings = []
        self.errors = []
        self.update_code_interface()

    def __getitem__(self, item):
        return self.data[item]

    def as_dict(self):
        return self.to_dict

    @abstractmethod
    def to_dict(self):
        """
        return a dictionary representation of self
        """

    def update_code_interface(self):
        self.code_interface = get_code_interface(self.get_code())

    def get_code(self):
        return self['code']

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
            print('Inputfile ', filename, ' not found exiting.')
            exit()
        self.update_code_interface()

    def update_interactive(self):
        """
        method to make changes to the GW input setting interactively
        """
        key = 'tmp'
        while len(key) != 0:
            print('Pseudos from: ', self.code_interface.read_ps_dir())
            print(self)
            key = raw_input('enter key to change (h for help): ')
            if key in self.data.keys():
                value = raw_input('enter new value: ')
                if isinstance(self.data[key], list):                        # list
                    if len(value) == 0:
                        print('removed', self.data[key].pop(-1))
                    else:
                        self.data[key].append(value)
                elif isinstance(self.data[key], bool):                      # logical
                    if value.lower() in ['true', 't']:
                        self.data[key] = True
                    elif value.lower() in ['false', 'f']:
                        self.data[key] = False
                    else:
                        print('undefined value, should be True or False')
                elif isinstance(self.data[key], int):                       # integer
                    self.data[key] = int(value)
                elif isinstance(self.data[key], float):                     # float
                    self.data[key] = float(value)
                elif isinstance(self.data[key], str):                       # string
                    self.data[key] = value
            elif key in ['help', 'h']:
                print(self.help)
            elif len(key) == 0:
                print('setup finished')
            else:
                print('undefined key')
        self.data['functional'] = self.data['functional'].upper()
        self.update_code_interface()

    def loop_structures(self, mode='i'):
        """
        reading the structures specified in spec, add special points, and excecute the specs
        mode:
        i: loop structures for input generation
        o: loop structures for output parsing
        w: print all results
        """
        print('loop structures mode ', mode)
        mp_key = os.environ['MP_KEY']

        mp_list_vasp = ['mp-149', 'mp-2534', 'mp-8062', 'mp-2469', 'mp-1550', 'mp-830', 'mp-1986', 'mp-10695', 'mp-66',
                        'mp-1639', 'mp-1265', 'mp-1138', 'mp-23155', 'mp-111']

        if self.data['source'] == 'mp-vasp':
            items_list = mp_list_vasp
        elif self.data['source'] == 'poscar':
            files = os.listdir('.')
            items_list = files
        elif self.data['source'] == 'mar_exp':
            items_list = []
            local_serv = pymongo.Connection("marilyn.pcpm.ucl.ac.be")
            local_db_gaps = local_serv.band_gaps
            pwd = os.environ['MAR_PAS']
            local_db_gaps.authenticate("setten", pwd)
            for c in local_db_gaps.exp.find():
                print(Structure.from_dict(c['icsd_data']['structure']).composition.reduced_formula, c['icsd_id'],\
                    c['MP_id'])
                items_list.append({'name': 'mp-' + c['MP_id'], 'icsd': c['icsd_id'], 'mp': c['MP_id']})
        else:
            items_list = [line.strip() for line in open(self.data['source'])]

        for item in items_list:
            print('\n')
            # special case, this should be encaptulated
            if self.data['source'] == 'mar_exp':
                print('structure from marilyn', item['name'], item['icsd'], item['mp'])
                exp = local_db_gaps.exp.find({'MP_id': item['mp']})[0]
                structure = Structure.from_dict(exp['icsd_data']['structure'])
                structure = refine_structure(structure)
                try:
                    kpts = local_db_gaps.GGA_BS.find({'transformations.history.0.id': item['icsd']})[0]['calculations']\
                    [-1]['band_structure']['kpoints']
                except (IndexError, KeyError):
                    kpts = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
                structure.kpts = kpts
                print('kpoints:', structure.kpts[0], structure.kpts[1])
                structure.item = item['name']
            else:
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
                        structure = add_gg_gap(structure)
                elif 'xyz' in item:
                    structure = pmg.read_structure(item)
                    raise NotImplementedError
                elif item.startswith('mp-'):
                    with MPRester(mp_key) as mp_database:
                        print('structure from mp database', item)
                        structure = mp_database.get_structure_by_material_id(item, final=True)
                        try:
                            bandstructure = mp_database.get_bandstructure_by_material_id(item)
                            structure.vbm_l = bandstructure.kpoints[bandstructure.get_vbm()['kpoint_index'][0]].label
                            structure.cbm_l = bandstructure.kpoints[bandstructure.get_cbm()['kpoint_index'][0]].label
                            structure.cbm = tuple(bandstructure.kpoints[bandstructure.get_cbm()['kpoint_index'][0]].frac_coords)
                            structure.vbm = tuple(bandstructure.kpoints[bandstructure.get_vbm()['kpoint_index'][0]].frac_coords)
                        except (MPRestError, IndexError, KeyError) as err:
                            print(err.message)
                            structure = add_gg_gap(structure)
                else:
                    continue
                structure.kpts = [list(structure.cbm), list(structure.vbm)]
                structure.item = item
            print(item, s_name(structure))
            if mode == 'i':
                self.excecute_flow(structure)
            elif mode == 'w':
                self.print_results(structure)
            elif mode == 's':
                self.insert_in_database(structure)
            elif mode == 'o':
                # if os.path.isdir(s_name(structure)) or os.path.isdir(s_name(structure)+'.conv'):
                self.process_data(structure)

        if 'ceci' in self.data['mode'] and mode == 'i':
            os.chmod("job_collection", stat.S_IRWXU)

    def reset_job_collection(self):
        if 'ceci' in self.data['mode']:
            if os.path.isfile('job_collection'):
                os.remove('job_collection')

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
    def print_results(self, structure):
        """
        method called in loopstructures in 'w' write mode, this method should print final results
        """

    @abstractmethod
    def excecute_flow(self, structure):
        """
        method called in loopstructures in 'i' input mode, this method should generate input, job script files etc
         or create fire_work workflows and put them in a database
        """

    @abstractmethod
    def process_data(self, structure):
        """
        method called in loopstructures in 'o' output mode, this method should take the results from the calcultations
        and perform analysis'
        """

    @abstractmethod
    def insert_in_database(self, structure):
        """
        method called in loopstructures in 's' store mode, this method should take the results from the calcultations
        and put them in a database'
        """


class GWSpecs(AbstractAbinitioSpec):
    """
    Class for GW specifications.
    """
    def __init__(self):
        super(GWSpecs, self).__init__()
        self.data.update({'jobs': ['prep', 'G0W0']})
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
                        '  converge     : %s' % (self.__class__.__name__, self.__doc__, self.get_code(),
                                                 self.data['source'], self['jobs'], self['mode'], self['functional'],
                                                 self['kp_grid_dens'], self['prec'], self['tol'], self['test'],
                                                 self['converge'])
        return self._message

    def __hash__(self):
        """
        this will work as long a there are only hashable items in the data dict
        """
        tmp = copy.deepcopy(self.data)
        tmp['jobs'] = tuple(tmp['jobs'])
        print(tmp)
        return hash(frozenset(tmp.items()))

    def hash_str(self):
        """
        old fist attempt ...
        """
        hashstr = 'code:%s;source:%sjobs:%s;mode:%s;functional:%s;kp_grid_dens:%sprec:%s;tol:%s;test:%s;converge:%s' % \
                  (self.get_code(), self.data['source'], self['jobs'], self['mode'], self['functional'],
                   self['kp_grid_dens'], self['prec'], self['tol'], self['test'], self['converge'])
        return hashstr

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

    def to_dict(self):
        return self.data

    def from_dict(self, data):
        self.data = data
        self.update_code_interface()
        self.test()

    def test(self):
        """
        test if there are inherent errors in the input
        """
        self.warnings, self.errors = self.code_interface.test(self.data)
        if self.data['mode'].lower() not in ['input', 'ceci', 'fw']:
            self.errors.append('unspecified mode')
        if self.data["source"] not in ['poscar', 'mp-vasp']:
            if not os.path.isfile(self.data['source']):
                self.warnings.append('no structures defined')
        if self.data["test"] and self.data["converge"]:
            self.errors.append('both "test" and "converge" are specified, only one can be done at the same time')
        if self.data["converge"] and not self.data['mode'] == 'fw':
            self.warnings.append('only real converging in fw mode, for other modes ALL convergence steps are created')
        if len(self.errors) > 0:
            print(str(len(self.errors)) + ' error(s) found:')
            for error in self.errors:
                print(' > ' + error)
            exit()
        if len(self.warnings) > 0:
            print(str(len(self.warnings)) + ' warning(s) found:')
            for warning in self.warnings:
                print(' > ' + warning)
        self.reset_job_collection()

    def excecute_flow(self, structure):
        """
        excecute spec prepare input/jobfiles or submit to fw for a given structure and the given code interface
        """
        # todo the mode should actually be handeled here... and not inside the code interface
        self.code_interface.excecute_flow(structure, self.data)

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
                    print('| no action needed al is done already')
            except (IOError, OSError, SyntaxError):
                pass

            data.set_type()
            while not done:
                if data.type['parm_scr']:
                    data.read()
                    #print data.data
                    # determine the parameters that give converged results
                    if len(data.data) == 0:
                        print('| parm_scr type calculation but no data found.')
                        break
                    if len(data.data) < 24:
                        print('| parm_scr type calculation but no complete data found,' \
                              ' check is all calculations are done.')
                        break

                    if data.find_conv_pars_scf('ecut', 'full_width', self['tol'])[0]:
                        print('| parm_scr type calculation, converged scf values found')
                        #print data.conv_res
                    else:
                        print('| parm_scr type calculation, no converged scf values found')
                        data.full_res.update({'remark': 'No converged SCf parameter found. '
                                                                             'Solution not implemented.'})
                        data.print_full_res()
                        data.conv_res['values'].update({'ecut': 40*eV_to_Ha})
                        data.conv_res['control'].update({'ecut': True})
                        done = True
                    # if converged ok, if not increase the grid parameter of the next set of calculations
                    extrapolated = data.find_conv_pars(self['tol'])
                    if data.conv_res['control']['nbands']:
                        print('| parm_scr type calculation, converged values found, extrapolated value: %s' %\
                              extrapolated[4])
                    else:
                        print('| parm_scr type calculation, no converged values found, increasing grid')
                        data.full_res['grid'] += 1
                    data.print_full_res()
                    data.print_conv_res()
                    # plot data:
                    print_gnuplot_header('plots', s_name(structure)+' tol = '+str(self['tol']), filetype=None)
                    data.print_gnuplot_line('plots')
                    data.print_plot_data()
                    done = True
                elif data.type['full']:
                    if not data.read_conv_res_from_file(s_name(structure)+'.conv_res'):
                        print('| Full type calculation but the conv_res file is not available, trying to reconstruct')
                        data.read()
                        data.find_conv_pars(self['tol'])
                        data.print_conv_res()
                    data.read(subset='.conv')
                    if len(data.data) == 0:
                        print('| Full type calculation but no data found.')
                        done = True
                    if len(data.data) < 4:
                        print('| Full type calculation but no complete data found.')
                        for item in data.data:
                            print(item)
                        done = True
                    if data.test_full_kp_results(tol_rel=1, tol_abs=0.001):
                        print('| Full type calculation and the full results agree with the parm_scr.' \
                              ' All_done for this compound.')
                        data.full_res.update({'all_done': True})
                        data.print_full_res()
                        done = True
                        #data.print_plot_data()
                        self.code_interface.store_results(name=s_name(structure))
                    else:
                        print('| Full type calculation but the full results do not agree with the parm_scr.')
                        print('|   Increase the tol to find beter converged parameters and test the full grid again.')
                        print('|   TODO')
                        # read the system specific tol for System.conv_res
                        # if it's not there create it from the global tol
                        # reduce tol
                        # set data.type to convergence
                        # loop
                        done = True
        elif self.data['test']:
            data.read()
            data.set_type()
            data.print_plot_data()
        else:
            data.read()
            data.set_type()
            data.print_plot_data()

    def print_results(self, structure, file_name='convergence_results'):
        """
        """
        data = GWConvergenceData(spec=self, structure=structure)
        if data.read_conv_res_from_file(os.path.join(s_name(structure)+'.res', s_name(structure)+'.conv_res')):
            s = '%s %s %s ' % (s_name(structure), str(data.conv_res['values']['ecuteps']), str(data.conv_res['values']['nscf_nbands']))
        else:
            s = '%s 0.0 0.0 ' % s_name(structure)
        con_dat = self.code_interface.read_convergence_data(s_name(structure)+'.res')
        if con_dat is not None:
            s += '%s ' % con_dat['gwgap']
        else:
            s += '0.0 '
        s += '\n'
        f = open(file_name, 'a')
        f.write(s)
        f.close()

    def insert_in_database(self, structure, clean_on_ok=False, db_name='GW_results', collection='general'):
        """
        insert the convergence data and the 'sigres' in a database
        """
        data = GWConvergenceData(spec=self, structure=structure)
        success = data.read_conv_res_from_file(os.path.join(s_name(structure)+'.res', s_name(structure)+'.conv_res'))
        con_dat = self.code_interface.read_convergence_data(s_name(structure)+'.res')
        try:
            f = open('extra_abivars', mode='r')
            extra = ast.literal_eval(f.read())
            f.close()
        except (OSError, IOError):
            extra = None
        ps = self.code_interface.read_ps_dir()
        results_file = os.path.join(s_name(structure)+'.res', self.code_interface.gw_data_file)
        data_file = os.path.join(s_name(structure)+'.res', s_name(structure)+'.data')
        if success and con_dat is not None:
            query = {'system': s_name(structure),
                     'item': structure.item,
                     'spec_hash': hash(self),
                     'extra_vars_hash': hash(None) if extra is None else hash(frozenset(extra.items())),
                     'ps': ps}
            print('query:', query)
            entry = copy.deepcopy(query)
            entry.update({'conv_res': data.conv_res,
                          'spec': self.to_dict(),
                          'extra_vars': extra,
                          'structure': structure.to_dict,
                          'gw_results': con_dat,
                          'results_file': results_file,
                          'data_file': data_file})

            # generic section that should go into the base class like
            #   insert_in_database(query, entry, db_name, collection, server="marilyn.pcpm.ucl.ac.be")
            local_serv = pymongo.Connection("marilyn.pcpm.ucl.ac.be")
            try:
                user = os.environ['MAR_USER']
            except KeyError:
                user = input('DataBase user name: ')
            try:
                pwd = os.environ['MAR_PAS']
            except KeyError:
                pwd = input('DataBase pwd: ')
            db = local_serv[db_name]
            db.authenticate(user, pwd)
            col = db[collection]
            print(col)
            gfs = gridfs.GridFS(db)
            count = col.find(query).count()
            if count == 0:
                try:
                    with open(entry['results_file'], 'r') as f:
                        entry['results_file'] = gfs.put(f.read())
                except IOError:
                    print(entry['results_file'], 'not found')
                try:
                    with open(entry['data_file'], 'r') as f:
                        entry['data_file'] = gfs.put(f.read())
                except IOError:
                    print(entry['data_file'], 'not found')
                col.insert(entry)
                print('inserted', s_name(structure))
            elif count == 1:
                new_entry = col.find_one(query)
                try:
                    print('removing file ', new_entry['results_file'], 'from db')
                    gfs.remove(new_entry['results_file'])
                except:
                    print('remove failed')
                try:
                    print('removing file ', new_entry['data_file'], 'from db')
                    gfs.remove(new_entry['data_file'])
                except:
                    print('remove failed')
                new_entry.update(entry)
                print('adding', new_entry['results_file'], new_entry['data_file'])
                try:
                    with open(new_entry['results_file'], 'r') as f:
                        new_entry['results_file'] = gfs.put(f.read())
                except IOError:
                    print(new_entry['results_file'], 'not found')
                try:
                    with open(new_entry['data_file'], 'r') as f:
                        new_entry['data_file'] = gfs.put(f.read())
                except IOError:
                    print(new_entry['data_file'], 'not found')
                print('as ', new_entry['results_file'], new_entry['data_file'])
                col.save(new_entry)
                print('updated', s_name(structure))
            else:
                print('duplicate entry ... ')
            local_serv.disconnect()
            # end generic section

            #todo remove the workfolders


class GWConvergenceData():
    """
    Class for GW data, reading, plotting and testing for convergence
    """
    def __init__(self, structure, spec):
        self.structure = structure
        self.spec = spec
        self.data = {}
        self.code_interface = get_code_interface(spec['code'])
        self.conv_res = {'control': {}, 'values': {}, 'derivatives': {}}
        self.full_res = {'all_done': False, 'grid': 0}
        self.name = s_name(structure)
        self.type = {'parm_scr': False, 'full': False, 'single': False, 'test': False}

    def read_conv_res_from_file(self, filename):
        """
        Read the results of a previous parameter screening set of calculations from file
        """
        print('reading')
        success = False
        try:
            f = open(filename, mode='r')
            self.conv_res = ast.literal_eval(f.read())
            f.close()
            success = True
            print(self.conv_res)
        except SyntaxError:
            print('Problems reading ', filename)
        except (OSError, IOError):
            print('Inputfile ', filename, ' not found.')
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
            print('Problems reading ', filename)
        except (OSError, IOError):
            print('Inputfile ', filename, ' not found.')
        return success

    def read(self, subset=''):
        """
        Read the G-G gaps from a GW calculation ABINIT / VASP.
        subset: read from System.subset (the . 'dot' should be included in subset)
        Data is placed in self.data, combined with the values of key parameters
        """
        n = 0
        self.data = {}
        tree = os.walk(self.name + subset)
        for dirs in tree:
            read = self.code_interface.read_convergence_data(dirs[0])
            if read:
                self.data.update({n: read})
                n += 1

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
        print(self.type)

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
        for z in zd:
            print(z)
        # print 'plot "'+self.name+'condat'+'"'
        for x in xs:
            zs = []
            for y in ys:
                try:
                    zs.append(zd[x][y])
                except KeyError:
                    pass
            conv_data = test_conv(ys, zs, name=self.name, tol=tol, extra='ecuteps at '+str(x))
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
            conv_data = test_conv(xs, z_conv, name=self.name, tol=tol, extra='nbands')
            if conv_data[0]:
                nbands_l = conv_data[0]
                nbands_c = conv_data[1]
                gap = conv_data[2]
                ecuteps_c = y_conv[conv_data[3]]
                nbands_d = conv_data[5]
                ecuteps_d = y_conv_der[conv_data[3]]
        self.conv_res['control'].update({'ecuteps': ecuteps_l, 'nbands': nbands_l})
        self.conv_res['values'].update({'ecuteps': ecuteps_c, 'nbands': nbands_c, 'gap': gap})
        self.conv_res['derivatives'].update({'ecuteps': ecuteps_d, 'nbands': nbands_d})
        return test_conv(xs, extrapolated, name=self.name, tol=-0.05, extra='nbands at extrapolated ecuteps')

    def find_conv_pars_scf(self, x_name, y_name, tol=0.0001):
        xs = self.get_var_range(x_name)
        ys = []
        #print self.get_data_array_2d(x_name, y_name)
        for x in xs:
            ys.append(self.get_data_array_2d(x_name, y_name)[x])
        conv_data = test_conv(xs, ys, name=self.name, tol=tol, extra=x_name)
        #print conv_data, {x_name: conv_data[0]}, {x_name: conv_data[1]}, {x_name: conv_data[5]}
        self.conv_res['control'].update({x_name: conv_data[0]})
        self.conv_res['values'].update({x_name: conv_data[1]})
        self.conv_res['derivatives'].update({x_name: conv_data[5]})
        return conv_data

    def test_full_kp_results(self, tol_rel=0.5, tol_abs=0.0001):
        """
        test if the slopes of the gap data at the full kp mesh are 'comparable' to those of the low kp mesh
        """
        print('test full kp results')
        self.read_conv_res_from_file(self.name+'.conv_res')
        nbs = self.get_var_range('nbands')
        ecs = self.get_var_range('ecuteps')
        zd = self.get_data_array()
        nb_slope = (zd[nbs[-1]][ecs[-1]] - zd[nbs[0]][ecs[-1]]) / (nbs[-1] - nbs[0])
        ec_slope = (zd[nbs[-1]][ecs[-1]] - zd[nbs[-1]][ecs[0]]) / (ecs[-1] - ecs[0])
        print('        parm_scan   full')
        lnb = abs(nb_slope) < (1 + tol_rel) * abs(self.conv_res['derivatives']['nbands']) or abs(nb_slope) < tol_abs
        print('nbands  %0.5f     %0.5f %r' % (abs(self.conv_res['derivatives']['nbands']), abs(nb_slope), lnb))
        lec = abs(ec_slope) < (1 + tol_rel) * abs(self.conv_res['derivatives']['ecuteps']) or abs(ec_slope) < tol_abs
        print('ecuteps %0.5f     %0.5f %r' % (abs(self.conv_res['derivatives']['ecuteps']), abs(ec_slope), lec))
        print('values: (nb, ec, gap)', nbs[0], ecs[0], zd[nbs[0]][ecs[0]])
        if lnb and lec:
            return True
        else:
            return False

    def get_sorted_data_list(self, data_type='gw_gap'):
        data_list = []
        for k in self.data:
            if data_type == 'gw_gap':
                if 'gwgap' in self.data[k].keys():
                    try:
                        data_list.append([self.data[k]['nbands'],
                                          self.data[k]['ecuteps'],
                                          self.data[k]['gwgap'],
                                          self.data[k]['nomega']])
                    except KeyError:
                        data_list.append([self.data[k]['nbands'],
                                          self.data[k]['ecuteps'],
                                          self.data[k]['gwgap']])
            elif data_type == 'ks_width':
                if 'ecut' in self.data[k].keys():
                    data_list.append([self.data[k]['ecut'],
                                      self.data[k]['min'],
                                      self.data[k]['max']])
        return sorted(data_list)

    def get_data_array(self):
        data_array = {}
        for k in self.data:
            try:
                try:
                    data_array[self.data[k]['nbands']].update({self.data[k]['ecuteps']: self.data[k]['gwgap']})
                except KeyError:
                    data_array.update({self.data[k]['nbands']: {self.data[k]['ecuteps']: self.data[k]['gwgap']}})
            except KeyError:
                pass
        return data_array

    def get_data_array_2d(self, x_name, y_name):
        data_array = {}
        for k in self.data:
            try:
                data_array.update({self.data[k][x_name]: self.data[k][y_name]})
            except KeyError:
                pass
        return data_array

    def get_var_range(self, var):
        var_range = []
        if self.data:
            for data_point in self.data.values():
                try:
                    if data_point[var] not in var_range:
                        var_range.append(data_point[var])
                except KeyError:
                    pass
            return sorted(var_range)
        else:
            return None

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
        with open(filename, mode='a') as f:
            f.write(string1)
            f.write(string2)

    def print_plot_data(self):
        """
        print the gap data in a way for 3d plotting using gnuplot to file
        """
        data_file = self.name + '.data'
        f = open(data_file, mode='w')
        f.write('\n')
        try:
            tmp = self.get_sorted_data_list()[0][0]
        except IndexError:
            tmp = 0
            pass
        for data in self.get_sorted_data_list():
            if tmp != data[0]:
                f.write('\n')
            tmp = data[0]
            try:
                f.write('%10.5f %10.5f %10.5f %10.5f \n' % (data[0], data[1], data[2], data[3]))
            except (KeyError, IndexError):
                f.write('%10.5f %10.5f %10.5f \n' % (data[0], data[1], data[2]))
        f.close()

    def print_conv_res(self):
        """
        print the results of a parameter screening calculation to file
        this file is later used in a subsequent 'full' calculation to perform the calculations at a higher kp-mesh
        """
        if self.conv_res['control']['nbands'] or True:
            filename = self.name + '.conv_res'
            f = open(filename, mode='w')
            string = "{'control': "+str(self.conv_res['control'])+", 'values': "
            string += self.code_interface.conv_res_string(self.conv_res)
            string += ", 'derivatives': "+str(self.conv_res['derivatives'])
            string += '}'
            f.write(string)
            f.close()

    def print_full_res(self):
        """
        print the results of full calculation to fiule
        """
        print(self.full_res)
        filename = self.name + '.full_res'
        f = open(filename, mode='w')
        string = str(self.full_res)
        f.write(string)
        f.close()


# API :


def get_spec(comp_type):
    """
    Main entry point
    """
    return {'GW': GWSpecs()}[comp_type]


def get_convergence_data_structure(comp_type, spec, structure):
    """
    Main entry point
    """
    return {'GW': GWConvergenceData(spec=spec, structure=structure)}[comp_type]
