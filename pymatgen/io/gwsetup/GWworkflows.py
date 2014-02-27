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
import os.path

from pymatgen.io.abinitio.abiobjects import asabistructure
from pymatgen.io.abinitio.calculations import g0w0_extended
from pymatgen.io.abinitio.flows import AbinitFlow
from pymatgen.io.abinitio.tasks import TaskManager
from pymatgen.io.abinitio.pseudos import PseudoTable
from pymatgen.io.gwsetup.GWtasks import *
from pymatgen.io.gwsetup.GWhelpers import now
from fireworks.core.firework import FireWork, Workflow
from fireworks.core.launchpad import LaunchPad

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


class VaspGWFWWorkFlow():
    """
    Object containing a VASP GW workflow for a single structure
    """
    def __init__(self, fw_specs):
        self.work_list = []
        self.connections = {}
        self.fw_id = 1
        self.prep_id = 1
        self.wf = []

    def add_work(self, parameters):
        tasks = []
        job = parameters['job']
        print 'adding job ' + job + ' to the workslist as ', self.fw_id
        if job == 'prep':
            launch_spec = {'task_type': 'Preparation job', '_category': 'cluster', '_queueadapter': 'qadapterdict'}
            task = VaspGWInputTask(parameters)
            tasks.append(task)
            task = VaspGWExecuteTask(parameters)
            tasks.append(task)
            task = VaspGWToDiagTask(parameters)
            tasks.append(task)
            task = VaspGWExecuteTask(parameters)
            tasks.append(task)
            fw = FireWork(tasks, spec=launch_spec, name=job, created_on=now(), fw_id=self.fw_id)
            self.connections[self.fw_id] = []
            self.prep_id = self.fw_id
            self.fw_id += 1
            print self.connections
        elif job in ['G0W0', 'GW0', 'scGW0']:
            launch_spec = {'task_type': 'GW job', '_category': 'cluster', '_queueadapter': 'qadapterdict'}
            task = VaspGWInputTask(parameters)
            tasks.append(task)
            task = VaspGWGetPrepResTask(parameters)
            tasks.append(task)
            task = VaspGWExecuteTask(parameters)
            tasks.append(task)
            if parameters['spec']['converge']:
                task = VaspGWWriteConDatTask(parameters)
                tasks.append(task)
                task = VaspGWTestConTask(parameters)
                tasks.append(task)
            fw = FireWork(tasks, spec=launch_spec, name=job, created_on=now(), fw_id=self.fw_id)
            self.connections[self.fw_id] = []
            self.connections[self.prep_id].append(self.fw_id)
            self.fw_id += 1
        else:
            fw = []
            connection = []
            print 'unspecified job, this should have been captured before !!'
            exit()

        self.work_list.append(fw)

    def create(self):
        self.wf = Workflow(self.work_list, self.connections, name='VaspGWFWWorkFlow', created_on=now())
        print 'creating workflow'

    def add_to_db(self):
        launchpad_file = os.path.join(os.environ['FW_CONFIG_DIR'], 'my_launchpad.yaml')
        lp = LaunchPad.from_file(launchpad_file)
        lp.add_wf(self.wf)


class SingleAbinitGWWorkFlow():
    """
    interface the
    """
    RESPONSE_MODELS = ["cd", "godby", "hybersten", "linden", "farid"]
    TESTS = {'ecuteps': {'test_range': (16, 24), 'method': 'direct', 'control': "gap", 'level': "sigma"},
             'nscf_nbands': {'test_range': (30, 40), 'method': 'set_bands', 'control': "gap", 'level': "nscf"},
             'response_model': {'test_range': RESPONSE_MODELS, 'method': 'direct', 'control': 'gap', 'level': 'screening'}}
    CONVS = {'ecuteps': {'test_range': (8, 16, 24, 32), 'method': 'direct', 'control': "gap", 'level': "sigma"},
             'nscf_nbands': {'test_range': (10, 20, 30, 40, 50, 60, 70), 'method': 'set_bands', 'control': "gap", 'level': "nscf"}}

    def __init__(self, structure, spec, option=None):
        self.structure = structure
        self.spec = spec
        self.option = option
        self.tests = self.__class__.get_defaults_tests()
        self.convs = self.__class__.get_defaults_convs()
        self.response_models = self.__class__.get_response_models()
        self.work_dir = self.structure.composition.reduced_formula
        abi_pseudo = '.GGA_PBE-JTH-paw.xml'
        abi_pseudo_dir = os.path.join(os.environ['ABINIT_PS'], 'GGA_PBE-JTH-paw')
        pseudos = []
        for element in self.structure.composition.element_composition:
            pseudo = os.path.join(abi_pseudo_dir, str(element) + abi_pseudo)
            pseudos.append(pseudo)
        self.pseudo_table = PseudoTable(pseudos)

    @classmethod
    def get_defaults_tests(cls):
        return cls.TESTS.copy()

    @classmethod
    def get_defaults_convs(cls):
        return cls.CONVS.copy()

    @classmethod
    def get_response_models(cls):
        return cls.RESPONSE_MODELS.copy()

    def get_electrons(self, structure):
        """
        Method for retrieving the number of valence electrons
        """
        electrons = 0

        for element in structure.species:
            entries = self.pseudo_table.pseudos_with_symbol(element.symbol)
            assert len(entries) == 1
            pseudo = entries[0]
            electrons += pseudo.Z_val
        return electrons

    def get_bands(self, structure):
        """
        Method for retrieving the standard number of bands
        """
        bands = self.get_electrons(structure) / 2 + len(structure)
        return int(bands)

    def get_work_dir(self):
            name = self.structure.composition.reduced_formula
            if self.option is not None:
                return str(name)+'_'+str(self.option['test'])+'_'+str(self.option['value'])
            else:
                return str(name)

    def create(self):
        """
        create single abinit G0W0 flow
        """
        manager = 'slurm' if 'ceci' in self.spec['mode'] else 'shell'

        abi_structure = asabistructure(self.structure)
        manager = TaskManager.from_user_config()
        # Initialize the flow.
        # FIXME
        # Don't know why protocol=-1 does not work here.
        flow = AbinitFlow(self.work_dir, manager, pickle_protocol=0)

        # kpoint grid defined over density 40 > ~ 3 3 3
        scf_kppa = self.spec.data['kp_grid_dens']
        gamma = True
        # alternatively:
        #nscf_ngkpt = [4,4,4]
        #nscf_shiftk = [0.0, 0.0, 0.0]

        # 100
        nscf_nband = [10 * self.get_bands(self.structure)]
        ecuteps = [8]
        ecutsigx = 8
        ecut = 16

        response_models = ['godby']

        if self.spec['test'] or self.spec['converge']:
            if self.spec['test']:
                tests = SingleAbinitGWWorkFlow(self.structure, self.spec).tests
                response_models = []
            else:
                tests = SingleAbinitGWWorkFlow(self.structure, self.spec).convs
            ecuteps = []
            nscf_nband = []
            for test in tests:
                for value in tests[test]['test_range']:
                    if test == 'nscf_nbands':
                        nscf_nband.append(value * self.get_bands(self.structure))
                        #scr_nband takes nscf_nbands if not specified
                        #sigma_nband takes scr_nbands if not specified
                    if test == 'ecuteps':
                        ecuteps.append(value)
                    if test == 'response_model':
                        response_models.append(value)

        extra_abivars = dict(
            ecut=[ecut],
            istwfk="*1",
            timopt=-1,
            pawecutdg=ecut*2,
            paral_kgb=0,
            nbdbuf=8
        )

        work = g0w0_extended(abi_structure, self.pseudo_table, scf_kppa, nscf_nband, ecuteps, ecutsigx,
                             accuracy="normal", spin_mode="unpolarized", smearing=None, response_models=response_models,
                             charge=0.0, inclvkb=2, sigma_nband=None, scr_nband=None, gamma=gamma,
                             **extra_abivars)

        flow.register_work(work)
        return flow.allocate()

    def create_job_file(self):
        job_file = open("job_collection", mode='a')
        job_file.write('nohup abirun.py ' + self.work_dir + ' scheduler > ' + self.work_dir + '.log & \n')
        job_file.close()
