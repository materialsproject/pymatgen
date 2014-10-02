# coding: utf-8

from __future__ import division, unicode_literals, print_function

"""
GW tasks for fireworks workflows, currently Vasp only
"""


__author__ = "Michiel van Setten"
__copyright__ = " "
__version__ = "0.9"
__maintainer__ = "Michiel van Setten"
__email__ = "mjvansetten@gmail.com"
__date__ = "May 2014"

import os
import os.path
import logging
import shutil
import subprocess
import socket

from pymatgen.core.structure import Structure
from pymatgen.io.gwwrapper.GWvaspinputsets import SingleVaspGWWork
from pymatgen.io.gwwrapper.GWvaspinputsets import GWG0W0VaspInputSet
from fireworks.core.firework import FireTaskBase, FWAction
from fireworks.utilities.fw_serializers import FWSerializable


MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

"""
MPGWVaspInputSet.joson contains the standards for GW calculations. This set contains all
parameters for the first sc dft calculation. The modifications for the subsequent
sub calculations are made below.
For many settings the number of cores on which the calculations will be run is needed, this
number is assumed to be on the environment variable NPARGWCALC.
"""


def get_mpi_runner_ntasks_keyword(mpi_runner):
    if mpi_runner == 'srun':
        return '-n'
    elif mpi_runner == 'mpirun':
        return '-np'
    else:
        raise EnvironmentError('Keyword for <number of tasks to run> is not found '
                               'for mpi_runner {}'.format(mpi_runner))


def get_vasp_environment(calculation_type='standard'):
    """
    Method to get the vasp environmental variables
    :param calculation_type: Type of calculation for the Vasp job ('standard' or 'hse' at the moment)
    :return: Dictionary with the necessary variables (vasp executable, vasp modules, vasp pseudos and mpi_runner)
    on a particular cluster for a particular type of calculation
    :raise ValueError: if the type of calculation is not allowed
    """
    if calculation_type == 'standard':
        try:
            vasp_exe = os.environ['VASP_EXE']
        except KeyError:
            hostname = socket.gethostname()
            raise EnvironmentError('VASP_EXE environment variable on '
                                   'cluster "{hn}" not found'.format(hn=hostname))
        try:
            vasp_modules = os.environ['VASP_MODULES'].split()
        except KeyError:
            hostname = socket.gethostname()
            raise EnvironmentError('VASP_MODULES environment variable on '
                                   'cluster "{hn}" not found'.format(hn=hostname))
        try:
            vasp_mpi_runner = os.environ['VASP_MPI_RUNNER']
        except KeyError:
            hostname = socket.gethostname()
            raise EnvironmentError('VASP_MPI_RUNNER environment variable on '
                                   'cluster "{hn}" not found'.format(hn=hostname))
        try:
            vasp_psp_dir = os.environ['VASP_PSP_DIR']
        except KeyError:
            hostname = socket.gethostname()
            raise EnvironmentError('VASP_PSP_DIR environment variable on '
                                   'cluster "{hn}" not found'.format(hn=hostname))
    elif calculation_type == 'hse':
        try:
            vasp_exe = os.environ['VASP_HSE_EXE']
        except KeyError:
            hostname = socket.gethostname()
            raise EnvironmentError('VASP_HSE_EXE environment variable on '
                                   'cluster "{hn}" not found'.format(hn=hostname))
        try:
            vasp_modules = os.environ['VASP_HSE_MODULES'].split()
        except KeyError:
            hostname = socket.gethostname()
            raise EnvironmentError('VASP_HSE_MODULES environment variable on '
                                   'cluster "{hn}" not found'.format(hn=hostname))
        try:
            vasp_mpi_runner = os.environ['VASP_HSE_MPI_RUNNER']
        except KeyError:
            hostname = socket.gethostname()
            raise EnvironmentError('VASP_HSE_MPI_RUNNER environment variable on '
                                   'cluster "{hn}" not found'.format(hn=hostname))
        try:
            vasp_psp_dir = os.environ['VASP_HSE_PSP_DIR']
        except KeyError:
            hostname = socket.gethostname()
            raise EnvironmentError('VASP_HSE_PSP_DIR environment variable on '
                                   'cluster "{hn}" not found'.format(hn=hostname))
    else:
        raise ValueError('Vasp calculation type "{}" is not allowed '
                         'in get_vasp_environment'.format(calculation_type))
    return {'vasp_exe': vasp_exe, 'vasp_modules': vasp_modules,
            'vasp_mpi_runner': vasp_mpi_runner, 'vasp_psp_dir': vasp_psp_dir,
            'mpi_runner_ntasks_option': get_mpi_runner_ntasks_keyword(vasp_mpi_runner)}


class VaspGWTask(FireTaskBase, FWSerializable):
    """
    Base class for vasp GW tasks
    """
    _fw_name = "Vasp GW task"

    def __init__(self, parameters):
        self.update(parameters)
        structure = Structure.from_dict(parameters['structure'])
        structure.vbm_l = parameters['band_structure']['vbm_l']
        structure.cbm_l = parameters['band_structure']['cbm_l']
        structure.vbm = (parameters['band_structure']['vbm_a'], parameters['band_structure']['vbm_b'], parameters['band_structure']['vbm_c'])
        structure.cbm = (parameters['band_structure']['cbm_a'], parameters['band_structure']['cbm_b'], parameters['band_structure']['cbm_c'])
        self.structure = structure
        self.job = parameters['job']
        self.spec = parameters['spec']
        self.option = parameters['option']

    def get_system(self):
        return self.structure.composition.reduced_formula

    def get_prep_dir(self):
        launch_dir = self.get_system()
        if self.option is not None:
            launch_dir = launch_dir + '.' + self.option['test_prep'] + str(self.option['value_prep'])
        return launch_dir

    def get_launch_dir(self):
        launch_dir = self.get_prep_dir()
        if self.job not in 'prep':
            if 'test' in self.option.keys():
                option_name = '.' + self.option['test'] + str(self.option['value'])
            else:
                option_name = ''
            launch_dir = os.path.join(launch_dir, self.job + option_name)
        return launch_dir



class VaspGWInputTask(VaspGWTask):
    _fw_name = "Vasp GW Input Task"
    """
    class for creating input for GW calculations
    """
    def __init__(self, parameters):
        VaspGWTask.__init__(self, parameters)

    def run_task(self, fw_spec):
        task = SingleVaspGWWork(structure=self.structure, job=self.job, spec=self.spec, option=self.option)
        task.create_input()
        return FWAction()


class VaspGWToDiagTask(VaspGWTask):
    _fw_name = 'Vasp GW To Diag Task'
    """
    class for copying the INCAR.DIAG to INCAR
    """
    def __init__(self, parameters):
        VaspGWTask.__init__(self, parameters)

    def run_task(self, fw_spec):
        shutil.copy2(os.path.join(self.get_prep_dir(), 'INCAR.DIAG'), os.path.join(self.get_prep_dir(), 'INCAR'))
        return FWAction()


class VaspGWGetPrepResTask(VaspGWTask):
    _fw_name = "Vasp GW Get Prep Res Task"
    """
    class for copying the results of the preparation job, the dft calculation, to the GW work folder
    """
    def __init__(self, parameters):
        VaspGWTask.__init__(self, parameters)

    def run_task(self, fw_spec):
        output_files = ['CHGCAR', 'WAVECAR', 'WAVEDER']
        for output_file in output_files:
            shutil.copy2(os.path.join(self.get_prep_dir(), output_file), self.get_launch_dir())
        return FWAction()


class VaspGWExecuteTask(VaspGWTask):
    _fw_name = 'Vasp GW Execute Task'
    """
    class for executing vasp
    """
    def __init__(self, parameters):
        VaspGWTask.__init__(self, parameters)

    def run_task(self, fw_spec):
        name = self.get_system() + self.job
        logging.basicConfig(filename=os.path.join(self.get_launch_dir(), name + '.log'), level=logging.DEBUG)
        vasp_exe = get_vasp_environment()['vasp_exe']
        vasp_mpi_executer = get_vasp_environment()['vasp_mpi_executer']
        n_tasks = GWG0W0VaspInputSet(self.structure).get_npar(self.structure)

        frontend_serial = True
        custodian = False

        if frontend_serial:
            #fake run, no actual vasp execution
            base = os.getcwdu()
            abs_work_dir = os.path.join(base, self.get_launch_dir())
            os.chdir(abs_work_dir)
            cmd = ["/home/setten/scripts/vasp"]
            subprocess.call(cmd)
            os.chdir(base)
        elif custodian:
            mpirunvasp = [vasp_mpi_executer, '-np', n_tasks, vasp_exe]
            pass

        return FWAction()


class VaspGWWriteConDatTask(VaspGWTask):
    fw_name = "Vasp GW Write Convergence Data"
    """
    Write the data needed to study the convergence

    """
    def __init__(self, parameters):
        VaspGWTask.__init__(self, parameters)

    def run_task(self, fw_spec):
        return FWAction()


class VaspGWTestConTask(VaspGWTask):
    fw_name = "Vasp GW Test converge"
    """
    test if the current conv option is converged if not calculate the nex value

    """
    def __init__(self, parameters):
        VaspGWTask.__init__(self, parameters)

    def get_n_pev_runs(self):
        """
        return the number of values that have already been calculated for this convergence parameter
        """
        n = 0
        return n

    def run_task(self, fw_spec):
        task = SingleVaspGWWork(structure=self.structure, job=self.job, spec=self.spec, option=self.option)
        task.create_input()
        return FWAction()
