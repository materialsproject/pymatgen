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
import logging
import shutil
import subprocess

from pymatgen.io.gwwrapper.GWvaspinputsets import SingleVaspGWWork
from pymatgen.io.gwwrapper.GWvaspinputsets import MPGWG0W0VaspInputSet
from fireworks.core.firework import FWAction
from uclworks.utils.clusters import get_vasp_environment
from uclworks.firetasks.vasptasks import VaspGWTask

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

"""
MPGWVaspInputSet.joson contains the standards for GW calculations. This set contains all
parameters for the first sc dft calculation. The modifications for the subsequent
sub calculations are made below.
For many settings the number of cores on which the calculations will be run is needed, this
number is assumed to be on the environment variable NPARGWCALC.
"""


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
        n_tasks = MPGWG0W0VaspInputSet(self.structure).get_npar(self.structure)

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

