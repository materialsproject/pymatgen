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
import json
import os.path
import pymatgen as pmg

from pymatgen.io.vaspio.vasp_input import Kpoints
from pymatgen.io.vaspio_set import DictVaspInputSet
from pymatgen.matproj.rest import MPRester
from pymatgen.io.abinitio.abiobjects import asabistructure
from pymatgen.io.abinitio.calculations import g0w0_with_ppmodel
from pymatgen.io.abinitio.flows import AbinitFlow

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

"""
MPGWVaspInputSet.joson contains the standards for GW calculations. This set contains all
parameters for the first sc dft calculation. The modifications for the subsequent
sub calculations are made below.
For many settings the number of cores on which the calculations will be run is needed, this
number is assumed to be on the environment variable NPARGWCALC.
"""


class MPGWscDFTPrepVaspInputSet(DictVaspInputSet):
    """
    Implementation of VaspInputSet overriding MaterialsProjectVaspInputSet
    for static calculations preparing for a GW calculation.
    """
    TESTS = {}

    def __init__(self, structure, **kwargs):
        """
        Supports the same kwargs as :class:`JSONVaspInputSet`.
        """
        self.structure = structure
        with open(os.path.join(MODULE_DIR, "MPGWVaspInputSet.json")) as f:
            DictVaspInputSet.__init__(
                self, "MP Static Self consistent run for GW", json.load(f), **kwargs)
        self.tests = self.__class__.get_defaults_tests()

    #  todo update the fromdict and todict ot include the new atributes

    @classmethod
    def get_defaults_tests(cls):
        return cls.TESTS.copy()

    def get_npar(self, structure):
        """
        get 'optimally' useful number of parallelism
        """
        npar = int(self.get_bands(structure)/2)
        npar = min(max(npar, 2), 52)
        return npar

    def set_test(self, _test_):
        """
        Method to switch a specific test on
        """
        all_tests = MPGWscDFTPrepVaspInputSet.get_defaults_tests()
        all_tests.update(MPGWDFTDiagVaspInputSet.get_defaults_tests())
        all_tests.update(MPGWG0W0VaspInputSet.get_defaults_tests())
        test_type = all_tests[_test_.keys()[0]]['method']
        npar = self.get_npar(self.structure)
        #npar int(os.environ['NPARGWCALC'])
        if test_type == 'incar_settings':
            self.incar_settings.update(_test_)
        if test_type == 'set_nomega':
            nomega = npar * int(_test_['NOMEGA'] / npar)
            self.incar_settings.update({"NOMEGA": int(nomega)})
        if test_type == 'set_nbands':
            #npar = int(os.environ['NPARGWCALC'])
            nbands = _test_['NBANDS'] * self.get_bands(self.structure)
            nbands = npar * int(nbands / npar + 1)
            self.incar_settings.update({"NBANDS": int(nbands)})
        if test_type == 'kpoint_grid':
            pass

    def get_kpoints(self, structure):
        """
        Writes out a KPOINTS file using the automated gamma grid method.
        VASP crashes GW calculations on none gamma centered meshes.
        """
        if self.sort_structure:
            structure = structure.get_sorted_structure()
        dens = int(self.kpoints_settings['grid_density'])
        return Kpoints.automatic_gamma_density(structure, dens)

    def get_bands(self, structure):
        """
        Method for retrieving the standard number of bands
        """
        valence_list = {}
        electrons = 0
        potcar = self.get_potcar(structure)
        for pot_single in potcar:
            valence_list.update({pot_single.element: pot_single.nelectrons})
        #for element in structure.species:
        #    electrons += valence_list[element.symbol]
        electrons = sum([valence_list[element.symbol] for element in structure.species])
        bands = electrons / 2 + len(structure)
        return int(bands)

    def set_test_calc(self):
        """
        absolute minimal setting for testing
        """
        self.incar_settings.update({"PREC": "low", "ENCUT": 250})
        self.kpoints_settings['grid_density'] = 1


class MPGWDFTDiagVaspInputSet(MPGWscDFTPrepVaspInputSet):
    """
    Implementation of VaspInputSet overriding MaterialsProjectVaspInputSet
    for static non self-consistend exact diagonalization step preparing for
    a GW calculation.
    """
    TESTS = {'NBANDS': {'test_range': (10, 15, 20), 'method': 'set_nbands', 'control': "gap"}}

    def __init__(self, structure, **kwargs):
        self.structure = structure
        """
        Supports the same kwargs as :class:`JSONVaspInputSet`.
        """
        with open(os.path.join(MODULE_DIR, "MPGWVaspInputSet.json")) as f:
            DictVaspInputSet.__init__(
                self, "MP Static exact diagonalization", json.load(f), **kwargs)
        self.tests = self.__class__.get_defaults_tests()
        #single step exact diagonalization, output WAVEDER
        self.incar_settings.update({"ALGO": "Exact", "NELM": 1, "LOPTICS": "TRUE"})
        #set the number of bands and npar
        gw_bands = self.get_bands(self.structure)
        #npar = int(os.environ['NPARGWCALC'])
        npar = self.get_npar(self.structure)
        gw_bands = npar * int((10 * gw_bands) / npar + 1)
        # for large systems exact diagonalization consumes too much memory
        if gw_bands > 800:
            self.incar_settings.update({"ALGO": 'fast'})
        self.incar_settings.update({"NBANDS": gw_bands})
        self.incar_settings.update({"NPAR": npar})


class MPGWG0W0VaspInputSet(MPGWDFTDiagVaspInputSet):
    """
    Implementation of VaspInputSet overriding MaterialsProjectVaspInputSet
    for static G0W0 calculation
    """
    TESTS = {'ENCUTGW': {'test_range': (150, 200, 250, 300), 'method': 'incar_settings', 'control': "gap"},
            'NOMEGA': {'test_range': (60, 80, 100), 'method': 'set_nomega', 'control': "gap"}}

    def __init__(self, structure, **kwargs):
        self.structure = structure
        """
        Supports the same kwargs as :class:`JSONVaspInputSet`.
        """
        with open(os.path.join(MODULE_DIR, "MPGWVaspInputSet.json")) as f:
            DictVaspInputSet.__init__(
                self, "MP Static G0W0", json.load(f), **kwargs)
        # G0W0 calculation with reduced cutoff for the response function
        self.incar_settings.update({"ALGO": "GW0", "ENCUTGW": 200, "LWAVE": "FALSE", "NELM": 1})
        # npar = int(os.environ['NPARGWCALC'])
        npar = self.get_npar(structure)
        # set nomega to smallest integer times npar smaller 100
        nomega = npar * int(104 / npar)
        gw_bands = self.get_bands(structure)
        gw_bands = npar * int((10 * gw_bands) / npar + 1)
        self.incar_settings.update({"NBANDS": gw_bands})
        self.incar_settings.update({"NPAR": npar})
        self.incar_settings.update({"NOMEGA": nomega})
        self.tests = self.__class__.get_defaults_tests()

    def spectral_off(self):
        """
        Method to switch the use of the spectral decomposition of the response function of
        this may be used to reduce memory demands if the calculation crashes due to memory shortage
        """
        self.incar_settings.update({"LSPECTRAL": ".False."})

    def gw0_on(self, niter=4, gwbandsfac=4, qpsc=False):
        """
        Method to switch to gw0 calculation with standard 4 iterations
        """
        # set the number of iterations of GW0
        self.incar_settings.update({"NELM": niter})
        # set the number of bands to update in the iteration of G
        npar = self.get_npar(self.structure)
        nbandsgw = self.get_bands(self.structure)*gwbandsfac
        nbandsgw = npar * int(nbandsgw / npar)
        self.incar_settings.update({"NBANDSGW": nbandsgw})
        # if set also updat the orbitals 'quasi particle self-consistency'
        if qpsc:
            self.incar_settings.update({"ALGO": "scGW0"})
        # todo update tests ....


def folder_name(option):
    """
    method to return the sub folder name
    """
    option_prep_name = str(option[0])
    option_name = str(option[1])
    for char in ["'", ":", " ", ",", "{", "}"]:
        option_prep_name = option_prep_name.replace(char, "")
        option_name = option_name.replace(char, "")
    if len(option_prep_name) > 0:
        option_prep_name = "." + option_prep_name
    if len(option_name) > 0:
        option_name = "." + option_name
    return [option_prep_name, option_name]


def create_single_vasp_gw_task(structure, task, spec, option=[{}, {}]):
    """
    Create VASP input for a single standard G0W0 and GW0 calculations
    """
    # create vasp input
    option_prep_name = str(option[0])
    option_name = str(option[1])
    for char in ["'", ":", " ", ",", "{", "}"]:
        option_prep_name = option_prep_name.replace(char, "")
        option_name = option_name.replace(char, "")
    if len(option_prep_name) > 0:
        option_prep_name = "." + option_prep_name
    if len(option_name) > 0:
        option_name = "." + option_name
    if task == 'prep':
        inpset = MPGWscDFTPrepVaspInputSet()
        if spec['test']:
            if option[0].keys()[0] in MPGWscDFTPrepVaspInputSet().tests.keys():
                inpset.set_test(option[0])
        inpset.write_input(structure, structure.composition.reduced_formula+option_prep_name)
        inpset = MPGWDFTDiagVaspInputSet()
        if spec['test']:
            inpset.set_test(option[0])
        inpset.get_incar(structure).write_file(structure.composition.reduced_formula+option_prep_name+'/INCAR.DIAG')
    if task == 'G0W0':
        inpset = MPGWG0W0VaspInputSet()
        if spec['test']:
            inpset.set_test(option[0])
            inpset.set_test(option[1])
        inpset.write_input(structure, structure.composition.reduced_formula+option_prep_name+'/G0W0'+option_name)
    if task == 'GW0':
        inpset = MPGWG0W0VaspInputSet()
        inpset.gw0_on()
        if spec['test']:
            inpset.set_test(option[0])
            inpset.set_test(option[1])
        inpset.write_input(structure, structure.composition.reduced_formula+option_prep_name+'/GW0'+option_name)
    if task == 'scGW0':
        inpset = MPGWG0W0VaspInputSet()
        inpset.gw0_on(qpsc=True)
        if spec['test']:
            inpset.set_test(option[0])
            inpset.set_test(option[1])
        inpset.write_input(structure, structure.composition.reduced_formula+option_prep_name+'/scGW0'+option_name)


def create_single_job_script(structure, task, option=[{}, {}]):
    """
    Create job script for ceci.
    """
    npar = MPGWscDFTPrepVaspInputSet().get_npar(structure)
    option_prep_name = str(option[0])
    option_name = str(option[1])
    for char in ["'", ":", " ", ",", "{", "}"]:
        option_prep_name = option_prep_name.replace(char, "")
        option_name = option_name.replace(char, "")
    if len(option_prep_name) > 0:
        option_prep_name = "." + option_prep_name
    if len(option_name) > 0:
        option_name = "." + option_name
    # npar = int(os.environ['NPARGWCALC'])
    header = ("#!/bin/bash \n"
              "## standard header for Ceci clusters ## \n"
              "#SBATCH --mail-user=michiel.vansetten@uclouvain.be \n"
              "#SBATCH --mail-type=ALL"
              "#SBATCH --time=2-24:0:0 \n"
              "#SBATCH --cpus-per-task=1 \n"
              "#SBATCH --mem-per-cpu=4000 \n")
    if task == 'prep':
        path = structure.composition.reduced_formula + option_prep_name
        # create this job
        job_file = open(name=path+'/job', mode='w')
        job_file.write(header)
        job_file.write('#SBATCH --job-name='+structure.composition.reduced_formula+task+'\n')
        job_file.write('#SBATCH --ntasks='+str(npar)+'\n')
        job_file.write('module load vasp \n')
        job_file.write('mpirun vasp \n')
        job_file.write('cp INCAR.DIAG INCAR \n')
        job_file.write('mpirun vasp \n')
        job_file.close()
        os.chmod(path+'/job', stat.S_IRWXU)
        print 'add this job to the set of jobs'
        job_file = open("job_collection", mode='a')
        job_file.write('cd ' + path + ' \n')
        job_file.write('sbatch job \n')
        job_file.write('cd .. \n')
        job_file.close()
        os.chmod("job_collection", stat.S_IRWXU)
    if task in ['G0W0', 'GW0', 'scGW0']:
        path = structure.composition.reduced_formula + option_prep_name + '/' + task + option_name
        # create this job
        job_file = open(name=path+'/job', mode='w')
        job_file.write(header)
        job_file.write('#SBATCH --job-name='+structure.composition.reduced_formula+task+'\n')
        job_file.write('#SBATCH --ntasks='+str(npar)+'\n')
        job_file.write('module load vasp \n')
        job_file.write('cp ../WAVECAR ../WAVEDER . \n')
        job_file.write('mpirun vasp \n')
        job_file.write('workon pymatgen-GW; get_gap > gap; deactivate')
        job_file.write('echo '+path+'`get_gap` >> ../../gaps.dat')
        job_file.close()
        os.chmod(path+'/job', stat.S_IRWXU)
        path = structure.composition.reduced_formula + option_prep_name
        # 'append submission of this job script to that of prep for this structure'
        job_file = open(name=path+'/job', mode='a')
        job_file.write('cd ' + task + option_name + ' \n')
        job_file.write('sbatch job \n')
        job_file.write('cd .. \n')
        job_file.close()


def create_single_ABINIT_gw_flow(structure, pseudos, manager, work_dir):
    """
    create single abinit G0W0 flow
    """
    abi_structure = asabistructure(structure)

    # Initialize the flow.
    # FIXME
    # Don't know why protocol=-1 does not work here.
    flow = AbinitFlow(work_dir, manager, pickle_protocol=0)

    scf_kppa = 40
    nscf_nband = 100
    #nscf_ngkpt = [4,4,4]
    #nscf_shiftk = [0.0, 0.0, 0.0]
    ecuteps, ecutsigx = 6, 8
    #scr_nband = 50
    #sigma_nband = 50

    extra_abivars = dict(
        ecut=8,
        istwfk="*1",
        timopt=-1,
    )

    work = g0w0_with_ppmodel(abi_structure, pseudos, scf_kppa, nscf_nband, ecuteps, ecutsigx,
                             accuracy="normal", spin_mode="unpolarized", smearing=None,
                             ppmodel="godby", charge=0.0, inclvkb=2, sigma_nband=None, scr_nband=None,
                             **extra_abivars)

    flow.register_work(work)
    return flow.allocate()


def main():
    """
    section for testing. The classes Should eventually be moved to vaspio_set.py
    source: where the structures come from
        poscar: Loop over all files of that match POSCAR_* in the current folder.
        mp: or use structures from hte mp database
    mode:
        input: create all input files locally
        ceci: create tar with input and job submission script
        fw: sumit jobs to fw database
    tasks: prep, G0W0, GW0, scGW0
    code: VASP, ABINIT
    """

    mp_list = ['mp-23818', 'mp-28069', 'mp-23037', 'mp-10627', 'mp-3807', 'mp-252', 'mp-504488', 'mp-1541', 'mp-650026']

    abi_pseudo_dir = os.environ('ABINIT_PS')+'GGA_PBE_FHI'

    spec = {'mode': 'ceci', 'jobs': ['prep'], 'test': False, 'source': 'mp', 'code': 'VASP'}

    if 'ceci' in spec['mode']:
       if os.path.isfile('job_collection'):
           os.remove('job_collection')

    if spec['source'] == 'mp':
        items_list = mp_list
    elif spec['source'] == 'poscar':
        files = os.listdir('.')
        items_list = files

    for item in items_list:
        if item.startswith('POSCAR_'):
            structure = pmg.read_structure(item)
        elif item.startswith('mp-'):
            print item
            with MPRester("wvfSQSO6T3wVdvB9") as mp_database:
                structure = mp_database.get_structure_by_material_id(item, final=False)
        print structure.composition.reduced_formula
        if 'input' or 'ceci' in spec['mode'] and spec['code'] == 'VASP':
            # create all input files locally
            if spec['test']:
                tests_prep = MPGWscDFTPrepVaspInputSet(structure).tests
                tests_prep.update(MPGWDFTDiagVaspInputSet(structure).tests)
                for test_prep in tests_prep:
                    print 'setting up test for: ' + test_prep
                    for value_prep in tests_prep[test_prep]['test_range']:
                        task = 'prep'
                        print "**" + str(value_prep) + "**"
                        create_single_vasp_gw_task(structure, 'prep', option=[{test_prep: value_prep}, {}])
                        if 'ceci' in spec['mode']:
                            create_single_job_script(structure, option=[{test_prep: value_prep}, {}])
                        for task in spec['jobs'][1:]:
                            if task == 'G0W0':
                                tests = MPGWG0W0VaspInputSet(structure).tests
                            if task in ['GW0', 'scGW0']:
                                input_set = MPGWG0W0VaspInputSet(structure)
                                input_set.gw0_on()
                                tests = input_set.tests
                            for test in tests:
                                print '    setting up test for: ' + test
                                for value in tests[test]['test_range']:
                                    print "    **" + str(value) + "**"
                                    create_single_vasp_gw_task(structure, task, spec, option=[{test_prep: value_prep}, {test: value}])
                                    if 'ceci' in spec['mode']:
                                        create_single_job_script(structure, task, option=[{test_prep: value_prep}, {test: value}])
            else:
                for task in spec['jobs']:
                    create_single_vasp_gw_task(structure, task, option=[{}, {}])
                    if 'ceci' in spec['mode']:
                        create_single_job_script()
        if 'ceci' in spec['mode']:
            # create the tar
            pass
        if 'fw' in spec['mode']:
            if spec['test']:
                pass
            else:
                for task in spec['jobs']:
                    pass
                    # submit job to FireWorks database
        if spec['code'] == 'ABINIT':
            # todo based on spec['mode'] set the manager
            manager = 'slurm' if 'ceci' in spec['mode'] else 'shell'
            pseudos = []
            for element in structure.species:
                pseudo = os.path.join(abi_pseudo_dir, element + '.ncps')
                pseudos.append(pseudo)
            work_dir = structure.structure.composition.reduced_formula
            flow = create_single_ABINIT_gw_flow(structure=structure, pseudos=pseudos, work_dir=work_dir, manager=manager)
            flow.build_and_pickle_dump()

if __name__ == "__main__":
    main()