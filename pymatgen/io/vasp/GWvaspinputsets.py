# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

"""
Input sets for VASP GW calculations
Single vasp GW work: Creates input and jobscripts from the input sets for a specific job
"""


__author__ = "Michiel van Setten"
__copyright__ = " "
__version__ = "0.9"
__maintainer__ = "Michiel van Setten"
__email__ = "mjvansetten@gmail.com"
__date__ = "May 2014"

import os
import json
import os.path
import stat

from monty.serialization import loadfn

from pymatgen.io.vasp.inputs import Kpoints, Potcar
from pymatgen.io.vasp.sets import DictVaspInputSet
from pymatgen.io.abinit.helpers import s_name

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

GWVaspInputSet = os.path.join(MODULE_DIR, "GWVaspInputSet.yaml")


"""
MPGWVaspInputSet.joson contains the standards for GW calculations. This set contains all
parameters for the first sc dft calculation. The modifications for the subsequent
sub calculations are made below.
For many settings the number of cores on which the calculations will be run is needed, this
number is assumed to be on the environment variable NPARGWCALC.
"""


class GWscDFTPrepVaspInputSet(DictVaspInputSet):
    """
    Implementation of VaspInputSet overriding MaterialsProjectVaspInputSet
    for static calculations preparing for a GW calculation.
    """
    TESTS = {}
    CONVS = {}

    def __init__(self, structure, spec, functional='PBE', sym_prec=0.01,
                 **kwargs):
        """
        Supports the same kwargs as :class:`JSONVaspInputSet`.
        """
        super(GWscDFTPrepVaspInputSet, self).__init__(
            "MP Static Self consistent run for GW",
            loadfn(GWVaspInputSet), **kwargs)
        self.structure = structure
        self.tests = self.__class__.get_defaults_tests()
        self.convs = self.__class__.get_defaults_convs()
        self.functional = functional
        self.set_dens(spec)
        self.sym_prec = sym_prec
    #  todo update the fromdict and todict ot include the new atributes

    @classmethod
    def get_defaults_tests(cls):
        return cls.TESTS.copy()

    @classmethod
    def get_defaults_convs(cls):
        return cls.CONVS.copy()

    def get_npar(self, structure):
        """
        get 'optimally' useful number of parallelism
        """
        npar = int(self.get_bands(structure) ** 2 * structure.volume / 200)
        npar = min(max(npar, 1), 52)
        return npar

    def set_test(self, test, value):
        """
        Method to switch a specific test on
        """
        all_tests = GWscDFTPrepVaspInputSet.get_defaults_tests()
        all_tests.update(GWDFTDiagVaspInputSet.get_defaults_tests())
        all_tests.update(GWG0W0VaspInputSet.get_defaults_tests())
        test_type = all_tests[test]['method']
        npar = self.get_npar(self.structure)
        if test_type == 'incar_settings':
            self.incar_settings.update({test: value})
        if test_type == 'set_nomega':
            nomega = npar * int(value / npar)
            self.incar_settings.update({"NOMEGA": int(nomega)})
        if test_type == 'set_nbands':
            nbands = value * self.get_bands(self.structure)
            nbands = npar * int(nbands / npar + 1)
            self.incar_settings.update({"NBANDS": int(nbands)})
        if test_type == 'kpoint_grid':
            pass

    def get_potcar(self, structure):
        """
        Method for getting LDA potcars
        """
        if self.sort_structure:
            structure = structure.get_sorted_structure()
        return Potcar(self.get_potcar_symbols(structure), functional=self.functional)

    def get_kpoints(self, structure):
        """
        Writes out a KPOINTS file using the automated gamma grid method.
        VASP crashes GW calculations on none gamma centered meshes.
        """
        if self.sort_structure:
            structure = structure.get_sorted_structure()
        dens = int(self.kpoints_settings['grid_density'])
        if dens == 1:
            return Kpoints.gamma_automatic()
        else:
            return Kpoints.automatic_gamma_density(structure, dens)

    def set_dens(self, spec):
        """
        sets the grid_density to the value specified in spec
        """
        self.kpoints_settings['grid_density'] = spec['kp_grid_dens']
        if spec['kp_grid_dens'] < 100:
            self.incar_settings.update({'ISMEAR': 0})

    def get_electrons(self, structure):
        """
        Method for retrieving the number of valence electrons
        """
        valence_list = {}
        potcar = self.get_potcar(structure)
        for pot_single in potcar:
            valence_list.update({pot_single.element: pot_single.nelectrons})
        electrons = sum([valence_list[element.symbol] for element in structure.species])
        return int(electrons)

    def get_bands(self, structure):
        """
        Method for retrieving the standard number of bands
        """
        bands = self.get_electrons(structure) / 2 + len(structure)
        return int(bands)

    def set_test_calc(self):
        """
        absolute minimal setting for testing
        """
        self.incar_settings.update({"PREC": "low", "ENCUT": 250})
        self.kpoints_settings['grid_density'] = 1

    def set_prec_high(self):
        self.incar_settings.update({"PREC": "Accurate", "ENCUT": 400})


class GWDFTDiagVaspInputSet(GWscDFTPrepVaspInputSet):
    """
    Implementation of VaspInputSet overriding MaterialsProjectVaspInputSet
    for static non self-consistent exact diagonalization step preparing for
    a GW calculation.
    """
    TESTS = {'NBANDS': {'test_range': (10, 20, 30), 'method': 'set_nbands', 'control': "gap"}}
    CONVS = {'NBANDS': {'test_range': (10, 20, 30, 40, 50, 60, 70), 'method': 'set_nbands', 'control': "gap"}}

    def __init__(self, structure, spec, functional='PBE', sym_prec=0.01, **kwargs):
        """
        Supports the same kwargs as :class:`JSONVaspInputSet`.
        """
        super(GWDFTDiagVaspInputSet, self).__init__(
            structure, spec, functional=functional, sym_prec=sym_prec, **kwargs)
        self.structure = structure
        self.tests = self.__class__.get_defaults_tests()
        self.convs = self.__class__.get_defaults_convs()
        self.functional = functional
        self.sym_prec = sym_prec
        self.set_dens(spec)
        npar = self.get_npar(self.structure)
        #single step exact diagonalization, output WAVEDER
        self.incar_settings.update({"ALGO": "Exact", "NELM": 1, "LOPTICS": "TRUE"})
        # for large systems exact diagonalization consumes too much memory
        self.set_gw_bands(15)
        self.incar_settings.update({"NPAR": npar})

    def set_gw_bands(self, factor=15):
        """
        method to set the number of bands for GW
        """
        gw_bands = self.get_bands(self.structure)
        gw_bands = self.get_npar(self.structure) * int((factor * gw_bands) / self.get_npar(self.structure) + 1)
        self.incar_settings.update({"NBANDS": gw_bands})
        if gw_bands > 800:
            self.incar_settings.update({"ALGO": 'fast'})

    def set_prec_high(self):
        super(GWDFTDiagVaspInputSet, self).set_prec_high()
        self.set_gw_bands(30)


class GWG0W0VaspInputSet(GWDFTDiagVaspInputSet):
    """
    Should go to Pymatgen vaspinputsets
    Implementation of VaspInputSet overriding MaterialsProjectVaspInputSet
    for static G0W0 calculation
    """
    TESTS = {'ENCUTGW': {'test_range': (200, 300, 400),
                         'method': 'incar_settings', 'control': "gap"},
             'NOMEGA': {'test_range': (80, 100, 120),
                        'method': 'set_nomega', 'control': "gap"}}
    CONVS = {'ENCUTGW': {'test_range': (200, 400, 600, 800),
                         'method': 'incar_settings', 'control': "gap"}}

    def __init__(self, structure, spec, functional='PBE', sym_prec=0.01,
                 **kwargs):
        """
        Supports the same kwargs as :class:`JSONVaspInputSet`.
        """
        super(GWG0W0VaspInputSet, self).__init__(
            structure, spec, functional=functional, sym_prec=sym_prec, **kwargs)
        self.structure = structure
        self.tests = self.__class__.get_defaults_tests()
        self.convs = self.__class__.get_defaults_convs()
        self.functional = functional
        self.sym_prec = sym_prec
        npar = self.get_npar(structure)
        # G0W0 calculation with reduced cutoff for the response function
        self.incar_settings.update({"ALGO": "GW0", "ENCUTGW": 250,
                                    "LWAVE": "FALSE", "NELM": 1})
        self.set_dens(spec)
        self.nomega_max = 2 * self.get_kpoints(structure).kpts[0][0]**3
        nomega = npar * int(self.nomega_max / npar)
        self.set_gw_bands(15)
        self.incar_settings.update({"NPAR": npar})
        self.incar_settings.update({"NOMEGA": nomega})
        self.tests = self.__class__.get_defaults_tests()

    def wannier_on(self):
        self.incar_settings.update({"LWANNIER90_RUN": ".TRUE."})
        self.incar_settings.update({"LWRITE_MMN_AMN": ".TRUE."})

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

    def set_prec_high(self):
        super(GWG0W0VaspInputSet, self).set_prec_high()
        self.incar_settings.update({"ENCUTGW": 400, "NOMEGA": int(self.incar_settings["NOMEGA"]*1.5)})
        self.incar_settings.update({"PRECFOCK": "accurate"})


class Wannier90InputSet():
    """
    class containing the input parameters for the wannier90.win file
    """
    def __init__(self, spec):
        self.file_name = "wannier90.win"
        self.settings = {"bands_plot": "true", "num_wann": 2, "num_bands": 4}
        self.parameters = {"n_include_bands": 1}
        self.spec = spec

    def make_kpoint_path(self, structure, f):
        f.write("\nbegin kpoint_path\n")
        line = str(structure.vbm_l) + " " + str(structure.vbm[0]) + " " + str(structure.vbm[1]) + " " + str(structure.vbm[2])
        line = line + " " + str(structure.cbm_l) + " " + str(structure.cbm[0]) + " " + str(structure.cbm[1]) + " " + str(structure.cbm[2])
        f.write(line)
        f.write("\nend kpoint_path\n\n")
        pass

    def make_exclude_bands(self, structure, f):
        nocc = GWscDFTPrepVaspInputSet(structure, self.spec).get_electrons(structure) / 2
        n1 = str(int(1))
        n2 = str(int(nocc - self.parameters["n_include_bands"]))
        n3 = str(int(nocc + 1 + self.parameters["n_include_bands"]))
        n4 = str(int(GWG0W0VaspInputSet(structure, self.spec).incar_settings["NBANDS"]))
        line = "exclude_bands : " + n1 + "-" + n2 + ", " + n3 + "-" + n4 + "\n"
        f.write(line)
        #todo there is still a bug here...
        pass

    def write_file(self, structure, path):
        f = open(os.path.join(path, self.file_name), mode='w')
        f.write("bands_plot = ")
        f.write(self.settings["bands_plot"])
        f.write("\n")
        self.make_kpoint_path(structure, f)
        f.write("num_wann  = ")
        f.write(str(self.settings["num_wann"]))
        f.write("\n")
        f.write("num_bands = ")
        f.write(str(self.settings["num_bands"]))
        f.write("\n")
        self.make_exclude_bands(structure, f)
        f.close()


class SingleVaspGWWork():
    """
    Create VASP input for a single standard G0W0 and GW0 calculation step
    the combination of job and option specifies what needs to be created
    """
    def __init__(self, structure, job, spec, option=None, converged=False):
        self.structure = structure
        self.job = job
        self.spec = spec
        self.option = option
        self.converged = converged

    def create_input(self):
        """
        create vasp input
        """
        option_name = ''
        path_add = ''
        if self.spec['converge'] and self.converged:
            path_add = '.conv'
        if self.option is None:
            path = s_name(self.structure)
        else:
            path = os.path.join(s_name(self.structure) + path_add,
                                str(self.option['test_prep'])+str(self.option['value_prep']))
            if 'test' in self.option.keys():
                option_name = '.'+str(self.option['test'])+str(self.option['value'])
        if self.job == 'prep':

            inpset = GWscDFTPrepVaspInputSet(self.structure, self.spec, functional=self.spec['functional'])
            if self.spec['converge'] and not self.converged:
                spec_tmp = self.spec.copy()
                spec_tmp.update({'kp_grid_dens': 2})
                inpset = GWscDFTPrepVaspInputSet(self.structure, spec_tmp, functional=self.spec['functional'])
                inpset.incar_settings.update({"ENCUT": 800})
            if self.spec['test'] or self.spec['converge']:
                if self.option['test_prep'] in GWscDFTPrepVaspInputSet.get_defaults_convs().keys() or self.option['test_prep'] in GWscDFTPrepVaspInputSet.get_defaults_tests().keys():
                    inpset.set_test(self.option['test_prep'], self.option['value_prep'])
            if self.spec["prec"] == "h":
                inpset.set_prec_high()
            inpset.write_input(self.structure, path)

            inpset = GWDFTDiagVaspInputSet(self.structure, self.spec, functional=self.spec['functional'])
            if self.spec["prec"] == "h":
                inpset.set_prec_high()
            if self.spec['converge'] and not self.converged:
                spec_tmp = self.spec.copy()
                spec_tmp.update({'kp_grid_dens': 2})
                inpset = GWDFTDiagVaspInputSet(self.structure, spec_tmp, functional=self.spec['functional'])
                inpset.incar_settings.update({"ENCUT": 800})
            if self.spec['test'] or self.spec['converge']:
                inpset.set_test(self.option['test_prep'], self.option['value_prep'])
            inpset.get_incar(self.structure).write_file(os.path.join(path, 'INCAR.DIAG'))

        if self.job == 'G0W0':

            inpset = GWG0W0VaspInputSet(self.structure, self.spec, functional=self.spec['functional'])
            if self.spec['converge'] and not self.converged:
                spec_tmp = self.spec.copy()
                spec_tmp.update({'kp_grid_dens': 2})
                inpset = GWG0W0VaspInputSet(self.structure, spec_tmp, functional=self.spec['functional'])
                inpset.incar_settings.update({"ENCUT": 800})
            if self.spec['test'] or self.spec['converge']:
                inpset.set_test(self.option['test_prep'], self.option['value_prep'])
                inpset.set_test(self.option['test'], self.option['value'])
            if self.spec["prec"] == "h":
                inpset.set_prec_high()
            if self.spec['kp_grid_dens'] > 20:
                #inpset.wannier_on()
                inpset.write_input(self.structure, os.path.join(path, 'G0W0'+option_name))
                #w_inpset = Wannier90InputSet(self.spec)
                #w_inpset.write_file(self.structure, os.path.join(path, 'G0W0'+option_name))
            else:
                inpset.write_input(self.structure, os.path.join(path, 'G0W0'+option_name))

        if self.job == 'GW0':

            inpset = GWG0W0VaspInputSet(self.structure, self.spec, functional=self.spec['functional'])
            if self.spec['converge'] and not self.converged:
                spec_tmp = self.spec.copy()
                spec_tmp.update({'kp_grid_dens': 2})
                inpset = GWG0W0VaspInputSet(self.structure, spec_tmp, functional=self.spec['functional'])
                inpset.incar_settings.update({"ENCUT": 800})
            if self.spec['test'] or self.spec['converge']:
                inpset.set_test(self.option['test_prep'], self.option['value_prep'])
                inpset.set_test(self.option['test'], self.option['value'])
            if self.spec["prec"] == "h":
                inpset.set_prec_high()
            inpset.gw0_on()
            if self.spec['kp_grid_dens'] > 20:
                #inpset.wannier_on()
                inpset.write_input(self.structure, os.path.join(path, 'GW0'+option_name))
                #w_inpset = Wannier90InputSet(self.spec)
                #w_inpset.write_file(self.structure, os.path.join(path, 'GW0'+option_name))
            else:
                inpset.write_input(self.structure, os.path.join(path, 'GW0'+option_name))

        if self.job == 'scGW0':

            inpset = GWG0W0VaspInputSet(self.structure, self.spec, functional=self.spec['functional'])
            if self.spec['converge'] and not self.converged:
                spec_tmp = self.spec.copy()
                spec_tmp.update({'kp_grid_dens': 2})
                inpset = GWG0W0VaspInputSet(self.structure, spec_tmp, functional=self.spec['functional'])
                inpset.incar_settings.update({"ENCUT": 800})
            if self.spec['test'] or self.spec['converge']:
                inpset.set_test(self.option['test_prep'], self.option['value_prep'])
                inpset.set_test(self.option['test'], self.option['value'])
            if self.spec["prec"] == "h":
                inpset.set_prec_high()
            inpset.gw0_on(qpsc=True)
            if self.spec['kp_grid_dens'] > 20:
                inpset.wannier_on()
                inpset.write_input(self.structure, os.path.join(path, 'scGW0'+option_name))
                w_inpset = Wannier90InputSet(self.spec)
                w_inpset.write_file(self.structure, os.path.join(path, 'scGW0'+option_name))
            else:
                inpset.write_input(self.structure, os.path.join(path, 'scGW0'+option_name))

    def create_job_script(self, add_to_collection=True, mode='pbspro'):
        if mode == 'slurm':
            """
            Create job script for ceci.
            """
            npar = GWscDFTPrepVaspInputSet(self.structure, self.spec,
                                             functional=self.spec['functional']).get_npar(self.structure)
            if self.option is not None:
                option_prep_name = str(self.option['test_prep']) + str(self.option['value_prep'])
                if 'test' in self.option.keys():
                    option_name = str('.') + str(self.option['test']) + str(self.option['value'])
            else:
                option_prep_name = option_name = ''
            # npar = int(os.environ['NPARGWCALC'])
            header = ("#!/bin/bash \n"
                      "## standard header for Ceci clusters ## \n"
                      "#SBATCH --mail-user=michiel.vansetten@uclouvain.be \n"
                      "#SBATCH --mail-type=ALL\n"
                      "#SBATCH --time=2-24:0:0 \n"
                      "#SBATCH --cpus-per-task=1 \n"
                      "#SBATCH --mem-per-cpu=4000 \n")
            path_add = ''
            if self.spec['converge'] and self.converged:
                path_add = '.conv'
            if self.job == 'prep':
                path = os.path.join(s_name(self.structure) + path_add, option_prep_name)
                # create this job
                job_file = open(name=os.path.join(path, 'job'), mode='w')
                job_file.write(header)
                job_file.write('#SBATCH --job-name='+s_name(self.structure)+self.job+'\n')
                job_file.write('#SBATCH --ntasks='+str(npar)+'\n')
                job_file.write('module load vasp \n')
                job_file.write('mpirun vasp \n')
                job_file.write('cp OUTCAR OUTCAR.sc \n')
                job_file.write('cp INCAR.DIAG INCAR \n')
                job_file.write('mpirun vasp \n')
                job_file.write('cp OUTCAR OUTCAR.diag \n')
                job_file.close()
                os.chmod(os.path.join(path, 'job'), stat.S_IRWXU)
                if add_to_collection:
                    job_file = open("job_collection", mode='a')
                    job_file.write('cd ' + path + ' \n')
                    job_file.write('sbatch job \n')
                    job_file.write('cd .. \n')
                    job_file.close()
                    os.chmod("job_collection", stat.S_IRWXU)
            if self.job in ['G0W0', 'GW0', 'scGW0']:
                path = os.path.join(s_name(self.structure) + path_add, option_prep_name, self.job + option_name)
                # create this job
                job_file = open(name=path+'/job', mode='w')
                job_file.write(header)
                job_file.write('#SBATCH --job-name='+s_name(self.structure)+self.job+'\n')
                job_file.write('#SBATCH --ntasks='+str(npar)+'\n')
                job_file.write('module load vasp/5.2_par_wannier90 \n')
                job_file.write('cp ../CHGCAR ../WAVECAR ../WAVEDER . \n')
                job_file.write('mpirun vasp \n')
                job_file.write('rm W* \n')
                #job_file.write('workon pymatgen-GW; get_gap > gap; deactivate')
                #job_file.write('echo '+path+'`get_gap` >> ../../gaps.dat')
                job_file.close()
                os.chmod(path+'/job', stat.S_IRWXU)
                path = os.path.join(s_name(self.structure) + path_add, option_prep_name)
                # 'append submission of this job script to that of prep for this structure'
                if add_to_collection:
                    job_file = open(name=os.path.join(path, 'job'), mode='a')
                    job_file.write('cd ' + self.job + option_name + ' \n')
                    job_file.write('sbatch job \n')
                    job_file.write('cd .. \n')
                    job_file.close()
        elif mode == 'pbspro':
            """
            Create job script for pbse pro Zenobe.
            """
            npar = GWscDFTPrepVaspInputSet(self.structure, self.spec,
                                           functional=self.spec['functional']).get_npar(self.structure)
            #npar = 96
            if self.option is not None:
                option_prep_name = str(self.option['test_prep']) + str(self.option['value_prep'])
                if 'test' in self.option.keys():
                    option_name = str('.') + str(self.option['test']) + str(self.option['value'])
            else:
                option_prep_name = option_name = ''
            # npar = int(os.environ['NPARGWCALC'])
            header = str("#!/bin/bash \n" +
                         "## standard header for zenobe ## \n" +
                         "#!/bin/bash \n" +
                         "#PBS -q main\n" +
                         "#PBS -l walltime=24:0:00\n" +
                         "#PBS -r y \n" +
                         "#PBS -m abe\n" +
                         "#PBS -M michiel.vansetten@uclouvain.be\n" +
                         "#PBS -W group_list=naps\n" +
                         "#PBS -l pvmem=1900mb\n")
            path_add = ''
            if self.spec['converge'] and self.converged:
                path_add = '.conv'
            if self.job == 'prep':
                path = os.path.join(s_name(self.structure) + path_add, option_prep_name)
                abs_path = os.path.abspath(path)
                # create this job
                job_file = open(name=os.path.join(path, 'job'), mode='w')
                job_file.write(header)
                job_file.write("#PBS -l select=%s:ncpus=1:vmem=1900mb:mpiprocs=1:ompthreads=1\n" % str(npar))
                job_file.write('#PBS -o %s/queue.qout\n#PBS -e %s/queue.qerr\ncd %s\n' % (abs_path, abs_path, abs_path))
                job_file.write('mpirun -n %s vasp \n' % str(npar))
                job_file.write('cp OUTCAR OUTCAR.sc \n')
                job_file.write('cp INCAR.DIAG INCAR \n')
                job_file.write('mpirun -n %s vasp \n' % str(npar))
                job_file.write('cp OUTCAR OUTCAR.diag \n')
                job_file.close()
                os.chmod(os.path.join(path, 'job'), stat.S_IRWXU)
                if add_to_collection:
                    job_file = open("job_collection", mode='a')
                    job_file.write('cd ' + path + ' \n')
                    job_file.write('qsub job \n')
                    job_file.write('cd ../.. \n')
                    job_file.close()
                    os.chmod("job_collection", stat.S_IRWXU)
            if self.job in ['G0W0', 'GW0', 'scGW0']:
                path = os.path.join(s_name(self.structure) + path_add, option_prep_name, self.job + option_name)
                abs_path = os.path.abspath(path)
                # create this job
                job_file = open(name=path+'/job', mode='w')
                job_file.write(header)
                job_file.write("#PBS -l select=%s:ncpus=1:vmem=1000mb:mpiprocs=1:ompthreads=1\n" % str(npar))
                job_file.write('#PBS -o %s/queue.qout\n#PBS -e %s/queue.qerr\ncd %s\n' % (abs_path, abs_path, abs_path))
                job_file.write('cp ../CHGCAR ../WAVECAR ../WAVEDER . \n')
                job_file.write('mpirun -n %s vasp \n' % str(npar))
                job_file.write('rm W* \n')
                #job_file.write('workon pymatgen-GW; get_gap > gap; deactivate')
                #job_file.write('echo '+path+'`get_gap` >> ../../gaps.dat')
                job_file.close()
                os.chmod(path+'/job', stat.S_IRWXU)
                path = os.path.join(s_name(self.structure) + path_add, option_prep_name)
                # 'append submission of this job script to that of prep for this structure'
                if add_to_collection:
                    job_file = open(name=os.path.join(path, 'job'), mode='a')
                    job_file.write('cd ' + self.job + option_name + ' \n')
                    job_file.write('qsub job \n')
                    job_file.write('cd .. \n')
                    job_file.close()
