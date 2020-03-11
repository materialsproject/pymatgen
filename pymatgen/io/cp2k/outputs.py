import json
import glob
import itertools
import logging
import math
import os
import re
import numpy as np
import pandas as pd
import warnings

from monty.io import zopen, reverse_readfile
from monty.json import MSONable
from monty.json import jsanitize
from monty.re import regrep
from monty.os.path import zpath

from pymatgen.core.sites import Site
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure, Molecule
from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.electronic_structure.dos import Dos, CompleteDos, add_densities
from pymatgen.io.xyz import XYZ

from pymatgen.io.cp2k.inputs import Cp2kInput

"""
Classes for reading/manipulating/writing CP2K output files. 
"""

__author__ = "Nicholas Winner"
__version__ = "0.1"
__status__ = "Development"

logger = logging.getLogger(__name__)

_hartree_to_ev_ = 2.72113838565563E+01
_static_run_names_ = ['ENERGY', 'ENERGY_FORCE', 'WAVEFUNCTION_OPTIMIZATION', 'WFN_OPT']
_bohr_to_angstrom_ = 5.29177208590000E-01


# TODO This might need more testing
def _postprocessor(s):
    """
    Helper function to post process the results of the pattern matching functions in Cp2kOutput and turn them to
    python types.
    """
    s = s.rstrip()  # Remove leading/trailing whitespace
    s = s.lower()   # Turn to lower case for convenience
    s = s.replace(' ', '_')  # Remove whitespaces

    if s == 'no' or s == 'none':
        return False
    elif s == 'yes':
        return True
    elif re.match(r'^-?\d+\.\d+', s) or re.match(r'^-?\.\d+', s):
        try:
            return float(s)
        except ValueError:
            raise IOError("Error in parsing CP2K output file.")
    elif re.match(r'^-?\d+', s):
        try:
            return int(s)
        except ValueError:
            raise IOError("Error in parsing CP2K output file.")
    elif re.match(r'\*+', s):
        try:
            return np.NaN
        except ValueError:
            raise IOError("Error in parsing CP2K output file.")
    else:
        return s


# TODO: Partition so that the constructor calls all default parsing methods -NW
class Cp2kOutput:

    """
    Class for parsing output file from CP2K. The CP2K Input file is very flexible in the way that it is returned.
    This class will automatically parse parameters that should always be present, but other parsing features may be
    called depending on the run type.

    Current RUN_TYPES supported:
        ENERGY
        ENERGY_FORCE
        GEO_OPT
        MD

    Args:
        filename: (str) Name of the CP2K output file to parse.
        verbose: (bool) Whether or not to parse with verbosity (will parse lots of data that may not be useful)

    """

    def __init__(self, filename, verbose=False, auto_load=False):
        # IO Info
        self.filename = filename
        self.dir = os.path.dirname(filename)

        self.data = {}

        # Material properties/results
        self.input = None
        self.initial_structure = None
        self.lattice = None
        self.final_structure = None
        self.composition = None
        self.efermi = None
        self.vbm = None
        self.cbm = None
        self.band_gap = None
        self.structures = []
        self.ionic_steps = []

        # parse the basic run parameters always
        self.parse_cp2k_params()
        self.parse_input()  # parse the input file
        self.parse_global_params()  # Always present, parse the global parameters, most important is what run type
        self.parse_dft_params()  # Present so long as a DFT calculation was performed
        self.parse_scf_params()
        self.parse_atomic_kind_info()

        # Auto-load will load the most crucial data into the data attribute
        if auto_load:
            self.ran_successfully()  # Only if job completed. No info about convergence etc.
            self.convergence()  # Checks to see if job converged

            self.parse_energies() # get total energy for each ionic step
            self.parse_forces() # get forces on all atoms (in order), if available
            self.parse_stresses() # get stress tensor and total stress at each ionic step, if available
            self.parse_ionic_steps() # collect energy, forces, and total stress into ionic steps variable

            self.parse_initial_structure()  # Get the initial structure by parsing lattice and then parsing coords
            self.parse_structures() # collect all structurs from the run

            self.parse_mo_eigenvalues()  # Get the eigenvalues of the MOs (for finding gaps, VBM, CBM)
            self.parse_timing()  # Get timing info (includes total CPU time consumed, but also much more)

            # TODO: Is this the best way to implement? Should there just be the option to select each individually?
            if verbose:
                self.parse_scf_opt()
                self.parse_opt_steps()
                self.parse_total_numbers()
                self.parse_mulliken()
                self.parse_hirshfeld()

    @property
    def cp2k_version(self):
        return self.data.get('cp2k_version', None)

    @property
    def completed(self):
        try:
            return self.data.get('completed', False)[0][0]
        except:
            return False

    @property
    def num_warnings(self):
        return self.data.get('num_warnings', 0)

    @property
    def run_type(self):
        return self.data.get('global').get('run_type')

    @property
    def project_name(self):
        return self.data.get('global').get("project_name")

    @property
    def spin_polarized(self):
        if ('UKS' or 'UNRESTRICTED_KOHN_SHAM' or 'LSD' or 'SPIN_POLARIZED') in \
                self.data['dft'].values():
            return True
        return False

    @property
    def is_metal(self):
        if self.band_gap <= 0:
            return True

    # TODO Maybe I should create a parse_files function that globs to get the file names instead of putting
        # it in each function seperate? -NW
    def parse_structures(self, trajectory_file=None, lattice_file=None):
        """
        Parses the structures from a cp2k calculation. Static calculations simply use the initial structure.
        For calculations with ionic motion, the function will look for the appropriate trajectory and lattice
        files based on naming convention. If no file is given, and no file is found, it is assumed
        that the lattice/structure remained constant, and the initial lattice/structure is used.
        Cp2k does not output the trajectory in the main output file by default, so non static calculations have to
        reference the trajectory file.
        """
        if lattice_file is None:
            lattice = "{}-1.cell".format(self.project_name)
            lattice = glob.glob(os.path.join(self.dir, lattice+'*'))
            if len(lattice) == 0:
                lattice = self.parse_initial_structure().lattice
            elif len(lattice) == 1:
                latfile = np.loadtxt(lattice[0])
                lattice = [l[2:11].reshape(3,3) for l in latfile] \
                    if len(latfile.shape)>1 else latfile[2:11].reshape(3, 3)
            else:
                raise FileNotFoundError("Unable to automatically determine lattice file. More than one exist.")
        else:
            latfile = np.loadtxt(lattice_file)
            lattice = [l[2:].reshape(3, 3) for l in latfile]

        if trajectory_file is None:
            trajectory_file = "{}-pos-1.xyz".format(self.project_name)
            trajectory_file = glob.glob(os.path.join(self.dir, trajectory_file+'*'))
            if len(trajectory_file) == 0:
                self.structures = []
                self.structures.append(self.parse_initial_structure())
                self.final_structure = self.structures[-1]
            elif len(trajectory_file) == 1:
                mols = XYZ.from_file(trajectory_file[0]).all_molecules
                self.structures = []
                for m, l in zip(mols, lattice):
                    self.structures.append(Structure(lattice=l, coords=[s.coords for s in m.sites],
                                                     species=[s.specie for s in m.sites], coords_are_cartesian=True))
                self.final_structure = self.structures[-1]
            else:
                raise FileNotFoundError("Unable to automatically determine trajectory file. More than one exist.")
        else:
            mols = XYZ.from_file(trajectory_file).all_molecules
            self.structures = []
            for m, l in zip(mols, lattice):
                self.structures.append(Structure(lattice=l, coords=[s.coords for s in m.sites],
                                                 species=[s.specie for s in m.sites], coords_are_cartesian=True))
            self.final_structure = self.structures[-1]

    # TODO: CP2K Seems to only output initial struc here in output file. If so this can turn into list of structures
    def parse_initial_structure(self):
        header = r"Atom\s+Kind\s+Element\s+X\s+Y\s+Z\s+Z\(eff\)\s+Mass"
        row = r"\s+(\d+)\s+(\d+)\s+(\w+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)"
        footer = r"^$"

        self.parse_cell_params()
        coord_table = self.read_table_pattern(header_pattern=header, row_pattern=row, footer_pattern=footer)
        self.initial_structure = Structure(self.lattice,
                                   species=[i[2] for i in coord_table],
                                   coords=[[float(i[4]), float(i[5]), float(i[6])] for i in coord_table])

        self.composition = self.initial_structure.composition
        return self.initial_structure

    def ran_successfully(self):
        """
        Sanity checks that the program ran successfully. Looks at the bottom of the CP2K output file
        for the "PROGRAM ENDED" line, which is printed when successfully ran. Also grabs the number
        of warnings issued.
        """
        program_ended_at = re.compile(r"PROGRAM ENDED AT\s+(\w+)")
        num_warnings = re.compile(r"The number of warnings for this run is : (\d+)")
        self.read_pattern(patterns={'completed': program_ended_at},
                          reverse=True, terminate_on_match=True, postprocess=bool)
        self.read_pattern(patterns={'num_warnings': num_warnings},
                          reverse=True, terminate_on_match=True, postprocess=int)

        if not self.completed:
            raise ValueError("The provided CP2K job did not finish running! Cannot parse the file reliably.")

    def convergence(self):
        # SCF Loops
        uncoverged_inner_loop = re.compile(r"(Leaving inner SCF loop)")
        scf_converged = re.compile(r"(SCF run converged)")
        scf_not_converged = re.compile(r"(SCF run NOT converged)")
        self.read_pattern(patterns={'uncoverged_inner_loop': uncoverged_inner_loop,
                                    'scf_converged': scf_converged,
                                    'scf_not_converged': scf_not_converged},
                          reverse=True,
                          terminate_on_match=False, postprocess=bool)

        # GEO_OPT
        geo_opt_not_converged = re.compile(r"(MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED)")
        geo_opt_converged = re.compile(r"(GEOMETRY OPTIMIZATION COMPLETED)")
        self.read_pattern(patterns={'geo_opt_converged': geo_opt_converged,
                                    'geo_opt_not_converged': geo_opt_not_converged},
                          reverse=True, terminate_on_match=True, postprocess=bool)

        if any(self.data['scf_not_converged']):
            warnings.warn("There is at least one unconverged SCF cycle in the provided cp2k calculation",
                          UserWarning)
        if any(self.data['geo_opt_not_converged']):
            warnings.warn("Geometry optimization did not converge",
                          UserWarning)

    def parse_energies(self):
        toten_pattern = re.compile(r"Total FORCE_EVAL.*\s(-?\d+.\d+)")
        self.read_pattern({'total_energy': toten_pattern},
                          terminate_on_match=True, postprocess=float, reverse=False)
        self.final_energy = self.data.get('total_energy', [])[-1][-1]

    def parse_forces(self):
        header_pattern = r"ATOMIC FORCES.+Z"
        row_pattern = r"\s+\d+\s+\d+\s+\w+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
        footer_pattern = r"SUM OF ATOMIC FORCES"

        self.data['forces'] = self.read_table_pattern(header_pattern=header_pattern, row_pattern=row_pattern,
                                                      footer_pattern=footer_pattern, postprocess=_postprocessor,
                                                      last_one_only=False)

    def parse_stresses(self):
        header_pattern = r"STRESS TENSOR.+Z"
        row_pattern = r"\s+\w+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
        footer_pattern = r"^$"

        self.data['stress_tensor'] = self.read_table_pattern(header_pattern=header_pattern, row_pattern=row_pattern,
                                                             footer_pattern=footer_pattern, postprocess=_postprocessor,
                                                             last_one_only=False)

        trace_pattern = re.compile(r"Trace\(stress tensor.+(-?\d+\.\d+E?-?\d+)")
        self.read_pattern({'stress': trace_pattern}, terminate_on_match=False, postprocess=float, reverse=False)

    def parse_ionic_steps(self):
        self.ionic_steps = []
        for i in range(len(self.data['total_energy'])):
            self.ionic_steps.append({})
            try:
                self.ionic_steps[i]['E'] = self.data['total_energy'][i][0]
            except:
                raise warnings.warn('No total energies idenfied! Check output file')
            try:
                self.ionic_steps[i]['forces'] = self.data['forces'][i]
            except:
                pass
            try:
                self.ionic_steps[i]['stress'] = self.data['stress'][i][0]
            except:
                pass

    def parse_cp2k_params(self):
        """
        Parse the CP2K general parameters from CP2K output file into a dictionary.
        """
        version = re.compile(r"\s+CP2K\|.+(\d\.\d)")
        input_file = re.compile(r"\s+CP2K\|\s+Input file name\s+(.+)$")
        self.read_pattern({'cp2k_version': version, 'input_filename': input_file},
                           terminate_on_match=True, reverse=False)

    def parse_input(self):
        input_filename = self.data['input_filename'][0][0]
        for ext in ['', '.gz', '.GZ', '.z', '.Z', '.bz2', '.BZ2']:
            if os.path.exists(os.path.join(self.dir, input_filename+ext)):
                self.input = Cp2kInput.from_file(os.path.join(self.dir, input_filename+ext))
                return
        raise warnings.warn("Original input file not found. Some info may be lost.")

    def parse_global_params(self):
        """
        Parse the GLOBAL section parameters from CP2K output file into a dictionary.
        """
        pat = re.compile(r"\s+GLOBAL\|\s+([\w+\s]*)\s+(\w+)")
        self.read_pattern({'global': pat}, terminate_on_match=False, reverse=False)
        for d in self.data['global']:
            d[0], d[1] = _postprocessor(d[0]), str(d[1])
        self.data['global'] = dict(self.data['global'])

    def parse_dft_params(self):
        """
        Parse the DFT parameters (as well as functional, HF, vdW params)
        """
        pat = re.compile(r"\s+DFT\|\s+(\w.*)\s\s\s(.*)$")
        self.read_pattern({'dft': pat}, terminate_on_match=False, postprocess=_postprocessor, reverse=False)
        self.data['dft'] = dict(self.data['dft'])

        self.data['dft']['cutoffs'] = {}
        self.data['dft']['cutoffs']['density'] = self.data['dft'].pop('cutoffs:_density')
        self.data['dft']['cutoffs']['gradient'] = self.data['dft'].pop('gradient')
        self.data['dft']['cutoffs']['tau'] = self.data['dft'].pop('tau')

        # Functional
        functional = re.compile(r"\s+FUNCTIONAL\|\s+(.+):")
        self.read_pattern({'functional': functional}, terminate_on_match=False, postprocess=_postprocessor, reverse=False)
        self.data['dft']['functional'] = [item for sublist in self.data.pop('functional') for item in sublist]

        # HF exchange info
        hfx = re.compile(r"\s+HFX_INFO\|\s+(.+):\s+(.*)$")
        self.read_pattern({'hfx': hfx}, terminate_on_match=False, postprocess=_postprocessor, reverse=False)
        if len(self.data['hfx']) > 0:
            self.data['dft']['hfx'] = dict(self.data.pop('hfx'))

        # Van der waals correction
        vdw = re.compile(r"\s+vdW POTENTIAL\|\s+(DFT-D.)\s")
        self.read_pattern({'vdw': vdw}, terminate_on_match=False, postprocess=_postprocessor, reverse=False)
        if len(self.data['vdw']) > 0:
            self.data['dft']['vdw'] = self.data.pop('vdw')[0][0]

    def parse_scf_params(self):
        """
        Retrieve the most import SCF parameters: the max number of scf cycles (max_scf),
        the convergence cutoff for scf (eps_scf),
        :return:
        """
        max_scf = re.compile(r"max_scf:\s+(\d+)")
        eps_scf = re.compile(r"eps_scf:\s+(\d+)")
        self.read_pattern({'max_scf': max_scf, 'eps_scf': eps_scf},
                          terminate_on_match=True, reverse=False)
        self.data['scf'] = {}
        self.data['scf']['max_scf'] = self.data.pop('max_scf')[0][0]
        self.data['scf']['eps_scf'] = self.data.pop('eps_scf')[0][0]

    def parse_cell_params(self):
        cell_volume = re.compile(r"\s+CELL\|\sVolume.*\s(\d+\.\d+)")
        vectors = re.compile(r"\s+CELL\| Vector.*\s(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)")
        angles = re.compile(r"\s+CELL\| Angle.*\s(\d+\.\d+)")
        self.read_pattern({'cell_volume': cell_volume, 'lattice': vectors, 'angles': angles},
                          terminate_on_match=False, postprocess=float, reverse=False)
        self.lattice = Lattice(self.data['lattice'])
        return self.lattice

    def parse_atomic_kind_info(self):
        kinds = re.compile(r"Atomic kind: (\w+)")
        orbital_basis_set = re.compile(r"Orbital Basis Set\s+(.+$)")
        potential_information = re.compile(r"Potential information for\s+(.+$)")
        auxiliary_basis_set = re.compile(r"Auxiliary Fit Basis Set\s+(.+$)")
        core_electrons = re.compile(r'Total number of core electrons\s+(\d+)')
        valence_electrons = re.compile(r'Total number of valence electrons\s+(\d+)')
        self.read_pattern({'kinds': kinds, 'orbital_basis_set': orbital_basis_set, 'potential_info':
                           potential_information, 'auxiliary_basis_set': auxiliary_basis_set,
                           "core_electrons": core_electrons, "valence_electrons": valence_electrons},
                          terminate_on_match=True, postprocess=str, reverse=False)
        atomic_kind_info = {}
        for i, kind in enumerate(self.data['kinds']):
            atomic_kind_info[kind[0]] = {'orbital_basis_set': self.data.get('orbital_basis_set')[i][0],
                                         'pseudo_potential': self.data.get('potential_info')[i][0]}
            try:
                atomic_kind_info[kind[0]]['valence_electrons'] = self.data.get('valence_electrons')[i][0]
            except:
                atomic_kind_info[kind[0]]['valence_electrons'] = None
            try:
                atomic_kind_info[kind[0]]['core_electrons'] = self.data.get('core_electrons')[i][0]
            except:
                atomic_kind_info[kind[0]]['core_electrons'] = None
            try:
                atomic_kind_info[kind[0]]['auxiliary_basis_set'] = self.data.get('auxiliary_basis_set')[i]
            except:
                atomic_kind_info[kind[0]]['auxiliary_basis_set'] = None

        self.data['atomic_kind_info'] = atomic_kind_info

    def parse_total_numbers(self):
        atomic_kinds = r'- Atomic kinds:\s+(\d+)'
        atoms = r'- Atoms:\s+(\d+)'
        shell_sets = r'- Shell sets:\s+(\d+)'
        shells = r'- Shells:\s+(\d+)'
        primitive_funcs = r'- Primitive Cartesian functions:\s+(\d+)'
        cart_base_funcs = r'- Cartesian basis functions:\s+(\d+)'
        spher_base_funcs = r'- Spherical basis functions:\s+(\d+)'

        self.read_pattern({'atomic_kinds': atomic_kinds, 'atoms': atoms, 'shell_sets': shell_sets,
                           'shells': shells, 'primitive_cartesian_functions': primitive_funcs,
                           'cartesian_basis_functions': cart_base_funcs, 'spherical_basis_functions':
                           spher_base_funcs},
                          terminate_on_match=True)

    def parse_scf_opt(self):
        header = r"Step\s+Update method\s+Time\s+Convergence\s+Total energy\s+Change" + \
            r"\s+\-+"
        row = r"(\d+)\s+(\w+\s?\w+)\s+(\d+\.\d+E\+\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+E[\+\-]?\d+)"
        footer = r"^$"

        scfs = self.read_table_pattern(header_pattern=header,
                                       row_pattern=row, footer_pattern=footer,
                                       last_one_only=False)

        self.data['electronic_steps'] = scfs

        #self.data['electronic_steps'] = []
        #for i in scfs:
        #    self.data['electronic_steps'].append([float(j[-2]) for j in i])

    def parse_timing(self):
        header = r"SUBROUTINE\s+CALLS\s+ASD\s+SELF TIME\s+TOTAL TIME" + \
                 r"\s+MAXIMUM\s+AVERAGE\s+MAXIMUM\s+AVERAGE\s+MAXIMUM"
        row = r"(\w+)\s+(.+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)"
        footer = r"\-+"

        timing = self.read_table_pattern(header_pattern=header, row_pattern=row, footer_pattern=footer,
                                         last_one_only=True, postprocess=_postprocessor)
        self.timing = {}
        for t in timing:
            self.timing[t[0]] = {'calls': {'max': t[1]},
                                'asd':t[2],
                                'self_time': {'average':  t[3], 'maximum':  t[4]},
                                'total_time': {'average':  t[5], 'maximum':  t[6]}}

    def parse_opt_steps(self):
        # "Informations at step =" Summary block (floating point terms)
        total_energy = re.compile(r"\s+Total Energy\s+=\s+(-?\d+.\d+)")
        real_energy_change = re.compile(r"\s+Real energy change\s+=\s+(-?\d+.\d+)")
        prediced_change_in_energy = re.compile(r"\s+Predicted change in energy\s+=\s+(-?\d+.\d+)")
        scaling_factor = re.compile(r"\s+Scaling factor\s+=\s+(-?\d+.\d+)")
        step_size = re.compile(r"\s+Step size\s+=\s+(-?\d+.\d+)")
        trust_radius = re.compile(r"\s+Trust radius\s+=\s+(-?\d+.\d+)")
        used_time = re.compile(r"\s+Used time\s+=\s+(-?\d+.\d+)")

        # For RUN_TYPE=CELL_OPT
        pressure_deviation = re.compile(r"\s+Pressure Deviation.*=\s+(-?\d+.\d+)")
        pressure_tolerance = re.compile(r"\s+Pressure Tolerance.*=\s+(-?\d+.\d+)")

        self.read_pattern({'total_energy': total_energy, 'real_energy_change': real_energy_change,
                           'predicted_change_in_energy': prediced_change_in_energy, 'scaling_factor': scaling_factor,
                           'step_size': step_size, 'trust_radius': trust_radius, 'used_time': used_time,
                           'pressure_deviation': pressure_deviation, 'pressure_tolerance': pressure_tolerance},
                          terminate_on_match=False, postprocess=float)

        # "Informations at step =" Summary block (bool terms)
        decrease_in_energy = re.compile(r"\s+Decrease in energy\s+=\s+(\w+)")
        converged_step_size = re.compile(r"\s+Convergence in step size\s+=\s+(\w+)")
        converged_rms_step = re.compile(r"\s+Convergence in RMS step\s+=\s+(\w+)")
        converged_in_grad = re.compile(r"\s+Conv\. in gradients\s+=\s+(\w+)")
        converged_in_rms_grad = re.compile(r"\s+Conv\. in RMS gradients\s+=\s+(\w+)")
        pressure_converged = re.compile(r"\s+Conv\. for  PRESSURE\s+=\s+(\w+)")

        self.read_pattern({'decrease_in_energy': decrease_in_energy, 'converged_step_size': converged_step_size,
                           'converged_rms_step': converged_rms_step, 'converged_in_grad': converged_in_grad,
                           'converged_in_rms_grad': converged_in_rms_grad, 'pressure_converged': pressure_converged},
                          terminate_on_match=False, postprocess=_postprocessor)

    def parse_mulliken(self):
        header = r"Mulliken Population Analysis.+Net charge"
        pattern = r"\s+(\d)\s+(\w+)\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
        footer = r".+Total charge"

        d = self.read_table_pattern(header_pattern=header, row_pattern=pattern,
                                    footer_pattern=footer, last_one_only=False)

    def parse_hirshfeld(self):
        uks = self.spin_polarized
        header = r"Hirshfeld Charges.+Net charge"
        footer = r"^$"

        if not uks:
            pattern = r"\s+(\d)\s+(\w+)\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
            d = self.read_table_pattern(header_pattern=header, row_pattern=pattern,
                                        footer_pattern=footer, last_one_only=False)
            for i, ionic_step in enumerate(d):
                population = []
                net_charge = []
                for site in ionic_step:
                    population.append(site[4])
                    net_charge.append(site[5])
                hirshfeld = [{'population': population[j],
                              'net_charge': net_charge[j]}
                             for j in range(len(population))]
                self.structures[i].add_site_property('hirshfield', hirshfeld)
        else:
            pattern = r"\s+(\d)\s+(\w+)\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
            d = self.read_table_pattern(header_pattern=header, row_pattern=pattern,
                                        footer_pattern=footer, last_one_only=False)
            for i, ionic_step in enumerate(d):
                population = []
                net_charge = []
                spin_moment = []
                for site in ionic_step:
                    population.append(tuple(site[4:5]))
                    spin_moment.append(site[6])
                    net_charge.append(site[7])
                hirshfeld = [{'population': population[j],
                              'net_charge': net_charge[j],
                              'spin_moment': spin_moment[j]}
                             for j in range(len(population))]
                self.structures[i].add_site_property('hirshfield', hirshfeld)

    def parse_mo_eigenvalues(self):
        """
        Parse the MO eigenvalues from the cp2k output file. Will get the eigenvalues (and band gap)
        at each ionic step (if more than one exist).

        Everything is decomposed by spin channel. If calculation was performed without spin polarization,
        then only Spin.up will be present, which represents the average of up and down.
        """
        eigenvalues = []
        band_gap = []
        efermi = []

        with zopen(self.filename, 'rt') as f:
            lines = iter(f.readlines())
            for line in lines:
                if line.__contains__(' occupied subspace spin'):
                    eigenvalues.append({'occupied': {Spin.up: [], Spin.down: []},
                                        'unoccupied': {Spin.up: [], Spin.down: []}})
                    band_gap.append({Spin.up: None, Spin.down: None})
                    efermi.append({Spin.up: None, Spin.down: None})

                    next(lines)
                    while True:
                        line = next(lines)
                        if line.__contains__('Fermi'):
                            efermi[-1][Spin.up] = float(line.split()[-1])
                            break
                        eigenvalues[-1]['occupied'][Spin.up].extend([_hartree_to_ev_*float(l) for l in line.split()])
                    next(lines)
                    line = next(lines)
                    if line.__contains__(' occupied subspace spin'):
                        next(lines)
                        while True:
                            line = next(lines)
                            if line.__contains__('Fermi'):
                                efermi[-1][Spin.down] = float(line.split()[-1])
                                break
                            eigenvalues[-1]['occupied'][Spin.down].extend([_hartree_to_ev_*float(l)
                                                                           for l in line.split()])
                if line.__contains__(' unoccupied subspace spin'):
                    next(lines)
                    next(lines)
                    while True:
                        line = next(lines)
                        if line.__contains__('Eigenvalues'):
                            break
                        elif line.__contains__('HOMO'):
                            band_gap[-1][Spin.up] = float(line.split()[-1])
                            break
                        eigenvalues[-1]['unoccupied'][Spin.up].extend([_hartree_to_ev_*float(l)
                                                                       for l in line.split()])
                    next(lines)
                    next(lines)
                    if line.__contains__(' unoccupied subspace spin'):
                        while True:
                            line = next(lines)
                            if line.__contains__('HOMO'):
                                band_gap[-1][Spin.up] = float(line.split()[-1])
                                line = next(lines)
                                band_gap[-1][Spin.down] = float(line.split()[-1])
                                break
                            eigenvalues[-1]['unoccupied'][Spin.down].extend([_hartree_to_ev_*float(l)
                                                                             for l in line.split()])

        self.data['eigenvalues'] = eigenvalues
        self.data['band_gap'] = band_gap

        # self.data will always contained the eigenvalues resolved by spin channel. The average vbm, cbm, gap,
        # and fermi are saved as class attributes, as there is (usually) no assymmetry in these values for
        # common materials
        if self.spin_polarized:
            self.data['vbm'] = {Spin.up: np.max(eigenvalues[-1]['occupied'][Spin.up]),
                                Spin.down: np.max(eigenvalues[-1]['occupied'][Spin.down])}
            self.data['cbm'] = {Spin.up: np.min(eigenvalues[-1]['unoccupied'][Spin.up]),
                                Spin.down: np.min(eigenvalues[-1]['unoccupied'][Spin.down])}
            self.vbm = (self.data['vbm'][Spin.up] + self.data['vbm'][Spin.down])/2
            self.cbm = (self.data['cbm'][Spin.up] + self.data['cbm'][Spin.down])/2
            self.efermi = (efermi[-1][Spin.up]+efermi[-1][Spin.down])/2
            self.band_gap = (band_gap[-1][Spin.up]+band_gap[-1][Spin.down])/2
        else:
            self.data['vbm'] = {Spin.up: np.max(eigenvalues[-1]['occupied'][Spin.up]),
                                Spin.down: None}
            self.data['cbm'] = {Spin.up: np.min(eigenvalues[-1]['unoccupied'][Spin.up]),
                                Spin.down: None}
            self.vbm = self.data['vbm'][Spin.up]
            self.cbm = self.data['cbm'][Spin.up]
            self.efermi = efermi[-1][Spin.up]
            self.band_gap = band_gap[-1][Spin.up]

    def parse_homo_lumo(self):
        """
        Find the HOMO - LUMO gap in [eV]. Returns the last value. For gaps/eigenvalues decomposed by
        spin up/spin down channel and over many ionic steps, see parse_mo_eigenvalues()
        """
        pattern = re.compile(r"HOMO - LUMO gap.*\s(-?\d+.\d+)")
        self.read_pattern(patterns={'band_gap': pattern}, reverse=True, terminate_on_match=True, postprocess=float)

    def parse_pdos(self, pdos_files=None):
        """
        Parse the pdos files created by cp2k, and assimilate them into a CompleteDos object.
        Either provide a list of PDOS file paths, or use glob to find the .pdos extension in
        the calculation directory.

        Args:
            pdos_files (list): list of pdos file paths

        Returns:
            CompleteDos
        """
        if pdos_files is None:
            pdos_files = glob.glob(os.path.join(self.dir, '*.pdos'))
            if len(pdos_files) == 0:
                raise FileNotFoundError("Unable to automatically located PDOS files in the calculation directory")
        pdoss = {}
        for pdos_file in pdos_files:
            if os.path.split(pdos_file)[-1].__contains__('ALPHA'):
                spin = 1
            else:
                spin = -1
            with zopen(pdos_file, 'rt') as f:
                lines = f.readlines()
                kind = re.search(r"atomic kind\s(.*)\sat iter", lines[0]).groups()[0]
                if kind not in pdoss.keys():
                    pdoss[kind] = {}
                efermi = float(lines[0].split()[-2]) * _hartree_to_ev_
                header = re.split(r'\s{2,}', lines[1].replace('#', '').strip())
                dat = np.loadtxt(pdos_file)
                dat[:, 1] = dat[:, 1] * _hartree_to_ev_
                for i in range(len(header)):
                    if header[i] == 'd-2':
                        header[i] = 'dxy'
                    elif header[i] == 'd-1':
                        header[i] = 'dyz'
                    elif header[i] == 'd0':
                        header[i] = 'dz2'
                    elif header[i] == 'd+1':
                        header[i] = 'dxz'
                    elif header[i] == 'd+2':
                        header[i] = 'dx2'
                    elif header[i] == 'f-3':
                        header[i] = 'f_3'
                    elif header[i] == 'f-2':
                        header[i] = 'f_2'
                    elif header[i] == 'f-1':
                        header[i] = 'f_1'
                    elif header[i] == 'f0':
                        header[i] = 'f0'
                    elif header[i] == 'f+1':
                        header[i] = 'f1'
                    elif header[i] == 'f+2':
                        header[i] = 'f2'
                    elif header[i] == 'f+3':
                        header[i] = 'f3'
                densities = []
                energies = []
                for i in range(2, len(header)):
                    if pdoss[kind].get(getattr(Orbital, header[i])):
                        continue
                    else:
                        pdoss[kind][getattr(Orbital, header[i])] = {Spin.up: [], Spin.down: []}
                for r in dat:
                    energies.append(float(r[1]))
                    for i in range(3, len(r)):
                        pdoss[kind][getattr(Orbital, header[i - 1])][Spin(spin)].append(r[i])

        print(pdoss)
        tdos_densities = {Spin.up: np.zeros(len(energies)), Spin.down: np.zeros(len(energies))}
        for el, orbitals in pdoss.items():
            for orbital, d in orbitals.items():
                tdos_densities = add_densities(tdos_densities, d)

        print(tdos_densities)

        return Dos(efermi=efermi, energies=energies, densities=tdos_densities)

        # sort and assimilate the total dos
        sort = np.argsort(tdos_energies)
        tdos_energies = list(np.array(tdos_energies)[sort])
        tdos_densities[Spin.up] = list(np.array(tdos_densities[Spin.up])[sort])
        tdos_densities[Spin.down] = list(np.array(tdos_densities[Spin.down])[sort])

        # TODO confim this works
        combined_energies = [tdos_energies[0]]
        combined_tdos = {Spin.up: [tdos_densities[Spin.up][0]],
                         Spin.down: [tdos_densities[Spin.down][0]]}
        for i in range(len(tdos_energies)-1):
            if (tdos_energies[i] == tdos_energies[i+1]) and (tdos_densities[Spin.up][i] != 0.0):
                combined_tdos[Spin.up][-1] += tdos_densities[Spin.up][i+1]
                combined_tdos[Spin.down][-1] += tdos_densities[Spin.down][i+1]
            else:
                combined_energies.append(tdos_energies[i+1])
                combined_tdos[Spin.up].append(tdos_densities[Spin.up][i+1])
                combined_tdos[Spin.down].append(tdos_densities[Spin.down][i+1])

        #print(tdos_energies)
        #print(tdos_densities[Spin.down])
        #print()
        #print(combined_energies)
        #print(combined_tdos[Spin.down])
        # TODO confirm combined dos works
        #tdos = Dos(efermi, tdos_densities, combined_tdos)
        tdos = Dos(efermi, tdos_energies, tdos_densities)
        return CompleteDos(self.initial_structure, total_dos=tdos, pdoss=pdoss)

    def read_pattern(self, patterns, reverse=False, terminate_on_match=False,
                     postprocess=str):
        """
        This function originally comes from pymatgen.io.vasp.outputs Outcar class

        General pattern reading. Uses monty's regrep method. Takes the same
        arguments.

        Args:
            patterns (dict): A dict of patterns, e.g.,
                {"energy": r"energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)"}.
            reverse (bool): Read files in reverse. Defaults to false. Useful for
                large files, esp OUTCARs, especially when used with
                terminate_on_match.
            terminate_on_match (bool): Whether to terminate when there is at
                least one match in each key in pattern.
            postprocess (callable): A post processing function to convert all
                matches. Defaults to str, i.e., no change.

        Renders accessible:
            Any attribute in patterns. For example,
            {"energy": r"energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)"} will set the
            value of self.data["energy"] = [[-1234], [-3453], ...], to the
            results from regex and postprocess. Note that the returned values
            are lists of lists, because you can grep multiple items on one line.
        """
        matches = regrep(self.filename, patterns, reverse=reverse,
                         terminate_on_match=terminate_on_match,
                         postprocess=postprocess)
        for k in patterns.keys():
            self.data[k] = [i[0] for i in matches.get(k, [])]

    def read_table_pattern(self, header_pattern, row_pattern, footer_pattern,
                           postprocess=str, attribute_name=None,
                           last_one_only=True):
        """
        This function originally comes from pymatgen.io.vasp.outputs Outcar class

        Parse table-like data. A table composes of three parts: header,
        main body, footer. All the data matches "row pattern" in the main body
        will be returned.

        Args:
            header_pattern (str): The regular expression pattern matches the
                table header. This pattern should match all the text
                immediately before the main body of the table. For multiple
                sections table match the text until the section of
                interest. MULTILINE and DOTALL options are enforced, as a
                result, the "." meta-character will also match "\n" in this
                section.
            row_pattern (str): The regular expression matches a single line in
                the table. Capture interested field using regular expression
                groups.
            footer_pattern (str): The regular expression matches the end of the
                table. E.g. a long dash line.
            postprocess (callable): A post processing function to convert all
                matches. Defaults to str, i.e., no change.
            attribute_name (str): Name of this table. If present the parsed data
                will be attached to "data. e.g. self.data["efg"] = [...]
            last_one_only (bool): All the tables will be parsed, if this option
                is set to True, only the last table will be returned. The
                enclosing list will be removed. i.e. Only a single table will
                be returned. Default to be True.

        Returns:
            List of tables. 1) A table is a list of rows. 2) A row if either a list of
            attribute values in case the the capturing group is defined without name in
            row_pattern, or a dict in case that named capturing groups are defined by
            row_pattern.
        """
        with zopen(self.filename, 'rt') as f:
            text = f.read()
        table_pattern_text = header_pattern + r"\s*^(?P<table_body>(?:\s+" + \
            row_pattern + r")+)\s+" + footer_pattern
        table_pattern = re.compile(table_pattern_text, re.MULTILINE | re.DOTALL)
        rp = re.compile(row_pattern)
        tables = []
        for mt in table_pattern.finditer(text):
            table_body_text = mt.group("table_body")
            table_contents = []
            for line in table_body_text.split("\n"):
                ml = rp.search(line)
                d = ml.groupdict()
                if len(d) > 0:
                    processed_line = {k: postprocess(v) for k, v in d.items()}
                else:
                    processed_line = [postprocess(v) for v in ml.groups()]
                table_contents.append(processed_line)
            tables.append(table_contents)
        if last_one_only:
            retained_data = tables[-1]
        else:
            retained_data = tables
        if attribute_name is not None:
            self.data[attribute_name] = retained_data
        return retained_data

    # TODO: for modularity, maybe this should return a raw dict instead of data dict
    # and the data dict can be created by the drone
    def as_dict(self):
        d = {'input': {}, 'output': {}}
        d['total_time'] = self.timing['cp2k']['total_time']['maximum']
        d['run_type'] = self.run_type
        d['input']['global'] = self.data.get('global')
        d['input']['dft'] = self.data.get('dft', None)
        d['input']['scf'] = self.data.get('scf', None)
        d['input']['structure'] = self.initial_structure.as_dict()
        d['input']['atomic_kind_info'] = self.data.get('atomic_kind_info', None)
        d['ran_successfully'] = self.completed
        d['cp2k_version'] = self.cp2k_version
        d['output']['structure'] = self.final_structure.as_dict()
        d['output']['ionic_steps'] = self.ionic_steps
        d['composition'] = self.composition.as_dict()
        d['output']['energy'] = self.final_energy
        d['output']['energy_per_atom'] = self.final_energy / self.composition.num_atoms
        d['output']['bandgap'] = self.band_gap
        d['output']['cbm'] = self.cbm
        d['output']['vbm'] = self.vbm
        d['output']['efermi'] = self.efermi
        d['output']['is_metal'] = self.is_metal
        return d


# TODO Use pymatgen's new "trajectory" object instead of a list of structures
def parse_structures(trajectory_file, lattice, final=False, step_skip=1):
    mols = XYZ.from_file(trajectory_file, step_skip=step_skip).all_molecules
    if isinstance(lattice, Lattice):
        lattice = [lattice]
    elif os.path.isfile(lattice):
        latfile = np.loadtxt(lattice)
        lattice = [l[2:].reshape(3,3) for l in latfile]
    else:
        raise ValueError("Lattice format not regonized! We take Lattice objects",
                         "or paths to a cp2k cell file.")
    structures = []
    for m, l in zip(mols, lattice):
        structures.append(Structure(lattice=l, coords=[s.coords for s in m.sites],
                                    species=[s.specie for s in m.sites],
                                    coords_are_cartesian=True))
    if final:
        return structures[-1]
    return structures


def parse_energy_file(energy_file):
    """
    Parses energy file for calculations with multiple ionic steps.
    """
    columns = [
        'step', 'kinetic_energy', 'temp', 'potential_energy', 'conserved_quantity', 'used_time'
    ]
    df = pd.read_table(energy_file, skiprows=1, names=columns, sep='\s+')
    df['kinetic_energy'] = df['kinetic_energy'] * _hartree_to_ev_
    df['potential_energy'] = df['potential_energy'] * _hartree_to_ev_
    df['conserved_quantity'] = df['conserved_quantity'] * _hartree_to_ev_
    df.astype(float)
    d = {c: df[c].values for c in columns}
    return d


def parse_pdos(pdos_file):
    with zopen(pdos_file, 'rt') as f:
        lines = f.readlines()
        efermi = float(lines[0].split()[-2])*_hartree_to_ev_
        header = re.split(r'\s{2,}', lines[1].replace('#', '').strip())
        dat = np.loadtxt(pdos_file)
        dat[:, 1] = dat[:, 1]*_hartree_to_ev_
        for i in range(len(header)):
            if header[i] == 'd-2':
                header[i] = 'dxy'
            elif header[i] == 'd-1':
                header[i] = 'dyz'
            elif header[i] == 'd0':
                header[i] = 'dz2'
            elif header[i] == 'd+1':
                header[i] = 'dxz'
            elif header[i] == 'd+2':
                header[i] = 'dx2'
            elif header[i] == 'f-3':
                header[i] = 'f_3'
            elif header[i] == 'f-2':
                header[i] = 'f_2'
            elif header[i] == 'f-1':
                header[i] = 'f_1'
            elif header[i] == 'f0':
                header[i] = 'f0'
            elif header[i] == 'f+1':
                header[i] = 'f1'
            elif header[i] == 'f+2':
                header[i] = 'f2'
            elif header[i] == 'f+3':
                header[i] = 'f3'
        densities = {}
        energies = []
        for i in range(2, len(header)):
            densities[getattr(Orbital, header[i])] = {Spin.up: []}
        for r in dat:
            energies.append(float(r[1]))
            for i in range(3, len(r)):
                densities[getattr(Orbital, header[i-1])][Spin.up].append(r[i])
        dos = {}
        for k,v in densities.items():
            dos[k] = Dos(efermi=efermi, energies=energies, densities=v)
        return dos


class Cube:
        """
        From ERG Research Group with minor modifications.
        """
        def __init__(self, fname):
            f = zopen(fname, 'rt')

            # skip header lines
            for i in range(2): f.readline()

            # number of atoms included in the file followed by the position of the origin of the volumetric data
            line = f.readline().split()
            self.natoms = int(line[0])
            self.origin = np.array(np.array(list(map(float, line[1:]))))

            # The next three lines give the number of voxels along each axis (x, y, z) followed by the axis vector.
            line = f.readline().split()
            self.NX = int(line[0])
            self.X = np.array([0.529177*float(l) for l in line[1:]])

            line = f.readline().split()
            self.NY = int(line[0])
            self.Y = np.array([0.529177*float(l) for l in line[1:]])

            line = f.readline().split()
            self.NZ = int(line[0])
            self.Z = np.array([0.529177*float(l) for l in line[1:]])

            self.volume = abs(np.dot(np.cross(self.X, self.Y), self.Z))

            # The last section in the header is one line for each atom consisting of 5 numbers,
            # the first is the atom number, second (?), the last three are the x,y,z coordinates of the atom center.
            self.sites = []
            for i in range(self.natoms):
                line = f.readline().split()
                self.sites.append(Site(line[0], list(map(float, line[2:]))))

            # Volumetric data
            self.data = np.zeros((self.NX, self.NY, self.NZ))
            i = 0
            for s in f:
                for v in s.split():
                    self.data[int(i / (self.NY * self.NZ)), int((i / self.NZ) % self.NY), int(i % self.NZ)] = float(v)
                    i += 1

        def mask_sphere(self, r, coord):
            # produce spheric volume mask with radius R and center @ [Cx,Cy,Cz]
            # can be used for integration over spherical part of the volume
            m = 0 * self.data
            Cx, Cy, Cz = coord
            for ix in range(int(math.ceil((Cx - r) / self.X[0])), int(math.floor((Cx + r) / self.X[0]))):
                ryz = math.sqrt(r ** 2 - (ix * self.X[0] - Cx) ** 2)
                for iy in range(int(math.ceil((Cy - ryz) / self.Y[1])), int(math.floor((Cy + ryz) / self.Y[1]))):
                    rz = math.sqrt(ryz ** 2 - (iy * self.Y[1] - Cy) ** 2)
                    for iz in range(int(math.ceil((Cz - rz) / self.Z[2])), int(math.floor((Cz + rz) / self.Z[2]))):
                        m[ix, iy, iz] = 1
            return m

        def dump(self, f):
            # output Gaussian cube into file descriptor "f".
            # Usage pattern: f=open('filename.cube'); cube.dump(f); f.close()
            print >> f, "CUBE file\ngenerated by piton _at_ erg.biophys.msu.ru"
            print >> f, "%4d %.6f %.6f %.6f" % (self.natoms, self.origin[0], self.origin[1], self.origin[2])
            print >> f, "%4d %.6f %.6f %.6f" % (self.NX, self.X[0], self.X[1], self.X[2])
            print >> f, "%4d %.6f %.6f %.6f" % (self.NY, self.Y[0], self.Y[1], self.Y[2])
            print >> f, "%4d %.6f %.6f %.6f" % (self.NZ, self.Z[0], self.Z[1], self.Z[2])
            for atom in self.atoms:
                print >> f, "%s %d %s %s %s" % (atom[0], 0, atom[1], atom[2], atom[3])
            for ix in range(self.NX):
                for iy in range(self.NY):
                    for iz in range(self.NZ):
                        print >> f, "%.5e " % self.data[ix, iy, iz],
                        if (iz % 6 == 5): print >> f, ''
                    print >> f, ""

        def integrate(self):
            return ((self.data)**2).sum()*self.volume

        def integrate_sphere(self, r, coord):
            mask = self.mask_sphere(r, coord)
            return ((mask*self.data)**2).sum()*self.volume