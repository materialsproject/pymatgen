import json
import glob
import itertools
import logging
import math
import os
import re
import warnings
from pathlib import Path
import xml.etree.cElementTree as ET
from collections import defaultdict
from io import StringIO

import numpy as np
from monty.io import zopen, reverse_readfile
from monty.json import MSONable
from monty.json import jsanitize
from monty.re import regrep
from monty.os.path import zpath

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure, Molecule
from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.electronic_structure.dos import Dos
from pymatgen.io.xyz import XYZ
"""
Classes for reading/manipulating/writing CP2K ouput files. 
Adapted from pymatgen.io.vasp.outputs 
"""

__author__ = "Nicholas Winner"
__version__ = "0.1"
__status__ = "Development"

logger = logging.getLogger(__name__)

_hartree_to_ev_ = 2.72113838565563E+01
_bohr_to_angstrom_ = 5.29177208590000E-01

# TODO This might need more testing
def _postprocessor(s):
    """
    Helper function to post process the results of the pattern matching functions in Cp2kOuput and turn them to
    python types.
    """
    s = s.rstrip()  # Remove leading/trailing whitespace
    s = s.lower()  # Turn to lower case for convenience

    if s == 'NO' or s == 'NONE':
        return False
    elif s == 'YES':
        return True
    elif re.search(s, '-?\d+\.\d+') or re.search(s, '-?.\d+'):
        try:
            return float(s)
        except ValueError:
            raise IOError("Error in parsing CP2K output file.")
    elif re.search(s, '-?\d+'):
        try:
            return int(s)
        except ValueError:
            raise IOError("Error in parsing CP2K output file.")
    else:
        return s


# TODO: Partition so that the constructor calls all default parsing methods -NW
class Cp2kOuput:

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

    def __init__(self, filename, verbose=False, auto_load=True):
        self.filename = filename
        self.is_stopped = False
        self.data = {}

        if auto_load:
            self._parse_global_params()  # Always present, parse the global parameters, most important is what run type
            self._parse_dft_params()  # TODO: assumes DFT, should be more general to support QMMM, MC, etc.
            self._ran_successfully()  # Only if job completed. No info about convergence etc.
            self._convergence()

            toten_pattern = re.compile(r"Total FORCE_EVAL.*\s(-?\d+.\d+)")
            self.read_pattern({'total_energy': toten_pattern},
                              terminate_on_match=True, postprocess=float, reverse=True)
            self.final_energy = self.data.get('total_energy', [])[-1][-1]

            self._parse_cell_params()  # Get the initial lattice
            self._parse_initial_structure()  # Get the initial structure by parsing lattice and then parsing coords
            self._parse_timing()  # Get timing info (includes total CPU time consumed, but also much more)

            # TODO: Is this the best way to implement? Should there just be the option to select each individually?
            if verbose:
                self._parse_scf_opt()
                self._parse_opt_steps()
                self._parse_total_numbers()

            self._parse_mulliken()
            self._parse_hirshfeld()

    @property
    def completed(self):
        return self.data.get('completed', False)

    @property
    def num_warnings(self):
        return self.data.get('num_warnings', 0)

    @property
    def run_type(self):
        return self.data.get('global').get('Run type')

    @property
    def spin_polarized(self):
        if ('UKS' or 'UNRESTRICTED_KOHN_SHAM' or 'LSD' or 'SPIN_POLARIZED') in \
                self.data['dft'].values():
            return True
        return False

    def _ran_successfully(self):
        """
        Sanity checks that the program ran successfully. Looks at the bottom of the CP2K ouput file
        for the "PROGRAM ENDED" line, which is printed when successfully ran. Also grabs the number
        of warnings issued.
        """
        program_ended_at = re.compile(r"PROGRAM ENDED AT\s+(\w+)")
        num_warnings = re.compile(r"The number of warnings for this run is : (\d+)")
        self.read_pattern(patterns={'completed': program_ended_at},
                          reverse=True, terminate_on_match=True, postprocess=bool)
        self.read_pattern(patterns={'num_warnings': num_warnings},
                          reverse=True, terminate_on_match=True, postprocess=int)

    def _convergence(self):
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

    def _parse_global_params(self):
        """
        Parse the GLOBAL section parameters from CP2K output file into a dictionary.
        """
        pat = re.compile(r"\s+GLOBAL\|\s+([\w+\s]*)\s+(\w+)")
        self.read_pattern({'global': pat}, terminate_on_match=False, postprocess=_postprocessor, reverse=False)
        self.data['global'] = dict(self.data['global'])

    def _parse_dft_params(self):
        """
        Parse the GLOBAL section parameters from CP2K output file into a dictionary.
        """
        pat = re.compile(r"\s+DFT\|\s+(\w.*)\s\s\s(.*)$")
        self.read_pattern({'dft': pat}, terminate_on_match=False, postprocess=_postprocessor, reverse=False)
        self.data['dft'] = dict(self.data['dft'])

        self.data['dft']['cutoffs'] = {}
        self.data['dft']['cutoffs']['density'] = self.data['dft']['cutoffs: density']
        self.data['dft']['cutoffs']['gradient'] = self.data['dft']['gradient']
        self.data['dft']['cutoffs']['tau'] = self.data['dft']['tau']

        self.data['dft'].pop('cutoffs: density')
        self.data['dft'].pop('gradient')
        self.data['dft'].pop('tau')

    def _parse_cell_params(self):
        cell_volume = re.compile(r"\s+CELL\|\sVolume.*\s(\d+\.\d+)")
        vectors = re.compile(r"\s+CELL\|.*(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)")
        angles = re.compile(r"\s+CELL\| Angle.*\s(\d+\.\d+)")
        self.read_pattern({'cell_volume': cell_volume, 'lattice': vectors, 'angles': angles},
                          terminate_on_match=False, postprocess=float, reverse=False)
        self.lattice = Lattice(self.data['lattice'])

    # TODO: CP2K Seems to only output initial struc here in output file. If so this can turn into list of structures
    def _parse_initial_structure(self):
        header = r"Atom\s+Kind\s+Element\s+X\s+Y\s+Z\s+Z\(eff\)\s+Mass"
        row = r"\s+(\d+)\s+(\d+)\s+(\w+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)"
        footer = r"^$"

        coord_table = self.read_table_pattern(header_pattern=header, row_pattern=row, footer_pattern=footer)

        self.structure = Structure(self.lattice,
                                   species=[i[2] for i in coord_table],
                                   coords=[[float(i[4]), float(i[5]), float(i[6])] for i in coord_table])

    def _parse_total_numbers(self):
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

    def _parse_scf_opt(self):
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

    def _parse_timing(self):
        header = r"SUBROUTINE\s+CALLS\s+ASD\s+SELF TIME\s+TOTAL TIME" + \
                 r"\s+MAXIMUM\s+AVERAGE\s+MAXIMUM\s+AVERAGE\s+MAXIMUM"
        row = r"(\w+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)"
        footer = r"\-+"

        timing = self.read_table_pattern(header_pattern=header, row_pattern=row, footer_pattern=footer,
                                         last_one_only=True)
        self.timing = {}
        for t in timing:
            self.timing[t[0]] = {'CALLS': {'MAX': int(t[1])},
                                'ASD': float(t[2]),
                                'SELF_TIME': {'AVERAGE':  float(t[3]), 'MAXIMUM':  float(t[4])},
                                'TOTAL_TIME': {'AVERAGE':  float(t[5]), 'MAXIMUM':  float(t[6])}}

    def _parse_opt_steps(self):
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

    def _parse_mulliken(self):
        header = r"Mulliken Population Analysis.+Net charge"
        pattern = r"\s+(\d)\s+(\w+)\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
        footer = r".+Total charge"

        d = self.read_table_pattern(header_pattern=header, row_pattern=pattern,
                                    footer_pattern=footer, last_one_only=False)

    def _parse_hirshfeld(self):
        header = r"Hirshfeld Charges.+Net charge"
        pattern = r"\s+(\d)\s+(\w+)\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
        footer = r"^$"

        d = self.read_table_pattern(header_pattern=header, row_pattern=pattern,
                                    footer_pattern=footer, last_one_only=False)

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


# TODO Use pymatgen's new "trajectory" object instead of a list of structures
def parse_structures(trajectory_file, lattice_file):
    mols = XYZ.from_file(trajectory_file).all_molecules
    lattice = np.loadtxt(lattice_file)
    structures = []
    for m, l in zip(mols, lattice):
        structures.append(Structure(lattice=l[2:].reshape(3,3),
                                    coords=[s.coords for s in m.sites],
                                    species=[s.specie for s in m.sites],
                                    coords_are_cartesian=True))


def parse_energies(energy_file):
    pass


def parse_pdos(pdos_file):

    with zopen(pdos_file) as f:
        lines = f.readlines()
        efermi = float(lines[0].split()[-2])*_hartree_to_ev_
        header = re.split(r'\s{2,}', lines[1].replace('#','').strip())

        dat = np.loadtxt('test_files/pdos')
        dat[:,1] = dat[:,1]*_hartree_to_ev_

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



