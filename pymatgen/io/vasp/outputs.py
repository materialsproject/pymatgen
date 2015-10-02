# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

"""
Classes for reading/manipulating/writing VASP ouput files.
"""


__author__ = "Shyue Ping Ong, Geoffroy Hautier, Rickard Armiento, " + \
    "Vincent L Chevrier, Ioannis Petousis, Stephen Dacek"
__credits__ = "Anubhav Jain"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.2"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Nov 30, 2012"

import os
import glob
import re
import math
import itertools
import warnings
from io import StringIO
import logging
from collections import defaultdict
from xml.etree.cElementTree import iterparse

from six.moves import map, zip
from six import string_types

import numpy as np

from monty.io import zopen, reverse_readfile
from monty.re import regrep
from monty.json import jsanitize

from pymatgen.util.io_utils import clean_lines, micro_pyawk
from pymatgen.core.structure import Structure
from pymatgen.core.units import unitized
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.electronic_structure.dos import CompleteDos, Dos
from pymatgen.electronic_structure.bandstructure import BandStructure, \
    BandStructureSymmLine, get_reconstructed_band_structure
from pymatgen.core.lattice import Lattice
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, Potcar
from pymatgen.entries.computed_entries import \
    ComputedEntry, ComputedStructureEntry
from pymatgen.serializers.json_coders import PMGSONable

logger = logging.getLogger(__name__)


def _parse_parameters(val_type, val):
    """
    Helper function to convert a Vasprun parameter into the proper type.
    Boolean, int and float types are converted.

    Args:
        val_type: Value type parsed from vasprun.xml.
        val: Actual string value parsed for vasprun.xml.
    """
    if val_type == "logical":
        return val == "T"
    elif val_type == "int":
        return int(val)
    elif val_type == "string":
        return val.strip()
    else:
        return float(val)


def _parse_v_parameters(val_type, val, filename, param_name):
    """
    Helper function to convert a Vasprun array-type parameter into the proper
    type. Boolean, int and float types are converted.

    Args:
        val_type: Value type parsed from vasprun.xml.
        val: Actual string value parsed for vasprun.xml.
        filename: Fullpath of vasprun.xml. Used for robust error handling.
            E.g., if vasprun.xml contains \*\*\* for some Incar parameters,
            the code will try to read from an INCAR file present in the same
            directory.
        param_name: Name of parameter.

    Returns:
        Parsed value.
    """
    if val_type == "logical":
        val = [i == "T" for i in val.split()]
    elif val_type == "int":
        try:
            val = [int(i) for i in val.split()]
        except ValueError:
            # Fix for stupid error in vasprun sometimes which displays
            # LDAUL/J as 2****
            val = _parse_from_incar(filename, param_name)
            if val is None:
                raise IOError("Error in parsing vasprun.xml")
    elif val_type == "string":
        val = val.split()
    else:
        try:
            val = [float(i) for i in val.split()]
        except ValueError:
            # Fix for stupid error in vasprun sometimes which displays
            # MAGMOM as 2****
            val = _parse_from_incar(filename, param_name)
            if val is None:
                raise IOError("Error in parsing vasprun.xml")
    return val


def _parse_varray(elem):
    return [[float(i) for i in v.text.split()] for v in elem]


def _parse_from_incar(filename, key):
    """
    Helper function to parse a parameter from the INCAR.
    """
    dirname = os.path.dirname(filename)
    for f in os.listdir(dirname):
        if re.search("INCAR", f):
            warnings.warn("INCAR found. Using " + key + " from INCAR.")
            incar = Incar.from_file(os.path.join(dirname, f))
            if key in incar:
                return incar[key]
            else:
                return None
    return None


def _vasprun_float(f):
    """
    Large numbers are often represented as ********* in the vasprun.
    This function parses these values as np.nan
    """
    try:
        return float(f)
    except ValueError as e:
        f = f.strip()
        if f == '*' * len(f):
            warnings.warn('Float overflow (*******) encountered in vasprun')
            return np.nan
        raise e


class Vasprun(PMGSONable):
    """
    Vastly improved cElementTree-based parser for vasprun.xml files. Uses
    iterparse to support incremental parsing of large files.
    Speedup over Dom is at least 2x for smallish files (~1Mb) to orders of
    magnitude for larger files (~10Mb).

    Args:
        filename (str): Filename to parse
        ionic_step_skip (int): If ionic_step_skip is a number > 1,
            only every ionic_step_skip ionic steps will be read for
            structure and energies. This is very useful if you are parsing
            very large vasprun.xml files and you are not interested in every
            single ionic step. Note that the final energies may not be the
            actual final energy in the vasprun.
        ionic_step_offset (int): Used together with ionic_step_skip. If set,
            the first ionic step read will be offset by the amount of
            ionic_step_offset. For example, if you want to start reading
            every 10th structure but only from the 3rd structure onwards,
            set ionic_step_skip to 10 and ionic_step_offset to 3. Main use
            case is when doing statistical structure analysis with
            extremely long time scale multiple VASP calculations of
            varying numbers of steps.
        parse_dos (bool): Whether to parse the dos. Defaults to True. Set
            to False to shave off significant time from the parsing if you
            are not interested in getting those data.
        parse_eigen (bool): Whether to parse the eigenvalues. Defaults to
            True. Set to False to shave off significant time from the
            parsing if you are not interested in getting those data.
        parse_projected_eigen (bool): Whether to parse the projected
            eigenvalues. Defaults to False. Set to True to obtain projected
            eigenvalues. **Note that this can take an extreme amount of time
            and memory.** So use this wisely.
        parse_potcar_file (bool/str): Whether to parse the potcar file to read the
            potcar hashes for the potcar_spec attribute. Defaults to True,
            where no hashes will be determined and the potcar_spec dictionaries
            will read {"symbol": ElSymbol, "hash": None}. By Default, looks in
            the same directory as the vasprun.xml, with same extensions as
             Vasprun.xml. If a string is provided, looks at that filepath

    **Vasp results**

    .. attribute:: ionic_steps

        All ionic steps in the run as a list of
        {"structure": structure at end of run,
        "electronic_steps": {All electronic step data in vasprun file},
        "stresses": stress matrix}

    .. attribute:: structures

        List of Structure objects for the structure at each ionic step.

    .. attribute:: tdos

        Total dos calculated at the end of run.

    .. attribute:: idos

        Integrated dos calculated at the end of run.

    .. attribute:: pdos

        List of list of PDos objects. Access as pdos[atomindex][orbitalindex]

    .. attribute:: efermi

        Fermi energy

    .. attribute:: eigenvalues

        Available only if parse_eigen=True. Final eigenvalues as a dict of
        {(spin, kpoint index):[[eigenvalue, occu]]}.
        This representation is based on actual ordering in VASP and is meant as
        an intermediate representation to be converted into proper objects. The
        kpoint index is 0-based (unlike the 1-based indexing in VASP).

    .. attribute:: projected_eigenvalues

        Final projected eigenvalues as a dict of
        {(atom index, band index, kpoint index, Orbital, Spin):float}
        This representation is based on actual ordering in VASP and is meant as
        an intermediate representation to be converted into proper objects. The
        kpoint, band and atom indices are 0-based (unlike the 1-based indexing
        in VASP).

    .. attribute:: dielectric

        The real and imaginary part of the dielectric constant (e.g., computed
        by RPA) in function of the energy (frequency). Optical properties (e.g.
        absorption coefficient) can be obtained through this.
        The data is given as a tuple of 3 values containing each of them
        the energy, the real part tensor, and the imaginary part tensor
        ([energies],[[real_partxx,real_partyy,real_partzz,real_partxy,
        real_partyz,real_partxz]],[[imag_partxx,imag_partyy,imag_partzz,
        imag_partxy, imag_partyz, imag_partxz]])

    .. attribute:: other_dielectric

        Dictionary, with the tag comment as key, containing other variants of
        the real and imaginary part of the dielectric constant (e.g., computed
        by RPA) in function of the energy (frequency). Optical properties (e.g.
        absorption coefficient) can be obtained through this.
        The data is given as a tuple of 3 values containing each of them
        the energy, the real part tensor, and the imaginary part tensor
        ([energies],[[real_partxx,real_partyy,real_partzz,real_partxy,
        real_partyz,real_partxz]],[[imag_partxx,imag_partyy,imag_partzz,
        imag_partxy, imag_partyz, imag_partxz]])

    .. attribute:: epsilon_static

        The static part of the dielectric constant. Present when it's a DFPT run
        (LEPSILON=TRUE)

    .. attribute:: epsilon_static_wolfe

        The static part of the dielectric constant without any local field effects.
        Present when it's a DFPT run (LEPSILON=TRUE)

    .. attribute:: epsilon_ionic

        The ionic part of the static dielectric constant. Present when it's a DFPT run
        (LEPSILON=TRUE) and IBRION=5, 6, 7 or 8

    .. attribute:: nionic_steps

        The total number of ionic steps. This number is always equal
        to the total number of steps in the actual run even if
        ionic_step_skip is used.

    **Vasp inputs**

    .. attribute:: incar

        Incar object for parameters specified in INCAR file.

    .. attribute:: parameters

        Incar object with parameters that vasp actually used, including all
        defaults.

    .. attribute:: kpoints

        Kpoints object for KPOINTS specified in run.

    .. attribute:: actual_kpoints

        List of actual kpoints, e.g.,
        [[0.25, 0.125, 0.08333333], [-0.25, 0.125, 0.08333333],
        [0.25, 0.375, 0.08333333], ....]

    .. attribute:: actual_kpoints_weights

        List of kpoint weights, E.g.,
        [0.04166667, 0.04166667, 0.04166667, 0.04166667, 0.04166667, ....]

    .. attribute:: atomic_symbols

        List of atomic symbols, e.g., ["Li", "Fe", "Fe", "P", "P", "P"]

    .. attribute:: potcar_symbols

        List of POTCAR symbols. e.g.,
        ["PAW_PBE Li 17Jan2003", "PAW_PBE Fe 06Sep2000", ..]

    Author: Shyue Ping Ong
    """

    def __init__(self, filename, ionic_step_skip=None,
                 ionic_step_offset=0, parse_dos=True,
                 parse_eigen=True, parse_projected_eigen=False,
                 parse_potcar_file=True):
        self.filename = filename
        self.ionic_step_skip = ionic_step_skip
        self.ionic_step_offset = ionic_step_offset

        with zopen(filename, "rt") as f:
            if ionic_step_skip or ionic_step_offset:
                # remove parts of the xml file and parse the string
                run = f.read()
                steps = run.split("<calculation>")
                #The text before the first <calculation> is the preamble!
                preamble = steps.pop(0)
                self.nionic_steps = len(steps)
                new_steps = steps[ionic_step_offset::int(ionic_step_skip)]
                #add the tailing informat in the last step from the run
                to_parse = "<calculation>".join(new_steps)
                if steps[-1] != new_steps[-1]:
                    to_parse = "{}<calculation>{}{}".format(
                        preamble, to_parse,
                        steps[-1].split("</calculation>")[-1])
                else:
                    to_parse = "{}<calculation>{}".format(preamble, to_parse)
                self._parse(StringIO(to_parse), parse_dos=parse_dos,
                            parse_eigen=parse_eigen,
                            parse_projected_eigen=parse_projected_eigen)
            else:
                self._parse(f, parse_dos=parse_dos, parse_eigen=parse_eigen,
                            parse_projected_eigen=parse_projected_eigen)
                self.nionic_steps = len(self.ionic_steps)

            if parse_potcar_file:
                self.update_potcar_spec(parse_potcar_file)

    def _parse(self, stream, parse_dos, parse_eigen, parse_projected_eigen):
        self.efermi = None
        self.eigenvalues = None
        self.projected_eigenvalues = None
        self.other_dielectric = {}
        ionic_steps = []
        parsed_header = False
        for event, elem in iterparse(stream):
            tag = elem.tag
            if not parsed_header:
                if tag == "generator":
                    self.generator = self._parse_params(elem)
                elif tag == "incar":
                    self.incar = self._parse_params(elem)
                elif tag == "kpoints":
                    self.kpoints, self.actual_kpoints, \
                        self.actual_kpoints_weights = self._parse_kpoints(elem)
                elif tag == "parameters":
                    self.parameters = self._parse_params(elem)
                elif tag == "structure" and elem.attrib.get("name") == \
                        "initialpos":
                    self.initial_structure = self._parse_structure(elem)
                elif tag == "atominfo":
                    self.atomic_symbols, self.potcar_symbols = \
                        self._parse_atominfo(elem)
                    self.potcar_spec = [{"titel": p,
                                         "hash": None} for
                                        p in self.potcar_symbols]
            if tag == "calculation":
                parsed_header = True
                ionic_steps.append(self._parse_calculation(elem))
            if tag == "dielectricfunction":
                if ("comment" not in elem.attrib) or \
                   elem.attrib["comment"] == "INVERSE MACROSCOPIC DIELECTRIC TENSOR (including local field effects in RPA (Hartree))":
                    self.dielectric = self._parse_diel(elem)
                else:
                    self.other_dielectric[elem.attrib["comment"]] = self._parse_diel(elem)
            elif parse_dos and tag == "dos":
                try:
                    self.tdos, self.idos, self.pdos = self._parse_dos(elem)
                    self.efermi = self.tdos.efermi
                    self.dos_has_errors = False
                except Exception as ex:
                    self.dos_has_errors = True
            elif parse_eigen and tag == "eigenvalues":
                self.eigenvalues = self._parse_eigen(elem)
            elif parse_projected_eigen and tag == "projected":
                self.projected_eigenvalues = self._parse_projected_eigen(elem)
            elif tag == "structure" and elem.attrib.get("name") == \
                    "finalpos":
                self.final_structure = self._parse_structure(elem)
        self.ionic_steps = ionic_steps
        self.vasp_version = self.generator["version"]

    @property
    def structures(self):
        return [step["structure"] for step in self.ionic_steps]

    @property
    def epsilon_static(self):
        """
        Property only available for DFPT calculations.
        """
        return self.ionic_steps[-1].get("epsilon", [])

    @property
    def epsilon_static_wolfe(self):
        """
        Property only available for DFPT calculations.
        """
        return self.ionic_steps[-1].get("epsilon_rpa", [])

    @property
    def epsilon_ionic(self):
        """
        Property only available for DFPT calculations and when IBRION=5, 6, 7 or 8.
        """
        return self.ionic_steps[-1].get("epsilon_ion", [])

    @property
    def lattice(self):
        return self.final_structure.lattice

    @property
    def lattice_rec(self):
        return self.final_structure.lattice.reciprocal_lattice

    @property
    def converged_electronic(self):
        """
        Checks that electronic step convergence has been reached in the final
        ionic step
        """
        final_esteps = self.ionic_steps[-1]["electronic_steps"]
        if 'LEPSILON' in self.incar and self.incar['LEPSILON']:
            i = 1
            to_check = set(['e_wo_entrp', 'e_fr_energy', 'e_0_energy'])
            while set(final_esteps[i].keys()) == to_check:
                i += 1
            return i + 1 != self.parameters["NELM"]
        return len(final_esteps) < self.parameters["NELM"]

    @property
    def converged_ionic(self):
        """
        Checks that ionic step convergence has been reached, i.e. that vasp
        exited before reaching the max ionic steps for a relaxation run
        """
        nsw = self.parameters.get("NSW", 0)
        return nsw <= 1 or len(self.ionic_steps) < nsw

    @property
    def converged(self):
        """
        Returns true if a relaxation run is converged.
        """
        return self.converged_electronic and self.converged_ionic

    @property
    @unitized("eV")
    def final_energy(self):
        """
        Final energy from the vasp run.
        """
        try:
            return self.ionic_steps[-1]["electronic_steps"][-1]["e_0_energy"]
        except (IndexError, KeyError):
            # not all calculations have a total energy, i.e. GW
            return np.inf

    @property
    def complete_dos(self):
        """
        A complete dos object which incorporates the total dos and all
        projected dos.
        """
        final_struct = self.final_structure
        pdoss = {final_struct[i]: pdos for i, pdos in enumerate(self.pdos)}
        return CompleteDos(self.final_structure, self.tdos, pdoss)

    @property
    def hubbards(self):
        """
        Hubbard U values used if a vasprun is a GGA+U run. {} otherwise.
        """
        symbols = [s.split()[1] for s in self.potcar_symbols]
        symbols = [re.split("_", s)[0] for s in symbols]
        if not self.incar.get("LDAU", False):
            return {}
        us = self.incar.get("LDAUU", self.parameters.get("LDAUU"))
        js = self.incar.get("LDAUJ", self.parameters.get("LDAUJ"))
        if len(us) == len(symbols):
            return {symbols[i]: us[i] - js[i] for i in range(len(symbols))}
        elif sum(us) == 0 and sum(js) == 0:
            return {}
        else:
            raise VaspParserError("Length of U value parameters and atomic "
                                  "symbols are mismatched")

    @property
    def run_type(self):
        """
        Returns the run type. Currently supports only GGA and HF calcs.

        TODO: Fix for other functional types like LDA, PW91, etc.
        """
        if self.is_hubbard:
            return "GGA+U"
        elif self.parameters.get("LHFCALC", False):
            return "HF"
        else:
            return "GGA"

    @property
    def is_hubbard(self):
        """
        True if run is a DFT+U run.
        """
        if len(self.hubbards) == 0:
            return False
        return sum(self.hubbards.values()) > 1e-8

    @property
    def is_spin(self):
        """
        True if run is spin-polarized.
        """
        return self.parameters.get("ISPIN", 1) == 2

    def get_computed_entry(self, inc_structure=False, parameters=None,
                           data=None):
        """
        Returns a ComputedStructureEntry from the vasprun.

        Args:
            inc_structure (bool): Set to True if you want
                ComputedStructureEntries to be returned instead of
                ComputedEntries.
            parameters (list): Input parameters to include. It has to be one of
                the properties supported by the Vasprun object. If
                parameters == None, a default set of parameters that are
                necessary for typical post-processing will be set.
            data (list): Output data to include. Has to be one of the properties
                supported by the Vasprun object.

        Returns:
            ComputedStructureEntry/ComputedEntry
        """
        param_names = {"is_hubbard", "hubbards", "potcar_symbols",
                       "potcar_spec", "run_type"}
        if parameters:
            param_names.update(parameters)
        params = {p: getattr(self, p) for p in param_names}
        data = {p: getattr(self, p) for p in data} if data is not None else {}

        if inc_structure:
            return ComputedStructureEntry(self.final_structure,
                                          self.final_energy, parameters=params,
                                          data=data)
        else:
            return ComputedEntry(self.final_structure.composition,
                                 self.final_energy, parameters=params,
                                 data=data)

    def get_band_structure(self, kpoints_filename=None, efermi=None,
                           line_mode=False):
        """
        Returns the band structure as a BandStructure object

        Args:
            kpoints_filename (str): Full path of the KPOINTS file from which
                the band structure is generated.
                If none is provided, the code will try to intelligently
                determine the appropriate KPOINTS file by substituting the
                filename of the vasprun.xml with KPOINTS.
                The latter is the default behavior.
            efermi (float): If you want to specify manually the fermi energy
                this is where you should do it. By default, the None value
                means the code will get it from the vasprun.
            line_mode (bool): Force the band structure to be considered as
                a run along symmetry lines.

        Returns:
            a BandStructure object (or more specifically a
            BandStructureSymmLine object if the run is detected to be a run
            along symmetry lines)

            Two types of runs along symmetry lines are accepted: non-sc with
            Line-Mode in the KPOINT file or hybrid, self-consistent with a
            uniform grid+a few kpoints along symmetry lines (explicit KPOINTS
            file) (it's not possible to run a non-sc band structure with hybrid
            functionals). The explicit KPOINTS file needs to have data on the
            kpoint label as commentary.
        """

        if not kpoints_filename:
            kpoints_filename = self.filename.replace('vasprun.xml', 'KPOINTS')
        if not os.path.exists(kpoints_filename) and line_mode is True:
            raise VaspParserError('KPOINTS needed to obtain band structure '
                                  'along symmetry lines.')

        if efermi is None:
            efermi = self.efermi

        kpoint_file = None
        if os.path.exists(kpoints_filename):
            kpoint_file = Kpoints.from_file(kpoints_filename)
        lattice_new = Lattice(self.lattice_rec.matrix)

        kpoints = [np.array(self.actual_kpoints[i])
                   for i in range(len(self.actual_kpoints))]
        dict_eigen = self.as_dict()['output']['eigenvalues']
        dict_p_eigen = {}
        if 'projected_eigenvalues' in self.as_dict()['output']:
            dict_p_eigen = self.as_dict()['output']['projected_eigenvalues']

        p_eigenvals = {}
        if "1" in dict_eigen["0"] and "-1" in dict_eigen["0"] \
                and self.incar['ISPIN'] == 2:
            eigenvals = {Spin.up: [], Spin.down: []}
            if len(dict_p_eigen) != 0:
                p_eigenvals = {Spin.up: [], Spin.down: []}
        else:
            eigenvals = {Spin.up: []}
            if len(dict_p_eigen) != 0:
                p_eigenvals = {Spin.up: []}

        neigenvalues = [len(v['1']) for v in dict_eigen.values()]
        min_eigenvalues = min(neigenvalues)
        get_orb = Orbital.from_string
        for i in range(min_eigenvalues):
            eigenvals[Spin.up].append([dict_eigen[str(j)]['1'][i][0]
                                       for j in range(len(kpoints))])
            if len(dict_p_eigen) != 0:
                p_eigenvals[Spin.up].append(
                    [{get_orb(orb): dict_p_eigen[j]['1'][i][orb]
                      for orb in dict_p_eigen[j]['1'][i]}
                     for j in range(len(kpoints))])
        if Spin.down in eigenvals:
            for i in range(min_eigenvalues):
                eigenvals[Spin.down].append([dict_eigen[str(j)]['-1'][i][0]
                                             for j in range(len(kpoints))])
                if len(dict_p_eigen) != 0:
                    p_eigenvals[Spin.down].append(
                        [{get_orb(orb): dict_p_eigen[j]['-1'][i][orb]
                          for orb in dict_p_eigen[j]['-1'][i]}
                         for j in range(len(kpoints))]
                    )

        # check if we have an hybrid band structure computation
        #for this we look at the presence of the LHFCALC tag
        hybrid_band = False
        if self.parameters.get('LHFCALC', False):
            hybrid_band = True

        if kpoint_file is not None:
            if kpoint_file.style == "Line_mode":
                line_mode = True

        if line_mode:
            labels_dict = {}
            if hybrid_band:
                start_bs_index = 0
                for i in range(len(self.actual_kpoints)):
                    if self.actual_kpoints_weights[i] == 0.0:
                        start_bs_index = i
                        break
                for i in range(len(kpoint_file.kpts)):
                    if kpoint_file.labels[i] is not None:
                        labels_dict[kpoint_file.labels[i]] = \
                            kpoint_file.kpts[i]
                #remake the data only considering line band structure k-points
                #(weight = 0.0 kpoints)
                kpoints = kpoints[start_bs_index:len(kpoints)]
                up_eigen = [eigenvals[Spin.up][i][
                            start_bs_index:len(eigenvals[Spin.up][i])]
                            for i in range(len(eigenvals[Spin.up]))]
                if self.is_spin:
                    down_eigen = [eigenvals[Spin.down][i]
                                  [start_bs_index:
                                  len(eigenvals[Spin.down][i])]
                                  for i in range(len(eigenvals[Spin.down]))]
                    eigenvals = {Spin.up: up_eigen,
                                 Spin.down: down_eigen}
                else:
                    eigenvals = {Spin.up: up_eigen}
            else:
                if '' in kpoint_file.labels:
                    raise Exception("A band structure along symmetry lines "
                                    "requires a label for each kpoint. "
                                    "Check your KPOINTS file")
                labels_dict = dict(zip(kpoint_file.labels, kpoint_file.kpts))
                labels_dict.pop(None, None)
            return BandStructureSymmLine(kpoints, eigenvals, lattice_new,
                                         efermi, labels_dict,
                                         structure=self.final_structure,
                                         projections=p_eigenvals)
        else:
            return BandStructure(kpoints, eigenvals, lattice_new, efermi,
                                 structure=self.final_structure,
                                 projections=p_eigenvals)

    @property
    def eigenvalue_band_properties(self):
        """
        Band properties from the eigenvalues as a tuple,
        (band gap, cbm, vbm, is_band_gap_direct).
        """
        vbm = -float("inf")
        vbm_kpoint = None
        cbm = float("inf")
        cbm_kpoint = None
        for k, val in self.eigenvalues.items():
            for (eigenval, occu) in val:
                if occu > 1e-8 and eigenval > vbm:
                    vbm = eigenval
                    vbm_kpoint = k[0]
                elif occu <= 1e-8 and eigenval < cbm:
                    cbm = eigenval
                    cbm_kpoint = k[0]
        return max(cbm - vbm, 0), cbm, vbm, vbm_kpoint == cbm_kpoint

    def update_potcar_spec(self, path):
        def get_potcar_in_path(p):
            for fn in os.listdir(os.path.abspath(p)):
                if 'POTCAR' in fn:
                    pc = Potcar.from_file(os.path.join(p, fn))
                    if {d.header for d in pc} == \
                            {sym for sym in self.potcar_symbols}:
                        return pc
            warnings.warn("No POTCAR file with matching TITEL fields"
                          " was found in {}".format(os.path.abspath(p)))

        if isinstance(path, string_types):
            if "POTCAR" in path:
                potcar = Potcar.from_file(path)
                if {d.TITEL for d in potcar} != \
                        {sym for sym in self.potcar_symbols}:
                    raise ValueError("Potcar TITELs do not match Vasprun")
            else:
                potcar = get_potcar_in_path(path)
        elif isinstance(path, bool) and path:
            potcar = get_potcar_in_path(os.path.split(self.filename)[0])
        else:
            potcar = None

        if potcar:
            self.potcar_spec = [{"titel": sym, "hash": ps.get_potcar_hash()}
                                for sym in self.potcar_symbols
                                for ps in potcar if
                                ps.symbol == sym.split()[1]]

    def as_dict(self):
        """
        Json-serializable dict representation.
        """
        d = {"vasp_version": self.vasp_version,
             "has_vasp_completed": self.converged,
             "nsites": len(self.final_structure)}
        comp = self.final_structure.composition
        d["unit_cell_formula"] = comp.as_dict()
        d["reduced_cell_formula"] = Composition(comp.reduced_formula).as_dict()
        d["pretty_formula"] = comp.reduced_formula
        symbols = [s.split()[1] for s in self.potcar_symbols]
        symbols = [re.split("_", s)[0] for s in symbols]
        d["is_hubbard"] = self.is_hubbard
        d["hubbards"] = {}
        if d["is_hubbard"]:
            us = self.incar.get("LDAUU", self.parameters.get("LDAUU"))
            js = self.incar.get("LDAUJ", self.parameters.get("LDAUJ"))
            if len(us) == len(symbols):
                d["hubbards"] = {symbols[i]: us[i] - js[i]
                                 for i in range(len(symbols))}
            else:
                raise VaspParserError("Length of U value parameters and atomic"
                                      " symbols are mismatched.")

        unique_symbols = sorted(list(set(self.atomic_symbols)))
        d["elements"] = unique_symbols
        d["nelements"] = len(unique_symbols)

        d["run_type"] = self.run_type

        vin = {"incar": {k: v for k, v in self.incar.items()},
               "crystal": self.initial_structure.as_dict(),
               "kpoints": self.kpoints.as_dict()}
        actual_kpts = [{"abc": list(self.actual_kpoints[i]),
                        "weight": self.actual_kpoints_weights[i]}
                       for i in range(len(self.actual_kpoints))]
        vin["kpoints"]["actual_points"] = actual_kpts
        vin["potcar"] = [s.split(" ")[1] for s in self.potcar_symbols]
        vin["potcar_spec"] = self.potcar_spec
        vin["potcar_type"] = [s.split(" ")[0] for s in self.potcar_symbols]
        vin["parameters"] = {k: v for k, v in self.parameters.items()}
        vin["lattice_rec"] = self.lattice_rec.as_dict()
        d["input"] = vin

        nsites = len(self.final_structure)

        try:
            vout = {"ionic_steps": self.ionic_steps,
                    "final_energy": self.final_energy,
                    "final_energy_per_atom": self.final_energy / nsites,
                    "crystal": self.final_structure.as_dict(),
                    "efermi": self.efermi}
        except (ArithmeticError, TypeError):
            vout = {"ionic_steps": self.ionic_steps,
                    "final_energy": self.final_energy,
                    "final_energy_per_atom": None,
                    "crystal": self.final_structure.as_dict(),
                    "efermi": self.efermi}

        if self.eigenvalues:
            eigen = defaultdict(dict)
            for (spin, index), values in self.eigenvalues.items():
                eigen[index][str(spin)] = values
            vout["eigenvalues"] = eigen
            (gap, cbm, vbm, is_direct) = self.eigenvalue_band_properties
            vout.update(dict(bandgap=gap, cbm=cbm, vbm=vbm,
                             is_gap_direct=is_direct))

            if self.projected_eigenvalues:
                peigen = []
                for i in range(len(eigen)):
                    peigen.append({})
                    for spin in eigen[i].keys():
                        peigen[i][spin] = []
                        for j in range(len(eigen[i][spin])):
                            peigen[i][spin].append({})
                for (spin, kpoint_index, band_index, ion_index, orbital), \
                        value in self.projected_eigenvalues.items():
                    beigen = peigen[kpoint_index][str(spin)][band_index]
                    if orbital not in beigen:
                        beigen[orbital] = [0.0] * nsites
                    beigen[orbital][ion_index] = value
                vout['projected_eigenvalues'] = peigen

        vout['epsilon_static'] = self.epsilon_static
        vout['epsilon_static_wolfe'] = self.epsilon_static_wolfe
        vout['epsilon_ionic'] = self.epsilon_ionic
        d['output'] = vout
        return jsanitize(d, strict=True)

    def _parse_params(self, elem):
        params = {}
        for c in elem:
            name = c.attrib.get("name")
            if c.tag not in ("i", "v"):
                p = self._parse_params(c)
                if name == "response functions":
                    # Delete duplicate fields from "response functions",
                    # which overrides the values in the root params.
                    p = {k: v for k, v in p.items() if k not in params}
                params.update(p)
            else:
                ptype = c.attrib.get("type")
                val = c.text.strip() if c.text else ""
                if c.tag == "i":
                    params[name] = _parse_parameters(ptype, val)
                else:
                    params[name] = _parse_v_parameters(ptype, val,
                                                       self.filename, name)
        elem.clear()
        return Incar(params)

    def _parse_atominfo(self, elem):
        for a in elem.findall("array"):
            if a.attrib["name"] == "atoms":
                atomic_symbols = [rc.find("c").text.strip()
                                  for rc in a.find("set")]
            elif a.attrib["name"] == "atomtypes":
                potcar_symbols = [rc.findall("c")[4].text.strip()
                                  for rc in a.find("set")]

        # ensure atomic symbols are valid elements
        def parse_atomic_symbol(symbol):
            try:
                return str(Element(symbol))
            # vasprun.xml uses X instead of Xe for xenon
            except KeyError as e:
                if symbol == "X":
                    return "Xe"
                raise e

        elem.clear()
        return [parse_atomic_symbol(sym) for
                sym in atomic_symbols], potcar_symbols

    def _parse_kpoints(self, elem):
        e = elem
        if elem.find("generation"):
            e = elem.find("generation")
        k = Kpoints("Kpoints from vasprun.xml")
        k.style = e.attrib["param"] if "param" in e.attrib else "Reciprocal"
        for v in e.findall("v"):
            name = v.attrib.get("name")
            toks = v.text.split()
            if name == "divisions":
                k.kpts = [[int(i) for i in toks]]
            elif name == "usershift":
                k.kpts_shift = [float(i) for i in toks]
            elif name in {"genvec1", "genvec2", "genvec3", "shift"}:
                setattr(k, name, [float(i) for i in toks])
        for va in elem.findall("varray"):
            name = va.attrib["name"]
            if name == "kpointlist":
                actual_kpoints = _parse_varray(va)
            elif name == "weights":
                weights = [i[0] for i in _parse_varray(va)]
        elem.clear()
        if k.style == "Reciprocal":
            k = Kpoints(comment="Kpoints from vasprun.xml",
                    style="Reciprocal", num_kpts=len(k.kpts),
                    kpts=actual_kpoints, kpts_weights=weights)
        return k, actual_kpoints, weights

    def _parse_structure(self, elem):
        latt = _parse_varray(elem.find("crystal").find("varray"))
        pos = _parse_varray(elem.find("varray"))
        return Structure(latt, self.atomic_symbols, pos)

    def _parse_diel(self, elem):
        imag = [[float(l) for l in r.text.split()] for r in elem.find("imag").find("array").find("set").findall("r")]
        real = [[float(l) for l in r.text.split()] for r in elem.find("real").find("array").find("set").findall("r")]
        return [e[0] for e in imag], [e[1:] for e in real], [e[1:] for e in imag]

    def _parse_calculation(self, elem):
        try:
            istep = {i.attrib["name"]: float(i.text)
                     for i in elem.find("energy").findall("i")}
        except AttributeError:  # not all calculations have an energy
            istep = {}
            pass
        esteps = []
        for scstep in elem.findall("scstep"):
            try:
                d = {i.attrib["name"]: _vasprun_float(i.text)
                     for i in scstep.find("energy").findall("i")}
                esteps.append(d)
            except AttributeError:  # not all calculations have an energy
                pass
        try:
            s = self._parse_structure(elem.find("structure"))
        except AttributeError:  # not all calculations have a structure
            s = None
            pass
        for va in elem.findall("varray"):
            istep[va.attrib["name"]] = _parse_varray(va)
        istep["electronic_steps"] = esteps
        istep["structure"] = s
        elem.clear()
        return istep

    def _parse_dos(self, elem):
        efermi = float(elem.find("i").text)
        energies = None
        tdensities = {}
        idensities = {}

        for s in elem.find("total").find("array").find("set").findall("set"):
            data = np.array(_parse_varray(s))
            energies = data[:, 0]
            spin = Spin.up if s.attrib["comment"] == "spin 1" else Spin.down
            tdensities[spin] = data[:, 1]
            idensities[spin] = data[:, 2]

        pdoss = []
        partial = elem.find("partial")
        if partial is not None:
            orbs = [ss.text for ss in partial.find("array").findall("field")]
            orbs.pop(0)
            lm = any(["x" in s for s in orbs])
            for s in partial.find("array").find("set").findall("set"):
                pdos = defaultdict(dict)

                for ss in s.findall("set"):
                    spin = Spin.up if ss.attrib["comment"] == "spin 1" else \
                        Spin.down
                    data = np.array(_parse_varray(ss))
                    nrow, ncol = data.shape
                    for j in range(1, ncol):
                        if lm:
                            orb = Orbital.from_vasp_index(j - 1)
                        else:
                            orb = orbs[j - 1].strip().upper()
                        pdos[orb][spin] = data[:, j]
                pdoss.append(pdos)
        elem.clear()
        return Dos(efermi, energies, tdensities), \
               Dos(efermi, energies, idensities), pdoss

    def _parse_eigen(self, elem):
        eigenvalues = {}
        for s in elem.find("array").find("set").findall("set"):
            spin = Spin.up if s.attrib["comment"] == "spin 1" else \
                Spin.down
            for i, ss in enumerate(s.findall("set")):
                eigenvalues[(spin, i)] = _parse_varray(ss)
        elem.clear()
        return eigenvalues

    def _parse_projected_eigen(self, elem):
        root = elem.find("array").find("set")
        proj_eigen = {}
        for s in root.findall("set"):
            spin = Spin.up if s.attrib["comment"] == "spin1" else \
                Spin.down
            for kpt, ss in enumerate(s.findall("set")):
                for band, sss in enumerate(ss.findall("set")):
                    for atom, data in enumerate(_parse_varray(sss)):
                        for i, v in enumerate(data):
                            orb = Orbital.from_vasp_index(i)
                            proj_eigen[(spin, kpt, band, atom, orb)] = v
        elem.clear()
        return proj_eigen


class Outcar(PMGSONable):
    """
    Parser for data in OUTCAR that is not available in Vasprun.xml

    Note, this class works a bit differently than most of the other
    VaspObjects, since the OUTCAR can be very different depending on which
    "type of run" performed.

    Creating the OUTCAR class with a filename reads "regular parameters" that
    are always present.

    Args:
        filename (str): OUTCAR filename to parse.

    .. attribute:: magnetization

        Magnetization on each ion as a tuple of dict, e.g.,
        ({"d": 0.0, "p": 0.003, "s": 0.002, "tot": 0.005}, ... )
        Note that this data is not always present.  LORBIT must be set to some
        other value than the default.

    .. attribute:: charge

        Charge on each ion as a tuple of dict, e.g.,
        ({"p": 0.154, "s": 0.078, "d": 0.0, "tot": 0.232}, ...)
        Note that this data is not always present.  LORBIT must be set to some
        other value than the default.

    .. attribute:: is_stopped

        True if OUTCAR is from a stopped run (using STOPCAR, see Vasp Manual).

    .. attribute:: run_stats

        Various useful run stats as a dict including "System time (sec)",
        "Total CPU time used (sec)", "Elapsed time (sec)",
        "Maximum memory used (kb)", "Average memory used (kb)",
        "User time (sec)".

    .. attribute:: elastic_tensor
        Total elastic moduli (Kbar) is given in a 6x6 array matrix.

    One can then call a specific reader depending on the type of run being
    performed. These are currently: read_igpar(), read_lepsilon() and
    read_lcalcpol(), read_core_state_eign().

    See the documentation of those methods for more documentation.

    Authors: Rickard Armiento, Shyue Ping Ong
    """
    def __init__(self, filename):
        self.filename = filename
        self.is_stopped = False

        # data from end of OUTCAR
        charge = []
        mag = []
        header = []
        run_stats = {}
        total_mag = None
        nelect = None
        efermi = None
        elastic_tensor = None

        time_patt = re.compile("\((sec|kb)\)")
        efermi_patt = re.compile("E-fermi\s*:\s*(\S+)")
        nelect_patt = re.compile("number of electron\s+(\S+)\s+"
                                 "magnetization\s+(\S+)")
        etensor_patt = re.compile("[X-Z][X-Z]+\s+-?\d+")

        all_lines = []
        for line in reverse_readfile(self.filename):
            clean = line.strip()
            all_lines.append(clean)
            if clean.find("soft stop encountered!  aborting job") != -1:
                self.is_stopped = True
            else:
                if time_patt.search(line):
                    tok = line.strip().split(":")
                    run_stats[tok[0].strip()] = float(tok[1].strip())
                    continue
                m = efermi_patt.search(clean)
                if m:
                    try:
                        #try-catch because VASP sometimes prints
                        #'E-fermi: ********     XC(G=0):  -6.1327
                        #alpha+bet : -1.8238'
                        efermi = float(m.group(1))
                        continue
                    except ValueError:
                        efermi = None
                        continue
                m = nelect_patt.search(clean)
                if m:
                    nelect = float(m.group(1))
                    total_mag = float(m.group(2))
            if all([nelect, total_mag is not None, efermi is not None,
                    run_stats]):
                break

        # For single atom systems, VASP doesn't print a total line, so
        # reverse parsing is very difficult
        read_charge = False
        read_mag = False
        all_lines.reverse()
        for clean in all_lines:
            if read_charge or read_mag:
                if clean.startswith("# of ion"):
                    header = re.split("\s{2,}", clean.strip())
                    header.pop(0)
                else:
                    m = re.match("\s*(\d+)\s+(([\d\.\-]+)\s+)+", clean)
                    if m:
                        toks = [float(i) for i in re.findall("[\d\.\-]+", clean)]
                        toks.pop(0)
                        if read_charge:
                            charge.append(dict(zip(header, toks)))
                        else:
                            mag.append(dict(zip(header, toks)))
                    elif clean.startswith('tot'):
                        read_charge = False
                        read_mag = False
            if clean == "total charge":
                charge = []
                read_charge = True
                read_mag = False
            elif clean == "magnetization (x)":
                mag = []
                read_mag = True
                read_charge = False

        # data from beginning of OUTCAR
        run_stats['cores'] = 0
        with zopen(filename, "rt") as f:
            for line in f:
                if "running" in line:
                    run_stats['cores'] = line.split()[2]
                    break

        # 6x6 tensor matrix for TOTAL ELASTIC MODULI
        tensor_matrix = []
        tag = "TOTAL ELASTIC MODULI (kBar)"
        if tag in all_lines:
            for clean in all_lines:
                if etensor_patt.search(clean):
                    tok = clean.strip().split()
                    tok.pop(0)
                    tok = [float(i) for i in tok]
                    tensor_matrix.append(tok)
            total_elm = [tensor_matrix[i] for i in range(18, 24)]
            elastic_tensor = np.asarray(total_elm).reshape(6, 6)
        else:
            pass

        self.run_stats = run_stats
        self.magnetization = tuple(mag)
        self.charge = tuple(charge)
        self.efermi = efermi
        self.nelect = nelect
        self.total_mag = total_mag
        self.elastic_tensor = elastic_tensor
        self.data = {}

    def read_pattern(self, patterns, reverse=False, terminate_on_match=False,
                     postprocess=str):
        """
        General pattern reading. Uses monty's regrep method. Takes the same
        arguments.

        Args:
            patterns (dict): A dict of patterns, e.g.,
                {"energy": "energy\(sigma->0\)\s+=\s+([\d\-\.]+)"}.
            reverse (bool): Read files in reverse. Defaults to false. Useful for
                large files, esp OUTCARs, especially when used with
                terminate_on_match.
            terminate_on_match (bool): Whether to terminate when there is at
                least one match in each key in pattern.
            postprocess (callable): A post processing function to convert all
                matches. Defaults to str, i.e., no change.

        Renders accessible:
            Any attribute in patterns. For example,
            {"energy": "energy\(sigma->0\)\s+=\s+([\d\-\.]+)"} will set the
            value of self.data["energy"] = [[-1234], [-3453], ...], to the
            results from regex and postprocess. Note that the returned values
            are lists of lists, because you can grep multiple items on one line.
        """
        matches = regrep(self.filename, patterns, reverse=reverse,
                         terminate_on_match=terminate_on_match,
                         postprocess=postprocess)
        for k in patterns.keys():
            self.data[k] = [i[0] for i in matches.get(k, [])]

    def read_neb(self, reverse=True, terminate_on_match=True):
        """
        Reads NEB data. This only works with OUTCARs from both normal
        VASP NEB calculations or from the CI NEB method implemented by
        Henkelman et al.

        Args:
            reverse (bool): Read files in reverse. Defaults to false. Useful for
                large files, esp OUTCARs, especially when used with
                terminate_on_match. Defaults to True here since we usually
                want only the final value.
            terminate_on_match (bool): Whether to terminate when there is at
                least one match in each key in pattern. Defaults to True here
                since we usually want only the final value.

        Renders accessible:
            tangent_force - Final tangent force.
            energy - Final energy.
            These can be accessed under Outcar.data[key]
        """
        patterns = {
            "energy": "energy\(sigma->0\)\s+=\s+([\d\-\.]+)",
            "tangent_force": "(NEB: projections on to tangent \(" \
                "spring, REAL\)\s+\S+|tangential force \(eV/A\))\s+(["
                                   "\d\-\.]+)"
        }
        self.read_pattern(patterns, reverse=reverse,
                          terminate_on_match=terminate_on_match,
                          postprocess=str)
        self.data["energy"] = float(self.data["energy"][0][0])
        if self.data.get("tangent_force"):
            self.data["tangent_force"] = float(self.data["tangent_force"][0][1])

    def read_igpar(self):
        """
        Renders accessible:
            er_ev = e<r>_ev (dictionary with Spin.up/Spin.down as keys)
            er_bp = e<r>_bp (dictionary with Spin.up/Spin.down as keys)
            er_ev_tot = spin up + spin down summed
            er_bp_tot = spin up + spin down summed
            p_elc = spin up + spin down summed
            p_ion = spin up + spin down summed

        (See VASP section "LBERRY,  IGPAR,  NPPSTR,  DIPOL tags" for info on
        what these are).
        """

        # variables to be filled
        self.er_ev = {}  # will  be  dict (Spin.up/down) of array(3*float)
        self.er_bp = {}  # will  be  dics (Spin.up/down) of array(3*float)
        self.er_ev_tot = None  # will be array(3*float)
        self.er_bp_tot = None  # will be array(3*float)
        self.p_elec = None
        self.p_ion = None
        try:
            search = []

            # Nonspin cases
            def er_ev(results, match):
                results.er_ev[Spin.up] = np.array(map(float,
                                                      match.groups()[1:4])) / 2
                results.er_ev[Spin.down] = results.er_ev[Spin.up]
                results.context = 2

            search.append(["^ *e<r>_ev=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           "*([-0-9.Ee+]*) *\)",
                           None, er_ev])

            def er_bp(results, match):
                results.er_bp[Spin.up] = np.array([float(match.group(i))
                                                   for i in range(1, 4)]) / 2
                results.er_bp[Spin.down] = results.er_bp[Spin.up]

            search.append(["^ *e<r>_bp=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           "*([-0-9.Ee+]*) *\)",
                           lambda results, line: results.context == 2, er_bp])

            # Spin cases
            def er_ev_up(results, match):
                results.er_ev[Spin.up] = np.array([float(match.group(i))
                                                   for i in range(1, 4)])
                results.context = Spin.up

            search.append(["^.*Spin component 1 *e<r>_ev=\( *([-0-9.Ee+]*) "
                           "*([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)",
                           None, er_ev_up])

            def er_bp_up(results, match):
                results.er_bp[Spin.up] = np.array([float(match.group(1)),
                                                   float(match.group(2)),
                                                   float(match.group(3))])

            search.append(["^ *e<r>_bp=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           "*([-0-9.Ee+]*) *\)",
                           lambda results,
                           line: results.context == Spin.up, er_bp_up])

            def er_ev_dn(results, match):
                results.er_ev[Spin.down] = np.array([float(match.group(1)),
                                                     float(match.group(2)),
                                                     float(match.group(3))])
                results.context = Spin.down
            search.append(["^.*Spin component 2 *e<r>_ev=\( *([-0-9.Ee+]*) "
                           "*([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)",
                           None, er_ev_dn])

            def er_bp_dn(results, match):
                results.er_bp[Spin.down] = np.array([float(match.group(i))
                                                     for i in range(1, 4)])
            search.append(["^ *e<r>_bp=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           "*([-0-9.Ee+]*) *\)",
                           lambda results,
                           line: results.context == Spin.down, er_bp_dn])

            # Always present spin/non-spin
            def p_elc(results, match):
                results.p_elc = np.array([float(match.group(i))
                                          for i in range(1, 4)])

            search.append(["^.*Total electronic dipole moment: "
                           "*p\[elc\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           "*([-0-9.Ee+]*) *\)", None, p_elc])

            def p_ion(results, match):
                results.p_ion = np.array([float(match.group(i))
                                          for i in range(1, 4)])

            search.append(["^.*ionic dipole moment: "
                           "*p\[ion\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           "*([-0-9.Ee+]*) *\)", None, p_ion])

            self.context = None
            self.er_ev = {Spin.up: None, Spin.down: None}
            self.er_bp = {Spin.up: None, Spin.down: None}

            micro_pyawk(self.filename, search, self)

            if self.er_ev[Spin.up] is not None and \
                    self.er_ev[Spin.down] is not None:
                self.er_ev_tot = self.er_ev[Spin.up] + self.er_ev[Spin.down]

            if self.er_bp[Spin.up] is not None and \
                    self.er_bp[Spin.down] is not None:
                self.er_bp_tot = self.er_bp[Spin.up] + self.er_bp[Spin.down]

        except:
            self.er_ev_tot = None
            self.er_bp_tot = None
            raise Exception("IGPAR OUTCAR could not be parsed.")

    def read_lepsilon(self):
        # variables to be filled
        try:
            search = []

            def dielectric_section_start(results, match):
                results.dielectric_index = -1

            search.append(["MACROSCOPIC STATIC DIELECTRIC TENSOR \(", None,
                           dielectric_section_start])

            def dielectric_section_start2(results, match):
                results.dielectric_index = 0

            search.append(
                ["-------------------------------------",
                lambda results, line: results.dielectric_index == -1,
                dielectric_section_start2])

            def dielectric_data(results, match):
                results.dielectric_tensor[results.dielectric_index, :] = \
                    np.array([float(match.group(i)) for i in range(1, 4)])
                results.dielectric_index += 1

            search.append(
                ["^ *([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+) *$",
                lambda results, line: results.dielectric_index >= 0
                                      if results.dielectric_index is not None
                                      else None,
                dielectric_data])

            def dielectric_section_stop(results, match):
                results.dielectric_index = None

            search.append(
                ["-------------------------------------",
                lambda results, line: results.dielectric_index >= 1
                                      if results.dielectric_index is not None
                                      else None,
                dielectric_section_stop])

            self.dielectric_index = None
            self.dielectric_tensor = np.zeros((3, 3))

            def piezo_section_start(results, match):
                results.piezo_index = 0

            search.append(["PIEZOELECTRIC TENSOR  for field in x, y, z        "
                           "\(C/m\^2\)",
                           None, piezo_section_start])

            def piezo_data(results, match):
                results.piezo_tensor[results.piezo_index, :] = \
                    np.array([float(match.group(i)) for i in range(1, 7)])
                results.piezo_index += 1

            search.append(
                ["^ *[xyz] +([-0-9.Ee+]+) +([-0-9.Ee+]+)" +
                 " +([-0-9.Ee+]+) *([-0-9.Ee+]+) +([-0-9.Ee+]+)" +
                 " +([-0-9.Ee+]+)*$",
                 lambda results, line: results.piezo_index >= 0
                                       if results.piezo_index is not None
                                       else None,
                 piezo_data])

            def piezo_section_stop(results, match):
                results.piezo_index = None

            search.append(
                ["-------------------------------------",
                lambda results, line: results.piezo_index >= 1
                                      if results.piezo_index is not None
                                      else None,
                piezo_section_stop])

            self.piezo_index = None
            self.piezo_tensor = np.zeros((3, 6))

            def born_section_start(results, match):
                results.born_ion = -1

            search.append(["BORN EFFECTIVE CHARGES " +
                           "\(in e, cummulative output\)",
                           None, born_section_start])

            def born_ion(results, match):
                results.born_ion = int(match.group(1)) - 1
                results.born[results.born_ion] = np.zeros((3, 3))

            search.append(["ion +([0-9]+)", lambda results,
                           line: results.born_ion is not None, born_ion])

            def born_data(results, match):
                results.born[results.born_ion][int(match.group(1)) - 1, :] = \
                    np.array([float(match.group(i)) for i in range(2, 5)])

            search.append(
                ["^ *([1-3]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+)$",
                lambda results, line: results.born_ion >= 0
                                      if results.born_ion is not None
                                      else results.born_ion,
                born_data])

            def born_section_stop(results, match):
                results.born_index = None

            search.append(
                ["-------------------------------------",
                lambda results, line: results.born_ion >= 1
                                      if results.born_ion is not None
                                      else results.born_ion,
                born_section_stop])

            self.born_ion = None
            self.born = {}

            micro_pyawk(self.filename, search, self)

            self.dielectric_tensor = self.dielectric_tensor.tolist()
            self.piezo_tensor = self.piezo_tensor.tolist()

        except:
            raise Exception("LEPSILON OUTCAR could not be parsed.")

    def read_lepsilon_ionic(self):
        # variables to be filled
        try:
            search = []

            def dielectric_section_start(results, match):
                results.dielectric_ionic_index = -1

            search.append(["MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC", None,
                           dielectric_section_start])

            def dielectric_section_start2(results, match):
                results.dielectric_ionic_index = 0

            search.append(
                ["-------------------------------------",
                lambda results, line: results.dielectric_ionic_index == -1
                                      if results.dielectric_ionic_index is not None
                                      else results.dielectric_ionic_index,
                dielectric_section_start2])

            def dielectric_data(results, match):
                results.dielectric_ionic_tensor[results.dielectric_ionic_index, :] = \
                    np.array([float(match.group(i)) for i in range(1, 4)])
                results.dielectric_ionic_index += 1

            search.append(
                ["^ *([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+) *$",
                lambda results, line: results.dielectric_ionic_index >= 0
                                      if results.dielectric_ionic_index is not None
                                      else results.dielectric_ionic_index,
                dielectric_data])

            def dielectric_section_stop(results, match):
                results.dielectric_ionic_index = None

            search.append(
                ["-------------------------------------",
                lambda results, line: results.dielectric_ionic_index >= 1
                                      if results.dielectric_ionic_index is not None
                                      else results.dielectric_ionic_index,
                dielectric_section_stop])

            self.dielectric_ionic_index = None
            self.dielectric_ionic_tensor = np.zeros((3, 3))

            def piezo_section_start(results, match):
                results.piezo_ionic_index = 0

            search.append(["PIEZOELECTRIC TENSOR IONIC CONTR  for field in x, y, z        ",
                           None, piezo_section_start])

            def piezo_data(results, match):
                results.piezo_ionic_tensor[results.piezo_ionic_index, :] = \
                    np.array([float(match.group(i)) for i in range(1, 7)])
                results.piezo_ionic_index += 1

            search.append(
                ["^ *[xyz] +([-0-9.Ee+]+) +([-0-9.Ee+]+)" +
                 " +([-0-9.Ee+]+) *([-0-9.Ee+]+) +([-0-9.Ee+]+)" +
                 " +([-0-9.Ee+]+)*$",
                 lambda results, line: results.piezo_ionic_index >= 0
                                       if results.piezo_ionic_index is not None
                                       else results.piezo_ionic_index,
                 piezo_data])

            def piezo_section_stop(results, match):
                results.piezo_ionic_index = None

            search.append(
                ["-------------------------------------",
                 lambda results, line: results.piezo_ionic_index >= 1
                                       if results.piezo_ionic_index is not None
                                       else results.piezo_ionic_index,
                 piezo_section_stop])

            self.piezo_ionic_index = None
            self.piezo_ionic_tensor = np.zeros((3, 6))

            micro_pyawk(self.filename, search, self)

            self.dielectric_ionic_tensor = self.dielectric_ionic_tensor.tolist()
            self.piezo_ionic_tensor = self.piezo_ionic_tensor.tolist()

        except:
            raise Exception("ionic part of LEPSILON OUTCAR could not be parsed.")

    def read_lcalcpol(self):
        # variables to be filled
        self.p_elec = None
        self.p_ion = None
        try:
            search = []

            # Always present spin/non-spin
            def p_elc(results, match):
                results.p_elc = np.array([float(match.group(1)),
                                          float(match.group(2)),
                                          float(match.group(3))])

            search.append(["^.*Total electronic dipole moment: "
                           "*p\[elc\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           "*([-0-9.Ee+]*) *\)",
                           None, p_elc])

            def p_ion(results, match):
                results.p_ion = np.array([float(match.group(1)),
                                          float(match.group(2)),
                                          float(match.group(3))])
            search.append(["^.*Ionic dipole moment: *p\[ion\]="
                           "\( *([-0-9.Ee+]*)"
                           " *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)",
                           None, p_ion])

            micro_pyawk(self.filename, search, self)

        except:
            raise Exception("CLACLCPOL OUTCAR could not be parsed.")

    def read_core_state_eigen(self):
        """
        Read the core state eigenenergies at each ionic step.

        Returns:
            A list of dict over the atom such as [{"AO":[core state eig]}].
            The core state eigenenergie list for each AO is over all ionic
            step.

        Example:
            The core state eigenenergie of the 2s AO of the 6th atom of the
            structure at the last ionic step is [5]["2s"][-1]
        """

        with zopen(self.filename, "rt") as foutcar:
            line = foutcar.readline()
            while line != "":
                line = foutcar.readline()
                if "NIONS =" in line:
                    natom = int(line.split("NIONS =")[1])
                    cl = [defaultdict(list) for i in range(natom)]
                if "the core state eigen" in line:
                    for iat in range(natom):
                        line = foutcar.readline()
                        data = line.split()[1:]
                        for i in range(0, len(data), 2):
                            cl[iat][data[i]].append(float(data[i+1]))
        return cl

    def as_dict(self):
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__, "efermi": self.efermi,
             "run_stats": self.run_stats, "magnetization": self.magnetization,
             "charge": self.charge, "total_magnetization": self.total_mag,
             "nelect": self.nelect, "is_stopped": self.is_stopped}
        return d


class VolumetricData(object):
    """
    Simple volumetric object for reading LOCPOT and CHGCAR type files.

    .. attribute:: structure

        Structure associated with the Volumetric Data object

    ..attribute:: is_spin_polarized

        True if run is spin polarized

    ..attribute:: dim

        Tuple of dimensions of volumetric grid in each direction (nx, ny, nz).

    ..attribute:: data

        Actual data as a dict of {string: np.array}. The string are "total"
        and "diff", in accordance to the output format of vasp LOCPOT and
        CHGCAR files where the total spin density is written first, followed
        by the difference spin density.

    .. attribute:: ngridpts

        Total number of grid points in volumetric data.
    """
    def __init__(self, structure, data, distance_matrix=None):
        """
        Typically, this constructor is not used directly and the static
        from_file constructor is used. This constructor is designed to allow
        summation and other operations between VolumetricData objects.

        Args:
            structure: Structure associated with the volumetric data
            data: Actual volumetric data.
            distance_matrix: A pre-computed distance matrix if available.
                Useful so pass distance_matrices between sums,
                shortcircuiting an otherwise expensive operation.
        """
        self.structure = structure
        self.is_spin_polarized = len(data) == 2
        self.dim = data["total"].shape
        self.data = data
        self.ngridpts = self.dim[0] * self.dim[1] * self.dim[2]
        #lazy init the spin data since this is not always needed.
        self._spin_data = {}
        self._distance_matrix = {} if not distance_matrix else distance_matrix

    @property
    def spin_data(self):
        """
        The data decomposed into actual spin data as {spin: data}.
        Essentially, this provides the actual Spin.up and Spin.down data
        instead of the total and diff.  Note that by definition, a
        non-spin-polarized run would have Spin.up data == Spin.down data.
        """
        if not self._spin_data:
            spin_data = dict()
            spin_data[Spin.up] = 0.5 * (self.data["total"] +
                                        self.data.get("diff", 0))
            spin_data[Spin.down] = 0.5 * (self.data["total"] -
                                          self.data.get("diff", 0))
            self._spin_data = spin_data
        return self._spin_data

    def get_axis_grid(self, ind):
        """
        Returns the grid for a particular axis.

        Args:
            ind (int): Axis index.
        """
        ng = self.dim
        num_pts = ng[ind]
        lengths = self.structure.lattice.abc
        return [i / num_pts * lengths[ind] for i in range(num_pts)]

    def __add__(self, other):
        return self.linear_add(other, 1.0)

    def __sub__(self, other):
        return self.linear_add(other, -1.0)

    def linear_add(self, other, scale_factor=1.0):
        """
        Method to do a linear sum of volumetric objects. Used by + and -
        operators as well. Returns a VolumetricData object containing the
        linear sum.

        Args:
            other (VolumetricData): Another VolumetricData object
            scale_factor (float): Factor to scale the other data by.

        Returns:
            VolumetricData corresponding to self + scale_factor * other.
        """
        if self.structure != other.structure:
            raise ValueError("Adding or subtraction operations can only be "
                             "performed for volumetric data with the exact "
                             "same structure.")
        #To add checks
        data = {}
        for k in self.data.keys():
            data[k] = self.data[k] + scale_factor * other.data[k]
        return VolumetricData(self.structure, data, self._distance_matrix)

    @staticmethod
    def parse_file(filename):
        """
        Convenience method to parse a generic volumetric data file in the vasp
        like format. Used by subclasses for parsing file.

        Args:
            filename (str): Path of file to parse

        Returns:
            (poscar, data)
        """
        poscar_read = False
        poscar_string = []
        dataset = []
        all_dataset = []
        dim = None
        dimline = None
        read_dataset = False
        ngrid_pts = 0
        data_count = 0
        poscar = None
        with zopen(filename) as f:
            for line in f:
                line = line.strip()
                if read_dataset:
                    toks = line.split()
                    for tok in toks:
                        if data_count < ngrid_pts:
                            #This complicated procedure is necessary because
                            #vasp outputs x as the fastest index, followed by y
                            #then z.
                            x = data_count % dim[0]
                            y = int(math.floor(data_count / dim[0])) % dim[1]
                            z = int(math.floor(data_count / dim[0] / dim[1]))
                            dataset[x, y, z] = float(tok)
                            data_count += 1
                    if data_count >= ngrid_pts:
                        read_dataset = False
                        data_count = 0
                        all_dataset.append(dataset)
                elif not poscar_read:
                    if line != "" or len(poscar_string) == 0:
                        poscar_string.append(line)
                    elif line == "":
                        poscar = Poscar.from_string("\n".join(poscar_string))
                        poscar_read = True
                elif not dim:
                    dim = [int(i) for i in line.split()]
                    ngrid_pts = dim[0] * dim[1] * dim[2]
                    dimline = line
                    read_dataset = True
                    dataset = np.zeros(dim)
                elif line == dimline:
                    read_dataset = True
                    dataset = np.zeros(dim)
            if len(all_dataset) == 2:
                data = {"total": all_dataset[0], "diff": all_dataset[1]}
            else:
                data = {"total": all_dataset[0]}
            return poscar, data

    def write_file(self, file_name, vasp4_compatible=False):
        """
        Write the VolumetricData object to a vasp compatible file.

        Args:
            file_name (str): Path to a file
            vasp4_compatible (bool): True if the format is vasp4 compatible
        """

        with zopen(file_name, "wt") as f:
            p = Poscar(self.structure)

            lines = p.comment + "\n"
            lines += "   1.00000000000000\n"
            latt = self.structure.lattice.matrix
            lines += " %12.6f%12.6f%12.6f\n" % tuple(latt[0,:])
            lines += " %12.6f%12.6f%12.6f\n" % tuple(latt[1,:])
            lines += " %12.6f%12.6f%12.6f\n" % tuple(latt[2,:])
            if not vasp4_compatible:
                lines += "".join(["%5s" % s for s in p.site_symbols]) + "\n"
            lines += "".join(["%6d" % x for x in p.natoms]) + "\n"
            lines += "Direct\n"
            for site in self.structure:
                lines += "%10.6f%10.6f%10.6f\n" % tuple(site.frac_coords)
            lines += "\n"
            f.write(lines)
            a = self.dim

            def write_spin(data_type):
                lines = []
                count = 0
                f.write("{} {} {}\n".format(a[0], a[1], a[2]))
                for (k, j, i) in itertools.product(list(range(a[2])), list(range(a[1])),
                                                   list(range(a[0]))):
                    lines.append("%0.11e" % self.data[data_type][i, j, k])
                    count += 1
                    if count % 5 == 0:
                        f.write("".join(lines) + "\n")
                        lines = []
                    else:
                        lines.append(" ")
                f.write("".join(lines) + "\n")

            write_spin("total")
            if self.is_spin_polarized:
                f.write("\n")
                write_spin("diff")

    def get_integrated_diff(self, ind, radius, nbins=1):
        """
        Get integrated difference of atom index ind up to radius. This can be
        an extremely computationally intensive process, depending on how many
        grid points are in the VolumetricData.

        Args:
            ind (int): Index of atom.
            radius (float): Radius of integration.
            nbins (int): Number of bins. Defaults to 1. This allows one to
                obtain the charge integration up to a list of the cumulative
                charge integration values for radii for [radius/nbins,
                2 * radius/nbins, ....].

        Returns:
            Differential integrated charge as a np array of [[radius, value],
            ...]. Format is for ease of plotting. E.g., plt.plot(data[:,0],
            data[:,1])
        """
        #For non-spin-polarized runs, this is zero by definition.
        if not self.is_spin_polarized:
            radii = [radius / nbins * (i + 1) for i in range(nbins)]
            data = np.zeros((nbins, 2))
            data[:, 0] = radii
            return data

        struct = self.structure
        a = self.dim
        if ind not in self._distance_matrix or\
                self._distance_matrix[ind]["max_radius"] < radius:
            coords = []
            for (x, y, z) in itertools.product(*[list(range(i)) for i in a]):
                coords.append([x / a[0], y / a[1], z / a[2]])
            sites_dist = struct.lattice.get_points_in_sphere(
                coords, struct[ind].coords, radius)
            self._distance_matrix[ind] = {"max_radius": radius,
                                          "data": np.array(sites_dist)}

        data = self._distance_matrix[ind]["data"]

        #Use boolean indexing to find all charges within the desired distance.
        inds = data[:, 1] <= radius
        dists = data[inds, 1]
        data_inds = np.rint(np.mod(list(data[inds, 0]), 1) *
                            np.tile(a, (len(dists), 1))).astype(int)
        vals = [self.data["diff"][x, y, z] for x, y, z in data_inds]

        hist, edges = np.histogram(dists, bins=nbins,
                                   range=[0, radius],
                                   weights=vals)
        data = np.zeros((nbins, 2))
        data[:, 0] = edges[1:]
        data[:, 1] = [sum(hist[0:i + 1]) / self.ngridpts
                      for i in range(nbins)]
        return data

    def get_average_along_axis(self, ind):
        """
        Get the averaged total of the volumetric data a certain axis direction.
        For example, useful for visualizing Hartree Potentials from a LOCPOT
        file.

        Args:
            ind (int): Index of axis.

        Returns:
            Average total along axis
        """
        m = self.data["total"]
        ng = self.dim
        if ind == 0:
            total = np.sum(np.sum(m, axis=1), 1)
        elif ind == 1:
            total = np.sum(np.sum(m, axis=0), 1)
        else:
            total = np.sum(np.sum(m, axis=0), 0)
        return total / ng[(ind + 1) % 3] / ng[(ind + 2) % 3]


class Locpot(VolumetricData):
    """
    Simple object for reading a LOCPOT file.

    Args:
        poscar (Poscar): Poscar object containing structure.
        data: Actual data.
    """

    def __init__(self, poscar, data):
        super(Locpot, self).__init__(poscar.structure, data)
        self.name = poscar.comment

    @staticmethod
    def from_file(filename):
        (poscar, data) = VolumetricData.parse_file(filename)
        return Locpot(poscar, data)


class Chgcar(VolumetricData):
    """
    Simple object for reading a CHGCAR file.

    Args:
        poscar (Poscar): Poscar object containing structure.
        data: Actual data.
    """

    def __init__(self, poscar, data):
        super(Chgcar, self).__init__(poscar.structure, data)
        self.poscar = poscar
        self.name = poscar.comment
        self._distance_matrix = {}

    @staticmethod
    def from_file(filename):
        (poscar, data) = VolumetricData.parse_file(filename)
        return Chgcar(poscar, data)


class Procar(object):
    """
    Object for reading a PROCAR file.

    Args:
        filename: Name of file containing PROCAR.

    .. attribute:: data

        A nested dict containing the PROCAR data of the form below. It should
        be noted that VASP uses 1-based indexing for atoms, but this is
        converted to zero-based indexing in this parser to be consistent with
        representation of structures in pymatgen::

            {
                atom_index: {
                    kpoint_index: {
                        "bands": {
                            band_index: {
                                "p": 0.002,
                                "s": 0.025,
                                "d": 0.0
                            },
                            ...
                        },
                        "weight": 0.03125
                    },
                    ...
            }
    """
    def __init__(self, filename):
        data = defaultdict(dict)
        headers = None
        with zopen(filename, "rt") as f:
            lines = list(clean_lines(f.readlines()))
            self.name = lines[0]
            kpointexpr = re.compile("^\s*k-point\s+(\d+).*weight = ([0-9\.]+)")
            bandexpr = re.compile("^\s*band\s+(\d+)")
            ionexpr = re.compile("^ion.*")
            expr = re.compile("^\s*([0-9]+)\s+")
            dataexpr = re.compile("[\.0-9]+")
            weight = 0
            current_kpoint = 0
            current_band = 0
            for l in lines:
                if bandexpr.match(l):
                    m = bandexpr.match(l)
                    current_band = int(m.group(1))
                elif kpointexpr.match(l):
                    m = kpointexpr.match(l)
                    current_kpoint = int(m.group(1))
                    weight = float(m.group(2))
                elif headers is None and ionexpr.match(l):
                    headers = l.split()
                    headers.pop(0)
                    headers.pop(-1)
                elif expr.match(l):
                    linedata = dataexpr.findall(l)
                    num_data = [float(i) for i in linedata]
                    #Convert to zero-based indexing for atoms.
                    index = int(num_data.pop(0)) - 1
                    num_data.pop(-1)
                    if current_kpoint not in data[index]:
                        data[index][current_kpoint] = {"weight": weight,
                                                       "bands": {}}
                    data[index][current_kpoint]["bands"][current_band] = \
                        dict(zip(headers, num_data))
            self.data = data
            self._nb_kpoints = len(data[0].keys())
            self._nb_bands = len(data[0][1]["bands"].keys())

    @property
    def nb_bands(self):
        """
        returns the number of bands in the band structure
        """
        return self._nb_bands

    @property
    def nb_kpoints(self):
        """
        Returns the number of k-points in the band structure calculation
        """
        return self._nb_kpoints

    def get_projection_on_elements(self, structure):
        """
        Method returning a dictionary of projections on elements.
        Spin polarized calculation are not supported.

        Args:
            structure (Structure): Input structure.

        Returns:
            a dictionary in the {Spin.up:[k index][b index][{Element:values}]]
        """
        dico = {Spin.up: []}
        dico[Spin.up] = [[defaultdict(float)
                          for i in range(self._nb_kpoints)]
                         for j in range(self.nb_bands)]

        for iat in self.data:
            name = structure.species[iat].symbol
            for k in self.data[iat]:
                for b in self.data[iat][k]["bands"]:
                    dico[Spin.up][b-1][k-1][name] = \
                        sum(self.data[iat][k]["bands"][b].values())

        return dico

    def get_occupation(self, atom_index, orbital):
        """
        Returns the occupation for a particular orbital of a particular atom.

        Args:
            atom_num (int): Index of atom in the PROCAR. It should be noted
                that VASP uses 1-based indexing for atoms, but this is
                converted to 0-based indexing in this parser to be
                consistent with representation of structures in pymatgen.
            orbital (str): An orbital. If it is a single character, e.g., s,
                p, d or f, the sum of all s-type, p-type, d-type or f-type
                orbitals occupations are returned respectively. If it is a
                specific orbital, e.g., px, dxy, etc., only the occupation
                of that orbital is returned.

        Returns:
            Sum occupation of orbital of atom.
        """
        total = 0
        found = False
        for kpoint, d in self.data[atom_index].items():
            wt = d["weight"]
            for band, dd in d["bands"].items():
                for orb, v in dd.items():
                    if orb.startswith(orbital):
                        found = True
                        total += v * wt
        if not found:
            raise ValueError("Invalid orbital {}".format(orbital))
        return total


class Oszicar(object):
    """
    A basic parser for an OSZICAR output from VASP.  In general, while the
    OSZICAR is useful for a quick look at the output from a VASP run, we
    recommend that you use the Vasprun parser instead, which gives far richer
    information about a run.

    Args:
        filename (str): Filename of file to parse

    .. attribute:: electronic_steps

            All electronic steps as a list of list of dict. e.g.,
            [[{"rms": 160.0, "E": 4507.24605593, "dE": 4507.2, "N": 1,
            "deps": -17777.0, "ncg": 16576}, ...], [....]
            where electronic_steps[index] refers the list of electronic steps
            in one ionic_step, electronic_steps[index][subindex] refers to a
            particular electronic step at subindex in ionic step at index. The
            dict of properties depends on the type of VASP run, but in general,
            "E", "dE" and "rms" should be present in almost all runs.

    .. attribute:: ionic_steps:

            All ionic_steps as a list of dict, e.g.,
            [{"dE": -526.36, "E0": -526.36024, "mag": 0.0, "F": -526.36024},
            ...]
            This is the typical output from VASP at the end of each ionic step.
    """

    def __init__(self, filename):
        electronic_steps = []
        ionic_steps = []
        ionic_pattern = re.compile("(\d+)\s+F=\s*([\d\-\.E\+]+)\s+"
                                   "E0=\s*([\d\-\.E\+]+)\s+"
                                   "d\s*E\s*=\s*([\d\-\.E\+]+)$")
        ionic_mag_pattern = re.compile("(\d+)\s+F=\s*([\d\-\.E\+]+)\s+"
                                       "E0=\s*([\d\-\.E\+]+)\s+"
                                       "d\s*E\s*=\s*([\d\-\.E\+]+)\s+"
                                       "mag=\s*([\d\-\.E\+]+)")
        ionic_MD_pattern = re.compile("(\d+)\s+T=\s*([\d\-\.E\+]+)\s+"
                                      "E=\s*([\d\-\.E\+]+)\s+"
                                      "F=\s*([\d\-\.E\+]+)\s+"
                                      "E0=\s*([\d\-\.E\+]+)\s+"
                                      "EK=\s*([\d\-\.E\+]+)\s+"
                                      "SP=\s*([\d\-\.E\+]+)\s+"
                                      "SK=\s*([\d\-\.E\+]+)")
        electronic_pattern = re.compile("\s*\w+\s*:(.*)")

        def smart_convert(header, num):
            try:
                if header == "N" or header == "ncg":
                    v = int(num)
                    return v
                v = float(num)
                return v
            except ValueError:
                return "--"

        header = []
        with zopen(filename, "rt") as fid:
            for line in fid:
                line = line.strip()
                m = electronic_pattern.match(line)
                if m:
                    toks = m.group(1).split()
                    data = {header[i]: smart_convert(header[i], toks[i])
                            for i in range(len(toks))}
                    if toks[0] == "1":
                        electronic_steps.append([data])
                    else:
                        electronic_steps[-1].append(data)
                elif ionic_pattern.match(line.strip()):
                    m = ionic_pattern.match(line.strip())
                    ionic_steps.append({"F": float(m.group(2)),
                                        "E0": float(m.group(3)),
                                        "dE": float(m.group(4))})
                elif ionic_mag_pattern.match(line.strip()):
                    m = ionic_mag_pattern.match(line.strip())
                    ionic_steps.append({"F": float(m.group(2)),
                                        "E0": float(m.group(3)),
                                        "dE": float(m.group(4)),
                                        "mag": float(m.group(5))})
                elif ionic_MD_pattern.match(line.strip()):
                    m = ionic_MD_pattern.match(line.strip())
                    ionic_steps.append({"T": float(m.group(2)),
                                        "E": float(m.group(3)),
                                        "F": float(m.group(4)),
                                        "E0": float(m.group(5)),
                                        "EK": float(m.group(6)),
                                        "SP": float(m.group(7)),
                                        "SK": float(m.group(8))})
                elif re.match("^\s*N\s+E\s*", line):
                    header = line.strip().replace("d eps", "deps").split()
        self.electronic_steps = electronic_steps
        self.ionic_steps = ionic_steps

    @property
    def all_energies(self):
        """
        Compilation of all energies from all electronic steps and ionic steps
        as a tuple of list of energies, e.g.,
        ((4507.24605593, 143.824705755, -512.073149912, ...), ...)
        """
        all_energies = []
        for i in range(len(self.electronic_steps)):
            energies = [step["E"] for step in self.electronic_steps[i]]
            energies.append(self.ionic_steps[i]["F"])
            all_energies.append(tuple(energies))
        return tuple(all_energies)

    @property
    @unitized("eV")
    def final_energy(self):
        """
        Final energy from run.
        """
        return self.ionic_steps[-1]["E0"]

    def as_dict(self):
        return {"electronic_steps": self.electronic_steps,
                "ionic_steps": self.ionic_steps}


class VaspParserError(Exception):
    """
    Exception class for VASP parsing.
    """
    pass


def get_band_structure_from_vasp_multiple_branches(dir_name, efermi=None,
                                                   projections=False):
    """
    This method is used to get band structure info from a VASP directory. It
    takes into account that the run can be divided in several branches named
    "branch_x". If the run has not been divided in branches the method will
    turn to parsing vasprun.xml directly.

    The method returns None is there"s a parsing error

    Args:
        dir_name: Directory containing all bandstructure runs.
        efermi: Efermi for bandstructure.
        projections: True if you want to get the data on site projections if
            any. Note that this is sometimes very large

    Returns:
        A BandStructure Object
    """
    #ToDo: Add better error handling!!!
    if os.path.exists(os.path.join(dir_name, "branch_0")):
        #get all branch dir names
        branch_dir_names = [os.path.abspath(d)
                            for d in glob.glob("{i}/branch_*"
                                               .format(i=dir_name))
                            if os.path.isdir(d)]

        #sort by the directory name (e.g, branch_10)
        sort_by = lambda x: int(x.split("_")[-1])
        sorted_branch_dir_names = sorted(branch_dir_names, key=sort_by)

        # populate branches with Bandstructure instances
        branches = []
        for dir_name in sorted_branch_dir_names:
            xml_file = os.path.join(dir_name, "vasprun.xml")
            if os.path.exists(xml_file):
                run = Vasprun(xml_file, parse_projected_eigen=projections)
                branches.append(run.get_band_structure(efermi=efermi))
            else:
                # It might be better to throw an exception
                warnings.warn("Skipping {}. Unable to find {}"
                              .format(d=dir_name, f=xml_file))

        return get_reconstructed_band_structure(branches, efermi)
    else:
        xml_file = os.path.join(dir_name, "vasprun.xml")
        #Better handling of Errors
        if os.path.exists(xml_file):
            return Vasprun(xml_file, parse_projected_eigen=projections)\
                .get_band_structure(kpoints_filename=None, efermi=efermi)
        else:
            return None


class Xdatcar(object):
    """
    Class representing an XDATCAR file. Only tested with VASP 5.x files.

    .. attribute:: structures

        List of structures parsed from XDATCAR.
    """

    def __init__(self, filename):
        """
        Init a Xdatcar.

        Args:
            filename (str): Filename of XDATCAR file.
        """
        preamble = None
        coords_str = []
        structures = []
        preamble_done = False
        with zopen(filename, "rt") as f:
            for l in f:
                l = l.strip()
                if preamble is None:
                    preamble = [l]
                elif not preamble_done:
                    if l == "" or "Direct configuration=" in l:
                        preamble_done = True
                    else:
                        preamble.append(l)
                elif l == "" or "Direct configuration=" in l:
                    p = Poscar.from_string("\n".join(preamble +
                                                     ["Direct"] + coords_str))
                    structures.append(p.structure)
                    coords_str = []
                else:
                    coords_str.append(l)
        self.structures = structures

class Dynmat(object):
    """
    Object for reading a DYNMAT file.

    Args:
        filename: Name of file containing DYNMAT.

    .. attribute:: data

        A nested dict containing the DYNMAT data of the form::
        [atom <int>][disp <int>]['dispvec'] =
            displacement vector (part of first line in dynmat block, e.g. "0.01 0 0")
        [atom <int>][disp <int>]['dynmat'] =
                <list> list of dynmat lines for this atom and this displacement

    Authors: Patrick Huck
    """
    def __init__(self, filename):
        with zopen(filename, "rt") as f:
            lines = list(clean_lines(f.readlines()))
            self._nspecs, self._natoms, self._ndisps = map(int, lines[0].split())
            self._masses = map(float, lines[1].split())
            self.data = defaultdict(dict)
            atom, disp = None, None
            for i,l in enumerate(lines[2:]):
                v = list(map(float, l.split()))
                if not i % (self._natoms+1):
                    atom, disp = map(int, v[:2])
                    if atom not in self.data: self.data[atom] = {}
                    if disp not in self.data[atom]: self.data[atom][disp] = {}
                    self.data[atom][disp]['dispvec'] = v[2:]
                else:
                    if 'dynmat' not in self.data[atom][disp]:
                        self.data[atom][disp]['dynmat'] = []
                    self.data[atom][disp]['dynmat'].append(v)

    def get_phonon_frequencies(self):
        """calculate phonon frequencies"""
        # TODO: the following is most likely not correct or suboptimal
        # hence for demonstration purposes only
        frequencies = []
        for k,v0 in self.data.iteritems():
            for v1 in v0.itervalues():
                vec = map(abs, v1['dynmat'][k-1])
                frequency = math.sqrt(sum(vec)) * 2.*math.pi*15.633302 # THz
                frequencies.append(frequency)
        return frequencies

    @property
    def nspecs(self):
        """returns the number of species"""
        return self._nspecs

    @property
    def natoms(self):
        """returns the number of atoms"""
        return self._natoms

    @property
    def ndisps(self):
        """returns the number of displacements"""
        return self._ndisps

    @property
    def masses(self):
        """returns the list of atomic masses"""
        return list(self._masses)

def get_adjusted_fermi_level(efermi, cbm, band_structure):
    """
    When running a band structure computations the fermi level needs to be
    take from the static run that gave the charge density used for the non-self
    consistent band structure run. Sometimes this fermi level is however a
    little too low because of the mismatch between the uniform grid used in
    the static run and the band structure k-points (e.g., the VBM is on Gamma
    and the Gamma point is not in the uniform mesh). Here we use a procedure
    consisting in looking for energy levels higher than the static fermi level
    (but lower than the LUMO) if any of these levels make the band structure
    appears insulating and not metallic anymore, we keep this adjusted fermi
    level. This procedure has shown to detect correctly most insulators.

    Args:
        efermi (float): Fermi energy of the static run
        cbm (float): Conduction band minimum of the static run
        run_bandstructure: a band_structure object

    Returns:
        a new adjusted fermi level
    """
    #make a working copy of band_structure
    bs_working = BandStructureSymmLine.from_dict(band_structure.as_dict())
    if bs_working.is_metal():
        e = efermi
        while e < cbm:
            e += 0.01
            bs_working._efermi = e
            if not bs_working.is_metal():
                return e
    return efermi

class Wavederf(object):
    """
    Object for reading a WAVEDERF file.

    Note: This file is only produced when LOPTICS is true AND vasp has been recompiled
    after uncommenting the line that calls WRT_CDER_BETWEEN_STATES_FORMATTED in
    linear_optics.F


    Args:
        filename: Name of file containing WAVEDERF.

    .. attribute:: data

        A numpy array containing the WAVEDERF data of the form below. It should
        be noted that VASP uses 1-based indexing for bands, but this is
        converted to 0-based numpy array indexing.

        For each kpoint (in the same order as in IBZKPT), and for each pair of bands:

            [ #kpoint index
             [ #band 1 index
              [ #band 2 index
               [cdum_x_real, cdum_x_imag, cdum_y_real, cdum_y_imag, cdum_z_real, cdum_z_imag]
              ]
             ]
            ]

        This structure follows the file format. Numpy array methods can be used to fetch data
        in a more useful way (e.g., get matrix elements between wo specific bands at each kpoint,
        fetch x/y/z components, real/imaginary parts, abs/phase, etc. )

    Author: Miguel Dias Costa
    """
    def __init__(self, filename):
        with zopen(filename, "rt") as f:
            header = f.readline().split()
            ispin = int(header[0])
            nb_kpoints = int(header[1])
            nb_bands = int(header[2])
            data = np.zeros((nb_kpoints,nb_bands,nb_bands,6))
            for ik in range(nb_kpoints):
                for ib1 in range(nb_bands):
                    for ib2 in range(nb_bands):
                        # each line in the file includes besides the band indexes, which are redundant,
                        # each band's energy and occupation, which are already available elsewhere,
                        # so we store only the 6 matrix elements after this 6 redundant values
                        data[ik][ib1][ib2] = [ float(element) for element in f.readline().split()[6:] ]

            self.data = data
            self._nb_kpoints = nb_kpoints
            self._nb_bands = nb_bands

    @property
    def nb_bands(self):
        """
        returns the number of bands in the band structure
        """
        return self._nb_bands

    @property
    def nb_kpoints(self):
        """
        Returns the number of k-points in the band structure calculation
        """
        return self._nb_kpoints

    def get_elements_between_bands(self, band_i, band_j):
        """
        Method returning a numpy array with elements

        [cdum_x_real, cdum_x_imag, cdum_y_real, cdum_y_imag, cdum_z_real, cdum_z_imag]

        between bands band_i and band_j (vasp 1-based indexing) for all kpoints.

        Args:
            band_i (Integer): Index of band i
            band_j (Integer): Index of band j

        Returns:
            a numpy list of elements for each kpoint
        """
        if band_i < 1 or band_i > self.nb_bands or band_j < 1 or band_j > self.nb_bands:
            raise ValueError("Band index out of bounds")

        return self.data[:,band_i-1,band_j-1,:] # using numpy array multidimensional slicing

