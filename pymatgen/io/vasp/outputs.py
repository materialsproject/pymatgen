# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Classes for reading/manipulating/writing VASP ouput files.
"""

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
import collections
from typing import Optional, Tuple, List

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from monty.io import zopen, reverse_readfile
from monty.json import MSONable
from monty.json import jsanitize
from monty.re import regrep
from monty.os.path import zpath
from monty.dev import deprecated

from pymatgen.core.composition import Composition
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.core.units import unitized
from pymatgen.electronic_structure.bandstructure import BandStructure, \
    BandStructureSymmLine, get_reconstructed_band_structure
from pymatgen.electronic_structure.core import Spin, Orbital, OrbitalType, Magmom
from pymatgen.electronic_structure.dos import CompleteDos, Dos
from pymatgen.entries.computed_entries import \
    ComputedEntry, ComputedStructureEntry
from pymatgen.io.wannier90 import Unk
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, Potcar
from pymatgen.util.io_utils import clean_lines, micro_pyawk
from pymatgen.util.num import make_symmetric_matrix_from_upper_tri


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
    if val_type == "int":
        return int(val)
    if val_type == "string":
        return val.strip()
    return float(val)


def _parse_v_parameters(val_type, val, filename, param_name):
    r"""
    Helper function to convert a Vasprun array-type parameter into the proper
    type. Boolean, int and float types are converted.

    Args:
        val_type: Value type parsed from vasprun.xml.
        val: Actual string value parsed for vasprun.xml.
        filename: Fullpath of vasprun.xml. Used for robust error handling.
            E.g., if vasprun.xml contains *** for some Incar parameters,
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
    if elem.get("type", None) == 'logical':
        m = [[i == 'T' for i in v.text.split()] for v in elem]
    else:
        m = [[_vasprun_float(i) for i in v.text.split()] for v in elem]
    return m


def _parse_from_incar(filename, key):
    """
    Helper function to parse a parameter from the INCAR.
    """
    dirname = os.path.dirname(filename)
    for f in os.listdir(dirname):
        if re.search(r"INCAR", f):
            warnings.warn("INCAR found. Using " + key + " from INCAR.")
            incar = Incar.from_file(os.path.join(dirname, f))
            if key in incar:
                return incar[key]
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


class Vasprun(MSONable):
    """
    Vastly improved cElementTree-based parser for vasprun.xml files. Uses
    iterparse to support incremental parsing of large files.
    Speedup over Dom is at least 2x for smallish files (~1Mb) to orders of
    magnitude for larger files (~10Mb).

    **Vasp results**

    .. attribute:: ionic_steps

        All ionic steps in the run as a list of
        {"structure": structure at end of run,
        "electronic_steps": {All electronic step data in vasprun file},
        "stresses": stress matrix}

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

        Final projected eigenvalues as a dict of {spin: nd-array}. To access
        a particular value, you need to do
        Vasprun.projected_eigenvalues[spin][kpoint index][band index][atom index][orbital_index]
        This representation is based on actual ordering in VASP and is meant as
        an intermediate representation to be converted into proper objects. The
        kpoint, band and atom indices are 0-based (unlike the 1-based indexing
        in VASP).


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

    .. attribute:: nionic_steps

        The total number of ionic steps. This number is always equal
        to the total number of steps in the actual run even if
        ionic_step_skip is used.

    .. attribute:: force_constants

        Force constants computed in phonon DFPT run(IBRION = 8).
        The data is a 4D numpy array of shape (natoms, natoms, 3, 3).

    .. attribute:: normalmode_eigenvals

        Normal mode frequencies.
        1D numpy array of size 3*natoms.

    .. attribute:: normalmode_eigenvecs

        Normal mode eigen vectors.
        3D numpy array of shape (3*natoms, natoms, 3).

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
                 parse_potcar_file=True, occu_tol=1e-8,
                 exception_on_bad_xml=True):
        """
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
            parse_potcar_file (bool/str): Whether to parse the potcar file to read
                the potcar hashes for the potcar_spec attribute. Defaults to True,
                where no hashes will be determined and the potcar_spec dictionaries
                will read {"symbol": ElSymbol, "hash": None}. By Default, looks in
                the same directory as the vasprun.xml, with same extensions as
                 Vasprun.xml. If a string is provided, looks at that filepath.
            occu_tol (float): Sets the minimum tol for the determination of the
                vbm and cbm. Usually the default of 1e-8 works well enough,
                but there may be pathological cases.
            exception_on_bad_xml (bool): Whether to throw a ParseException if a
                malformed XML is detected. Default to True, which ensures only
                proper vasprun.xml are parsed. You can set to False if you want
                partial results (e.g., if you are monitoring a calculation during a
                run), but use the results with care. A warning is issued.
        """
        self.filename = filename
        self.ionic_step_skip = ionic_step_skip
        self.ionic_step_offset = ionic_step_offset
        self.occu_tol = occu_tol
        self.exception_on_bad_xml = exception_on_bad_xml

        with zopen(filename, "rt") as f:
            if ionic_step_skip or ionic_step_offset:
                # remove parts of the xml file and parse the string
                run = f.read()
                steps = run.split("<calculation>")
                # The text before the first <calculation> is the preamble!
                preamble = steps.pop(0)
                self.nionic_steps = len(steps)
                new_steps = steps[ionic_step_offset::int(ionic_step_skip)]
                # add the tailing informat in the last step from the run
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
                self.update_charge_from_potcar(parse_potcar_file)

        if self.incar.get("ALGO", "") != "BSE" and (not self.converged):
            msg = "%s is an unconverged VASP run.\n" % filename
            msg += "Electronic convergence reached: %s.\n" % \
                   self.converged_electronic
            msg += "Ionic convergence reached: %s." % self.converged_ionic
            warnings.warn(msg, UnconvergedVASPWarning)

    def _parse(self, stream, parse_dos, parse_eigen, parse_projected_eigen):
        self.efermi = None
        self.eigenvalues = None
        self.projected_eigenvalues = None
        self.dielectric_data = {}
        self.other_dielectric = {}
        ionic_steps = []
        parsed_header = False
        try:
            for event, elem in ET.iterparse(stream):
                tag = elem.tag
                if not parsed_header:
                    if tag == "generator":
                        self.generator = self._parse_params(elem)
                    elif tag == "incar":
                        self.incar = self._parse_params(elem)
                    elif tag == "kpoints":
                        if not hasattr(self, 'kpoints'):
                            self.kpoints, self.actual_kpoints, self.actual_kpoints_weights = self._parse_kpoints(elem)
                    elif tag == "parameters":
                        self.parameters = self._parse_params(elem)
                    elif tag == "structure" and elem.attrib.get("name") == "initialpos":
                        self.initial_structure = self._parse_structure(elem)
                    elif tag == "atominfo":
                        self.atomic_symbols, self.potcar_symbols = self._parse_atominfo(elem)
                        self.potcar_spec = [{"titel": p,
                                             "hash": None} for
                                            p in self.potcar_symbols]
                if tag == "calculation":
                    parsed_header = True
                    if not self.parameters.get("LCHIMAG", False):
                        ionic_steps.append(self._parse_calculation(elem))
                    else:
                        ionic_steps.extend(self._parse_chemical_shielding_calculation(elem))
                elif parse_dos and tag == "dos":
                    try:
                        self.tdos, self.idos, self.pdos = self._parse_dos(elem)
                        self.efermi = self.tdos.efermi
                        self.dos_has_errors = False
                    except Exception:
                        self.dos_has_errors = True
                elif parse_eigen and tag == "eigenvalues":
                    self.eigenvalues = self._parse_eigen(elem)
                elif parse_projected_eigen and tag == "projected":
                    self.projected_eigenvalues = self._parse_projected_eigen(
                        elem)
                elif tag == "dielectricfunction":
                    if ("comment" not in elem.attrib or
                            elem.attrib["comment"] ==
                            "INVERSE MACROSCOPIC DIELECTRIC TENSOR (including "
                            "local field effects in RPA (Hartree))"):
                        if 'density' not in self.dielectric_data:
                            self.dielectric_data['density'] = self._parse_diel(
                                elem)
                        elif 'velocity' not in self.dielectric_data:
                            # "velocity-velocity" is also named
                            # "current-current" in OUTCAR
                            self.dielectric_data['velocity'] = self._parse_diel(
                                elem)
                        else:
                            raise NotImplementedError(
                                'This vasprun.xml has >2 unlabelled dielectric '
                                'functions')
                    else:
                        comment = elem.attrib["comment"]
                        # VASP 6+ has labels for the density and current
                        # derived dielectric constants
                        if comment == "density-density":
                            self.dielectric_data["density"] = self._parse_diel(
                                elem)
                        elif comment == "current-current":
                            self.dielectric_data["velocity"] = self._parse_diel(
                                elem)
                        else:
                            self.other_dielectric[comment] = self._parse_diel(
                                elem)

                elif tag == "varray" and elem.attrib.get("name") == 'opticaltransitions':
                    self.optical_transition = np.array(_parse_varray(elem))
                elif tag == "structure" and elem.attrib.get("name") == \
                        "finalpos":
                    self.final_structure = self._parse_structure(elem)
                elif tag == "dynmat":
                    hessian, eigenvalues, eigenvectors = self._parse_dynmat(elem)
                    natoms = len(self.atomic_symbols)
                    hessian = np.array(hessian)
                    self.force_constants = np.zeros((natoms, natoms, 3, 3), dtype='double')
                    for i in range(natoms):
                        for j in range(natoms):
                            self.force_constants[i, j] = hessian[i * 3:(i + 1) * 3, j * 3:(j + 1) * 3]
                    phonon_eigenvectors = []
                    for ev in eigenvectors:
                        phonon_eigenvectors.append(np.array(ev).reshape(natoms, 3))
                    self.normalmode_eigenvals = np.array(eigenvalues)
                    self.normalmode_eigenvecs = np.array(phonon_eigenvectors)
        except ET.ParseError as ex:
            if self.exception_on_bad_xml:
                raise ex
            warnings.warn(
                "XML is malformed. Parsing has stopped but partial data"
                "is available.", UserWarning)
        self.ionic_steps = ionic_steps
        self.vasp_version = self.generator["version"]

    @property
    def structures(self):
        """
        Returns:
             List of Structure objects for the structure at each ionic step.
        """
        return [step["structure"] for step in self.ionic_steps]

    @property
    def epsilon_static(self):
        """
        Property only available for DFPT calculations.

        Returns:
            The static part of the dielectric constant. Present when it's a DFPT run
            (LEPSILON=TRUE)
        """
        return self.ionic_steps[-1].get("epsilon", [])

    @property
    def epsilon_static_wolfe(self):
        """
        Property only available for DFPT calculations.

        Returns:
            The static part of the dielectric constant without any local field
            effects. Present when it's a DFPT run (LEPSILON=TRUE)
        """
        return self.ionic_steps[-1].get("epsilon_rpa", [])

    @property
    def epsilon_ionic(self):
        """
        Property only available for DFPT calculations and when IBRION=5, 6, 7 or 8.

        Returns:
            The ionic part of the static dielectric constant. Present when it's a
            DFPT run (LEPSILON=TRUE) and IBRION=5, 6, 7 or 8

        """
        return self.ionic_steps[-1].get("epsilon_ion", [])

    @property
    def dielectric(self):
        """
        Returns:
            The real and imaginary part of the dielectric constant (e.g., computed
            by RPA) in function of the energy (frequency). Optical properties (e.g.
            absorption coefficient) can be obtained through this.
            The data is given as a tuple of 3 values containing each of them
            the energy, the real part tensor, and the imaginary part tensor
            ([energies],[[real_partxx,real_partyy,real_partzz,real_partxy,
            real_partyz,real_partxz]],[[imag_partxx,imag_partyy,imag_partzz,
            imag_partxy, imag_partyz, imag_partxz]])
        """
        return self.dielectric_data['density']

    @property
    def optical_absorption_coeff(self):
        """
        Calculate the optical absorption coefficient
        from the dielectric constants. Note that this method is only
        implemented for optical properties calculated with GGA and BSE.
        Returns:
        optical absorption coefficient in list
        """
        if self.dielectric_data["density"]:
            real_avg = [sum(self.dielectric_data["density"][1][i][0:3]) / 3
                        for i in range(len(self.dielectric_data["density"][0]))]
            imag_avg = [sum(self.dielectric_data["density"][2][i][0:3]) / 3
                        for i in range(len(self.dielectric_data["density"][0]))]

            def f(freq, real, imag):
                """
                The optical absorption coefficient calculated in terms of
                equation
                """
                hbar = 6.582119514e-16  # eV/K
                coeff = np.sqrt(np.sqrt(real ** 2 + imag ** 2) - real) * np.sqrt(2) / hbar * freq
                return coeff

            absorption_coeff = [f(freq, real, imag) for freq, real, imag in
                                zip(self.dielectric_data["density"][0], real_avg, imag_avg)]
        return absorption_coeff

    @property
    def converged_electronic(self):
        """
        Returns:
            True if electronic step convergence has been reached in the final
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
        Returns:
            True if ionic step convergence has been reached, i.e. that vasp
            exited before reaching the max ionic steps for a relaxation run
        """
        nsw = self.parameters.get("NSW", 0)
        return nsw <= 1 or len(self.ionic_steps) < nsw

    @property
    def converged(self):
        """
        Returns:
            True if a relaxation run is converged both ionically and
            electronically.
        """
        return self.converged_electronic and self.converged_ionic

    @property  # type: ignore
    @unitized("eV")
    def final_energy(self):
        """
        Final energy from the vasp run.
        """
        try:
            final_istep = self.ionic_steps[-1]
            if final_istep["e_wo_entrp"] != final_istep['electronic_steps'][-1]["e_0_energy"]:
                warnings.warn("Final e_wo_entrp differs from the final "
                              "electronic step. VASP may have included some "
                              "corrections, e.g., vdw. Vasprun will return "
                              "the final e_wo_entrp, i.e., including "
                              "corrections in such instances.")
                return final_istep["e_wo_entrp"]
            return final_istep['electronic_steps'][-1]["e_0_energy"]
        except (IndexError, KeyError):
            warnings.warn("Calculation does not have a total energy. "
                          "Possibly a GW or similar kind of run. A value of "
                          "infinity is returned.")
            return float('inf')

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
        symbols = [re.split(r"_", s)[0] for s in symbols]
        if not self.incar.get("LDAU", False):
            return {}
        us = self.incar.get("LDAUU", self.parameters.get("LDAUU"))
        js = self.incar.get("LDAUJ", self.parameters.get("LDAUJ"))
        if len(js) != len(us):
            js = [0] * len(us)
        if len(us) == len(symbols):
            return {symbols[i]: us[i] - js[i] for i in range(len(symbols))}
        if sum(us) == 0 and sum(js) == 0:
            return {}
        raise VaspParserError("Length of U value parameters and atomic "
                              "symbols are mismatched")

    @property
    def run_type(self):
        """
        Returns the run type. Currently supports LDA, GGA, vdW-DF and HF calcs.

        TODO: Fix for other functional types like PW91, other vdW types, etc.
        """
        GGA_TYPES = {"RE": "revPBE", "PE": "PBE", "PS": "PBESol", "RP": "RevPBE+PADE", "AM": "AM05", "OR": "optPBE",
                     "BO": "optB88", "MK": "optB86b", "--": "GGA"}

        METAGGA_TYPES = {"TPSS": "TPSS", "RTPSS": "revTPSS", "M06L": "M06-L", "MBJ": "modified Becke-Johnson",
                         "SCAN": "SCAN", "MS0": "MadeSimple0", "MS1": "MadeSimple1", "MS2": "MadeSimple2"}

        if self.parameters.get("AEXX", 1.00) == 1.00:
            rt = "HF"
        elif self.parameters.get("HFSCREEN", 0.30) == 0.30:
            rt = "HSE03"
        elif self.parameters.get("HFSCREEN", 0.20) == 0.20:
            rt = "HSE06"
        elif self.parameters.get("AEXX", 0.20) == 0.20:
            rt = "B3LYP"
        elif self.parameters.get("LHFCALC", True):
            rt = "PBEO or other Hybrid Functional"
        elif self.parameters.get("LUSE_VDW", False):
            if self.incar.get("METAGGA", "").strip().upper() in METAGGA_TYPES:
                rt = METAGGA_TYPES[self.incar.get("METAGGA", "").strip().upper()] + "+rVV10"
            else:
                rt = GGA_TYPES[self.parameters.get("GGA", "").strip().upper()] + "+rVV10"
        elif self.incar.get("METAGGA", "").strip().upper() in METAGGA_TYPES:
            rt = METAGGA_TYPES[self.incar.get("METAGGA", "").strip().upper()]
            if self.is_hubbard or self.parameters.get("LDAU", True):
                rt += "+U"
        elif self.potcar_symbols[0].split()[0] == 'PAW':
            rt = "LDA"
        elif self.parameters.get("GGA", "").strip().upper() in GGA_TYPES:
            rt = GGA_TYPES[self.parameters.get("GGA", "").strip().upper()]
            if self.is_hubbard or self.parameters.get("LDAU", True):
                rt += "+U"
        return rt

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

    def get_computed_entry(self, inc_structure=True, parameters=None,
                           data=None):
        """
        Returns a ComputedStructureEntry from the vasprun.

        Args:
            inc_structure (bool): Set to True if you want
                ComputedStructureEntries to be returned instead of
                ComputedEntries.
            parameters (list): Input parameters to include. It has to be one of
                the properties supported by the Vasprun object. If
                parameters is None, a default set of parameters that are
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
            return ComputedStructureEntry(self.final_structure, self.final_energy, parameters=params, data=data)
        return ComputedEntry(self.final_structure.composition, self.final_energy, parameters=params, data=data)

    def get_band_structure(self, kpoints_filename=None, efermi=None,
                           line_mode=False, force_hybrid_mode=False):
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
            force_hybrid_mode (bool): Makes it possible to read in self-consistent band structure calculations for
                every type of functional

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
            kpoints_filename = zpath(
                os.path.join(os.path.dirname(self.filename), 'KPOINTS'))
        if not os.path.exists(kpoints_filename) and line_mode is True:
            raise VaspParserError('KPOINTS needed to obtain band structure '
                                  'along symmetry lines.')

        if efermi is None:
            efermi = self.efermi

        kpoint_file = None
        if os.path.exists(kpoints_filename):
            kpoint_file = Kpoints.from_file(kpoints_filename)
        lattice_new = Lattice(self.final_structure.lattice.reciprocal_lattice.matrix)

        kpoints = [np.array(self.actual_kpoints[i])
                   for i in range(len(self.actual_kpoints))]

        p_eigenvals = defaultdict(list)
        eigenvals = defaultdict(list)

        nkpts = len(kpoints)

        for spin, v in self.eigenvalues.items():
            v = np.swapaxes(v, 0, 1)
            eigenvals[spin] = v[:, :, 0]

            if self.projected_eigenvalues:
                peigen = self.projected_eigenvalues[spin]
                # Original axes for self.projected_eigenvalues are kpoints,
                # band, ion, orb.
                # For BS input, we need band, kpoints, orb, ion.
                peigen = np.swapaxes(peigen, 0, 1)  # Swap kpoint and band axes
                peigen = np.swapaxes(peigen, 2, 3)  # Swap ion and orb axes

                p_eigenvals[spin] = peigen
                # for b in range(min_eigenvalues):
                #     p_eigenvals[spin].append(
                #         [{Orbital(orb): v for orb, v in enumerate(peigen[b, k])}
                #          for k in range(nkpts)])

        # check if we have an hybrid band structure computation
        # for this we look at the presence of the LHFCALC tag
        hybrid_band = False
        if self.parameters.get('LHFCALC', False) or \
                0. in self.actual_kpoints_weights:
            hybrid_band = True

        if kpoint_file is not None:
            if kpoint_file.style == Kpoints.supported_modes.Line_mode:
                line_mode = True

        if line_mode:
            labels_dict = {}
            if hybrid_band or force_hybrid_mode:
                start_bs_index = 0
                for i in range(len(self.actual_kpoints)):
                    if self.actual_kpoints_weights[i] == 0.0:
                        start_bs_index = i
                        break
                for i in range(start_bs_index, len(kpoint_file.kpts)):
                    if kpoint_file.labels[i] is not None:
                        labels_dict[kpoint_file.labels[i]] = \
                            kpoint_file.kpts[i]
                # remake the data only considering line band structure k-points
                # (weight = 0.0 kpoints)
                nbands = len(eigenvals[Spin.up])
                kpoints = kpoints[start_bs_index:nkpts]
                up_eigen = [eigenvals[Spin.up][i][start_bs_index:nkpts]
                            for i in range(nbands)]
                if self.projected_eigenvalues:
                    p_eigenvals[Spin.up] = [p_eigenvals[Spin.up][i][
                                            start_bs_index:nkpts]
                                            for i in range(nbands)]
                if self.is_spin:
                    down_eigen = [eigenvals[Spin.down][i][start_bs_index:nkpts]
                                  for i in range(nbands)]
                    eigenvals = {Spin.up: up_eigen, Spin.down: down_eigen}
                    if self.projected_eigenvalues:
                        p_eigenvals[Spin.down] = [p_eigenvals[Spin.down][i][
                                                  start_bs_index:nkpts]
                                                  for i in range(nbands)]
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
        for spin, d in self.eigenvalues.items():
            for k, val in enumerate(d):
                for (eigenval, occu) in val:
                    if occu > self.occu_tol and eigenval > vbm:
                        vbm = eigenval
                        vbm_kpoint = k
                    elif occu <= self.occu_tol and eigenval < cbm:
                        cbm = eigenval
                        cbm_kpoint = k
        return max(cbm - vbm, 0), cbm, vbm, vbm_kpoint == cbm_kpoint

    def get_potcars(self, path):
        """
        :param path: Path to search for POTCARs
        :return: Potcar from path.
        """
        def get_potcar_in_path(p):
            for fn in os.listdir(os.path.abspath(p)):
                if fn.startswith('POTCAR') and ".spec" not in fn:
                    pc = Potcar.from_file(os.path.join(p, fn))
                    if {d.header for d in pc} == set(self.potcar_symbols):
                        return pc
            warnings.warn("No POTCAR file with matching TITEL fields"
                          " was found in {}".format(os.path.abspath(p)))
            return None

        if isinstance(path, (str, Path)):
            path = str(path)
            if "POTCAR" in path:
                potcar = Potcar.from_file(path)
                if {d.TITEL for d in potcar} != set(self.potcar_symbols):
                    raise ValueError("Potcar TITELs do not match Vasprun")
            else:
                potcar = get_potcar_in_path(path)
        elif isinstance(path, bool) and path:
            potcar = get_potcar_in_path(os.path.split(self.filename)[0])
        else:
            potcar = None

        return potcar

    def get_trajectory(self):
        """
        This method returns a Trajectory object, which is an alternative
        representation of self.structures into a single object. Forces are
        added to the Trajectory as site properties.

        Returns: a Trajectory
        """
        # required due to circular imports
        # TODO: fix pymatgen.core.trajectory so it does not load from io.vasp(!)
        from pymatgen.core.trajectory import Trajectory

        structs = []
        for step in self.ionic_steps:
            struct = step['structure'].copy()
            struct.add_site_property('forces', step['forces'])
            structs.append(struct)
        return Trajectory.from_structures(structs, constant_lattice=False)

    def update_potcar_spec(self, path):
        """
        :param path: Path to search for POTCARs
        :return: Potcar spec from path.
        """
        potcar = self.get_potcars(path)
        if potcar:
            self.potcar_spec = [{"titel": sym, "hash": ps.get_potcar_hash()}
                                for sym in self.potcar_symbols
                                for ps in potcar if
                                ps.symbol == sym.split()[1]]

    def update_charge_from_potcar(self, path):
        """
        Sets the charge of a structure based on the POTCARs found.

        :param path: Path to search for POTCARs
        """
        potcar = self.get_potcars(path)

        if potcar and self.incar.get("ALGO", "") not in ["GW0", "G0W0", "GW", "BSE"]:
            nelect = self.parameters["NELECT"]
            if len(potcar) == len(self.initial_structure.composition.element_composition):
                potcar_nelect = sum([
                    self.initial_structure.composition.element_composition[ps.element] * ps.ZVAL
                    for ps in potcar])
            else:
                nums = [len(list(g)) for _, g in
                        itertools.groupby(self.atomic_symbols)]
                potcar_nelect = sum(ps.ZVAL * num for ps, num in
                                    zip(potcar, nums))
            charge = nelect - potcar_nelect

            if charge:
                for s in self.structures:
                    s._charge = charge

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
        symbols = [re.split(r"_", s)[0] for s in symbols]
        d["is_hubbard"] = self.is_hubbard
        d["hubbards"] = self.hubbards

        unique_symbols = sorted(list(set(self.atomic_symbols)))
        d["elements"] = unique_symbols
        d["nelements"] = len(unique_symbols)

        d["run_type"] = self.run_type

        vin = {"incar": dict(self.incar.items()),
               "crystal": self.initial_structure.as_dict(),
               "kpoints": self.kpoints.as_dict()}
        actual_kpts = [{"abc": list(self.actual_kpoints[i]),
                        "weight": self.actual_kpoints_weights[i]}
                       for i in range(len(self.actual_kpoints))]
        vin["kpoints"]["actual_points"] = actual_kpts
        vin["nkpoints"] = len(actual_kpts)
        vin["potcar"] = [s.split(" ")[1] for s in self.potcar_symbols]
        vin["potcar_spec"] = self.potcar_spec
        vin["potcar_type"] = [s.split(" ")[0] for s in self.potcar_symbols]
        vin["parameters"] = dict(self.parameters.items())
        vin["lattice_rec"] = self.final_structure.lattice.reciprocal_lattice.as_dict()
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
            eigen = {str(spin): v.tolist()
                     for spin, v in self.eigenvalues.items()}
            vout["eigenvalues"] = eigen
            (gap, cbm, vbm, is_direct) = self.eigenvalue_band_properties
            vout.update(dict(bandgap=gap, cbm=cbm, vbm=vbm,
                             is_gap_direct=is_direct))

            if self.projected_eigenvalues:
                vout['projected_eigenvalues'] = {
                    str(spin): v.tolist()
                    for spin, v in self.projected_eigenvalues.items()}

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

    @staticmethod
    def _parse_atominfo(elem):
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
            except ValueError as e:
                if symbol == "X":
                    return "Xe"
                if symbol == "r":
                    return "Zr"
                raise e

        elem.clear()
        return [parse_atomic_symbol(sym) for sym in atomic_symbols], potcar_symbols

    @staticmethod
    def _parse_kpoints(elem):
        e = elem
        if elem.find("generation"):
            e = elem.find("generation")
        k = Kpoints("Kpoints from vasprun.xml")
        k.style = Kpoints.supported_modes.from_string(
            e.attrib["param"] if "param" in e.attrib else "Reciprocal")
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
        if k.style == Kpoints.supported_modes.Reciprocal:
            k = Kpoints(comment="Kpoints from vasprun.xml",
                        style=Kpoints.supported_modes.Reciprocal,
                        num_kpts=len(k.kpts),
                        kpts=actual_kpoints, kpts_weights=weights)
        return k, actual_kpoints, weights

    def _parse_structure(self, elem):
        latt = _parse_varray(elem.find("crystal").find("varray"))
        pos = _parse_varray(elem.find("varray"))
        struct = Structure(latt, self.atomic_symbols, pos)
        sdyn = elem.find("varray/[@name='selective']")
        if sdyn:
            struct.add_site_property('selective_dynamics',
                                     _parse_varray(sdyn))
        return struct

    @staticmethod
    def _parse_diel(elem):
        imag = [[_vasprun_float(l) for l in r.text.split()]
                for r in elem.find("imag").find("array").find("set").findall("r")]
        real = [[_vasprun_float(l) for l in r.text.split()]
                for r in elem.find("real").find("array").find("set").findall("r")]
        elem.clear()
        return [e[0] for e in imag], [e[1:] for e in real], [e[1:] for e in imag]

    @staticmethod
    def _parse_optical_transition(elem):
        for va in elem.findall("varray"):
            if va.attrib.get("name") == "opticaltransitions":
                # opticaltransitions array contains oscillator strength and probability of transition
                oscillator_strength = np.array(_parse_varray(va))[0:, ]
                probability_transition = np.array(_parse_varray(va))[0:, 1]
        return oscillator_strength, probability_transition

    def _parse_chemical_shielding_calculation(self, elem):
        calculation = []
        istep = {}
        try:
            s = self._parse_structure(elem.find("structure"))
        except AttributeError:  # not all calculations have a structure
            s = None
            pass
        for va in elem.findall("varray"):
            istep[va.attrib["name"]] = _parse_varray(va)
        istep["structure"] = s
        istep["electronic_steps"] = []
        calculation.append(istep)
        for scstep in elem.findall("scstep"):
            try:
                d = {i.attrib["name"]: _vasprun_float(i.text)
                     for i in scstep.find("energy").findall("i")}
                cur_ene = d['e_fr_energy']
                min_steps = 1 if len(calculation) >= 1 else self.parameters.get("NELMIN", 5)
                if len(calculation[-1]["electronic_steps"]) <= min_steps:
                    calculation[-1]["electronic_steps"].append(d)
                else:
                    last_ene = calculation[-1]["electronic_steps"][-1]["e_fr_energy"]
                    if abs(cur_ene - last_ene) < 1.0:
                        calculation[-1]["electronic_steps"].append(d)
                    else:
                        calculation.append({"electronic_steps": [d]})
            except AttributeError:  # not all calculations have an energy
                pass
        calculation[-1].update(calculation[-1]["electronic_steps"][-1])
        return calculation

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

    @staticmethod
    def _parse_dos(elem):
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
                            orb = Orbital(j - 1)
                        else:
                            orb = OrbitalType(j - 1)
                        pdos[orb][spin] = data[:, j]
                pdoss.append(pdos)
        elem.clear()
        return Dos(efermi, energies, tdensities), Dos(efermi, energies, idensities), pdoss

    @staticmethod
    def _parse_eigen(elem):
        eigenvalues = defaultdict(list)
        for s in elem.find("array").find("set").findall("set"):
            spin = Spin.up if s.attrib["comment"] == "spin 1" else Spin.down
            for ss in s.findall("set"):
                eigenvalues[spin].append(_parse_varray(ss))
        eigenvalues = {spin: np.array(v) for spin, v in eigenvalues.items()}
        elem.clear()
        return eigenvalues

    @staticmethod
    def _parse_projected_eigen(elem):
        root = elem.find("array").find("set")
        proj_eigen = defaultdict(list)
        for s in root.findall("set"):
            spin = int(re.match(r"spin(\d+)", s.attrib["comment"]).group(1))

            # Force spin to be +1 or -1
            spin = Spin.up if spin == 1 else Spin.down
            for kpt, ss in enumerate(s.findall("set")):
                dk = []
                for band, sss in enumerate(ss.findall("set")):
                    db = _parse_varray(sss)
                    dk.append(db)
                proj_eigen[spin].append(dk)
        proj_eigen = {spin: np.array(v) for spin, v in proj_eigen.items()}
        elem.clear()
        return proj_eigen

    @staticmethod
    def _parse_dynmat(elem):
        hessian = []
        eigenvalues = []
        eigenvectors = []
        for v in elem.findall("v"):
            if v.attrib["name"] == "eigenvalues":
                eigenvalues = [float(i) for i in v.text.split()]
        for va in elem.findall("varray"):
            if va.attrib["name"] == "hessian":
                for v in va.findall("v"):
                    hessian.append([float(i) for i in v.text.split()])
            elif va.attrib["name"] == "eigenvectors":
                for v in va.findall("v"):
                    eigenvectors.append([float(i) for i in v.text.split()])
        return hessian, eigenvalues, eigenvectors


class BSVasprun(Vasprun):
    """
    A highly optimized version of Vasprun that parses only eigenvalues for
    bandstructures. All other properties like structures, parameters,
    etc. are ignored.
    """

    def __init__(self, filename, parse_projected_eigen=False,
                 parse_potcar_file=False, occu_tol=1e-8):
        """
        Args:
            filename (str): Filename to parse
            parse_projected_eigen (bool): Whether to parse the projected
                eigenvalues. Defaults to False. Set to True to obtain projected
                eigenvalues. **Note that this can take an extreme amount of time
                and memory.** So use this wisely.
            parse_potcar_file (bool/str): Whether to parse the potcar file to read
                the potcar hashes for the potcar_spec attribute. Defaults to True,
                where no hashes will be determined and the potcar_spec dictionaries
                will read {"symbol": ElSymbol, "hash": None}. By Default, looks in
                the same directory as the vasprun.xml, with same extensions as
                 Vasprun.xml. If a string is provided, looks at that filepath.
            occu_tol (float): Sets the minimum tol for the determination of the
                vbm and cbm. Usually the default of 1e-8 works well enough,
                but there may be pathological cases.
        """
        self.filename = filename
        self.occu_tol = occu_tol

        with zopen(filename, "rt") as f:
            self.efermi = None
            parsed_header = False
            self.eigenvalues = None
            self.projected_eigenvalues = None
            for event, elem in ET.iterparse(f):
                tag = elem.tag
                if not parsed_header:
                    if tag == "generator":
                        self.generator = self._parse_params(elem)
                    elif tag == "incar":
                        self.incar = self._parse_params(elem)
                    elif tag == "kpoints":
                        self.kpoints, self.actual_kpoints, self.actual_kpoints_weights = self._parse_kpoints(elem)
                    elif tag == "parameters":
                        self.parameters = self._parse_params(elem)
                    elif tag == "atominfo":
                        self.atomic_symbols, self.potcar_symbols = self._parse_atominfo(elem)
                        self.potcar_spec = [{"titel": p, "hash": None} for p in self.potcar_symbols]
                        parsed_header = True
                elif tag == "i" and elem.attrib.get("name") == "efermi":
                    self.efermi = float(elem.text)
                elif tag == "eigenvalues":
                    self.eigenvalues = self._parse_eigen(elem)
                elif parse_projected_eigen and tag == "projected":
                    self.projected_eigenvalues = self._parse_projected_eigen(
                        elem)
                elif tag == "structure" and elem.attrib.get("name") == \
                        "finalpos":
                    self.final_structure = self._parse_structure(elem)
        self.vasp_version = self.generator["version"]
        if parse_potcar_file:
            self.update_potcar_spec(parse_potcar_file)

    def as_dict(self):
        """
        Json-serializable dict representation.
        """
        d = {"vasp_version": self.vasp_version,
             "has_vasp_completed": True,
             "nsites": len(self.final_structure)}
        comp = self.final_structure.composition
        d["unit_cell_formula"] = comp.as_dict()
        d["reduced_cell_formula"] = Composition(comp.reduced_formula).as_dict()
        d["pretty_formula"] = comp.reduced_formula
        d["is_hubbard"] = self.is_hubbard
        d["hubbards"] = self.hubbards

        unique_symbols = sorted(list(set(self.atomic_symbols)))
        d["elements"] = unique_symbols
        d["nelements"] = len(unique_symbols)

        d["run_type"] = self.run_type

        vin = {"incar": dict(self.incar),
               "crystal": self.final_structure.as_dict(),
               "kpoints": self.kpoints.as_dict()}
        actual_kpts = [{"abc": list(self.actual_kpoints[i]),
                        "weight": self.actual_kpoints_weights[i]}
                       for i in range(len(self.actual_kpoints))]
        vin["kpoints"]["actual_points"] = actual_kpts
        vin["potcar"] = [s.split(" ")[1] for s in self.potcar_symbols]
        vin["potcar_spec"] = self.potcar_spec
        vin["potcar_type"] = [s.split(" ")[0] for s in self.potcar_symbols]
        vin["parameters"] = dict(self.parameters)
        vin["lattice_rec"] = self.final_structure.lattice.reciprocal_lattice.as_dict()
        d["input"] = vin

        vout = {"crystal": self.final_structure.as_dict(),
                "efermi": self.efermi}

        if self.eigenvalues:
            eigen = defaultdict(dict)
            for spin, values in self.eigenvalues.items():
                for i, v in enumerate(values):
                    eigen[i][str(spin)] = v
            vout["eigenvalues"] = eigen
            (gap, cbm, vbm, is_direct) = self.eigenvalue_band_properties
            vout.update(dict(bandgap=gap, cbm=cbm, vbm=vbm,
                             is_gap_direct=is_direct))

            if self.projected_eigenvalues:
                peigen = []
                for i in range(len(eigen)):
                    peigen.append({})
                for spin, v in self.projected_eigenvalues.items():
                    for kpoint_index, vv in enumerate(v):
                        if str(spin) not in peigen[kpoint_index]:
                            peigen[kpoint_index][str(spin)] = vv
                vout['projected_eigenvalues'] = peigen

        d['output'] = vout
        return jsanitize(d, strict=True)


class Outcar:
    """
    Parser for data in OUTCAR that is not available in Vasprun.xml

    Note, this class works a bit differently than most of the other
    VaspObjects, since the OUTCAR can be very different depending on which
    "type of run" performed.

    Creating the OUTCAR class with a filename reads "regular parameters" that
    are always present.

    .. attribute:: magnetization

        Magnetization on each ion as a tuple of dict, e.g.,
        ({"d": 0.0, "p": 0.003, "s": 0.002, "tot": 0.005}, ... )
        Note that this data is not always present.  LORBIT must be set to some
        other value than the default.

    .. attribute:: chemical_shielding

        chemical shielding on each ion as a dictionary with core and valence contributions

    .. attribute:: unsym_cs_tensor

        Unsymmetrized chemical shielding tensor matrixes on each ion as a list.
        e.g.,
        [[[sigma11, sigma12, sigma13],
          [sigma21, sigma22, sigma23],
          [sigma31, sigma32, sigma33]],
          ...
         [[sigma11, sigma12, sigma13],
          [sigma21, sigma22, sigma23],
          [sigma31, sigma32, sigma33]]]

    .. attribute:: cs_g0_contribution

        G=0 contribution to chemical shielding. 2D rank 3 matrix

    .. attribute:: cs_core_contribution

        Core contribution to chemical shielding. dict. e.g.,
        {'Mg': -412.8, 'C': -200.5, 'O': -271.1}

    .. attribute:: efg

        Electric Field Gradient (EFG) tensor on each ion as a tuple of dict, e.g.,
        ({"cq": 0.1, "eta", 0.2, "nuclear_quadrupole_moment": 0.3},
         {"cq": 0.7, "eta", 0.8, "nuclear_quadrupole_moment": 0.9},
         ...)

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

    .. attribute:: drift
        Total drift for each step in eV/Atom

    .. attribute:: ngf
        Dimensions for the Augementation grid

    .. attribute: sampling_radii
        Size of the sampling radii in VASP for the test charges for
        the electrostatic potential at each atom. Total array size is the number
        of elements present in the calculation

    .. attribute: electrostatic_potential
        Average electrostatic potential at each atomic position in order
        of the atoms in POSCAR.

    ..attribute: final_energy_contribs
        Individual contributions to the total final energy as a dictionary.
        Include contirbutions from keys, e.g.:
        {'DENC': -505778.5184347, 'EATOM': 15561.06492564, 'EBANDS': -804.53201231,
        'EENTRO': -0.08932659, 'EXHF': 0.0, 'Ediel_sol': 0.0,
        'PAW double counting': 664.6726974100002, 'PSCENC': 742.48691646,
        'TEWEN': 489742.86847338, 'XCENC': -169.64189814}

    One can then call a specific reader depending on the type of run being
    performed. These are currently: read_igpar(), read_lepsilon() and
    read_lcalcpol(), read_core_state_eign(), read_avg_core_pot().

    See the documentation of those methods for more documentation.

    Authors: Rickard Armiento, Shyue Ping Ong
    """

    def __init__(self, filename):
        """
        Args:
            filename (str): OUTCAR filename to parse.
        """
        self.filename = filename
        self.is_stopped = False

        # data from end of OUTCAR
        charge = []
        mag_x = []
        mag_y = []
        mag_z = []
        header = []
        run_stats = {}
        total_mag = None
        nelect = None
        efermi = None
        total_energy = None

        time_patt = re.compile(r"\((sec|kb)\)")
        efermi_patt = re.compile(r"E-fermi\s*:\s*(\S+)")
        nelect_patt = re.compile(r"number of electron\s+(\S+)\s+magnetization")
        mag_patt = re.compile(r"number of electron\s+\S+\s+magnetization\s+("
                              r"\S+)")
        toten_pattern = re.compile(r"free  energy   TOTEN\s+=\s+([\d\-\.]+)")

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
                        # try-catch because VASP sometimes prints
                        # 'E-fermi: ********     XC(G=0):  -6.1327
                        # alpha+bet : -1.8238'
                        efermi = float(m.group(1))
                        continue
                    except ValueError:
                        efermi = None
                        continue
                m = nelect_patt.search(clean)
                if m:
                    nelect = float(m.group(1))
                m = mag_patt.search(clean)
                if m:
                    total_mag = float(m.group(1))
                if total_energy is None:
                    m = toten_pattern.search(clean)
                    if m:
                        total_energy = float(m.group(1))
            if all([nelect, total_mag is not None, efermi is not None,
                    run_stats]):
                break

        # For single atom systems, VASP doesn't print a total line, so
        # reverse parsing is very difficult
        read_charge = False
        read_mag_x = False
        read_mag_y = False  # for SOC calculations only
        read_mag_z = False
        all_lines.reverse()
        for clean in all_lines:
            if read_charge or read_mag_x or read_mag_y or read_mag_z:
                if clean.startswith("# of ion"):
                    header = re.split(r"\s{2,}", clean.strip())
                    header.pop(0)
                else:
                    m = re.match(r"\s*(\d+)\s+(([\d\.\-]+)\s+)+", clean)
                    if m:
                        toks = [float(i)
                                for i in re.findall(r"[\d\.\-]+", clean)]
                        toks.pop(0)
                        if read_charge:
                            charge.append(dict(zip(header, toks)))
                        elif read_mag_x:
                            mag_x.append(dict(zip(header, toks)))
                        elif read_mag_y:
                            mag_y.append(dict(zip(header, toks)))
                        elif read_mag_z:
                            mag_z.append(dict(zip(header, toks)))
                    elif clean.startswith('tot'):
                        read_charge = False
                        read_mag_x = False
                        read_mag_y = False
                        read_mag_z = False
            if clean == "total charge":
                charge = []
                read_charge = True
                read_mag_x, read_mag_y, read_mag_z = False, False, False
            elif clean == "magnetization (x)":
                mag_x = []
                read_mag_x = True
                read_charge, read_mag_y, read_mag_z = False, False, False
            elif clean == "magnetization (y)":
                mag_y = []
                read_mag_y = True
                read_charge, read_mag_x, read_mag_z = False, False, False
            elif clean == "magnetization (z)":
                mag_z = []
                read_mag_z = True
                read_charge, read_mag_x, read_mag_y = False, False, False
            elif re.search("electrostatic", clean):
                read_charge, read_mag_x, read_mag_y, read_mag_z = False, False, False, False

        # merge x, y and z components of magmoms if present (SOC calculation)
        if mag_y and mag_z:
            # TODO: detect spin axis
            mag = []
            for idx in range(len(mag_x)):
                mag.append({
                    key: Magmom([mag_x[idx][key], mag_y[idx][key], mag_z[idx][key]])
                    for key in mag_x[0].keys()
                })
        else:
            mag = mag_x

        # data from beginning of OUTCAR
        run_stats['cores'] = 0
        with zopen(filename, "rt") as f:
            for line in f:
                if "running" in line:
                    run_stats['cores'] = line.split()[2]
                    break

        self.run_stats = run_stats
        self.magnetization = tuple(mag)
        self.charge = tuple(charge)
        self.efermi = efermi
        self.nelect = nelect
        self.total_mag = total_mag
        self.final_energy = total_energy
        self.data = {}

        # Read "total number of plane waves", NPLWV:
        self.read_pattern(
            {"nplwv": r"total plane-waves  NPLWV =\s+(\*{6}|\d+)"},
            terminate_on_match=True
        )
        try:
            self.data["nplwv"] = [[int(self.data["nplwv"][0][0])]]
        except ValueError:
            self.data["nplwv"] = [[None]]

        nplwvs_at_kpoints = [
            n for [n] in self.read_table_pattern(
                r"\n{3}-{104}\n{3}",
                r".+plane waves:\s+(\*{6,}|\d+)",
                r"maximum and minimum number of plane-waves"
            )
        ]
        self.data["nplwvs_at_kpoints"] = [None for n in nplwvs_at_kpoints]
        for (n, nplwv) in enumerate(nplwvs_at_kpoints):
            try:
                self.data["nplwvs_at_kpoints"][n] = int(nplwv)
            except ValueError:
                pass

        # Read the drift:
        self.read_pattern({
            "drift": r"total drift:\s+([\.\-\d]+)\s+([\.\-\d]+)\s+([\.\-\d]+)"},
            terminate_on_match=False,
            postprocess=float)
        self.drift = self.data.get('drift', [])

        # Check if calculation is spin polarized
        self.spin = False
        self.read_pattern({'spin': 'ISPIN  =      2'})
        if self.data.get('spin', []):
            self.spin = True

        # Check if calculation is noncollinear
        self.noncollinear = False
        self.read_pattern({'noncollinear': 'LNONCOLLINEAR =      T'})
        if self.data.get('noncollinear', []):
            self.noncollinear = False

        # Check if the calculation type is DFPT
        self.dfpt = False
        self.read_pattern({'ibrion': r"IBRION =\s+([\-\d]+)"},
                          terminate_on_match=True,
                          postprocess=int)
        if self.data.get("ibrion", [[0]])[0][0] > 6:
            self.dfpt = True
            self.read_internal_strain_tensor()

        # Check to see if LEPSILON is true and read piezo data if so
        self.lepsilon = False
        self.read_pattern({'epsilon': 'LEPSILON=     T'})
        if self.data.get('epsilon', []):
            self.lepsilon = True
            self.read_lepsilon()
            # only read ionic contribution if DFPT is turned on
            if self.dfpt:
                self.read_lepsilon_ionic()

        # Check to see if LCALCPOL is true and read polarization data if so
        self.lcalcpol = False
        self.read_pattern({'calcpol': 'LCALCPOL   =     T'})
        if self.data.get('calcpol', []):
            self.lcalcpol = True
            self.read_lcalcpol()
            self.read_pseudo_zval()

        # Read electrostatic potential
        self.read_pattern({
            'electrostatic': r"average \(electrostatic\) potential at core"})
        if self.data.get('electrostatic', []):
            self.read_electrostatic_potential()

        self.nmr_cs = False
        self.read_pattern({"nmr_cs": r"LCHIMAG   =     (T)"})
        if self.data.get("nmr_cs", None):
            self.nmr_cs = True
            self.read_chemical_shielding()
            self.read_cs_g0_contribution()
            self.read_cs_core_contribution()
            self.read_cs_raw_symmetrized_tensors()

        self.nmr_efg = False
        self.read_pattern({"nmr_efg": r"NMR quadrupolar parameters"})
        if self.data.get("nmr_efg", None):
            self.nmr_efg = True
            self.read_nmr_efg()
            self.read_nmr_efg_tensor()

        self.has_onsite_density_matrices = False
        self.read_pattern({"has_onsite_density_matrices": r"onsite density matrix"},
                          terminate_on_match=True)
        if "has_onsite_density_matrices" in self.data:
            self.has_onsite_density_matrices = True
            self.read_onsite_density_matrices()

        # Store the individual contributions to the final total energy
        final_energy_contribs = {}
        for k in ["PSCENC", "TEWEN", "DENC", "EXHF", "XCENC", "PAW double counting",
                  "EENTRO", "EBANDS", "EATOM", "Ediel_sol"]:
            if k == "PAW double counting":
                self.read_pattern({k: r"%s\s+=\s+([\.\-\d]+)\s+([\.\-\d]+)" % (k)})
            else:
                self.read_pattern({k: r"%s\s+=\s+([\d\-\.]+)" % (k)})
            if not self.data[k]:
                continue
            final_energy_contribs[k] = sum([float(f) for f in self.data[k][-1]])
        self.final_energy_contribs = final_energy_contribs

    def read_pattern(self, patterns, reverse=False, terminate_on_match=False,
                     postprocess=str):
        r"""
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
        r"""
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
        table_pattern_text = header_pattern + r"\s*^(?P<table_body>(?:\s+" + row_pattern + r")+)\s+" + footer_pattern
        table_pattern = re.compile(table_pattern_text, re.MULTILINE | re.DOTALL)
        rp = re.compile(row_pattern)
        tables = []
        for mt in table_pattern.finditer(text):
            table_body_text = mt.group("table_body")
            table_contents = []
            for line in table_body_text.split("\n"):
                ml = rp.search(line)
                # skip empty lines
                if not ml:
                    continue
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

    def read_electrostatic_potential(self):
        """
        Parses the eletrostatic potential for the last ionic step
        """
        pattern = {"ngf": r"\s+dimension x,y,z NGXF=\s+([\.\-\d]+)\sNGYF=\s+([\.\-\d]+)\sNGZF=\s+([\.\-\d]+)"}
        self.read_pattern(pattern, postprocess=int)
        self.ngf = self.data.get("ngf", [[]])[0]

        pattern = {"radii": r"the test charge radii are((?:\s+[\.\-\d]+)+)"}
        self.read_pattern(pattern, reverse=True, terminate_on_match=True, postprocess=str)
        self.sampling_radii = [float(f) for f in self.data["radii"][0][0].split()]

        header_pattern = r"\(the norm of the test charge is\s+[\.\-\d]+\)"
        table_pattern = r"((?:\s+\d+\s*[\.\-\d]+)+)"
        footer_pattern = r"\s+E-fermi :"

        pots = self.read_table_pattern(header_pattern, table_pattern, footer_pattern)
        pots = "".join(itertools.chain.from_iterable(pots))

        pots = re.findall(r"\s+\d+\s*([\.\-\d]+)+", pots)

        self.electrostatic_potential = [float(f) for f in pots]

    @staticmethod
    def _parse_sci_notation(line):
        """
        Method to parse lines with values in scientific notation and potentially
        without spaces in between the values. This assumes that the scientific
        notation always lists two digits for the exponent, e.g. 3.535E-02
        Args:
            line: line to parse

        Returns: an array of numbers if found, or empty array if not

        """
        m = re.findall(r"[\.\-\d]+E[\+\-]\d{2}", line)
        if m:
            return [float(t) for t in m]
        return []

    def read_freq_dielectric(self):
        """
        Parses the frequency dependent dielectric function (obtained with
        LOPTICS). Frequencies (in eV) are in self.frequencies, and dielectric
        tensor function is given as self.dielectric_tensor_function.
        """

        plasma_pattern = r"plasma frequency squared.*"
        dielectric_pattern = r"frequency dependent\s+IMAGINARY " \
                             r"DIELECTRIC FUNCTION \(independent particle, " \
                             r"no local field effects\)(\sdensity-density)*$"
        row_pattern = r"\s+".join([r"([\.\-\d]+)"] * 3)
        plasma_frequencies = collections.defaultdict(list)
        read_plasma = False
        read_dielectric = False
        energies = []
        data = {"REAL": [], "IMAGINARY": []}
        count = 0
        component = "IMAGINARY"
        with zopen(self.filename, "rt") as f:
            for l in f:
                l = l.strip()
                if re.match(plasma_pattern, l):
                    read_plasma = "intraband" if "intraband" in l else "interband"
                elif re.match(dielectric_pattern, l):
                    read_plasma = False
                    read_dielectric = True
                    row_pattern = r"\s+".join([r"([\.\-\d]+)"] * 7)

                if read_plasma and re.match(row_pattern, l):
                    plasma_frequencies[read_plasma].append(
                        [float(t) for t in l.strip().split()])
                elif read_plasma and Outcar._parse_sci_notation(l):
                    plasma_frequencies[read_plasma].append(
                        Outcar._parse_sci_notation(l))
                elif read_dielectric:
                    toks = None
                    if re.match(row_pattern, l.strip()):
                        toks = l.strip().split()
                    elif Outcar._parse_sci_notation(l.strip()):
                        toks = Outcar._parse_sci_notation(l.strip())
                    elif re.match(r"\s*-+\s*", l):
                        count += 1

                    if toks:
                        if component == "IMAGINARY":
                            energies.append(float(toks[0]))
                        xx, yy, zz, xy, yz, xz = [float(t) for t in toks[1:]]
                        matrix = [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]
                        data[component].append(matrix)

                    if count == 2:
                        component = "REAL"
                    elif count == 3:
                        break

        self.plasma_frequencies = {k: np.array(v[:3])
                                   for k, v in plasma_frequencies.items()}
        self.dielectric_energies = np.array(energies)
        self.dielectric_tensor_function = np.array(data["REAL"]) + 1j * np.array(data["IMAGINARY"])

    @property  # type: ignore
    @deprecated(message="frequencies has been renamed to dielectric_energies.")
    def frequencies(self):
        """
        Renamed to dielectric energies.
        """
        return self.dielectric_energies

    def read_chemical_shielding(self):
        """
        Parse the NMR chemical shieldings data. Only the second part "absolute, valence and core"
        will be parsed. And only the three right most field (ISO_SHIELDING, SPAN, SKEW) will be retrieved.

        Returns:
            List of chemical shieldings in the order of atoms from the OUTCAR. Maryland notation is adopted.
        """
        header_pattern = r"\s+CSA tensor \(J\. Mason, Solid State Nucl\. Magn\. Reson\. 2, " \
                         r"285 \(1993\)\)\s+" \
                         r"\s+-{50,}\s+" \
                         r"\s+EXCLUDING G=0 CONTRIBUTION\s+INCLUDING G=0 CONTRIBUTION\s+" \
                         r"\s+-{20,}\s+-{20,}\s+" \
                         r"\s+ATOM\s+ISO_SHIFT\s+SPAN\s+SKEW\s+ISO_SHIFT\s+SPAN\s+SKEW\s+" \
                         r"-{50,}\s*$"
        first_part_pattern = r"\s+\(absolute, valence only\)\s+$"
        swallon_valence_body_pattern = r".+?\(absolute, valence and core\)\s+$"
        row_pattern = r"\d+(?:\s+[-]?\d+\.\d+){3}\s+" + r'\s+'.join(
            [r"([-]?\d+\.\d+)"] * 3)
        footer_pattern = r"-{50,}\s*$"
        h1 = header_pattern + first_part_pattern
        cs_valence_only = self.read_table_pattern(
            h1, row_pattern, footer_pattern, postprocess=float,
            last_one_only=True)
        h2 = header_pattern + swallon_valence_body_pattern
        cs_valence_and_core = self.read_table_pattern(
            h2, row_pattern, footer_pattern, postprocess=float,
            last_one_only=True)
        all_cs = {}
        for name, cs_table in [["valence_only", cs_valence_only],
                               ["valence_and_core", cs_valence_and_core]]:
            all_cs[name] = cs_table
        self.data["chemical_shielding"] = all_cs

    def read_cs_g0_contribution(self):
        """
        Parse the  G0 contribution of NMR chemical shielding.

        Returns:
            G0 contribution matrix as list of list.
        """
        header_pattern = r'^\s+G\=0 CONTRIBUTION TO CHEMICAL SHIFT \(field along BDIR\)\s+$\n' \
                         r'^\s+-{50,}$\n' \
                         r'^\s+BDIR\s+X\s+Y\s+Z\s*$\n' \
                         r'^\s+-{50,}\s*$\n'
        row_pattern = r'(?:\d+)\s+' + r'\s+'.join([r'([-]?\d+\.\d+)'] * 3)
        footer_pattern = r'\s+-{50,}\s*$'
        self.read_table_pattern(header_pattern, row_pattern, footer_pattern, postprocess=float,
                                last_one_only=True, attribute_name="cs_g0_contribution")

    def read_cs_core_contribution(self):
        """
            Parse the core contribution of NMR chemical shielding.

            Returns:
            G0 contribution matrix as list of list.
        """
        header_pattern = r'^\s+Core NMR properties\s*$\n' \
                         r'\n' \
                         r'^\s+typ\s+El\s+Core shift \(ppm\)\s*$\n' \
                         r'^\s+-{20,}$\n'
        row_pattern = r'\d+\s+(?P<element>[A-Z][a-z]?\w?)\s+(?P<shift>[-]?\d+\.\d+)'
        footer_pattern = r'\s+-{20,}\s*$'
        self.read_table_pattern(header_pattern, row_pattern, footer_pattern, postprocess=str,
                                last_one_only=True, attribute_name="cs_core_contribution")
        core_contrib = {d['element']: float(d['shift'])
                        for d in self.data["cs_core_contribution"]}
        self.data["cs_core_contribution"] = core_contrib

    def read_cs_raw_symmetrized_tensors(self):
        """
        Parse the matrix form of NMR tensor before corrected to table.

        Returns:
            nsymmetrized tensors list in the order of atoms.
        """
        header_pattern = r"\s+-{50,}\s+" \
                         r"\s+Absolute Chemical Shift tensors\s+" \
                         r"\s+-{50,}$"
        first_part_pattern = r"\s+UNSYMMETRIZED TENSORS\s+$"
        row_pattern = r"\s+".join([r"([-]?\d+\.\d+)"] * 3)
        unsym_footer_pattern = r"^\s+SYMMETRIZED TENSORS\s+$"

        with zopen(self.filename, 'rt') as f:
            text = f.read()
        unsym_table_pattern_text = header_pattern + first_part_pattern + r"(?P<table_body>.+)" + unsym_footer_pattern
        table_pattern = re.compile(unsym_table_pattern_text, re.MULTILINE | re.DOTALL)
        rp = re.compile(row_pattern)
        m = table_pattern.search(text)
        if m:
            table_text = m.group("table_body")
            micro_header_pattern = r"ion\s+\d+"
            micro_table_pattern_text = micro_header_pattern + r"\s*^(?P<table_body>(?:\s*" + row_pattern + r")+)\s+"
            micro_table_pattern = re.compile(micro_table_pattern_text,
                                             re.MULTILINE | re.DOTALL)
            unsym_tensors = []
            for mt in micro_table_pattern.finditer(table_text):
                table_body_text = mt.group("table_body")
                tensor_matrix = []
                for line in table_body_text.rstrip().split("\n"):
                    ml = rp.search(line)
                    processed_line = [float(v) for v in ml.groups()]
                    tensor_matrix.append(processed_line)
                unsym_tensors.append(tensor_matrix)
            self.data["unsym_cs_tensor"] = unsym_tensors
        else:
            raise ValueError("NMR UNSYMMETRIZED TENSORS is not found")

    def read_nmr_efg_tensor(self):
        """
        Parses the NMR Electric Field Gradient Raw Tensors

        Returns:
            A list of Electric Field Gradient Tensors in the order of Atoms from OUTCAR
        """

        header_pattern = r'Electric field gradients \(V/A\^2\)\n' \
                         r'-*\n' \
                         r' ion\s+V_xx\s+V_yy\s+V_zz\s+V_xy\s+V_xz\s+V_yz\n' \
                         r'-*\n'

        row_pattern = r'\d+\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)'
        footer_pattern = r'-*\n'

        data = self.read_table_pattern(header_pattern, row_pattern, footer_pattern, postprocess=float)
        tensors = [make_symmetric_matrix_from_upper_tri(d) for d in data]
        self.data["unsym_efg_tensor"] = tensors
        return tensors

    def read_nmr_efg(self):
        """
        Parse the NMR Electric Field Gradient interpretted values.

        Returns:
            Electric Field Gradient tensors as a list of dict in the order of atoms from OUTCAR.
            Each dict key/value pair corresponds to a component of the tensors.
        """
        header_pattern = r'^\s+NMR quadrupolar parameters\s+$\n' \
                         r'^\s+Cq : quadrupolar parameter\s+Cq=e[*]Q[*]V_zz/h$\n' \
                         r'^\s+eta: asymmetry parameters\s+\(V_yy - V_xx\)/ V_zz$\n' \
                         r'^\s+Q  : nuclear electric quadrupole moment in mb \(millibarn\)$\n' \
                         r'^-{50,}$\n' \
                         r'^\s+ion\s+Cq\(MHz\)\s+eta\s+Q \(mb\)\s+$\n' \
                         r'^-{50,}\s*$\n'
        row_pattern = r'\d+\s+(?P<cq>[-]?\d+\.\d+)\s+(?P<eta>[-]?\d+\.\d+)\s+' \
                      r'(?P<nuclear_quadrupole_moment>[-]?\d+\.\d+)'
        footer_pattern = r'-{50,}\s*$'
        self.read_table_pattern(header_pattern, row_pattern, footer_pattern, postprocess=float,
                                last_one_only=True, attribute_name="efg")

    def read_elastic_tensor(self):
        """
        Parse the elastic tensor data.

        Returns:
            6x6 array corresponding to the elastic tensor from the OUTCAR.
        """
        header_pattern = r"TOTAL ELASTIC MODULI \(kBar\)\s+" \
                         r"Direction\s+([X-Z][X-Z]\s+)+" \
                         r"\-+"
        row_pattern = r"[X-Z][X-Z]\s+" + r"\s+".join([r"(\-*[\.\d]+)"] * 6)
        footer_pattern = r"\-+"
        et_table = self.read_table_pattern(header_pattern, row_pattern,
                                           footer_pattern, postprocess=float)
        self.data["elastic_tensor"] = et_table

    def read_piezo_tensor(self):
        """
        Parse the piezo tensor data
        """
        header_pattern = r"PIEZOELECTRIC TENSOR  for field in x, y, " \
                         r"z\s+\(C/m\^2\)\s+([X-Z][X-Z]\s+)+\-+"
        row_pattern = r"[x-z]\s+" + r"\s+".join([r"(\-*[\.\d]+)"] * 6)
        footer_pattern = r"BORN EFFECTIVE"
        pt_table = self.read_table_pattern(header_pattern, row_pattern,
                                           footer_pattern, postprocess=float)
        self.data["piezo_tensor"] = pt_table

    def read_onsite_density_matrices(self):
        """
        Parse the onsite density matrices, returns list with index corresponding
        to atom index in Structure.
        """

        # matrix size will vary depending on if d or f orbitals are present
        # therefore regex assumes f, but filter out None values if d

        header_pattern = r"spin component  1\n"
        row_pattern = r'[^\S\r\n]*(?:(-?[\d.]+))' + r'(?:[^\S\r\n]*(-?[\d.]+)[^\S\r\n]*)?' * 6 + r'.*?'
        footer_pattern = r"\nspin component  2"
        spin1_component = self.read_table_pattern(header_pattern, row_pattern,
                                                  footer_pattern, postprocess=lambda x: float(x) if x else None,
                                                  last_one_only=False)

        # filter out None values
        spin1_component = [[[e for e in row if e is not None] for row in matrix] for matrix in spin1_component]

        # and repeat for Spin.down

        header_pattern = r"spin component  2\n"
        row_pattern = r'[^\S\r\n]*(?:([\d.-]+))' + r'(?:[^\S\r\n]*(-?[\d.]+)[^\S\r\n]*)?' * 6 + r'.*?'
        footer_pattern = r"\n occupancies and eigenvectors"
        spin2_component = self.read_table_pattern(header_pattern, row_pattern,
                                                  footer_pattern, postprocess=lambda x: float(x) if x else None,
                                                  last_one_only=False)

        spin2_component = [[[e for e in row if e is not None] for row in matrix] for matrix in spin2_component]

        self.data["onsite_density_matrices"] = [
            {
                Spin.up: spin1_component[idx],
                Spin.down: spin2_component[idx]
            }
            for idx in range(len(spin1_component))
        ]

    def read_corrections(self, reverse=True, terminate_on_match=True):
        """
        Reads the dipol qudropol corrections into the
        Outcar.data["dipol_quadrupol_correction"].

        :param reverse: Whether to start from end of OUTCAR.
        :param terminate_on_match: Whether to terminate once match is found.
        """
        patterns = {
            "dipol_quadrupol_correction": r"dipol\+quadrupol energy "
                                          r"correction\s+([\d\-\.]+)"
        }
        self.read_pattern(patterns, reverse=reverse,
                          terminate_on_match=terminate_on_match,
                          postprocess=float)
        self.data["dipol_quadrupol_correction"] = self.data["dipol_quadrupol_correction"][0][0]

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
            "energy": r"energy\(sigma->0\)\s+=\s+([\d\-\.]+)",
            "tangent_force": r"(NEB: projections on to tangent \(spring, REAL\)\s+\S+|tangential force \(eV/A\))\s+"
                             r"([\d\-\.]+)"
        }
        self.read_pattern(patterns, reverse=reverse,
                          terminate_on_match=terminate_on_match,
                          postprocess=str)
        self.data["energy"] = float(self.data["energy"][0][0])
        if self.data.get("tangent_force"):
            self.data["tangent_force"] = float(
                self.data["tangent_force"][0][1])

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

            search.append([r"^ *e<r>_ev=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           r"*([-0-9.Ee+]*) *\)",
                           None, er_ev])

            def er_bp(results, match):
                results.er_bp[Spin.up] = np.array([float(match.group(i))
                                                   for i in range(1, 4)]) / 2
                results.er_bp[Spin.down] = results.er_bp[Spin.up]

            search.append([r"^ *e<r>_bp=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           r"*([-0-9.Ee+]*) *\)",
                           lambda results, line: results.context == 2, er_bp])

            # Spin cases
            def er_ev_up(results, match):
                results.er_ev[Spin.up] = np.array([float(match.group(i))
                                                   for i in range(1, 4)])
                results.context = Spin.up

            search.append([r"^.*Spin component 1 *e<r>_ev=\( *([-0-9.Ee+]*) "
                           r"*([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)",
                           None, er_ev_up])

            def er_bp_up(results, match):
                results.er_bp[Spin.up] = np.array([float(match.group(1)),
                                                   float(match.group(2)),
                                                   float(match.group(3))])

            search.append([r"^ *e<r>_bp=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           r"*([-0-9.Ee+]*) *\)",
                           lambda results, line: results.context == Spin.up, er_bp_up])

            def er_ev_dn(results, match):
                results.er_ev[Spin.down] = np.array([float(match.group(1)),
                                                     float(match.group(2)),
                                                     float(match.group(3))])
                results.context = Spin.down

            search.append([r"^.*Spin component 2 *e<r>_ev=\( *([-0-9.Ee+]*) "
                           r"*([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)",
                           None, er_ev_dn])

            def er_bp_dn(results, match):
                results.er_bp[Spin.down] = np.array([float(match.group(i))
                                                     for i in range(1, 4)])

            search.append([r"^ *e<r>_bp=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           r"*([-0-9.Ee+]*) *\)",
                           lambda results, line: results.context == Spin.down, er_bp_dn])

            # Always present spin/non-spin
            def p_elc(results, match):
                results.p_elc = np.array([float(match.group(i)) for i in range(1, 4)])

            search.append([r"^.*Total electronic dipole moment: "
                           r"*p\[elc\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           r"*([-0-9.Ee+]*) *\)", None, p_elc])

            def p_ion(results, match):
                results.p_ion = np.array([float(match.group(i)) for i in range(1, 4)])

            search.append([r"^.*ionic dipole moment: "
                           r"*p\[ion\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           r"*([-0-9.Ee+]*) *\)", None, p_ion])

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

        except Exception:
            self.er_ev_tot = None
            self.er_bp_tot = None
            raise Exception("IGPAR OUTCAR could not be parsed.")

    def read_internal_strain_tensor(self):
        """
        Reads the internal strain tensor and populates self.internal_strain_tensor with an array of voigt notation
            tensors for each site.
        """
        search = []

        def internal_strain_start(results, match):
            results.internal_strain_ion = int(match.group(1)) - 1
            results.internal_strain_tensor.append(np.zeros((3, 6)))

        search.append([r"INTERNAL STRAIN TENSOR FOR ION\s+(\d+)\s+for displacements in x,y,z  \(eV/Angst\):",
                       None, internal_strain_start])

        def internal_strain_data(results, match):
            if match.group(1).lower() == "x":
                index = 0
            elif match.group(1).lower() == "y":
                index = 1
            elif match.group(1).lower() == "z":
                index = 2
            else:
                raise Exception(
                    "Couldn't parse row index from symbol for internal strain tensor: {}".format(match.group(1)))
            results.internal_strain_tensor[results.internal_strain_ion][index] = np.array([float(match.group(i))
                                                                                           for i in range(2, 8)])
            if index == 2:
                results.internal_strain_ion = None

        search.append([r"^\s+([x,y,z])\s+" + r"([-]?\d+\.\d+)\s+" * 6,
                       lambda results, line: results.internal_strain_ion is not None,
                       internal_strain_data])

        self.internal_strain_ion = None
        self.internal_strain_tensor = []
        micro_pyawk(self.filename, search, self)

    def read_lepsilon(self):
        """
        Reads an LEPSILON run.

        # TODO: Document the actual variables.
        """
        try:
            search = []

            def dielectric_section_start(results, match):
                results.dielectric_index = -1

            search.append([r"MACROSCOPIC STATIC DIELECTRIC TENSOR \(", None,
                           dielectric_section_start])

            def dielectric_section_start2(results, match):
                results.dielectric_index = 0

            search.append(
                [r"-------------------------------------",
                 lambda results, line: results.dielectric_index == -1,
                 dielectric_section_start2])

            def dielectric_data(results, match):
                results.dielectric_tensor[results.dielectric_index, :] = \
                    np.array([float(match.group(i)) for i in range(1, 4)])
                results.dielectric_index += 1

            search.append(
                [r"^ *([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+) *$",
                 lambda results, line: results.dielectric_index >= 0
                 if results.dielectric_index is not None
                 else None,
                 dielectric_data])

            def dielectric_section_stop(results, match):
                results.dielectric_index = None

            search.append(
                [r"-------------------------------------",
                 lambda results, line: results.dielectric_index >= 1
                 if results.dielectric_index is not None
                 else None,
                 dielectric_section_stop])

            self.dielectric_index = None
            self.dielectric_tensor = np.zeros((3, 3))

            def piezo_section_start(results, match):
                results.piezo_index = 0

            search.append([r"PIEZOELECTRIC TENSOR  for field in x, y, z        "
                           r"\(C/m\^2\)",
                           None, piezo_section_start])

            def piezo_data(results, match):
                results.piezo_tensor[results.piezo_index, :] = \
                    np.array([float(match.group(i)) for i in range(1, 7)])
                results.piezo_index += 1

            search.append(
                [r"^ *[xyz] +([-0-9.Ee+]+) +([-0-9.Ee+]+)" +
                 r" +([-0-9.Ee+]+) *([-0-9.Ee+]+) +([-0-9.Ee+]+)" +
                 r" +([-0-9.Ee+]+)*$",
                 lambda results, line: results.piezo_index >= 0
                 if results.piezo_index is not None
                 else None,
                 piezo_data])

            def piezo_section_stop(results, match):
                results.piezo_index = None

            search.append(
                [r"-------------------------------------",
                 lambda results, line: results.piezo_index >= 1
                 if results.piezo_index is not None
                 else None,
                 piezo_section_stop])

            self.piezo_index = None
            self.piezo_tensor = np.zeros((3, 6))

            def born_section_start(results, match):
                results.born_ion = -1

            search.append([r"BORN EFFECTIVE CHARGES ",
                           None, born_section_start])

            def born_ion(results, match):
                results.born_ion = int(match.group(1)) - 1
                results.born.append(np.zeros((3, 3)))

            search.append([r"ion +([0-9]+)", lambda results, line: results.born_ion is not None, born_ion])

            def born_data(results, match):
                results.born[results.born_ion][int(match.group(1)) - 1, :] = \
                    np.array([float(match.group(i)) for i in range(2, 5)])

            search.append(
                [r"^ *([1-3]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+)$",
                 lambda results, line: results.born_ion >= 0
                 if results.born_ion is not None
                 else results.born_ion,
                 born_data])

            def born_section_stop(results, match):
                results.born_ion = None

            search.append(
                [r"-------------------------------------",
                 lambda results, line: results.born_ion >= 1
                 if results.born_ion is not None
                 else results.born_ion,
                 born_section_stop])

            self.born_ion = None
            self.born = []

            micro_pyawk(self.filename, search, self)

            self.born = np.array(self.born)

            self.dielectric_tensor = self.dielectric_tensor.tolist()
            self.piezo_tensor = self.piezo_tensor.tolist()

        except Exception:
            raise Exception("LEPSILON OUTCAR could not be parsed.")

    def read_lepsilon_ionic(self):
        """
        Reads an LEPSILON run, the ionic component.

        # TODO: Document the actual variables.
        """
        try:
            search = []

            def dielectric_section_start(results, match):
                results.dielectric_ionic_index = -1

            search.append([r"MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC", None,
                           dielectric_section_start])

            def dielectric_section_start2(results, match):
                results.dielectric_ionic_index = 0

            search.append(
                [r"-------------------------------------",
                 lambda results, line: results.dielectric_ionic_index == -1
                 if results.dielectric_ionic_index is not None
                 else results.dielectric_ionic_index,
                 dielectric_section_start2])

            def dielectric_data(results, match):
                results.dielectric_ionic_tensor[results.dielectric_ionic_index, :] = \
                    np.array([float(match.group(i)) for i in range(1, 4)])
                results.dielectric_ionic_index += 1

            search.append(
                [r"^ *([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+) *$",
                 lambda results, line: results.dielectric_ionic_index >= 0
                 if results.dielectric_ionic_index is not None
                 else results.dielectric_ionic_index,
                 dielectric_data])

            def dielectric_section_stop(results, match):
                results.dielectric_ionic_index = None

            search.append(
                [r"-------------------------------------",
                 lambda results, line: results.dielectric_ionic_index >= 1
                 if results.dielectric_ionic_index is not None
                 else results.dielectric_ionic_index,
                 dielectric_section_stop])

            self.dielectric_ionic_index = None
            self.dielectric_ionic_tensor = np.zeros((3, 3))

            def piezo_section_start(results, match):
                results.piezo_ionic_index = 0

            search.append([r"PIEZOELECTRIC TENSOR IONIC CONTR  for field in "
                           r"x, y, z        ",
                           None, piezo_section_start])

            def piezo_data(results, match):
                results.piezo_ionic_tensor[results.piezo_ionic_index, :] = \
                    np.array([float(match.group(i)) for i in range(1, 7)])
                results.piezo_ionic_index += 1

            search.append(
                [r"^ *[xyz] +([-0-9.Ee+]+) +([-0-9.Ee+]+)" +
                 r" +([-0-9.Ee+]+) *([-0-9.Ee+]+) +([-0-9.Ee+]+)" +
                 r" +([-0-9.Ee+]+)*$",
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

        except Exception:
            raise Exception(
                "ionic part of LEPSILON OUTCAR could not be parsed.")

    def read_lcalcpol(self):
        """
        Reads the lcalpol.

        # TODO: Document the actual variables.
        """
        self.p_elec = None
        self.p_sp1 = None
        self.p_sp2 = None
        self.p_ion = None
        try:
            search = []

            # Always present spin/non-spin
            def p_elec(results, match):
                results.p_elec = np.array([float(match.group(1)),
                                           float(match.group(2)),
                                           float(match.group(3))])

            search.append([r"^.*Total electronic dipole moment: "
                           r"*p\[elc\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                           r"*([-0-9.Ee+]*) *\)",
                           None, p_elec])

            # If spin-polarized (and not noncollinear)
            # save spin-polarized electronic values
            if self.spin and not self.noncollinear:
                def p_sp1(results, match):
                    results.p_sp1 = np.array([float(match.group(1)),
                                              float(match.group(2)),
                                              float(match.group(3))])

                search.append([r"^.*p\[sp1\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                               r"*([-0-9.Ee+]*) *\)",
                               None, p_sp1])

                def p_sp2(results, match):
                    results.p_sp2 = np.array([float(match.group(1)),
                                              float(match.group(2)),
                                              float(match.group(3))])

                search.append([r"^.*p\[sp2\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                               r"*([-0-9.Ee+]*) *\)",
                               None, p_sp2])

            def p_ion(results, match):
                results.p_ion = np.array([float(match.group(1)),
                                          float(match.group(2)),
                                          float(match.group(3))])

            search.append([r"^.*Ionic dipole moment: *p\[ion\]="
                           r"\( *([-0-9.Ee+]*)"
                           r" *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)",
                           None, p_ion])

            micro_pyawk(self.filename, search, self)

        except Exception:
            raise Exception("LCALCPOL OUTCAR could not be parsed.")

    def read_pseudo_zval(self):
        """
        Create pseudopotential ZVAL dictionary.
        """
        # pylint: disable=E1101
        try:
            def atom_symbols(results, match):
                element_symbol = match.group(1)
                if not hasattr(results, 'atom_symbols'):
                    results.atom_symbols = []
                results.atom_symbols.append(element_symbol.strip())

            def zvals(results, match):
                zvals = match.group(1)
                results.zvals = map(float, re.findall(r'-?\d+\.\d*', zvals))

            search = []
            search.append([r'(?<=VRHFIN =)(.*)(?=:)', None, atom_symbols])
            search.append([r'^\s+ZVAL.*=(.*)', None, zvals])

            micro_pyawk(self.filename, search, self)

            zval_dict = {}
            for x, y in zip(self.atom_symbols, self.zvals):
                zval_dict.update({x: y})
            self.zval_dict = zval_dict

            # Clean-up
            del self.atom_symbols
            del self.zvals
        except Exception:
            raise Exception("ZVAL dict could not be parsed.")

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
                    iat = -1
                    while line != "":
                        line = foutcar.readline()
                        # don't know number of lines to parse without knowing
                        # specific species, so stop parsing when we reach
                        # "E-fermi" instead
                        if "E-fermi" in line:
                            break
                        data = line.split()
                        # data will contain odd number of elements if it is
                        # the start of a new entry, or even number of elements
                        # if it continues the previous entry
                        if len(data) % 2 == 1:
                            iat += 1  # started parsing a new ion
                            data = data[1:]  # remove element with ion number
                        for i in range(0, len(data), 2):
                            cl[iat][data[i]].append(float(data[i + 1]))
        return cl

    def read_avg_core_poten(self):
        """
        Read the core potential at each ionic step.

        Returns:
            A list for each ionic step containing a list of the average core
            potentials for each atom: [[avg core pot]].

        Example:
            The average core potential of the 2nd atom of the structure at the
            last ionic step is: [-1][1]
        """

        def pairwise(iterable):
            """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
            a = iter(iterable)
            return zip(a, a)

        with zopen(self.filename, "rt") as foutcar:
            line = foutcar.readline()
            aps = []
            while line != "":
                line = foutcar.readline()
                if "the norm of the test charge is" in line:
                    ap = []
                    while line != "":
                        line = foutcar.readline()
                        # don't know number of lines to parse without knowing
                        # specific species, so stop parsing when we reach
                        # "E-fermi" instead
                        if "E-fermi" in line:
                            aps.append(ap)
                            break
                        data = line.split()
                        # the average core potentials of up to 5 elements are
                        # given per line
                        for i, pot in pairwise(data):
                            ap.append(float(pot))
        return aps

    def as_dict(self):
        """
        :return: MSONAble dict.
        """
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__, "efermi": self.efermi,
             "run_stats": self.run_stats, "magnetization": self.magnetization,
             "charge": self.charge, "total_magnetization": self.total_mag,
             "nelect": self.nelect, "is_stopped": self.is_stopped,
             "drift": self.drift, "ngf": self.ngf,
             "sampling_radii": self.sampling_radii,
             "electrostatic_potential": self.electrostatic_potential}

        if self.lepsilon:
            d.update({"piezo_tensor": self.piezo_tensor,
                      "dielectric_tensor": self.dielectric_tensor,
                      "born": self.born})

        if self.dfpt:
            d.update({"internal_strain_tensor": self.internal_strain_tensor})

        if self.dfpt and self.lepsilon:
            d.update({"piezo_ionic_tensor": self.piezo_ionic_tensor,
                      "dielectric_ionic_tensor": self.dielectric_ionic_tensor})

        if self.lcalcpol:
            d.update({'p_elec': self.p_elec,
                      'p_ion': self.p_ion})
            if self.spin and not self.noncollinear:
                d.update({'p_sp1': self.p_sp1,
                          'p_sp2': self.p_sp2})
            d.update({'zval_dict': self.zval_dict})

        if self.nmr_cs:
            d.update({"nmr_cs": {"valence and core": self.data["chemical_shielding"]["valence_and_core"],
                                 "valence_only": self.data["chemical_shielding"]["valence_only"],
                                 "g0": self.data["cs_g0_contribution"],
                                 "core": self.data["cs_core_contribution"],
                                 "raw": self.data["unsym_cs_tensor"]}})

        if self.nmr_efg:
            d.update({"nmr_efg": {"raw": self.data["unsym_efg_tensor"],
                                  "parameters": self.data["efg"]}})

        if self.has_onsite_density_matrices:
            # cast Spin to str for consistency with electronic_structure
            # TODO: improve handling of Enum (de)serialization in monty
            onsite_density_matrices = [{str(k): v for k, v in d.items()}
                                       for d in self.data["onsite_density_matrices"]]
            d.update({"onsite_density_matrices": onsite_density_matrices})

        return d

    def read_fermi_contact_shift(self):
        """
        output example:
        Fermi contact (isotropic) hyperfine coupling parameter (MHz)
        -------------------------------------------------------------
        ion      A_pw      A_1PS     A_1AE     A_1c      A_tot
        -------------------------------------------------------------
         1      -0.002    -0.002    -0.051     0.000    -0.052
         2      -0.002    -0.002    -0.051     0.000    -0.052
         3       0.056     0.056     0.321    -0.048     0.321
        -------------------------------------------------------------
        , which corresponds to
        [[-0.002, -0.002, -0.051, 0.0, -0.052],
         [-0.002, -0.002, -0.051, 0.0, -0.052],
         [0.056, 0.056, 0.321, -0.048, 0.321]] from 'fch' data
        """

        # Fermi contact (isotropic) hyperfine coupling parameter (MHz)
        header_pattern1 = r"\s*Fermi contact \(isotropic\) hyperfine coupling parameter \(MHz\)\s+" \
                          r"\s*\-+" \
                          r"\s*ion\s+A_pw\s+A_1PS\s+A_1AE\s+A_1c\s+A_tot\s+" \
                          r"\s*\-+"
        row_pattern1 = r'(?:\d+)\s+' + r'\s+'.join([r'([-]?\d+\.\d+)'] * 5)
        footer_pattern = r"\-+"
        fch_table = self.read_table_pattern(header_pattern1, row_pattern1,
                                            footer_pattern, postprocess=float,
                                            last_one_only=True)

        # Dipolar hyperfine coupling parameters (MHz)
        header_pattern2 = r"\s*Dipolar hyperfine coupling parameters \(MHz\)\s+" \
                          r"\s*\-+" \
                          r"\s*ion\s+A_xx\s+A_yy\s+A_zz\s+A_xy\s+A_xz\s+A_yz\s+" \
                          r"\s*\-+"
        row_pattern2 = r'(?:\d+)\s+' + r'\s+'.join([r'([-]?\d+\.\d+)'] * 6)
        dh_table = self.read_table_pattern(header_pattern2, row_pattern2,
                                           footer_pattern, postprocess=float,
                                           last_one_only=True)

        # Total hyperfine coupling parameters after diagonalization (MHz)
        header_pattern3 = r"\s*Total hyperfine coupling parameters after diagonalization \(MHz\)\s+" \
                          r"\s*\(convention: \|A_zz\| > \|A_xx\| > \|A_yy\|\)\s+" \
                          r"\s*\-+" \
                          r"\s*ion\s+A_xx\s+A_yy\s+A_zz\s+asymmetry \(A_yy - A_xx\)/ A_zz\s+" \
                          r"\s*\-+"
        row_pattern3 = r'(?:\d+)\s+' + r'\s+'.join([r'([-]?\d+\.\d+)'] * 4)
        th_table = self.read_table_pattern(header_pattern3, row_pattern3,
                                           footer_pattern, postprocess=float,
                                           last_one_only=True)

        fc_shift_table = {'fch': fch_table, 'dh': dh_table, 'th': th_table}

        self.data["fermi_contact_shift"] = fc_shift_table


class VolumetricData(MSONable):
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

    def __init__(self, structure, data, distance_matrix=None, data_aug=None):
        """
        Typically, this constructor is not used directly and the static
        from_file constructor is used. This constructor is designed to allow
        summation and other operations between VolumetricData objects.

        Args:
            structure: Structure associated with the volumetric data
            data: Actual volumetric data.
            data_aug: Any extra information associated with volumetric data
                (typically augmentation charges)
            distance_matrix: A pre-computed distance matrix if available.
                Useful so pass distance_matrices between sums,
                shortcircuiting an otherwise expensive operation.
        """
        self.structure = structure
        self.is_spin_polarized = len(data) >= 2
        self.is_soc = len(data) >= 4
        self.dim = data["total"].shape
        self.data = data
        self.data_aug = data_aug if data_aug else {}
        self.ngridpts = self.dim[0] * self.dim[1] * self.dim[2]
        # lazy init the spin data since this is not always needed.
        self._spin_data = {}
        self._distance_matrix = {} if not distance_matrix else distance_matrix
        self.xpoints = np.linspace(0.0, 1.0, num=self.dim[0])
        self.ypoints = np.linspace(0.0, 1.0, num=self.dim[1])
        self.zpoints = np.linspace(0.0, 1.0, num=self.dim[2])
        self.interpolator = RegularGridInterpolator((self.xpoints, self.ypoints, self.zpoints),
                                                    self.data['total'], bounds_error=True)
        self.name = "VolumetricData"

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

    def copy(self):
        """
        :return: Copy of Volumetric object
        """
        return VolumetricData(
            self.structure,
            {k: v.copy() for k, v in self.data.items()},
            distance_matrix=self._distance_matrix,
            data_aug=self.data_aug
        )

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
            warnings.warn("Structures are different. Make sure you know what "
                          "you are doing...")
        if self.data.keys() != other.data.keys():
            raise ValueError("Data have different keys! Maybe one is spin-"
                             "polarized and the other is not?")

        # To add checks
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
        # pylint: disable=E1136,E1126
        poscar_read = False
        poscar_string = []
        dataset = []
        all_dataset = []
        # for holding any strings in input that are not Poscar
        # or VolumetricData (typically augmentation charges)
        all_dataset_aug = {}
        dim = None
        dimline = None
        read_dataset = False
        ngrid_pts = 0
        data_count = 0
        poscar = None
        with zopen(filename, "rt") as f:
            for line in f:
                original_line = line
                line = line.strip()
                if read_dataset:
                    for tok in line.split():
                        if data_count < ngrid_pts:
                            # This complicated procedure is necessary because
                            # vasp outputs x as the fastest index, followed by y
                            # then z.
                            no_x = data_count // dim[0]
                            dataset[
                                data_count % dim[0],
                                no_x % dim[1],
                                no_x // dim[1]
                            ] = float(tok)
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
                    # when line == dimline, expect volumetric data to follow
                    # so set read_dataset to True
                    read_dataset = True
                    dataset = np.zeros(dim)
                else:
                    # store any extra lines that were not part of the
                    # volumetric data so we know which set of data the extra
                    # lines are associated with
                    key = len(all_dataset) - 1
                    if key not in all_dataset_aug:
                        all_dataset_aug[key] = []
                    all_dataset_aug[key].append(original_line)
            if len(all_dataset) == 4:

                data = {"total": all_dataset[0], "diff_x": all_dataset[1],
                        "diff_y": all_dataset[2], "diff_z": all_dataset[3]}
                data_aug = {"total": all_dataset_aug.get(0, None),
                            "diff_x": all_dataset_aug.get(1, None),
                            "diff_y": all_dataset_aug.get(2, None),
                            "diff_z": all_dataset_aug.get(3, None)}

                # construct a "diff" dict for scalar-like magnetization density,
                # referenced to an arbitrary direction (using same method as
                # pymatgen.electronic_structure.core.Magmom, see
                # Magmom documentation for justification for this)
                # TODO: re-examine this, and also similar behavior in
                # Magmom - @mkhorton
                # TODO: does CHGCAR change with different SAXIS?
                diff_xyz = np.array([data["diff_x"], data["diff_y"],
                                     data["diff_z"]])
                diff_xyz = diff_xyz.reshape((3, dim[0] * dim[1] * dim[2]))
                ref_direction = np.array([1.01, 1.02, 1.03])
                ref_sign = np.sign(np.dot(ref_direction, diff_xyz))
                diff = np.multiply(np.linalg.norm(diff_xyz, axis=0), ref_sign)
                data["diff"] = diff.reshape((dim[0], dim[1], dim[2]))

            elif len(all_dataset) == 2:
                data = {"total": all_dataset[0], "diff": all_dataset[1]}
                data_aug = {"total": all_dataset_aug.get(0, None),
                            "diff": all_dataset_aug.get(1, None)}
            else:
                data = {"total": all_dataset[0]}
                data_aug = {"total": all_dataset_aug.get(0, None)}
            return poscar, data, data_aug

    def write_file(self, file_name, vasp4_compatible=False):
        """
        Write the VolumetricData object to a vasp compatible file.

        Args:
            file_name (str): Path to a file
            vasp4_compatible (bool): True if the format is vasp4 compatible
        """

        def _print_fortran_float(f):
            """
            Fortran codes print floats with a leading zero in scientific
            notation. When writing CHGCAR files, we adopt this convention
            to ensure written CHGCAR files are byte-to-byte identical to
            their input files as far as possible.
            :param f: float
            :return: str
            """
            s = "{:.10E}".format(f)
            if f >= 0:
                return "0." + s[0] + s[2:12] + 'E' + "{:+03}".format(int(s[13:]) + 1)
            return "-." + s[1] + s[3:13] + 'E' + "{:+03}".format(int(s[14:]) + 1)

        with zopen(file_name, "wt") as f:
            p = Poscar(self.structure)

            # use original name if it's been set (e.g. from Chgcar)
            comment = getattr(self, 'name', p.comment)

            lines = comment + "\n"
            lines += "   1.00000000000000\n"
            latt = self.structure.lattice.matrix
            lines += " %12.6f%12.6f%12.6f\n" % tuple(latt[0, :])
            lines += " %12.6f%12.6f%12.6f\n" % tuple(latt[1, :])
            lines += " %12.6f%12.6f%12.6f\n" % tuple(latt[2, :])
            if not vasp4_compatible:
                lines += "".join(["%5s" % s for s in p.site_symbols]) + "\n"
            lines += "".join(["%6d" % x for x in p.natoms]) + "\n"
            lines += "Direct\n"
            for site in self.structure:
                lines += "%10.6f%10.6f%10.6f\n" % tuple(site.frac_coords)
            lines += " \n"
            f.write(lines)
            a = self.dim

            def write_spin(data_type):
                lines = []
                count = 0
                f.write("   {}   {}   {}\n".format(a[0], a[1], a[2]))
                for (k, j, i) in itertools.product(list(range(a[2])),
                                                   list(range(a[1])),
                                                   list(range(a[0]))):
                    lines.append(_print_fortran_float(self.data[data_type][i, j, k]))
                    count += 1
                    if count % 5 == 0:
                        f.write(" " + "".join(lines) + "\n")
                        lines = []
                    else:
                        lines.append(" ")
                if count % 5 != 0:
                    f.write(" " + "".join(lines) + " \n")
                f.write("".join(self.data_aug.get(data_type, [])))

            write_spin("total")
            if self.is_spin_polarized and self.is_soc:
                write_spin("diff_x")
                write_spin("diff_y")
                write_spin("diff_z")
            elif self.is_spin_polarized:
                write_spin("diff")

    def value_at(self, x, y, z):
        """
        Get a data value from self.data at a given point (x, y, z) in terms
        of fractional lattice parameters. Will be interpolated using a
        RegularGridInterpolator on self.data if (x, y, z) is not in the original
        set of data points.

        Args:
            x (float): Fraction of lattice vector a.
            y (float): Fraction of lattice vector b.
            z (float): Fraction of lattice vector c.

        Returns:
            Value from self.data (potentially interpolated) correspondisng to
            the point (x, y, z).
        """
        return self.interpolator([x, y, z])[0]

    def linear_slice(self, p1, p2, n=100):
        """
        Get a linear slice of the volumetric data with n data points from
        point p1 to point p2, in the form of a list.

        Args:
            p1 (list): 3-element list containing fractional coordinates of the first point.
            p2 (list): 3-element list containing fractional coordinates of the second point.
            n (int): Number of data points to collect, defaults to 100.

        Returns:
            List of n data points (mostly interpolated) representing a linear slice of the
            data from point p1 to point p2.
        """
        assert type(p1) in [list, np.ndarray] and type(p2) in [list, np.ndarray]
        assert len(p1) == 3 and len(p2) == 3
        xpts = np.linspace(p1[0], p2[0], num=n)
        ypts = np.linspace(p1[1], p2[1], num=n)
        zpts = np.linspace(p1[2], p2[2], num=n)
        return [self.value_at(xpts[i], ypts[i], zpts[i]) for i in range(n)]

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
        # For non-spin-polarized runs, this is zero by definition.
        if not self.is_spin_polarized:
            radii = [radius / nbins * (i + 1) for i in range(nbins)]
            data = np.zeros((nbins, 2))
            data[:, 0] = radii
            return data

        struct = self.structure
        a = self.dim
        if ind not in self._distance_matrix or \
                self._distance_matrix[ind]["max_radius"] < radius:
            coords = []
            for (x, y, z) in itertools.product(*[list(range(i)) for i in a]):
                coords.append([x / a[0], y / a[1], z / a[2]])
            sites_dist = struct.lattice.get_points_in_sphere(
                coords, struct[ind].coords, radius)
            self._distance_matrix[ind] = {"max_radius": radius,
                                          "data": np.array(sites_dist)}

        data = self._distance_matrix[ind]["data"]

        # Use boolean indexing to find all charges within the desired distance.
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

    def to_hdf5(self, filename):
        """
        Writes the VolumetricData to a HDF5 format, which is a highly optimized
        format for reading storing large data. The mapping of the VolumetricData
        to this file format is as follows:

        VolumetricData.data -> f["vdata"]
        VolumetricData.structure ->
            f["Z"]: Sequence of atomic numbers
            f["fcoords"]: Fractional coords
            f["lattice"]: Lattice in the pymatgen.core.lattice.Lattice matrix
                format
            f.attrs["structure_json"]: String of json representation

        Args:
            filename (str): Filename to output to.
        """
        import h5py
        with h5py.File(filename, "w") as f:
            ds = f.create_dataset("lattice", (3, 3), dtype='float')
            ds[...] = self.structure.lattice.matrix
            ds = f.create_dataset("Z", (len(self.structure.species),),
                                  dtype="i")
            ds[...] = np.array([sp.Z for sp in self.structure.species])
            ds = f.create_dataset("fcoords", self.structure.frac_coords.shape,
                                  dtype='float')
            ds[...] = self.structure.frac_coords
            dt = h5py.special_dtype(vlen=str)
            ds = f.create_dataset("species", (len(self.structure.species),),
                                  dtype=dt)
            ds[...] = [str(sp) for sp in self.structure.species]
            grp = f.create_group("vdata")
            for k, v in self.data.items():
                ds = grp.create_dataset(k, self.data[k].shape, dtype='float')
                ds[...] = self.data[k]
            f.attrs["name"] = self.name
            f.attrs["structure_json"] = json.dumps(self.structure.as_dict())

    @classmethod
    def from_hdf5(cls, filename, **kwargs):
        """
        Reads VolumetricData from HDF5 file.

        :param filename: Filename
        :return: VolumetricData
        """
        import h5py
        with h5py.File(filename, "r") as f:
            data = {k: np.array(v) for k, v in f["vdata"].items()}
            data_aug = None
            if 'vdata_aug' in f:
                data_aug = {k: np.array(v) for k, v in f["vdata_aug"].items()}
            structure = Structure.from_dict(json.loads(f.attrs["structure_json"]))
            return cls(structure, data=data, data_aug=data_aug, **kwargs)


class Locpot(VolumetricData):
    """
    Simple object for reading a LOCPOT file.
    """

    def __init__(self, poscar, data):
        """
        Args:
            poscar (Poscar): Poscar object containing structure.
            data: Actual data.
        """
        super().__init__(poscar.structure, data)
        self.name = poscar.comment

    @classmethod
    def from_file(cls, filename, **kwargs):
        """
        Reads a LOCPOT file.

        :param filename: Filename
        :return: Locpot
        """
        (poscar, data, data_aug) = VolumetricData.parse_file(filename)
        return cls(poscar, data, **kwargs)


class Chgcar(VolumetricData):
    """
    Simple object for reading a CHGCAR file.
    """

    def __init__(self, poscar, data, data_aug=None):
        """
        Args:
            poscar (Poscar or Structure): Object containing structure.
            data: Actual data.
            data_aug: Augmentation charge data
        """
        # allow for poscar or structure files to be passed
        if isinstance(poscar, Poscar):
            tmp_struct = poscar.structure
            self.poscar = poscar
            self.name = poscar.comment
        elif isinstance(poscar, Structure):
            tmp_struct = poscar
            self.poscar = Poscar(poscar)
            self.name = None

        super().__init__(tmp_struct, data, data_aug=data_aug)
        self._distance_matrix = {}

    @staticmethod
    def from_file(filename):
        """
        Reads a CHGCAR file.

        :param filename: Filename
        :return: Chgcar
        """
        (poscar, data, data_aug) = VolumetricData.parse_file(filename)
        return Chgcar(poscar, data, data_aug=data_aug)

    @property
    def net_magnetization(self):
        """
        :return: Net magnetization from Chgcar
        """
        if self.is_spin_polarized:
            return np.sum(self.data['diff'])
        return None


class Elfcar(VolumetricData):
    """
    Read an ELFCAR file which contains the Electron Localization Function (ELF)
    as calculated by VASP.

    For ELF, "total" key refers to Spin.up, and "diff" refers to Spin.down.

    This also contains information on the kinetic energy density.
    """

    def __init__(self, poscar, data):
        """
        Args:
            poscar (Poscar or Structure): Object containing structure.
            data: Actual data.
        """
        # allow for poscar or structure files to be passed
        if isinstance(poscar, Poscar):
            tmp_struct = poscar.structure
            self.poscar = poscar
        elif isinstance(poscar, Structure):
            tmp_struct = poscar
            self.poscar = Poscar(poscar)

        super().__init__(tmp_struct, data)
        # TODO: modify VolumetricData so that the correct keys can be used.
        # for ELF, instead of "total" and "diff" keys we have
        # "Spin.up" and "Spin.down" keys
        # I believe this is correct, but there's not much documentation -mkhorton
        self.data = data

    @classmethod
    def from_file(cls, filename):
        """
        Reads a ELFCAR file.

        :param filename: Filename
        :return: Elfcar
        """
        (poscar, data, data_aug) = VolumetricData.parse_file(filename)
        return cls(poscar, data)

    def get_alpha(self):
        """
        Get the parameter alpha where ELF = 1/(1+alpha^2).
        """
        alpha_data = {}
        for k, v in self.data.items():
            alpha = 1 / v
            alpha = alpha - 1
            alpha = np.sqrt(alpha)
            alpha_data[k] = alpha
        return VolumetricData(self.structure, alpha_data)


class Procar:
    """
    Object for reading a PROCAR file.


    .. attribute:: data

        The PROCAR data of the form below. It should VASP uses 1-based indexing,
        but all indices are converted to 0-based here.::

            {
                spin: nd.array accessed with (k-point index, band index,
                                              ion index, orbital index)
            }

    .. attribute:: weights

        The weights associated with each k-point as an nd.array of lenght
        nkpoints.

    ..attribute:: phase_factors

        Phase factors, where present (e.g. LORBIT = 12). A dict of the form:
        {
            spin: complex nd.array accessed with (k-point index, band index,
                                                  ion index, orbital index)
        }

    ..attribute:: nbands

        Number of bands

    ..attribute:: nkpoints

        Number of k-points

    ..attribute:: nions

        Number of ions
    """

    def __init__(self, filename):
        """
        Args:
            filename: Name of file containing PROCAR.
        """
        headers = None

        with zopen(filename, "rt") as f:
            preambleexpr = re.compile(
                r"# of k-points:\s*(\d+)\s+# of bands:\s*(\d+)\s+# of "
                r"ions:\s*(\d+)")
            kpointexpr = re.compile(r"^k-point\s+(\d+).*weight = ([0-9\.]+)")
            bandexpr = re.compile(r"^band\s+(\d+)")
            ionexpr = re.compile(r"^ion.*")
            expr = re.compile(r"^([0-9]+)\s+")
            current_kpoint = 0
            current_band = 0
            done = False
            spin = Spin.down
            weights = None
            # pylint: disable=E1137
            for l in f:
                l = l.strip()
                if bandexpr.match(l):
                    m = bandexpr.match(l)
                    current_band = int(m.group(1)) - 1
                    done = False
                elif kpointexpr.match(l):
                    m = kpointexpr.match(l)
                    current_kpoint = int(m.group(1)) - 1
                    weights[current_kpoint] = float(m.group(2))
                    if current_kpoint == 0:
                        spin = Spin.up if spin == Spin.down else Spin.down
                    done = False
                elif headers is None and ionexpr.match(l):
                    headers = l.split()
                    headers.pop(0)
                    headers.pop(-1)

                    def ff():
                        return np.zeros((nkpoints, nbands, nions, len(headers)))

                    data = defaultdict(ff)

                    def f2():
                        return np.full((nkpoints, nbands, nions, len(headers)),
                                       np.NaN, dtype=np.complex128)

                    phase_factors = defaultdict(f2)
                elif expr.match(l):
                    toks = l.split()
                    index = int(toks.pop(0)) - 1
                    num_data = np.array([float(t) for t in toks[:len(headers)]])
                    if not done:
                        data[spin][current_kpoint, current_band, index, :] = num_data
                    else:
                        if len(toks) > len(headers):
                            # new format of PROCAR (vasp 5.4.4)
                            num_data = np.array([float(t)
                                                 for t in toks[:2 * len(headers)]])
                            for orb in range(len(headers)):
                                phase_factors[spin][current_kpoint, current_band,
                                                    index, orb] = complex(num_data[2 * orb], num_data[2 * orb + 1])
                        else:
                            # old format of PROCAR (vasp 5.4.1 and before)
                            if np.isnan(phase_factors[spin][current_kpoint, current_band, index, 0]):
                                phase_factors[spin][current_kpoint, current_band, index, :] = num_data
                            else:
                                phase_factors[spin][current_kpoint, current_band, index, :] += 1j * num_data
                elif l.startswith("tot"):
                    done = True
                elif preambleexpr.match(l):
                    m = preambleexpr.match(l)
                    nkpoints = int(m.group(1))
                    nbands = int(m.group(2))
                    nions = int(m.group(3))
                    weights = np.zeros(nkpoints)

            self.nkpoints = nkpoints
            self.nbands = nbands
            self.nions = nions
            self.weights = weights
            self.orbitals = headers
            self.data = data
            self.phase_factors = phase_factors

    def get_projection_on_elements(self, structure):
        """
        Method returning a dictionary of projections on elements.

        Args:
            structure (Structure): Input structure.

        Returns:
            a dictionary in the {Spin.up:[k index][b index][{Element:values}]]
        """
        dico = {}
        for spin in self.data.keys():
            dico[spin] = [[defaultdict(float)
                           for i in range(self.nkpoints)]
                          for j in range(self.nbands)]

        for iat in range(self.nions):
            name = structure.species[iat].symbol
            for spin, d in self.data.items():
                for k, b in itertools.product(range(self.nkpoints),
                                              range(self.nbands)):
                    dico[spin][b][k][name] = np.sum(d[k, b, iat, :])

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

        orbital_index = self.orbitals.index(orbital)
        return {spin: np.sum(d[:, :, atom_index, orbital_index] * self.weights[:, None])
                for spin, d in self.data.items()}


class Oszicar:
    """
    A basic parser for an OSZICAR output from VASP.  In general, while the
    OSZICAR is useful for a quick look at the output from a VASP run, we
    recommend that you use the Vasprun parser instead, which gives far richer
    information about a run.

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
        """
        Args:
            filename (str): Filename of file to parse
        """
        electronic_steps = []
        ionic_steps = []
        ionic_pattern = re.compile(r"(\d+)\s+F=\s*([\d\-\.E\+]+)\s+"
                                   r"E0=\s*([\d\-\.E\+]+)\s+"
                                   r"d\s*E\s*=\s*([\d\-\.E\+]+)$")
        ionic_mag_pattern = re.compile(r"(\d+)\s+F=\s*([\d\-\.E\+]+)\s+"
                                       r"E0=\s*([\d\-\.E\+]+)\s+"
                                       r"d\s*E\s*=\s*([\d\-\.E\+]+)\s+"
                                       r"mag=\s*([\d\-\.E\+]+)")
        ionic_MD_pattern = re.compile(r"(\d+)\s+T=\s*([\d\-\.E\+]+)\s+"
                                      r"E=\s*([\d\-\.E\+]+)\s+"
                                      r"F=\s*([\d\-\.E\+]+)\s+"
                                      r"E0=\s*([\d\-\.E\+]+)\s+"
                                      r"EK=\s*([\d\-\.E\+]+)\s+"
                                      r"SP=\s*([\d\-\.E\+]+)\s+"
                                      r"SK=\s*([\d\-\.E\+]+)")
        electronic_pattern = re.compile(r"\s*\w+\s*:(.*)")

        def smart_convert(header, num):
            try:
                if header in ("N", "ncg"):
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
                elif re.match(r"^\s*N\s+E\s*", line):
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

    @property  # type: ignore
    @unitized("eV")
    def final_energy(self):
        """
        Final energy from run.
        """
        return self.ionic_steps[-1]["E0"]

    def as_dict(self):
        """
        :return: MSONable dict
        """
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
    # TODO: Add better error handling!!!
    if os.path.exists(os.path.join(dir_name, "branch_0")):
        # get all branch dir names
        branch_dir_names = [os.path.abspath(d)
                            for d in glob.glob("{i}/branch_*"
                                               .format(i=dir_name))
                            if os.path.isdir(d)]

        # sort by the directory name (e.g, branch_10)
        sorted_branch_dir_names = sorted(branch_dir_names, key=lambda x: int(x.split("_")[-1]))

        # populate branches with Bandstructure instances
        branches = []
        for dname in sorted_branch_dir_names:
            xml_file = os.path.join(dname, "vasprun.xml")
            if os.path.exists(xml_file):
                run = Vasprun(xml_file, parse_projected_eigen=projections)
                branches.append(run.get_band_structure(efermi=efermi))
            else:
                # It might be better to throw an exception
                warnings.warn("Skipping {}. Unable to find {}".format(dname, xml_file))

        return get_reconstructed_band_structure(branches, efermi)

    xml_file = os.path.join(dir_name, "vasprun.xml")
    # Better handling of Errors
    if os.path.exists(xml_file):
        return Vasprun(xml_file, parse_projected_eigen=projections) \
            .get_band_structure(kpoints_filename=None, efermi=efermi)

    return None


class Xdatcar:
    """
    Class representing an XDATCAR file. Only tested with VASP 5.x files.

    .. attribute:: structures

        List of structures parsed from XDATCAR.
    .. attribute:: comment

        Optional comment string.
    Authors: Ram Balachandran
    """

    def __init__(self, filename, ionicstep_start=1,
                 ionicstep_end=None, comment=None):
        """
        Init a Xdatcar.

        Args:
            filename (str): Filename of input XDATCAR file.
            ionicstep_start (int): Starting number of ionic step.
            ionicstep_end (int): Ending number of ionic step.
        """
        preamble = None
        coords_str = []
        structures = []
        preamble_done = False
        if ionicstep_start < 1:
            raise Exception('Start ionic step cannot be less than 1')
        if ionicstep_end is not None and ionicstep_start < 1:
            raise Exception('End ionic step cannot be less than 1')

        # pylint: disable=E1136
        ionicstep_cnt = 1
        with zopen(filename, "rt") as f:
            for l in f:
                l = l.strip()
                if preamble is None:
                    preamble = [l]
                elif not preamble_done:
                    if l == "" or "Direct configuration=" in l:
                        preamble_done = True
                        tmp_preamble = [preamble[0]]
                        for i in range(1, len(preamble)):
                            if preamble[0] != preamble[i]:
                                tmp_preamble.append(preamble[i])
                            else:
                                break
                        preamble = tmp_preamble
                    else:
                        preamble.append(l)
                elif l == "" or "Direct configuration=" in l:
                    p = Poscar.from_string("\n".join(preamble +
                                                     ["Direct"] + coords_str))
                    if ionicstep_end is None:
                        if ionicstep_cnt >= ionicstep_start:
                            structures.append(p.structure)
                    else:
                        if ionicstep_start <= ionicstep_cnt < ionicstep_end:
                            structures.append(p.structure)
                        if ionicstep_cnt >= ionicstep_end:
                            break
                    ionicstep_cnt += 1
                    coords_str = []
                else:
                    coords_str.append(l)
            p = Poscar.from_string("\n".join(preamble +
                                             ["Direct"] + coords_str))
            if ionicstep_end is None:
                if ionicstep_cnt >= ionicstep_start:
                    structures.append(p.structure)
            else:
                if ionicstep_start <= ionicstep_cnt < ionicstep_end:
                    structures.append(p.structure)
        self.structures = structures
        self.comment = comment or self.structures[0].formula

    @property
    def site_symbols(self):
        """
        Sequence of symbols associated with the Xdatcar. Similar to 6th line in
        vasp 5+ Xdatcar.
        """
        syms = [site.specie.symbol for site in self.structures[0]]
        return [a[0] for a in itertools.groupby(syms)]

    @property
    def natoms(self):
        """
        Sequence of number of sites of each type associated with the Poscar.
        Similar to 7th line in vasp 5+ Xdatcar.
        """
        syms = [site.specie.symbol for site in self.structures[0]]
        return [len(tuple(a[1])) for a in itertools.groupby(syms)]

    def concatenate(self, filename, ionicstep_start=1,
                    ionicstep_end=None):
        """
        Concatenate structures in file to Xdatcar.

        Args:
            filename (str): Filename of XDATCAR file to be concatenated.
            ionicstep_start (int): Starting number of ionic step.
            ionicstep_end (int): Ending number of ionic step.
        TODO(rambalachandran):
           Requires a check to ensure if the new concatenating file has the
           same lattice structure and atoms as the Xdatcar class.
        """
        preamble = None
        coords_str = []
        structures = self.structures
        preamble_done = False
        if ionicstep_start < 1:
            raise Exception('Start ionic step cannot be less than 1')
        if (ionicstep_end is not None and
                ionicstep_start < 1):
            raise Exception('End ionic step cannot be less than 1')

        # pylint: disable=E1136
        ionicstep_cnt = 1
        with zopen(filename, "rt") as f:
            for l in f:
                l = l.strip()
                if preamble is None:
                    preamble = [l]
                elif not preamble_done:
                    if l == "" or "Direct configuration=" in l:
                        preamble_done = True
                        tmp_preamble = [preamble[0]]
                        for i in range(1, len(preamble)):
                            if preamble[0] != preamble[i]:
                                tmp_preamble.append(preamble[i])
                            else:
                                break
                        preamble = tmp_preamble
                    else:
                        preamble.append(l)
                elif l == "" or "Direct configuration=" in l:
                    p = Poscar.from_string("\n".join(preamble +
                                                     ["Direct"] + coords_str))
                    if ionicstep_end is None:
                        if ionicstep_cnt >= ionicstep_start:
                            structures.append(p.structure)
                    else:
                        if ionicstep_start <= ionicstep_cnt < ionicstep_end:
                            structures.append(p.structure)
                    ionicstep_cnt += 1
                    coords_str = []
                else:
                    coords_str.append(l)
            p = Poscar.from_string("\n".join(preamble +
                                             ["Direct"] + coords_str))
            if ionicstep_end is None:
                if ionicstep_cnt >= ionicstep_start:
                    structures.append(p.structure)
            else:
                if ionicstep_start <= ionicstep_cnt < ionicstep_end:
                    structures.append(p.structure)
        self.structures = structures

    def get_string(self, ionicstep_start=1,
                   ionicstep_end=None,
                   significant_figures=8):
        """
        Write  Xdatcar class to a string.

        Args:
            ionicstep_start (int): Starting number of ionic step.
            ionicstep_end (int): Ending number of ionic step.
            significant_figures (int): Number of significant figures.
        """
        if ionicstep_start < 1:
            raise Exception('Start ionic step cannot be less than 1')
        if ionicstep_end is not None and ionicstep_end < 1:
            raise Exception('End ionic step cannot be less than 1')
        latt = self.structures[0].lattice
        if np.linalg.det(latt.matrix) < 0:
            latt = Lattice(-latt.matrix)
        lines = [self.comment, "1.0", str(latt)]
        lines.append(" ".join(self.site_symbols))
        lines.append(" ".join([str(x) for x in self.natoms]))
        format_str = "{{:.{0}f}}".format(significant_figures)
        ionicstep_cnt = 1
        output_cnt = 1
        for cnt, structure in enumerate(self.structures):
            ionicstep_cnt = cnt + 1
            if ionicstep_end is None:
                if ionicstep_cnt >= ionicstep_start:
                    lines.append("Direct configuration=" +
                                 ' ' * (7 - len(str(output_cnt))) + str(output_cnt))
                    for (i, site) in enumerate(structure):
                        coords = site.frac_coords
                        line = " ".join([format_str.format(c) for c in coords])
                        lines.append(line)
                    output_cnt += 1
            else:
                if ionicstep_start <= ionicstep_cnt < ionicstep_end:
                    lines.append("Direct configuration=" +
                                 ' ' * (7 - len(str(output_cnt))) + str(output_cnt))
                    for (i, site) in enumerate(structure):
                        coords = site.frac_coords
                        line = " ".join([format_str.format(c) for c in coords])
                        lines.append(line)
                    output_cnt += 1
        return "\n".join(lines) + "\n"

    def write_file(self, filename, **kwargs):
        """
        Write  Xdatcar class into a file.

        Args:
            filename (str): Filename of output XDATCAR file.
            The supported kwargs are the same as those for the
            Xdatcar.get_string method and are passed through directly.
        """
        with zopen(filename, "wt") as f:
            f.write(self.get_string(**kwargs))

    def __str__(self):
        return self.get_string()


class Dynmat:
    """
    Object for reading a DYNMAT file.

    .. attribute:: data

        A nested dict containing the DYNMAT data of the form::
        [atom <int>][disp <int>]['dispvec'] =
            displacement vector (part of first line in dynmat block, e.g. "0.01 0 0")
        [atom <int>][disp <int>]['dynmat'] =
                <list> list of dynmat lines for this atom and this displacement

    Authors: Patrick Huck
    """

    def __init__(self, filename):
        """
        Args:
            filename: Name of file containing DYNMAT
        """
        with zopen(filename, "rt") as f:
            lines = list(clean_lines(f.readlines()))
            self._nspecs, self._natoms, self._ndisps = map(int, lines[
                0].split())
            self._masses = map(float, lines[1].split())
            self.data = defaultdict(dict)
            atom, disp = None, None
            for i, l in enumerate(lines[2:]):
                v = list(map(float, l.split()))
                if not i % (self._natoms + 1):
                    atom, disp = map(int, v[:2])
                    if atom not in self.data:
                        self.data[atom] = {}
                    if disp not in self.data[atom]:
                        self.data[atom][disp] = {}
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
        for k, v0 in self.data.items():
            for v1 in v0.itervalues():
                vec = map(abs, v1['dynmat'][k - 1])
                frequency = math.sqrt(sum(vec)) * 2. * math.pi * 15.633302  # THz
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
    # make a working copy of band_structure
    bs_working = BandStructureSymmLine.from_dict(band_structure.as_dict())
    if bs_working.is_metal():
        e = efermi
        while e < cbm:
            e += 0.01
            bs_working._efermi = e
            if not bs_working.is_metal():
                return e
    return efermi


# a note to future confused people (i.e. myself):
# I use numpy.fromfile instead of scipy.io.FortranFile here because the records
# are of fixed length, so the record length is only written once. In fortran,
# this amounts to using open(..., form='unformatted', recl=recl_len). In
# constrast when you write UNK files, the record length is written at the
# beginning of each record. This allows you to use scipy.io.FortranFile. In
# fortran, this amounts to using open(..., form='unformatted') [i.e. no recl=].
class Wavecar:
    """
    This is a class that contains the (pseudo-) wavefunctions from VASP.

    Coefficients are read from the given WAVECAR file and the corresponding
    G-vectors are generated using the algorithm developed in WaveTrans (see
    acknowledgments below). To understand how the wavefunctions are evaluated,
    please see the evaluate_wavefunc docstring.

    It should be noted that the pseudopotential augmentation is not included in
    the WAVECAR file. As a result, some caution should be exercised when
    deriving value from this information.

    The usefulness of this class is to allow the user to do projections or band
    unfolding style manipulations of the wavefunction. An example of this can
    be seen in the work of Shen et al. 2017
    (https://doi.org/10.1103/PhysRevMaterials.1.065001).

    .. attribute:: filename

        String of the input file (usually WAVECAR)

    .. attribute:: vasp_type

        String that determines VASP type the WAVECAR was generated with (either
        'std', 'gam', or 'ncl')

    .. attribute:: nk

        Number of k-points from the WAVECAR

    .. attribute:: nb

        Number of bands per k-point

    .. attribute:: encut

        Energy cutoff (used to define G_{cut})

    .. attribute:: efermi

        Fermi energy

    .. attribute:: a

        Primitive lattice vectors of the cell (e.g. a_1 = self.a[0, :])

    .. attribute:: b

        Reciprocal lattice vectors of the cell (e.g. b_1 = self.b[0, :])

    .. attribute:: vol

        The volume of the unit cell in real space

    .. attribute:: kpoints

        The list of k-points read from the WAVECAR file

    .. attribute:: band_energy

        The list of band eigenenergies (and corresponding occupancies) for
        each kpoint, where the first index corresponds to the index of the
        k-point (e.g. self.band_energy[kp])

    .. attribute:: Gpoints

        The list of generated G-points for each k-point (a double list), which
        are used with the coefficients for each k-point and band to recreate
        the wavefunction (e.g. self.Gpoints[kp] is the list of G-points for
        k-point kp). The G-points depend on the k-point and reciprocal lattice
        and therefore are identical for each band at the same k-point. Each
        G-point is represented by integer multipliers (e.g. assuming
        Gpoints[kp][n] == [n_1, n_2, n_3], then
        G_n = n_1*b_1 + n_2*b_2 + n_3*b_3)

    .. attribute:: coeffs

        The list of coefficients for each k-point and band for reconstructing
        the wavefunction. For non-spin-polarized, the first index corresponds
        to the kpoint and the second corresponds to the band (e.g.
        self.coeffs[kp][b] corresponds to k-point kp and band b). For
        spin-polarized calculations, the first index is for the spin.
        If the calculation was non-collinear, then self.coeffs[kp][b] will have
        two columns (one for each component of the spinor).

    Acknowledgments:
        This code is based upon the Fortran program, WaveTrans, written by
        R. M. Feenstra and M. Widom from the Dept. of Physics at Carnegie
        Mellon University. To see the original work, please visit:
        https://www.andrew.cmu.edu/user/feenstra/wavetrans/

    Author: Mark Turiansky
    """

    def __init__(self, filename='WAVECAR', verbose=False, precision='normal',
                 vasp_type=None):
        """
        Information is extracted from the given WAVECAR

        Args:
            filename (str): input file (default: WAVECAR)
            verbose (bool): determines whether processing information is shown
            precision (str): determines how fine the fft mesh is (normal or
                             accurate), only the first letter matters
            vasp_type (str): determines the VASP type that is used, allowed
                             values are ['std', 'gam', 'ncl']
                             (only first letter is required)
        """
        self.filename = filename
        if not (vasp_type is None or vasp_type.lower()[0] in ['s', 'g', 'n']):
            raise ValueError(f'invalid vasp_type {vasp_type}')
        self.vasp_type = vasp_type

        # c = 0.26246582250210965422
        # 2m/hbar^2 in agreement with VASP
        self._C = 0.262465831
        with open(self.filename, 'rb') as f:
            # read the header information
            recl, spin, rtag = np.fromfile(f, dtype=np.float64, count=3) \
                .astype(np.int)
            if verbose:
                print('recl={}, spin={}, rtag={}'.format(recl, spin, rtag))
            recl8 = int(recl / 8)
            self.spin = spin

            # check to make sure we have precision correct
            if rtag not in (45200, 45210, 53300, 53310):
                # note that rtag=45200 and 45210 may not work if file was actually
                # generated by old version of VASP, since that would write eigenvalues
                # and occupations in way that does not span FORTRAN records, but
                # reader below appears to assume that record boundaries can be ignored
                # (see OUTWAV vs. OUTWAV_4 in vasp fileio.F)
                raise ValueError('invalid rtag of {}'.format(rtag))

            # padding to end of fortran REC=1
            np.fromfile(f, dtype=np.float64, count=(recl8 - 3))

            # extract kpoint, bands, energy, and lattice information
            self.nk, self.nb, self.encut = \
                np.fromfile(f, dtype=np.float64, count=3).astype(np.int)
            self.a = np.fromfile(f, dtype=np.float64, count=9).reshape((3, 3))
            self.efermi = np.fromfile(f, dtype=np.float64, count=1)[0]
            if verbose:
                print('kpoints = {}, bands = {}, energy cutoff = {}, fermi '
                      'energy= {:.04f}\n'.format(self.nk, self.nb, self.encut,
                                                 self.efermi))
                print('primitive lattice vectors = \n{}'.format(self.a))

            self.vol = np.dot(self.a[0, :],
                              np.cross(self.a[1, :], self.a[2, :]))
            if verbose:
                print('volume = {}\n'.format(self.vol))

            # calculate reciprocal lattice
            b = np.array([np.cross(self.a[1, :], self.a[2, :]),
                          np.cross(self.a[2, :], self.a[0, :]),
                          np.cross(self.a[0, :], self.a[1, :])])
            b = 2 * np.pi * b / self.vol
            self.b = b
            if verbose:
                print('reciprocal lattice vectors = \n{}'.format(b))
                print('reciprocal lattice vector magnitudes = \n{}\n'
                      .format(np.linalg.norm(b, axis=1)))

            # calculate maximum number of b vectors in each direction
            self._generate_nbmax()
            if verbose:
                print('max number of G values = {}\n\n'.format(self._nbmax))
            self.ng = self._nbmax * 3 if precision.lower()[0] == 'n' else \
                self._nbmax * 4

            # padding to end of fortran REC=2
            np.fromfile(f, dtype=np.float64, count=recl8 - 13)

            # reading records
            self.Gpoints = [None for _ in range(self.nk)]
            self.kpoints = []
            if spin == 2:
                self.coeffs = [[[None for i in range(self.nb)]
                                for j in range(self.nk)] for _ in range(spin)]
                self.band_energy = [[] for _ in range(spin)]
            else:
                self.coeffs = [[None for i in range(self.nb)]
                               for j in range(self.nk)]
                self.band_energy = []

            for ispin in range(spin):
                if verbose:
                    print('reading spin {}'.format(ispin))

                for ink in range(self.nk):
                    # information for this kpoint
                    nplane = int(np.fromfile(f, dtype=np.float64, count=1)[0])
                    kpoint = np.fromfile(f, dtype=np.float64, count=3)

                    if ispin == 0:
                        self.kpoints.append(kpoint)
                    else:
                        assert np.allclose(self.kpoints[ink], kpoint)

                    if verbose:
                        print('kpoint {: 4} with {: 5} plane waves at {}'
                              .format(ink, nplane, kpoint))

                    # energy and occupation information
                    enocc = np.fromfile(f, dtype=np.float64,
                                        count=3 * self.nb).reshape((self.nb, 3))
                    if spin == 2:
                        self.band_energy[ispin].append(enocc)
                    else:
                        self.band_energy.append(enocc)

                    if verbose:
                        print('enocc =\n', enocc[:, [0, 2]])

                    # padding to end of record that contains nplane, kpoints, evals and occs
                    np.fromfile(f, dtype=np.float64, count=(recl8 - 4 - 3 * self.nb) % recl8)

                    if self.vasp_type is None:
                        (self.Gpoints[ink], extra_gpoints, extra_coeff_inds) = \
                            self._generate_G_points(kpoint, gamma=True)
                        if len(self.Gpoints[ink]) == nplane:
                            self.vasp_type = 'gam'
                        else:
                            (self.Gpoints[ink], extra_gpoints, extra_coeff_inds) = \
                                self._generate_G_points(kpoint, gamma=False)
                            self.vasp_type = \
                                'std' if len(self.Gpoints[ink]) == nplane else 'ncl'

                        if verbose:
                            print('\ndetermined vasp_type =', self.vasp_type, '\n')
                    else:
                        (self.Gpoints[ink], extra_gpoints, extra_coeff_inds) = \
                            self._generate_G_points(kpoint, gamma=(self.vasp_type.lower()[0] == 'g'))

                    if len(self.Gpoints[ink]) != nplane and 2*len(self.Gpoints[ink]) != nplane:
                        raise ValueError(f'Incorrect value of vasp_type given ({vasp_type}).'
                                         ' Please open an issue if you are certain this WAVECAR'
                                         ' was generated with the given vasp_type.')

                    self.Gpoints[ink] = \
                        np.array(self.Gpoints[ink] + extra_gpoints, dtype=np.float64)

                    # extract coefficients
                    for inb in range(self.nb):
                        if rtag in (45200, 53300):
                            data = np.fromfile(f, dtype=np.complex64, count=nplane)
                            np.fromfile(f, dtype=np.float64, count=recl8 - nplane)
                        elif rtag in (45210, 53310):
                            # this should handle double precision coefficients
                            # but I don't have a WAVECAR to test it with
                            data = np.fromfile(f, dtype=np.complex128, count=nplane)
                            np.fromfile(f, dtype=np.float64, count=recl8 - 2 * nplane)

                        extra_coeffs = []
                        if len(extra_coeff_inds) > 0:
                            # reconstruct extra coefficients missing from gamma-only executable WAVECAR
                            for G_ind in extra_coeff_inds:
                                # no idea where this factor of sqrt(2) comes from, but empirically
                                # it appears to be necessary
                                data[G_ind] /= np.sqrt(2)
                                extra_coeffs.append(np.conj(data[G_ind]))

                        if spin == 2:
                            self.coeffs[ispin][ink][inb] = \
                                np.array(list(data) + extra_coeffs, dtype=np.complex64)
                        else:
                            self.coeffs[ink][inb] = \
                                np.array(list(data) + extra_coeffs, dtype=np.complex128)

                        if self.vasp_type.lower()[0] == 'n':
                            self.coeffs[ink][inb].shape = (2, nplane//2)

    def _generate_nbmax(self) -> None:
        """
        Helper function that determines maximum number of b vectors for
        each direction.

        This algorithm is adapted from WaveTrans (see Class docstring). There
        should be no reason for this function to be called outside of
        initialization.
        """
        bmag = np.linalg.norm(self.b, axis=1)
        b = self.b

        # calculate maximum integers in each direction for G
        phi12 = np.arccos(np.dot(b[0, :], b[1, :]) / (bmag[0] * bmag[1]))
        sphi123 = np.dot(b[2, :], np.cross(b[0, :], b[1, :])) \
            / (bmag[2] * np.linalg.norm(np.cross(b[0, :], b[1, :])))
        nbmaxA = np.sqrt(self.encut * self._C) / bmag
        nbmaxA[0] /= np.abs(np.sin(phi12))
        nbmaxA[1] /= np.abs(np.sin(phi12))
        nbmaxA[2] /= np.abs(sphi123)
        nbmaxA += 1

        phi13 = np.arccos(np.dot(b[0, :], b[2, :]) / (bmag[0] * bmag[2]))
        sphi123 = np.dot(b[1, :], np.cross(b[0, :], b[2, :])) \
            / (bmag[1] * np.linalg.norm(np.cross(b[0, :], b[2, :])))
        nbmaxB = np.sqrt(self.encut * self._C) / bmag
        nbmaxB[0] /= np.abs(np.sin(phi13))
        nbmaxB[1] /= np.abs(sphi123)
        nbmaxB[2] /= np.abs(np.sin(phi13))
        nbmaxB += 1

        phi23 = np.arccos(np.dot(b[1, :], b[2, :]) / (bmag[1] * bmag[2]))
        sphi123 = np.dot(b[0, :], np.cross(b[1, :], b[2, :])) \
            / (bmag[0] * np.linalg.norm(np.cross(b[1, :], b[2, :])))
        nbmaxC = np.sqrt(self.encut * self._C) / bmag
        nbmaxC[0] /= np.abs(sphi123)
        nbmaxC[1] /= np.abs(np.sin(phi23))
        nbmaxC[2] /= np.abs(np.sin(phi23))
        nbmaxC += 1

        self._nbmax = np.max([nbmaxA, nbmaxB, nbmaxC], axis=0).astype(np.int)

    def _generate_G_points(self, kpoint: np.ndarray, gamma: bool = False) -> Tuple[List, List, List]:
        """
        Helper function to generate G-points based on nbmax.

        This function iterates over possible G-point values and determines
        if the energy is less than G_{cut}. Valid values are appended to
        the output array. This function should not be called outside of
        initialization.

        Args:
            kpoint (np.array): the array containing the current k-point value
            gamma (bool): determines if G points for gamma-point only executable
                          should be generated

        Returns:
            a list containing valid G-points
        """

        if gamma:
            kmax = self._nbmax[0] + 1
        else:
            kmax = 2 * self._nbmax[0] + 1

        gpoints = []
        extra_gpoints = []
        extra_coeff_inds = []
        G_ind = 0
        for i in range(2 * self._nbmax[2] + 1):
            i3 = i - 2 * self._nbmax[2] - 1 if i > self._nbmax[2] else i
            for j in range(2 * self._nbmax[1] + 1):
                j2 = j - 2 * self._nbmax[1] - 1 if j > self._nbmax[1] else j
                for k in range(kmax):
                    k1 = k - 2 * self._nbmax[0] - 1 if k > self._nbmax[0] else k
                    if gamma and ((k1 == 0 and j2 < 0) or (k1 == 0 and j2 == 0 and i3 < 0)):
                        continue
                    G = np.array([k1, j2, i3])
                    v = kpoint + G
                    g = np.linalg.norm(np.dot(v, self.b))
                    E = g ** 2 / self._C
                    if E < self.encut:
                        gpoints.append(G)
                        if gamma and (k1, j2, i3) != (0, 0, 0):
                            extra_gpoints.append(-G)
                            extra_coeff_inds.append(G_ind)
                        G_ind += 1
        return (gpoints, extra_gpoints, extra_coeff_inds)

    def evaluate_wavefunc(self, kpoint: int, band: int, r: np.ndarray,
                          spin: int = 0, spinor: int = 0) -> np.complex64:
        r"""
        Evaluates the wavefunction for a given position, r.

        The wavefunction is given by the k-point and band. It is evaluated
        at the given position by summing over the components. Formally,

        \psi_n^k (r) = \sum_{i=1}^N c_i^{n,k} \exp (i (k + G_i^{n,k}) \cdot r)

        where \psi_n^k is the wavefunction for the nth band at k-point k, N is
        the number of plane waves, c_i^{n,k} is the ith coefficient that
        corresponds to the nth band and k-point k, and G_i^{n,k} is the ith
        G-point corresponding to k-point k.

        NOTE: This function is very slow; a discrete fourier transform is the
        preferred method of evaluation (see Wavecar.fft_mesh).

        Args:
            kpoint (int): the index of the kpoint where the wavefunction
                            will be evaluated
            band (int): the index of the band where the wavefunction will be
                            evaluated
            r (np.array): the position where the wavefunction will be evaluated
            spin (int):  spin index for the desired wavefunction (only for
                            ISPIN = 2, default = 0)
            spinor (int): component of the spinor that is evaluated (only used
                            if vasp_type == 'ncl')
        Returns:
            a complex value corresponding to the evaluation of the wavefunction
        """
        v = self.Gpoints[kpoint] + self.kpoints[kpoint]
        u = np.dot(np.dot(v, self.b), r)
        if self.vasp_type.lower()[0] == 'n':
            c = self.coeffs[kpoint][band][spinor, :]
        elif self.spin == 2:
            c = self.coeffs[spin][kpoint][band]
        else:
            c = self.coeffs[kpoint][band]
        return np.sum(np.dot(c, np.exp(1j * u, dtype=np.complex64))) / np.sqrt(self.vol)

    def fft_mesh(self, kpoint: int, band: int, spin: int = 0, spinor: int = 0,
                 shift: bool = True) -> np.ndarray:
        """
        Places the coefficients of a wavefunction onto an fft mesh.

        Once the mesh has been obtained, a discrete fourier transform can be
        used to obtain real-space evaluation of the wavefunction. The output
        of this function can be passed directly to numpy's fft function. For
        example:

            mesh = Wavecar('WAVECAR').fft_mesh(kpoint, band)
            evals = np.fft.ifftn(mesh)

        Args:
            kpoint (int): the index of the kpoint where the wavefunction
                            will be evaluated
            band (int): the index of the band where the wavefunction will be
                            evaluated
            spin (int):  the spin of the wavefunction for the desired
                            wavefunction (only for ISPIN = 2, default = 0)
            spinor (int): component of the spinor that is evaluated (only used
                            if vasp_type == 'ncl')
            shift (bool): determines if the zero frequency coefficient is
                            placed at index (0, 0, 0) or centered
        Returns:
            a numpy ndarray representing the 3D mesh of coefficients
        """
        if self.vasp_type.lower()[0] == 'n':
            tcoeffs = self.coeffs[kpoint][band][spinor, :]
        elif self.spin == 2:
            tcoeffs = self.coeffs[spin][kpoint][band]
        else:
            tcoeffs = self.coeffs[kpoint][band]

        mesh = np.zeros(tuple(self.ng), dtype=np.complex)
        for gp, coeff in zip(self.Gpoints[kpoint], tcoeffs):
            t = tuple(gp.astype(np.int) + (self.ng / 2).astype(np.int))
            mesh[t] = coeff

        if shift:
            return np.fft.ifftshift(mesh)
        return mesh

    def get_parchg(self, poscar: Poscar, kpoint: int, band: int,
                   spin: Optional[int] = None, spinor: Optional[int] = None,
                   phase: bool = False, scale: int = 2) -> Chgcar:
        """
        Generates a Chgcar object, which is the charge density of the specified
        wavefunction.

        This function generates a Chgcar object with the charge density of the
        wavefunction specified by band and kpoint (and spin, if the WAVECAR
        corresponds to a spin-polarized calculation). The phase tag is a
        feature that is not present in VASP. For a real wavefunction, the phase
        tag being turned on means that the charge density is multiplied by the
        sign of the wavefunction at that point in space. A warning is generated
        if the phase tag is on and the chosen kpoint is not Gamma.

        Note: Augmentation from the PAWs is NOT included in this function. The
        maximal charge density will differ from the PARCHG from VASP, but the
        qualitative shape of the charge density will match.

        Args:
            poscar (pymatgen.io.vasp.inputs.Poscar): Poscar object that has the
                structure associated with the WAVECAR file
            kpoint (int): the index of the kpoint for the wavefunction
            band (int): the index of the band for the wavefunction
            spin (int): optional argument to specify the spin. If the Wavecar
                has ISPIN = 2, spin is None generates a Chgcar with total spin
                and magnetization, and spin == {0, 1} specifies just the spin
                up or down component.
            spinor (int): optional argument to specify the spinor component
                for noncollinear data wavefunctions (allowed values of None,
                0, or 1)
            phase (bool): flag to determine if the charge density is multiplied
                by the sign of the wavefunction. Only valid for real
                wavefunctions.
            scale (int): scaling for the FFT grid. The default value of 2 is at
                least as fine as the VASP default.
        Returns:
            a pymatgen.io.vasp.outputs.Chgcar object
        """
        if phase and not np.all(self.kpoints[kpoint] == 0.):
            warnings.warn('phase == True should only be used for the Gamma '
                          'kpoint! I hope you know what you\'re doing!')

        # scaling of ng for the fft grid, need to restore value at the end
        temp_ng = self.ng
        self.ng = self.ng * scale
        N = np.prod(self.ng)

        data = {}
        if self.spin == 2:
            if spin is not None:
                wfr = np.fft.ifftn(self.fft_mesh(kpoint, band, spin=spin)) * N
                den = np.abs(np.conj(wfr) * wfr)
                if phase:
                    den = np.sign(np.real(wfr)) * den
                data['total'] = den
            else:
                wfr = np.fft.ifftn(self.fft_mesh(kpoint, band, spin=0)) * N
                denup = np.abs(np.conj(wfr) * wfr)
                wfr = np.fft.ifftn(self.fft_mesh(kpoint, band, spin=1)) * N
                dendn = np.abs(np.conj(wfr) * wfr)
                data['total'] = denup + dendn
                data['diff'] = denup - dendn
        else:
            if spinor is not None:
                wfr = np.fft.ifftn(self.fft_mesh(kpoint, band, spinor=spinor)) * N
                den = np.abs(np.conj(wfr) * wfr)
            else:
                wfr = np.fft.ifftn(self.fft_mesh(kpoint, band, spinor=0)) * N
                wfr_t = np.fft.ifftn(self.fft_mesh(kpoint, band, spinor=1)) * N
                den = np.abs(np.conj(wfr) * wfr)
                den += np.abs(np.conj(wfr_t) * wfr_t)

            if phase and not (self.vasp_type.lower()[0] == 'n' and spinor is None):
                den = np.sign(np.real(wfr)) * den
            data['total'] = den

        self.ng = temp_ng
        return Chgcar(poscar, data)

    def write_unks(self, directory: str) -> None:
        """
        Write the UNK files to the given directory.

        Writes the cell-periodic part of the bloch wavefunctions from the
        WAVECAR file to each of the UNK files. There will be one UNK file for
        each of the kpoints in the WAVECAR file.

        Note:
            wannier90 expects the full kpoint grid instead of the symmetry-
            reduced one that VASP stores the wavefunctions on. You should run
            a nscf calculation with ISYM=0 to obtain the correct grid.

        Args:
            directory (str): directory where the UNK files are written
        """
        out_dir = Path(directory).expanduser()
        if not out_dir.exists():
            out_dir.mkdir(parents=False)
        elif not out_dir.is_dir():
            raise ValueError('invalid directory')

        N = np.prod(self.ng)
        for ik in range(self.nk):
            fname = f'UNK{ik+1:05d}.'
            if self.vasp_type.lower()[0] == 'n':
                data = np.empty((self.nb, 2, *self.ng), dtype=np.complex128)
                for ib in range(self.nb):
                    data[ib, 0, :, :, :] = \
                        np.fft.ifftn(self.fft_mesh(ik, ib, spinor=0)) * N
                    data[ib, 1, :, :, :] = \
                        np.fft.ifftn(self.fft_mesh(ik, ib, spinor=1)) * N
                Unk(ik+1, data).write_file(str(out_dir / (fname + 'NC')))
            else:
                data = np.empty((self.nb, *self.ng), dtype=np.complex128)
                for ispin in range(self.spin):
                    for ib in range(self.nb):
                        data[ib, :, :, :] = \
                            np.fft.ifftn(self.fft_mesh(ik, ib, spin=ispin)) * N
                    Unk(ik+1, data).write_file(str(out_dir / (fname + f'{ispin+1}')))


class Eigenval:
    """
    Object for reading EIGENVAL file.

    .. attribute:: filename

        string containing input filename

    .. attribute:: occu_tol

        tolerance for determining occupation in band properties

    .. attribute:: ispin

        spin polarization tag (int)

    .. attribute:: nelect

        number of electrons

    .. attribute:: nkpt

        number of kpoints

    .. attribute:: nbands

        number of bands

    .. attribute:: kpoints

        list of kpoints

    .. attribute:: kpoints_weights

        weights of each kpoint in the BZ, should sum to 1.

    .. attribute:: eigenvalues

        Eigenvalues as a dict of {(spin): np.ndarray(shape=(nkpt, nbands, 2))}.
        This representation is based on actual ordering in VASP and is meant as
        an intermediate representation to be converted into proper objects. The
        kpoint index is 0-based (unlike the 1-based indexing in VASP).
    """

    def __init__(self, filename, occu_tol=1e-8):
        """
        Reads input from filename to construct Eigenval object

        Args:
            filename (str):     filename of EIGENVAL to read in
            occu_tol (float):   tolerance for determining band gap

        Returns:
            a pymatgen.io.vasp.outputs.Eigenval object
        """

        self.filename = filename
        self.occu_tol = occu_tol

        with zopen(filename, 'r') as f:
            self.ispin = int(f.readline().split()[-1])

            # useless header information
            for _ in range(4):
                f.readline()

            self.nelect, self.nkpt, self.nbands = \
                list(map(int, f.readline().split()))

            self.kpoints = []
            self.kpoints_weights = []
            if self.ispin == 2:
                self.eigenvalues = \
                    {Spin.up: np.zeros((self.nkpt, self.nbands, 2)),
                     Spin.down: np.zeros((self.nkpt, self.nbands, 2))}
            else:
                self.eigenvalues = \
                    {Spin.up: np.zeros((self.nkpt, self.nbands, 2))}

            ikpt = -1
            for line in f:
                if re.search(r'(\s+[\-+0-9eE.]+){4}', str(line)):
                    ikpt += 1
                    kpt = list(map(float, line.split()))
                    self.kpoints.append(kpt[:-1])
                    self.kpoints_weights.append(kpt[-1])
                    for i in range(self.nbands):
                        sl = list(map(float, f.readline().split()))
                        if len(sl) == 3:
                            self.eigenvalues[Spin.up][ikpt, i, 0] = sl[1]
                            self.eigenvalues[Spin.up][ikpt, i, 1] = sl[2]
                        elif len(sl) == 5:
                            self.eigenvalues[Spin.up][ikpt, i, 0] = sl[1]
                            self.eigenvalues[Spin.up][ikpt, i, 1] = sl[3]
                            self.eigenvalues[Spin.down][ikpt, i, 0] = sl[2]
                            self.eigenvalues[Spin.down][ikpt, i, 1] = sl[4]

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
        for spin, d in self.eigenvalues.items():
            for k, val in enumerate(d):
                for (eigenval, occu) in val:
                    if occu > self.occu_tol and eigenval > vbm:
                        vbm = eigenval
                        vbm_kpoint = k
                    elif occu <= self.occu_tol and eigenval < cbm:
                        cbm = eigenval
                        cbm_kpoint = k
        return max(cbm - vbm, 0), cbm, vbm, vbm_kpoint == cbm_kpoint


class Wavederf:
    """
    Object for reading a WAVEDERF file.

    Note: This file is only produced when LOPTICS is true AND vasp has been
    recompiled after uncommenting the line that calls
    WRT_CDER_BETWEEN_STATES_FORMATTED in linear_optics.F

    .. attribute:: data

        A numpy array containing the WAVEDERF data of the form below. It should
        be noted that VASP uses 1-based indexing for bands, but this is
        converted to 0-based numpy array indexing.

        For each kpoint (in the same order as in IBZKPT), and for each pair of
        bands:

            [ #kpoint index
             [ #band 1 index
              [ #band 2 index
               [cdum_x_real, cdum_x_imag, cdum_y_real, cdum_y_imag, cdum_z_real, cdum_z_imag]
              ]
             ]
            ]

        This structure follows the file format. Numpy array methods can be used
        to fetch data in a more useful way (e.g., get matrix elements between
        wo specific bands at each kpoint, fetch x/y/z components,
        real/imaginary parts, abs/phase, etc. )

    Author: Miguel Dias Costa
    """

    def __init__(self, filename):
        """
        Args:
            filename: Name of file containing WAVEDERF.
        """
        with zopen(filename, "rt") as f:
            header = f.readline().split()
            nb_kpoints = int(header[1])
            nb_bands = int(header[2])
            data = np.zeros((nb_kpoints, nb_bands, nb_bands, 6))
            for ik in range(nb_kpoints):
                for ib1 in range(nb_bands):
                    for ib2 in range(nb_bands):
                        # each line in the file includes besides the band
                        # indexes, which are redundant, each band's energy
                        # and occupation, which are already available elsewhere,
                        # so we store only the 6 matrix elements after this 6
                        # redundant values
                        data[ik][ib1][ib2] = [float(element)
                                              for element in f.readline().split()[6:]]

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

        return self.data[:, band_i - 1, band_j - 1, :]


class Waveder:
    """
    Class for reading a WAVEDER file.
    The LOPTICS tag produces a WAVEDER file.
    The WAVEDER contains the derivative of the orbitals with respect to k.
    Author: Kamal Choudhary, NIST
    """

    def __init__(self, filename, gamma_only=False):
        """
        Args:
            filename: Name of file containing WAVEDER.
        """
        with open(filename, 'rb') as fp:
            def readData(dtype):
                """ Read records from Fortran binary file and convert to
                np.array of given dtype. """
                data = b''
                while 1:
                    prefix = np.fromfile(fp, dtype=np.int32, count=1)[0]
                    data += fp.read(abs(prefix))
                    suffix = np.fromfile(fp, dtype=np.int32, count=1)[0]
                    if abs(prefix) - abs(suffix):
                        raise RuntimeError("Read wrong amount of bytes.\n"
                                           "Expected: %d, read: %d, suffix: %d." % (prefix, len(data), suffix))
                    if prefix > 0:
                        break
                return np.frombuffer(data, dtype=dtype)

            nbands, nelect, nk, ispin = readData(np.int32)
            _ = readData(np.float)  # nodes_in_dielectric_function
            _ = readData(np.float)  # wplasmon
            if gamma_only:
                cder = readData(np.float)
            else:
                cder = readData(np.complex64)

            cder_data = cder.reshape((3, ispin, nk, nelect, nbands)).T

            self._cder_data = cder_data
            self._nkpoints = nk
            self._ispin = ispin
            self._nelect = nelect
            self._nbands = nbands

    @property
    def cder_data(self):
        """
        Returns the orbital derivative between states
        """
        return self._cder_data

    @property
    def nbands(self):
        """
        Returns the number of bands in the calculation
        """
        return self._nbands

    @property
    def nkpoints(self):
        """
        Returns the number of k-points in the calculation
        """
        return self._nkpoints

    @property
    def nelect(self):
        """
        Returns the number of electrons in the calculation
        """
        return self._nelect

    def get_orbital_derivative_between_states(self, band_i, band_j, kpoint, spin, cart_dir):
        """
        Method returning a value
        between bands band_i and band_j for k-point index, spin-channel and cartesian direction.
        Args:
            band_i (Integer): Index of band i
            band_j (Integer): Index of band j
            kpoint (Integer): Index of k-point
            spin   (Integer): Index of spin-channel (0 or 1)
            cart_dir (Integer): Index of cartesian direction (0,1,2)

        Returns:
            a float value
        """
        if band_i < 0 or band_i > self.nbands - 1 or band_j < 0 or band_j > self.nelect - 1:
            raise ValueError("Band index out of bounds")
        if kpoint > self.nkpoints:
            raise ValueError("K-point index out of bounds")
        if cart_dir > 2 or cart_dir < 0:
            raise ValueError("cart_dir index out of bounds")

        return self._cder_data[band_i, band_j, kpoint, spin, cart_dir]


class UnconvergedVASPWarning(Warning):
    """
    Warning for unconverged vasp run.
    """
    pass
