"""Classes for reading/manipulating/writing VASP output files."""

from __future__ import annotations

import datetime
import itertools
import logging
import math
import os
import re
import warnings
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass
from glob import glob
from io import StringIO
from pathlib import Path
from typing import Literal

import numpy as np
from monty.io import reverse_readfile, zopen
from monty.json import MSONable, jsanitize
from monty.os.path import zpath
from monty.re import regrep
from numpy.testing import assert_allclose

from pymatgen.core.composition import Composition
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.core.units import unitized
from pymatgen.electronic_structure.bandstructure import (
    BandStructure,
    BandStructureSymmLine,
    get_reconstructed_band_structure,
)
from pymatgen.electronic_structure.core import Magmom, Orbital, OrbitalType, Spin
from pymatgen.electronic_structure.dos import CompleteDos, Dos
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
from pymatgen.io.common import VolumetricData as BaseVolumetricData
from pymatgen.io.core import ParseError
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, Potcar
from pymatgen.io.wannier90 import Unk
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
    """
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
    err = ValueError("Error in parsing vasprun.xml")
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
                raise err
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
                raise err
    return val


def _parse_varray(elem):
    if elem.get("type") == "logical":
        m = [[i == "T" for i in v.text.split()] for v in elem]
    else:
        m = [[_vasprun_float(i) for i in v.text.split()] for v in elem]
    return m


def _parse_from_incar(filename, key):
    """Helper function to parse a parameter from the INCAR."""
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
    This function parses these values as np.nan.
    """
    try:
        return float(f)
    except ValueError as e:
        f = f.strip()
        if f == "*" * len(f):
            warnings.warn("Float overflow (*******) encountered in vasprun")
            return np.nan
        raise e


class Vasprun(MSONable):
    """
    Vastly improved cElementTree-based parser for vasprun.xml files. Uses
    iterparse to support incremental parsing of large files.
    Speedup over Dom is at least 2x for smallish files (~1Mb) to orders of
    magnitude for larger files (~10Mb).

    **VASP results**

    Attributes:
        ionic_steps (list): All ionic steps in the run as a list of {"structure": structure at end of run,
            "electronic_steps": {All electronic step data in vasprun file}, "stresses": stress matrix}.
        tdos (Dos): Total dos calculated at the end of run.
        idos (Dos): Integrated dos calculated at the end of run.
        pdos (list): List of list of PDos objects. Access as pdos[atomindex][orbitalindex].
        efermi (float): Fermi energy.
        eigenvalues (dict): Final eigenvalues as a dict of {(spin, kpoint index):[[eigenvalue, occu]]}.
            The kpoint index is 0-based (unlike the 1-based indexing in VASP).
        projected_eigenvalues (dict): Final projected eigenvalues as a dict of {spin: nd-array}.
            To access a particular value, you need to do
            Vasprun.projected_eigenvalues[spin][kpoint index][band index][atom index][orbital_index].
            The kpoint, band and atom indices are 0-based (unlike the 1-based indexing in VASP).
        projected_magnetisation (numpy.ndarray): Final projected magnetization as a numpy array with the
            shape (nkpoints, nbands, natoms, norbitals, 3). Where the last axis is the contribution in the
            3 Cartesian directions. This attribute is only set if spin-orbit coupling (LSORBIT = True) or
            non-collinear magnetism (LNONCOLLINEAR = True) is turned on in the INCAR.
        other_dielectric (dict): Dictionary, with the tag comment as key, containing other variants of
            the real and imaginary part of the dielectric constant (e.g., computed by RPA) in function of
            the energy (frequency). Optical properties (e.g. absorption coefficient) can be obtained through this.
            The data is given as a tuple of 3 values containing each of them the energy, the real part tensor,
            and the imaginary part tensor ([energies],[[real_partxx,real_partyy,real_partzz,real_partxy,
            real_partyz,real_partxz]],[[imag_partxx,imag_partyy,imag_partzz,imag_partxy, imag_partyz, imag_partxz]]).
        nionic_steps (int): The total number of ionic steps. This number is always equal to the total number
            of steps in the actual run even if ionic_step_skip is used.
        force_constants (numpy.ndarray): Force constants computed in phonon DFPT run(IBRION = 8).
            The data is a 4D numpy array of shape (natoms, natoms, 3, 3).
        normalmode_eigenvals (numpy.ndarray): Normal mode frequencies. 1D numpy array of size 3*natoms.
        normalmode_eigenvecs (numpy.ndarray): Normal mode eigen vectors. 3D numpy array of shape (3*natoms, natoms, 3).
        md_data (list): Available only for ML MD runs, i.e., INCAR with ML_LMLFF = .TRUE. md_data is a list of
            dict with the following format: [{'energy': {'e_0_energy': -525.07195568, 'e_fr_energy': -525.07195568,
            'e_wo_entrp': -525.07195568, 'kinetic': 3.17809233, 'lattice kinetic': 0.0, 'nosekinetic': 1.323e-5,
            'nosepot': 0.0, 'total': -521.89385012}, 'forces': [[0.17677989, 0.48309874, 1.85806696], ...],
            'structure': Structure object}].
        incar (Incar): Incar object for parameters specified in INCAR file.
        parameters (Incar): Incar object with parameters that vasp actually used, including all defaults.
        kpoints (Kpoints): Kpoints object for KPOINTS specified in run.
        actual_kpoints (list): List of actual kpoints, e.g., [[0.25, 0.125, 0.08333333], [-0.25, 0.125, 0.08333333],
            [0.25, 0.375, 0.08333333], ....].
        actual_kpoints_weights (list): List of kpoint weights, E.g., [0.04166667, 0.04166667, 0.04166667, 0.04166667,
            0.04166667, ....].
        atomic_symbols (list): List of atomic symbols, e.g., ["Li", "Fe", "Fe", "P", "P", "P"].
        potcar_symbols (list): List of POTCAR symbols. e.g., ["PAW_PBE Li 17Jan2003", "PAW_PBE Fe 06Sep2000", ..].

    Author: Shyue Ping Ong
    """

    def __init__(
        self,
        filename,
        ionic_step_skip=None,
        ionic_step_offset=0,
        parse_dos=True,
        parse_eigen=True,
        parse_projected_eigen=False,
        parse_potcar_file=True,
        occu_tol=1e-8,
        separate_spins=False,
        exception_on_bad_xml=True,
    ):
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
                eigenvalues and magnetization. Defaults to False. Set to True to obtain
                projected eigenvalues and magnetization. **Note that this can take an
                extreme amount of time and memory.** So use this wisely.
            parse_potcar_file (bool/str): Whether to parse the potcar file to read
                the potcar hashes for the potcar_spec attribute. Defaults to True,
                where no hashes will be determined and the potcar_spec dictionaries
                will read {"symbol": ElSymbol, "hash": None}. By Default, looks in
                the same directory as the vasprun.xml, with same extensions as
                 Vasprun.xml. If a string is provided, looks at that filepath.
            occu_tol (float): Sets the minimum tol for the determination of the
                vbm and cbm. Usually the default of 1e-8 works well enough,
                but there may be pathological cases.
            separate_spins (bool): Whether the band gap, CBM, and VBM should be
                reported for each individual spin channel. Defaults to False,
                which computes the eigenvalue band properties independent of
                the spin orientation. If True, the calculation must be spin-polarized.
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
        self.separate_spins = separate_spins
        self.exception_on_bad_xml = exception_on_bad_xml

        with zopen(filename, "rt") as f:
            if ionic_step_skip or ionic_step_offset:
                # remove parts of the xml file and parse the string
                run = f.read()
                steps = run.split("<calculation>")
                # The text before the first <calculation> is the preamble!
                preamble = steps.pop(0)
                self.nionic_steps = len(steps)
                new_steps = steps[ionic_step_offset :: int(ionic_step_skip)]
                # add the tailing information in the last step from the run
                to_parse = "<calculation>".join(new_steps)
                if steps[-1] != new_steps[-1]:
                    to_parse = f"{preamble}<calculation>{to_parse}{steps[-1].split('</calculation>')[-1]}"
                else:
                    to_parse = f"{preamble}<calculation>{to_parse}"
                self._parse(
                    StringIO(to_parse),
                    parse_dos=parse_dos,
                    parse_eigen=parse_eigen,
                    parse_projected_eigen=parse_projected_eigen,
                )
            else:
                self._parse(
                    f,
                    parse_dos=parse_dos,
                    parse_eigen=parse_eigen,
                    parse_projected_eigen=parse_projected_eigen,
                )
                self.nionic_steps = len(self.ionic_steps)

            if parse_potcar_file:
                self.update_potcar_spec(parse_potcar_file)
                self.update_charge_from_potcar(parse_potcar_file)

        if self.incar.get("ALGO") not in ["CHI", "BSE"] and not self.converged and self.parameters.get("IBRION") != 0:
            msg = f"{filename} is an unconverged VASP run.\n"
            msg += f"Electronic convergence reached: {self.converged_electronic}.\n"
            msg += f"Ionic convergence reached: {self.converged_ionic}."
            warnings.warn(msg, UnconvergedVASPWarning)

    def _parse(self, stream, parse_dos, parse_eigen, parse_projected_eigen):
        self.efermi = self.eigenvalues = self.projected_eigenvalues = self.projected_magnetisation = None
        self.dielectric_data = {}
        self.other_dielectric = {}
        self.incar = {}
        ionic_steps = []

        md_data = []
        parsed_header = False
        try:
            for _, elem in ET.iterparse(stream):
                tag = elem.tag
                if not parsed_header:
                    if tag == "generator":
                        self.generator = self._parse_params(elem)
                    elif tag == "incar":
                        self.incar = self._parse_params(elem)
                    elif tag == "kpoints":
                        if not hasattr(self, "kpoints"):
                            self.kpoints, self.actual_kpoints, self.actual_kpoints_weights = self._parse_kpoints(elem)
                    elif tag == "parameters":
                        self.parameters = self._parse_params(elem)
                    elif tag == "structure" and elem.attrib.get("name") == "initialpos":
                        self.initial_structure = self._parse_structure(elem)
                        self.final_structure = self.initial_structure
                    elif tag == "atominfo":
                        self.atomic_symbols, self.potcar_symbols = self._parse_atominfo(elem)
                        self.potcar_spec = [{"titel": p, "hash": None} for p in self.potcar_symbols]
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
                    self.projected_eigenvalues, self.projected_magnetisation = self._parse_projected_eigen(elem)
                elif tag == "dielectricfunction":
                    if (
                        "comment" not in elem.attrib
                        or elem.attrib["comment"] == "INVERSE MACROSCOPIC DIELECTRIC TENSOR (including "
                        "local field effects in RPA (Hartree))"
                    ):
                        if "density" not in self.dielectric_data:
                            self.dielectric_data["density"] = self._parse_diel(elem)
                        elif "velocity" not in self.dielectric_data:
                            # "velocity-velocity" is also named
                            # "current-current" in OUTCAR
                            self.dielectric_data["velocity"] = self._parse_diel(elem)
                        else:
                            raise NotImplementedError("This vasprun.xml has >2 unlabelled dielectric functions")
                    else:
                        comment = elem.attrib["comment"]
                        # VASP 6+ has labels for the density and current
                        # derived dielectric constants
                        if comment == "density-density":
                            self.dielectric_data["density"] = self._parse_diel(elem)
                        elif comment == "current-current":
                            self.dielectric_data["velocity"] = self._parse_diel(elem)
                        else:
                            self.other_dielectric[comment] = self._parse_diel(elem)

                elif tag == "varray" and elem.attrib.get("name") == "opticaltransitions":
                    self.optical_transition = np.array(_parse_varray(elem))
                elif tag == "structure" and elem.attrib.get("name") == "finalpos":
                    self.final_structure = self._parse_structure(elem)
                elif tag == "dynmat":
                    hessian, eigenvalues, eigenvectors = self._parse_dynmat(elem)
                    # n_atoms is not the total number of atoms, only those for which force constants were calculated
                    # https://github.com/materialsproject/pymatgen/issues/3084
                    n_atoms = len(hessian) // 3
                    hessian = np.array(hessian)
                    self.force_constants = np.zeros((n_atoms, n_atoms, 3, 3), dtype="double")
                    for ii in range(n_atoms):
                        for jj in range(n_atoms):
                            self.force_constants[ii, jj] = hessian[ii * 3 : (ii + 1) * 3, jj * 3 : (jj + 1) * 3]
                    phonon_eigenvectors = []
                    for ev in eigenvectors:
                        phonon_eigenvectors.append(np.array(ev).reshape(n_atoms, 3))
                    self.normalmode_eigenvals = np.array(eigenvalues)
                    self.normalmode_eigenvecs = np.array(phonon_eigenvectors)
                elif self.incar.get("ML_LMLFF"):
                    if tag == "structure" and elem.attrib.get("name") is None:
                        md_data.append({})
                        md_data[-1]["structure"] = self._parse_structure(elem)
                    elif tag == "varray" and elem.attrib.get("name") == "forces":
                        md_data[-1]["forces"] = _parse_varray(elem)
                    elif tag == "energy":
                        d = {i.attrib["name"]: float(i.text) for i in elem.findall("i")}
                        if "kinetic" in d:
                            md_data[-1]["energy"] = {i.attrib["name"]: float(i.text) for i in elem.findall("i")}
        except ET.ParseError as exc:
            if self.exception_on_bad_xml:
                raise exc
            warnings.warn(
                "XML is malformed. Parsing has stopped but partial data is available.",
                UserWarning,
            )
        self.ionic_steps = ionic_steps
        self.md_data = md_data
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
            imag_partxy, imag_partyz, imag_partxz]]).
        """
        return self.dielectric_data["density"]

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
            real_avg = [
                sum(self.dielectric_data["density"][1][i][0:3]) / 3
                for i in range(len(self.dielectric_data["density"][0]))
            ]
            imag_avg = [
                sum(self.dielectric_data["density"][2][i][0:3]) / 3
                for i in range(len(self.dielectric_data["density"][0]))
            ]

            def f(freq, real, imag):
                """
                The optical absorption coefficient calculated in terms of
                equation, the unit is in cm-1.
                """
                hc = 1.23984 * 1e-4  # plank constant times speed of light, in the unit of eV*cm
                return 2 * 3.14159 * np.sqrt(np.sqrt(real**2 + imag**2) - real) * np.sqrt(2) / hc * freq

            absorption_coeff = [
                f(freq, real, imag) for freq, real, imag in zip(self.dielectric_data["density"][0], real_avg, imag_avg)
            ]
        return absorption_coeff

    @property
    def converged_electronic(self):
        """
        Returns:
            bool: True if electronic step convergence has been reached in the final ionic step.
        """
        final_elec_steps = (
            self.ionic_steps[-1]["electronic_steps"] if self.incar.get("ALGO", "").lower() != "chi" else 0
        )
        # In a response function run there is no ionic steps, there is no scf step
        if self.incar.get("LEPSILON"):
            idx = 1
            to_check = {"e_wo_entrp", "e_fr_energy", "e_0_energy"}
            while set(final_elec_steps[idx]) == to_check:
                idx += 1
            return idx + 1 != self.parameters["NELM"]
        return len(final_elec_steps) < self.parameters["NELM"]

    @property
    def converged_ionic(self):
        """
        Returns:
            bool: True if ionic step convergence has been reached, i.e. that vasp
                exited before reaching the max ionic steps for a relaxation run.
        """
        nsw = self.parameters.get("NSW", 0)
        return nsw <= 1 or len(self.ionic_steps) < nsw

    @property
    def converged(self):
        """
        Returns:
            bool: True if a relaxation run is both ionically and electronically converged.
        """
        return self.converged_electronic and self.converged_ionic

    @property  # type: ignore
    @unitized("eV")
    def final_energy(self):
        """Final energy from the vasp run."""
        try:
            final_istep = self.ionic_steps[-1]
            total_energy = final_istep["e_0_energy"]

            # Addresses a bug in vasprun.xml. See https://www.vasp.at/forum/viewtopic.php?f=3&t=16942
            final_estep = final_istep["electronic_steps"][-1]
            electronic_energy_diff = final_estep["e_0_energy"] - final_estep["e_fr_energy"]
            total_energy_bugfix = np.round(electronic_energy_diff + final_istep["e_fr_energy"], 8)
            if np.abs(total_energy - total_energy_bugfix) > 1e-7:
                return total_energy_bugfix

            return total_energy
        except (IndexError, KeyError):
            warnings.warn(
                "Calculation does not have a total energy. "
                "Possibly a GW or similar kind of run. A value of "
                "infinity is returned."
            )
            return float("inf")

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
    def complete_dos_normalized(self) -> CompleteDos:
        """
        A CompleteDos object which incorporates the total DOS and all
        projected DOS. Normalized by the volume of the unit cell with
        units of states/eV/unit cell volume.
        """
        final_struct = self.final_structure
        pdoss = {final_struct[i]: pdos for i, pdos in enumerate(self.pdos)}
        return CompleteDos(self.final_structure, self.tdos, pdoss, normalize=True)

    @property
    def hubbards(self):
        """Hubbard U values used if a vasprun is a GGA+U run. {} otherwise."""
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
        raise VaspParseError("Length of U value parameters and atomic symbols are mismatched")

    @property
    def run_type(self):
        """
        Returns the run type. Currently detects GGA, metaGGA, HF, HSE, B3LYP,
        and hybrid functionals based on relevant INCAR tags. LDA is assigned if
        PAW POTCARs are used and no other functional is detected.

        Hubbard U terms and vdW corrections are detected automatically as well.
        """
        GGA_TYPES = {
            "RE": "revPBE",
            "PE": "PBE",
            "PS": "PBEsol",
            "RP": "revPBE+Padé",
            "AM": "AM05",
            "OR": "optPBE",
            "BO": "optB88",
            "MK": "optB86b",
            "--": "GGA",
        }

        METAGGA_TYPES = {
            "TPSS": "TPSS",
            "RTPSS": "revTPSS",
            "M06L": "M06-L",
            "MBJ": "modified Becke-Johnson",
            "SCAN": "SCAN",
            "R2SCAN": "R2SCAN",
            "RSCAN": "RSCAN",
            "MS0": "MadeSimple0",
            "MS1": "MadeSimple1",
            "MS2": "MadeSimple2",
        }

        IVDW_TYPES = {
            1: "DFT-D2",
            10: "DFT-D2",
            11: "DFT-D3",
            12: "DFT-D3-BJ",
            2: "TS",
            20: "TS",
            21: "TS-H",
            202: "MBD",
            4: "dDsC",
        }

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
        elif self.incar.get("METAGGA") and self.incar.get("METAGGA") not in [
            "--",
            "None",
        ]:
            incar_tag = self.incar.get("METAGGA", "").strip().upper()
            rt = METAGGA_TYPES.get(incar_tag, incar_tag)
        elif self.parameters.get("GGA"):
            incar_tag = self.parameters.get("GGA", "").strip().upper()
            rt = GGA_TYPES.get(incar_tag, incar_tag)
        elif self.potcar_symbols[0].split()[0] == "PAW":
            rt = "LDA"
        else:
            rt = "unknown"
            warnings.warn("Unknown run type!")

        if self.is_hubbard or self.parameters.get("LDAU", True):
            rt += "+U"

        if self.parameters.get("LUSE_VDW", False):
            rt += "+rVV10"
        elif self.incar.get("IVDW") in IVDW_TYPES:
            rt += "+vdW-" + IVDW_TYPES[self.incar.get("IVDW")]
        elif self.incar.get("IVDW"):
            rt += "+vdW-unknown"

        return rt

    @property
    def is_hubbard(self) -> bool:
        """True if run is a DFT+U run."""
        if len(self.hubbards) == 0:
            return False
        return sum(self.hubbards.values()) > 1e-8

    @property
    def is_spin(self) -> bool:
        """True if run is spin-polarized."""
        return self.parameters.get("ISPIN", 1) == 2

    def get_computed_entry(self, inc_structure=True, parameters=None, data=None, entry_id: str | None = None):
        """
        Returns a ComputedEntry or ComputedStructureEntry from the Vasprun.

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
            entry_id (str): Specify an entry id for the ComputedEntry. Defaults to
                "vasprun-{current datetime}"

        Returns:
            ComputedStructureEntry/ComputedEntry
        """
        if entry_id is None:
            entry_id = f"vasprun-{datetime.datetime.now()}"
        param_names = {
            "is_hubbard",
            "hubbards",
            "potcar_symbols",
            "potcar_spec",
            "run_type",
        }
        if parameters:
            param_names.update(parameters)
        params = {p: getattr(self, p) for p in param_names}
        data = {p: getattr(self, p) for p in data} if data is not None else {}

        if inc_structure:
            return ComputedStructureEntry(
                self.final_structure, self.final_energy, parameters=params, data=data, entry_id=entry_id
            )
        return ComputedEntry(
            self.final_structure.composition, self.final_energy, parameters=params, data=data, entry_id=entry_id
        )

    def get_band_structure(
        self,
        kpoints_filename: str | None = None,
        efermi: float | Literal["smart"] | None = None,
        line_mode: bool = False,
        force_hybrid_mode: bool = False,
    ) -> BandStructureSymmLine | BandStructure:
        """Get the band structure as a BandStructure object.

        Args:
            kpoints_filename: Full path of the KPOINTS file from which
                the band structure is generated.
                If none is provided, the code will try to intelligently
                determine the appropriate KPOINTS file by substituting the
                filename of the vasprun.xml with KPOINTS.
                The latter is the default behavior.
            efermi: The Fermi energy associated with the bandstructure, in eV. By
                default (None), uses the value reported by VASP in vasprun.xml. To
                manually set the Fermi energy, pass a float. Pass 'smart' to use the
                `calculate_efermi()` method, which calculates the Fermi level by first
                checking whether it lies within a small tolerance (by default 0.001 eV)
                of a band edge) If it does, the Fermi level is placed in the center of
                the bandgap. Otherwise, the value is identical to the value reported by
                VASP.
            line_mode: Force the band structure to be considered as
                a run along symmetry lines. (Default: False)
            force_hybrid_mode: Makes it possible to read in self-consistent band
                structure calculations for every type of functional. (Default: False)

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
            kpoints_filename = zpath(os.path.join(os.path.dirname(self.filename), "KPOINTS"))
        if kpoints_filename and not os.path.exists(kpoints_filename) and line_mode is True:
            raise VaspParseError("KPOINTS not found but needed to obtain band structure along symmetry lines.")

        if efermi == "smart":
            e_fermi = self.calculate_efermi()
        elif efermi is None:
            e_fermi = self.efermi
        else:
            e_fermi = efermi

        kpoint_file: Kpoints = None  # type: ignore
        if kpoints_filename and os.path.exists(kpoints_filename):
            kpoint_file = Kpoints.from_file(kpoints_filename)
        lattice_new = Lattice(self.final_structure.lattice.reciprocal_lattice.matrix)

        kpoints = [np.array(kpt) for kpt in self.actual_kpoints]

        p_eigenvals: defaultdict[Spin, list] = defaultdict(list)
        eigenvals: defaultdict[Spin, list] = defaultdict(list)

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
        if self.parameters.get("LHFCALC", False) or 0.0 in self.actual_kpoints_weights:
            hybrid_band = True

        if kpoint_file is not None and kpoint_file.style == Kpoints.supported_modes.Line_mode:
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
                        labels_dict[kpoint_file.labels[i]] = kpoint_file.kpts[i]
                # remake the data only considering line band structure k-points
                # (weight = 0.0 kpoints)
                nbands = len(eigenvals[Spin.up])
                kpoints = kpoints[start_bs_index:nkpts]
                up_eigen = [eigenvals[Spin.up][i][start_bs_index:nkpts] for i in range(nbands)]
                if self.projected_eigenvalues:
                    p_eigenvals[Spin.up] = [p_eigenvals[Spin.up][i][start_bs_index:nkpts] for i in range(nbands)]
                if self.is_spin:
                    down_eigen = [eigenvals[Spin.down][i][start_bs_index:nkpts] for i in range(nbands)]
                    eigenvals[Spin.up] = up_eigen
                    eigenvals[Spin.down] = down_eigen
                    if self.projected_eigenvalues:
                        p_eigenvals[Spin.down] = [
                            p_eigenvals[Spin.down][i][start_bs_index:nkpts] for i in range(nbands)
                        ]
                else:
                    eigenvals[Spin.up] = up_eigen
            else:
                if "" in kpoint_file.labels:
                    raise Exception(
                        "A band structure along symmetry lines "
                        "requires a label for each kpoint. "
                        "Check your KPOINTS file"
                    )
                labels_dict = dict(zip(kpoint_file.labels, kpoint_file.kpts))
                labels_dict.pop(None, None)
            return BandStructureSymmLine(
                kpoints,
                eigenvals,
                lattice_new,
                e_fermi,
                labels_dict,
                structure=self.final_structure,
                projections=p_eigenvals,
            )
        return BandStructure(
            kpoints,  # type: ignore[arg-type]
            eigenvals,  # type: ignore[arg-type]
            lattice_new,
            e_fermi,
            structure=self.final_structure,
            projections=p_eigenvals,  # type: ignore[arg-type]
        )

    @property
    def eigenvalue_band_properties(self):
        """
        Band properties from the eigenvalues as a tuple,
        (band gap, cbm, vbm, is_band_gap_direct). In the case of separate_spins=True,
        the band gap, cbm, vbm, and is_band_gap_direct are each lists of length 2,
        with index 0 representing the spin-up channel and index 1 representing
        the spin-down channel.
        """
        vbm = -float("inf")
        vbm_kpoint = None
        cbm = float("inf")
        cbm_kpoint = None
        vbm_spins = []
        vbm_spins_kpoints = []
        cbm_spins = []
        cbm_spins_kpoints = []
        if self.separate_spins and len(self.eigenvalues) != 2:
            raise ValueError("The separate_spins flag can only be True if ISPIN = 2")

        for d in self.eigenvalues.values():
            if self.separate_spins:
                vbm = -float("inf")
                cbm = float("inf")
            for k, val in enumerate(d):
                for eigenval, occu in val:
                    if occu > self.occu_tol and eigenval > vbm:
                        vbm = eigenval
                        vbm_kpoint = k
                    elif occu <= self.occu_tol and eigenval < cbm:
                        cbm = eigenval
                        cbm_kpoint = k
            if self.separate_spins:
                vbm_spins.append(vbm)
                vbm_spins_kpoints.append(vbm_kpoint)
                cbm_spins.append(cbm)
                cbm_spins_kpoints.append(cbm_kpoint)
        if self.separate_spins:
            return (
                [max(cbm_spins[0] - vbm_spins[0], 0), max(cbm_spins[1] - vbm_spins[1], 0)],
                [cbm_spins[0], cbm_spins[1]],
                [vbm_spins[0], vbm_spins[1]],
                [vbm_spins_kpoints[0] == cbm_spins_kpoints[0], vbm_spins_kpoints[1] == cbm_spins_kpoints[1]],
            )
        return max(cbm - vbm, 0), cbm, vbm, vbm_kpoint == cbm_kpoint

    def calculate_efermi(self, tol: float = 0.001):
        """
        Calculate the Fermi level using a robust algorithm.

        Sometimes VASP can put the Fermi level just inside of a band due to issues in
        the way band occupancies are handled. This algorithm tries to detect and correct
        for this bug.

        Slightly more details are provided here: https://www.vasp.at/forum/viewtopic.php?f=4&t=17981
        """
        # drop weights and set shape nbands, nkpoints
        all_eigs = np.concatenate([eigs[:, :, 0].transpose(1, 0) for eigs in self.eigenvalues.values()])

        def crosses_band(fermi):
            eigs_below = np.any(all_eigs < fermi, axis=1)
            eigs_above = np.any(all_eigs > fermi, axis=1)
            return np.any(eigs_above & eigs_below)

        def get_vbm_cbm(fermi):
            return np.max(all_eigs[all_eigs < fermi]), np.min(all_eigs[all_eigs > fermi])

        if not crosses_band(self.efermi):
            # Fermi doesn't cross a band; safe to use VASP Fermi level
            return self.efermi

        # if the Fermi level crosses a band, check if we are very close to band gap;
        # if so, then likely this is a VASP tetrahedron bug
        if not crosses_band(self.efermi + tol):
            # efermi placed slightly in the valence band
            # set Fermi level half way between valence and conduction bands
            vbm, cbm = get_vbm_cbm(self.efermi + tol)
            return (cbm + vbm) / 2

        if not crosses_band(self.efermi - tol):
            # efermi placed slightly in the conduction band
            # set Fermi level half way between valence and conduction bands
            vbm, cbm = get_vbm_cbm(self.efermi - tol)
            return (cbm + vbm) / 2

        # it is actually a metal
        return self.efermi

    def get_potcars(self, path: str | Path) -> Potcar | None:
        """
        Returns the POTCAR from the specified path.

        Args:
            path (str): The path to search for POTCARs.

        Returns:
            Potcar | None: The POTCAR from the specified path.
        """

        def get_potcar_in_path(p):
            for fn in os.listdir(os.path.abspath(p)):
                if fn.startswith("POTCAR") and ".spec" not in fn:
                    pc = Potcar.from_file(os.path.join(p, fn))
                    if {d.header for d in pc} == set(self.potcar_symbols):
                        return pc
            warnings.warn(f"No POTCAR file with matching TITEL fields was found in {os.path.abspath(p)}")
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

        Returns:
            Trajectory: from pymatgen.core.trajectory
        """
        # required due to circular imports
        # TODO: fix pymatgen.core.trajectory so it does not load from io.vasp(!)
        from pymatgen.core.trajectory import Trajectory

        structs = []
        for step in self.ionic_steps:
            struct = step["structure"].copy()
            struct.add_site_property("forces", step["forces"])
            structs.append(struct)
        return Trajectory.from_structures(structs, constant_lattice=False)

    def update_potcar_spec(self, path):
        """
        Args:
            path: Path to search for POTCARs

        Returns:
            Potcar spec from path.
        """
        if potcar := self.get_potcars(path):
            self.potcar_spec = [
                {"titel": sym, "hash": ps.get_potcar_hash()}
                for sym in self.potcar_symbols
                for ps in potcar
                if ps.symbol == sym.split()[1]
            ]

    def update_charge_from_potcar(self, path):
        """
        Sets the charge of a structure based on the POTCARs found.

        :param path: Path to search for POTCARs
        """
        potcar = self.get_potcars(path)

        if potcar and self.incar.get("ALGO", "") not in ["GW0", "G0W0", "GW", "BSE"]:
            nelect = self.parameters["NELECT"]
            if len(potcar) == len(self.initial_structure.composition.element_composition):
                potcar_nelect = sum(
                    self.initial_structure.composition.element_composition[ps.element] * ps.ZVAL for ps in potcar
                )
            else:
                nums = [len(list(g)) for _, g in itertools.groupby(self.atomic_symbols)]
                potcar_nelect = sum(ps.ZVAL * num for ps, num in zip(potcar, nums))
            charge = potcar_nelect - nelect

            for s in self.structures:
                s._charge = charge
            if hasattr(self, "initial_structure"):
                self.initial_structure._charge = charge
            if hasattr(self, "final_structure"):
                self.final_structure._charge = charge

    def as_dict(self):
        """JSON-serializable dict representation."""
        dct = {
            "vasp_version": self.vasp_version,
            "has_vasp_completed": self.converged,
            "nsites": len(self.final_structure),
        }
        comp = self.final_structure.composition
        dct["unit_cell_formula"] = comp.as_dict()
        dct["reduced_cell_formula"] = Composition(comp.reduced_formula).as_dict()
        dct["pretty_formula"] = comp.reduced_formula
        symbols = [s.split()[1] for s in self.potcar_symbols]
        symbols = [re.split(r"_", s)[0] for s in symbols]
        dct["is_hubbard"] = self.is_hubbard
        dct["hubbards"] = self.hubbards

        unique_symbols = sorted(set(self.atomic_symbols))
        dct["elements"] = unique_symbols
        dct["nelements"] = len(unique_symbols)

        dct["run_type"] = self.run_type

        vin = {
            "incar": dict(self.incar.items()),
            "crystal": self.initial_structure.as_dict(),
            "kpoints": self.kpoints.as_dict(),
        }
        actual_kpts = [
            {
                "abc": list(self.actual_kpoints[i]),
                "weight": self.actual_kpoints_weights[i],
            }
            for i in range(len(self.actual_kpoints))
        ]
        vin["kpoints"]["actual_points"] = actual_kpts
        vin["nkpoints"] = len(actual_kpts)
        vin["potcar"] = [s.split(" ")[1] for s in self.potcar_symbols]
        vin["potcar_spec"] = self.potcar_spec
        vin["potcar_type"] = [s.split(" ")[0] for s in self.potcar_symbols]
        vin["parameters"] = dict(self.parameters.items())
        vin["lattice_rec"] = self.final_structure.lattice.reciprocal_lattice.as_dict()
        dct["input"] = vin

        n_sites = len(self.final_structure)

        try:
            vout = {
                "ionic_steps": self.ionic_steps,
                "final_energy": self.final_energy,
                "final_energy_per_atom": self.final_energy / n_sites,
                "crystal": self.final_structure.as_dict(),
                "efermi": self.efermi,
            }
        except (ArithmeticError, TypeError):
            vout = {
                "ionic_steps": self.ionic_steps,
                "final_energy": self.final_energy,
                "final_energy_per_atom": None,
                "crystal": self.final_structure.as_dict(),
                "efermi": self.efermi,
            }

        if self.eigenvalues:
            eigen = {str(spin): v.tolist() for spin, v in self.eigenvalues.items()}
            vout["eigenvalues"] = eigen
            (gap, cbm, vbm, is_direct) = self.eigenvalue_band_properties
            vout.update({"bandgap": gap, "cbm": cbm, "vbm": vbm, "is_gap_direct": is_direct})

            if self.projected_eigenvalues:
                vout["projected_eigenvalues"] = {
                    str(spin): v.tolist() for spin, v in self.projected_eigenvalues.items()
                }

            if self.projected_magnetisation is not None:
                vout["projected_magnetisation"] = self.projected_magnetisation.tolist()

        vout["epsilon_static"] = self.epsilon_static
        vout["epsilon_static_wolfe"] = self.epsilon_static_wolfe
        vout["epsilon_ionic"] = self.epsilon_ionic
        dct["output"] = vout
        return jsanitize(dct, strict=True)

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
                try:
                    if c.tag == "i":
                        params[name] = _parse_parameters(ptype, val)
                    else:
                        params[name] = _parse_v_parameters(ptype, val, self.filename, name)
                except Exception as exc:
                    if name == "RANDOM_SEED":
                        # Handles the case where RANDOM SEED > 99999, which results in *****
                        params[name] = None
                    else:
                        raise exc
        elem.clear()
        return Incar(params)

    @staticmethod
    def _parse_atominfo(elem):
        for a in elem.findall("array"):
            if a.attrib["name"] == "atoms":
                atomic_symbols = [rc.find("c").text.strip() for rc in a.find("set")]
            elif a.attrib["name"] == "atomtypes":
                potcar_symbols = [rc.findall("c")[4].text.strip() for rc in a.find("set")]

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
        k.style = Kpoints.supported_modes.from_str(e.attrib["param"] if "param" in e.attrib else "Reciprocal")
        for v in e.findall("v"):
            name = v.attrib.get("name")
            tokens = v.text.split()
            if name == "divisions":
                k.kpts = [[int(i) for i in tokens]]
            elif name == "usershift":
                k.kpts_shift = [float(i) for i in tokens]
            elif name in {"genvec1", "genvec2", "genvec3", "shift"}:
                setattr(k, name, [float(i) for i in tokens])
        for va in elem.findall("varray"):
            name = va.attrib["name"]
            if name == "kpointlist":
                actual_kpoints = _parse_varray(va)
            elif name == "weights":
                weights = [i[0] for i in _parse_varray(va)]
        elem.clear()
        if k.style == Kpoints.supported_modes.Reciprocal:
            k = Kpoints(
                comment="Kpoints from vasprun.xml",
                style=Kpoints.supported_modes.Reciprocal,
                num_kpts=len(k.kpts),
                kpts=actual_kpoints,
                kpts_weights=weights,
            )
        return k, actual_kpoints, weights

    def _parse_structure(self, elem):
        latt = _parse_varray(elem.find("crystal").find("varray"))
        pos = _parse_varray(elem.find("varray"))
        struct = Structure(latt, self.atomic_symbols, pos)
        sdyn = elem.find("varray/[@name='selective']")
        if sdyn:
            struct.add_site_property("selective_dynamics", _parse_varray(sdyn))
        return struct

    @staticmethod
    def _parse_diel(elem):
        imag = [
            [_vasprun_float(line) for line in r.text.split()]
            for r in elem.find("imag").find("array").find("set").findall("r")
        ]
        real = [
            [_vasprun_float(line) for line in r.text.split()]
            for r in elem.find("real").find("array").find("set").findall("r")
        ]
        elem.clear()
        return [e[0] for e in imag], [e[1:] for e in real], [e[1:] for e in imag]

    @staticmethod
    def _parse_optical_transition(elem):
        for va in elem.findall("varray"):
            if va.attrib.get("name") == "opticaltransitions":
                # opticaltransitions array contains oscillator strength and probability of transition
                oscillator_strength = np.array(_parse_varray(va))[0:]
                probability_transition = np.array(_parse_varray(va))[0:, 1]
        return oscillator_strength, probability_transition

    def _parse_chemical_shielding_calculation(self, elem):
        calculation = []
        istep = {}
        try:
            struct = self._parse_structure(elem.find("structure"))
        except AttributeError:  # not all calculations have a structure
            struct = None
        for va in elem.findall("varray"):
            istep[va.attrib["name"]] = _parse_varray(va)
        istep["structure"] = struct
        istep["electronic_steps"] = []
        calculation.append(istep)
        for scstep in elem.findall("scstep"):
            try:
                d = {i.attrib["name"]: _vasprun_float(i.text) for i in scstep.find("energy").findall("i")}
                cur_ene = d["e_fr_energy"]
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
            istep = {i.attrib["name"]: float(i.text) for i in elem.find("energy").findall("i")}
        except AttributeError:  # not all calculations have an energy
            istep = {}
        esteps = []
        for scstep in elem.findall("scstep"):
            try:
                d = {i.attrib["name"]: _vasprun_float(i.text) for i in scstep.find("energy").findall("i")}
                esteps.append(d)
            except AttributeError:  # not all calculations have an energy
                pass
        try:
            struct = self._parse_structure(elem.find("structure"))
        except AttributeError:  # not all calculations have a structure
            struct = None
        for va in elem.findall("varray"):
            istep[va.attrib["name"]] = _parse_varray(va)
        istep["electronic_steps"] = esteps
        istep["structure"] = struct
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
            lm = any("x" in s for s in orbs)
            for s in partial.find("array").find("set").findall("set"):
                pdos = defaultdict(dict)

                for ss in s.findall("set"):
                    spin = Spin.up if ss.attrib["comment"] == "spin 1" else Spin.down
                    data = np.array(_parse_varray(ss))
                    nrow, ncol = data.shape
                    for j in range(1, ncol):
                        orb = Orbital(j - 1) if lm else OrbitalType(j - 1)
                        pdos[orb][spin] = data[:, j]
                pdoss.append(pdos)
        elem.clear()
        return (
            Dos(efermi, energies, tdensities),
            Dos(efermi, energies, idensities),
            pdoss,
        )

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
            for ss in s.findall("set"):
                dk = []
                for sss in ss.findall("set"):
                    db = _parse_varray(sss)
                    dk.append(db)
                proj_eigen[spin].append(dk)
        proj_eigen = {spin: np.array(v) for spin, v in proj_eigen.items()}

        if len(proj_eigen) > 2:
            # non-collinear magentism (also spin-orbit coupling) enabled, last three
            # "spin channels" are the projected magnetization of the orbitals in the
            # x, y, and z Cartesian coordinates
            proj_mag = np.stack([proj_eigen.pop(i) for i in range(2, 5)], axis=-1)
            proj_eigen = {Spin.up: proj_eigen[1]}
        else:
            proj_eigen = {Spin.up if k == 1 else Spin.down: v for k, v in proj_eigen.items()}
            proj_mag = None

        elem.clear()
        return proj_eigen, proj_mag

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

    def __init__(
        self,
        filename: str,
        parse_projected_eigen: bool | str = False,
        parse_potcar_file: bool | str = False,
        occu_tol: float = 1e-8,
        separate_spins: bool = False,
    ):
        """
        Args:
            filename: Filename to parse
            parse_projected_eigen: Whether to parse the projected
                eigenvalues. Defaults to False. Set to True to obtain projected
                eigenvalues. **Note that this can take an extreme amount of time
                and memory.** So use this wisely.
            parse_potcar_file: Whether to parse the potcar file to read
                the potcar hashes for the potcar_spec attribute. Defaults to True,
                where no hashes will be determined and the potcar_spec dictionaries
                will read {"symbol": ElSymbol, "hash": None}. By Default, looks in
                the same directory as the vasprun.xml, with same extensions as
                 Vasprun.xml. If a string is provided, looks at that filepath.
            occu_tol: Sets the minimum tol for the determination of the
                vbm and cbm. Usually the default of 1e-8 works well enough,
                but there may be pathological cases.
            separate_spins (bool): Whether the band gap, CBM, and VBM should be
                reported for each individual spin channel. Defaults to False,
                which computes the eigenvalue band properties independent of
                the spin orientation. If True, the calculation must be spin-polarized.
        """
        self.filename = filename
        self.occu_tol = occu_tol
        self.separate_spins = separate_spins

        with zopen(filename, "rt") as f:
            self.efermi = None
            parsed_header = False
            self.eigenvalues = self.projected_eigenvalues = None
            for _, elem in ET.iterparse(f):
                tag = elem.tag
                if not parsed_header:
                    if tag == "generator":
                        self.generator = self._parse_params(elem)
                    elif tag == "incar":
                        self.incar = self._parse_params(elem)
                    elif tag == "kpoints":
                        (
                            self.kpoints,
                            self.actual_kpoints,
                            self.actual_kpoints_weights,
                        ) = self._parse_kpoints(elem)
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
                    self.projected_eigenvalues, self.projected_magnetisation = self._parse_projected_eigen(elem)
                elif tag == "structure" and elem.attrib.get("name") == "finalpos":
                    self.final_structure = self._parse_structure(elem)
        self.vasp_version = self.generator["version"]
        if parse_potcar_file:
            self.update_potcar_spec(parse_potcar_file)

    def as_dict(self):
        """JSON-serializable dict representation."""
        dct = {
            "vasp_version": self.vasp_version,
            "has_vasp_completed": True,
            "nsites": len(self.final_structure),
        }
        comp = self.final_structure.composition
        dct["unit_cell_formula"] = comp.as_dict()
        dct["reduced_cell_formula"] = Composition(comp.reduced_formula).as_dict()
        dct["pretty_formula"] = comp.reduced_formula
        dct["is_hubbard"] = self.is_hubbard
        dct["hubbards"] = self.hubbards

        unique_symbols = sorted(set(self.atomic_symbols))
        dct["elements"] = unique_symbols
        dct["nelements"] = len(unique_symbols)

        dct["run_type"] = self.run_type

        vin = {
            "incar": dict(self.incar),
            "crystal": self.final_structure.as_dict(),
            "kpoints": self.kpoints.as_dict(),
        }
        actual_kpts = [
            {
                "abc": list(self.actual_kpoints[i]),
                "weight": self.actual_kpoints_weights[i],
            }
            for i in range(len(self.actual_kpoints))
        ]
        vin["kpoints"]["actual_points"] = actual_kpts
        vin["potcar"] = [s.split(" ")[1] for s in self.potcar_symbols]
        vin["potcar_spec"] = self.potcar_spec
        vin["potcar_type"] = [s.split(" ")[0] for s in self.potcar_symbols]
        vin["parameters"] = dict(self.parameters)
        vin["lattice_rec"] = self.final_structure.lattice.reciprocal_lattice.as_dict()
        dct["input"] = vin

        vout = {"crystal": self.final_structure.as_dict(), "efermi": self.efermi}

        if self.eigenvalues:
            eigen = defaultdict(dict)
            for spin, values in self.eigenvalues.items():
                for i, v in enumerate(values):
                    eigen[i][str(spin)] = v
            vout["eigenvalues"] = eigen
            (gap, cbm, vbm, is_direct) = self.eigenvalue_band_properties
            vout.update({"bandgap": gap, "cbm": cbm, "vbm": vbm, "is_gap_direct": is_direct})

            if self.projected_eigenvalues:
                peigen = [{} for _ in eigen]
                for spin, v in self.projected_eigenvalues.items():
                    for kpoint_index, vv in enumerate(v):
                        if str(spin) not in peigen[kpoint_index]:
                            peigen[kpoint_index][str(spin)] = vv
                vout["projected_eigenvalues"] = peigen

        dct["output"] = vout
        return jsanitize(dct, strict=True)


class Outcar:
    """
    Parser for data in OUTCAR that is not available in Vasprun.xml.

    Note, this class works a bit differently than most of the other
    VaspObjects, since the OUTCAR can be very different depending on which
    "type of run" performed.

    Creating the OUTCAR class with a filename reads "regular parameters" that
    are always present.

    Attributes:
        magnetization (tuple): Magnetization on each ion as a tuple of dict, e.g.,
            ({"d": 0.0, "p": 0.003, "s": 0.002, "tot": 0.005}, ... )
        chemical_shielding (dict): Chemical shielding on each ion as a dictionary with core and valence contributions.
        unsym_cs_tensor (list): Unsymmetrized chemical shielding tensor matrixes on each ion as a list.
            e.g., [[[sigma11, sigma12, sigma13], [sigma21, sigma22, sigma23], [sigma31, sigma32, sigma33]], ...]
        cs_g0_contribution (numpy.ndarray): G=0 contribution to chemical shielding. 2D rank 3 matrix.
        cs_core_contribution (dict): Core contribution to chemical shielding. dict. e.g.,
            {'Mg': -412.8, 'C': -200.5, 'O': -271.1}
        efg (tuple): Electric Field Gradient (EFG) tensor on each ion as a tuple of dict, e.g.,
            ({"cq": 0.1, "eta", 0.2, "nuclear_quadrupole_moment": 0.3}, {"cq": 0.7, "eta", 0.8,
            "nuclear_quadrupole_moment": 0.9}, ...)
        charge (tuple): Charge on each ion as a tuple of dict, e.g.,
            ({"p": 0.154, "s": 0.078, "d": 0.0, "tot": 0.232}, ...)
        is_stopped (bool): True if OUTCAR is from a stopped run (using STOPCAR, see VASP Manual).
        run_stats (dict): Various useful run stats as a dict including "System time (sec)", "Total CPU time used (sec)",
            "Elapsed time (sec)", "Maximum memory used (kb)", "Average memory used (kb)", "User time (sec)", "cores".
        elastic_tensor (numpy.ndarray): Total elastic moduli (Kbar) is given in a 6x6 array matrix.
        drift (numpy.ndarray): Total drift for each step in eV/Atom.
        ngf (tuple): Dimensions for the Augmentation grid.
        sampling_radii (numpy.ndarray): Size of the sampling radii in VASP for the test charges for the electrostatic
            potential at each atom. Total array size is the number of elements present in the calculation.
        electrostatic_potential (numpy.ndarray): Average electrostatic potential at each atomic position in order of
            the atoms in POSCAR.
        final_energy_contribs (dict): Individual contributions to the total final energy as a dictionary.
            Include contributions from keys, e.g.:
            {'DENC': -505778.5184347, 'EATOM': 15561.06492564, 'EBANDS': -804.53201231, 'EENTRO': -0.08932659,
            'EXHF': 0.0, 'Ediel_sol': 0.0, 'PAW double counting': 664.6726974100002, 'PSCENC': 742.48691646,
            'TEWEN': 489742.86847338, 'XCENC': -169.64189814}
        efermi (float): Fermi energy.
        filename (str): Filename.
        final_energy (float): Final energy after extrapolation of sigma back to 0, i.e. energy(sigma->0).
        final_energy_wo_entrp (float): Final energy before extrapolation of sigma, i.e. energy without entropy.
        final_fr_energy (float): Final "free energy", i.e. free energy TOTEN.
        has_onsite_density_matrices (bool): Boolean for if onsite density matrices have been set.
        lcalcpol (bool): If LCALCPOL has been set.
        lepsilon (bool): If LEPSILON has been set.
        nelect (float): Returns the number of electrons in the calculation.
        spin (bool): If spin-polarization was enabled via ISPIN.
        total_mag (float): Total magnetization (in terms of the number of unpaired electrons).

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

        # Assume a compilation with parallelization enabled.
        # Will be checked later.
        # If VASP is compiled in serial, the OUTCAR is written slightly differently.
        serial_compilation = False

        # data from end of OUTCAR
        charge = []
        mag_x = []
        mag_y = []
        mag_z = []
        header = []
        run_stats = {}
        total_mag = nelect = efermi = e_fr_energy = e_wo_entrp = e0 = None

        time_patt = re.compile(r"\((sec|kb)\)")
        efermi_patt = re.compile(r"E-fermi\s*:\s*(\S+)")
        nelect_patt = re.compile(r"number of electron\s+(\S+)\s+magnetization")
        mag_patt = re.compile(r"number of electron\s+\S+\s+magnetization\s+(\S+)")
        e_fr_energy_pattern = re.compile(r"free  energy   TOTEN\s+=\s+([\d\-\.]+)")
        e_wo_entrp_pattern = re.compile(r"energy  without entropy\s*=\s+([\d\-\.]+)")
        e0_pattern = re.compile(r"energy\(sigma->0\)\s*=\s+([\d\-\.]+)")

        all_lines = []
        for line in reverse_readfile(self.filename):
            clean = line.strip()
            all_lines.append(clean)
            if clean.find("soft stop encountered!  aborting job") != -1:
                self.is_stopped = True
            else:
                if time_patt.search(line):
                    tok = line.strip().split(":")
                    try:
                        # try-catch because VASP 6.2.0 may print
                        # Average memory used (kb):          N/A
                        # which cannot be parsed as float
                        run_stats[tok[0].strip()] = float(tok[1].strip())
                    except ValueError:
                        run_stats[tok[0].strip()] = None
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

                if e_fr_energy is None:
                    m = e_fr_energy_pattern.search(clean)
                    if m:
                        e_fr_energy = float(m.group(1))
                if e_wo_entrp is None:
                    m = e_wo_entrp_pattern.search(clean)
                    if m:
                        e_wo_entrp = float(m.group(1))
                if e0 is None:
                    m = e0_pattern.search(clean)
                    if m:
                        e0 = float(m.group(1))
            if all([nelect, total_mag is not None, efermi is not None, run_stats]):
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
                        tokens = [float(i) for i in re.findall(r"[\d\.\-]+", clean)]
                        tokens.pop(0)
                        if read_charge:
                            charge.append(dict(zip(header, tokens)))
                        elif read_mag_x:
                            mag_x.append(dict(zip(header, tokens)))
                        elif read_mag_y:
                            mag_y.append(dict(zip(header, tokens)))
                        elif read_mag_z:
                            mag_z.append(dict(zip(header, tokens)))
                    elif clean.startswith("tot"):
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
                read_charge, read_mag_x, read_mag_y, read_mag_z = (
                    False,
                    False,
                    False,
                    False,
                )

        # merge x, y and z components of magmoms if present (SOC calculation)
        if mag_y and mag_z:
            # TODO: detect spin axis
            mag = []
            for idx in range(len(mag_x)):
                mag.append({key: Magmom([mag_x[idx][key], mag_y[idx][key], mag_z[idx][key]]) for key in mag_x[0]})
        else:
            mag = mag_x

        # data from beginning of OUTCAR
        run_stats["cores"] = None
        with zopen(filename, "rt") as f:
            for line in f:
                if "serial" in line:
                    # activate the serial parallelization
                    run_stats["cores"] = 1
                    serial_compilation = True
                    break
                if "running" in line:
                    if line.split()[1] == "on":
                        run_stats["cores"] = int(line.split()[2])
                    else:
                        run_stats["cores"] = int(line.split()[1])
                    break

        self.run_stats = run_stats
        self.magnetization = tuple(mag)
        self.charge = tuple(charge)
        self.efermi = efermi
        self.nelect = nelect
        self.total_mag = total_mag
        self.final_energy = e0
        self.final_energy_wo_entrp = e_wo_entrp
        self.final_fr_energy = e_fr_energy
        self.data = {}

        # Read "total number of plane waves", NPLWV:
        self.read_pattern(
            {"nplwv": r"total plane-waves  NPLWV =\s+(\*{6}|\d+)"},
            terminate_on_match=True,
        )
        try:
            self.data["nplwv"] = [[int(self.data["nplwv"][0][0])]]
        except ValueError:
            self.data["nplwv"] = [[None]]

        nplwvs_at_kpoints = [
            n
            for [n] in self.read_table_pattern(
                r"\n{3}-{104}\n{3}",
                r".+plane waves:\s+(\*{6,}|\d+)",
                r"maximum number of plane-waves"
                if serial_compilation
                else r"maximum and minimum number of plane-waves",
                last_one_only=False,
                first_one_only=True,
            )
        ]
        self.data["nplwvs_at_kpoints"] = [None for n in nplwvs_at_kpoints]
        for n, nplwv in enumerate(nplwvs_at_kpoints):
            try:
                self.data["nplwvs_at_kpoints"][n] = int(nplwv)
            except ValueError:
                pass

        # Read the drift:
        self.read_pattern(
            {"drift": r"total drift:\s+([\.\-\d]+)\s+([\.\-\d]+)\s+([\.\-\d]+)"},
            terminate_on_match=False,
            postprocess=float,
        )
        self.drift = self.data.get("drift", [])

        # Check if calculation is spin polarized
        self.spin = False
        self.read_pattern({"spin": "ISPIN  =      2"})
        if self.data.get("spin", []):
            self.spin = True

        # Check if calculation is noncollinear
        self.noncollinear = False
        self.read_pattern({"noncollinear": "LNONCOLLINEAR =      T"})
        if self.data.get("noncollinear", []):
            self.noncollinear = False

        # Check if the calculation type is DFPT
        self.dfpt = False
        self.read_pattern(
            {"ibrion": r"IBRION =\s+([\-\d]+)"},
            terminate_on_match=True,
            postprocess=int,
        )
        if self.data.get("ibrion", [[0]])[0][0] > 6:
            self.dfpt = True
            self.read_internal_strain_tensor()

        # Check to see if LEPSILON is true and read piezo data if so
        self.lepsilon = False
        self.read_pattern({"epsilon": "LEPSILON=     T"})
        if self.data.get("epsilon", []):
            self.lepsilon = True
            self.read_lepsilon()
            # only read ionic contribution if DFPT is turned on
            if self.dfpt:
                self.read_lepsilon_ionic()

        # Check to see if LCALCPOL is true and read polarization data if so
        self.lcalcpol = False
        self.read_pattern({"calcpol": "LCALCPOL   =     T"})
        if self.data.get("calcpol", []):
            self.lcalcpol = True
            self.read_lcalcpol()
            self.read_pseudo_zval()

        # Read electrostatic potential
        self.electrostatic_potential = self.ngf = self.sampling_radii = None
        self.read_pattern({"electrostatic": r"average \(electrostatic\) potential at core"})
        if self.data.get("electrostatic", []):
            self.read_electrostatic_potential()

        self.nmr_cs = False
        self.read_pattern({"nmr_cs": r"LCHIMAG   =     (T)"})
        if self.data.get("nmr_cs"):
            self.nmr_cs = True
            self.read_chemical_shielding()
            self.read_cs_g0_contribution()
            self.read_cs_core_contribution()
            self.read_cs_raw_symmetrized_tensors()

        self.nmr_efg = False
        self.read_pattern({"nmr_efg": r"NMR quadrupolar parameters"})
        if self.data.get("nmr_efg"):
            self.nmr_efg = True
            self.read_nmr_efg()
            self.read_nmr_efg_tensor()

        self.has_onsite_density_matrices = False
        self.read_pattern(
            {"has_onsite_density_matrices": r"onsite density matrix"},
            terminate_on_match=True,
        )
        if "has_onsite_density_matrices" in self.data:
            self.has_onsite_density_matrices = True
            self.read_onsite_density_matrices()

        # Store the individual contributions to the final total energy
        final_energy_contribs = {}
        for key in [
            "PSCENC",
            "TEWEN",
            "DENC",
            "EXHF",
            "XCENC",
            "PAW double counting",
            "EENTRO",
            "EBANDS",
            "EATOM",
            "Ediel_sol",
        ]:
            if key == "PAW double counting":
                self.read_pattern({key: rf"{key}\s+=\s+([\.\-\d]+)\s+([\.\-\d]+)"})
            else:
                self.read_pattern({key: rf"{key}\s+=\s+([\d\-\.]+)"})
            if not self.data[key]:
                continue
            final_energy_contribs[key] = sum(float(f) for f in self.data[key][-1])
        self.final_energy_contribs = final_energy_contribs

    def read_pattern(self, patterns, reverse=False, terminate_on_match=False, postprocess=str):
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
        matches = regrep(
            self.filename,
            patterns,
            reverse=reverse,
            terminate_on_match=terminate_on_match,
            postprocess=postprocess,
        )
        for k in patterns:
            self.data[k] = [i[0] for i in matches.get(k, [])]

    def read_table_pattern(
        self,
        header_pattern,
        row_pattern,
        footer_pattern,
        postprocess=str,
        attribute_name=None,
        last_one_only=True,
        first_one_only=False,
    ):
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
                be returned. Default to be True. Incompatible with first_one_only.
            first_one_only (bool): Only the first occurrence of the table will be
                parsed and the parsing procedure will stop. The enclosing list
                will be removed. i.e. Only a single table will be returned.
                Incompatible with last_one_only.

        Returns:
            List of tables. 1) A table is a list of rows. 2) A row if either a list of
            attribute values in case the capturing group is defined without name in
            row_pattern, or a dict in case that named capturing groups are defined by
            row_pattern.
        """
        if last_one_only and first_one_only:
            raise ValueError("last_one_only and first_one_only options are incompatible")

        with zopen(self.filename, "rt") as f:
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
            if first_one_only:
                break
        retained_data = tables[-1] if last_one_only or first_one_only else tables
        if attribute_name is not None:
            self.data[attribute_name] = retained_data
        return retained_data

    def read_electrostatic_potential(self):
        """Parses the eletrostatic potential for the last ionic step."""
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
            line: line to parse.

        Returns:
            list[float]: numbers if found, empty ist if not
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
        dielectric_pattern = (
            r"frequency dependent\s+IMAGINARY "
            r"DIELECTRIC FUNCTION \(independent particle, "
            r"no local field effects\)(\sdensity-density)*$"
        )
        row_pattern = r"\s+".join([r"([\.\-\d]+)"] * 3)
        plasma_frequencies = defaultdict(list)
        read_plasma = False
        read_dielectric = False
        energies = []
        data = {"REAL": [], "IMAGINARY": []}
        count = 0
        component = "IMAGINARY"
        with zopen(self.filename, "rt") as file:
            for line in file:
                line = line.strip()
                if re.match(plasma_pattern, line):
                    read_plasma = "intraband" if "intraband" in line else "interband"
                elif re.match(dielectric_pattern, line):
                    read_plasma = False
                    read_dielectric = True
                    row_pattern = r"\s+".join([r"([\.\-\d]+)"] * 7)

                if read_plasma and re.match(row_pattern, line):
                    plasma_frequencies[read_plasma].append([float(t) for t in line.strip().split()])
                elif read_plasma and Outcar._parse_sci_notation(line):
                    plasma_frequencies[read_plasma].append(Outcar._parse_sci_notation(line))
                elif read_dielectric:
                    tokens = None
                    if re.match(row_pattern, line.strip()):
                        tokens = line.strip().split()
                    elif Outcar._parse_sci_notation(line.strip()):
                        tokens = Outcar._parse_sci_notation(line.strip())
                    elif re.match(r"\s*-+\s*", line):
                        count += 1

                    if tokens:
                        if component == "IMAGINARY":
                            energies.append(float(tokens[0]))
                        xx, yy, zz, xy, yz, xz = (float(t) for t in tokens[1:])
                        matrix = [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]
                        data[component].append(matrix)

                    if count == 2:
                        component = "REAL"
                    elif count == 3:
                        break

        self.plasma_frequencies = {k: np.array(v[:3]) for k, v in plasma_frequencies.items()}
        self.dielectric_energies = np.array(energies)
        self.dielectric_tensor_function = np.array(data["REAL"]) + 1j * np.array(data["IMAGINARY"])

    def read_chemical_shielding(self):
        """
        Parse the NMR chemical shieldings data. Only the second part "absolute, valence and core"
        will be parsed. And only the three right most field (ISO_SHIELDING, SPAN, SKEW) will be retrieved.

        Returns:
            List of chemical shieldings in the order of atoms from the OUTCAR. Maryland notation is adopted.
        """
        header_pattern = (
            r"\s+CSA tensor \(J\. Mason, Solid State Nucl\. Magn\. Reson\. 2, "
            r"285 \(1993\)\)\s+"
            r"\s+-{50,}\s+"
            r"\s+EXCLUDING G=0 CONTRIBUTION\s+INCLUDING G=0 CONTRIBUTION\s+"
            r"\s+-{20,}\s+-{20,}\s+"
            r"\s+ATOM\s+ISO_SHIFT\s+SPAN\s+SKEW\s+ISO_SHIFT\s+SPAN\s+SKEW\s+"
            r"-{50,}\s*$"
        )
        first_part_pattern = r"\s+\(absolute, valence only\)\s+$"
        swallon_valence_body_pattern = r".+?\(absolute, valence and core\)\s+$"
        row_pattern = r"\d+(?:\s+[-]?\d+\.\d+){3}\s+" + r"\s+".join([r"([-]?\d+\.\d+)"] * 3)
        footer_pattern = r"-{50,}\s*$"
        h1 = header_pattern + first_part_pattern
        cs_valence_only = self.read_table_pattern(
            h1, row_pattern, footer_pattern, postprocess=float, last_one_only=True
        )
        h2 = header_pattern + swallon_valence_body_pattern
        cs_valence_and_core = self.read_table_pattern(
            h2, row_pattern, footer_pattern, postprocess=float, last_one_only=True
        )
        all_cs = {}
        for name, cs_table in [
            ["valence_only", cs_valence_only],
            ["valence_and_core", cs_valence_and_core],
        ]:
            all_cs[name] = cs_table
        self.data["chemical_shielding"] = all_cs

    def read_cs_g0_contribution(self):
        """
        Parse the  G0 contribution of NMR chemical shielding.

        Returns:
            G0 contribution matrix as list of list.
        """
        header_pattern = (
            r"^\s+G\=0 CONTRIBUTION TO CHEMICAL SHIFT \(field along BDIR\)\s+$\n"
            r"^\s+-{50,}$\n"
            r"^\s+BDIR\s+X\s+Y\s+Z\s*$\n"
            r"^\s+-{50,}\s*$\n"
        )
        row_pattern = r"(?:\d+)\s+" + r"\s+".join([r"([-]?\d+\.\d+)"] * 3)
        footer_pattern = r"\s+-{50,}\s*$"
        self.read_table_pattern(
            header_pattern,
            row_pattern,
            footer_pattern,
            postprocess=float,
            last_one_only=True,
            attribute_name="cs_g0_contribution",
        )

    def read_cs_core_contribution(self):
        """
        Parse the core contribution of NMR chemical shielding.

        Returns:
            list[list]: G0 contribution matrix.
        """
        header_pattern = r"^\s+Core NMR properties\s*$\n\n^\s+typ\s+El\s+Core shift \(ppm\)\s*$\n^\s+-{20,}$\n"
        row_pattern = r"\d+\s+(?P<element>[A-Z][a-z]?\w?)\s+(?P<shift>[-]?\d+\.\d+)"
        footer_pattern = r"\s+-{20,}\s*$"
        self.read_table_pattern(
            header_pattern,
            row_pattern,
            footer_pattern,
            postprocess=str,
            last_one_only=True,
            attribute_name="cs_core_contribution",
        )
        core_contrib = {d["element"]: float(d["shift"]) for d in self.data["cs_core_contribution"]}
        self.data["cs_core_contribution"] = core_contrib

    def read_cs_raw_symmetrized_tensors(self):
        """
        Parse the matrix form of NMR tensor before corrected to table.

        Returns:
            nsymmetrized tensors list in the order of atoms.
        """
        header_pattern = r"\s+-{50,}\s+\s+Absolute Chemical Shift tensors\s+\s+-{50,}$"
        first_part_pattern = r"\s+UNSYMMETRIZED TENSORS\s+$"
        row_pattern = r"\s+".join([r"([-]?\d+\.\d+)"] * 3)
        unsym_footer_pattern = r"^\s+SYMMETRIZED TENSORS\s+$"

        with zopen(self.filename, "rt") as f:
            text = f.read()
        unsym_table_pattern_text = header_pattern + first_part_pattern + r"(?P<table_body>.+)" + unsym_footer_pattern
        table_pattern = re.compile(unsym_table_pattern_text, re.MULTILINE | re.DOTALL)
        rp = re.compile(row_pattern)
        m = table_pattern.search(text)
        if m:
            table_text = m.group("table_body")
            micro_header_pattern = r"ion\s+\d+"
            micro_table_pattern_text = micro_header_pattern + r"\s*^(?P<table_body>(?:\s*" + row_pattern + r")+)\s+"
            micro_table_pattern = re.compile(micro_table_pattern_text, re.MULTILINE | re.DOTALL)
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
        Parses the NMR Electric Field Gradient Raw Tensors.

        Returns:
            A list of Electric Field Gradient Tensors in the order of Atoms from OUTCAR
        """
        header_pattern = (
            r"Electric field gradients \(V/A\^2\)\n-*\n ion\s+V_xx\s+V_yy\s+V_zz\s+V_xy\s+V_xz\s+V_yz\n-*\n"
        )

        row_pattern = r"\d+\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)"
        footer_pattern = r"-*\n"

        data = self.read_table_pattern(header_pattern, row_pattern, footer_pattern, postprocess=float)
        tensors = [make_symmetric_matrix_from_upper_tri(d) for d in data]
        self.data["unsym_efg_tensor"] = tensors
        return tensors

    def read_nmr_efg(self):
        """
        Parse the NMR Electric Field Gradient interpreted values.

        Returns:
            Electric Field Gradient tensors as a list of dict in the order of atoms from OUTCAR.
            Each dict key/value pair corresponds to a component of the tensors.
        """
        header_pattern = (
            r"^\s+NMR quadrupolar parameters\s+$\n"
            r"^\s+Cq : quadrupolar parameter\s+Cq=e[*]Q[*]V_zz/h$\n"
            r"^\s+eta: asymmetry parameters\s+\(V_yy - V_xx\)/ V_zz$\n"
            r"^\s+Q  : nuclear electric quadrupole moment in mb \(millibarn\)$\n"
            r"^-{50,}$\n"
            r"^\s+ion\s+Cq\(MHz\)\s+eta\s+Q \(mb\)\s+$\n"
            r"^-{50,}\s*$\n"
        )
        row_pattern = (
            r"\d+\s+(?P<cq>[-]?\d+\.\d+)\s+(?P<eta>[-]?\d+\.\d+)\s+(?P<nuclear_quadrupole_moment>[-]?\d+\.\d+)"
        )
        footer_pattern = r"-{50,}\s*$"
        self.read_table_pattern(
            header_pattern,
            row_pattern,
            footer_pattern,
            postprocess=float,
            last_one_only=True,
            attribute_name="efg",
        )

    def read_elastic_tensor(self):
        """
        Parse the elastic tensor data.

        Returns:
            6x6 array corresponding to the elastic tensor from the OUTCAR.
        """
        header_pattern = r"TOTAL ELASTIC MODULI \(kBar\)\s+Direction\s+([X-Z][X-Z]\s+)+\-+"
        row_pattern = r"[X-Z][X-Z]\s+" + r"\s+".join([r"(\-*[\.\d]+)"] * 6)
        footer_pattern = r"\-+"
        et_table = self.read_table_pattern(header_pattern, row_pattern, footer_pattern, postprocess=float)
        self.data["elastic_tensor"] = et_table

    def read_piezo_tensor(self):
        """Parse the piezo tensor data."""
        header_pattern = r"PIEZOELECTRIC TENSOR  for field in x, y, z\s+\(C/m\^2\)\s+([X-Z][X-Z]\s+)+\-+"
        row_pattern = r"[x-z]\s+" + r"\s+".join([r"(\-*[\.\d]+)"] * 6)
        footer_pattern = r"BORN EFFECTIVE"
        pt_table = self.read_table_pattern(header_pattern, row_pattern, footer_pattern, postprocess=float)
        self.data["piezo_tensor"] = pt_table

    def read_onsite_density_matrices(self):
        """
        Parse the onsite density matrices, returns list with index corresponding
        to atom index in Structure.
        """
        # matrix size will vary depending on if d or f orbitals are present
        # therefore regex assumes f, but filter out None values if d

        header_pattern = r"spin component  1\n"
        row_pattern = r"[^\S\r\n]*(?:(-?[\d.]+))" + r"(?:[^\S\r\n]*(-?[\d.]+)[^\S\r\n]*)?" * 6 + r".*?"
        footer_pattern = r"\nspin component  2"
        spin1_component = self.read_table_pattern(
            header_pattern,
            row_pattern,
            footer_pattern,
            postprocess=lambda x: float(x) if x else None,
            last_one_only=False,
        )

        # filter out None values
        spin1_component = [[[e for e in row if e is not None] for row in matrix] for matrix in spin1_component]

        # and repeat for Spin.down

        header_pattern = r"spin component  2\n"
        row_pattern = r"[^\S\r\n]*(?:([\d.-]+))" + r"(?:[^\S\r\n]*(-?[\d.]+)[^\S\r\n]*)?" * 6 + r".*?"
        footer_pattern = r"\n occupancies and eigenvectors"
        spin2_component = self.read_table_pattern(
            header_pattern,
            row_pattern,
            footer_pattern,
            postprocess=lambda x: float(x) if x else None,
            last_one_only=False,
        )

        spin2_component = [[[e for e in row if e is not None] for row in matrix] for matrix in spin2_component]

        self.data["onsite_density_matrices"] = [
            {Spin.up: spin1_component[idx], Spin.down: spin2_component[idx]} for idx in range(len(spin1_component))
        ]

    def read_corrections(self, reverse=True, terminate_on_match=True):
        """
        Reads the dipol qudropol corrections into the
        Outcar.data["dipol_quadrupol_correction"].

        Args:
            reverse (bool): Whether to start from end of OUTCAR. Defaults to True.
            terminate_on_match (bool): Whether to terminate once match is found. Defaults to True.
        """
        patterns = {"dipol_quadrupol_correction": r"dipol\+quadrupol energy correction\s+([\d\-\.]+)"}
        self.read_pattern(
            patterns,
            reverse=reverse,
            terminate_on_match=terminate_on_match,
            postprocess=float,
        )
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
            r"([\d\-\.]+)",
        }
        self.read_pattern(
            patterns,
            reverse=reverse,
            terminate_on_match=terminate_on_match,
            postprocess=str,
        )
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
            p_ion = spin up + spin down summed.

        (See VASP section "LBERRY,  IGPAR,  NPPSTR,  DIPOL tags" for info on
        what these are).
        """
        # variables to be filled
        self.er_ev = {}  # will  be  dict (Spin.up/down) of array(3*float)
        self.er_bp = {}  # will  be  dics (Spin.up/down) of array(3*float)
        self.er_ev_tot = None  # will be array(3*float)
        self.er_bp_tot = None  # will be array(3*float)
        self.p_elec = self.p_ion = None
        try:
            search = []

            # Nonspin cases
            def er_ev(results, match):
                results.er_ev[Spin.up] = np.array(map(float, match.groups()[1:4])) / 2
                results.er_ev[Spin.down] = results.er_ev[Spin.up]
                results.context = 2

            search.append(
                [
                    r"^ *e<r>_ev=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)",
                    None,
                    er_ev,
                ]
            )

            def er_bp(results, match):
                results.er_bp[Spin.up] = np.array([float(match.group(i)) for i in range(1, 4)]) / 2
                results.er_bp[Spin.down] = results.er_bp[Spin.up]

            search.append(
                [
                    r"^ *e<r>_bp=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)",
                    lambda results, _line: results.context == 2,
                    er_bp,
                ]
            )

            # Spin cases
            def er_ev_up(results, match):
                results.er_ev[Spin.up] = np.array([float(match.group(i)) for i in range(1, 4)])
                results.context = Spin.up

            search.append(
                [
                    r"^.*Spin component 1 *e<r>_ev=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)",
                    None,
                    er_ev_up,
                ]
            )

            def er_bp_up(results, match):
                results.er_bp[Spin.up] = np.array(
                    [
                        float(match.group(1)),
                        float(match.group(2)),
                        float(match.group(3)),
                    ]
                )

            search.append(
                [
                    r"^ *e<r>_bp=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)",
                    lambda results, _line: results.context == Spin.up,
                    er_bp_up,
                ]
            )

            def er_ev_dn(results, match):
                results.er_ev[Spin.down] = np.array(
                    [
                        float(match.group(1)),
                        float(match.group(2)),
                        float(match.group(3)),
                    ]
                )
                results.context = Spin.down

            search.append(
                [
                    r"^.*Spin component 2 *e<r>_ev=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)",
                    None,
                    er_ev_dn,
                ]
            )

            def er_bp_dn(results, match):
                results.er_bp[Spin.down] = np.array([float(match.group(i)) for i in range(1, 4)])

            search.append(
                [
                    r"^ *e<r>_bp=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)",
                    lambda results, _line: results.context == Spin.down,
                    er_bp_dn,
                ]
            )

            # Always present spin/non-spin
            def p_elc(results, match):
                results.p_elc = np.array([float(match.group(i)) for i in range(1, 4)])

            search.append(
                [
                    r"^.*Total electronic dipole moment: "
                    r"*p\[elc\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                    r"*([-0-9.Ee+]*) *\)",
                    None,
                    p_elc,
                ]
            )

            def p_ion(results, match):
                results.p_ion = np.array([float(match.group(i)) for i in range(1, 4)])

            search.append(
                [
                    r"^.*ionic dipole moment: *p\[ion\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)",
                    None,
                    p_ion,
                ]
            )

            self.context = None
            self.er_ev = {Spin.up: None, Spin.down: None}
            self.er_bp = {Spin.up: None, Spin.down: None}

            micro_pyawk(self.filename, search, self)

            if self.er_ev[Spin.up] is not None and self.er_ev[Spin.down] is not None:
                self.er_ev_tot = self.er_ev[Spin.up] + self.er_ev[Spin.down]

            if self.er_bp[Spin.up] is not None and self.er_bp[Spin.down] is not None:
                self.er_bp_tot = self.er_bp[Spin.up] + self.er_bp[Spin.down]

        except Exception:
            self.er_ev_tot = self.er_bp_tot = None
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

        search.append(
            [
                r"INTERNAL STRAIN TENSOR FOR ION\s+(\d+)\s+for displacements in x,y,z  \(eV/Angst\):",
                None,
                internal_strain_start,
            ]
        )

        def internal_strain_data(results, match):
            if match.group(1).lower() == "x":
                index = 0
            elif match.group(1).lower() == "y":
                index = 1
            elif match.group(1).lower() == "z":
                index = 2
            else:
                raise Exception(f"Couldn't parse row index from symbol for internal strain tensor: {match.group(1)}")
            results.internal_strain_tensor[results.internal_strain_ion][index] = np.array(
                [float(match.group(i)) for i in range(2, 8)]
            )
            if index == 2:
                results.internal_strain_ion = None

        search.append(
            [
                r"^\s+([x,y,z])\s+" + r"([-]?\d+\.\d+)\s+" * 6,
                lambda results, _line: results.internal_strain_ion is not None,
                internal_strain_data,
            ]
        )

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

            search.append(
                [
                    r"MACROSCOPIC STATIC DIELECTRIC TENSOR \(",
                    None,
                    dielectric_section_start,
                ]
            )

            def dielectric_section_start2(results, match):
                results.dielectric_index = 0

            search.append(
                [
                    r"-------------------------------------",
                    lambda results, _line: results.dielectric_index == -1,
                    dielectric_section_start2,
                ]
            )

            def dielectric_data(results, match):
                results.dielectric_tensor[results.dielectric_index, :] = np.array(
                    [float(match.group(i)) for i in range(1, 4)]
                )
                results.dielectric_index += 1

            search.append(
                [
                    r"^ *([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+) *$",
                    lambda results, _line: results.dielectric_index >= 0
                    if results.dielectric_index is not None
                    else None,
                    dielectric_data,
                ]
            )

            def dielectric_section_stop(results, match):
                results.dielectric_index = None

            search.append(
                [
                    r"-------------------------------------",
                    lambda results, _line: results.dielectric_index >= 1
                    if results.dielectric_index is not None
                    else None,
                    dielectric_section_stop,
                ]
            )

            self.dielectric_index = None
            self.dielectric_tensor = np.zeros((3, 3))

            def piezo_section_start(results, _match):
                results.piezo_index = 0

            search.append(
                [
                    r"PIEZOELECTRIC TENSOR  for field in x, y, z        \(C/m\^2\)",
                    None,
                    piezo_section_start,
                ]
            )

            def piezo_data(results, match):
                results.piezo_tensor[results.piezo_index, :] = np.array([float(match.group(i)) for i in range(1, 7)])
                results.piezo_index += 1

            search.append(
                [
                    r"^ *[xyz] +([-0-9.Ee+]+) +([-0-9.Ee+]+)"
                    r" +([-0-9.Ee+]+) *([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+)*$",
                    lambda results, _line: results.piezo_index >= 0 if results.piezo_index is not None else None,
                    piezo_data,
                ]
            )

            def piezo_section_stop(results, _match):
                results.piezo_index = None

            search.append(
                [
                    r"-------------------------------------",
                    lambda results, _line: results.piezo_index >= 1 if results.piezo_index is not None else None,
                    piezo_section_stop,
                ]
            )

            self.piezo_index = None
            self.piezo_tensor = np.zeros((3, 6))

            def born_section_start(results, _match):
                results.born_ion = -1

            search.append([r"BORN EFFECTIVE CHARGES ", None, born_section_start])

            def born_ion(results, match):
                results.born_ion = int(match.group(1)) - 1
                results.born.append(np.zeros((3, 3)))

            search.append(
                [
                    r"ion +([0-9]+)",
                    lambda results, _line: results.born_ion is not None,
                    born_ion,
                ]
            )

            def born_data(results, match):
                results.born[results.born_ion][int(match.group(1)) - 1, :] = np.array(
                    [float(match.group(i)) for i in range(2, 5)]
                )

            search.append(
                [
                    r"^ *([1-3]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+)$",
                    lambda results, _line: results.born_ion >= 0 if results.born_ion is not None else results.born_ion,
                    born_data,
                ]
            )

            def born_section_stop(results, _match):
                results.born_ion = None

            search.append(
                [
                    r"-------------------------------------",
                    lambda results, _line: results.born_ion >= 1 if results.born_ion is not None else results.born_ion,
                    born_section_stop,
                ]
            )

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

            def dielectric_section_start(results, _match):
                results.dielectric_ionic_index = -1

            search.append(
                [
                    r"MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC",
                    None,
                    dielectric_section_start,
                ]
            )

            def dielectric_section_start2(results, _match):
                results.dielectric_ionic_index = 0

            search.append(
                [
                    r"-------------------------------------",
                    lambda results, _line: results.dielectric_ionic_index == -1
                    if results.dielectric_ionic_index is not None
                    else results.dielectric_ionic_index,
                    dielectric_section_start2,
                ]
            )

            def dielectric_data(results, match):
                results.dielectric_ionic_tensor[results.dielectric_ionic_index, :] = np.array(
                    [float(match.group(i)) for i in range(1, 4)]
                )
                results.dielectric_ionic_index += 1

            search.append(
                [
                    r"^ *([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+) *$",
                    lambda results, _line: results.dielectric_ionic_index >= 0
                    if results.dielectric_ionic_index is not None
                    else results.dielectric_ionic_index,
                    dielectric_data,
                ]
            )

            def dielectric_section_stop(results, _match):
                results.dielectric_ionic_index = None

            search.append(
                [
                    r"-------------------------------------",
                    lambda results, _line: results.dielectric_ionic_index >= 1
                    if results.dielectric_ionic_index is not None
                    else results.dielectric_ionic_index,
                    dielectric_section_stop,
                ]
            )

            self.dielectric_ionic_index = None
            self.dielectric_ionic_tensor = np.zeros((3, 3))

            def piezo_section_start(results, _match):
                results.piezo_ionic_index = 0

            search.append(
                [
                    r"PIEZOELECTRIC TENSOR IONIC CONTR  for field in x, y, z        ",
                    None,
                    piezo_section_start,
                ]
            )

            def piezo_data(results, match):
                results.piezo_ionic_tensor[results.piezo_ionic_index, :] = np.array(
                    [float(match.group(i)) for i in range(1, 7)]
                )
                results.piezo_ionic_index += 1

            search.append(
                [
                    r"^ *[xyz] +([-0-9.Ee+]+) +([-0-9.Ee+]+)"
                    r" +([-0-9.Ee+]+) *([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+)*$",
                    lambda results, _line: results.piezo_ionic_index >= 0
                    if results.piezo_ionic_index is not None
                    else results.piezo_ionic_index,
                    piezo_data,
                ]
            )

            def piezo_section_stop(results, _match):
                results.piezo_ionic_index = None

            search.append(
                [
                    "-------------------------------------",
                    lambda results, _line: results.piezo_ionic_index >= 1
                    if results.piezo_ionic_index is not None
                    else results.piezo_ionic_index,
                    piezo_section_stop,
                ]
            )

            self.piezo_ionic_index = None
            self.piezo_ionic_tensor = np.zeros((3, 6))

            micro_pyawk(self.filename, search, self)

            self.dielectric_ionic_tensor = self.dielectric_ionic_tensor.tolist()
            self.piezo_ionic_tensor = self.piezo_ionic_tensor.tolist()

        except Exception:
            raise Exception("ionic part of LEPSILON OUTCAR could not be parsed.")

    def read_lcalcpol(self):
        """
        Reads the lcalpol.

        # TODO: Document the actual variables.
        """
        self.p_elec = self.p_sp1 = self.p_sp2 = self.p_ion = None
        try:
            search = []

            # Always present spin/non-spin
            def p_elec(results, match):
                results.p_elec = np.array(
                    [
                        float(match.group(1)),
                        float(match.group(2)),
                        float(match.group(3)),
                    ]
                )

            search.append(
                [
                    r"^.*Total electronic dipole moment: "
                    r"*p\[elc\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) "
                    r"*([-0-9.Ee+]*) *\)",
                    None,
                    p_elec,
                ]
            )

            # If spin-polarized (and not noncollinear)
            # save spin-polarized electronic values
            if self.spin and not self.noncollinear:

                def p_sp1(results, match):
                    results.p_sp1 = np.array(
                        [
                            float(match.group(1)),
                            float(match.group(2)),
                            float(match.group(3)),
                        ]
                    )

                search.append(
                    [
                        r"^.*p\[sp1\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)",
                        None,
                        p_sp1,
                    ]
                )

                def p_sp2(results, match):
                    results.p_sp2 = np.array(
                        [
                            float(match.group(1)),
                            float(match.group(2)),
                            float(match.group(3)),
                        ]
                    )

                search.append(
                    [
                        r"^.*p\[sp2\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)",
                        None,
                        p_sp2,
                    ]
                )

            def p_ion(results, match):
                results.p_ion = np.array(
                    [
                        float(match.group(1)),
                        float(match.group(2)),
                        float(match.group(3)),
                    ]
                )

            search.append(
                [
                    r"^.*Ionic dipole moment: *p\[ion\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)",
                    None,
                    p_ion,
                ]
            )

            micro_pyawk(self.filename, search, self)

            # fix polarization units in new versions of vasp
            regex = r"^.*Ionic dipole moment: .*"
            search = [[regex, None, lambda x, y: x.append(y.group(0))]]
            r = micro_pyawk(self.filename, search, [])

            if "|e|" in r[0]:
                self.p_elec *= -1
                self.p_ion *= -1
                if self.spin and not self.noncollinear:
                    self.p_sp1 *= -1
                    self.p_sp2 *= -1

        except Exception as exc:
            print(exc.args)
            raise Exception("LCALCPOL OUTCAR could not be parsed.") from exc

    def read_pseudo_zval(self):
        """Create pseudopotential ZVAL dictionary."""
        # pylint: disable=E1101
        try:

            def atom_symbols(results, match):
                element_symbol = match.group(1)
                if not hasattr(results, "atom_symbols"):
                    results.atom_symbols = []
                results.atom_symbols.append(element_symbol.strip())

            def zvals(results, match):
                zvals = match.group(1)
                results.zvals = map(float, re.findall(r"-?\d+\.\d*", zvals))

            search = []
            search.append([r"(?<=VRHFIN =)(.*)(?=:)", None, atom_symbols])
            search.append([r"^\s+ZVAL.*=(.*)", None, zvals])

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

                        # the average core potentials of up to 5 elements are
                        # given per line; the potentials are separated by several
                        # spaces and numbered from 1 to natoms; the potentials are
                        # parsed in a fixed width format
                        npots = int((len(line) - 1) / 17)
                        for i in range(npots):
                            start = i * 17
                            ap.append(float(line[start + 8 : start + 17]))

        return aps

    def as_dict(self):
        """MSONable dict."""
        dct = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "efermi": self.efermi,
            "run_stats": self.run_stats,
            "magnetization": self.magnetization,
            "charge": self.charge,
            "total_magnetization": self.total_mag,
            "nelect": self.nelect,
            "is_stopped": self.is_stopped,
            "drift": self.drift,
            "ngf": self.ngf,
            "sampling_radii": self.sampling_radii,
            "electrostatic_potential": self.electrostatic_potential,
        }

        if self.lepsilon:
            dct.update(
                {
                    "piezo_tensor": self.piezo_tensor,
                    "dielectric_tensor": self.dielectric_tensor,
                    "born": self.born,
                }
            )

        if self.dfpt:
            dct["internal_strain_tensor"] = self.internal_strain_tensor

        if self.dfpt and self.lepsilon:
            dct.update(
                {
                    "piezo_ionic_tensor": self.piezo_ionic_tensor,
                    "dielectric_ionic_tensor": self.dielectric_ionic_tensor,
                }
            )

        if self.lcalcpol:
            dct.update({"p_elec": self.p_elec, "p_ion": self.p_ion})
            if self.spin and not self.noncollinear:
                dct.update({"p_sp1": self.p_sp1, "p_sp2": self.p_sp2})
            dct["zval_dict"] = self.zval_dict

        if self.nmr_cs:
            dct.update(
                {
                    "nmr_cs": {
                        "valence and core": self.data["chemical_shielding"]["valence_and_core"],
                        "valence_only": self.data["chemical_shielding"]["valence_only"],
                        "g0": self.data["cs_g0_contribution"],
                        "core": self.data["cs_core_contribution"],
                        "raw": self.data["unsym_cs_tensor"],
                    }
                }
            )

        if self.nmr_efg:
            dct.update(
                {
                    "nmr_efg": {
                        "raw": self.data["unsym_efg_tensor"],
                        "parameters": self.data["efg"],
                    }
                }
            )

        if self.has_onsite_density_matrices:
            # cast Spin to str for consistency with electronic_structure
            # TODO: improve handling of Enum (de)serialization in monty
            onsite_density_matrices = [{str(k): v for k, v in d.items()} for d in self.data["onsite_density_matrices"]]
            dct["onsite_density_matrices"] = onsite_density_matrices

        return dct

    def read_fermi_contact_shift(self):
        """
        Output example:
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
         [0.056, 0.056, 0.321, -0.048, 0.321]] from 'fch' data.
        """
        # Fermi contact (isotropic) hyperfine coupling parameter (MHz)
        header_pattern1 = (
            r"\s*Fermi contact \(isotropic\) hyperfine coupling parameter \(MHz\)\s+"
            r"\s*\-+"
            r"\s*ion\s+A_pw\s+A_1PS\s+A_1AE\s+A_1c\s+A_tot\s+"
            r"\s*\-+"
        )
        row_pattern1 = r"(?:\d+)\s+" + r"\s+".join([r"([-]?\d+\.\d+)"] * 5)
        footer_pattern = r"\-+"
        fch_table = self.read_table_pattern(
            header_pattern1,
            row_pattern1,
            footer_pattern,
            postprocess=float,
            last_one_only=True,
        )

        # Dipolar hyperfine coupling parameters (MHz)
        header_pattern2 = (
            r"\s*Dipolar hyperfine coupling parameters \(MHz\)\s+"
            r"\s*\-+"
            r"\s*ion\s+A_xx\s+A_yy\s+A_zz\s+A_xy\s+A_xz\s+A_yz\s+"
            r"\s*\-+"
        )
        row_pattern2 = r"(?:\d+)\s+" + r"\s+".join([r"([-]?\d+\.\d+)"] * 6)
        dh_table = self.read_table_pattern(
            header_pattern2,
            row_pattern2,
            footer_pattern,
            postprocess=float,
            last_one_only=True,
        )

        # Total hyperfine coupling parameters after diagonalization (MHz)
        header_pattern3 = (
            r"\s*Total hyperfine coupling parameters after diagonalization \(MHz\)\s+"
            r"\s*\(convention: \|A_zz\| > \|A_xx\| > \|A_yy\|\)\s+"
            r"\s*\-+"
            r"\s*ion\s+A_xx\s+A_yy\s+A_zz\s+asymmetry \(A_yy - A_xx\)/ A_zz\s+"
            r"\s*\-+"
        )
        row_pattern3 = r"(?:\d+)\s+" + r"\s+".join([r"([-]?\d+\.\d+)"] * 4)
        th_table = self.read_table_pattern(
            header_pattern3,
            row_pattern3,
            footer_pattern,
            postprocess=float,
            last_one_only=True,
        )

        fc_shift_table = {"fch": fch_table, "dh": dh_table, "th": th_table}

        self.data["fermi_contact_shift"] = fc_shift_table


class VolumetricData(BaseVolumetricData):
    """
    Container for volumetric data that allows
    for reading/writing with Poscar-type data.
    """

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
        dim = dimline = None
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
                            dataset[data_count % dim[0], no_x % dim[1], no_x // dim[1]] = float(tok)
                            data_count += 1
                    if data_count >= ngrid_pts:
                        read_dataset = False
                        data_count = 0
                        all_dataset.append(dataset)
                elif not poscar_read:
                    if line != "" or len(poscar_string) == 0:
                        poscar_string.append(line)
                    elif line == "":
                        poscar = Poscar.from_str("\n".join(poscar_string))
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
                data = {
                    "total": all_dataset[0],
                    "diff_x": all_dataset[1],
                    "diff_y": all_dataset[2],
                    "diff_z": all_dataset[3],
                }
                data_aug = {
                    "total": all_dataset_aug.get(0),
                    "diff_x": all_dataset_aug.get(1),
                    "diff_y": all_dataset_aug.get(2),
                    "diff_z": all_dataset_aug.get(3),
                }

                # construct a "diff" dict for scalar-like magnetization density,
                # referenced to an arbitrary direction (using same method as
                # pymatgen.electronic_structure.core.Magmom, see
                # Magmom documentation for justification for this)
                # TODO: re-examine this, and also similar behavior in
                # Magmom - @mkhorton
                # TODO: does CHGCAR change with different SAXIS?
                diff_xyz = np.array([data["diff_x"], data["diff_y"], data["diff_z"]])
                diff_xyz = diff_xyz.reshape((3, dim[0] * dim[1] * dim[2]))
                ref_direction = np.array([1.01, 1.02, 1.03])
                ref_sign = np.sign(np.dot(ref_direction, diff_xyz))
                diff = np.multiply(np.linalg.norm(diff_xyz, axis=0), ref_sign)
                data["diff"] = diff.reshape((dim[0], dim[1], dim[2]))

            elif len(all_dataset) == 2:
                data = {"total": all_dataset[0], "diff": all_dataset[1]}
                data_aug = {
                    "total": all_dataset_aug.get(0),
                    "diff": all_dataset_aug.get(1),
                }
            else:
                data = {"total": all_dataset[0]}
                data_aug = {"total": all_dataset_aug.get(0)}
            return poscar, data, data_aug

    def write_file(self, file_name, vasp4_compatible=False):
        """
        Write the VolumetricData object to a vasp compatible file.

        Args:
            file_name (str): Path to a file
            vasp4_compatible (bool): True if the format is vasp4 compatible
        """

        def _print_fortran_float(flt):
            """
            Fortran codes print floats with a leading zero in scientific
            notation. When writing CHGCAR files, we adopt this convention
            to ensure written CHGCAR files are byte-to-byte identical to
            their input files as far as possible.

            Args:
                flt (float): Float to print.

            Returns:
                str: String representation of float in Fortran format.
            """
            s = f"{flt:.10E}"
            if flt >= 0:
                return f"0.{s[0]}{s[2:12]}E{int(s[13:]) + 1:+03}"
            return f"-.{s[1]}{s[3:13]}E{int(s[14:]) + 1:+03}"

        with zopen(file_name, "wt") as f:
            p = Poscar(self.structure)

            # use original name if it's been set (e.g. from Chgcar)
            comment = getattr(self, "name", p.comment)

            lines = comment + "\n"
            lines += "   1.00000000000000\n"
            for vec in self.structure.lattice.matrix:
                lines += f" {vec[0]:12.6f}{vec[1]:12.6f}{vec[2]:12.6f}\n"
            if not vasp4_compatible:
                lines += "".join(f"{s:5}" for s in p.site_symbols) + "\n"
            lines += "".join(f"{x:6}" for x in p.natoms) + "\n"
            lines += "Direct\n"
            for site in self.structure:
                a, b, c = site.frac_coords
                lines += f"{a:10.6f}{b:10.6f}{c:10.6f}\n"
            lines += " \n"
            f.write(lines)
            a = self.dim

            def write_spin(data_type):
                lines = []
                count = 0
                f.write(f"   {a[0]}   {a[1]}   {a[2]}\n")
                for k, j, i in itertools.product(list(range(a[2])), list(range(a[1])), list(range(a[0]))):
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


class Locpot(VolumetricData):
    """Simple object for reading a LOCPOT file."""

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
        """Read a LOCPOT file.

        Args:
            filename (str): Path to LOCPOT file.

        Returns:
            Locpot
        """
        (poscar, data, data_aug) = VolumetricData.parse_file(filename)
        return cls(poscar, data, **kwargs)


class Chgcar(VolumetricData):
    """Simple object for reading a CHGCAR file."""

    def __init__(self, poscar, data, data_aug=None):
        """
        Args:
            poscar (Poscar or Structure): Object containing structure.
            data: Actual data.
            data_aug: Augmentation charge data.
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
    def from_file(filename: str):
        """Read a CHGCAR file.

        Args:
            filename (str): Path to CHGCAR file.

        Returns:
            Chgcar
        """
        poscar, data, data_aug = VolumetricData.parse_file(filename)
        return Chgcar(poscar, data, data_aug=data_aug)

    @property
    def net_magnetization(self):
        """Net magnetization from Chgcar"""
        if self.is_spin_polarized:
            return np.sum(self.data["diff"])
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

        Returns:
            Elfcar
        """
        (poscar, data, data_aug) = VolumetricData.parse_file(filename)
        return cls(poscar, data)

    def get_alpha(self):
        """Get the parameter alpha where ELF = 1/(1+alpha^2)."""
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

    Attributes:
        data (dict): The PROCAR data of the form below. It should VASP uses 1-based indexing,
            but all indices are converted to 0-based here.
            { spin: nd.array accessed with (k-point index, band index, ion index, orbital index) }
        weights (numpy.ndarray): The weights associated with each k-point as an nd.array of length nkpoints.
        phase_factors (dict): Phase factors, where present (e.g. LORBIT = 12). A dict of the form:
            { spin: complex nd.array accessed with (k-point index, band index, ion index, orbital index) }
        nbands (int): Number of bands.
        nkpoints (int): Number of k-points.
        nions (int): Number of ions.
    """

    def __init__(self, filename):
        """
        Args:
            filename: Name of file containing PROCAR.
        """
        headers = None

        with zopen(filename, "rt") as file_handle:
            preambleexpr = re.compile(r"# of k-points:\s*(\d+)\s+# of bands:\s*(\d+)\s+# of ions:\s*(\d+)")
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
            for line in file_handle:
                line = line.strip()
                if bandexpr.match(line):
                    m = bandexpr.match(line)
                    current_band = int(m.group(1)) - 1
                    done = False
                elif kpointexpr.match(line):
                    m = kpointexpr.match(line)
                    current_kpoint = int(m.group(1)) - 1
                    weights[current_kpoint] = float(m.group(2))
                    if current_kpoint == 0:
                        spin = Spin.up if spin == Spin.down else Spin.down
                    done = False
                elif headers is None and ionexpr.match(line):
                    headers = line.split()
                    headers.pop(0)
                    headers.pop(-1)

                    data = defaultdict(lambda: np.zeros((n_kpoints, n_bands, n_ions, len(headers))))

                    phase_factors = defaultdict(
                        lambda: np.full(
                            (n_kpoints, n_bands, n_ions, len(headers)),
                            np.NaN,
                            dtype=np.complex128,
                        )
                    )
                elif expr.match(line):
                    tokens = line.split()
                    index = int(tokens.pop(0)) - 1
                    num_data = np.array([float(t) for t in tokens[: len(headers)]])
                    if not done:
                        data[spin][current_kpoint, current_band, index, :] = num_data
                    elif len(tokens) > len(headers):
                        # new format of PROCAR (vasp 5.4.4)
                        num_data = np.array([float(t) for t in tokens[: 2 * len(headers)]])
                        for orb in range(len(headers)):
                            phase_factors[spin][current_kpoint, current_band, index, orb] = complex(
                                num_data[2 * orb], num_data[2 * orb + 1]
                            )
                    elif np.isnan(phase_factors[spin][current_kpoint, current_band, index, 0]):
                        # old format of PROCAR (vasp 5.4.1 and before)
                        phase_factors[spin][current_kpoint, current_band, index, :] = num_data
                    else:
                        phase_factors[spin][current_kpoint, current_band, index, :] += 1j * num_data
                elif line.startswith("tot"):
                    done = True
                elif preambleexpr.match(line):
                    m = preambleexpr.match(line)
                    n_kpoints = int(m.group(1))
                    n_bands = int(m.group(2))
                    n_ions = int(m.group(3))
                    weights = np.zeros(n_kpoints)

            self.nkpoints = n_kpoints
            self.nbands = n_bands
            self.nions = n_ions
            self.weights = weights
            self.orbitals = headers
            self.data = data
            self.phase_factors = phase_factors

    def get_projection_on_elements(self, structure: Structure):
        """
        Method returning a dictionary of projections on elements.

        Args:
            structure (Structure): Input structure.

        Returns:
            a dictionary in the {Spin.up:[k index][b index][{Element:values}]]
        """
        dico: dict[Spin, list] = {}
        for spin in self.data:
            dico[spin] = [[defaultdict(float) for i in range(self.nkpoints)] for j in range(self.nbands)]

        for iat in range(self.nions):
            name = structure.species[iat].symbol
            for spin, d in self.data.items():
                for k, b in itertools.product(range(self.nkpoints), range(self.nbands)):
                    dico[spin][b][k][name] += np.sum(d[k, b, iat, :])

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
        return {
            spin: np.sum(d[:, :, atom_index, orbital_index] * self.weights[:, None]) for spin, d in self.data.items()
        }


class Oszicar:
    """
    A basic parser for an OSZICAR output from VASP. In general, while the
    OSZICAR is useful for a quick look at the output from a VASP run, we
    recommend that you use the Vasprun parser instead, which gives far richer
    information about a run.

    Attributes:
        electronic_steps (list): All electronic steps as a list of list of dict. e.g.,
            [[{"rms": 160.0, "E": 4507.24605593, "dE": 4507.2, "N": 1, "deps": -17777.0, "ncg": 16576}, ...], [....]
            where electronic_steps[index] refers the list of electronic steps in one ionic_step,
            electronic_steps[index][subindex] refers to a particular electronic step at subindex in ionic step at
            index. The dict of properties depends on the type of VASP run, but in general, "E", "dE" and "rms" should
            be present in almost all runs.
        ionic_steps (list): All ionic_steps as a list of dict, e.g.,
            [{"dE": -526.36, "E0": -526.36024, "mag": 0.0, "F": -526.36024}, ...]
            This is the typical output from VASP at the end of each ionic step. The stored dict might be different
            depending on the type of VASP run.
    """

    def __init__(self, filename):
        """
        Args:
            filename (str): Filename of file to parse.
        """
        electronic_steps = []
        ionic_steps = []
        ionic_general_pattern = re.compile(r"(\w+)=\s*(\S+)")
        electronic_pattern = re.compile(r"\s*\w+\s*:(.*)")

        def smart_convert(header, num):
            try:
                if header in ("N", "ncg"):
                    return int(num)
                return float(num)
            except ValueError:
                return "--"

        header = []
        with zopen(filename, "rt") as fid:
            for line in fid:
                m = electronic_pattern.match(line.strip())
                if m:
                    tokens = m.group(1).split()
                    data = {header[i]: smart_convert(header[i], tokens[i]) for i in range(len(tokens))}
                    if tokens[0] == "1":
                        electronic_steps.append([data])
                    else:
                        electronic_steps[-1].append(data)
                elif re.match(r"^\s*N\s+E\s*", line.strip()):
                    header = line.strip().replace("d eps", "deps").split()
                elif line.strip() != "":
                    # remove space first and apply field agnostic extraction
                    matches = re.findall(ionic_general_pattern, re.sub(r"d E ", "dE", line))
                    ionic_steps.append({key: float(value) for key, value in matches})

        self.electronic_steps = electronic_steps
        self.ionic_steps = ionic_steps

    @property
    def all_energies(self):
        """
        Compilation of all energies from all electronic steps and ionic steps
        as a tuple of list of energies, e.g.,
        ((4507.24605593, 143.824705755, -512.073149912, ...), ...).
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
        """Final energy from run."""
        return self.ionic_steps[-1]["E0"]

    def as_dict(self):
        """MSONable dict"""
        return {
            "electronic_steps": self.electronic_steps,
            "ionic_steps": self.ionic_steps,
        }


class VaspParseError(ParseError):
    """Exception class for VASP parsing."""


def get_band_structure_from_vasp_multiple_branches(dir_name, efermi=None, projections=False):
    """
    This method is used to get band structure info from a VASP directory. It
    takes into account that the run can be divided in several branches named
    "branch_x". If the run has not been divided in branches the method will
    turn to parsing vasprun.xml directly.

    The method returns None is there's a parsing error

    Args:
        dir_name: Directory containing all bandstructure runs.
        efermi: Efermi for bandstructure.
        projections: True if you want to get the data on site projections if
            any. Note that this is sometimes very large

    Returns:
        A BandStructure Object
    """
    # TODO: Add better error handling!!!
    if os.path.exists(f"{dir_name}/branch_0"):
        # get all branch dir names
        branch_dir_names = [os.path.abspath(d) for d in glob(f"{dir_name}/branch_*") if os.path.isdir(d)]

        # sort by the directory name (e.g, branch_10)
        sorted_branch_dir_names = sorted(branch_dir_names, key=lambda x: int(x.split("_")[-1]))

        # populate branches with Bandstructure instances
        branches = []
        for dname in sorted_branch_dir_names:
            xml_file = f"{dname}/vasprun.xml"
            if os.path.exists(xml_file):
                run = Vasprun(xml_file, parse_projected_eigen=projections)
                branches.append(run.get_band_structure(efermi=efermi))
            else:
                # It might be better to throw an exception
                warnings.warn(f"Skipping {dname}. Unable to find {xml_file}")

        return get_reconstructed_band_structure(branches, efermi)

    xml_file = f"{dir_name}/vasprun.xml"
    # Better handling of Errors
    if os.path.exists(xml_file):
        return Vasprun(xml_file, parse_projected_eigen=projections).get_band_structure(
            kpoints_filename=None, efermi=efermi
        )

    return None


class Xdatcar:
    """
    Class representing an XDATCAR file. Only tested with VASP 5.x files.

    Attributes:
        structures (list): List of structures parsed from XDATCAR.
        comment (str): Optional comment string.

    Authors: Ram Balachandran
    """

    def __init__(self, filename, ionicstep_start=1, ionicstep_end=None, comment=None):
        """
        Init a Xdatcar.

        Args:
            filename (str): Filename of input XDATCAR file.
            ionicstep_start (int): Starting number of ionic step.
            ionicstep_end (int): Ending number of ionic step.
            comment (str): Optional comment attached to this set of structures.
        """
        preamble = None
        coords_str = []
        structures = []
        preamble_done = False
        if ionicstep_start < 1:
            raise Exception("Start ionic step cannot be less than 1")
        if ionicstep_end is not None and ionicstep_start < 1:
            raise Exception("End ionic step cannot be less than 1")

        # pylint: disable=E1136
        ionicstep_cnt = 1
        with zopen(filename, "rt") as file:
            for line in file:
                line = line.strip()
                if preamble is None:
                    preamble = [line]
                    title = line
                elif title == line:
                    preamble_done = False
                    p = Poscar.from_str("\n".join([*preamble, "Direct", *coords_str]))
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
                    preamble = [line]
                elif not preamble_done:
                    if line == "" or "Direct configuration=" in line:
                        preamble_done = True
                        tmp_preamble = [preamble[0]]
                        for i in range(1, len(preamble)):
                            if preamble[0] != preamble[i]:
                                tmp_preamble.append(preamble[i])
                            else:
                                break
                        preamble = tmp_preamble
                    else:
                        preamble.append(line)
                elif line == "" or "Direct configuration=" in line:
                    p = Poscar.from_str("\n".join([*preamble, "Direct", *coords_str]))
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
                    coords_str.append(line)
            p = Poscar.from_str("\n".join([*preamble, "Direct", *coords_str]))
            if ionicstep_end is None:
                if ionicstep_cnt >= ionicstep_start:
                    structures.append(p.structure)
            elif ionicstep_start <= ionicstep_cnt < ionicstep_end:
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

    def concatenate(self, filename, ionicstep_start=1, ionicstep_end=None):
        """
        Concatenate structures in file to Xdatcar.

        Args:
            filename (str): Filename of XDATCAR file to be concatenated.
            ionicstep_start (int): Starting number of ionic step.
            ionicstep_end (int): Ending number of ionic step.
        TODO (rambalachandran): Requires a check to ensure if the new concatenating file
            has the same lattice structure and atoms as the Xdatcar class.
        """
        preamble = None
        coords_str = []
        structures = self.structures
        preamble_done = False
        if ionicstep_start < 1:
            raise Exception("Start ionic step cannot be less than 1")
        if ionicstep_end is not None and ionicstep_start < 1:
            raise Exception("End ionic step cannot be less than 1")

        # pylint: disable=E1136
        ionicstep_cnt = 1
        with zopen(filename, "rt") as f:
            for line in f:
                line = line.strip()
                if preamble is None:
                    preamble = [line]
                elif not preamble_done:
                    if line == "" or "Direct configuration=" in line:
                        preamble_done = True
                        tmp_preamble = [preamble[0]]
                        for i in range(1, len(preamble)):
                            if preamble[0] != preamble[i]:
                                tmp_preamble.append(preamble[i])
                            else:
                                break
                        preamble = tmp_preamble
                    else:
                        preamble.append(line)
                elif line == "" or "Direct configuration=" in line:
                    p = Poscar.from_str("\n".join([*preamble, "Direct", *coords_str]))
                    if ionicstep_end is None:
                        if ionicstep_cnt >= ionicstep_start:
                            structures.append(p.structure)
                    elif ionicstep_start <= ionicstep_cnt < ionicstep_end:
                        structures.append(p.structure)
                    ionicstep_cnt += 1
                    coords_str = []
                else:
                    coords_str.append(line)
            p = Poscar.from_str("\n".join([*preamble, "Direct", *coords_str]))
            if ionicstep_end is None:
                if ionicstep_cnt >= ionicstep_start:
                    structures.append(p.structure)
            elif ionicstep_start <= ionicstep_cnt < ionicstep_end:
                structures.append(p.structure)
        self.structures = structures

    @np.deprecate(message="Use get_str instead")
    def get_string(self, *args, **kwargs) -> str:
        return self.get_str(*args, **kwargs)

    def get_str(self, ionicstep_start: int = 1, ionicstep_end: int | None = None, significant_figures: int = 8) -> str:
        """
        Write  Xdatcar class to a string.

        Args:
            ionicstep_start (int): Starting number of ionic step.
            ionicstep_end (int): Ending number of ionic step.
            significant_figures (int): Number of significant figures.
        """
        if ionicstep_start < 1:
            raise Exception("Start ionic step cannot be less than 1")
        if ionicstep_end is not None and ionicstep_end < 1:
            raise Exception("End ionic step cannot be less than 1")
        latt = self.structures[0].lattice
        if np.linalg.det(latt.matrix) < 0:
            latt = Lattice(-latt.matrix)
        lines = [self.comment, "1.0", str(latt)]
        lines.append(" ".join(self.site_symbols))
        lines.append(" ".join(str(x) for x in self.natoms))
        format_str = f"{{:.{significant_figures}f}}"
        ionicstep_cnt = 1
        output_cnt = 1
        for cnt, structure in enumerate(self.structures):
            ionicstep_cnt = cnt + 1
            if ionicstep_end is None:
                if ionicstep_cnt >= ionicstep_start:
                    lines.append(f"Direct configuration={' '*(7-len(str(output_cnt)))}{output_cnt}")
                    for site in structure:
                        coords = site.frac_coords
                        line = " ".join(format_str.format(c) for c in coords)
                        lines.append(line)
                    output_cnt += 1
            elif ionicstep_start <= ionicstep_cnt < ionicstep_end:
                lines.append(f"Direct configuration={' '*(7-len(str(output_cnt)))}{output_cnt}")
                for site in structure:
                    coords = site.frac_coords
                    line = " ".join(format_str.format(c) for c in coords)
                    lines.append(line)
                output_cnt += 1
        return "\n".join(lines) + "\n"

    def write_file(self, filename, **kwargs):
        """
        Write Xdatcar class into a file.

        Args:
            filename (str): Filename of output XDATCAR file.
            **kwargs: Supported kwargs are the same as those for the
                Xdatcar.get_string method and are passed through directly.
        """
        with zopen(filename, "wt") as f:
            f.write(self.get_str(**kwargs))

    def __str__(self):
        return self.get_str()


class Dynmat:
    """
    Object for reading a DYNMAT file.

    Attributes:
        data (dict): A nested dict containing the DYNMAT data of the form:
            [atom <int>][disp <int>]['dispvec'] =
                displacement vector (part of first line in dynmat block, e.g. "0.01 0 0")
            [atom <int>][disp <int>]['dynmat'] =
                    <list> list of dynmat lines for this atom and this displacement

    Authors: Patrick Huck
    """

    def __init__(self, filename):
        """
        Args:
            filename: Name of file containing DYNMAT.
        """
        with zopen(filename, "rt") as f:
            lines = list(clean_lines(f.readlines()))
            self._nspecs, self._natoms, self._ndisps = map(int, lines[0].split())
            self._masses = map(float, lines[1].split())
            self.data = defaultdict(dict)
            atom, disp = None, None
            for idx, line in enumerate(lines[2:]):
                v = list(map(float, line.split()))
                if not idx % (self._natoms + 1):
                    atom, disp = map(int, v[:2])
                    if atom not in self.data:
                        self.data[atom] = {}
                    if disp not in self.data[atom]:
                        self.data[atom][disp] = {}
                    self.data[atom][disp]["dispvec"] = v[2:]
                else:
                    if "dynmat" not in self.data[atom][disp]:
                        self.data[atom][disp]["dynmat"] = []
                    self.data[atom][disp]["dynmat"].append(v)

    def get_phonon_frequencies(self):
        """Calculate phonon frequencies."""
        # TODO: the following is most likely not correct or suboptimal
        # hence for demonstration purposes only
        frequencies = []
        for k, v0 in self.data.items():
            for v1 in v0.itervalues():
                vec = map(abs, v1["dynmat"][k - 1])
                frequency = math.sqrt(sum(vec)) * 2.0 * math.pi * 15.633302  # THz
                frequencies.append(frequency)
        return frequencies

    @property
    def nspecs(self):
        """Returns the number of species."""
        return self._nspecs

    @property
    def natoms(self):
        """Returns the number of atoms."""
        return self._natoms

    @property
    def ndisps(self):
        """Returns the number of displacements."""
        return self._ndisps

    @property
    def masses(self):
        """Returns the list of atomic masses."""
        return list(self._masses)


def get_adjusted_fermi_level(efermi, cbm, band_structure):
    """
    When running a band structure computations the Fermi level needs to be
    take from the static run that gave the charge density used for the non-self
    consistent band structure run. Sometimes this Fermi level is however a
    little too low because of the mismatch between the uniform grid used in
    the static run and the band structure k-points (e.g., the VBM is on Gamma
    and the Gamma point is not in the uniform mesh). Here we use a procedure
    consisting in looking for energy levels higher than the static Fermi level
    (but lower than the LUMO) if any of these levels make the band structure
    appears insulating and not metallic anymore, we keep this adjusted fermi
    level. This procedure has shown to detect correctly most insulators.

    Args:
        efermi (float): The Fermi energy of the static run.
        cbm (float): The conduction band minimum of the static run.
        band_structure (BandStructureSymmLine): A band structure object.

    Returns:
        float: A new adjusted Fermi level.
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
# contrast when you write UNK files, the record length is written at the
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

    Attributes:
        filename (str): String of the input file (usually WAVECAR).
        vasp_type (str): String that determines VASP type the WAVECAR was generated with.
            One of 'std', 'gam', 'ncl'.
        nk (int): Number of k-points from the WAVECAR.
        nb (int): Number of bands per k-point.
        encut (float): Energy cutoff (used to define G_{cut}).
        efermi (float): Fermi energy.
        a (numpy.ndarray): Primitive lattice vectors of the cell (e.g. a_1 = self.a[0, :]).
        b (numpy.ndarray): Reciprocal lattice vectors of the cell (e.g. b_1 = self.b[0, :]).
        vol (float): The volume of the unit cell in real space.
        kpoints (numpy.ndarray): The list of k-points read from the WAVECAR file.
        band_energy (list): The list of band eigenenergies (and corresponding occupancies) for each kpoint,
            where the first index corresponds to the index of the k-point (e.g. self.band_energy[kp]).
        Gpoints (list): The list of generated G-points for each k-point (a double list), which
            are used with the coefficients for each k-point and band to recreate
            the wavefunction (e.g. self.Gpoints[kp] is the list of G-points for
            k-point kp). The G-points depend on the k-point and reciprocal lattice
            and therefore are identical for each band at the same k-point. Each
            G-point is represented by integer multipliers (e.g. assuming
            Gpoints[kp][n] == [n_1, n_2, n_3], then
            G_n = n_1*b_1 + n_2*b_2 + n_3*b_3)
        coeffs (list): The list of coefficients for each k-point and band for reconstructing the wavefunction.
            For non-spin-polarized, the first index corresponds to the kpoint and the second corresponds to the band
            (e.g. self.coeffs[kp][b] corresponds to k-point kp and band b). For spin-polarized calculations,
            the first index is for the spin. If the calculation was non-collinear, then self.coeffs[kp][b] will have
            two columns (one for each component of the spinor).

    Acknowledgments:
        This code is based upon the Fortran program, WaveTrans, written by
        R. M. Feenstra and M. Widom from the Dept. of Physics at Carnegie
        Mellon University. To see the original work, please visit:
        https://www.andrew.cmu.edu/user/feenstra/wavetrans/

    Author: Mark Turiansky
    """

    def __init__(self, filename="WAVECAR", verbose=False, precision="normal", vasp_type=None):
        """
        Information is extracted from the given WAVECAR.

        Args:
            filename (str): input file (default: WAVECAR)
            verbose (bool): determines whether processing information is shown
            precision (str): determines how fine the fft mesh is (normal or
                accurate), only the first letter matters
            vasp_type (str): determines the VASP type that is used, allowed
                values are ['std', 'gam', 'ncl'] (only first letter is required)
        """
        self.filename = filename
        valid_types = ["std", "gam", "ncl"]
        initials = {x[0] for x in valid_types}
        if not (vasp_type is None or vasp_type.lower()[0] in initials):
            raise ValueError(
                f"invalid {vasp_type=}, must be one of {valid_types} (we only check the first letter {initials})"
            )
        self.vasp_type = vasp_type

        # c = 0.26246582250210965422
        # 2m/hbar^2 in agreement with VASP
        self._C = 0.262465831
        with open(self.filename, "rb") as f:
            # read the header information
            recl, spin, rtag = np.fromfile(f, dtype=np.float64, count=3).astype(int)
            if verbose:
                print(f"{recl=}, {spin=}, {rtag=}")
            recl8 = int(recl / 8)
            self.spin = spin

            # check to make sure we have precision correct
            valid_rtags = {45200, 45210, 53300, 53310}
            if rtag not in valid_rtags:
                # note that rtag=45200 and 45210 may not work if file was actually
                # generated by old version of VASP, since that would write eigenvalues
                # and occupations in way that does not span FORTRAN records, but
                # reader below appears to assume that record boundaries can be ignored
                # (see OUTWAV vs. OUTWAV_4 in vasp fileio.F)
                raise ValueError(f"invalid {rtag=}, must be one of {valid_rtags}")

            # padding to end of fortran REC=1
            np.fromfile(f, dtype=np.float64, count=recl8 - 3)

            # extract kpoint, bands, energy, and lattice information
            self.nk, self.nb = np.fromfile(f, dtype=np.float64, count=2).astype(int)
            self.encut = np.fromfile(f, dtype=np.float64, count=1)[0]
            self.a = np.fromfile(f, dtype=np.float64, count=9).reshape((3, 3))
            self.efermi = np.fromfile(f, dtype=np.float64, count=1)[0]
            if verbose:
                print(
                    f"kpoints = {self.nk}, bands = {self.nb}, energy cutoff = {self.encut}, fermi "
                    f"energy= {self.efermi:.04f}\n"
                )
                print(f"primitive lattice vectors = \n{self.a}")

            self.vol = np.dot(self.a[0, :], np.cross(self.a[1, :], self.a[2, :]))
            if verbose:
                print(f"volume = {self.vol}\n")

            # calculate reciprocal lattice
            b = np.array(
                [
                    np.cross(self.a[1, :], self.a[2, :]),
                    np.cross(self.a[2, :], self.a[0, :]),
                    np.cross(self.a[0, :], self.a[1, :]),
                ]
            )
            b = 2 * np.pi * b / self.vol
            self.b = b
            if verbose:
                print(f"reciprocal lattice vectors = \n{b}")
                print(f"reciprocal lattice vector magnitudes = \n{np.linalg.norm(b, axis=1)}\n")

            # calculate maximum number of b vectors in each direction
            self._generate_nbmax()
            if verbose:
                print(f"max number of G values = {self._nbmax}\n\n")
            self.ng = self._nbmax * 3 if precision.lower()[0] == "n" else self._nbmax * 4

            # padding to end of fortran REC=2
            np.fromfile(f, dtype=np.float64, count=recl8 - 13)

            # reading records
            self.Gpoints = [None for _ in range(self.nk)]
            self.kpoints = []
            if spin == 2:
                self.coeffs = [[[None for i in range(self.nb)] for j in range(self.nk)] for _ in range(spin)]
                self.band_energy = [[] for _ in range(spin)]
            else:
                self.coeffs = [[None for i in range(self.nb)] for j in range(self.nk)]
                self.band_energy = []

            for ispin in range(spin):
                if verbose:
                    print(f"reading spin {ispin}")

                for ink in range(self.nk):
                    # information for this kpoint
                    nplane = int(np.fromfile(f, dtype=np.float64, count=1)[0])
                    kpoint = np.fromfile(f, dtype=np.float64, count=3)

                    if ispin == 0:
                        self.kpoints.append(kpoint)
                    else:
                        assert_allclose(self.kpoints[ink], kpoint)

                    if verbose:
                        print(f"kpoint {ink: 4} with {nplane: 5} plane waves at {kpoint}")

                    # energy and occupation information
                    enocc = np.fromfile(f, dtype=np.float64, count=3 * self.nb).reshape((self.nb, 3))
                    if spin == 2:
                        self.band_energy[ispin].append(enocc)
                    else:
                        self.band_energy.append(enocc)

                    if verbose:
                        print("enocc =\n", enocc[:, [0, 2]])

                    # padding to end of record that contains nplane, kpoints, evals and occs
                    np.fromfile(f, dtype=np.float64, count=(recl8 - 4 - 3 * self.nb) % recl8)

                    if self.vasp_type is None:
                        self.Gpoints[ink], extra_gpoints, extra_coeff_inds = self._generate_G_points(kpoint, gamma=True)
                        if len(self.Gpoints[ink]) == nplane:
                            self.vasp_type = "gam"
                        else:
                            self.Gpoints[ink], extra_gpoints, extra_coeff_inds = self._generate_G_points(
                                kpoint, gamma=False
                            )
                            self.vasp_type = "std" if len(self.Gpoints[ink]) == nplane else "ncl"

                        if verbose:
                            print(f"\ndetermined {self.vasp_type = }\n")
                    else:
                        self.Gpoints[ink], extra_gpoints, extra_coeff_inds = self._generate_G_points(
                            kpoint, gamma=self.vasp_type.lower()[0] == "g"
                        )

                    if len(self.Gpoints[ink]) != nplane and 2 * len(self.Gpoints[ink]) != nplane:
                        raise ValueError(
                            f"Incorrect {vasp_type=}. Please open an issue if you are certain this WAVECAR"
                            " was generated with the given vasp_type."
                        )

                    self.Gpoints[ink] = np.array(self.Gpoints[ink] + extra_gpoints, dtype=np.float64)

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
                            self.coeffs[ispin][ink][inb] = np.array(list(data) + extra_coeffs, dtype=np.complex64)
                        else:
                            self.coeffs[ink][inb] = np.array(list(data) + extra_coeffs, dtype=np.complex128)

                        if self.vasp_type.lower()[0] == "n":
                            self.coeffs[ink][inb].shape = (2, nplane // 2)

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
        sphi123 = np.dot(b[2, :], np.cross(b[0, :], b[1, :])) / (bmag[2] * np.linalg.norm(np.cross(b[0, :], b[1, :])))
        nbmaxA = np.sqrt(self.encut * self._C) / bmag
        nbmaxA[0] /= np.abs(np.sin(phi12))
        nbmaxA[1] /= np.abs(np.sin(phi12))
        nbmaxA[2] /= np.abs(sphi123)
        nbmaxA += 1

        phi13 = np.arccos(np.dot(b[0, :], b[2, :]) / (bmag[0] * bmag[2]))
        sphi123 = np.dot(b[1, :], np.cross(b[0, :], b[2, :])) / (bmag[1] * np.linalg.norm(np.cross(b[0, :], b[2, :])))
        nbmaxB = np.sqrt(self.encut * self._C) / bmag
        nbmaxB[0] /= np.abs(np.sin(phi13))
        nbmaxB[1] /= np.abs(sphi123)
        nbmaxB[2] /= np.abs(np.sin(phi13))
        nbmaxB += 1

        phi23 = np.arccos(np.dot(b[1, :], b[2, :]) / (bmag[1] * bmag[2]))
        sphi123 = np.dot(b[0, :], np.cross(b[1, :], b[2, :])) / (bmag[0] * np.linalg.norm(np.cross(b[1, :], b[2, :])))
        nbmaxC = np.sqrt(self.encut * self._C) / bmag
        nbmaxC[0] /= np.abs(sphi123)
        nbmaxC[1] /= np.abs(np.sin(phi23))
        nbmaxC[2] /= np.abs(np.sin(phi23))
        nbmaxC += 1

        self._nbmax = np.max([nbmaxA, nbmaxB, nbmaxC], axis=0).astype(int)

    def _generate_G_points(self, kpoint: np.ndarray, gamma: bool = False) -> tuple[list, list, list]:
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
        kmax = self._nbmax[0] + 1 if gamma else 2 * self._nbmax[0] + 1

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
                    E = g**2 / self._C
                    if self.encut > E:
                        gpoints.append(G)
                        if gamma and (k1, j2, i3) != (0, 0, 0):
                            extra_gpoints.append(-G)
                            extra_coeff_inds.append(G_ind)
                        G_ind += 1
        return gpoints, extra_gpoints, extra_coeff_inds

    def evaluate_wavefunc(self, kpoint: int, band: int, r: np.ndarray, spin: int = 0, spinor: int = 0) -> np.complex64:
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
            kpoint (int): the index of the kpoint where the wavefunction will be evaluated
            band (int): the index of the band where the wavefunction will be evaluated
            r (np.array): the position where the wavefunction will be evaluated
            spin (int): spin index for the desired wavefunction (only for
                ISPIN = 2, default = 0)
            spinor (int): component of the spinor that is evaluated (only used
                if vasp_type == 'ncl')

        Returns:
            a complex value corresponding to the evaluation of the wavefunction
        """
        v = self.Gpoints[kpoint] + self.kpoints[kpoint]
        u = np.dot(np.dot(v, self.b), r)
        if self.vasp_type.lower()[0] == "n":
            c = self.coeffs[kpoint][band][spinor, :]
        elif self.spin == 2:
            c = self.coeffs[spin][kpoint][band]
        else:
            c = self.coeffs[kpoint][band]
        return np.sum(np.dot(c, np.exp(1j * u, dtype=np.complex64))) / np.sqrt(self.vol)

    def fft_mesh(self, kpoint: int, band: int, spin: int = 0, spinor: int = 0, shift: bool = True) -> np.ndarray:
        """
        Places the coefficients of a wavefunction onto an fft mesh.

        Once the mesh has been obtained, a discrete fourier transform can be
        used to obtain real-space evaluation of the wavefunction. The output
        of this function can be passed directly to numpy's fft function. For
        example:

            mesh = Wavecar('WAVECAR').fft_mesh(kpoint, band)
            evals = np.fft.ifftn(mesh)

        Args:
            kpoint (int): the index of the kpoint where the wavefunction will be evaluated
            band (int): the index of the band where the wavefunction will be evaluated
            spin (int): the spin of the wavefunction for the desired
                wavefunction (only for ISPIN = 2, default = 0)
            spinor (int): component of the spinor that is evaluated (only used
                if vasp_type == 'ncl')
            shift (bool): determines if the zero frequency coefficient is
                placed at index (0, 0, 0) or centered

        Returns:
            a numpy ndarray representing the 3D mesh of coefficients
        """
        if self.vasp_type.lower()[0] == "n":
            tcoeffs = self.coeffs[kpoint][band][spinor, :]
        elif self.spin == 2:
            tcoeffs = self.coeffs[spin][kpoint][band]
        else:
            tcoeffs = self.coeffs[kpoint][band]

        mesh = np.zeros(tuple(self.ng), dtype=np.complex_)
        for gp, coeff in zip(self.Gpoints[kpoint], tcoeffs):
            t = tuple(gp.astype(int) + (self.ng / 2).astype(int))
            mesh[t] = coeff

        if shift:
            return np.fft.ifftshift(mesh)
        return mesh

    def get_parchg(
        self,
        poscar: Poscar,
        kpoint: int,
        band: int,
        spin: int | None = None,
        spinor: int | None = None,
        phase: bool = False,
        scale: int = 2,
    ) -> Chgcar:
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
        if phase and not np.all(self.kpoints[kpoint] == 0.0):
            warnings.warn("phase is True should only be used for the Gamma kpoint! I hope you know what you're doing!")

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
                data["total"] = den
            else:
                wfr = np.fft.ifftn(self.fft_mesh(kpoint, band, spin=0)) * N
                denup = np.abs(np.conj(wfr) * wfr)
                wfr = np.fft.ifftn(self.fft_mesh(kpoint, band, spin=1)) * N
                dendn = np.abs(np.conj(wfr) * wfr)
                data["total"] = denup + dendn
                data["diff"] = denup - dendn
        else:
            if spinor is not None:
                wfr = np.fft.ifftn(self.fft_mesh(kpoint, band, spinor=spinor)) * N
                den = np.abs(np.conj(wfr) * wfr)
            else:
                wfr = np.fft.ifftn(self.fft_mesh(kpoint, band, spinor=0)) * N
                wfr_t = np.fft.ifftn(self.fft_mesh(kpoint, band, spinor=1)) * N
                den = np.abs(np.conj(wfr) * wfr)
                den += np.abs(np.conj(wfr_t) * wfr_t)

            if phase and not (self.vasp_type.lower()[0] == "n" and spinor is None):
                den = np.sign(np.real(wfr)) * den
            data["total"] = den

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
            raise ValueError("invalid directory")

        N = np.prod(self.ng)
        for ik in range(self.nk):
            fname = f"UNK{ik+1:05d}."
            if self.vasp_type.lower()[0] == "n":
                data = np.empty((self.nb, 2, *self.ng), dtype=np.complex128)
                for ib in range(self.nb):
                    data[ib, 0, :, :, :] = np.fft.ifftn(self.fft_mesh(ik, ib, spinor=0)) * N
                    data[ib, 1, :, :, :] = np.fft.ifftn(self.fft_mesh(ik, ib, spinor=1)) * N
                Unk(ik + 1, data).write_file(str(out_dir / (fname + "NC")))
            else:
                data = np.empty((self.nb, *self.ng), dtype=np.complex128)
                for ispin in range(self.spin):
                    for ib in range(self.nb):
                        data[ib, :, :, :] = np.fft.ifftn(self.fft_mesh(ik, ib, spin=ispin)) * N
                    Unk(ik + 1, data).write_file(str(out_dir / f"{fname}{ispin+1}"))


class Eigenval:
    """
    Object for reading EIGENVAL file.

    Attributes:
        filename (str): String containing input filename.
        occu_tol (float): Tolerance for determining occupation in band properties.
        ispin (int): Spin polarization tag.
        nelect (int): Number of electrons.
        nkpt (int): Number of kpoints.
        nbands (int): Number of bands.
        kpoints (list): List of kpoints.
        kpoints_weights (list): Weights of each kpoint in the BZ, should sum to 1.
        eigenvalues (dict): Eigenvalues as a dict of {(spin): np.ndarray(shape=(nkpt, nbands, 2))}.
            This representation is based on actual ordering in VASP and is meant as an intermediate representation
            to be converted into proper objects. The kpoint index is 0-based (unlike the 1-based indexing in VASP).
    """

    def __init__(self, filename, occu_tol=1e-8, separate_spins=False):
        """
        Reads input from filename to construct Eigenval object.

        Args:
            filename (str):     filename of EIGENVAL to read in
            occu_tol (float):   tolerance for determining band gap
            separate_spins (bool):   whether the band gap, CBM, and VBM should be
                reported for each individual spin channel. Defaults to False,
                which computes the eigenvalue band properties independent of
                the spin orientation. If True, the calculation must be spin-polarized.

        Returns:
            a pymatgen.io.vasp.outputs.Eigenval object
        """
        self.filename = filename
        self.occu_tol = occu_tol
        self.separate_spins = separate_spins

        with zopen(filename, "r") as f:
            self.ispin = int(f.readline().split()[-1])

            # useless header information
            for _ in range(4):
                f.readline()

            self.nelect, self.nkpt, self.nbands = list(map(int, f.readline().split()))

            self.kpoints = []
            self.kpoints_weights = []
            if self.ispin == 2:
                self.eigenvalues = {
                    Spin.up: np.zeros((self.nkpt, self.nbands, 2)),
                    Spin.down: np.zeros((self.nkpt, self.nbands, 2)),
                }
            else:
                self.eigenvalues = {Spin.up: np.zeros((self.nkpt, self.nbands, 2))}

            ikpt = -1
            for line in f:
                if re.search(r"(\s+[\-+0-9eE.]+){4}", str(line)):
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
        (band gap, cbm, vbm, is_band_gap_direct). In the case of separate_spins=True,
        the band gap, cbm, vbm, and is_band_gap_direct are each lists of length 2,
        with index 0 representing the spin-up channel and index 1 representing
        the spin-down channel.
        """
        vbm = -float("inf")
        vbm_kpoint = None
        cbm = float("inf")
        cbm_kpoint = None
        vbm_spins = []
        vbm_spins_kpoints = []
        cbm_spins = []
        cbm_spins_kpoints = []
        if self.separate_spins and len(self.eigenvalues) != 2:
            raise ValueError("The separate_spins flag can only be True if ISPIN = 2")

        for d in self.eigenvalues.values():
            if self.separate_spins:
                vbm = -float("inf")
                cbm = float("inf")
            for k, val in enumerate(d):
                for eigenval, occu in val:
                    if occu > self.occu_tol and eigenval > vbm:
                        vbm = eigenval
                        vbm_kpoint = k
                    elif occu <= self.occu_tol and eigenval < cbm:
                        cbm = eigenval
                        cbm_kpoint = k
            if self.separate_spins:
                vbm_spins.append(vbm)
                vbm_spins_kpoints.append(vbm_kpoint)
                cbm_spins.append(cbm)
                cbm_spins_kpoints.append(cbm_kpoint)
        if self.separate_spins:
            return (
                [max(cbm_spins[0] - vbm_spins[0], 0), max(cbm_spins[1] - vbm_spins[1], 0)],
                [cbm_spins[0], cbm_spins[1]],
                [vbm_spins[0], vbm_spins[1]],
                [vbm_spins_kpoints[0] == cbm_spins_kpoints[0], vbm_spins_kpoints[1] == cbm_spins_kpoints[1]],
            )
        return max(cbm - vbm, 0), cbm, vbm, vbm_kpoint == cbm_kpoint


@dataclass
class Waveder(MSONable):
    """Representation of the WAVEDER file.

    The LOPTICS tag produces a WAVEDER file which contains the derivative of the orbitals with respect to k.
    Since the data is complex, we need to split it into the real and imaginary parts for JSON serialization.

    Note:
        The way that VASP writes the WAVEDER and WAVEDERF has slightly different logic when indexing the bands.
        This results in the formatted WAVDERF only indexing between filled bands. (i.e. all the matrix elements
        are between the states i=1:8 and j=1:8 in a two atom Si calculation, which is likely a VASP bug).
        As such, it is recommended to used the hidden ``LVEL=.True.`` flag in VASP which will force indexing over
        all bands.

    The order of the indices of the data are:
    [
        band index1,
        band index2,
        kpoint index,
        spin index,
        cartesian direction,
    ]

    Attributes:
        cder_real: Real part of the derivative of the orbitals with respect to k.
        cder_imag: Imaginary part of the derivative of the orbitals with respect to k.

    Author: Miguel Dias Costa, Kamal Choudhary, Jimmy-Xuan Shen
    """

    cder_real: np.ndarray
    cder_imag: np.ndarray

    @classmethod
    def from_formatted(cls, filename):
        """Reads the WAVEDERF file and returns a Waveder object.

        Note: This file is only produced when LOPTICS is true AND vasp has been
        recompiled after uncommenting the line that calls
        WRT_CDER_BETWEEN_STATES_FORMATTED in linear_optics.F
        It is recommended to use `from_binary` instead since the binary file is
        much smaller and contains the same information.

        Args:
            filename (str): The name of the WAVEDER file.

        Returns:
            A Waveder object.
        """
        with zopen(filename, "rt") as f:
            nspin, nkpts, nbands = f.readline().split()
        # 1 and 4 are the eigenvalues of the bands (this data is missing in the WAVEDER file)
        # 6:12 are the complex matrix elements in each cartesian direction.
        data = np.loadtxt(filename, skiprows=1, usecols=(1, 4, 6, 7, 8, 9, 10, 11))
        data = data.reshape(int(nspin), int(nkpts), int(nbands), int(nbands), 8)  # slowest to fastest
        cder_real = data[:, :, :, :, 2::2]
        cder_imag = data[:, :, :, :, 3::2]
        # change to [band1, band2, kpt, spin, cartesian]
        cder_real = np.transpose(cder_real, (2, 3, 1, 0, 4))
        cder_imag = np.transpose(cder_imag, (2, 3, 1, 0, 4))
        # TODO add eigenvalues?
        return cls(cder_real, cder_imag)

    @classmethod
    def from_binary(cls, filename, data_type="complex64"):
        """Read the WAVEDER file and returns a Waveder object.

        Args:
            filename: Name of file containing WAVEDER.
            data_type: Data type of the WAVEDER file. Default is complex64.
                If the file was generated with the "gamma" version of VASP,
                the data type can be either "float64" or "float32".

        Returns:
            Waveder object.
        """
        with open(filename, "rb") as fp:

            def readData(dtype):
                """Read records from Fortran binary file and convert to np.array of given dtype."""
                data = b""
                while True:
                    prefix = np.fromfile(fp, dtype=np.int32, count=1)[0]
                    data += fp.read(abs(prefix))
                    suffix = np.fromfile(fp, dtype=np.int32, count=1)[0]
                    if abs(prefix) - abs(suffix):
                        raise RuntimeError(
                            f"Read wrong amount of bytes.\nExpected: {prefix}, read: {len(data)}, suffix: {suffix}."
                        )
                    if prefix > 0:
                        break
                return np.frombuffer(data, dtype=dtype)

            nbands, nelect, nk, ispin = readData(np.int32)
            _ = readData(np.float_)  # nodes_in_dielectric_function
            _ = readData(np.float_)  # wplasmon
            me_datatype = np.dtype(data_type)
            cder = readData(me_datatype)

            cder_data = cder.reshape((3, ispin, nk, nelect, nbands)).T
            return cls(cder_data.real, cder_data.imag)

    @property
    def cder(self):
        """Return the complex derivative of the orbitals with respect to k."""
        if self.cder_real.shape[0] != self.cder_real.shape[1]:  # pragma: no cover
            warnings.warn(
                "Not all band pairs are present in the WAVEDER file."
                "If you want to get all the matrix elements set LVEL=.True. in the INCAR."
            )
        return self.cder_real + 1j * self.cder_imag

    @property
    def nspin(self):
        """Returns the number of spin channels."""
        return self.cder_real.shape[3]

    @property
    def nkpoints(self):
        """Returns the number of k-points."""
        return self.cder_real.shape[2]

    @property
    def nbands(self):
        """Returns the number of bands."""
        return self.cder_real.shape[0]

    def get_orbital_derivative_between_states(self, band_i, band_j, kpoint, spin, cart_dir):
        """
        Method returning a value
        between bands band_i and band_j for k-point index, spin-channel and Cartesian direction.

        Args:
            band_i (int): Index of band i
            band_j (int): Index of band j
            kpoint (int): Index of k-point
            spin   (int): Index of spin-channel (0 or 1)
            cart_dir (int): Index of Cartesian direction (0,1,2)

        Returns:
            a float value
        """
        return self.cder[band_i, band_j, kpoint, spin, cart_dir]


@dataclass
class WSWQ(MSONable):
    r"""
    Class for reading a WSWQ file.
    The WSWQ file is used to calculation the wave function overlaps between
        - W: Wavefunctions in the currenct directory's WAVECAR file
        - WQ: Wavefunctions stored in a filed named the WAVECAR.qqq.

    The overlap is computed using the overlap operator S
    which make the PAW wavefunctions orthogonormal:
        <W_k,m| S | W_k,n> = \delta_{mn}

    The WSWQ file contains matrix elements of the overlap operator S evaluated
    between the planewave wavefunctions W and WQ:
        COVL_k,mn = < W_s,k,m | S | WQ_s,k,n >

    The indices of WSWQ.data are:
        [spin][kpoint][band_i][band_j]

    Attributes:
        nspin: Number of spin channels
        nkpoints: Number of k-points
        nbands: Number of bands
        me_real: Real part of the overlap matrix elements
        me_imag: Imaginary part of the overlap matrix elements
    """

    nspin: int
    nkpoints: int
    nbands: int
    me_real: np.ndarray
    me_imag: np.ndarray

    @property
    def data(self):
        """Complex overlap matrix."""
        return self.me_real + 1j * self.me_imag

    @classmethod
    def from_file(cls, filename: str) -> WSWQ:
        """Constructs a WSWQ object from a file.

        Args:
            filename (str): Name of WSWQ file.

        Returns:
            WSWQ object.
        """
        # Grab the spin and kpoint values from lines containing "spin"
        spin_res = regrep(
            filename,
            {"spin": r"spin\s*=\s*(\d+)\s?\,\s?kpoint\s*=\s*(\d+)"},
            reverse=True,
            terminate_on_match=True,
            postprocess=int,
        )["spin"]
        (nspin, nkpoints), _ = spin_res[0]
        # Grab the index values from the parts with "i    =" and "j    ="
        ij_res = regrep(
            filename,
            {"ij": r"i\s*=\s*(\d+)\s?\,\s?j\s*=\s*(\d+)"},
            reverse=True,
            terminate_on_match=True,
            postprocess=int,
        )["ij"]
        (nbands, _), _ = ij_res[0]
        # Grab the matrix elements from the parts after the ":"
        data_res = regrep(
            filename,
            {"data": r"\:\s*([-+]?\d*\.\d+)\s+([-+]?\d*\.\d+)"},
            reverse=False,
            terminate_on_match=False,
            postprocess=float,
        )["data"]
        assert len(data_res) == nspin * nkpoints * nbands * nbands

        data = np.array([complex(real_part, img_part) for (real_part, img_part), _ in data_res])

        # NOTE: loop order (slow->fast) spin -> kpoint -> j -> i
        data = data.reshape((nspin, nkpoints, nbands, nbands))
        data = np.swapaxes(data, 2, 3)  # swap i and j
        return cls(nspin=nspin, nkpoints=nkpoints, nbands=nbands, me_real=np.real(data), me_imag=np.imag(data))


class UnconvergedVASPWarning(Warning):
    """Warning for unconverged vasp run."""
