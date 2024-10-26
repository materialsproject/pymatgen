"""Classes for reading/manipulating/writing VASP output files."""

from __future__ import annotations

import collections
import itertools
import math
import os
import re
import typing
import warnings
import xml.etree.ElementTree as ET
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import datetime, timezone
from glob import glob
from io import StringIO
from pathlib import Path
from typing import TYPE_CHECKING, cast

import numpy as np
from monty.io import reverse_readfile, zopen
from monty.json import MSONable, jsanitize
from monty.os.path import zpath
from monty.re import regrep
from numpy.testing import assert_allclose
from tqdm import tqdm

from pymatgen.core import Composition, Element, Lattice, Structure
from pymatgen.core.trajectory import Trajectory
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
from pymatgen.util.typing import Kpoint, Tuple3Floats, Vector3D

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, Literal

    # Avoid name conflict with pymatgen.core.Element
    from xml.etree.ElementTree import Element as XML_Element

    from numpy.typing import NDArray
    from typing_extensions import Self

    from pymatgen.util.typing import PathLike


def _parse_parameters(val_type: str, val: str) -> bool | str | float | int:
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


def _parse_v_parameters(
    val_type: str,
    val: str,
    filename: PathLike,
    param_name: str,
) -> list[bool] | list[float] | list[int] | list[str]:
    """
    Helper function to convert a Vasprun array-type parameter into
    the proper type. Boolean, int and float types are converted.

    Args:
        val_type: Value type parsed from vasprun.xml.
        val: Actual string value parsed for vasprun.xml.
        filename: Fullpath of vasprun.xml. Used for robust error handling.
            e.g. if vasprun.xml contains *** for some Incar parameters,
            the code will try to read from an INCAR file present in the same
            directory.
        param_name: Name of parameter.

    Returns:
        Parsed value.
    """
    err = ValueError("Error in parsing vasprun.xml")

    if val_type == "logical":
        return [i == "T" for i in val.split()]

    if val_type == "string":
        return val.split()

    if val_type == "int":
        try:
            ints = [int(i) for i in val.split()]
        except ValueError as exc:
            # Fix an error in vasprun where
            # sometimes LDAUL/J is displayed as 2****
            ints = _parse_from_incar(filename, param_name)
            if ints is None:
                raise err from exc

        return ints

    try:
        floats = [float(i) for i in val.split()]
    except ValueError as exc:
        # Fix an error in vasprun where
        # sometimes MAGMOM is displayed as 2****
        floats = _parse_from_incar(filename, param_name)
        if floats is None:
            raise err from exc

    return floats


def _parse_vasp_array(elem) -> list[list[float]]:
    if elem.get("type") == "logical":
        return [[i == "T" for i in v.text.split()] for v in elem]
    return [[_vasprun_float(i) for i in v.text.split()] for v in elem]


def _parse_from_incar(filename: PathLike, key: str) -> Any:
    """Helper function to parse a parameter from the INCAR."""
    dirname = os.path.dirname(filename)
    for fn in os.listdir(dirname):
        if re.search("INCAR", fn):
            warnings.warn(f"INCAR found. Using {key} from INCAR.", stacklevel=2)
            incar = Incar.from_file(os.path.join(dirname, fn))
            return incar.get(key)
    return None


def _vasprun_float(flt: float | str) -> float:
    """
    Large numbers are often represented as ********* in the vasprun.
    This function parses these values as np.nan.
    """
    try:
        return float(flt)

    except ValueError:
        flt = cast(str, flt)
        _flt: str = flt.strip()
        if _flt == "*" * len(_flt):
            warnings.warn("Float overflow (*******) encountered in vasprun")
            return np.nan
        raise


@dataclass
class KpointOptProps:
    """Simple container class to store KPOINTS_OPT data in a separate namespace. Used by Vasprun."""

    tdos: Dos | None = None
    idos: Dos | None = None
    pdos: list | None = None
    efermi: float | None = None
    eigenvalues: dict | None = None
    projected_eigenvalues: dict | None = None
    projected_magnetisation: np.ndarray | None = None
    kpoints: Kpoints | None = None
    actual_kpoints: list | None = None
    actual_kpoints_weights: list | None = None
    dos_has_errors: bool | None = None


class Vasprun(MSONable):
    """
    Vastly improved cElementTree-based parser for vasprun.xml files. Uses
    iterparse to support incremental parsing of large files.
    Speedup over Dom is at least 2x for smallish files (~1 Mb) to orders of
    magnitude for larger files (~10 Mb).

    **VASP results**

    Attributes:
        ionic_steps (list): All ionic steps in the run as a list of {"structure": structure at end of run,
            "electronic_steps": {All electronic step data in vasprun file}, "stresses": stress matrix}.
        tdos (Dos): Total dos calculated at the end of run. Note that this is rounded to 4 decimal
            places by VASP.
        idos (Dos): Integrated dos calculated at the end of run. Rounded to 4 decimal places by VASP.
        pdos (list): List of list of PDos objects. Access as pdos[atomindex][orbitalindex].
        efermi (float): Fermi energy.
        eigenvalues (dict): Final eigenvalues as a dict of {(spin, kpoint index):[[eigenvalue, occu]]}.
            The kpoint index is 0-based (unlike the 1-based indexing in VASP).
        projected_eigenvalues (dict): Final projected eigenvalues as a dict of {spin: nd-array}.
            To access a particular value, you need to do
            Vasprun.projected_eigenvalues[spin][kpoint index][band index][atom index][orbital_index].
            The kpoint, band and atom indices are 0-based (unlike the 1-based indexing in VASP).
        projected_magnetisation (np.array): Final projected magnetization as a numpy array with the
            shape (nkpoints, nbands, natoms, norbitals, 3). Where the last axis is the contribution in the
            3 Cartesian directions. This attribute is only set if spin-orbit coupling (LSORBIT = True) or
            non-collinear magnetism (LNONCOLLINEAR = True) is turned on in the INCAR.
        dielectric_data (dict): Dictionary, with the tag comment as key, containing other variants of
            the real and imaginary part of the dielectric constant (e.g., computed by RPA) in function of
            the energy (frequency). Optical properties (e.g. absorption coefficient) can be obtained through this.
            The data is given as a tuple of 3 values containing each of them the energy, the real part tensor,
            and the imaginary part tensor ([energies],[[real_partxx,real_partyy,real_partzz,real_partxy,
            real_partyz,real_partxz]],[[imag_partxx,imag_partyy,imag_partzz,imag_partxy, imag_partyz, imag_partxz]]).
            The data can be the current, density or freq_dependent (BSE) dielectric data.
        nionic_steps (int): The total number of ionic steps. This number is always equal to the total number
            of steps in the actual run even if ionic_step_skip is used.
        force_constants (np.array): Force constants computed in phonon DFPT run(IBRION = 8).
            The data is a 4D numpy array of shape (natoms, natoms, 3, 3).
        normalmode_eigenvals (np.array): Normal mode frequencies. 1D numpy array of size 3*natoms.
        normalmode_eigenvecs (np.array): Normal mode eigen vectors. 3D numpy array of shape (3*natoms, natoms, 3).
        md_data (list): Available only for ML MD runs, i.e., INCAR with ML_LMLFF = .TRUE. md_data is a list of
            dict with the following format: [{'energy': {'e_0_energy': -525.07195568, 'e_fr_energy': -525.07195568,
            'e_wo_entrp': -525.07195568, 'kinetic': 3.17809233, 'lattice kinetic': 0.0, 'nosekinetic': 1.323e-5,
            'nosepot': 0.0, 'total': -521.89385012}, 'forces': [[0.17677989, 0.48309874, 1.85806696], ...],
            'structure': Structure object}].
        incar (Incar): Incar object for parameters specified in INCAR file.
        parameters (Incar): Incar object with parameters that VASP actually used, including all defaults.
        kpoints (Kpoints): Kpoints object for KPOINTS specified in run.
        actual_kpoints (list): List of actual kpoints, e.g. [[0.25, 0.125, 0.08333333], [-0.25, 0.125, 0.08333333],
            [0.25, 0.375, 0.08333333], ....].
        actual_kpoints_weights (list): List of kpoint weights, e.g. [0.04166667, 0.04166667, 0.04166667, 0.04166667,
            0.04166667, ....].
        atomic_symbols (list): List of atomic symbols, e.g. ["Li", "Fe", "Fe", "P", "P", "P"].
        potcar_symbols (list): List of POTCAR symbols. e.g. ["PAW_PBE Li 17Jan2003", "PAW_PBE Fe 06Sep2000", ..].
        kpoints_opt_props (object): Object whose attributes are the data from KPOINTS_OPT (if present,
            else None). Attributes of the same name have the same format and meaning as Vasprun (or they are
            None if absent). Attributes are: tdos, idos, pdos, efermi, eigenvalues, projected_eigenvalues,
            projected magnetisation, kpoints, actual_kpoints, actual_kpoints_weights, dos_has_errors.

    Author: Shyue Ping Ong
    """

    def __init__(
        self,
        filename: PathLike,
        ionic_step_skip: int | None = None,
        ionic_step_offset: int = 0,
        parse_dos: bool = True,
        parse_eigen: bool = True,
        parse_projected_eigen: bool = False,
        parse_potcar_file: PathLike | bool = True,
        occu_tol: float = 1e-8,
        separate_spins: bool = False,
        exception_on_bad_xml: bool = True,
        ignore_dielectric: bool = False,
    ) -> None:
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
                Note that the DOS output from VASP is rounded to 4 decimal places,
                which can give some slight inaccuracies.
            parse_eigen (bool): Whether to parse the eigenvalues. Defaults to
                True. Set to False to shave off significant time from the
                parsing if you are not interested in getting those data.
            parse_projected_eigen (bool): Whether to parse the projected
                eigenvalues and magnetization. Defaults to False. Set to True to obtain
                projected eigenvalues and magnetization. **Note that this can take an
                extreme amount of time and memory.** So use this wisely.
            parse_potcar_file (bool | PathLike): Whether to parse the potcar file to read
                the potcar hashes for the potcar_spec attribute. Defaults to True,
                where no hashes will be determined and the potcar_spec dictionaries
                will read {"symbol": ElSymbol, "hash": None}. By Default, looks in
                the same directory as the vasprun.xml, with same extensions as
                Vasprun.xml. If a path is provided, look at that path.
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
            ignore_dielectric (bool): Whether to ignore the parsing errors related to the dielectric function. This
                is because the dielectric function can usually be parsed from the OUTCAR instead.
        """
        self.filename = filename
        self.ionic_step_skip = ionic_step_skip
        self.ionic_step_offset = ionic_step_offset
        self.occu_tol = occu_tol
        self.separate_spins = separate_spins
        self.exception_on_bad_xml = exception_on_bad_xml
        self.ignore_dielectric = ignore_dielectric

        with zopen(filename, mode="rt") as file:
            if ionic_step_skip or ionic_step_offset:
                # Remove parts of the xml file and parse the string
                content: str = file.read()
                steps: list[str] = content.split("<calculation>")

                # The text before the first <calculation> is the preamble!
                preamble: str = steps.pop(0)
                self.nionic_steps: int = len(steps)
                new_steps = steps[ionic_step_offset :: int(ionic_step_skip or 1)]

                # Add the tailing information in the last step from the run
                to_parse: str = "<calculation>".join(new_steps)
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
                    file,
                    parse_dos=parse_dos,
                    parse_eigen=parse_eigen,
                    parse_projected_eigen=parse_projected_eigen,
                )
                self.nionic_steps = len(self.ionic_steps)

            if parse_potcar_file:
                self.update_potcar_spec(parse_potcar_file)
                self.update_charge_from_potcar(parse_potcar_file)

        if self.incar.get("ALGO") not in {"Chi", "Bse"} and not self.converged and self.parameters.get("IBRION") != 0:
            msg = f"{filename} is an unconverged VASP run.\n"
            msg += f"Electronic convergence reached: {self.converged_electronic}.\n"
            msg += f"Ionic convergence reached: {self.converged_ionic}."
            warnings.warn(msg, UnconvergedVASPWarning)

    def _parse(
        self,
        stream,
        parse_dos: bool,
        parse_eigen: bool,
        parse_projected_eigen: bool,
    ) -> None:
        self.efermi: float | None = None
        self.eigenvalues: dict[Any, NDArray] | None = None
        self.projected_eigenvalues: dict[Any, NDArray] | None = None
        self.projected_magnetisation: NDArray | None = None
        self.dielectric_data: dict[str, tuple] = {}
        self.other_dielectric: dict[str, tuple] = {}
        self.incar: Incar = {}
        self.kpoints_opt_props: KpointOptProps | None = None
        ionic_steps: list[dict[str, Any]] = []

        md_data: list[dict] = []
        parsed_header: bool = False
        in_kpoints_opt: bool = False
        try:
            # When parsing XML, start tags tell us when we have entered a block
            # while end tags are when we have actually read the data.
            # To know if a particular tag is nested within another block,
            # we have to read the start tags (otherwise we will only learn of
            # the nesting after we have left the data behind).
            # When parsing KPOINTS_OPT data, some of the tags in vasprun.xml
            # can only be distinguished from their regular counterparts by
            # whether they are nested within another block. This is why we
            # must read both start and end tags and have flags to tell us
            # when we have entered or left a block. (2024-01-26)
            for event, elem in ET.iterparse(stream, events=["start", "end"]):
                tag = elem.tag
                if event == "start":
                    # The start event tells us when we have entered blocks
                    if tag == "calculation":
                        parsed_header = True
                    elif tag in ("eigenvalues_kpoints_opt", "projected_kpoints_opt"):
                        in_kpoints_opt = True

                else:  # event == "end":
                    # The end event happens when we have read a block, so have
                    # its data.
                    if not parsed_header:
                        if tag == "generator":
                            self.generator = self._parse_params(elem)
                        elif tag == "incar":
                            self.incar = self._parse_params(elem)
                        elif tag == "kpoints":
                            if not hasattr(self, "kpoints"):
                                (
                                    self.kpoints,
                                    self.actual_kpoints,
                                    self.actual_kpoints_weights,
                                ) = self._parse_kpoints(elem)
                        elif tag == "parameters":
                            self.parameters = self._parse_params(elem)
                        elif tag == "structure" and elem.attrib.get("name") == "initialpos":
                            self.initial_structure = self._parse_structure(elem)
                            self.final_structure = self.initial_structure
                        elif tag == "atominfo":
                            self.atomic_symbols, self.potcar_symbols = self._parse_atominfo(elem)
                            self.potcar_spec: list[dict] = [
                                {"titel": titel, "hash": None, "summary_stats": {}} for titel in self.potcar_symbols
                            ]

                    if tag == "calculation":
                        parsed_header = True
                        if not self.parameters.get("LCHIMAG", False):
                            ionic_steps.append(self._parse_ionic_step(elem))
                        else:
                            ionic_steps.extend(self._parse_chemical_shielding(elem))

                    elif parse_dos and tag == "dos":
                        if elem.get("comment") == "kpoints_opt":
                            kpoints_opt_props = self.kpoints_opt_props = self.kpoints_opt_props or KpointOptProps()
                            try:
                                (
                                    kpoints_opt_props.tdos,
                                    kpoints_opt_props.idos,
                                    kpoints_opt_props.pdos,
                                ) = self._parse_dos(elem)
                                kpoints_opt_props.efermi = kpoints_opt_props.tdos.efermi
                                kpoints_opt_props.dos_has_errors = False
                            except Exception:
                                kpoints_opt_props.dos_has_errors = True
                        else:
                            try:
                                self.tdos, self.idos, self.pdos = self._parse_dos(elem)
                                self.efermi = self.tdos.efermi
                                self.dos_has_errors = False
                            except Exception:
                                self.dos_has_errors = True

                    elif parse_eigen and tag == "eigenvalues" and not in_kpoints_opt:
                        self.eigenvalues = self._parse_eigen(elem)

                    elif parse_projected_eigen and tag == "projected" and not in_kpoints_opt:
                        self.projected_eigenvalues, self.projected_magnetisation = self._parse_projected_eigen(elem)

                    elif tag in ("eigenvalues_kpoints_opt", "projected_kpoints_opt"):
                        in_kpoints_opt = False
                        if self.kpoints_opt_props is None:
                            self.kpoints_opt_props = KpointOptProps()
                        if parse_eigen:
                            # projected_kpoints_opt includes occupation information whereas
                            # eigenvalues_kpoints_opt doesn't.
                            self.kpoints_opt_props.eigenvalues = self._parse_eigen(elem.find("eigenvalues"))
                        if tag == "eigenvalues_kpoints_opt":
                            (
                                self.kpoints_opt_props.kpoints,
                                self.kpoints_opt_props.actual_kpoints,
                                self.kpoints_opt_props.actual_kpoints_weights,
                            ) = self._parse_kpoints(elem.find("kpoints"))
                        elif parse_projected_eigen:  # and tag == "projected_kpoints_opt": (implied)
                            (
                                self.kpoints_opt_props.projected_eigenvalues,
                                self.kpoints_opt_props.projected_magnetisation,
                            ) = self._parse_projected_eigen(elem)

                    elif tag == "dielectricfunction":
                        if (
                            "comment" not in elem.attrib
                            or elem.attrib["comment"]
                            == "INVERSE MACROSCOPIC DIELECTRIC TENSOR (including local field effects in RPA (Hartree))"
                        ):
                            if self.incar.get("ALGO", "Normal").upper() == "BSE":
                                self.dielectric_data["freq_dependent"] = self._parse_diel(elem)
                            elif "density" not in self.dielectric_data:
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
                        self.optical_transition = np.array(_parse_vasp_array(elem))

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
                                self.force_constants[ii, jj] = hessian[ii * 3 : (ii + 1) * 3, jj * 3 : (jj + 1) * 3]  # type: ignore[call-overload]
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
                            md_data[-1]["forces"] = _parse_vasp_array(elem)
                        elif tag == "varray" and elem.attrib.get("name") == "stress":
                            md_data[-1]["stress"] = _parse_vasp_array(elem)
                        elif tag == "energy":
                            d = {i.attrib["name"]: float(i.text) for i in elem.findall("i")}
                            if "kinetic" in d:
                                md_data[-1]["energy"] = {i.attrib["name"]: float(i.text) for i in elem.findall("i")}

        except ET.ParseError:
            if self.exception_on_bad_xml:
                raise
            warnings.warn(
                "XML is malformed. Parsing has stopped but partial data is available.",
                UserWarning,
                stacklevel=2,
            )

        self.ionic_steps = ionic_steps
        self.md_data = md_data
        self.vasp_version = self.generator["version"]

    @property
    def structures(self) -> list[Structure]:
        """List of Structures for each ionic step."""
        return [step["structure"] for step in self.ionic_steps]

    @property
    def epsilon_static(self) -> list[float]:
        """The static part of the dielectric constant.
        Present only when it's a DFPT run (LEPSILON=TRUE).
        """
        return self.ionic_steps[-1].get("epsilon", [])

    @property
    def epsilon_static_wolfe(self) -> list[float]:
        """The static part of the dielectric constant without any local
        field effects. Present only when it's a DFPT run (LEPSILON=TRUE).
        """
        return self.ionic_steps[-1].get("epsilon_rpa", [])

    @property
    def epsilon_ionic(self) -> list[float]:
        """The ionic part of the static dielectric constant.
        Present when it's a DFPT run (LEPSILON=TRUE) and IBRION=5, 6, 7 or 8.
        """
        return self.ionic_steps[-1].get("epsilon_ion", [])

    @property
    def dielectric(self) -> tuple[list, list, list]:
        """The real and imaginary part of the dielectric constant (e.g.,
        computed by RPA) in function of the energy (frequency).
        Optical properties (e.g. absorption coefficient) can be obtained through this.

        Returns:
            The data is given as a tuple of 3 for the energy, the real part
            tensor, and the imaginary part tensor:
            ([energies], [[real_partxx, real_partyy, real_partzz, real_partxy,
            real_partyz, real_partxz]], [[imag_partxx, imag_partyy, imag_partzz,
            imag_partxy, imag_partyz, imag_partxz]]).
        """
        return self.dielectric_data["density"]

    @property
    def optical_absorption_coeff(self) -> list[float] | None:
        """The optical absorption coefficient from the dielectric constants.
        Note that this method is only implemented for optical properties
        calculated with GGA and BSE.
        """
        diel_data = self.dielectric_data.get("freq_dependent") or self.dielectric_data["density"]
        if diel_data:
            real_avg = [sum(diel_data[1][i][:3]) / 3 for i in range(len(diel_data[0]))]
            imag_avg = [sum(diel_data[2][i][:3]) / 3 for i in range(len(diel_data[0]))]

            def optical_absorb_coeff(freq: float, real: float, imag: float) -> float:
                """Calculate optical absorption coefficient,
                the unit is cm^-1.
                """
                hc = 1.23984 * 1e-4  # plank constant times speed of light, in the unit of eV*cm
                return 2 * 3.14159 * np.sqrt(np.sqrt(real**2 + imag**2) - real) * np.sqrt(2) / hc * freq

            return list(
                itertools.starmap(
                    optical_absorb_coeff,
                    zip(
                        diel_data[0],
                        real_avg,
                        imag_avg,
                        strict=True,
                    ),
                )
            )
        return None

    @property
    def converged_electronic(self) -> bool:
        """Whether electronic step converged in the final ionic step."""
        final_elec_steps: list[dict[str, Any]] | Literal[0] = (
            self.ionic_steps[-1]["electronic_steps"] if self.incar.get("ALGO", "").lower() != "chi" else 0
        )
        # In a response function run there is no ionic steps, there is no SCF step
        if final_elec_steps == 0:
            raise ValueError("there is no ionic step in response function ALGO=CHI.")

        if self.incar.get("LEPSILON"):
            idx = 1
            to_check = {"e_wo_entrp", "e_fr_energy", "e_0_energy"}
            while set(final_elec_steps[idx]) == to_check:
                idx += 1
            return idx + 1 != self.parameters["NELM"]
        if self.incar.get("ALGO", "").upper() == "EXACT" and self.incar.get("NELM") == 1:
            return True
        return len(final_elec_steps) < self.parameters["NELM"]

    @property
    def converged_ionic(self) -> bool:
        """Whether ionic step convergence has been reached, i.e. VASP
        exited before reaching the max ionic steps for a relaxation run.
        In case IBRION=0 (MD) or EDIFFG=0, returns True if the max ionic
        steps are reached.
        """
        nsw = self.parameters.get("NSW", 0)
        ibrion = self.parameters.get("IBRION", -1 if nsw in (-1, 0) else 0)
        if ibrion == 0:
            return nsw <= 1 or self.md_n_steps == nsw

        # Context re EDIFFG: the use case for EDIFFG=0 is to ensure a relaxation runs for
        # NSW steps (the non-AIMD way to generate a relaxation trajectory with DFT). In
        # that case, user isn't worried about convergence w.r.t. forces or energy. The
        # next if statement prevents custodian from trying to correct the calc because
        # Vasprun.converged_ionic = False.
        ediffg = self.parameters.get("EDIFFG", 1)
        if ibrion in {1, 2} and ediffg == 0:
            return nsw <= 1 or nsw == len(self.ionic_steps)

        return nsw <= 1 or len(self.ionic_steps) < nsw

    @property
    def converged(self) -> bool:
        """Whether a relaxation run has both ionically and
        electronically converged.
        """
        return self.converged_electronic and self.converged_ionic

    @property
    @unitized("eV")
    def final_energy(self) -> float:
        """Final energy from the VASP run."""
        try:
            final_istep = self.ionic_steps[-1]
            total_energy = final_istep["e_0_energy"]

            # Fix a bug in vasprun.xml.
            # See https://www.vasp.at/forum/viewtopic.php?f=3&t=16942
            final_estep = final_istep["electronic_steps"][-1]
            electronic_energy_diff = final_estep["e_0_energy"] - final_estep["e_fr_energy"]
            total_energy_bugfix = np.round(electronic_energy_diff + final_istep["e_fr_energy"], 8)
            if np.abs(total_energy - total_energy_bugfix) > 1e-7:
                return total_energy_bugfix

            return total_energy

        except (IndexError, KeyError):
            warnings.warn(
                "Calculation does not have a total energy. "
                "Possibly a GW or similar kind of run. "
                "Infinity is returned."
            )
            return float("inf")

    @property
    def complete_dos(self) -> CompleteDos:
        """A CompleteDos object which incorporates the total DOS and all projected DOS."""
        final_struct = self.final_structure
        pdoss = {final_struct[i]: pdos for i, pdos in enumerate(self.pdos)}
        return CompleteDos(self.final_structure, self.tdos, pdoss)

    @property
    def complete_dos_normalized(self) -> CompleteDos:
        """A CompleteDos object which incorporates the total DOS and all projected DOS.
        Normalized by the volume of the unit cell with units of states/eV/unit cell
        volume.
        """
        final_struct = self.final_structure
        pdoss = {final_struct[i]: pdos for i, pdos in enumerate(self.pdos)}
        return CompleteDos(self.final_structure, self.tdos, pdoss, normalize=True)

    @property
    def hubbards(self) -> dict[str, float]:
        """Hubbard U values used for a GGA+U run, otherwise an empty dict."""
        symbols = [s.split()[1] for s in self.potcar_symbols]
        symbols = [re.split(r"_", s)[0] for s in symbols]
        if not self.incar.get("LDAU", False):
            return {}

        us = self.incar.get("LDAUU", self.parameters.get("LDAUU"))
        js = self.incar.get("LDAUJ", self.parameters.get("LDAUJ"))
        if len(js) != len(us):
            js = [0] * len(us)
        if len(us) == len(symbols):
            return {symbols[idx]: us[idx] - js[idx] for idx in range(len(symbols))}
        if sum(us) == 0 and sum(js) == 0:
            return {}

        raise VaspParseError("Length of U value parameters and atomic symbols are mismatched")

    @property
    def run_type(self) -> str:
        """The run type. Currently detects GGA, metaGGA, HF, HSE, B3LYP,
        and hybrid functionals based on relevant INCAR tags. LDA is assigned if
        PAW POTCARs are used and no other functional is detected.

        Hubbard U terms and vdW corrections are detected automatically as well.
        """
        # Care should be taken: if a GGA tag is not specified, VASP will default
        # to the functional specified by the POTCAR. It is not clear how
        # VASP handles the "--" value.
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
            0: "no-correction",
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
            run_type = "HF"
        elif self.parameters.get("HFSCREEN", 0.30) == 0.30:
            run_type = "HSE03"
        elif self.parameters.get("HFSCREEN", 0.20) == 0.20:
            run_type = "HSE06"
        elif self.parameters.get("AEXX", 0.20) == 0.20:
            run_type = "B3LYP"
        elif self.parameters.get("LHFCALC", True):
            run_type = "PBEO or other Hybrid Functional"
        elif self.incar.get("METAGGA") and self.incar.get("METAGGA") not in {
            "--",
            "None",
        }:
            incar_tag = self.incar.get("METAGGA", "").strip().upper()
            run_type = METAGGA_TYPES.get(incar_tag, incar_tag)
        elif self.parameters.get("GGA"):
            incar_tag = self.parameters.get("GGA", "").strip().upper()
            run_type = GGA_TYPES.get(incar_tag, incar_tag)
        elif self.potcar_symbols[0].split()[0] == "PAW":
            run_type = "LDA"
        else:
            run_type = "unknown"
            warnings.warn("Unknown run type!")

        if self.is_hubbard or self.parameters.get("LDAU", True):
            run_type += "+U"

        if self.parameters.get("LUSE_VDW", False):
            run_type += "+rVV10"
        elif self.incar.get("IVDW") in IVDW_TYPES:
            run_type += f"+vdW-{IVDW_TYPES[self.incar.get('IVDW', 0)]}"
        elif self.incar.get("IVDW"):
            run_type += "+vdW-unknown"

        return run_type

    @property
    def is_hubbard(self) -> bool:
        """Whether is a DFT+U run."""
        if len(self.hubbards) == 0:
            return False
        return sum(self.hubbards.values()) > 1e-8

    @property
    def is_spin(self) -> bool:
        """Whether is spin-polarized."""
        return self.parameters.get("ISPIN", 1) == 2

    @property
    def md_n_steps(self) -> int:
        """Number of steps for MD runs.

        Count all the actual MD steps if ML enabled.
        """
        return len(self.md_data) if self.md_data else self.nionic_steps

    def get_computed_entry(
        self,
        inc_structure: bool = True,
        parameters: list[str] | None = None,
        data: dict | None = None,
        entry_id: str | None = None,
    ) -> ComputedStructureEntry | ComputedEntry:
        """Get a ComputedEntry or ComputedStructureEntry from the Vasprun.

        Args:
            inc_structure (bool): Whether to return ComputedStructureEntries
                instead of ComputedEntries.
            parameters (list): Input parameters to include. It has to be one of
                the properties supported by the Vasprun object. If is None,
                a default set of parameters that are
                necessary for typical post-processing will be set.
            data (dict): Output data to include. Have to be the properties
                supported by the Vasprun object.
            entry_id (str): An entry id for the ComputedEntry.
                Defaults to "vasprun-{current datetime}"

        Returns:
            ComputedStructureEntry/ComputedEntry
        """
        if entry_id is None:
            entry_id = f"vasprun-{datetime.now(tz=timezone.utc)}"
        param_names = {
            "is_hubbard",
            "hubbards",
            "potcar_symbols",
            "potcar_spec",
            "run_type",
        }
        if parameters is not None and len(parameters) > 0:
            param_names.update(parameters)
        params = {param: getattr(self, param) for param in param_names}
        data = {} if data is None else {param: getattr(self, param) for param in data}

        if inc_structure:
            return ComputedStructureEntry(
                self.final_structure,
                self.final_energy,
                parameters=params,
                data=data,
                entry_id=entry_id,
            )
        return ComputedEntry(
            self.final_structure.composition,
            self.final_energy,
            parameters=params,
            data=data,
            entry_id=entry_id,
        )

    def get_band_structure(
        self,
        kpoints_filename: str | None = None,
        efermi: float | Literal["smart"] | None = None,
        line_mode: bool = False,
        force_hybrid_mode: bool = False,
        ignore_kpoints_opt: bool = False,
    ) -> BandStructureSymmLine | BandStructure:
        """Get the band structure.

        Args:
            kpoints_filename: Full path of the KPOINTS file from which
                the band structure is generated.
                If None is provided, the code will try to intelligently
                determine the appropriate KPOINTS file by substituting the
                filename of the vasprun.xml with KPOINTS (or KPOINTS_OPT).
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
            ignore_kpoints_opt: Normally, if KPOINTS_OPT data exists, it has
                the band structure data. Set this flag to ignore it. (Default: False)

        Returns:
            BandStructure (or more specifically a BandStructureSymmLine object if the run
            is detected to be a run along symmetry lines)

            Two types of runs along symmetry lines are accepted: non-sc with
            Line-Mode in the KPOINT file or hybrid, self-consistent with a
            uniform grid+a few kpoints along symmetry lines (explicit KPOINTS
            file) (it's not possible to run a non-sc band structure with hybrid
            functionals). The explicit KPOINTS file needs to have data on the
            kpoint label as commentary.

            If VASP was run with KPOINTS_OPT, it reads the data from that
            file unless told otherwise. This overrides hybrid mode.
        """
        use_kpoints_opt = not ignore_kpoints_opt and (getattr(self, "kpoints_opt_props", None) is not None)
        if not kpoints_filename:
            kpts_path = os.path.join(
                os.path.dirname(self.filename),
                "KPOINTS_OPT" if use_kpoints_opt else "KPOINTS",
            )
            kpoints_filename = zpath(kpts_path)
        if kpoints_filename and not os.path.isfile(kpoints_filename) and line_mode:
            name = "KPOINTS_OPT" if use_kpoints_opt else "KPOINTS"
            raise VaspParseError(f"{name} not found but needed to obtain band structure along symmetry lines.")

        # Note that we're using the Fermi energy of the self-consistent grid
        # run even if we're reading bands from KPOINTS_OPT.
        if efermi == "smart":
            e_fermi: float | None = self.calculate_efermi()
        elif efermi is None:
            e_fermi = self.efermi
        else:
            e_fermi = efermi

        if e_fermi is None:
            raise ValueError("e_fermi is None.")

        kpoint_file: Kpoints | None = None
        if kpoints_filename and os.path.isfile(kpoints_filename):
            kpoint_file = Kpoints.from_file(kpoints_filename)
        lattice_new = Lattice(self.final_structure.lattice.reciprocal_lattice.matrix)

        if use_kpoints_opt:
            if self.kpoints_opt_props is None or self.kpoints_opt_props.actual_kpoints is None:
                raise RuntimeError("KPOINTS_opt or actual_kpoints is None.")
            kpoints = [np.array(kpt) for kpt in self.kpoints_opt_props.actual_kpoints]
        else:
            if self.actual_kpoints is None:
                raise RuntimeError("actual_kpoints is None.")
            kpoints = [np.array(kpt) for kpt in self.actual_kpoints]

        p_eig_vals: defaultdict[Spin, list] = defaultdict(list)
        eigenvals: defaultdict[Spin, list] = defaultdict(list)

        n_kpts = len(kpoints)

        if use_kpoints_opt:
            if self.kpoints_opt_props is None:
                raise RuntimeError("KPOINTS_opt is None.")
            eig_vals = self.kpoints_opt_props.eigenvalues
            projected_eig_vals = self.kpoints_opt_props.projected_eigenvalues
        else:
            eig_vals = self.eigenvalues
            projected_eig_vals = self.projected_eigenvalues

        if eig_vals is None:
            raise ValueError("Eigenvalues are None.")

        for spin, val in eig_vals.items():
            val = np.swapaxes(val, 0, 1)
            eigenvals[spin] = val[:, :, 0]

            if projected_eig_vals:
                proj_eig_vals = projected_eig_vals[spin]
                # Original axes for self.projected_eigenvalues are kpoints, band, ion, orb.
                # For BS input, we need band, kpoints, orb, ion.
                proj_eig_vals = np.swapaxes(proj_eig_vals, 0, 1)  # Swap kpoint and band axes
                proj_eig_vals = np.swapaxes(proj_eig_vals, 2, 3)  # Swap ion and orb axes

                p_eig_vals[spin] = proj_eig_vals

        # Check if we have an hybrid band structure computation
        # for this we look at the presence of the LHFCALC tag
        # (but hybrid mode is redundant if using kpoints_opt)
        hybrid_band = False
        if self.parameters.get("LHFCALC", False) or 0.0 in self.actual_kpoints_weights:
            hybrid_band = True

        if kpoint_file is not None and kpoint_file.style == Kpoints.supported_modes.Line_mode:
            line_mode = True

        if line_mode:
            labels_dict = {}
            if kpoint_file is None:
                raise RuntimeError("Kpoint file cannot be None for line mode.")

            if (hybrid_band or force_hybrid_mode) and not use_kpoints_opt:
                start_bs_index = 0
                for i in range(len(self.actual_kpoints)):
                    if self.actual_kpoints_weights[i] == 0.0:
                        start_bs_index = i
                        break
                for i in range(start_bs_index, len(kpoint_file.kpts)):
                    if kpoint_file.labels is not None and kpoint_file.labels[i] is not None:
                        labels_dict[kpoint_file.labels[i]] = kpoint_file.kpts[i]
                # Remake the data only considering line band structure k-points
                # (weight = 0.0 kpoints)
                n_bands = len(eigenvals[Spin.up])
                kpoints = kpoints[start_bs_index:n_kpts]
                up_eigen = [eigenvals[Spin.up][i][start_bs_index:n_kpts] for i in range(n_bands)]
                if self.projected_eigenvalues:
                    p_eig_vals[Spin.up] = [p_eig_vals[Spin.up][i][start_bs_index:n_kpts] for i in range(n_bands)]
                if self.is_spin:
                    down_eigen = [eigenvals[Spin.down][i][start_bs_index:n_kpts] for i in range(n_bands)]
                    eigenvals[Spin.up] = up_eigen
                    eigenvals[Spin.down] = down_eigen
                    if self.projected_eigenvalues:
                        p_eig_vals[Spin.down] = [
                            p_eig_vals[Spin.down][i][start_bs_index:n_kpts] for i in range(n_bands)
                        ]
                else:
                    eigenvals[Spin.up] = up_eigen
            else:
                if kpoint_file.labels is not None:
                    if "" in kpoint_file.labels:
                        raise ValueError(
                            "A band structure along symmetry lines requires a label "
                            "for each kpoint. Check your KPOINTS file"
                        )
                    labels_dict = dict(zip(kpoint_file.labels, kpoint_file.kpts, strict=True))
                labels_dict.pop(None, None)  # type: ignore[call-overload]

            return BandStructureSymmLine(
                kpoints,
                eigenvals,
                lattice_new,
                e_fermi,
                labels_dict,  # type: ignore[arg-type]
                structure=self.final_structure,
                projections=p_eig_vals,
            )

        return BandStructure(
            kpoints,
            eigenvals,
            lattice_new,
            e_fermi,
            structure=self.final_structure,
            projections=p_eig_vals,
        )

    @property
    def eigenvalue_band_properties(
        self,
    ) -> (
        tuple[float, float, float, bool]
        | tuple[
            tuple[float, float],
            tuple[float, float],
            tuple[float, float],
            tuple[bool, bool],
        ]
    ):
        """Band properties from the eigenvalues as a tuple,
        (band gap, cbm, vbm, is_band_gap_direct).
        In the case of separate_spins=True, the band gap, cbm, vbm, and is_band_gap_direct are each
        lists of length 2, with index 0 representing the spin-up channel and index 1 representing the
        spin-down channel.
        """
        vbm = -float("inf")
        vbm_kpoint = None
        cbm = float("inf")
        cbm_kpoint = None

        vbm_spins = []
        vbm_spins_kpoints = []
        cbm_spins = []
        cbm_spins_kpoints = []

        if self.eigenvalues is None:
            raise ValueError("eigenvalues is None.")

        if self.separate_spins and len(self.eigenvalues) != 2:
            raise ValueError("The separate_spins flag can only be True if ISPIN = 2")

        for eigenvalue in self.eigenvalues.values():
            if self.separate_spins:
                vbm = -float("inf")
                cbm = float("inf")
            for kpoint, val in enumerate(eigenvalue):
                for eigenval, occu in val:
                    if occu > self.occu_tol and eigenval > vbm:
                        vbm = eigenval
                        vbm_kpoint = kpoint
                    elif occu <= self.occu_tol and eigenval < cbm:
                        cbm = eigenval
                        cbm_kpoint = kpoint
            if self.separate_spins:
                vbm_spins.append(vbm)
                vbm_spins_kpoints.append(vbm_kpoint)
                cbm_spins.append(cbm)
                cbm_spins_kpoints.append(cbm_kpoint)
        if self.separate_spins:
            return (
                (
                    max(cbm_spins[0] - vbm_spins[0], 0),
                    max(cbm_spins[1] - vbm_spins[1], 0),
                ),
                (cbm_spins[0], cbm_spins[1]),
                (vbm_spins[0], vbm_spins[1]),
                (
                    vbm_spins_kpoints[0] == cbm_spins_kpoints[0],
                    vbm_spins_kpoints[1] == cbm_spins_kpoints[1],
                ),
            )
        return max(cbm - vbm, 0), cbm, vbm, vbm_kpoint == cbm_kpoint

    def calculate_efermi(self, tol: float = 0.001) -> float:
        """Calculate the Fermi level using a robust algorithm.

        Sometimes VASP can set the Fermi level inside a band
        due to issues in the way band occupancies are handled.
        This method tries to detect and correct this.

        More details: https://www.vasp.at/forum/viewtopic.php?f=4&t=17981.
        """
        if self.eigenvalues is None:
            raise ValueError("eigenvalues cannot be None.")
        if self.efermi is None:
            raise ValueError("efermi cannot be None.")

        # Drop weights and set shape n_bands, n_kpoints
        all_eigs = np.concatenate([eigs[:, :, 0].transpose(1, 0) for eigs in self.eigenvalues.values()])

        def crosses_band(fermi: float) -> bool:
            eigs_below = np.any(all_eigs < fermi, axis=1)
            eigs_above = np.any(all_eigs > fermi, axis=1)
            return np.any(eigs_above & eigs_below)

        def get_vbm_cbm(fermi):
            return np.max(all_eigs[all_eigs < fermi]), np.min(all_eigs[all_eigs > fermi])

        # Fermi level doesn't cross a band: safe to use VASP Fermi level
        if not crosses_band(self.efermi):
            return self.efermi

        # If the Fermi level crosses a band, check if we are very close to band gap.
        # If so, then likely this is a VASP tetrahedron bug
        if not crosses_band(self.efermi + tol):
            # efermi resides in the valence band
            # Set Fermi level as average of CBM and VBM
            vbm, cbm = get_vbm_cbm(self.efermi + tol)
            return (cbm + vbm) / 2

        if not crosses_band(self.efermi - tol):
            # efermi resides in the conduction band
            # Set Fermi level as average of CBM and VBM
            vbm, cbm = get_vbm_cbm(self.efermi - tol)
            return (cbm + vbm) / 2

        # If is a metal
        return self.efermi

    def get_potcars(self, path: PathLike | bool) -> Potcar | None:
        """Get the POTCAR from the specified path.

        Args:
            path (PathLike | bool): If a str or Path, the path to search for POTCARs.
                If a bool, whether to take the search path from the specified vasprun.xml

        Returns:
            Potcar | None: The POTCAR from the specified path or None if not found/no path specified.
        """
        if not path:
            return None

        if isinstance(path, str | Path) and "POTCAR" in str(path):
            potcar_paths = [str(path)]
        else:
            # the abspath is needed here in cases where no leading directory is specified,
            # e.g. Vasprun("vasprun.xml"). see gh-3586:
            search_path = os.path.dirname(os.path.abspath(self.filename)) if path is True else str(path)

            potcar_paths = [
                f"{search_path}/{fn}" for fn in os.listdir(search_path) if fn.startswith("POTCAR") and ".spec" not in fn
            ]

        for potcar_path in potcar_paths:
            try:
                potcar = Potcar.from_file(potcar_path)
                if {d.header for d in potcar} == set(self.potcar_symbols):
                    return potcar
            except Exception:
                continue

        warnings.warn("No POTCAR file with matching TITEL fields was found in\n" + "\n  ".join(potcar_paths))

        return None

    def get_trajectory(self) -> Trajectory:
        """
        Get a Trajectory, an alternative representation of self.structures
        as a single object. Forces are added as site properties.

        Returns:
            Trajectory
        """
        structs: list[Structure] = []
        steps = self.md_data or self.ionic_steps
        for step in steps:
            struct = step["structure"].copy()
            struct.add_site_property("forces", step["forces"])
            structs.append(struct)
        return Trajectory.from_structures(structs, constant_lattice=False)

    def update_potcar_spec(self, path: PathLike | bool) -> None:
        """Update the specs based on the POTCARs found.

        Args:
            path (PathLike | bool): Path to search for POTCARs.
        """
        if potcar := self.get_potcars(path):
            self.potcar_spec = [
                {
                    "titel": sym,
                    "hash": ps.md5_header_hash,
                    "summary_stats": ps._summary_stats,
                }
                for sym in self.potcar_symbols
                for ps in potcar
                if ps.symbol == sym.split()[1]
            ]

    def update_charge_from_potcar(self, path: PathLike | bool) -> None:
        """Update the charge of a structure based on the POTCARs found.

        Args:
            path (PathLike | bool): Path to search for POTCARs.
        """
        potcar = self.get_potcars(path)

        if potcar and self.incar.get("ALGO", "").upper() not in {
            "GW0",
            "G0W0",
            "GW",
            "BSE",
            # VASP renamed the GW tags in v6.
            "QPGW",
            "QPGW0",
            "EVGW",
            "EVGW0",
            "GWR",
            "GW0R",
        }:
            nelect = self.parameters["NELECT"]
            if len(potcar) == len(self.initial_structure.composition.element_composition):
                potcar_nelect = sum(
                    self.initial_structure.composition.element_composition[ps.element] * ps.ZVAL for ps in potcar
                )
            else:
                nums = [len(list(g)) for _, g in itertools.groupby(self.atomic_symbols)]
                potcar_nelect = sum(ps.ZVAL * num for ps, num in zip(potcar, nums, strict=False))
            charge = potcar_nelect - nelect

            for struct in self.structures:
                struct._charge = charge
            if hasattr(self, "initial_structure"):
                self.initial_structure._charge = charge
            if hasattr(self, "final_structure"):
                self.final_structure._charge = charge

    def as_dict(self) -> dict:
        """JSON-serializable dict representation."""
        comp = self.final_structure.composition
        unique_symbols = sorted(set(self.atomic_symbols))
        dct: dict[str, Any] = {
            "vasp_version": self.vasp_version,
            "has_vasp_completed": self.converged,
            "nsites": len(self.final_structure),
            "unit_cell_formula": comp.as_dict(),
            "reduced_cell_formula": Composition(comp.reduced_formula).as_dict(),
            "pretty_formula": comp.reduced_formula,
            "is_hubbard": self.is_hubbard,
            "hubbards": self.hubbards,
            "elements": unique_symbols,
            "nelements": len(unique_symbols),
            "run_type": self.run_type,
        }

        vin: dict[str, Any] = {
            "incar": dict(self.incar.items()),
            "crystal": self.initial_structure.as_dict(),
            "kpoints": self.kpoints.as_dict(),
        }
        actual_kpts = [
            {
                "abc": list(self.actual_kpoints[idx]),
                "weight": self.actual_kpoints_weights[idx],
            }
            for idx in range(len(self.actual_kpoints))
        ]
        vin["kpoints"]["actual_points"] = actual_kpts
        vin["nkpoints"] = len(actual_kpts)
        if kpt_opt_props := self.kpoints_opt_props:
            if (
                kpt_opt_props.kpoints is None
                or kpt_opt_props.actual_kpoints is None
                or kpt_opt_props.actual_kpoints_weights is None
            ):
                raise ValueError("kpt_opt_props.kpoints/actual_kpoints/actual_kpoints_weights cannot be None.")
            vin["kpoints_opt"] = kpt_opt_props.kpoints.as_dict()
            actual_kpts = [
                {
                    "abc": list(kpt_opt_props.actual_kpoints[idx]),
                    "weight": kpt_opt_props.actual_kpoints_weights[idx],
                }
                for idx in range(len(kpt_opt_props.actual_kpoints))
            ]
            vin["kpoints_opt"]["actual_kpoints"] = actual_kpts
            vin["nkpoints_opt"] = len(actual_kpts)
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
            gap, cbm, vbm, is_direct = self.eigenvalue_band_properties
            vout |= {"bandgap": gap, "cbm": cbm, "vbm": vbm, "is_gap_direct": is_direct}

            if self.projected_eigenvalues:
                vout["projected_eigenvalues"] = {
                    str(spin): v.tolist() for spin, v in self.projected_eigenvalues.items()
                }

            if self.projected_magnetisation is not None:
                vout["projected_magnetisation"] = self.projected_magnetisation.tolist()

        if kpt_opt_props and kpt_opt_props.eigenvalues:
            eigen = {str(spin): v.tolist() for spin, v in kpt_opt_props.eigenvalues.items()}
            vout["eigenvalues_kpoints_opt"] = eigen
            # TODO: implement kpoints_opt eigenvalue_band_proprties.
            # gap, cbm, vbm, is_direct = self.eigenvalue_band_properties
            # vout |= {"bandgap": gap, "cbm": cbm, "vbm": vbm, "is_gap_direct": is_direct}

            if kpt_opt_props.projected_eigenvalues:
                vout["projected_eigenvalues_kpoints_opt"] = {
                    str(spin): v.tolist() for spin, v in kpt_opt_props.projected_eigenvalues.items()
                }

            if kpt_opt_props.projected_magnetisation is not None:
                vout["projected_magnetisation_kpoints_opt"] = kpt_opt_props.projected_magnetisation.tolist()

        vout["epsilon_static"] = self.epsilon_static
        vout["epsilon_static_wolfe"] = self.epsilon_static_wolfe
        vout["epsilon_ionic"] = self.epsilon_ionic
        dct["output"] = vout
        return jsanitize(dct, strict=True)

    def _parse_params(self, elem: XML_Element) -> Incar[str, Any]:
        """Parse INCAR parameters and more."""
        params: dict[str, Any] = {}
        for c in elem:
            # VASP 6.4.3 can add trailing whitespace
            # for example, <i type="string" name="GGA    ">PE</i>
            name = c.attrib.get("name", "").strip()
            if c.tag not in {"i", "v"}:
                p = self._parse_params(c)
                if name == "response functions":
                    # Delete duplicate fields from "response functions",
                    # which overrides the values in the root params.
                    p = {k: v for k, v in p.items() if k not in params}
                params |= p

            else:
                ptype = c.attrib.get("type", "")
                val = c.text.strip() if c.text else ""
                try:
                    if c.tag == "i":
                        params[name] = _parse_parameters(ptype, val)
                    else:
                        params[name] = _parse_v_parameters(ptype, val, self.filename, name)

                except Exception:
                    if name == "RANDOM_SEED":
                        # Handle the case where RANDOM_SEED > 99999, which results in *****
                        params[name] = None
                    else:
                        raise
        elem.clear()
        return Incar(params)

    @staticmethod
    def _parse_atominfo(elem: XML_Element) -> tuple[list[str], list[str]]:
        """Parse atom symbols and POTCAR symbols."""

        def parse_atomic_symbol(symbol: str) -> str:
            """Parse and ensure atomic symbols are valid elements."""
            try:
                return str(Element(symbol))

            # vasprun.xml uses "X" instead of "Xe" for Xenon
            except ValueError:
                if symbol == "X":
                    return "Xe"
                if symbol == "r":
                    return "Zr"
                raise

        atomic_symbols = []
        potcar_symbols = []

        for a in elem.findall("array"):
            if a.attrib["name"] == "atoms":
                atomic_symbols = [rc.find("c").text.strip() for rc in a.find("set")]  # type: ignore[union-attr]
            elif a.attrib["name"] == "atomtypes":
                potcar_symbols = [rc.findall("c")[4].text.strip() for rc in a.find("set")]  # type: ignore[union-attr]

        elem.clear()
        return [parse_atomic_symbol(sym) for sym in atomic_symbols], potcar_symbols

    @staticmethod
    def _parse_kpoints(
        elem: XML_Element,
    ) -> tuple[Kpoints, list[Tuple3Floats], list[float]]:
        """Parse Kpoints."""
        e = elem if elem.find("generation") is None else elem.find("generation")
        kpoint = Kpoints("Kpoints from vasprun.xml")
        kpoint.style = Kpoints.supported_modes.from_str(e.attrib.get("param", "Reciprocal"))  # type: ignore[union-attr]

        for v in e.findall("v"):  # type: ignore[union-attr]
            name = v.attrib.get("name")  # type: ignore[union-attr]
            tokens = v.text.split()  # type: ignore[union-attr]

            if name == "divisions":
                kpoint.kpts = [
                    cast(Kpoint, tuple(int(i) for i in tokens)),
                ]
            elif name == "usershift":
                kpoint.kpts_shift = cast(Vector3D, tuple(float(i) for i in tokens))
            elif name in {"genvec1", "genvec2", "genvec3", "shift"}:
                setattr(kpoint, name, [float(i) for i in tokens])

        actual_kpoints = []
        weights = []
        for va in elem.findall("varray"):
            name = va.attrib["name"]
            if name == "kpointlist":
                actual_kpoints = cast(list[Tuple3Floats], list(map(tuple, _parse_vasp_array(va))))
            elif name == "weights":
                weights = [i[0] for i in _parse_vasp_array(va)]
        elem.clear()

        if kpoint.style == Kpoints.supported_modes.Reciprocal:
            kpoint = Kpoints(
                comment="Kpoints from vasprun.xml",
                style=Kpoints.supported_modes.Reciprocal,
                num_kpts=len(kpoint.kpts),
                kpts=actual_kpoints,
                kpts_weights=weights,
            )
        return kpoint, actual_kpoints, weights

    def _parse_structure(self, elem: XML_Element) -> Structure:
        """Parse Structure with lattice, positions and selective dynamics info."""
        lattice = _parse_vasp_array(elem.find("crystal").find("varray"))  # type: ignore[union-attr]
        pos = _parse_vasp_array(elem.find("varray"))
        struct = Structure(lattice, self.atomic_symbols, pos)
        selective_dyn = elem.find("varray/[@name='selective']")

        if selective_dyn is not None:
            struct.add_site_property("selective_dynamics", _parse_vasp_array(selective_dyn))
        return struct

    @staticmethod
    def _parse_diel(elem: XML_Element) -> tuple[list, list, list]:
        """Parse dielectric properties."""
        imag = [
            [_vasprun_float(line) for line in r.text.split()]  # type: ignore[union-attr]
            for r in elem.find("imag").find("array").find("set").findall("r")  # type: ignore[union-attr]
        ]
        real = [
            [_vasprun_float(line) for line in r.text.split()]  # type: ignore[union-attr]
            for r in elem.find("real").find("array").find("set").findall("r")  # type: ignore[union-attr]
        ]
        elem.clear()
        return [e[0] for e in imag], [e[1:] for e in real], [e[1:] for e in imag]

    @staticmethod
    def _parse_optical_transition(elem: XML_Element) -> tuple[NDArray, NDArray]:
        """Parse optical transitions."""
        for va in elem.findall("varray"):
            if va.attrib.get("name") == "opticaltransitions":
                # optical transitions array contains oscillator strength and probability of transition
                oscillator_strength = np.array(_parse_vasp_array(va))[:]
                probability_transition = np.array(_parse_vasp_array(va))[:, 1]

                return oscillator_strength, probability_transition

        raise RuntimeError("Failed to parse optical transitions.")

    def _parse_chemical_shielding(self, elem: XML_Element) -> list[dict[str, Any]]:
        """Parse NMR chemical shielding."""
        istep: dict[str, Any] = {}
        # not all calculations have a structure
        _struct = elem.find("structure")
        struct = None if _struct is None else self._parse_structure(_struct)

        for va in elem.findall("varray"):
            istep[va.attrib["name"]] = _parse_vasp_array(va)
        istep["structure"] = struct
        istep["electronic_steps"] = []
        calculation = [
            istep,
        ]
        for scstep in elem.findall("scstep"):
            try:
                e_steps_dict = {i.attrib["name"]: _vasprun_float(i.text) for i in scstep.find("energy").findall("i")}  # type: ignore[union-attr, arg-type]
                cur_ene = e_steps_dict["e_fr_energy"]
                min_steps = 1 if len(calculation) >= 1 else self.parameters.get("NELMIN", 5)
                if len(calculation[-1]["electronic_steps"]) <= min_steps:
                    calculation[-1]["electronic_steps"].append(e_steps_dict)
                else:
                    last_ene = calculation[-1]["electronic_steps"][-1]["e_fr_energy"]  # type: ignore[call-overload]
                    if abs(cur_ene - last_ene) < 1.0:
                        calculation[-1]["electronic_steps"].append(e_steps_dict)
                    else:
                        calculation.append({"electronic_steps": [e_steps_dict]})
            except AttributeError:  # not all calculations have an energy
                pass
        calculation[-1].update(calculation[-1]["electronic_steps"][-1])
        return calculation

    def _parse_ionic_step(self, elem: XML_Element) -> dict[str, float]:
        """Parse an ionic step."""
        try:
            ion_step: dict[str, Any] = {
                i.attrib["name"]: _vasprun_float(i.text)  # type: ignore[arg-type]
                for i in elem.find("energy").findall("i")  # type: ignore[union-attr]
            }
        # Not all calculations have an energy
        except AttributeError:
            ion_step = {}

        elec_steps = []
        for scstep in elem.findall("scstep"):
            try:
                e_step_dict = {i.attrib["name"]: _vasprun_float(i.text) for i in scstep.find("energy").findall("i")}  # type: ignore[union-attr, arg-type]
                elec_steps.append(e_step_dict)
            # Not all calculations have an energy
            except AttributeError:
                pass

        try:
            struct = self._parse_structure(elem.find("structure"))  # type: ignore[arg-type]
        except AttributeError:  # not all calculations have a structure
            struct = None

        for va in elem.findall("varray"):
            ion_step[va.attrib["name"]] = _parse_vasp_array(va)
        ion_step["electronic_steps"] = elec_steps
        ion_step["structure"] = struct
        elem.clear()
        return ion_step

    @staticmethod
    def _parse_dos(elem: XML_Element) -> tuple[Dos, Dos, list[dict]]:
        """Parse density of states (DOS)."""
        efermi = float(elem.find("i").text)  # type: ignore[union-attr, arg-type]
        energies: NDArray | None = None
        tdensities = {}
        idensities = {}

        for s in elem.find("total").find("array").find("set").findall("set"):  # type: ignore[union-attr]
            data = np.array(_parse_vasp_array(s))
            energies = data[:, 0]
            spin = Spin.up if s.attrib["comment"] == "spin 1" else Spin.down
            tdensities[spin] = data[:, 1]
            idensities[spin] = data[:, 2]

        pdoss = []
        partial = elem.find("partial")
        if partial is not None:
            orbs = [ss.text for ss in partial.find("array").findall("field")]  # type: ignore[union-attr]
            orbs.pop(0)
            lm = any("x" in s for s in orbs if s is not None)
            for s in partial.find("array").find("set").findall("set"):  # type: ignore[union-attr]
                pdos: dict[Orbital | OrbitalType, dict[Spin, np.ndarray]] = defaultdict(dict)

                for ss in s.findall("set"):
                    spin = Spin.up if ss.attrib["comment"] == "spin 1" else Spin.down
                    data = np.array(_parse_vasp_array(ss))
                    _n_row, n_col = data.shape
                    for col_idx in range(1, n_col):
                        orb = Orbital(col_idx - 1) if lm else OrbitalType(col_idx - 1)
                        pdos[orb][spin] = data[:, col_idx]  # type: ignore[index]
                pdoss.append(pdos)
        elem.clear()

        if energies is None:
            raise ValueError("energies is None")
        return (
            Dos(efermi, energies, tdensities),
            Dos(efermi, energies, idensities),
            pdoss,
        )

    @staticmethod
    def _parse_eigen(elem: XML_Element) -> dict[Spin, NDArray]:
        """Parse eigenvalues."""
        eigenvalues: dict[Spin, np.ndarray] = defaultdict(list)
        for s in elem.find("array").find("set").findall("set"):  # type: ignore[union-attr]
            spin = Spin.up if s.attrib["comment"] == "spin 1" else Spin.down
            for ss in s.findall("set"):
                eigenvalues[spin].append(_parse_vasp_array(ss))
        eigenvalues = {spin: np.array(v) for spin, v in eigenvalues.items()}
        elem.clear()
        return eigenvalues

    @staticmethod
    def _parse_projected_eigen(
        elem: XML_Element,
    ) -> tuple[dict[Spin, NDArray], NDArray | None]:
        """Parse projected eigenvalues."""
        root = elem.find("array").find("set")  # type: ignore[union-attr]
        _proj_eigen: dict[int, np.ndarray] = defaultdict(list)
        for s in root.findall("set"):  # type: ignore[union-attr]
            spin: int = int(re.match(r"spin(\d+)", s.attrib["comment"])[1])  # type: ignore[index]

            # Force spin to be +1 or -1
            for ss in s.findall("set"):
                dk = []
                for sss in ss.findall("set"):
                    db = _parse_vasp_array(sss)
                    dk.append(db)
                _proj_eigen[spin].append(dk)
        _proj_eigen = {spin: np.array(v) for spin, v in _proj_eigen.items()}

        if len(_proj_eigen) > 2:
            # non-collinear magentism (also spin-orbit coupling) enabled, last three
            # "spin channels" are the projected magnetization of the orbitals in the
            # x, y, and z Cartesian coordinates
            proj_mag = np.stack([_proj_eigen.pop(i) for i in range(2, 5)], axis=-1)  # type: ignore[call-overload]
            proj_eigen: dict[Spin, np.ndarray] = {Spin.up: _proj_eigen[1]}
        else:
            proj_eigen = {Spin.up if k == 1 else Spin.down: v for k, v in _proj_eigen.items()}
            proj_mag = None

        elem.clear()
        return proj_eigen, proj_mag

    @staticmethod
    def _parse_dynmat(elem: XML_Element) -> tuple[list, list, list]:
        """Parse dynamical matrix."""
        hessian = []
        eigenvalues = []
        eigenvectors = []

        for v in elem.findall("v"):
            if v.attrib["name"] == "eigenvalues":
                eigenvalues = [float(i) for i in v.text.split()]  # type: ignore[union-attr]

        for va in elem.findall("varray"):
            if va.attrib["name"] == "hessian":
                for v in va.findall("v"):
                    hessian.append([float(i) for i in v.text.split()])  # type: ignore[union-attr]

            elif va.attrib["name"] == "eigenvectors":
                for v in va.findall("v"):
                    eigenvectors.append([float(i) for i in v.text.split()])  # type: ignore[union-attr]

        return hessian, eigenvalues, eigenvectors


class BSVasprun(Vasprun):
    """
    A highly optimized version of Vasprun that parses only eigenvalues for
    bandstructures. All other properties like structures, parameters,
    etc. are ignored.
    """

    def __init__(
        self,
        filename: PathLike,
        parse_projected_eigen: bool | str = False,
        parse_potcar_file: bool | str = False,
        occu_tol: float = 1e-8,
        separate_spins: bool = False,
    ) -> None:
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

        with zopen(filename, mode="rt") as file:
            self.efermi = None
            parsed_header = False
            in_kpoints_opt = False
            self.eigenvalues = self.projected_eigenvalues = None
            self.kpoints_opt_props = None
            for event, elem in ET.iterparse(file, events=["start", "end"]):
                tag = elem.tag
                if event == "start":
                    # The start event tells us when we have entered blocks
                    if tag in {"eigenvalues_kpoints_opt", "projected_kpoints_opt"} or (
                        tag == "dos" and elem.attrib.get("comment") == "kpoints_opt"
                    ):
                        in_kpoints_opt = True
                elif not parsed_header:
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
                        self.potcar_spec = [
                            {"titel": p, "hash": None, "summary_stats": {}} for p in self.potcar_symbols
                        ]
                        parsed_header = True
                elif tag == "i" and elem.attrib.get("name") == "efermi":
                    if in_kpoints_opt:
                        if self.kpoints_opt_props is None:
                            self.kpoints_opt_props = KpointOptProps()
                        self.kpoints_opt_props.efermi = float(elem.text)
                        in_kpoints_opt = False
                    else:
                        self.efermi = float(elem.text)
                elif tag == "eigenvalues" and not in_kpoints_opt:
                    self.eigenvalues = self._parse_eigen(elem)
                elif parse_projected_eigen and tag == "projected" and not in_kpoints_opt:
                    self.projected_eigenvalues, self.projected_magnetisation = self._parse_projected_eigen(elem)
                elif tag in ("eigenvalues_kpoints_opt", "projected_kpoints_opt"):
                    if self.kpoints_opt_props is None:
                        self.kpoints_opt_props = KpointOptProps()
                    in_kpoints_opt = False
                    # projected_kpoints_opt includes occupation information whereas
                    # eigenvalues_kpoints_opt doesn't.
                    self.kpoints_opt_props.eigenvalues = self._parse_eigen(elem.find("eigenvalues"))
                    if tag == "eigenvalues_kpoints_opt":
                        (
                            self.kpoints_opt_props.kpoints,
                            self.kpoints_opt_props.actual_kpoints,
                            self.kpoints_opt_props.actual_kpoints_weights,
                        ) = self._parse_kpoints(elem.find("kpoints"))
                    elif parse_projected_eigen:  # and tag == "projected_kpoints_opt": (implied)
                        (
                            self.kpoints_opt_props.projected_eigenvalues,
                            self.kpoints_opt_props.projected_magnetisation,
                        ) = self._parse_projected_eigen(elem)
                elif tag == "structure" and elem.attrib.get("name") == "finalpos":
                    self.final_structure = self._parse_structure(elem)
        self.vasp_version = self.generator["version"]
        if parse_potcar_file:
            self.update_potcar_spec(parse_potcar_file)

    def as_dict(self) -> dict:
        """JSON-serializable dict representation."""
        comp = self.final_structure.composition
        unique_symbols = sorted(set(self.atomic_symbols))
        dct = {
            "vasp_version": self.vasp_version,
            "has_vasp_completed": True,
            "nsites": len(self.final_structure),
            "unit_cell_formula": comp.as_dict(),
            "reduced_cell_formula": Composition(comp.reduced_formula).as_dict(),
            "pretty_formula": comp.reduced_formula,
            "is_hubbard": self.is_hubbard,
            "hubbards": self.hubbards,
            "elements": unique_symbols,
            "nelements": len(unique_symbols),
            "run_type": self.run_type,
        }

        vin: dict[str, Any] = {
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
        if kpt_opt_props := getattr(self, "kpoints_opt_props", None):
            vin["kpoints_opt"] = kpt_opt_props.kpoints.as_dict()
            actual_kpts = [
                {
                    "abc": list(kpt_opt_props.actual_kpoints[idx]),
                    "weight": kpt_opt_props.actual_kpoints_weights[idx],
                }
                for idx in range(len(kpt_opt_props.actual_kpoints))
            ]
            vin["kpoints_opt"]["actual_kpoints"] = actual_kpts
            vin["nkpoints_opt"] = len(actual_kpts)
        vin["potcar"] = [s.split(" ")[1] for s in self.potcar_symbols]
        vin["potcar_spec"] = self.potcar_spec
        vin["potcar_type"] = [s.split(" ")[0] for s in self.potcar_symbols]
        vin["parameters"] = dict(self.parameters)
        vin["lattice_rec"] = self.final_structure.lattice.reciprocal_lattice.as_dict()
        dct["input"] = vin

        vout: dict[str, Any] = {
            "crystal": self.final_structure.as_dict(),
            "efermi": self.efermi,
        }

        if self.eigenvalues:
            eigen: dict[Any, Any] = defaultdict(dict)
            for spin, values in self.eigenvalues.items():
                for idx, val in enumerate(values):
                    eigen[idx][str(spin)] = val
            vout["eigenvalues"] = eigen
            gap, cbm, vbm, is_direct = self.eigenvalue_band_properties
            vout |= {"bandgap": gap, "cbm": cbm, "vbm": vbm, "is_gap_direct": is_direct}

            if self.projected_eigenvalues:
                peigen: list[dict] = [{} for _ in eigen]
                for spin, val in self.projected_eigenvalues.items():
                    for kpoint_index, vv in enumerate(val):
                        if str(spin) not in peigen[kpoint_index]:
                            peigen[kpoint_index][str(spin)] = vv
                vout["projected_eigenvalues"] = peigen

        if kpt_opt_props and kpt_opt_props.eigenvalues:
            eigen = {str(spin): v.tolist() for spin, v in kpt_opt_props.eigenvalues.items()}
            vout["eigenvalues_kpoints_opt"] = eigen
            # TODO: implement kpoints_opt eigenvalue_band_proprties.
            # gap, cbm, vbm, is_direct = self.eigenvalue_band_properties
            # vout |= {"bandgap": gap, "cbm": cbm, "vbm": vbm, "is_gap_direct": is_direct}

            if kpt_opt_props.projected_eigenvalues:
                vout["projected_eigenvalues_kpoints_opt"] = {
                    str(spin): v.tolist() for spin, v in kpt_opt_props.projected_eigenvalues.items()
                }

        dct["output"] = vout
        return jsanitize(dct, strict=True)


class Outcar:
    """Parser for data in OUTCAR that is not available in Vasprun.xml.

    Note, this class works a bit differently than most of the other
    VASP objects, since OUTCAR can be very different depending on which
    "type of run" performed.

    Create the OUTCAR class with a filename reads "regular parameters" that
    are always present.

    Attributes:
        magnetization (tuple): Magnetization on each ion as a tuple of dict, e.g.
            ({"d": 0.0, "p": 0.003, "s": 0.002, "tot": 0.005}, ... )
        chemical_shielding (dict): Chemical shielding on each ion as a dictionary with core and valence contributions.
        unsym_cs_tensor (list): Unsymmetrized chemical shielding tensor matrixes on each ion as a list.
            e.g. [[[sigma11, sigma12, sigma13], [sigma21, sigma22, sigma23], [sigma31, sigma32, sigma33]], ...]
        cs_g0_contribution (np.array): G=0 contribution to chemical shielding. 2D rank 3 matrix.
        cs_core_contribution (dict): Core contribution to chemical shielding. dict. e.g.
            {'Mg': -412.8, 'C': -200.5, 'O': -271.1}
        efg (tuple): Electric Field Gradient (EFG) tensor on each ion as a tuple of dict, e.g.
            ({"cq": 0.1, "eta", 0.2, "nuclear_quadrupole_moment": 0.3}, {"cq": 0.7, "eta", 0.8,
            "nuclear_quadrupole_moment": 0.9}, ...)
        charge (tuple): Charge on each ion as a tuple of dict, e.g.
            ({"p": 0.154, "s": 0.078, "d": 0.0, "tot": 0.232}, ...)
        is_stopped (bool): True if OUTCAR is from a stopped run (using STOPCAR, see VASP Manual).
        run_stats (dict): Various useful run stats as a dict including "System time (sec)", "Total CPU time used (sec)",
            "Elapsed time (sec)", "Maximum memory used (kb)", "Average memory used (kb)", "User time (sec)", "cores".
        elastic_tensor (np.array): Total elastic moduli (Kbar) is given in a 6x6 array matrix.
        drift (np.array): Total drift for each step in eV/Atom.
        ngf (tuple): Dimensions for the Augmentation grid.
        sampling_radii (np.array): Size of the sampling radii in VASP for the test charges for the electrostatic
            potential at each atom. Total array size is the number of elements present in the calculation.
        electrostatic_potential (np.array): Average electrostatic potential at each atomic position in order of
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
        has_onsite_density_matrices (bool): Whether onsite density matrices have been set.
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

    def __init__(self, filename: PathLike) -> None:
        """
        Args:
            filename (PathLike): OUTCAR file to parse.
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
        run_stats: dict[str, float | None] = {}
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

                if match := efermi_patt.search(clean):
                    try:
                        # try-catch because VASP sometimes prints
                        # 'E-fermi: ********     XC(G=0):  -6.1327
                        # alpha+bet : -1.8238'
                        efermi = float(match[1])
                        continue
                    except ValueError:
                        efermi = None
                        continue
                if match := nelect_patt.search(clean):
                    nelect = float(match[1])

                if match := mag_patt.search(clean):
                    total_mag = float(match[1])

                if e_fr_energy is None and (match := e_fr_energy_pattern.search(clean)):
                    e_fr_energy = float(match[1])
                if e_wo_entrp is None and (match := e_wo_entrp_pattern.search(clean)):
                    e_wo_entrp = float(match[1])
                if e0 is None and (match := e0_pattern.search(clean)):
                    e0 = float(match[1])
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
                elif match := re.match(r"\s*(\d+)\s+(([\d\.\-]+)\s+)+", clean):
                    tokens = [float(i) for i in re.findall(r"[\d\.\-]+", clean)]
                    tokens.pop(0)
                    if read_charge:
                        charge.append(dict(zip(header, tokens, strict=True)))
                    elif read_mag_x:
                        mag_x.append(dict(zip(header, tokens, strict=True)))
                    elif read_mag_y:
                        mag_y.append(dict(zip(header, tokens, strict=True)))
                    elif read_mag_z:
                        mag_z.append(dict(zip(header, tokens, strict=True)))
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

        # Merge x, y and z components of magmoms if present (SOC calculation)
        if mag_y and mag_z:
            # TODO: detect spin axis
            mag = []
            for idx in range(len(mag_x)):
                mag.append({key: Magmom([mag_x[idx][key], mag_y[idx][key], mag_z[idx][key]]) for key in mag_x[0]})
        else:
            mag = mag_x

        # Data from beginning of OUTCAR
        run_stats["cores"] = None
        with zopen(filename, mode="rt") as file:
            for line in file:
                if "serial" in line:
                    # Activate serial parallelization
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
        self.data: dict = {}

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
                (
                    r"maximum number of plane-waves"
                    if serial_compilation
                    else r"maximum and minimum number of plane-waves"
                ),
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

        # Read the drift
        self.read_pattern(
            {"drift": r"total drift:\s+([\.\-\d]+)\s+([\.\-\d]+)\s+([\.\-\d]+)"},
            terminate_on_match=False,
            postprocess=float,
        )
        self.drift = self.data.get("drift", [])

        # Check if calculation is spin polarized
        self.read_pattern({"spin": r"ISPIN\s*=\s*2"})
        self.spin = bool(self.data.get("spin", []))

        # Check if calculation is non-collinear
        self.read_pattern({"noncollinear": r"LNONCOLLINEAR\s*=\s*T"})
        self.noncollinear = bool(self.data.get("noncollinear", []))

        # Check if the calculation type is DFPT
        self.read_pattern(
            {"ibrion": r"IBRION =\s+([\-\d]+)"},
            terminate_on_match=True,
            postprocess=int,
        )
        if self.data.get("ibrion", [[0]])[0][0] > 6:
            self.dfpt = True
            self.read_internal_strain_tensor()
        else:
            self.dfpt = False

        # Check if LEPSILON is True and read piezo data if so
        self.read_pattern({"epsilon": r"LEPSILON\s*=\s*T"})
        if self.data.get("epsilon", []):
            self.lepsilon = True
            self.read_lepsilon()
            # Only read ionic contribution if DFPT is turned on
            if self.dfpt:
                self.read_lepsilon_ionic()
        else:
            self.lepsilon = False

        # Check if LCALCPOL is True and read polarization data if so
        self.read_pattern({"calcpol": r"LCALCPOL\s*=\s*T"})
        if self.data.get("calcpol", []):
            self.lcalcpol = True
            self.read_lcalcpol()
            self.read_pseudo_zval()
        else:
            self.lcalcpol = False

        # Read electrostatic potential
        self.electrostatic_potential: list[float] | None = None
        self.ngf = None
        self.sampling_radii: list[float] | None = None
        self.read_pattern({"electrostatic": r"average \(electrostatic\) potential at core"})
        if self.data.get("electrostatic", []):
            self.read_electrostatic_potential()

        self.read_pattern({"nmr_cs": r"LCHIMAG\s*=\s*(T)"})
        if self.data.get("nmr_cs"):
            self.nmr_cs = True
            self.read_chemical_shielding()
            self.read_cs_g0_contribution()
            self.read_cs_core_contribution()
            self.read_cs_raw_symmetrized_tensors()
        else:
            self.nmr_cs = False

        self.read_pattern({"nmr_efg": r"NMR quadrupolar parameters"})
        if self.data.get("nmr_efg"):
            self.nmr_efg = True
            self.read_nmr_efg()
            self.read_nmr_efg_tensor()
        else:
            self.nmr_efg = False

        self.read_pattern(
            {"has_onsite_density_matrices": r"onsite density matrix"},
            terminate_on_match=True,
        )
        if "has_onsite_density_matrices" in self.data:
            self.has_onsite_density_matrices = True
            self.read_onsite_density_matrices()
        else:
            self.has_onsite_density_matrices = False

        # Store the individual contributions to the final total energy
        final_energy_contribs = {}
        for key in (
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
        ):
            if key == "PAW double counting":
                self.read_pattern({key: rf"{key}\s+=\s+([\.\-\d]+)\s+([\.\-\d]+)"})
            else:
                self.read_pattern({key: rf"{key}\s+=\s+([\d\-\.]+)"})
            if not self.data[key]:
                continue
            final_energy_contribs[key] = sum(map(float, self.data[key][-1]))
        self.final_energy_contribs = final_energy_contribs

    def read_pattern(
        self,
        patterns: dict[str, str],
        reverse: bool = False,
        terminate_on_match: bool = False,
        postprocess: Callable = str,
    ) -> None:
        r"""
        General pattern reading. Use monty's regrep method and take the same
        arguments.

        Args:
            patterns (dict): A dict of patterns, e.g.
                {"energy": r"energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)"}.
            reverse (bool): Read files in reverse. Defaults to false. Useful for
                large files, esp OUTCARs, especially when used with
                terminate_on_match.
            terminate_on_match (bool): Whether to terminate when there is at
                least one match in each key in pattern.
            postprocess (Callable): A post processing function to convert all
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
        header_pattern: str,
        row_pattern: str,
        footer_pattern: str,
        postprocess: Callable = str,
        attribute_name: str | None = None,
        last_one_only: bool = True,
        first_one_only: bool = False,
    ) -> list:
        r"""Parse table-like data. A table composes of three parts: header,
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
            postprocess (Callable): A post processing function to convert all
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

        with zopen(self.filename, mode="rt") as file:
            text = file.read()
        table_pattern_text = header_pattern + r"\s*^(?P<table_body>(?:\s+" + row_pattern + r")+)\s+" + footer_pattern
        table_pattern = re.compile(table_pattern_text, re.MULTILINE | re.DOTALL)
        rp = re.compile(row_pattern)
        tables: list[list] = []
        for mt in table_pattern.finditer(text):
            table_body_text = mt.group("table_body")
            table_contents = []
            for line in table_body_text.split("\n"):
                ml = rp.search(line)
                # Skip empty lines
                if not ml:
                    continue
                d = ml.groupdict()
                if len(d) > 0:
                    processed_line: dict | list = {k: postprocess(v) for k, v in d.items()}
                else:
                    processed_line = [postprocess(v) for v in ml.groups()]
                table_contents.append(processed_line)
            tables.append(table_contents)
            if first_one_only:
                break
        retained_data: list = tables[-1] if last_one_only or first_one_only else tables
        if attribute_name is not None:
            self.data[attribute_name] = retained_data
        return retained_data

    def read_electrostatic_potential(self) -> None:
        """Parse the eletrostatic potential for the last ionic step."""
        pattern = {"ngf": r"\s+dimension x,y,z NGXF=\s+([\.\-\d]+)\sNGYF=\s+([\.\-\d]+)\sNGZF=\s+([\.\-\d]+)"}
        self.read_pattern(pattern, postprocess=int)
        self.ngf = self.data.get("ngf", [[]])[0]

        pattern = {"radii": r"the test charge radii are((?:\s+[\.\-\d]+)+)"}
        self.read_pattern(pattern, reverse=True, terminate_on_match=True, postprocess=str)
        self.sampling_radii = [*map(float, self.data["radii"][0][0].split())]

        header_pattern = r"\(the norm of the test charge is\s+[\.\-\d]+\)"
        table_pattern = r"((?:\s+\d+\s*[\.\-\d]+)+)"
        footer_pattern = r"\s+E-fermi :"

        pots: list = self.read_table_pattern(header_pattern, table_pattern, footer_pattern)
        _pots: str = "".join(itertools.chain.from_iterable(pots))

        pots = re.findall(r"\s+\d+\s*([\.\-\d]+)+", _pots)

        self.electrostatic_potential = [*map(float, pots)]

    @staticmethod
    def _parse_sci_notation(line: str) -> list[float]:
        """
        Parse lines with values in scientific notation and potentially
        without spaces in between the values. This assumes that the scientific
        notation always lists two digits for the exponent, e.g. 3.535E-02.

        Args:
            line: line to parse.

        Returns:
            list[float]: numbers if found, empty if not.
        """
        if match := re.findall(r"[\.\-\d]+E[\+\-]\d{2}", line):
            return [float(t) for t in match]
        return []

    def read_freq_dielectric(self) -> None:
        """
        Parse the frequency dependent dielectric function (obtained with
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
        read_plasma: str | bool = False
        read_dielectric = False
        energies = []
        data: dict[str, Any] = {"REAL": [], "IMAGINARY": []}
        count = 0
        component = "IMAGINARY"
        with zopen(self.filename, mode="rt") as file:
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
                elif read_plasma and type(self)._parse_sci_notation(line):
                    plasma_frequencies[read_plasma].append(type(self)._parse_sci_notation(line))
                elif read_dielectric:
                    tokens = None
                    if re.match(row_pattern, line.strip()):
                        tokens = line.strip().split()
                    elif type(self)._parse_sci_notation(line.strip()):
                        tokens = type(self)._parse_sci_notation(line.strip())
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

    def read_chemical_shielding(self) -> None:
        """Parse the NMR chemical shieldings data. Only the second part "absolute, valence and core"
        will be parsed. And only the three right most field (ISO_SHIELDING, SPAN, SKEW) will be retrieved.

        Set self.data["chemical_shielding"] as:
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
        self.data["chemical_shielding"] = {
            "valence_only": cs_valence_only,
            "valence_and_core": cs_valence_and_core,
        }

    def read_cs_g0_contribution(self) -> None:
        """Parse the G0 contribution of NMR chemical shielding.

        Set self.data["cs_g0_contribution"] as:
            list[list]: G0 contribution matrix.
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

    def read_cs_core_contribution(self) -> None:
        """Parse the core contribution of NMR chemical shielding.

        Set self.data["cs_core_contribution"] as:
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

    def read_cs_raw_symmetrized_tensors(self) -> None:
        """Parse the matrix form of NMR tensor before corrected to table.

        Returns:
            nsymmetrized tensors list in the order of atoms.
        """
        header_pattern = r"\s+-{50,}\s+\s+Absolute Chemical Shift tensors\s+\s+-{50,}$"
        first_part_pattern = r"\s+UNSYMMETRIZED TENSORS\s+$"
        row_pattern = r"\s+".join([r"([-]?\d+\.\d+)"] * 3)
        unsym_footer_pattern = r"^\s+SYMMETRIZED TENSORS\s+$"

        with zopen(self.filename, mode="rt") as file:
            text = file.read()
        unsym_table_pattern_text = header_pattern + first_part_pattern + r"(?P<table_body>.+)" + unsym_footer_pattern
        table_pattern = re.compile(unsym_table_pattern_text, re.MULTILINE | re.DOTALL)
        row_pat = re.compile(row_pattern)

        if match := table_pattern.search(text):
            table_text = match["table_body"]
            micro_header_pattern = r"ion\s+\d+"
            micro_table_pattern_text = micro_header_pattern + r"\s*^(?P<table_body>(?:\s*" + row_pattern + r")+)\s+"
            micro_table_pattern = re.compile(micro_table_pattern_text, re.MULTILINE | re.DOTALL)
            unsym_tensors = []
            for mt in micro_table_pattern.finditer(table_text):
                table_body_text = mt.group("table_body")
                tensor_matrix = []
                for line in table_body_text.rstrip().split("\n"):
                    ml = row_pat.search(line)
                    if ml is None:
                        raise RuntimeError(f"failure to find pattern, {ml=}")
                    processed_line = [float(v) for v in ml.groups()]
                    tensor_matrix.append(processed_line)
                unsym_tensors.append(tensor_matrix)
            self.data["unsym_cs_tensor"] = unsym_tensors
        else:
            raise ValueError("NMR UNSYMMETRIZED TENSORS is not found")

    def read_nmr_efg_tensor(self) -> list[NDArray]:
        """Parses the NMR Electric Field Gradient Raw Tensors.

        Returns:
            A list of Electric Field Gradient Tensors in the order of Atoms from OUTCAR.
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

    def read_nmr_efg(self) -> None:
        """Parse the NMR Electric Field Gradient interpreted values.

        Set self.data["efg"] as:
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

    def read_elastic_tensor(self) -> None:
        """
        Parse the elastic tensor data.

        Set self.data["elastic_tensor"] as:
            6x6 array corresponding to the elastic tensor from the OUTCAR.
        """
        header_pattern = r"TOTAL ELASTIC MODULI \(kBar\)\s+Direction\s+([X-Z][X-Z]\s+)+\-+"
        row_pattern = r"[X-Z][X-Z]\s+" + r"\s+".join([r"(\-*[\.\d]+)"] * 6)
        footer_pattern = r"\-+"
        et_table = self.read_table_pattern(header_pattern, row_pattern, footer_pattern, postprocess=float)
        self.data["elastic_tensor"] = et_table

    def read_piezo_tensor(self) -> None:
        """Parse the piezo tensor data."""
        header_pattern = r"PIEZOELECTRIC TENSOR  for field in x, y, z\s+\(C/m\^2\)\s+([X-Z][X-Z]\s+)+\-+"
        row_pattern = r"[x-z]\s+" + r"\s+".join([r"(\-*[\.\d]+)"] * 6)
        footer_pattern = r"BORN EFFECTIVE"
        pt_table = self.read_table_pattern(header_pattern, row_pattern, footer_pattern, postprocess=float)
        self.data["piezo_tensor"] = pt_table

    def read_onsite_density_matrices(self) -> None:
        """Parse the onsite density matrices.

        Set self.data["onsite_density_matrices"] as:
            List with index corresponding to atom index in Structure.
        """
        # Matrix size will vary depending on if d or f orbitals are present.
        # Therefore regex assumes f, but filter out None values if d.
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

        # Filter out None
        spin1_component = [[[e for e in row if e is not None] for row in matrix] for matrix in spin1_component]

        # And repeat for Spin.down
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

    def read_corrections(
        self,
        reverse: bool = True,
        terminate_on_match: bool = True,
    ) -> None:
        """Read the dipol qudropol corrections into
        self.data["dipol_quadrupol_correction"].

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

    def read_neb(
        self,
        reverse: bool = True,
        terminate_on_match: bool = True,
    ) -> None:
        """
        Read NEB data. This only works with OUTCARs from both normal
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

    def read_igpar(self) -> None:
        """Read IGPAR.

        See VASP sections "LBERRY, IGPAR, NPPSTR, DIPOL" for info on
        what these are.

        Renders accessible:
            er_ev = e<r>_ev (dictionary with Spin.up/Spin.down as keys)
            er_bp = e<r>_bp (dictionary with Spin.up/Spin.down as keys)
            er_ev_tot = spin up + spin down summed
            er_bp_tot = spin up + spin down summed
            p_elc = spin up + spin down summed
            p_ion = spin up + spin down summed.
        """
        # Variables to be filled
        self.er_ev = {}  # dict (Spin.up/down) of array(3*float)
        self.er_bp = {}  # dict (Spin.up/down) of array(3*float)
        self.er_ev_tot = None  # array(3*float)
        self.er_bp_tot = None  # array(3*float)
        self.p_elec: int | None = None
        self.p_ion: int | None = None
        try:
            search = []

            # Non-spin cases
            def er_ev(results, match):
                results.er_ev[Spin.up] = np.array(map(float, match.groups()[1:4])) / 2
                results.er_ev[Spin.down] = results.er_ev[Spin.up]
                results.context = 2

            er_ev_pattern = r"^ *e<r>_ev=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)"
            search.append([er_ev_pattern, None, er_ev])

            def er_bp(results, match):
                results.er_bp[Spin.up] = np.array([float(match[i]) for i in range(1, 4)]) / 2
                results.er_bp[Spin.down] = results.er_bp[Spin.up]

            er_bp_pattern = r"^ *e<r>_bp=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)"
            search.append([er_bp_pattern, lambda results, _line: results.context == 2, er_bp])

            # Spin cases
            def er_ev_up(results, match):
                results.er_ev[Spin.up] = np.array([float(match[i]) for i in range(1, 4)])
                results.context = Spin.up

            spin1_ev_pattern = r"^.*Spin component 1 *e<r>_ev=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)"
            search.append([spin1_ev_pattern, None, er_ev_up])

            def er_bp_up(results, match):
                results.er_bp[Spin.up] = np.array([float(match[1]), float(match[2]), float(match[3])])

            spin_bp_pattern = r"^ *e<r>_bp=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)"
            search.append(
                [
                    spin_bp_pattern,
                    lambda results, _line: results.context == Spin.up,
                    er_bp_up,
                ]
            )

            def er_ev_dn(results, match):
                results.er_ev[Spin.down] = np.array([float(match[1]), float(match[2]), float(match[3])])
                results.context = Spin.down

            spin2_pattern = r"^.*Spin component 2 *e<r>_ev=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)"
            search.append([spin2_pattern, None, er_ev_dn])

            def er_bp_dn(results, match):
                results.er_bp[Spin.down] = np.array([float(match[i]) for i in range(1, 4)])

            e_r_bp_pattern = r"^ *e<r>_bp=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)"
            search.append(
                [
                    e_r_bp_pattern,
                    lambda results, _line: results.context == Spin.down,
                    er_bp_dn,
                ]
            )

            # Always present spin/non-spin
            def p_elc(results, match):
                results.p_elc = np.array([float(match[i]) for i in range(1, 4)])

            elec_dipole_moment_pattern = (
                r"^.*Total electronic dipole moment: *p\[elc\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)"
            )
            search.append([elec_dipole_moment_pattern, None, p_elc])

            def p_ion(results, match):
                results.p_ion = np.array([float(match[i]) for i in range(1, 4)])

            ionic_dipole_moment_pattern = (
                r"^.*ionic dipole moment: *p\[ion\]=\( *([-0-9.Ee+]*) *([-0-9.Ee+]*) *([-0-9.Ee+]*) *\)"
            )
            search.append([ionic_dipole_moment_pattern, None, p_ion])

            self.context = None
            self.er_ev = {Spin.up: None, Spin.down: None}
            self.er_bp = {Spin.up: None, Spin.down: None}

            micro_pyawk(self.filename, search, self)

            if self.er_ev[Spin.up] is not None and self.er_ev[Spin.down] is not None:
                self.er_ev_tot = self.er_ev[Spin.up] + self.er_ev[Spin.down]  # type: ignore[operator]

            if self.er_bp[Spin.up] is not None and self.er_bp[Spin.down] is not None:
                self.er_bp_tot = self.er_bp[Spin.up] + self.er_bp[Spin.down]  # type: ignore[operator]

        except Exception as exc:
            raise RuntimeError("IGPAR OUTCAR could not be parsed.") from exc

    def read_internal_strain_tensor(self):
        """Read the internal strain tensor and populates
        self.internal_strain_tensor with an array of voigt notation
        tensors for each site.
        """
        search = []

        def internal_strain_start(results, match: str) -> None:
            results.internal_strain_ion = int(match[1]) - 1
            results.internal_strain_tensor.append(np.zeros((3, 6)))

        search.append(
            [
                r"INTERNAL STRAIN TENSOR FOR ION\s+(\d+)\s+for displacements in x,y,z  \(eV/Angst\):",
                None,
                internal_strain_start,
            ]
        )

        def internal_strain_data(results, match: str) -> None:
            if match[1].lower() == "x":
                index = 0
            elif match[1].lower() == "y":
                index = 1
            elif match[1].lower() == "z":
                index = 2
            else:
                raise IndexError(f"Couldn't parse row index from symbol for internal strain tensor: {match[1]}")
            results.internal_strain_tensor[results.internal_strain_ion][index] = np.array(
                [float(match[i]) for i in range(2, 8)]
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

    def read_lepsilon(self) -> None:
        """Read a LEPSILON run.

        TODO: Document the actual variables.
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
                    [float(match[i]) for i in range(1, 4)]
                )
                results.dielectric_index += 1

            search.append(
                [
                    r"^ *([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+) *$",
                    lambda results, _line: (
                        results.dielectric_index >= 0 if results.dielectric_index is not None else None
                    ),
                    dielectric_data,
                ]
            )

            def dielectric_section_stop(results, match):
                results.dielectric_index = None

            search.append(
                [
                    r"-------------------------------------",
                    lambda results, _line: (
                        results.dielectric_index >= 1 if results.dielectric_index is not None else None
                    ),
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
                results.piezo_tensor[results.piezo_index, :] = np.array([float(match[i]) for i in range(1, 7)])
                results.piezo_index += 1

            search.append(
                [
                    r"^ *[xyz] +([-0-9.Ee+]+) +([-0-9.Ee+]+)"
                    r" +([-0-9.Ee+]+) *([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+)*$",
                    lambda results, _line: (results.piezo_index >= 0 if results.piezo_index is not None else None),
                    piezo_data,
                ]
            )

            def piezo_section_stop(results, _match):
                results.piezo_index = None

            search.append(
                [
                    r"-------------------------------------",
                    lambda results, _line: (results.piezo_index >= 1 if results.piezo_index is not None else None),
                    piezo_section_stop,
                ]
            )

            self.piezo_index = None
            self.piezo_tensor = np.zeros((3, 6))

            def born_section_start(results, _match):
                results.born_ion = -1

            search.append([r"BORN EFFECTIVE CHARGES ", None, born_section_start])

            def born_ion(results, match):
                results.born_ion = int(match[1]) - 1
                results.born.append(np.zeros((3, 3)))

            search.append(
                [
                    r"ion +([0-9]+)",
                    lambda results, _line: results.born_ion is not None,
                    born_ion,
                ]
            )

            def born_data(results, match):
                results.born[results.born_ion][int(match[1]) - 1, :] = np.array([float(match[i]) for i in range(2, 5)])

            search.append(
                [
                    r"^ *([1-3]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+)$",
                    lambda results, _line: (
                        results.born_ion >= 0 if results.born_ion is not None else results.born_ion
                    ),
                    born_data,
                ]
            )

            def born_section_stop(results, _match):
                results.born_ion = None

            search.append(
                [
                    r"-------------------------------------",
                    lambda results, _line: (
                        results.born_ion >= 1 if results.born_ion is not None else results.born_ion
                    ),
                    born_section_stop,
                ]
            )

            self.born_ion = None
            self.born: list | np.ndarray = []

            micro_pyawk(self.filename, search, self)

            self.born = np.array(self.born)

            self.dielectric_tensor = self.dielectric_tensor.tolist()
            self.piezo_tensor = self.piezo_tensor.tolist()

        except Exception as exc:
            raise RuntimeError("LEPSILON OUTCAR could not be parsed.") from exc

    def read_lepsilon_ionic(self) -> None:
        """Read the ionic component of a LEPSILON run.

        TODO: Document the actual variables.
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
                    lambda results, _line: (
                        results.dielectric_ionic_index == -1
                        if results.dielectric_ionic_index is not None
                        else results.dielectric_ionic_index
                    ),
                    dielectric_section_start2,
                ]
            )

            def dielectric_data(results, match):
                results.dielectric_ionic_tensor[results.dielectric_ionic_index, :] = np.array(
                    [float(match[i]) for i in range(1, 4)]
                )
                results.dielectric_ionic_index += 1

            search.append(
                [
                    r"^ *([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+) *$",
                    lambda results, _line: (
                        results.dielectric_ionic_index >= 0
                        if results.dielectric_ionic_index is not None
                        else results.dielectric_ionic_index
                    ),
                    dielectric_data,
                ]
            )

            def dielectric_section_stop(results, _match):
                results.dielectric_ionic_index = None

            search.append(
                [
                    r"-------------------------------------",
                    lambda results, _line: (
                        results.dielectric_ionic_index >= 1
                        if results.dielectric_ionic_index is not None
                        else results.dielectric_ionic_index
                    ),
                    dielectric_section_stop,
                ]
            )

            self.dielectric_ionic_index = None
            self.dielectric_ionic_tensor = np.zeros((3, 3))

            def piezo_section_start(results, _match):
                results.piezo_ionic_index = 0

            search.append(
                [
                    "PIEZOELECTRIC TENSOR IONIC CONTR  for field in x, y, z        ",
                    None,
                    piezo_section_start,
                ]
            )

            def piezo_data(results, match):
                results.piezo_ionic_tensor[results.piezo_ionic_index, :] = np.array(
                    [float(match[i]) for i in range(1, 7)]
                )
                results.piezo_ionic_index += 1

            search.append(
                [
                    r"^ *[xyz] +([-0-9.Ee+]+) +([-0-9.Ee+]+)"
                    r" +([-0-9.Ee+]+) *([-0-9.Ee+]+) +([-0-9.Ee+]+) +([-0-9.Ee+]+)*$",
                    lambda results, _line: (
                        results.piezo_ionic_index >= 0
                        if results.piezo_ionic_index is not None
                        else results.piezo_ionic_index
                    ),
                    piezo_data,
                ]
            )

            def piezo_section_stop(results, _match):
                results.piezo_ionic_index = None

            search.append(
                [
                    "-------------------------------------",
                    lambda results, _line: (
                        results.piezo_ionic_index >= 1
                        if results.piezo_ionic_index is not None
                        else results.piezo_ionic_index
                    ),
                    piezo_section_stop,
                ]
            )

            self.piezo_ionic_index = None
            self.piezo_ionic_tensor = np.zeros((3, 6))

            micro_pyawk(self.filename, search, self)

            self.dielectric_ionic_tensor = self.dielectric_ionic_tensor.tolist()
            self.piezo_ionic_tensor = self.piezo_ionic_tensor.tolist()

        except Exception as exc:
            raise RuntimeError("ionic part of LEPSILON OUTCAR could not be parsed.") from exc

    def read_lcalcpol(self) -> None:
        """Read the LCALCPOL.

        TODO: Document the actual variables.
        """
        self.p_elec = None
        self.p_sp1: int | None = None
        self.p_sp2: int | None = None
        self.p_ion = None

        try:
            search = []

            # Always present spin/non-spin
            def p_elec(results, match):
                results.p_elec = np.array(
                    [
                        float(match[1]),
                        float(match[2]),
                        float(match[3]),
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
                            float(match[1]),
                            float(match[2]),
                            float(match[3]),
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
                            float(match[1]),
                            float(match[2]),
                            float(match[3]),
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
                        float(match[1]),
                        float(match[2]),
                        float(match[3]),
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

            # Fix polarization units in new versions of VASP
            regex = r"^.*Ionic dipole moment: .*"
            search = [[regex, None, lambda x, y: x.append(y.group(0))]]
            results = micro_pyawk(self.filename, search, [])

            if "|e|" in results[0]:
                self.p_elec *= -1  # type: ignore[operator]
                self.p_ion *= -1  # type: ignore[operator]
                if self.spin and not self.noncollinear:
                    self.p_sp1 *= -1  # type: ignore[operator]
                    self.p_sp2 *= -1  # type: ignore[operator]

        except Exception as exc:
            raise RuntimeError("LCALCPOL OUTCAR could not be parsed.") from exc

    def read_pseudo_zval(self) -> None:
        """Create a pseudopotential ZVAL dictionary."""
        try:

            def atom_symbols(results, match):
                element_symbol = match[1]
                if not hasattr(results, "atom_symbols"):
                    results.atom_symbols = []
                results.atom_symbols.append(element_symbol.strip())

            def zvals(results, match):
                zvals = match[1]
                results.zvals = map(float, re.findall(r"-?\d+\.\d*", zvals))

            search: list[list] = []
            search.extend(
                (
                    ["(?<=VRHFIN =)(.*)(?=:)", None, atom_symbols],
                    ["^\\s+ZVAL.*=(.*)", None, zvals],
                )
            )

            micro_pyawk(self.filename, search, self)

            self.zval_dict = dict(zip(self.atom_symbols, self.zvals, strict=True))  # type: ignore[attr-defined]

            # Clean up
            del self.atom_symbols  # type: ignore[attr-defined]
            del self.zvals  # type: ignore[attr-defined]

        except Exception as exc:
            raise RuntimeError("ZVAL dict could not be parsed.") from exc

    def read_core_state_eigen(self) -> list[dict]:
        """Read the core state eigenenergies at each ionic step.

        Returns:
            A list of dict over the atom such as [{"AO":[core state eig]}].
            The core state eigenenergie list for each AO is over all ionic
            step.

        Example:
            The core state eigenenergie of the 2s AO of the 6th atom of the
            structure at the last ionic step is [5]["2s"][-1].
        """
        with zopen(self.filename, mode="rt") as foutcar:
            line = foutcar.readline()
            cl: list[dict] = []

            while line != "":
                line = foutcar.readline()

                if "NIONS =" in line:
                    natom = int(line.split("NIONS =")[1])
                    cl = [defaultdict(list) for _ in range(natom)]

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

    def read_avg_core_poten(self) -> list[list]:
        """Read the core potential at each ionic step.

        Returns:
            A list for each ionic step containing a list of the average core
            potentials for each atom: [[avg core pot]].

        Example:
            The average core potential of the 2nd atom of the structure at the
            last ionic step is: [-1][1]
        """
        with zopen(self.filename, mode="rt") as foutcar:
            line = foutcar.readline()
            aps: list[list[float]] = []
            while line != "":
                line = foutcar.readline()
                if "the norm of the test charge is" in line:
                    ap: list[float] = []
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

    def as_dict(self) -> dict:
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
            dct |= {
                "piezo_tensor": self.piezo_tensor,
                "dielectric_tensor": self.dielectric_tensor,
                "born": self.born,
            }

        if self.dfpt:
            dct["internal_strain_tensor"] = self.internal_strain_tensor

        if self.dfpt and self.lepsilon:
            dct |= {
                "piezo_ionic_tensor": self.piezo_ionic_tensor,
                "dielectric_ionic_tensor": self.dielectric_ionic_tensor,
            }

        if self.lcalcpol:
            dct |= {"p_elec": self.p_elec, "p_ion": self.p_ion}
            if self.spin and not self.noncollinear:
                dct |= {"p_sp1": self.p_sp1, "p_sp2": self.p_sp2}
            dct["zval_dict"] = self.zval_dict

        if self.nmr_cs:
            dct.update(
                nmr_cs={
                    "valence and core": self.data["chemical_shielding"]["valence_and_core"],
                    "valence_only": self.data["chemical_shielding"]["valence_only"],
                    "g0": self.data["cs_g0_contribution"],
                    "core": self.data["cs_core_contribution"],
                    "raw": self.data["unsym_cs_tensor"],
                }
            )

        if self.nmr_efg:
            dct.update(
                nmr_efg={
                    "raw": self.data["unsym_efg_tensor"],
                    "parameters": self.data["efg"],
                }
            )

        if self.has_onsite_density_matrices:
            # Cast Spin to str for consistency with electronic_structure
            # TODO: improve handling of Enum (de)serialization in monty
            onsite_density_matrices = [{str(k): v for k, v in d.items()} for d in self.data["onsite_density_matrices"]]
            dct["onsite_density_matrices"] = onsite_density_matrices

        return dct

    def read_fermi_contact_shift(self) -> None:
        """Read Fermi contact (isotropic) hyperfine coupling parameter.

        Output example:
        Fermi contact (isotropic) hyperfine coupling parameter (MHz)
        -------------------------------------------------------------
        ion      A_pw      A_1PS     A_1AE     A_1c      A_tot
        -------------------------------------------------------------
         1      -0.002    -0.002    -0.051     0.000    -0.052
         2      -0.002    -0.002    -0.051     0.000    -0.052
         3       0.056     0.056     0.321    -0.048     0.321
        -------------------------------------------------------------
        which corresponds to:
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
    """Container for volumetric data that allows
    for reading/writing with Poscar-type data.
    """

    @staticmethod
    def parse_file(filename: PathLike) -> tuple[Poscar, dict, dict]:
        """
        Parse a generic volumetric data file in the VASP like format.
        Used by subclasses for parsing files.

        Args:
            filename (PathLike): Path of file to parse.

        Returns:
            tuple[Poscar, dict, dict]: Poscar object, data dict, data_aug dict
        """
        poscar_read = False
        poscar_string: list[str] = []
        dataset: np.ndarray = np.zeros((1, 1, 1))
        all_dataset: list[np.ndarray] = []
        # for holding any strings in input that are not Poscar
        # or VolumetricData (typically augmentation charges)
        all_dataset_aug: dict[int, list[str]] = {}
        dim: list[int] = []
        dimline = ""
        read_dataset = False
        ngrid_pts = 0
        data_count = 0
        poscar = None
        with zopen(filename, mode="rt") as file:
            for line in file:
                original_line = line
                line = line.strip()
                if read_dataset:
                    for tok in line.split():
                        if data_count < ngrid_pts:
                            # This complicated procedure is necessary because
                            # VASP outputs x as the fastest index, followed by y
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

                # Construct a "diff" dict for scalar-like magnetization density,
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
            return poscar, data, data_aug  # type: ignore[return-value]

    def write_file(
        self,
        file_name: PathLike,
        vasp4_compatible: bool = False,
    ) -> None:
        """Write the VolumetricData object to a VASP compatible file.

        Args:
            file_name (PathLike): The output file.
            vasp4_compatible (bool): True if the format is VASP4 compatible.
        """

        def format_fortran_float(flt: float) -> str:
            """Fortran code prints floats with a leading zero in scientific
            notation. When writing CHGCAR files, we adopt this convention
            to ensure written CHGCAR files are byte-to-byte identical to
            their input files as far as possible.

            Args:
                flt (float): Float to print.

            Returns:
                str: The float in Fortran format.
            """
            flt_str = f"{flt:.10E}"
            if flt >= 0:
                return f"0.{flt_str[0]}{flt_str[2:12]}E{int(flt_str[13:]) + 1:+03}"
            return f"-.{flt_str[1]}{flt_str[3:13]}E{int(flt_str[14:]) + 1:+03}"

        def write_spin(data_type: str) -> None:
            lines = []
            count = 0
            file.write(f"   {dim[0]}   {dim[1]}   {dim[2]}\n")
            for k, j, i in itertools.product(list(range(dim[2])), list(range(dim[1])), list(range(dim[0]))):
                lines.append(format_fortran_float(self.data[data_type][i, j, k]))
                count += 1
                if count % 5 == 0:
                    file.write(" " + "".join(lines) + "\n")
                    lines = []
                else:
                    lines.append(" ")
            if count % 5 != 0:
                file.write(" " + "".join(lines) + " \n")

            data = self.data_aug.get(data_type, [])
            if isinstance(data, Iterable):
                file.write("".join(data))

        with zopen(file_name, mode="wt") as file:
            poscar = Poscar(self.structure)

            # Use original name if it's been set (e.g. from Chgcar)
            comment = getattr(self, "name", poscar.comment)

            lines = f"{comment}\n"
            lines += "   1.00000000000000\n"
            for vec in self.structure.lattice.matrix:
                lines += f" {vec[0]:12.6f}{vec[1]:12.6f}{vec[2]:12.6f}\n"
            if not vasp4_compatible:
                lines += "".join(f"{s:5}" for s in poscar.site_symbols) + "\n"
            lines += "".join(f"{x:6}" for x in poscar.natoms) + "\n"
            lines += "Direct\n"
            for site in self.structure:
                dim, b, c = site.frac_coords
                lines += f"{dim:10.6f}{b:10.6f}{c:10.6f}\n"
            lines += " \n"
            file.write(lines)
            dim = self.dim

            write_spin("total")
            if self.is_spin_polarized:
                if self.is_soc:
                    write_spin("diff_x")
                    write_spin("diff_y")
                    write_spin("diff_z")
                else:
                    write_spin("diff")


class Locpot(VolumetricData):
    """LOCPOT file reader."""

    def __init__(self, poscar: Poscar, data: np.ndarray, **kwargs) -> None:
        """
        Args:
            poscar (Poscar): Poscar object containing structure.
            data (np.ndarray): Actual data.
        """
        super().__init__(poscar.structure, data, **kwargs)
        self.name = poscar.comment

    @classmethod
    def from_file(cls, filename: PathLike, **kwargs) -> Self:
        """Read a LOCPOT file.

        Args:
            filename (PathLike): Path to LOCPOT file.

        Returns:
            Locpot
        """
        poscar, data, _data_aug = VolumetricData.parse_file(filename)
        return cls(poscar, data, **kwargs)


class Chgcar(VolumetricData):
    """CHGCAR file reader."""

    def __init__(
        self,
        poscar: Poscar | Structure,
        data: dict[str, NDArray],
        data_aug: NDArray | None = None,
    ) -> None:
        """
        Args:
            poscar (Poscar | Structure): Object containing structure.
            data: Actual data.
            data_aug: Augmentation charge data.
        """
        # Allow Poscar or Structure to be passed
        if isinstance(poscar, Poscar):
            struct = poscar.structure
            self.poscar = poscar
            self.name: str | None = poscar.comment  # type: ignore[assignment]
        elif isinstance(poscar, Structure):
            struct = poscar
            self.poscar = Poscar(poscar)
            self.name = None
        else:
            raise TypeError("Unsupported POSCAR type.")

        super().__init__(struct, data, data_aug=data_aug)
        self._distance_matrix: dict = {}

    @classmethod
    def from_file(cls, filename: str) -> Self:
        """Read a CHGCAR file.

        Args:
            filename (str): Path to CHGCAR file.

        Returns:
            Chgcar
        """
        poscar, data, data_aug = VolumetricData.parse_file(filename)
        return cls(poscar, data, data_aug=data_aug)

    @property
    def net_magnetization(self) -> float | None:
        """Net magnetic moment from Chgcar."""
        return np.sum(self.data["diff"]) if self.is_spin_polarized else None


class Elfcar(VolumetricData):
    """Read an ELFCAR file which contains the Electron Localization Function (ELF).

    For ELF, "total" key refers to Spin.up, and "diff" refers to Spin.down.

    This also contains information on the kinetic energy density.
    """

    def __init__(
        self,
        poscar: Poscar | Structure,
        data: dict[str, NDArray],
    ) -> None:
        """
        Args:
            poscar (Poscar or Structure): Object containing structure.
            data: Actual data.
        """
        # Allow Poscar or Structure to be passed
        if isinstance(poscar, Poscar):
            tmp_struct = poscar.structure
            self.poscar = poscar
        elif isinstance(poscar, Structure):
            tmp_struct = poscar
            self.poscar = Poscar(poscar)
        else:
            raise TypeError("Unsupported POSCAR type.")

        super().__init__(tmp_struct, data)
        # TODO: (mkhorton) modify VolumetricData so that the correct keys can be used.
        # for ELF, instead of "total" and "diff" keys we have
        # "Spin.up" and "Spin.down" keys
        # I believe this is correct, but there's not much documentation.
        self.data = data

    @classmethod
    def from_file(cls, filename: str) -> Self:
        """
        Read a ELFCAR file.

        Args:
            filename: Filename

        Returns:
            Elfcar
        """
        poscar, data, _data_aug = VolumetricData.parse_file(filename)
        return cls(poscar, data)

    def get_alpha(self) -> VolumetricData:
        """Get the parameter alpha where ELF = 1/(1 + alpha^2)."""
        alpha_data = {key: np.sqrt((1 / val) - 1) for key, val in self.data.items()}
        return VolumetricData(self.structure, alpha_data)


class Procar(MSONable):
    """
    PROCAR file reader.

    Updated to use code from easyunfold (https://smtg-bham.github.io/easyunfold; band-structure
    unfolding package) to allow SOC PROCAR parsing, and parsing multiple PROCAR files together.
    easyunfold's PROCAR parser can be used if finer control over projections (k-point weighting,
    normalisation per band, quick orbital sub-selection etc) is needed.

    Attributes:
        data (dict): The PROCAR data of the form below. It should VASP uses 1-based indexing,
            but all indices are converted to 0-based here.
            { spin: nd.array accessed with (k-point index, band index, ion index, orbital index) }
        weights (np.array): The weights associated with each k-point as an nd.array of length nkpoints.
        phase_factors (dict): Phase factors, where present (e.g. LORBIT = 12). A dict of the form:
            { spin: complex nd.array accessed with (k-point index, band index, ion index, orbital index) }
        nbands (int): Number of bands.
        nkpoints (int): Number of k-points.
        nions (int): Number of ions.
        nspins (int): Number of spins.
        is_soc (bool): Whether the PROCAR contains spin-orbit coupling (LSORBIT = True) data.
        kpoints (np.array): The k-points as an nd.array of shape (nkpoints, 3).
        occupancies (dict): The occupancies of the bands as a dict of the form:
            { spin: nd.array accessed with (k-point index, band index) }
        eigenvalues (dict): The eigenvalues of the bands as a dict of the form:
            { spin: nd.array accessed with (k-point index, band index) }
        xyz_data (dict): The PROCAR projections data along the x,y and z magnetisation projection
            directions, with is_soc = True (see VASP wiki for more info).
            { 'x'/'y'/'z': nd.array accessed with (k-point index, band index, ion index, orbital index) }
    """

    def __init__(self, filename: PathLike | list[PathLike]):
        """
        Args:
            filename: The path to PROCAR(.gz) file to read, or list of paths.
        """
        # get PROCAR filenames list to parse:
        filenames = filename if isinstance(filename, list) else [filename]
        self.nions: int | None = None  # used to check for consistency in files later
        self.nspins: int | None = None  # used to check for consistency in files later
        self.is_soc: bool | None = None  # used to check for consistency in files later
        self.orbitals = None  # used to check for consistency in files later
        self.read(filenames)

    def read(self, filenames: list[PathLike]):
        """
        Read in PROCAR projections data, possibly from multiple files.

        Args:
            filenames: List of PROCAR files to read.
        """
        parsed_kpoints = None
        occupancies_list, kpoints_list, weights_list = [], [], []
        eigenvalues_list, data_list, xyz_data_list = [], [], []
        phase_factors_list = []
        for filename in tqdm(filenames, desc="Reading PROCARs", unit="file", disable=len(filenames) == 1):
            (
                kpoints,
                weights,
                eigenvalues,
                occupancies,
                data,
                phase_factors,
                xyz_data,
            ) = self._read(filename, parsed_kpoints=parsed_kpoints)

            # Append to respective lists
            occupancies_list.append(occupancies)
            kpoints_list.append(kpoints)
            weights_list.append(weights)
            eigenvalues_list.append(eigenvalues)
            data_list.append(data)
            xyz_data_list.append(xyz_data)
            phase_factors_list.append(phase_factors)

        # Combine arrays along the kpoints axis:
        # nbands (axis = 2) could differ between arrays, so set missing values to zero:
        max_nbands = max(eig_dict[Spin.up].shape[1] for eig_dict in eigenvalues_list)
        for dict_array in itertools.chain(
            occupancies_list,
            eigenvalues_list,
            data_list,
            xyz_data_list,
            phase_factors_list,
        ):
            if dict_array:
                for key, array in dict_array.items():
                    if array.shape[1] < max_nbands:
                        if len(array.shape) == 2:  # occupancies, eigenvalues
                            dict_array[key] = np.pad(
                                array,
                                ((0, 0), (0, max_nbands - array.shape[2])),
                                mode="constant",
                            )
                        elif len(array.shape) == 4:  # data, phase_factors
                            dict_array[key] = np.pad(
                                array,
                                (
                                    (0, 0),
                                    (0, max_nbands - array.shape[2]),
                                    (0, 0),
                                    (0, 0),
                                ),
                                mode="constant",
                            )
                        elif len(array.shape) == 5:  # xyz_data
                            dict_array[key] = np.pad(
                                array,
                                (
                                    (0, 0),
                                    (0, max_nbands - array.shape[2]),
                                    (0, 0),
                                    (0, 0),
                                    (0, 0),
                                ),
                                mode="constant",
                            )
                        else:
                            raise ValueError("Unexpected array shape encountered!")

        # set nbands, nkpoints, and other attributes:
        self.nbands = max_nbands
        self.kpoints = np.concatenate(kpoints_list, axis=0)
        self.nkpoints = len(self.kpoints)
        self.occupancies = {
            spin: np.concatenate([occupancies[spin] for occupancies in occupancies_list], axis=0)
            for spin in occupancies_list[0]
        }
        self.eigenvalues = {
            spin: np.concatenate([eigenvalues[spin] for eigenvalues in eigenvalues_list], axis=0)
            for spin in eigenvalues_list[0]
        }
        self.weights = np.concatenate(weights_list, axis=0)
        self.data = {spin: np.concatenate([data[spin] for data in data_list], axis=0) for spin in data_list[0]}
        self.phase_factors = {
            spin: np.concatenate([phase_factors[spin] for phase_factors in phase_factors_list], axis=0)
            for spin in phase_factors_list[0]
        }
        if self.is_soc:
            self.xyz_data: dict | None = {
                key: np.concatenate([xyz_data[key] for xyz_data in xyz_data_list], axis=0) for key in xyz_data_list[0]
            }
        else:
            self.xyz_data = None

    def _parse_kpoint_line(self, line):
        """
        Parse k-point vector from a PROCAR line.

        Sometimes VASP outputs the kpoints joined together like
        '0.00000000-0.50000000-0.50000000' when there are negative signs,
        so need to be able to recognise and handle this.
        """
        fields = line.split()
        kpoint_fields = fields[3 : fields.index("weight")]
        kpoint_fields = [" -".join(field.split("-")).split() for field in kpoint_fields]
        kpoint_fields = [val for sublist in kpoint_fields for val in sublist]  # flatten

        return tuple(round(float(val), 5) for val in kpoint_fields)  # tuple to make it hashable,
        # rounded to 5 decimal places to ensure proper kpoint matching

    def _read(self, filename: PathLike, parsed_kpoints: set[tuple[Kpoint]] | None = None):
        """Main function for reading in the PROCAR projections data.

        Args:
            filename (PathLike): Path to PROCAR file to read.
            parsed_kpoints (set[tuple[Kpoint]]): Set of tuples of already-parsed kpoints (e.g. from multiple
                zero-weighted bandstructure calculations), to ensure redundant/duplicate parsing.
        """
        if parsed_kpoints is None:
            parsed_kpoints = set()

        with zopen(filename, mode="rt") as file_handle:
            preamble_expr = re.compile(r"# of k-points:\s*(\d+)\s+# of bands:\s*(\d+)\s+# of ions:\s*(\d+)")
            kpoint_expr = re.compile(r"^k-point\s+(\d+).*weight = ([0-9\.]+)")
            band_expr = re.compile(r"^band\s+(\d+)")
            ion_expr = re.compile(r"^ion.*")
            total_expr = re.compile(r"^tot.*")
            expr = re.compile(r"^([0-9]+)\s+")
            current_kpoint = 0
            current_band = 0
            spin = Spin.down  # switched to Spin.up for first block

            n_kpoints = None
            kpoints: list[tuple[float, float, float]] = []
            n_bands = None
            n_ions = None
            weights: np.ndarray[float] | None = None
            headers = None
            data: dict[Spin, np.ndarray] = {}
            eigenvalues: dict[Spin, np.ndarray] | None = None
            occupancies: dict[Spin, np.ndarray] | None = None
            phase_factors: dict[Spin, np.ndarray] | None = None
            xyz_data: dict[str, np.ndarray] | None = None  # 'x'/'y'/'z' as keys for SOC projections dict
            # keep track of parsed kpoints, to avoid redundant/duplicate parsing with multiple PROCARs:
            this_procar_parsed_kpoints = (
                set()
            )  # set of tuples of parsed (kvectors, 0/1 for Spin.up/down) for this PROCAR

            # first dynamically determine whether PROCAR is SOC or not; SOC PROCARs have 4 lists of projections (
            # total and x,y,z) for each band, while non-SOC have only 1 list of projections:
            tot_count = 0
            band_count = 0
            for line in file_handle:
                if total_expr.match(line):
                    tot_count += 1
                elif band_expr.match(line):
                    band_count += 1
                if band_count == 2:
                    break

            file_handle.seek(0)  # reset file handle to beginning
            if tot_count == 1:
                is_soc = False
            elif tot_count == 4:
                is_soc = True
            else:
                raise ValueError(
                    "Number of lines starting with 'tot' in PROCAR does not match expected values (4x or 1x number of "
                    "lines with 'band'), indicating a corrupted file!"
                )
            if self.is_soc is not None and self.is_soc != is_soc:
                raise ValueError("Mismatch in SOC setting (LSORBIT) in supplied PROCARs!")
            self.is_soc = is_soc

            skipping_kpoint = False  # true when skipping projections for a previously-parsed kpoint
            ion_line_count = 0  # printed twice when phase factors present
            proj_data_parsed_for_band = 0  # 0 for non-SOC, 1-4 for SOC/phase factors
            for line in file_handle:
                line = line.strip()
                if ion_expr.match(line):
                    ion_line_count += 1

                if kpoint_expr.match(line):
                    kvec = self._parse_kpoint_line(line)
                    match = kpoint_expr.match(line)
                    current_kpoint = int(match[1]) - 1  # type: ignore[index]
                    if current_kpoint == 0:
                        spin = Spin.up if spin == Spin.down else Spin.down

                    if (
                        kvec not in parsed_kpoints
                        and (kvec, {Spin.down: 0, Spin.up: 1}[spin]) not in this_procar_parsed_kpoints
                    ):
                        this_procar_parsed_kpoints.add((kvec, {Spin.down: 0, Spin.up: 1}[spin]))
                        skipping_kpoint = False
                        if spin == Spin.up:
                            kpoints.append(kvec)  # only add once
                    else:  # skip ahead to next kpoint:
                        skipping_kpoint = True
                        continue

                    if spin == Spin.up:  # record k-weight only once
                        weights[current_kpoint] = float(match[2])  # type: ignore[index]
                    proj_data_parsed_for_band = 0

                elif skipping_kpoint:
                    continue

                elif band_expr.match(line):
                    ion_line_count = 0  # printed a second time when phase factors present
                    match = band_expr.match(line)
                    current_band = int(match[1]) - 1  # type: ignore[index]
                    tokens = line.split()
                    eigenvalues[spin][current_kpoint, current_band] = float(tokens[4])  # type: ignore[index]
                    occupancies[spin][current_kpoint, current_band] = float(tokens[-1])  # type: ignore[index]
                    # keep track of parsed projections for each band (1x w/non-SOC, 4x w/SOC):
                    proj_data_parsed_for_band = 0

                elif headers is None and ion_expr.match(line):
                    headers = line.split()
                    headers.pop(0)
                    headers.pop(-1)

                    data = defaultdict(lambda: np.zeros((n_kpoints, n_bands, n_ions, len(headers))))
                    phase_factors = defaultdict(
                        lambda: np.full(
                            (n_kpoints, n_bands, n_ions, len(headers)),
                            np.nan,
                            dtype=np.complex128,
                        )
                    )
                    if self.is_soc:  # dict keys are now "x", "y", "z" rather than Spin.up/down
                        xyz_data = defaultdict(lambda: np.zeros((n_kpoints, n_bands, n_ions, len(headers))))

                elif expr.match(line):
                    tokens = line.split()
                    index = int(tokens.pop(0)) - 1
                    if headers is None:
                        raise ValueError("headers is None")
                    num_data = np.array([float(t) for t in tokens[: len(headers)]])
                    if phase_factors is None:
                        raise ValueError("phase_factors is None")

                    if proj_data_parsed_for_band == 0:
                        data[spin][current_kpoint, current_band, index, :] = num_data

                    elif self.is_soc and proj_data_parsed_for_band < 4:
                        proj_direction = {1: "x", 2: "y", 3: "z"}[proj_data_parsed_for_band]
                        if xyz_data is None:
                            raise ValueError(f"{xyz_data=}")
                        xyz_data[proj_direction][current_kpoint, current_band, index, :] = num_data

                    elif len(tokens) > len(headers):  # note no xyz projected phase factors with SOC
                        # New format of PROCAR (VASP 5.4.4)
                        num_data = np.array([float(t) for t in tokens[: 2 * len(headers)]])
                        for orb in range(len(headers)):
                            phase_factors[spin][current_kpoint, current_band, index, orb] = complex(
                                num_data[2 * orb], num_data[2 * orb + 1]
                            )
                    elif np.isnan(phase_factors[spin][current_kpoint, current_band, index, 0]):
                        # Old format of PROCAR (VASP 5.4.1 and before)
                        phase_factors[spin][current_kpoint, current_band, index, :] = num_data
                    else:
                        phase_factors[spin][current_kpoint, current_band, index, :] += 1j * num_data

                elif total_expr.match(line):
                    proj_data_parsed_for_band += 1

                elif preamble_expr.match(line):
                    match = preamble_expr.match(line)
                    if match is None:
                        raise RuntimeError(f"Failed to find preamable pattern, {match=}")
                    n_kpoints = int(match[1])
                    n_bands = int(match[2])
                    if eigenvalues is None:  # first spin
                        weights = np.zeros(n_kpoints)
                        eigenvalues = defaultdict(lambda: np.zeros((n_kpoints, n_bands)))
                        occupancies = defaultdict(lambda: np.zeros((n_kpoints, n_bands)))
                    n_ions = int(match[3])

                    if self.nions is not None and self.nions != n_ions:  # parsing multiple PROCARs but nions mismatch!
                        raise ValueError(f"Mismatch in number of ions in supplied PROCARs: ({n_ions} vs {self.nions})!")

            self.nions = n_ions  # attributes that should be consistent between multiple files are set here
            if self.orbitals is not None and self.orbitals != headers:  # multiple PROCARs but orbitals mismatch!
                raise ValueError(f"Mismatch in orbital headers in supplied PROCARs: {headers} vs {self.orbitals}!")
            self.orbitals = headers
            if self.nspins is not None and self.nspins != len(data):  # parsing multiple PROCARs but nspins mismatch!
                raise ValueError("Mismatch in number of spin channels in supplied PROCARs!")
            self.nspins = len(data)

            # chop off empty kpoints in arrays and redetermine nkpoints as we may have skipped previously-parsed kpoints
            nkpoints = current_kpoint + 1
            weights = np.array(weights[:nkpoints])  # type: ignore[index]
            data = {spin: data[spin][:nkpoints] for spin in data}  # type: ignore[index]
            eigenvalues = {spin: eigenvalues[spin][:nkpoints] for spin in eigenvalues}  # type: ignore[union-attr,index]
            occupancies = {spin: occupancies[spin][:nkpoints] for spin in occupancies}  # type: ignore[union-attr,index]
            phase_factors = {spin: phase_factors[spin][:nkpoints] for spin in phase_factors}  # type: ignore[union-attr,index]
            if self.is_soc:
                xyz_data = {spin: xyz_data[spin][:nkpoints] for spin in xyz_data}  # type: ignore[union-attr,index]

            # Update the parsed kpoints
            parsed_kpoints.update({kvec_spin_tuple[0] for kvec_spin_tuple in this_procar_parsed_kpoints})

            return (
                kpoints,
                weights,
                eigenvalues,
                occupancies,
                data,
                phase_factors,
                xyz_data,
            )

    def get_projection_on_elements(self, structure: Structure) -> dict[Spin, list[list[dict[str, float]]]]:
        """Get a dict of projections on elements.

        Args:
            structure (Structure): Input structure.

        Returns:
            A dict as {Spin: [band index][kpoint index][{Element: values}]].
        """
        if self.data is None:
            raise ValueError("data cannot be None.")
        if self.nkpoints is None:
            raise ValueError("nkpoints cannot be None.")
        if self.nbands is None:
            raise ValueError("nbands cannot be None.")
        if self.nions is None:
            raise ValueError("nions cannot be None.")

        elem_proj: dict[Spin, list] = {}
        for spin in self.data:
            elem_proj[spin] = [[defaultdict(float) for _ in range(self.nkpoints)] for _ in range(self.nbands)]

        for iat in range(self.nions):
            name = structure.species[iat].symbol
            for spin, data in self.data.items():
                for kpoint, band in itertools.product(range(self.nkpoints), range(self.nbands)):
                    elem_proj[spin][band][kpoint][name] += np.sum(data[kpoint, band, iat, :])

        return elem_proj

    def get_occupation(self, atom_index: int, orbital: str) -> dict:
        """Get the occupation for a particular orbital of a particular atom.

        Args:
            atom_num (int): Index of atom in the PROCAR. It should be noted
                that VASP uses 1-based indexing for atoms, but this is
                converted to 0-based indexing in this parser to be
                consistent with representation of structures in pymatgen.
            orbital (str): An orbital. If it is a single character, e.g. s,
                p, d or f, the sum of all s-type, p-type, d-type or f-type
                orbitals occupations are returned respectively. If it is a
                specific orbital, e.g. px, dxy, etc., only the occupation
                of that orbital is returned.

        Returns:
            Sum occupation of orbital of atom.
        """
        if self.orbitals is None:
            raise ValueError("orbitals is None")
        orbital_index = self.orbitals.index(orbital)

        if self.data is None:
            raise ValueError("data is None")
        return {
            spin: np.sum(data[:, :, atom_index, orbital_index] * self.weights[:, None])  # type: ignore[call-overload]
            for spin, data in self.data.items()
        }


class Oszicar:
    """OSZICAR parser for VASP.

    In general, while OSZICAR is useful for a quick look at the
    output from a VASP run, we recommend using the Vasprun parser
    instead, which gives far richer information.

    Attributes:
        electronic_steps (list): All electronic steps as a list of list of dict. e.g.
            [[{"rms": 160.0, "E": 4507.24605593, "dE": 4507.2, "N": 1, "deps": -17777.0, "ncg": 16576}, ...], [....]
            where electronic_steps[index] refers the list of electronic steps in one ionic_step,
            electronic_steps[index][subindex] refers to a particular electronic step at subindex in ionic step at
            index. The dict of properties depends on the type of VASP run, but in general, "E", "dE" and "rms" should
            be present in almost all runs.
        ionic_steps (list): All ionic_steps as a list of dict, e.g.
            [{"dE": -526.36, "E0": -526.36024, "mag": 0.0, "F": -526.36024}, ...]
            This is the typical output from VASP at the end of each ionic step. The stored dict might be different
            depending on the type of VASP run.
    """

    def __init__(self, filename: PathLike) -> None:
        """
        Args:
            filename (PathLike): The file to parse.
        """

        def smart_convert(header: str, num: float | str) -> float | str:
            try:
                return int(num) if header in {"N", "ncg"} else float(num)

            except ValueError:
                return "--"

        electronic_steps = []
        ionic_steps = []
        ionic_general_pattern = re.compile(r"(\w+)=\s*(\S+)")
        electronic_pattern = re.compile(r"\s*\w+\s*:(.*)")

        header: list = []
        with zopen(filename, mode="rt") as fid:
            for line in fid:
                if match := electronic_pattern.match(line.strip()):
                    tokens = match[1].split()
                    data = {header[idx]: smart_convert(header[idx], tokens[idx]) for idx in range(len(tokens))}
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
    def all_energies(self) -> tuple[tuple[float | str, ...], ...]:
        """Compilation of all energies from all electronic steps and ionic steps
        as a tuple of list of energies, e.g.
        ((4507.24605593, 143.824705755, -512.073149912, ...), ...).
        """
        all_energies: list[tuple] = []
        for i in range(len(self.electronic_steps)):
            energies: list[float | str] = [step["E"] for step in self.electronic_steps[i]]
            energies.append(self.ionic_steps[i]["F"])
            all_energies.append(tuple(energies))
        return tuple(all_energies)

    @property
    @unitized("eV")
    def final_energy(self) -> float:
        """Final energy from a run."""
        return self.ionic_steps[-1]["E0"]

    def as_dict(self) -> dict[str, list]:
        """MSONable dict."""
        return {
            "electronic_steps": self.electronic_steps,
            "ionic_steps": self.ionic_steps,
        }


class VaspParseError(ParseError):
    """Exception class for VASP parsing."""


def get_band_structure_from_vasp_multiple_branches(
    dir_name: str,
    efermi: float | None = None,
    projections: bool = False,
) -> BandStructureSymmLine | BandStructure | None:
    """Get band structure info from a VASP directory.

    It takes into account that a run can be divided in several branches named
    "branch_x". If the run has not been divided in branches the method will
    turn to parsing vasprun.xml directly.

    Args:
        dir_name: Directory containing all bandstructure runs.
        efermi: Efermi for bandstructure.
        projections: True if you want to get the data on site projections if
            any. Note that this is sometimes very large

    Returns:
        A BandStructure Object.
        None is there's a parsing error.
    """
    # TODO: Add better error handling!!!
    if os.path.isfile(f"{dir_name}/branch_0"):
        # Get all branch dir names
        branch_dir_names = [os.path.abspath(d) for d in glob(f"{dir_name}/branch_*") if os.path.isdir(d)]

        # Sort by the directory name (e.g, branch_10)
        sorted_branch_dir_names = sorted(branch_dir_names, key=lambda x: int(x.split("_")[-1]))

        # Populate branches with Bandstructure instances
        branches = []
        for dname in sorted_branch_dir_names:
            xml_file = f"{dname}/vasprun.xml"
            if os.path.isfile(xml_file):
                run = Vasprun(xml_file, parse_projected_eigen=projections)
                branches.append(run.get_band_structure(efermi=efermi))
            else:
                # TODO: It might be better to throw an exception
                warnings.warn(f"Skipping {dname}. Unable to find {xml_file}")

        return get_reconstructed_band_structure(branches, efermi)

    xml_file = f"{dir_name}/vasprun.xml"
    # Better handling of Errors
    if os.path.isfile(xml_file):
        return Vasprun(xml_file, parse_projected_eigen=projections).get_band_structure(
            kpoints_filename=None, efermi=efermi
        )

    return None


class Xdatcar:
    """XDATCAR parser. Only tested with VASP 5.x files.

    Attributes:
        structures (list): List of structures parsed from XDATCAR.
        comment (str): Optional comment string.

    Authors: Ram Balachandran
    """

    def __init__(
        self,
        filename: PathLike,
        ionicstep_start: int = 1,
        ionicstep_end: int | None = None,
        comment: str | None = None,
    ) -> None:
        """
        Init a Xdatcar.

        Args:
            filename (PathLike): The XDATCAR file.
            ionicstep_start (int): Starting number of ionic step.
            ionicstep_end (int): Ending number of ionic step.
            comment (str): Optional comment attached to this set of structures.
        """
        preamble = None
        coords_str: list = []
        structures: list = []
        preamble_done: bool = False
        parse_poscar: bool = False
        num_sites: int | None = None
        restart_preamble: bool = False
        if ionicstep_start < 1:
            raise ValueError("Start ionic step cannot be less than 1")
        if ionicstep_end is not None and ionicstep_end < 1:
            raise ValueError("End ionic step cannot be less than 1")

        file_len = sum(1 for _ in zopen(filename, mode="rt"))
        ionicstep_cnt = 1
        ionicstep_start = ionicstep_start or 0
        with zopen(filename, mode="rt") as file:
            title = None
            for iline, line in enumerate(file):
                line = line.strip()
                if preamble is None:
                    preamble = [line]
                    title = line

                elif title == line and len(coords_str) > 0:
                    # sometimes the title is the same as the only chemical species in the structure
                    # only enter this block if the coords have been read
                    parse_poscar = True
                    restart_preamble = True
                    preamble_done = False

                elif not preamble_done:
                    if line == "" or "Direct configuration=" in line:
                        preamble_done = True
                    else:
                        preamble.append(line)

                elif line == "" or "Direct configuration=" in line and len(coords_str) > 0:
                    parse_poscar = True
                    restart_preamble = False
                else:
                    coords_str.append(line)

                if (parse_poscar and (num_sites is None or len(coords_str) == num_sites)) or iline == file_len - 1:
                    if num_sites is None:
                        num_sites = len(coords_str)

                    poscar = Poscar.from_str("\n".join([*preamble, "Direct", *coords_str]))
                    if (ionicstep_end is None and ionicstep_cnt >= ionicstep_start) or (
                        ionicstep_end is not None and ionicstep_start <= ionicstep_cnt < ionicstep_end
                    ):
                        structures.append(poscar.structure)
                    elif (ionicstep_end is not None) and ionicstep_cnt >= ionicstep_end:
                        break

                    ionicstep_cnt += 1
                    coords_str = []
                    parse_poscar = False
                    if restart_preamble:
                        preamble = [line]

            if preamble is None:
                raise ValueError("preamble is None")

        self.structures = structures
        self.comment = comment or self.structures[0].formula

    def __str__(self) -> str:
        return self.get_str()

    @property
    def site_symbols(self) -> list[str]:
        """Sequence of symbols associated with the Xdatcar.
        Similar to 6th line in VASP 5+ Xdatcar.
        """
        syms = [site.specie.symbol for site in self.structures[0]]
        return [a[0] for a in itertools.groupby(syms)]

    @property
    def natoms(self) -> list[int]:
        """Sequence of number of sites of each type associated with the Poscar.
        Similar to 7th line in VASP 5+ Xdatcar.
        """
        syms = [site.specie.symbol for site in self.structures[0]]
        return [len(tuple(a[1])) for a in itertools.groupby(syms)]

    def concatenate(
        self,
        filename: PathLike,
        ionicstep_start: int = 1,
        ionicstep_end: int | None = None,
    ) -> None:
        """Concatenate structures in file to Xdatcar.

        Args:
            filename (PathLike): The XDATCAR file to be concatenated.
            ionicstep_start (int): Starting number of ionic step.
            ionicstep_end (int): Ending number of ionic step.

        TODO (rambalachandran): Check to ensure the new concatenating file
            has the same lattice structure and atoms as the Xdatcar class.
        """
        preamble = None
        coords_str: list[str] = []
        structures = self.structures
        preamble_done = False
        if ionicstep_start < 1:
            raise ValueError("Start ionic step cannot be less than 1")
        if ionicstep_end is not None and ionicstep_end < 1:
            raise ValueError("End ionic step cannot be less than 1")

        ionicstep_cnt = 1
        with zopen(filename, mode="rt") as file:
            for line in file:
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
                    poscar = Poscar.from_str("\n".join([*preamble, "Direct", *coords_str]))
                    if (
                        ionicstep_end is None and ionicstep_cnt >= ionicstep_start
                    ) or ionicstep_start <= ionicstep_cnt < ionicstep_end:
                        structures.append(poscar.structure)
                    ionicstep_cnt += 1
                    coords_str = []
                else:
                    coords_str.append(line)

            if preamble is None:
                raise ValueError("preamble is None")
            poscar = Poscar.from_str("\n".join([*preamble, "Direct", *coords_str]))

            if (
                (ionicstep_end is None and ionicstep_cnt >= ionicstep_start)
                or ionicstep_start <= ionicstep_cnt < ionicstep_end  # type: ignore[operator]
            ):
                structures.append(poscar.structure)
        self.structures = structures

    def get_str(
        self,
        ionicstep_start: int = 1,
        ionicstep_end: int | None = None,
        significant_figures: int = 8,
    ) -> str:
        """Write Xdatcar to a string.

        Args:
            ionicstep_start (int): Starting number of ionic step.
            ionicstep_end (int): Ending number of ionic step.
            significant_figures (int): Number of significant digits.
        """
        if ionicstep_start < 1:
            raise ValueError("Start ionic step cannot be less than 1")
        if ionicstep_end is not None and ionicstep_end < 1:
            raise ValueError("End ionic step cannot be less than 1")

        lattice = self.structures[0].lattice
        if np.linalg.det(lattice.matrix) < 0:
            lattice = Lattice(-lattice.matrix)
        lines = [self.comment, "1.0", str(lattice)]
        lines.extend((" ".join(self.site_symbols), " ".join(map(str, self.natoms))))

        format_str = f"{{:.{significant_figures}f}}"
        ionicstep_cnt = 1
        output_cnt = 1
        for cnt, structure in enumerate(self.structures, start=1):
            ionicstep_cnt = cnt
            if (
                (ionicstep_end is None and ionicstep_cnt >= ionicstep_start)
                or ionicstep_start <= ionicstep_cnt < ionicstep_end  # type: ignore[operator]
            ):
                lines.append(f"Direct configuration={' ' * (7 - len(str(output_cnt)))}{output_cnt}")
                for site in structure:
                    coords = site.frac_coords
                    line = " ".join(format_str.format(c) for c in coords)
                    lines.append(line)
                output_cnt += 1
        return "\n".join(lines) + "\n"

    def write_file(self, filename: PathLike, **kwargs) -> None:
        """Write Xdatcar into a file.

        Args:
            filename (PathLike): The output XDATCAR file.
            **kwargs: The same as those for the Xdatcar.get_str
                method and are passed through directly.
        """
        with zopen(filename, mode="wt") as file:
            file.write(self.get_str(**kwargs))


class Dynmat:
    """DYNMAT file reader.

    Attributes:
        data (dict): A nested dict containing the DYNMAT data of the form:
            [atom <int>][disp <int>]['dispvec'] =
                displacement vector (part of first line in dynmat block, e.g. "0.01 0 0")
            [atom <int>][disp <int>]['dynmat'] =
                    <list> list of dynmat lines for this atom and this displacement

    Authors: Patrick Huck
    """

    def __init__(self, filename: PathLike) -> None:
        """
        Args:
            filename: Name of file containing DYNMAT.
        """
        with zopen(filename, mode="rt") as file:
            lines = list(clean_lines(file.readlines()))
            self._nspecs, self._natoms, self._ndisps = map(int, lines[0].split())
            self._masses = map(float, lines[1].split())
            self.data: dict[int, dict] = defaultdict(dict)
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
                    if "dynmat" not in self.data[atom][disp]:  # type: ignore[index]
                        self.data[atom][disp]["dynmat"] = []  # type: ignore[index]
                    self.data[atom][disp]["dynmat"].append(v)  # type: ignore[index]

    def get_phonon_frequencies(self) -> list:
        """Calculate phonon frequencies.

        WARNING: This method is most likely incorrect or suboptimal,
        hence for demonstration purposes only.
        """
        frequencies = []
        for k, v0 in self.data.items():
            for v1 in v0.values():
                vec = map(abs, v1["dynmat"][k - 1])
                frequency = math.sqrt(sum(vec)) * 2.0 * math.pi * 15.633302  # THz
                frequencies.append(frequency)
        return frequencies

    @property
    def nspecs(self) -> int:
        """The number of species."""
        return self._nspecs

    @property
    def natoms(self) -> int:
        """The number of atoms."""
        return self._natoms

    @property
    def ndisps(self) -> int:
        """The number of displacements."""
        return self._ndisps

    @property
    def masses(self) -> list[float]:
        """The list of atomic masses."""
        return list(self._masses)


def get_adjusted_fermi_level(
    efermi: float,
    cbm: float,
    band_structure: BandStructureSymmLine,
    energy_step: float = 0.01,
) -> float:
    """
    When running a band structure computation, the Fermi level needs to be
    taken from the static run that gave the charge density used for the non-self
    consistent band structure run. Sometimes this Fermi level is too low
    because of the mismatch between the uniform grid used in the static run
    and the band structure k-points (e.g., the VBM is on Gamma and the Gamma
    point is not in the uniform mesh).

    Here we use a procedure looking for energy levels higher than the static
    Fermi level and lower than the LUMO. If any of these levels make the
    band structure appears insulating (not metallic anymore), we keep this
    adjusted fermi level.

    This procedure has shown to detect most insulators correctly.

    Args:
        efermi (float): The Fermi energy of the static run.
        cbm (float): The conduction band minimum of the static run.
        band_structure (BandStructureSymmLine): A band structure object.
        energy_step (float): The step length for increasing energy during search.

    Returns:
        float: A new adjusted Fermi level.
    """
    # Make a working copy of band_structure
    bs_working = BandStructureSymmLine.from_dict(band_structure.as_dict())
    if bs_working.is_metal():
        energy = efermi
        while energy < cbm:
            energy += energy_step
            bs_working._efermi = energy
            if not bs_working.is_metal():
                return energy
    return efermi


# A note to future confused people (i.e. myself):
# I use numpy.fromfile instead of scipy.io.FortranFile here because the records
# are of fixed length, so the record length is only written once. In fortran,
# this amounts to using open(..., form='unformatted', recl=recl_len). In
# contrast when you write UNK files, the record length is written at the
# beginning of each record. This allows you to use scipy.io.FortranFile. In
# fortran, this amounts to using open(..., form='unformatted') [i.e. no recl=].
class Wavecar:
    """
    Container for the (pseudo-) wavefunctions from VASP.

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
        vasp_type (str): String that determines VASP type the WAVECAR was generated with.
            One of 'std', 'gam', 'ncl'.
        nk (int): Number of k-points from the WAVECAR.
        nb (int): Number of bands per k-point.
        encut (float): Energy cutoff (used to define G_{cut}).
        efermi (float): Fermi energy.
        a (np.array): Primitive lattice vectors of the cell (e.g. a_1 = self.a[0, :]).
        b (np.array): Reciprocal lattice vectors of the cell (e.g. b_1 = self.b[0, :]).
        vol (float): The volume of the unit cell in real space.
        kpoints (np.array): The list of k-points read from the WAVECAR file.
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

    def __init__(
        self,
        filename: PathLike = "WAVECAR",
        verbose: bool = False,
        precision: Literal["normal", "accurate"] = "normal",
        vasp_type: Literal["std", "gam", "ncl"] | None = None,
    ) -> None:
        """Extract information from the given WAVECAR.

        Args:
            filename (PathLike): input file. Defaults to WAVECAR.
            verbose (bool): determines whether processing information is shown
            precision (str): determines how fine the fft mesh is (normal or
                accurate), only the first letter matters.
            vasp_type (str): determines the VASP type that is used, allowed
                values are {'std', 'gam', 'ncl'} (only first letter is required).
        """
        self.filename = filename
        valid_types = {"std", "gam", "ncl"}
        initials = {x[0] for x in valid_types}
        if vasp_type is not None and vasp_type.lower()[0] not in initials:
            raise ValueError(
                f"invalid {vasp_type=}, must be one of {valid_types} (we only check the first letter {initials})"
            )
        self.vasp_type = vasp_type

        # c = 0.26246582250210965422
        # 2m/hbar^2 in agreement with VASP
        self._C = 0.262465831
        with open(self.filename, "rb") as file:
            # Read the header information
            recl, spin, rtag = np.fromfile(file, dtype=np.float64, count=3).astype(int)
            if verbose:
                print(f"{recl=}, {spin=}, {rtag=}")
            recl8 = int(recl / 8)
            self.spin = spin

            # Make sure we have correct precision
            valid_rtags = {45200, 45210, 53300, 53310}
            if rtag not in valid_rtags:
                # note that rtag=45200 and 45210 may not work if file was actually
                # generated by old version of VASP, since that would write eigenvalues
                # and occupations in way that does not span FORTRAN records, but
                # reader below appears to assume that record boundaries can be ignored
                # (see OUTWAV vs. OUTWAV_4 in VASP fileio.F)
                raise ValueError(f"Invalid {rtag=}, must be one of {valid_rtags}")

            # Pad to end of fortran REC=1
            np.fromfile(file, dtype=np.float64, count=recl8 - 3)

            # Extract kpoint, bands, energy, and lattice information
            self.nk, self.nb = np.fromfile(file, dtype=np.float64, count=2).astype(int)
            self.encut = np.fromfile(file, dtype=np.float64, count=1)[0]
            self.a = np.fromfile(file, dtype=np.float64, count=9).reshape((3, 3))
            self.efermi = np.fromfile(file, dtype=np.float64, count=1)[0]
            if verbose:
                print(
                    f"kpoints = {self.nk}, bands = {self.nb}, energy cutoff = {self.encut}, fermi "
                    f"energy= {self.efermi:.04f}\n"
                )
                print(f"primitive lattice vectors = \n{self.a}")

            self.vol = np.dot(self.a[0, :], np.cross(self.a[1, :], self.a[2, :]))
            if verbose:
                print(f"volume = {self.vol}\n")

            # Calculate reciprocal lattice
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

            # Calculate maximum number of b vectors in each direction
            self._generate_nbmax()
            if verbose:
                print(f"max number of G values = {self._nbmax}\n\n")
            self.ng = self._nbmax * 3 if precision.lower()[0] == "n" else self._nbmax * 4

            # Pad to end of fortran REC=2
            np.fromfile(file, dtype=np.float64, count=recl8 - 13)

            # Read records
            self.Gpoints = [None for _ in range(self.nk)]
            self.kpoints = []
            if spin == 2:
                self.coeffs: list[list[list[None]]] | list[list[None]] = [
                    [[None for _ in range(self.nb)] for _ in range(self.nk)] for _ in range(spin)
                ]
                self.band_energy: list = [[] for _ in range(spin)]
            else:
                self.coeffs = [[None for _ in range(self.nb)] for _ in range(self.nk)]
                self.band_energy = []

            for i_spin in range(spin):
                if verbose:
                    print(f"Reading spin {i_spin}")

                for i_nk in range(self.nk):
                    # Information for this kpoint
                    nplane = int(np.fromfile(file, dtype=np.float64, count=1)[0])
                    kpoint = np.fromfile(file, dtype=np.float64, count=3)

                    if i_spin == 0:
                        self.kpoints.append(kpoint)
                    else:
                        assert_allclose(self.kpoints[i_nk], kpoint)

                    if verbose:
                        print(f"kpoint {i_nk: 4} with {nplane: 5} plane waves at {kpoint}")

                    # Energy and occupation information
                    enocc = np.fromfile(file, dtype=np.float64, count=3 * self.nb).reshape((self.nb, 3))
                    if spin == 2:
                        self.band_energy[i_spin].append(enocc)
                    else:
                        self.band_energy.append(enocc)

                    if verbose:
                        print("enocc =\n", enocc[:, [0, 2]])

                    # Pad the end of record that contains nplane, kpoints, evals and occs
                    np.fromfile(file, dtype=np.float64, count=(recl8 - 4 - 3 * self.nb) % recl8)

                    if self.vasp_type is None:
                        self.Gpoints[i_nk], extra_gpoints, extra_coeff_inds = self._generate_G_points(  # type: ignore[call-overload]
                            kpoint, gamma=True
                        )
                        if len(self.Gpoints[i_nk]) == nplane:  # type: ignore[arg-type]
                            self.vasp_type = "gam"
                        else:
                            self.Gpoints[i_nk], extra_gpoints, extra_coeff_inds = self._generate_G_points(  # type: ignore[call-overload]
                                kpoint, gamma=False
                            )
                            self.vasp_type = "std" if len(self.Gpoints[i_nk]) == nplane else "ncl"  # type: ignore[arg-type]

                        if verbose:
                            print(f"\ndetermined {self.vasp_type = }\n")
                    else:
                        self.Gpoints[i_nk], extra_gpoints, extra_coeff_inds = self._generate_G_points(  # type: ignore[call-overload]
                            kpoint, gamma=self.vasp_type.lower()[0] == "g"
                        )

                    if len(self.Gpoints[i_nk]) != nplane and 2 * len(self.Gpoints[i_nk]) != nplane:  # type: ignore[arg-type]
                        raise ValueError(
                            f"Incorrect {vasp_type=}. Please open an issue if you are certain this WAVECAR"
                            " was generated with the given vasp_type."
                        )

                    self.Gpoints[i_nk] = np.array(self.Gpoints[i_nk] + extra_gpoints, dtype=np.float64)  # type: ignore[arg-type, operator]

                    # Extract coefficients
                    for inb in range(self.nb):
                        if rtag in (45200, 53300):
                            data = np.fromfile(file, dtype=np.complex64, count=nplane)
                            np.fromfile(file, dtype=np.float64, count=recl8 - nplane)
                        elif rtag in (45210, 53310):
                            # TODO: This should handle double precision coefficients,
                            # but I don't have a WAVECAR to test it with
                            data = np.fromfile(file, dtype=np.complex128, count=nplane)
                            np.fromfile(file, dtype=np.float64, count=recl8 - 2 * nplane)
                        else:
                            raise RuntimeError("Invalid rtag value.")

                        extra_coeffs = []
                        if len(extra_coeff_inds) > 0:
                            # Reconstruct extra coefficients missing from gamma-only executable WAVECAR
                            for G_ind in extra_coeff_inds:
                                # No idea where this factor of sqrt(2) comes from,
                                # but empirically it appears to be necessary
                                data[G_ind] /= np.sqrt(2)
                                extra_coeffs.append(np.conj(data[G_ind]))

                        if spin == 2:
                            self.coeffs[i_spin][i_nk][inb] = np.array(list(data) + extra_coeffs, dtype=np.complex64)  # type: ignore[index]
                        else:
                            self.coeffs[i_nk][inb] = np.array(list(data) + extra_coeffs, dtype=np.complex128)

                        if self.vasp_type is not None and self.vasp_type.lower()[0] == "n":
                            self.coeffs[i_nk][inb].shape = (2, nplane // 2)  # type: ignore[union-attr]

    def _generate_nbmax(self) -> None:
        """Helper function to determine maximum number of b vectors for
        each direction.

        This algorithm is adapted from WaveTrans (see Class docstring).
        """
        bmag = np.linalg.norm(self.b, axis=1)
        b = self.b

        # Calculate maximum integers in each direction for G
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

    def _generate_G_points(
        self,
        kpoint: np.ndarray,
        gamma: bool = False,
    ) -> tuple[list, list, list]:
        """Helper method to generate G-points based on nbmax.

        This function iterates over possible G-point values and determines
        if the energy is less than G_{cut}. Valid values are appended to
        the output array. This function should not be called outside of
        initialization.

        Args:
            kpoint (np.array): the array containing the current k-point value
            gamma (bool): determines if G points for gamma-point only executable
                          should be generated

        Returns:
            A tuple containing valid G-points
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

    def evaluate_wavefunc(
        self,
        kpoint: int,
        band: int,
        r: np.ndarray,
        spin: int = 0,
        spinor: int = 0,
    ) -> np.complex64:
        r"""Evaluate the wavefunction for a given position, r.

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
            kpoint (int): the index of the kpoint where the wavefunction will be evaluated.
            band (int): the index of the band where the wavefunction will be evaluated.
            r (np.array): the position where the wavefunction will be evaluated.
            spin (int): spin index for the desired wavefunction (only for
                ISPIN = 2, default = 0).
            spinor (int): component of the spinor that is evaluated (only used
                if vasp_type == 'ncl').

        Returns:
            A complex value corresponding to the evaluation of the wavefunction.
        """
        if self.vasp_type is None:
            raise RuntimeError("vasp_type cannot be None.")

        v = self.Gpoints[kpoint] + self.kpoints[kpoint]
        u = np.dot(np.dot(v, self.b), r)

        if self.vasp_type.lower()[0] == "n":
            c = self.coeffs[kpoint][band][spinor, :]  # type: ignore[call-overload, index]
        elif self.spin == 2:
            c = self.coeffs[spin][kpoint][band]  # type: ignore[index]
        else:
            c = self.coeffs[kpoint][band]
        return np.sum(np.dot(c, np.exp(1j * u, dtype=np.complex64))) / np.sqrt(self.vol)

    def fft_mesh(
        self,
        kpoint: int,
        band: int,
        spin: int = 0,
        spinor: int = 0,
        shift: bool = True,
    ) -> np.ndarray:
        """Place the coefficients of a wavefunction onto an fft mesh.

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
        if self.vasp_type is None:
            raise RuntimeError("vasp_type cannot be None.")

        if self.vasp_type.lower()[0] == "n":
            tcoeffs = self.coeffs[kpoint][band][spinor, :]  # type: ignore[call-overload, index]
        elif self.spin == 2:
            tcoeffs = self.coeffs[spin][kpoint][band]  # type: ignore[index]
        else:
            tcoeffs = self.coeffs[kpoint][band]

        mesh = np.zeros(tuple(self.ng), dtype=np.complex128)
        for gp, coeff in zip(self.Gpoints[kpoint], tcoeffs, strict=False):  # type: ignore[call-overload]
            t = tuple(gp.astype(int) + (self.ng / 2).astype(int))
            mesh[t] = coeff

        return np.fft.ifftshift(mesh) if shift else mesh

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
        """Generate a Chgcar object, which is the charge density of the specified
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
            A Chgcar object.
        """
        if phase and not np.all(self.kpoints[kpoint] == 0.0):
            warnings.warn("phase is True should only be used for the Gamma kpoint! I hope you know what you're doing!")

        # Scaling of ng for the fft grid, need to restore value at the end
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

            if self.vasp_type is None:
                raise RuntimeError("vasp_type cannot be None.")

            if phase and (self.vasp_type.lower()[0] != "n" or spinor is not None):
                den = np.sign(np.real(wfr)) * den
            data["total"] = den

        self.ng = temp_ng
        return Chgcar(poscar, data)

    def write_unks(self, directory: PathLike) -> None:
        """Write the UNK files to the given directory.

        Write the cell-periodic part of the Bloch wavefunctions from the
        WAVECAR file to each of the UNK files. There will be one UNK file for
        each of the kpoints in the WAVECAR file.

        Note:
            Wannier90 expects the full kpoint grid instead of the symmetry-
            reduced one that VASP stores the wavefunctions on. You should run
            a nscf calculation with ISYM=0 to obtain the correct grid.

        Args:
            directory (PathLike): directory to write the UNK files.
        """
        out_dir = Path(directory).expanduser()
        if not out_dir.exists():
            out_dir.mkdir(parents=False)
        elif not out_dir.is_dir():
            raise ValueError("invalid directory")

        if self.vasp_type is None:
            raise RuntimeError("vasp_type cannot be None.")

        N = np.prod(self.ng)
        for ik in range(self.nk):
            fname = f"UNK{ik + 1:05d}."
            if self.vasp_type.lower()[0] == "n":
                data = np.empty((self.nb, 2, *self.ng), dtype=np.complex128)
                for ib in range(self.nb):
                    data[ib, 0, :, :, :] = np.fft.ifftn(self.fft_mesh(ik, ib, spinor=0)) * N
                    data[ib, 1, :, :, :] = np.fft.ifftn(self.fft_mesh(ik, ib, spinor=1)) * N
                Unk(ik + 1, data).write_file(str(out_dir / f"{fname}NC"))
            else:
                data = np.empty((self.nb, *self.ng), dtype=np.complex128)
                for ispin in range(self.spin):
                    for ib in range(self.nb):
                        data[ib, :, :, :] = np.fft.ifftn(self.fft_mesh(ik, ib, spin=ispin)) * N
                    Unk(ik + 1, data).write_file(str(out_dir / f"{fname}{ispin + 1}"))


class Eigenval:
    """EIGENVAL file reader.

    Attributes:
        filename (PathLike): The input file.
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

    def __init__(
        self,
        filename: PathLike,
        occu_tol: float = 1e-8,
        separate_spins: bool = False,
    ) -> None:
        """Read input from filename to construct Eigenval object.

        Args:
            filename (PathLike): filename of EIGENVAL to read.
            occu_tol (float): tolerance for determining band gap.
            separate_spins (bool): whether the band gap, CBM, and VBM should be
                reported for each individual spin channel. Defaults to False,
                which computes the eigenvalue band properties independent of
                the spin orientation. If True, the calculation must be spin-polarized.
        """
        self.filename = filename
        self.occu_tol = occu_tol
        self.separate_spins = separate_spins

        with zopen(filename, mode="r") as file:
            self.ispin = int(file.readline().split()[-1])

            # Remove useless header information
            for _ in range(4):
                file.readline()

            self.nelect, self.nkpt, self.nbands = list(map(int, file.readline().split()))

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
            for line in file:
                if re.search(r"(\s+[\-+0-9eE.]+){4}", str(line)):
                    ikpt += 1
                    kpt = list(map(float, line.split()))
                    self.kpoints.append(kpt[:-1])
                    self.kpoints_weights.append(kpt[-1])
                    for i in range(self.nbands):
                        sl = list(map(float, file.readline().split()))
                        if len(sl) == 3:
                            self.eigenvalues[Spin.up][ikpt, i, 0] = sl[1]
                            self.eigenvalues[Spin.up][ikpt, i, 1] = sl[2]
                        elif len(sl) == 5:
                            self.eigenvalues[Spin.up][ikpt, i, 0] = sl[1]
                            self.eigenvalues[Spin.up][ikpt, i, 1] = sl[3]
                            self.eigenvalues[Spin.down][ikpt, i, 0] = sl[2]
                            self.eigenvalues[Spin.down][ikpt, i, 1] = sl[4]

    @property
    def eigenvalue_band_properties(
        self,
    ) -> (
        tuple[float, float, float, bool]
        | tuple[
            tuple[float, float],
            tuple[float, float],
            tuple[float, float],
            tuple[bool, bool],
        ]
    ):
        """Band properties from the eigenvalues as a tuple of
        (band gap, cbm, vbm, is_band_gap_direct).
        In the case of separate_spins=True,
        the band gap, cbm, vbm, and is_band_gap_direct are each tuples of 2,
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
                (
                    max(cbm_spins[0] - vbm_spins[0], 0),
                    max(cbm_spins[1] - vbm_spins[1], 0),
                ),
                (cbm_spins[0], cbm_spins[1]),
                (vbm_spins[0], vbm_spins[1]),
                (
                    vbm_spins_kpoints[0] == cbm_spins_kpoints[0],
                    vbm_spins_kpoints[1] == cbm_spins_kpoints[1],
                ),
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
    def from_formatted(cls, filename: PathLike) -> Self:
        """Read the WAVEDERF file.

        Note: This file is only produced when LOPTICS is true and VASP has been
        recompiled after uncommenting the line that calls
        WRT_CDER_BETWEEN_STATES_FORMATTED in linear_optics.F.

        It is recommended to use `from_binary` instead since the binary file is
        much smaller and contains the same information.

        Args:
            filename (PathLike): The name of the WAVEDER file.

        Returns:
            A Waveder object.
        """
        with zopen(filename, mode="rt") as file:
            nspin, nkpts, nbands = file.readline().split()
        # 1 and 4 are the eigenvalues of the bands (this data is missing in the WAVEDER file)
        # 6:12 are the complex matrix elements in each cartesian direction.
        data = np.loadtxt(filename, skiprows=1, usecols=(1, 4, 6, 7, 8, 9, 10, 11))
        data = data.reshape(int(nspin), int(nkpts), int(nbands), int(nbands), 8)  # slowest to fastest
        cder_real = data[:, :, :, :, 2::2]
        cder_imag = data[:, :, :, :, 3::2]
        # Change to [band1, band2, kpt, spin, cartesian]
        cder_real = np.transpose(cder_real, (2, 3, 1, 0, 4))
        cder_imag = np.transpose(cder_imag, (2, 3, 1, 0, 4))
        # TODO: add eigenvalues?
        return cls(cder_real, cder_imag)

    @classmethod
    def from_binary(
        cls,
        filename: PathLike,
        data_type: Literal["complex64", "float64", "float32"] = "complex64",
    ) -> Self:
        """Read the WAVEDER file.

        Args:
            filename: Name of file containing WAVEDER.
            data_type: Data type of the WAVEDER file. Default is complex64.
                If the file was generated with the "gamma" version of VASP,
                the data type can be either "float64" or "float32".

        Returns:
            Waveder object.
        """

        def read_data(dtype):
            """Read records from Fortran binary file and convert to np.array of given dtype."""
            data = b""
            while True:
                prefix = np.fromfile(file, dtype=np.int32, count=1)[0]
                data += file.read(abs(prefix))
                suffix = np.fromfile(file, dtype=np.int32, count=1)[0]
                if abs(prefix) - abs(suffix):
                    raise RuntimeError(
                        f"Read wrong amount of bytes.\nExpected: {prefix}, read: {len(data)}, suffix: {suffix}."
                    )
                if prefix > 0:
                    break
            return np.frombuffer(data, dtype=dtype)

        with open(filename, "rb") as file:
            nbands, nelect, nk, ispin = read_data(np.int32)
            _ = read_data(np.float64)  # nodes_in_dielectric_function
            _ = read_data(np.float64)  # wplasmon
            me_datatype = np.dtype(data_type)
            cder = read_data(me_datatype)

            cder_data = cder.reshape((3, ispin, nk, nelect, nbands)).T
            return cls(cder_data.real, cder_data.imag)

    @property
    def cder(self) -> np.ndarray:
        """The complex derivative of the orbitals with respect to k."""
        if self.cder_real.shape[0] != self.cder_real.shape[1]:  # pragma: no cover
            warnings.warn(
                "Not all band pairs are present in the WAVEDER file."
                "If you want to get all the matrix elements set LVEL=.True. in the INCAR."
            )
        return self.cder_real + 1j * self.cder_imag

    @property
    def nspin(self) -> int:
        """The number of spin channels."""
        return self.cder_real.shape[3]

    @property
    def nkpoints(self) -> int:
        """The number of k-points."""
        return self.cder_real.shape[2]

    @property
    def nbands(self) -> int:
        """The number of bands."""
        return self.cder_real.shape[0]

    def get_orbital_derivative_between_states(
        self,
        band_i: int,
        band_j: int,
        kpoint: int,
        spin: Literal[0, 1],
        cart_dir: Literal[0, 1, 2],
    ) -> float:
        """
        Get a float between bands band_i and band_j for k-point index,
        spin-channel and Cartesian direction.

        Args:
            band_i (int): Index of band i
            band_j (int): Index of band j
            kpoint (int): Index of k-point
            spin (int): Spin-channel (0 or 1)
            cart_dir (int): Index of Cartesian direction (0, 1, 2)

        Returns:
            a float value
        """
        return self.cder[band_i, band_j, kpoint, spin, cart_dir]


@dataclass
class WSWQ(MSONable):
    r"""Read a WSWQ file.
    The WSWQ file is used to calculation the wave function overlaps between:
        - W: Wavefunctions in the current directory's WAVECAR file.
        - WQ: Wavefunctions stored in the WAVECAR.qqq file.

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
    def data(self) -> np.ndarray:
        """Complex overlap matrix."""
        return self.me_real + 1j * self.me_imag

    @classmethod
    def from_file(cls, filename: str) -> Self:
        """Construct a WSWQ object from a file.

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
        if len(data_res) != nspin * nkpoints * nbands * nbands:
            raise ValueError("incorrect length of data_res")

        data = np.array([complex(real_part, img_part) for (real_part, img_part), _ in data_res])

        # NOTE: loop order (slow->fast) spin -> kpoint -> j -> i
        data = data.reshape((nspin, nkpoints, nbands, nbands))
        data = np.swapaxes(data, 2, 3)  # swap i and j
        return cls(
            nspin=nspin,
            nkpoints=nkpoints,
            nbands=nbands,
            me_real=np.real(data),
            me_imag=np.imag(data),
        )


class UnconvergedVASPWarning(Warning):
    """Warning for unconverged VASP run."""


class VaspDir(collections.abc.Mapping):
    """
    User-friendly class to access all files in a VASP calculation directory as pymatgen objects in a dict.
    Note that the files are lazily parsed to minimize initialization costs since not all files will be needed by all
    users.

    Example:

    ```
    d = VaspDir(".")
    print(d["INCAR"]["NELM"])
    print(d["vasprun.xml"].parameters)
    ```
    """

    FILE_MAPPINGS: typing.ClassVar = {
        n: globals()[n.capitalize()]
        for n in [
            "INCAR",
            "POSCAR",
            "KPOINTS",
            "POTCAR",
            "vasprun",
            "OUTCAR",
            "OSZICAR",
            "CHGCAR",
            "WAVECAR",
            "WAVEDER",
            "LOCPOT",
            "XDATCAR",
            "EIGENVAL",
            "PROCAR",
            "ELFCAR",
            "DYNMAT",
        ]
    } | {
        "CONTCAR": Poscar,
        "IBZKPT": Kpoints,
        "WSWQ": WSWQ,
    }

    def __init__(self, dirname: str | Path):
        """
        Args:
            dirname: The directory containing the VASP calculation as a string or Path.
        """
        self.path = Path(dirname).absolute()
        self.reset()

    def reset(self):
        """
        Reset all loaded files and recheck the directory for files. Use this when the contents of the directory has
        changed.
        """
        # Note that py3.12 has Path.walk(). But we need to use os.walk to ensure backwards compatibility for now.
        self.files = [str(Path(d) / f).lstrip(str(self.path)) for d, _, fnames in os.walk(self.path) for f in fnames]
        self._parsed_files: dict[str, Any] = {}

    def __len__(self):
        return len(self.files)

    def __iter__(self):
        return iter(self.files)

    def __getitem__(self, item):
        if item in self._parsed_files:
            return self._parsed_files[item]
        fpath = self.path / item

        if not (self.path / item).exists():
            raise ValueError(f"{item} not found in {self.path}. List of files are {self.files}.")

        for k, cls_ in VaspDir.FILE_MAPPINGS.items():
            if k in item:
                try:
                    self._parsed_files[item] = cls_.from_file(fpath)
                except AttributeError:
                    self._parsed_files[item] = cls_(fpath)

                return self._parsed_files[item]

        warnings.warn(
            f"No parser defined for {item}. Contents are returned as a string.",
            UserWarning,
        )
        with zopen(fpath, "rt") as f:
            return f.read()

    def get_files_by_name(self, name: str) -> dict[str, Any]:
        """
        Returns all files with a given name. E.g., if you want all the OUTCAR files, set name="OUTCAR".

        Returns:
            {filename: object from VaspDir[filename]}
        """
        return {f: self[f] for f in self.files if name in f}
