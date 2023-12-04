"""
Classes for reading/manipulating/writing VASP input files.
All major VASP input files.
"""

from __future__ import annotations

import codecs
import hashlib
import itertools
import json
import logging
import math
import os
import re
import subprocess
import warnings
from collections import namedtuple
from enum import Enum
from glob import glob
from hashlib import sha256
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
import scipy.constants as const
from monty.dev import deprecated
from monty.io import zopen
from monty.json import MontyDecoder, MSONable
from monty.os import cd
from monty.os.path import zpath
from monty.serialization import dumpfn, loadfn
from tabulate import tabulate

from pymatgen.core import SETTINGS, Element, Lattice, Structure, get_el_sp
from pymatgen.electronic_structure.core import Magmom
from pymatgen.util.io_utils import clean_lines
from pymatgen.util.string import str_delimited

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence

    from numpy.typing import ArrayLike

    from pymatgen.core.trajectory import Vector3D
    from pymatgen.util.typing import PathLike


__author__ = "Shyue Ping Ong, Geoffroy Hautier, Rickard Armiento, Vincent L Chevrier, Stephen Dacek"
__copyright__ = "Copyright 2011, The Materials Project"

logger = logging.getLogger(__name__)
module_dir = os.path.dirname(os.path.abspath(__file__))

# hashes computed from the full POTCAR file contents by pymatgen (not 1st-party VASP hashes)
PYMATGEN_POTCAR_HASHES = loadfn(f"{module_dir}/vasp_potcar_pymatgen_hashes.json")
# written to some newer POTCARs by VASP
VASP_POTCAR_HASHES = loadfn(f"{module_dir}/vasp_potcar_file_hashes.json")
POTCAR_STATS_PATH = os.path.join(module_dir, "potcar-summary-stats.json.bz2")


class Poscar(MSONable):
    """
    Object for representing the data in a POSCAR or CONTCAR file.

    Attributes:
        structure: Associated Structure.
        comment: Optional comment string.
        true_names: Boolean indication whether Poscar contains actual real names parsed
            from either a POTCAR or the POSCAR itself.
        selective_dynamics: Selective dynamics attribute for each site if available.
            A Nx3 array of booleans.
        velocities: Velocities for each site (typically read in from a CONTCAR).
            A Nx3 array of floats.
        predictor_corrector: Predictor corrector coordinates and derivatives for each site;
            i.e. a list of three 1x3 arrays for each site (typically read in from an MD CONTCAR).
        predictor_corrector_preamble: Predictor corrector preamble contains the predictor-corrector key,
            POTIM, and thermostat parameters that precede the site-specific predictor corrector data in MD CONTCAR.
        lattice_velocities: Lattice velocities and current lattice (typically read
            in from an MD CONTCAR). A 6x3 array of floats.
        temperature: Temperature of velocity Maxwell-Boltzmann initialization.
            Initialized to -1 (MB hasn't been performed).
    """

    def __init__(
        self,
        structure: Structure,
        comment: str | None = None,
        selective_dynamics: ArrayLike | None = None,
        true_names: bool = True,
        velocities: ArrayLike | None = None,
        predictor_corrector: ArrayLike | None = None,
        predictor_corrector_preamble: str | None = None,
        lattice_velocities: ArrayLike | None = None,
        sort_structure: bool = False,
    ):
        """
        Args:
            structure (Structure): Structure object.
            comment (str | None, optional): Optional comment line for POSCAR. Defaults to unit
                cell formula of structure. Defaults to None.
            selective_dynamics (ArrayLike | None, optional): Bool values for selective dynamics,
                where N is the number of sites. Defaults to None.
            true_names (bool, optional): Set to False if the names in the POSCAR are not
                well-defined and ambiguous. This situation arises commonly in
                VASP < 5 where the POSCAR sometimes does not contain element
                symbols. Defaults to True.
            velocities (ArrayLike | None, optional): Velocities for the POSCAR. Typically parsed
                in MD runs or can be used to initialize velocities. Defaults to None.
            predictor_corrector (ArrayLike | None, optional): Predictor corrector for the POSCAR.
                Typically parsed in MD runs. Defaults to None.
            predictor_corrector_preamble (str | None, optional): Preamble to the predictor
                corrector. Defaults to None.
            lattice_velocities (ArrayLike | None, optional): Lattice velocities and current
                lattice for the POSCAR. Available in MD runs with variable cell. Defaults to None.
            sort_structure (bool, optional): Whether to sort the structure. Useful if species
                are not grouped properly together. Defaults to False.
        """
        if structure.is_ordered:
            site_properties = {}

            if selective_dynamics is not None:
                selective_dynamics = np.array(selective_dynamics)
                if not selective_dynamics.all():
                    site_properties["selective_dynamics"] = selective_dynamics

            if velocities is not None:
                velocities = np.array(velocities)
                if velocities.any():
                    site_properties["velocities"] = velocities

            if predictor_corrector is not None:
                predictor_corrector = np.array(predictor_corrector)
                if predictor_corrector.any():
                    site_properties["predictor_corrector"] = predictor_corrector

            structure = Structure.from_sites(structure)
            self.structure = structure.copy(site_properties=site_properties)
            if sort_structure:
                self.structure = self.structure.get_sorted_structure()
            self.true_names = true_names
            self.comment = structure.formula if comment is None else comment
            if predictor_corrector_preamble:
                self.structure.properties["predictor_corrector_preamble"] = predictor_corrector_preamble

            if lattice_velocities is not None and np.any(lattice_velocities):
                self.structure.properties["lattice_velocities"] = np.asarray(lattice_velocities)
        else:
            raise ValueError("Disordered structure with partial occupancies cannot be converted into POSCAR!")

        self.temperature = -1.0

    @property
    def velocities(self):
        """Velocities in Poscar."""
        return self.structure.site_properties.get("velocities")

    @property
    def selective_dynamics(self):
        """Selective dynamics in Poscar."""
        return self.structure.site_properties.get("selective_dynamics")

    @property
    def predictor_corrector(self):
        """Predictor corrector in Poscar."""
        return self.structure.site_properties.get("predictor_corrector")

    @property
    def predictor_corrector_preamble(self):
        """Predictor corrector preamble in Poscar."""
        return self.structure.properties.get("predictor_corrector_preamble")

    @property
    def lattice_velocities(self):
        """Lattice velocities in Poscar (including the current lattice vectors)."""
        return self.structure.properties.get("lattice_velocities")

    @velocities.setter  # type: ignore
    def velocities(self, velocities):
        """Setter for Poscar.velocities."""
        self.structure.add_site_property("velocities", velocities)

    @selective_dynamics.setter  # type: ignore
    def selective_dynamics(self, selective_dynamics):
        """Setter for Poscar.selective_dynamics."""
        self.structure.add_site_property("selective_dynamics", selective_dynamics)

    @predictor_corrector.setter  # type: ignore
    def predictor_corrector(self, predictor_corrector):
        """Setter for Poscar.predictor_corrector."""
        self.structure.add_site_property("predictor_corrector", predictor_corrector)

    @predictor_corrector_preamble.setter  # type: ignore
    def predictor_corrector_preamble(self, predictor_corrector_preamble):
        """Setter for Poscar.predictor_corrector."""
        self.structure.properties["predictor_corrector"] = predictor_corrector_preamble

    @lattice_velocities.setter  # type: ignore
    def lattice_velocities(self, lattice_velocities: ArrayLike) -> None:
        """Setter for Poscar.lattice_velocities."""
        self.structure.properties["lattice_velocities"] = np.asarray(lattice_velocities)

    @property
    def site_symbols(self) -> list[str]:
        """
        Sequence of symbols associated with the Poscar. Similar to 6th line in VASP 5+ POSCAR.
        """
        syms = [site.specie.symbol for site in self.structure]
        return [a[0] for a in itertools.groupby(syms)]

    @property
    def natoms(self) -> list[int]:
        """
        Sequence of number of sites of each type associated with the Poscar.
        Similar to 7th line in vasp 5+ POSCAR or the 6th line in vasp 4 POSCAR.
        """
        syms = [site.specie.symbol for site in self.structure]
        return [len(tuple(a[1])) for a in itertools.groupby(syms)]

    def __setattr__(self, name, value):
        if name in ("selective_dynamics", "velocities") and value is not None and len(value) > 0:
            value = np.array(value)
            dim = value.shape
            if dim[1] != 3 or dim[0] != len(self.structure):
                raise ValueError(name + " array must be same length as the structure.")
            value = value.tolist()
        super().__setattr__(name, value)

    @classmethod
    def from_file(cls, filename, check_for_potcar=True, read_velocities=True, **kwargs) -> Poscar:
        """
        Reads a Poscar from a file.

        The code will try its best to determine the elements in the POSCAR in
        the following order:

        1. If check_for_potcar is True, the code will try to check if a POTCAR
        is in the same directory as the POSCAR and use elements from that by
        default. (This is the VASP default sequence of priority).
        2. If the input file is VASP5-like and contains element symbols in the
        6th line, the code will use that if check_for_potcar is False or there
        is no POTCAR found.
        3. Failing (2), the code will check if a symbol is provided at the end
        of each coordinate.

        If all else fails, the code will just assign the first n elements in
        increasing atomic number, where n is the number of species, to the
        Poscar. For example, H, He, Li, .... This will ensure at least a
        unique element is assigned to each site and any analysis that does not
        require specific elemental properties should work fine.

        Args:
            filename (str): File name containing Poscar data.
            check_for_potcar (bool): Whether to check if a POTCAR is present
                in the same directory as the POSCAR. Defaults to True.
            read_velocities (bool): Whether to read or not velocities if they
                are present in the POSCAR. Default is True.

        Returns:
            Poscar object.
        """
        if "check_for_POTCAR" in kwargs:
            warnings.warn("check_for_POTCAR is deprecated. Use check_for_potcar instead.", DeprecationWarning)
            check_for_potcar = kwargs.pop("check_for_POTCAR")

        dirname = os.path.dirname(os.path.abspath(filename))
        names = None
        if check_for_potcar and SETTINGS.get("PMG_POTCAR_CHECKS") is not False:
            potcars = glob(f"{dirname}/*POTCAR*")
            if potcars:
                try:
                    potcar = Potcar.from_file(sorted(potcars)[0])
                    names = [sym.split("_")[0] for sym in potcar.symbols]
                    [get_el_sp(n) for n in names]  # ensure valid names
                except Exception:
                    names = None
        with zopen(filename, "rt") as f:
            return cls.from_str(f.read(), names, read_velocities=read_velocities)

    @classmethod
    @deprecated(message="Use from_str instead")
    def from_string(cls, *args, **kwargs):
        return cls.from_str(*args, **kwargs)

    @classmethod
    def from_str(cls, data, default_names=None, read_velocities=True):
        """
        Reads a Poscar from a string.

        The code will try its best to determine the elements in the POSCAR in
        the following order:

        1. If default_names are supplied and valid, it will use those. Usually,
        default names comes from an external source, such as a POTCAR in the
        same directory.

        2. If there are no valid default names but the input file is VASP5-like
        and contains element symbols in the 6th line, the code will use that.

        3. Failing (2), the code will check if a symbol is provided at the end
        of each coordinate.

        If all else fails, the code will just assign the first n elements in
        increasing atomic number, where n is the number of species, to the
        Poscar. For example, H, He, Li, .... This will ensure at least a
        unique element is assigned to each site and any analysis that does not
        require specific elemental properties should work fine.

        Args:
            data (str): String containing Poscar data.
            default_names ([str]): Default symbols for the POSCAR file,
                usually coming from a POTCAR in the same directory.
            read_velocities (bool): Whether to read or not velocities if they
                are present in the POSCAR. Default is True.

        Returns:
            Poscar object.
        """
        # "^\s*$" doesn't match lines with no whitespace
        chunks = re.split(r"\n\s*\n", data.rstrip(), flags=re.MULTILINE)
        try:
            if chunks[0] == "":
                chunks.pop(0)
                chunks[0] = "\n" + chunks[0]
        except IndexError:
            raise ValueError("Empty POSCAR")

        # Parse positions
        lines = tuple(clean_lines(chunks[0].split("\n"), remove_empty_lines=False))
        comment = lines[0]
        scale = float(lines[1])
        lattice = np.array([[float(i) for i in line.split()] for line in lines[2:5]])
        if scale < 0:
            # In vasp, a negative scale factor is treated as a volume. We need
            # to translate this to a proper lattice vector scaling.
            vol = abs(np.linalg.det(lattice))
            lattice *= (-scale / vol) ** (1 / 3)
        else:
            lattice *= scale

        vasp5_symbols = False
        try:
            n_atoms = [int(i) for i in lines[5].split()]
            ipos = 6
        except ValueError:
            vasp5_symbols = True
            symbols = lines[5].split()

            """
            Atoms and number of atoms in POSCAR written with vasp appear on
            multiple lines when atoms of the same type are not grouped together
            and more than 20 groups are then defined ...
            Example :
            Cr16 Fe35 Ni2
               1.00000000000000
                 8.5415010000000002   -0.0077670000000000   -0.0007960000000000
                -0.0077730000000000    8.5224019999999996    0.0105580000000000
                -0.0007970000000000    0.0105720000000000    8.5356889999999996
               Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Ni   Fe   Cr   Fe   Cr
               Fe   Ni   Fe   Cr   Fe
                 1   1   2   4   2   1   1   1     2     1     1     1     4     1     1     1     5     3     6     1
                 2   1   3   2   5
            Direct
              ...
            """
            n_lines_symbols = 1
            for n_lines_symbols in range(1, 11):
                try:
                    int(lines[5 + n_lines_symbols].split()[0])
                    break
                except ValueError:
                    pass
            for i_line_symbols in range(6, 5 + n_lines_symbols):
                symbols.extend(lines[i_line_symbols].split())
            n_atoms = []
            iline_natoms_start = 5 + n_lines_symbols
            for iline_natoms in range(iline_natoms_start, iline_natoms_start + n_lines_symbols):
                n_atoms.extend([int(i) for i in lines[iline_natoms].split()])
            atomic_symbols = []
            for i, nat in enumerate(n_atoms):
                atomic_symbols.extend([symbols[i]] * nat)
            ipos = 5 + 2 * n_lines_symbols

        pos_type = lines[ipos].split()[0]

        has_selective_dynamics = False
        # Selective dynamics
        if pos_type[0] in "sS":
            has_selective_dynamics = True
            ipos += 1
            pos_type = lines[ipos].split()[0]

        cart = pos_type[0] in "cCkK"
        n_sites = sum(n_atoms)

        # If default_names is specified (usually coming from a POTCAR), use
        # them. This is in line with VASP's parsing order that the POTCAR
        # specified is the default used.
        if default_names:
            try:
                atomic_symbols = []
                for i, nat in enumerate(n_atoms):
                    atomic_symbols.extend([default_names[i]] * nat)
                vasp5_symbols = True
            except IndexError:
                pass

        if not vasp5_symbols:
            ind = 3 if not has_selective_dynamics else 6
            try:
                # Check if names are appended at the end of the coordinates.
                atomic_symbols = [line.split()[ind] for line in lines[ipos + 1 : ipos + 1 + n_sites]]
                # Ensure symbols are valid elements
                if not all(Element.is_valid_symbol(sym) for sym in atomic_symbols):
                    raise ValueError("Non-valid symbols detected.")
                vasp5_symbols = True
            except (ValueError, IndexError):
                # Defaulting to false names.
                atomic_symbols = []
                for i, nat in enumerate(n_atoms):
                    sym = Element.from_Z(i + 1).symbol
                    atomic_symbols.extend([sym] * nat)
                warnings.warn(f"Elements in POSCAR cannot be determined. Defaulting to false names {atomic_symbols}.")
        # read the atomic coordinates
        coords = []
        selective_dynamics = [] if has_selective_dynamics else None
        for i in range(n_sites):
            tokens = lines[ipos + 1 + i].split()
            crd_scale = scale if cart else 1
            coords.append([float(j) * crd_scale for j in tokens[:3]])
            if has_selective_dynamics:
                selective_dynamics.append([tok.upper()[0] == "T" for tok in tokens[3:6]])
        struct = Structure(
            lattice,
            atomic_symbols,
            coords,
            to_unit_cell=False,
            validate_proximity=False,
            coords_are_cartesian=cart,
        )

        if read_velocities:
            # Parse the lattice velocities and current lattice, if present.
            # The header line should contain "Lattice velocities and vectors"
            # There is no space between the coordinates and this section, so
            # it appears in the lines of the first chunk
            lattice_velocities = []
            if len(lines) > ipos + n_sites + 1 and lines[ipos + n_sites + 1].lower().startswith("l"):
                for line in lines[ipos + n_sites + 3 : ipos + n_sites + 9]:
                    lattice_velocities.append([float(tok) for tok in line.split()])

            # Parse velocities if any
            velocities = []
            if len(chunks) > 1:
                for line in chunks[1].strip().split("\n"):
                    velocities.append([float(tok) for tok in line.split()])

            # Parse the predictor-corrector data
            predictor_corrector = []
            predictor_corrector_preamble = None

            if len(chunks) > 2:
                lines = chunks[2].strip().split("\n")
                # There are 3 sets of 3xN Predictor corrector parameters
                # So can't be stored as a single set of "site_property"

                # First line in chunk is a key in CONTCAR
                # Second line is POTIM
                # Third line is the thermostat parameters
                predictor_corrector_preamble = lines[0] + "\n" + lines[1] + "\n" + lines[2]
                # Rest is three sets of parameters, each set contains
                # x, y, z predictor-corrector parameters for every atom in order
                lines = lines[3:]
                for st in range(n_sites):
                    d1 = [float(tok) for tok in lines[st].split()]
                    d2 = [float(tok) for tok in lines[st + n_sites].split()]
                    d3 = [float(tok) for tok in lines[st + 2 * n_sites].split()]
                    predictor_corrector.append([d1, d2, d3])
        else:
            velocities = predictor_corrector = predictor_corrector_preamble = lattice_velocities = None

        return cls(
            struct,
            comment,
            selective_dynamics,
            vasp5_symbols,
            velocities=velocities,
            predictor_corrector=predictor_corrector,
            predictor_corrector_preamble=predictor_corrector_preamble,
            lattice_velocities=lattice_velocities,
        )

    @np.deprecate(message="Use get_str instead")
    def get_string(self, *args, **kwargs) -> str:
        return self.get_str(*args, **kwargs)

    def get_str(self, direct: bool = True, vasp4_compatible: bool = False, significant_figures: int = 16) -> str:
        """
        Returns a string to be written as a POSCAR file. By default, site
        symbols are written, which means compatibility is for vasp >= 5.

        Args:
            direct (bool): Whether coordinates are output in direct or
                Cartesian. Defaults to True.
            vasp4_compatible (bool): Set to True to omit site symbols on 6th
                line to maintain backward vasp 4.x compatibility. Defaults
                to False.
            significant_figures (int): No. of significant figures to
                output all quantities. Defaults to 16. Note that positions are
                output in fixed point, while velocities are output in
                scientific format.

        Returns:
            String representation of POSCAR.
        """
        # This corrects for VASP really annoying bug of crashing on lattices
        # which have triple product < 0. We will just invert the lattice
        # vectors.
        latt = self.structure.lattice
        if np.linalg.det(latt.matrix) < 0:
            latt = Lattice(-latt.matrix)

        format_str = f"{{:{significant_figures + 5}.{significant_figures}f}}"
        lines = [self.comment, "1.0"]
        for vec in latt.matrix:
            lines.append(" ".join(format_str.format(c) for c in vec))

        if self.true_names and not vasp4_compatible:
            lines.append(" ".join(self.site_symbols))
        lines.append(" ".join(map(str, self.natoms)))
        if self.selective_dynamics:
            lines.append("Selective dynamics")
        lines.append("direct" if direct else "cartesian")

        selective_dynamics = self.selective_dynamics
        for i, site in enumerate(self.structure):
            coords = site.frac_coords if direct else site.coords
            line = " ".join(format_str.format(c) for c in coords)
            if selective_dynamics is not None:
                sd = ["T" if j else "F" for j in selective_dynamics[i]]
                line += f" {sd[0]} {sd[1]} {sd[2]}"
            line += " " + site.species_string
            lines.append(line)

        if self.lattice_velocities is not None:
            try:
                lines.append("Lattice velocities and vectors")
                lines.append("  1")
                for velo in self.lattice_velocities:
                    # VASP is strict about the format when reading this quantity
                    lines.append(" ".join(f" {val: .7E}" for val in velo))
            except Exception:
                warnings.warn("Lattice velocities are missing or corrupted.")

        if self.velocities:
            try:
                lines.append("")
                for velo in self.velocities:
                    lines.append(" ".join(format_str.format(val) for val in velo))
            except Exception:
                warnings.warn("Velocities are missing or corrupted.")

        if self.predictor_corrector:
            lines.append("")
            if self.predictor_corrector_preamble:
                lines.append(self.predictor_corrector_preamble)
                pred = np.array(self.predictor_corrector)
                for col in range(3):
                    for z in pred[:, col]:
                        lines.append(" ".join(format_str.format(i) for i in z))
            else:
                warnings.warn(
                    "Preamble information missing or corrupt. Writing Poscar with no predictor corrector data."
                )

        return "\n".join(lines) + "\n"

    def __repr__(self):
        return self.get_str()

    def __str__(self):
        """String representation of Poscar file."""
        return self.get_str()

    def write_file(self, filename: PathLike, **kwargs):
        """
        Writes POSCAR to a file. The supported kwargs are the same as those for
        the Poscar.get_string method and are passed through directly.
        """
        with zopen(filename, "wt") as f:
            f.write(self.get_str(**kwargs))

    def as_dict(self) -> dict:
        """MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "structure": self.structure.as_dict(),
            "true_names": self.true_names,
            "selective_dynamics": np.array(self.selective_dynamics).tolist(),
            "velocities": self.velocities,
            "predictor_corrector": self.predictor_corrector,
            "comment": self.comment,
        }

    @classmethod
    def from_dict(cls, d: dict) -> Poscar:
        """
        :param d: Dict representation.

        Returns:
            Poscar
        """
        return Poscar(
            Structure.from_dict(d["structure"]),
            comment=d["comment"],
            selective_dynamics=d["selective_dynamics"],
            true_names=d["true_names"],
            velocities=d.get("velocities"),
            predictor_corrector=d.get("predictor_corrector"),
        )

    def set_temperature(self, temperature: float):
        """
        Initializes the velocities based on Maxwell-Boltzmann distribution.
        Removes linear, but not angular drift (same as VASP).

        Scales the energies to the exact temperature (microcanonical ensemble)
        Velocities are given in A/fs. This is the vasp default when
        direct/cartesian is not specified (even when positions are given in
        direct coordinates)

        Overwrites imported velocities, if any.

        Args:
            temperature (float): Temperature in Kelvin.
        """
        # mean 0 variance 1
        velocities = np.random.randn(len(self.structure), 3)

        # in AMU, (N,1) array
        atomic_masses = np.array([site.specie.atomic_mass.to("kg") for site in self.structure])
        dof = 3 * len(self.structure) - 3

        # remove linear drift (net momentum)
        velocities -= np.average(atomic_masses[:, np.newaxis] * velocities, axis=0) / np.average(atomic_masses)

        # scale velocities due to atomic masses
        # mean 0 std proportional to sqrt(1/m)
        velocities /= atomic_masses[:, np.newaxis] ** (1 / 2)

        # scale velocities to get correct temperature
        energy = np.sum(1 / 2 * atomic_masses * np.sum(velocities**2, axis=1))
        scale = (temperature * dof / (2 * energy / const.k)) ** (1 / 2)

        velocities *= scale * 1e-5  # these are in A/fs

        self.temperature = temperature
        self.structure.site_properties.pop("selective_dynamics", None)
        self.structure.site_properties.pop("predictor_corrector", None)
        # returns as a list of lists to be consistent with the other
        # initializations

        self.structure.add_site_property("velocities", velocities.tolist())


with open(f"{module_dir}/incar_parameters.json") as json_file:
    incar_params = json.loads(json_file.read())


class BadIncarWarning(UserWarning):
    """Warning class for bad Incar parameters."""


class Incar(dict, MSONable):
    """
    INCAR object for reading and writing INCAR files. Essentially consists of
    a dictionary with some helper functions.
    """

    def __init__(self, params: dict[str, Any] | None = None):
        """
        Creates an Incar object.

        Args:
            params (dict): A set of input parameters as a dictionary.
        """
        super().__init__()
        if params:
            # if Incar contains vector-like magmoms given as a list
            # of floats, convert to a list of lists
            if (params.get("MAGMOM") and isinstance(params["MAGMOM"][0], (int, float))) and (
                params.get("LSORBIT") or params.get("LNONCOLLINEAR")
            ):
                val = []
                for i in range(len(params["MAGMOM"]) // 3):
                    val.append(params["MAGMOM"][i * 3 : (i + 1) * 3])
                params["MAGMOM"] = val

            self.update(params)

    def __setitem__(self, key: str, val: Any):
        """
        Add parameter-val pair to Incar. Warns if parameter is not in list of
        valid INCAR tags. Also cleans the parameter and val by stripping
        leading and trailing white spaces.
        """
        super().__setitem__(
            key.strip(),
            Incar.proc_val(key.strip(), val.strip()) if isinstance(val, str) else val,
        )

    def as_dict(self) -> dict:
        """MSONable dict."""
        dct = dict(self)
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        return dct

    @classmethod
    def from_dict(cls, d) -> Incar:
        """
        :param d: Dict representation.

        Returns:
            Incar
        """
        if d.get("MAGMOM") and isinstance(d["MAGMOM"][0], dict):
            d["MAGMOM"] = [Magmom.from_dict(m) for m in d["MAGMOM"]]
        return Incar({k: v for k, v in d.items() if k not in ("@module", "@class")})

    @np.deprecate(message="Use get_str instead")
    def get_string(self, *args, **kwargs) -> str:
        return self.get_str(*args, **kwargs)

    def get_str(self, sort_keys: bool = False, pretty: bool = False) -> str:
        """
        Returns a string representation of the INCAR. The reason why this
        method is different from the __str__ method is to provide options for
        pretty printing.

        Args:
            sort_keys (bool): Set to True to sort the INCAR parameters
                alphabetically. Defaults to False.
            pretty (bool): Set to True for pretty aligned output. Defaults
                to False.
        """
        keys = list(self)
        if sort_keys:
            keys = sorted(keys)
        lines = []
        for k in keys:
            if k == "MAGMOM" and isinstance(self[k], list):
                value = []

                if isinstance(self[k][0], (list, Magmom)) and (self.get("LSORBIT") or self.get("LNONCOLLINEAR")):
                    value.append(" ".join(str(i) for j in self[k] for i in j))
                elif self.get("LSORBIT") or self.get("LNONCOLLINEAR"):
                    for m, g in itertools.groupby(self[k]):
                        value.append(f"3*{len(tuple(g))}*{m}")
                else:
                    # float() to ensure backwards compatibility between
                    # float magmoms and Magmom objects
                    for m, g in itertools.groupby(self[k], key=float):
                        value.append(f"{len(tuple(g))}*{m}")

                lines.append([k, " ".join(value)])
            elif isinstance(self[k], list):
                lines.append([k, " ".join(map(str, self[k]))])
            else:
                lines.append([k, self[k]])

        if pretty:
            return str(tabulate([[line[0], "=", line[1]] for line in lines], tablefmt="plain"))
        return str_delimited(lines, None, " = ") + "\n"

    def __str__(self):
        return self.get_str(sort_keys=True, pretty=False)

    def write_file(self, filename: PathLike):
        """Write Incar to a file.

        Args:
            filename (str): filename to write to.
        """
        with zopen(filename, "wt") as f:
            f.write(str(self))

    @classmethod
    def from_file(cls, filename: PathLike) -> Incar:
        """Reads an Incar object from a file.

        Args:
            filename (str): Filename for file

        Returns:
            Incar object
        """
        with zopen(filename, "rt") as f:
            return cls.from_str(f.read())

    @classmethod
    @np.deprecate(message="Use from_str instead")
    def from_string(cls, *args, **kwargs):
        return cls.from_str(*args, **kwargs)

    @classmethod
    def from_str(cls, string: str) -> Incar:
        """Reads an Incar object from a string.

        Args:
            string (str): Incar string

        Returns:
            Incar object
        """
        lines = list(clean_lines(string.splitlines()))
        params = {}
        for line in lines:
            for sline in line.split(";"):
                m = re.match(r"(\w+)\s*=\s*(.*)", sline.strip())
                if m:
                    key = m.group(1).strip()
                    val = m.group(2).strip()
                    val = Incar.proc_val(key, val)
                    params[key] = val
        return cls(params)

    @staticmethod
    def proc_val(key: str, val: Any):
        """Helper method to convert INCAR parameters to proper types like ints, floats, lists, etc.

        Args:
            key: INCAR parameter key
            val: Actual value of INCAR parameter.
        """
        list_keys = ("LDAUU", "LDAUL", "LDAUJ", "MAGMOM", "DIPOL", "LANGEVIN_GAMMA", "QUAD_EFG", "EINT")
        bool_keys = (
            "LDAU",
            "LWAVE",
            "LSCALU",
            "LCHARG",
            "LPLANE",
            "LUSE_VDW",
            "LHFCALC",
            "ADDGRID",
            "LSORBIT",
            "LNONCOLLINEAR",
        )
        float_keys = ("EDIFF", "SIGMA", "TIME", "ENCUTFOCK", "HFSCREEN", "POTIM", "EDIFFG", "AGGAC", "PARAM1", "PARAM2")
        int_keys = (
            "NSW",
            "NBANDS",
            "NELMIN",
            "ISIF",
            "IBRION",
            "ISPIN",
            "ISTART",
            "ICHARG",
            "NELM",
            "ISMEAR",
            "NPAR",
            "LDAUPRINT",
            "LMAXMIX",
            "ENCUT",
            "NSIM",
            "NKRED",
            "NUPDOWN",
            "ISPIND",
            "LDAUTYPE",
            "IVDW",
        )

        def smart_int_or_float(num_str):
            if num_str.find(".") != -1 or num_str.lower().find("e") != -1:
                return float(num_str)
            return int(num_str)

        try:
            if key in list_keys:
                output = []
                tokens = re.findall(r"(-?\d+\.?\d*)\*?(-?\d+\.?\d*)?\*?(-?\d+\.?\d*)?", val)
                for tok in tokens:
                    if tok[2] and "3" in tok[0]:
                        output.extend([smart_int_or_float(tok[2])] * int(tok[0]) * int(tok[1]))
                    elif tok[1]:
                        output.extend([smart_int_or_float(tok[1])] * int(tok[0]))
                    else:
                        output.append(smart_int_or_float(tok[0]))
                return output
            if key in bool_keys:
                m = re.match(r"^\.?([T|F|t|f])[A-Za-z]*\.?", val)
                if m:
                    return m.group(1).lower() == "t"

                raise ValueError(key + " should be a boolean type!")

            if key in float_keys:
                return float(re.search(r"^-?\d*\.?\d*[e|E]?-?\d*", val).group(0))  # type: ignore

            if key in int_keys:
                return int(re.match(r"^-?[0-9]+", val).group(0))  # type: ignore

        except ValueError:
            pass

        # Not in standard keys. We will try a hierarchy of conversions.
        try:
            return int(val)
        except ValueError:
            pass

        try:
            return float(val)
        except ValueError:
            pass

        if "true" in val.lower():
            return True

        if "false" in val.lower():
            return False

        return val.strip().capitalize()

    def diff(self, other: Incar) -> dict[str, dict[str, Any]]:
        """
        Diff function for Incar. Compares two Incars and indicates which
        parameters are the same and which are not. Useful for checking whether
        two runs were done using the same parameters.

        Args:
            other (Incar): The other Incar object to compare to.

        Returns:
            dict[str, dict]: of the following format:
                {"Same" : parameters_that_are_the_same, "Different": parameters_that_are_different}
                Note that the parameters are return as full dictionaries of values. E.g. {"ISIF":3}
        """
        similar_param = {}
        different_param = {}
        for k1, v1 in self.items():
            if k1 not in other:
                different_param[k1] = {"INCAR1": v1, "INCAR2": None}
            elif v1 != other[k1]:
                different_param[k1] = {"INCAR1": v1, "INCAR2": other[k1]}
            else:
                similar_param[k1] = v1
        for k2, v2 in other.items():
            if k2 not in similar_param and k2 not in different_param and k2 not in self:
                different_param[k2] = {"INCAR1": None, "INCAR2": v2}
        return {"Same": similar_param, "Different": different_param}

    def __add__(self, other):
        """
        Add all the values of another INCAR object to this object.
        Facilitates the use of "standard" INCARs.
        """
        params = dict(self.items())
        for key, val in other.items():
            if key in self and val != self[key]:
                raise ValueError(f"Incars have conflicting values for {key}: {self[key]} != {val}")
            params[key] = val
        return Incar(params)

    def check_params(self):
        """
        Raises a warning for nonsensical or non-existent INCAR tags and
        parameters. If a keyword doesn't exist (e.g. there's a typo in a
        keyword), your calculation will still run, however VASP will ignore the
        parameter without letting you know, hence why we have this Incar method.
        """
        for key, val in self.items():
            # First check if this parameter even exists
            if key not in incar_params:
                warnings.warn(f"Cannot find {key} in the list of INCAR flags", BadIncarWarning, stacklevel=2)

            if key in incar_params:
                if type(incar_params[key]).__name__ == "str":
                    # Now we check if this is an appropriate parameter type
                    if incar_params[key] == "float":
                        if not type(val) not in ["float", "int"]:
                            warnings.warn(f"{key}: {val} is not real", BadIncarWarning, stacklevel=2)
                    elif type(val).__name__ != incar_params[key]:
                        warnings.warn(f"{key}: {val} is not a {incar_params[key]}", BadIncarWarning, stacklevel=2)

                # if we have a list of possible parameters, check
                # if the user given parameter is in this list
                elif type(incar_params[key]).__name__ == "list" and val not in incar_params[key]:
                    warnings.warn(f"{key}: Cannot find {val} in the list of parameters", BadIncarWarning, stacklevel=2)


class KpointsSupportedModes(Enum):
    """Enum type of all supported modes for Kpoint generation."""

    Automatic = 0
    Gamma = 1
    Monkhorst = 2
    Line_mode = 3
    Cartesian = 4
    Reciprocal = 5

    def __str__(self):
        return str(self.name)

    @classmethod
    @np.deprecate(message="Use from_str instead")
    def from_string(cls, *args, **kwargs):
        return cls.from_str(*args, **kwargs)

    @classmethod
    def from_str(cls, mode: str) -> KpointsSupportedModes:
        """
        :param s: String

        Returns:
            Kpoints_supported_modes
        """
        initial = mode.lower()[0]
        for key in KpointsSupportedModes:
            if key.name.lower()[0] == initial:
                return key
        raise ValueError(f"Invalid Kpoint {mode=}")


class Kpoints(MSONable):
    """KPOINT reader/writer."""

    supported_modes = KpointsSupportedModes

    def __init__(
        self,
        comment: str = "Default gamma",
        num_kpts: int = 0,
        style: KpointsSupportedModes = supported_modes.Gamma,
        kpts: Sequence[float | Sequence] = ((1, 1, 1),),
        kpts_shift: Vector3D = (0, 0, 0),
        kpts_weights=None,
        coord_type=None,
        labels=None,
        tet_number: int = 0,
        tet_weight: float = 0,
        tet_connections=None,
    ):
        """
        Highly flexible constructor for Kpoints object. The flexibility comes
        at the cost of usability and in general, it is recommended that you use
        the default constructor only if you know exactly what you are doing and
        requires the flexibility. For most usage cases, the three automatic
        schemes can be constructed far more easily using the convenience static
        constructors (automatic, gamma_automatic, monkhorst_automatic) and it
        is recommended that you use those.

        Args:
            comment (str): String comment for Kpoints. Defaults to "Default gamma".
            num_kpts: Following VASP method of defining the KPOINTS file, this
                parameter is the number of kpoints specified. If set to 0
                (or negative), VASP automatically generates the KPOINTS.
            style: Style for generating KPOINTS. Use one of the
                Kpoints.supported_modes enum types.
            kpts (2D array): 2D array of kpoints. Even when only a single
                specification is required, e.g. in the automatic scheme,
                the kpts should still be specified as a 2D array. e.g.,
                [[20]] or [[2,2,2]].
            kpts_shift (3x1 array): Shift for Kpoints.
            kpts_weights: Optional weights for kpoints. Weights should be
                integers. For explicit kpoints.
            coord_type: In line-mode, this variable specifies whether the
                Kpoints were given in Cartesian or Reciprocal coordinates.
            labels: In line-mode, this should provide a list of labels for
                each kpt. It is optional in explicit kpoint mode as comments for
                k-points.
            tet_number: For explicit kpoints, specifies the number of
                tetrahedrons for the tetrahedron method.
            tet_weight: For explicit kpoints, specifies the weight for each
                tetrahedron for the tetrahedron method.
            tet_connections: For explicit kpoints, specifies the connections
                of the tetrahedrons for the tetrahedron method.
                Format is a list of tuples, [ (sym_weight, [tet_vertices]),
                ...]

        The default behavior of the constructor is for a Gamma centered,
        1x1x1 KPOINTS with no shift.
        """
        if num_kpts > 0 and not labels and not kpts_weights:
            raise ValueError("For explicit or line-mode kpoints, either the labels or kpts_weights must be specified.")

        self.comment = comment
        self.num_kpts = num_kpts
        self.kpts = kpts
        self.style = style
        self.coord_type = coord_type
        self.kpts_weights = kpts_weights
        self.kpts_shift = kpts_shift
        self.labels = labels
        self.tet_number = tet_number
        self.tet_weight = tet_weight
        self.tet_connections = tet_connections

    @property
    def style(self) -> KpointsSupportedModes:
        """
        Style for kpoint generation. One of Kpoints_supported_modes enum.
        """
        return self._style

    @style.setter
    def style(self, style):
        """
        :param style: Style

        Returns:
            Sets the style for the Kpoints. One of Kpoints_supported_modes
            enum.
        """
        if isinstance(style, str):
            style = Kpoints.supported_modes.from_str(style)

        if (
            style
            in (Kpoints.supported_modes.Automatic, Kpoints.supported_modes.Gamma, Kpoints.supported_modes.Monkhorst)
            and len(self.kpts) > 1
        ):
            raise ValueError(
                "For fully automatic or automatic gamma or monk "
                "kpoints, only a single line for the number of "
                "divisions is allowed."
            )

        self._style = style

    @staticmethod
    def automatic(subdivisions):
        """
        Convenient static constructor for a fully automatic Kpoint grid, with
        gamma centered Monkhorst-Pack grids and the number of subdivisions
        along each reciprocal lattice vector determined by the scheme in the
        VASP manual.

        Args:
            subdivisions: Parameter determining number of subdivisions along
                each reciprocal lattice vector.

        Returns:
            Kpoints object
        """
        return Kpoints(
            "Fully automatic kpoint scheme", 0, style=Kpoints.supported_modes.Automatic, kpts=[[subdivisions]]
        )

    @staticmethod
    def gamma_automatic(kpts: tuple[int, int, int] = (1, 1, 1), shift: Vector3D = (0, 0, 0)):
        """
        Convenient static constructor for an automatic Gamma centered Kpoint
        grid.

        Args:
            kpts: Subdivisions N_1, N_2 and N_3 along reciprocal lattice
                vectors. Defaults to (1,1,1)
            shift: Shift to be applied to the kpoints. Defaults to (0,0,0).

        Returns:
            Kpoints object
        """
        return Kpoints("Automatic kpoint scheme", 0, Kpoints.supported_modes.Gamma, kpts=[kpts], kpts_shift=shift)

    @staticmethod
    def monkhorst_automatic(kpts: tuple[int, int, int] = (2, 2, 2), shift: Vector3D = (0, 0, 0)):
        """
        Convenient static constructor for an automatic Monkhorst pack Kpoint
        grid.

        Args:
            kpts: Subdivisions N_1, N_2, N_3 along reciprocal lattice
                vectors. Defaults to (2,2,2)
            shift: Shift to be applied to the kpoints. Defaults to (0,0,0).

        Returns:
            Kpoints object
        """
        return Kpoints("Automatic kpoint scheme", 0, Kpoints.supported_modes.Monkhorst, kpts=[kpts], kpts_shift=shift)

    @staticmethod
    def automatic_density(structure: Structure, kppa: float, force_gamma: bool = False):
        """
        Returns an automatic Kpoint object based on a structure and a kpoint
        density. Uses Gamma centered meshes for hexagonal cells and face-centered cells,
        Monkhorst-Pack grids otherwise.

        Algorithm:
            Uses a simple approach scaling the number of divisions along each
            reciprocal lattice vector proportional to its length.

        Args:
            structure (Structure): Input structure
            kppa (float): Grid density
            force_gamma (bool): Force a gamma centered mesh (default is to
                use gamma only for hexagonal cells or odd meshes)

        Returns:
            Kpoints
        """
        comment = f"pymatgen with grid density = {kppa:.0f} / number of atoms"
        if math.fabs((math.floor(kppa ** (1 / 3) + 0.5)) ** 3 - kppa) < 1:
            kppa += kppa * 0.01
        lattice = structure.lattice
        lengths = lattice.abc
        ngrid = kppa / len(structure)
        mult = (ngrid * lengths[0] * lengths[1] * lengths[2]) ** (1 / 3)

        num_div = [int(math.floor(max(mult / length, 1))) for length in lengths]

        is_hexagonal = lattice.is_hexagonal()
        is_face_centered = structure.get_space_group_info()[0][0] == "F"
        has_odd = any(idx % 2 == 1 for idx in num_div)
        if has_odd or is_hexagonal or is_face_centered or force_gamma:
            style = Kpoints.supported_modes.Gamma
        else:
            style = Kpoints.supported_modes.Monkhorst

        return Kpoints(comment, 0, style, [num_div], (0, 0, 0))

    @staticmethod
    def automatic_gamma_density(structure: Structure, kppa: float):
        """
        Returns an automatic Kpoint object based on a structure and a kpoint
        density. Uses Gamma centered meshes always. For GW.

        Algorithm:
            Uses a simple approach scaling the number of divisions along each
            reciprocal lattice vector proportional to its length.

        Args:
            structure: Input structure
            kppa: Grid density
        """
        latt = structure.lattice
        a, b, c = latt.abc
        ngrid = kppa / len(structure)

        mult = (ngrid * a * b * c) ** (1 / 3)
        num_div = [int(round(mult / length)) for length in latt.abc]

        # ensure that all num_div[i] > 0
        num_div = [idx if idx > 0 else 1 for idx in num_div]

        # VASP documentation recommends to use even grids for n <= 8 and odd grids for n > 8.
        num_div = [idx + idx % 2 if idx <= 8 else idx - idx % 2 + 1 for idx in num_div]

        style = Kpoints.supported_modes.Gamma

        comment = f"pymatgen with grid density = {kppa:.0f} / number of atoms"

        num_kpts = 0
        return Kpoints(comment, num_kpts, style, [num_div], (0, 0, 0))

    @staticmethod
    def automatic_density_by_vol(structure: Structure, kppvol: int, force_gamma: bool = False) -> Kpoints:
        """
        Returns an automatic Kpoint object based on a structure and a kpoint
        density per inverse Angstrom^3 of reciprocal cell.

        Algorithm:
            Same as automatic_density()

        Args:
            structure (Structure): Input structure
            kppvol (int): Grid density per Angstrom^(-3) of reciprocal cell
            force_gamma (bool): Force a gamma centered mesh

        Returns:
            Kpoints
        """
        vol = structure.lattice.reciprocal_lattice.volume
        kppa = kppvol * vol * len(structure)
        return Kpoints.automatic_density(structure, kppa, force_gamma=force_gamma)

    @staticmethod
    def automatic_density_by_lengths(
        structure: Structure, length_densities: Sequence[float], force_gamma: bool = False
    ):
        """
        Returns an automatic Kpoint object based on a structure and a k-point
        density normalized by lattice constants.

        Algorithm:
            For a given dimension, the # of k-points is chosen as
            length_density = # of kpoints * lattice constant, e.g. [50.0, 50.0, 1.0] would
            have k-points of 50/a x 50/b x 1/c.

        Args:
            structure (Structure): Input structure
            length_densities (list[floats]): Defines the density of k-points in each
            dimension, e.g. [50.0, 50.0, 1.0].
            force_gamma (bool): Force a gamma centered mesh

        Returns:
            Kpoints
        """
        if len(length_densities) != 3:
            msg = f"The dimensions of length_densities must be 3, not {len(length_densities)}"
            raise ValueError(msg)
        comment = f"k-point density of {length_densities}/[a, b, c]"
        lattice = structure.lattice
        abc = lattice.abc
        num_div = [np.ceil(ld / abc[idx]) for idx, ld in enumerate(length_densities)]
        is_hexagonal = lattice.is_hexagonal()
        is_face_centered = structure.get_space_group_info()[0][0] == "F"
        has_odd = any(idx % 2 == 1 for idx in num_div)
        if has_odd or is_hexagonal or is_face_centered or force_gamma:
            style = Kpoints.supported_modes.Gamma
        else:
            style = Kpoints.supported_modes.Monkhorst

        return Kpoints(comment, 0, style, [num_div], (0, 0, 0))

    @staticmethod
    def automatic_linemode(divisions, ibz):
        """
        Convenient static constructor for a KPOINTS in mode line_mode.
        gamma centered Monkhorst-Pack grids and the number of subdivisions
        along each reciprocal lattice vector determined by the scheme in the
        VASP manual.

        Args:
            divisions: Parameter determining the number of k-points along each high symmetry line.
            ibz: HighSymmKpath object (pymatgen.symmetry.bandstructure)

        Returns:
            Kpoints object
        """
        kpoints = []
        labels = []
        for path in ibz.kpath["path"]:
            kpoints.append(ibz.kpath["kpoints"][path[0]])
            labels.append(path[0])
            for i in range(1, len(path) - 1):
                kpoints.append(ibz.kpath["kpoints"][path[i]])
                labels.append(path[i])
                kpoints.append(ibz.kpath["kpoints"][path[i]])
                labels.append(path[i])

            kpoints.append(ibz.kpath["kpoints"][path[-1]])
            labels.append(path[-1])

        return Kpoints(
            "Line_mode KPOINTS file",
            style=Kpoints.supported_modes.Line_mode,
            coord_type="Reciprocal",
            kpts=kpoints,
            labels=labels,
            num_kpts=int(divisions),
        )

    @classmethod
    def from_file(cls, filename):
        """
        Reads a Kpoints object from a KPOINTS file.

        Args:
            filename (str): filename to read from.

        Returns:
            Kpoints object
        """
        with zopen(filename, "rt") as f:
            return cls.from_str(f.read())

    @classmethod
    @np.deprecate(message="Use from_str instead")
    def from_string(cls, *args, **kwargs):
        return cls.from_str(*args, **kwargs)

    @classmethod
    def from_str(cls, string):
        """
        Reads a Kpoints object from a KPOINTS string.

        Args:
            string (str): KPOINTS string.

        Returns:
            Kpoints object
        """
        lines = [line.strip() for line in string.splitlines()]

        comment = lines[0]
        num_kpts = int(lines[1].split()[0].strip())
        style = lines[2].lower()[0]

        # Fully automatic KPOINTS
        if style == "a":
            return cls.automatic(int(lines[3].split()[0].strip()))

        coord_pattern = re.compile(r"^\s*([\d+.\-Ee]+)\s+([\d+.\-Ee]+)\s+([\d+.\-Ee]+)")

        # Automatic gamma and Monk KPOINTS, with optional shift
        if style in ["g", "m"]:
            kpts = [int(i) for i in lines[3].split()]
            kpts_shift = (0, 0, 0)
            if len(lines) > 4 and coord_pattern.match(lines[4]):
                try:
                    kpts_shift = [float(i) for i in lines[4].split()]
                except ValueError:
                    pass
            return cls.gamma_automatic(kpts, kpts_shift) if style == "g" else cls.monkhorst_automatic(kpts, kpts_shift)

        # Automatic kpoints with basis
        if num_kpts <= 0:
            style = cls.supported_modes.Cartesian if style in "ck" else cls.supported_modes.Reciprocal
            kpts = [[float(j) for j in lines[i].split()] for i in range(3, 6)]
            kpts_shift = [float(i) for i in lines[6].split()]
            return Kpoints(
                comment=comment,
                num_kpts=num_kpts,
                style=style,
                kpts=kpts,
                kpts_shift=kpts_shift,
            )

        # Line-mode KPOINTS, usually used with band structures
        if style == "l":
            coord_type = "Cartesian" if lines[3].lower()[0] in "ck" else "Reciprocal"
            style = cls.supported_modes.Line_mode
            kpts = []
            labels = []
            patt = re.compile(r"([e0-9.\-]+)\s+([e0-9.\-]+)\s+([e0-9.\-]+)\s*!*\s*(.*)")
            for i in range(4, len(lines)):
                line = lines[i]
                m = patt.match(line)
                if m:
                    kpts.append([float(m.group(1)), float(m.group(2)), float(m.group(3))])
                    labels.append(m.group(4).strip())
            return Kpoints(
                comment=comment,
                num_kpts=num_kpts,
                style=style,
                kpts=kpts,
                coord_type=coord_type,
                labels=labels,
            )

        # Assume explicit KPOINTS if all else fails.
        style = cls.supported_modes.Cartesian if style in "ck" else cls.supported_modes.Reciprocal
        kpts = []
        kpts_weights = []
        labels = []
        tet_number = 0
        tet_weight = 0
        tet_connections = None

        for i in range(3, 3 + num_kpts):
            tokens = lines[i].split()
            kpts.append([float(j) for j in tokens[0:3]])
            kpts_weights.append(float(tokens[3]))
            if len(tokens) > 4:
                labels.append(tokens[4])
            else:
                labels.append(None)
        try:
            # Deal with tetrahedron method
            if lines[3 + num_kpts].strip().lower()[0] == "t":
                tokens = lines[4 + num_kpts].split()
                tet_number = int(tokens[0])
                tet_weight = float(tokens[1])
                tet_connections = []
                for i in range(5 + num_kpts, 5 + num_kpts + tet_number):
                    tokens = lines[i].split()
                    tet_connections.append((int(tokens[0]), [int(tokens[j]) for j in range(1, 5)]))
        except IndexError:
            pass

        return cls(
            comment=comment,
            num_kpts=num_kpts,
            style=cls.supported_modes[str(style)],
            kpts=kpts,
            kpts_weights=kpts_weights,
            tet_number=tet_number,
            tet_weight=tet_weight,
            tet_connections=tet_connections,
            labels=labels,
        )

    def write_file(self, filename):
        """
        Write Kpoints to a file.

        Args:
            filename (str): Filename to write to.
        """
        with zopen(filename, "wt") as f:
            f.write(str(self))

    def __repr__(self):
        lines = [self.comment, str(self.num_kpts), self.style.name]
        style = self.style.name.lower()[0]
        if style == "l":
            lines.append(self.coord_type)
        for idx, kpt in enumerate(self.kpts):
            lines.append(" ".join(map(str, kpt)))
            if style == "l":
                lines[-1] += " ! " + self.labels[idx]
                if idx % 2 == 1:
                    lines[-1] += "\n"
            elif self.num_kpts > 0:
                if self.labels is not None:
                    lines[-1] += f" {int(self.kpts_weights[idx])} {self.labels[idx]}"
                else:
                    lines[-1] += f" {int(self.kpts_weights[idx])}"

        # Print tetrahedron parameters if the number of tetrahedrons > 0
        if style not in "lagm" and self.tet_number > 0:
            lines.extend(("Tetrahedron", f"{self.tet_number} {self.tet_weight:f}"))
            for sym_weight, vertices in self.tet_connections:
                a, b, c, d = vertices
                lines.append(f"{sym_weight} {a} {b} {c} {d}")

        # Print shifts for automatic kpoints types if not zero.
        if self.num_kpts <= 0 and tuple(self.kpts_shift) != (0, 0, 0):
            lines.append(" ".join(map(str, self.kpts_shift)))
        return "\n".join(lines) + "\n"

    def as_dict(self):
        """MSONable dict."""
        dct = {
            "comment": self.comment,
            "nkpoints": self.num_kpts,
            "generation_style": self.style.name,
            "kpoints": self.kpts,
            "usershift": self.kpts_shift,
            "kpts_weights": self.kpts_weights,
            "coord_type": self.coord_type,
            "labels": self.labels,
            "tet_number": self.tet_number,
            "tet_weight": self.tet_weight,
            "tet_connections": self.tet_connections,
        }

        optional_paras = ["genvec1", "genvec2", "genvec3", "shift"]
        for para in optional_paras:
            if para in self.__dict__:
                dct[para] = self.__dict__[para]
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        return dct

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation.

        Returns:
            Kpoints
        """
        comment = d.get("comment", "")
        generation_style = d.get("generation_style")
        kpts = d.get("kpoints", [[1, 1, 1]])
        kpts_shift = d.get("usershift", [0, 0, 0])
        num_kpts = d.get("nkpoints", 0)
        return cls(
            comment=comment,
            kpts=kpts,
            style=generation_style,
            kpts_shift=kpts_shift,
            num_kpts=num_kpts,
            kpts_weights=d.get("kpts_weights"),
            coord_type=d.get("coord_type"),
            labels=d.get("labels"),
            tet_number=d.get("tet_number", 0),
            tet_weight=d.get("tet_weight", 0),
            tet_connections=d.get("tet_connections"),
        )


def _parse_bool(s):
    m = re.match(r"^\.?([TFtf])[A-Za-z]*\.?", s)
    if m:
        return m.group(1) in ["T", "t"]
    raise ValueError(s + " should be a boolean type!")


def _parse_float(s):
    return float(re.search(r"^-?\d*\.?\d*[eE]?-?\d*", s).group(0))


def _parse_int(s):
    return int(re.match(r"^-?[0-9]+", s).group(0))


def _parse_list(s):
    return [float(y) for y in re.split(r"\s+", s.strip()) if not y.isalpha()]


Orbital = namedtuple("Orbital", ["n", "l", "j", "E", "occ"])
OrbitalDescription = namedtuple("OrbitalDescription", ["l", "E", "Type", "Rcut", "Type2", "Rcut2"])


class UnknownPotcarWarning(UserWarning):
    """Warning raised when POTCAR hashes do not pass validation."""


class PotcarSingle:
    """
    Object for a **single** POTCAR. The builder assumes the POTCAR contains
    the complete untouched data in "data" as a string and a dict of keywords.

    Attributes:
        data (str): POTCAR data as a string.
        keywords (dict): Keywords parsed from the POTCAR as a dict. All keywords are also
            accessible as attributes in themselves. E.g., potcar.enmax, potcar.encut, etc.

    md5 hashes of the entire POTCAR file and the actual data are validated
    against a database of known good hashes. Appropriate warnings or errors
    are raised if a POTCAR hash fails validation.
    """

    """
    NB: there are multiple releases of the {LDA,PBE} {52,54} POTCARs
        the original (univie) releases include no SHA256 hashes nor COPYR fields
        in the PSCTR/header field.
    We indicate the older release in `functional_dir` as PBE_52, PBE_54, LDA_52, LDA_54.
    The newer release is indicated as PBE_52_W_HASH, etc.
    """
    functional_dir = dict(
        PBE="POT_GGA_PAW_PBE",
        PBE_52="POT_GGA_PAW_PBE_52",
        PBE_52_W_HASH="POTPAW_PBE_52",
        PBE_54="POT_GGA_PAW_PBE_54",
        PBE_54_W_HASH="POTPAW_PBE_54",
        PBE_64="POT_PAW_PBE_64",
        LDA="POT_LDA_PAW",
        LDA_52="POT_LDA_PAW_52",
        LDA_52_W_HASH="POTPAW_LDA_52",
        LDA_54="POT_LDA_PAW_54",
        LDA_54_W_HASH="POTPAW_LDA_54",
        LDA_64="POT_LDA_PAW_64",
        PW91="POT_GGA_PAW_PW91",
        LDA_US="POT_LDA_US",
        PW91_US="POT_GGA_US_PW91",
        Perdew_Zunger81="POT_LDA_PAW",
    )

    functional_tags = {
        "pe": {"name": "PBE", "class": "GGA"},
        "91": {"name": "PW91", "class": "GGA"},
        "rp": {"name": "revPBE", "class": "GGA"},
        "am": {"name": "AM05", "class": "GGA"},
        "ps": {"name": "PBEsol", "class": "GGA"},
        "pw": {"name": "PW86", "class": "GGA"},
        "lm": {"name": "Langreth-Mehl-Hu", "class": "GGA"},
        "pb": {"name": "Perdew-Becke", "class": "GGA"},
        "ca": {"name": "Perdew-Zunger81", "class": "LDA"},
        "hl": {"name": "Hedin-Lundquist", "class": "LDA"},
        "wi": {"name": "Wigner Interpolation", "class": "LDA"},
    }

    parse_functions = dict(
        LULTRA=_parse_bool,
        LUNSCR=_parse_bool,
        LCOR=_parse_bool,
        LPAW=_parse_bool,
        EATOM=_parse_float,
        RPACOR=_parse_float,
        POMASS=_parse_float,
        ZVAL=_parse_float,
        RCORE=_parse_float,
        RWIGS=_parse_float,
        ENMAX=_parse_float,
        ENMIN=_parse_float,
        EMMIN=_parse_float,
        EAUG=_parse_float,
        DEXC=_parse_float,
        RMAX=_parse_float,
        RAUG=_parse_float,
        RDEP=_parse_float,
        RDEPT=_parse_float,
        QCUT=_parse_float,
        QGAM=_parse_float,
        RCLOC=_parse_float,
        IUNSCR=_parse_int,
        ICORE=_parse_int,
        NDATA=_parse_int,
        VRHFIN=str.strip,
        LEXCH=str.strip,
        TITEL=str.strip,
        STEP=_parse_list,
        RRKJ=_parse_list,
        GGA=_parse_list,
        SHA256=str.strip,
        COPYR=str.strip,
    )

    # used for POTCAR validation
    potcar_summary_stats = loadfn(POTCAR_STATS_PATH)

    def __init__(self, data: str, symbol: str | None = None) -> None:
        """
        Args:
            data (str): Complete and single POTCAR file as a string.
            symbol (str): POTCAR symbol corresponding to the filename suffix e.g. "Tm_3" for POTCAR.TM_3".
                If not given, pymatgen will attempt to extract the symbol from the file itself. This is
                not always reliable!
        """
        self.data = data  # raw POTCAR as a string

        # VASP parses header in vasprun.xml and this differs from the titel
        self.header = data.split("\n")[0].strip()

        match = re.search(r"(?s)(parameters from PSCTR are:.*?END of PSCTR-controll parameters)", data)
        search_lines = match.group(1) if match else ""

        keywords = {}
        for key, val in re.findall(r"(\S+)\s*=\s*(.*?)(?=;|$)", search_lines, flags=re.MULTILINE):
            try:
                keywords[key] = self.parse_functions[key](val)  # type: ignore
            except KeyError:
                warnings.warn(f"Ignoring unknown variable type {key}")

        PSCTR: dict[str, Any] = {}

        array_search = re.compile(r"(-*[0-9.]+)")
        orbitals = []
        descriptions = []
        atomic_config_match = re.search(r"(?s)Atomic configuration(.*?)Description", search_lines)
        if atomic_config_match:
            lines = atomic_config_match.group(1).splitlines()
            match = re.search(r"([0-9]+)", lines[1])
            num_entries = int(match.group(1)) if match else 0
            PSCTR["nentries"] = num_entries
            for line in lines[3:]:
                orbit = array_search.findall(line)
                if orbit:
                    orbitals.append(
                        Orbital(int(orbit[0]), int(orbit[1]), float(orbit[2]), float(orbit[3]), float(orbit[4]))
                    )
            PSCTR["Orbitals"] = tuple(orbitals)

        description_string = re.search(
            r"(?s)Description\s*\n(.*?)Error from kinetic energy argument \(eV\)",
            search_lines,
        )
        if description_string:
            for line in description_string.group(1).splitlines():
                description = array_search.findall(line)
                if description:
                    descriptions.append(
                        OrbitalDescription(
                            int(description[0]),
                            float(description[1]),
                            int(description[2]),
                            float(description[3]),
                            int(description[4]) if len(description) > 4 else None,
                            float(description[5]) if len(description) > 4 else None,
                        )
                    )

        if descriptions:
            PSCTR["OrbitalDescriptions"] = tuple(descriptions)

        rrkj_kinetic_energy_string = re.search(
            r"(?s)Error from kinetic energy argument \(eV\)\s*\n(.*?)END of PSCTR-controll parameters",
            search_lines,
        )
        rrkj_array = []
        if rrkj_kinetic_energy_string:
            for line in rrkj_kinetic_energy_string.group(1).splitlines():
                if "=" not in line:
                    rrkj_array += _parse_list(line.strip("\n"))
            if rrkj_array:
                PSCTR["RRKJ"] = tuple(rrkj_array)

        self.keywords = dict(sorted({**PSCTR, **keywords}.items()))

        if symbol:
            self._symbol = symbol
        else:
            try:
                self._symbol = keywords["TITEL"].split(" ")[1].strip()
            except IndexError:
                self._symbol = keywords["TITEL"].strip()

        # Compute the POTCAR meta to check them against the database of known metadata,
        # and possibly SHA256 hashes contained in the file itself.

        if not self.is_valid:
            warnings.warn(
                f"POTCAR data with symbol {self.symbol} is not known to pymatgen. Your "
                "POTCAR may be corrupted or pymatgen's POTCAR database is incomplete.",
                UnknownPotcarWarning,
            )

    def __str__(self) -> str:
        return self.data + "\n"

    @property
    def electron_configuration(self) -> list[tuple[int, str, int]] | None:
        """Electronic configuration of the PotcarSingle."""
        if not self.nelectrons.is_integer():
            warnings.warn("POTCAR has non-integer charge, electron configuration not well-defined.")
            return None
        el = Element.from_Z(self.atomic_no)
        full_config = el.full_electronic_structure
        nelect = self.nelectrons
        config = []
        while nelect > 0:
            e = full_config.pop(-1)
            config.append(e)
            nelect -= e[-1]
        return config

    def write_file(self, filename: str) -> None:
        """
        Write PotcarSingle to a file.

        Args:
            filename (str): Filename to write to.
        """
        with zopen(filename, "wt") as file:
            file.write(str(self))

    @classmethod
    def from_file(cls, filename: str) -> PotcarSingle:
        """
        Reads PotcarSingle from file.

        :param filename: Filename.

        Returns:
            PotcarSingle
        """
        match = re.search(r"(?<=POTCAR\.)(.*)(?=.gz)", str(filename))
        symbol = match.group(0) if match else ""

        try:
            with zopen(filename, "rt") as file:
                return cls(file.read(), symbol=symbol or None)
        except UnicodeDecodeError:
            warnings.warn("POTCAR contains invalid unicode errors. We will attempt to read it by ignoring errors.")

            with codecs.open(filename, "r", encoding="utf-8", errors="ignore") as file:
                return cls(file.read(), symbol=symbol or None)

    @classmethod
    def from_symbol_and_functional(cls, symbol: str, functional: str | None = None):
        """
        Makes a PotcarSingle from a symbol and functional.

        Args:
            symbol (str): Symbol, e.g., Li_sv
            functional (str): Functional, e.g., PBE

        Returns:
            PotcarSingle
        """
        functional = functional or SETTINGS.get("PMG_DEFAULT_FUNCTIONAL", "PBE")
        assert isinstance(functional, str)  # mypy type narrowing
        funcdir = cls.functional_dir[functional]
        PMG_VASP_PSP_DIR = SETTINGS.get("PMG_VASP_PSP_DIR")
        if PMG_VASP_PSP_DIR is None:
            raise ValueError(
                f"No POTCAR for {symbol} with {functional=} found. Please set the PMG_VASP_PSP_DIR in .pmgrc.yaml."
            )
        paths_to_try = [
            os.path.join(PMG_VASP_PSP_DIR, funcdir, f"POTCAR.{symbol}"),
            os.path.join(PMG_VASP_PSP_DIR, funcdir, symbol, "POTCAR"),
        ]
        for path in paths_to_try:
            path = os.path.expanduser(path)
            path = zpath(path)
            if os.path.isfile(path):
                return cls.from_file(path)
        raise OSError(
            f"You do not have the right POTCAR with {functional=} and {symbol=} "
            f"in your {PMG_VASP_PSP_DIR=}. Paths tried: {paths_to_try}"
        )

    @property
    def element(self) -> str:
        """Attempt to return the atomic symbol based on the VRHFIN keyword."""
        element = self.keywords["VRHFIN"].split(":")[0].strip()
        try:
            return Element(element).symbol
        except ValueError:
            # VASP incorrectly gives the element symbol for Xe as "X"
            # Some potentials, e.g., Zr_sv, gives the symbol as r.
            if element == "X":
                return "Xe"
            return Element(self.symbol.split("_")[0]).symbol

    @property
    def atomic_no(self) -> int:
        """Attempt to return the atomic number based on the VRHFIN keyword."""
        return Element(self.element).Z

    @property
    def nelectrons(self) -> float:
        """Number of electrons"""
        return self.zval

    @property
    def symbol(self) -> str:
        """The POTCAR symbol, e.g. W_pv"""
        return self._symbol

    @property
    def potential_type(self) -> Literal["NC", "PAW", "US"]:
        """Type of PSP. E.g., US, PAW, etc."""
        if self.lultra:
            return "US"
        if self.lpaw:
            return "PAW"
        return "NC"

    @property
    def functional(self) -> str | None:
        """Functional associated with PotcarSingle."""
        return self.functional_tags.get(self.LEXCH.lower(), {}).get("name")

    @property
    def functional_class(self):
        """Functional class associated with PotcarSingle."""
        return self.functional_tags.get(self.LEXCH.lower(), {}).get("class")

    def verify_potcar(self) -> tuple[bool, bool]:
        """
        Attempts to verify the integrity of the POTCAR data.

        This method checks the whole file (removing only the SHA256
        metadata) against the SHA256 hash in the header if this is found.
        If no SHA256 hash is found in the file, the file hash (md5 hash of the
        whole file) is checked against all POTCAR file hashes known to pymatgen.

        Returns:
            tuple[bool, bool]: has_sha256 and passed_hash_check are returned.
        """
        if self.hash_sha256_from_file:
            has_sha256 = True
            hash_is_valid = self.hash_sha256_from_file == self.sha256_computed_file_hash
        else:
            has_sha256 = False
            # if no sha256 hash is found in the POTCAR file, compare the whole
            # file with known potcar file hashes.
            md5_file_hash = self.md5_computed_file_hash
            hash_is_valid = md5_file_hash in VASP_POTCAR_HASHES
        return has_sha256, hash_is_valid

    def identify_potcar(
        self, mode: Literal["data", "file"] = "data", data_tol: float = 1e-6
    ) -> tuple[list[str], list[str]]:
        """
        Identify the symbol and compatible functionals associated with this PotcarSingle.

        This method checks the summary statistics of either the POTCAR metadadata
        (PotcarSingle._summary_stats[key]["header"] for key in ("keywords", "stats") )
        or the entire POTCAR file (PotcarSingle._summary_stats) against a database
        of hashes for POTCARs distributed with VASP 5.4.4.

        Args:
            mode ('data' | 'file'): 'data' mode checks the POTCAR header keywords and stats only
                while 'file' mode checks the entire summary stats.
            data_tol (float): Tolerance for comparing the summary statistics of the POTCAR
                with the reference statistics.

        Returns:
            symbol (list): List of symbols associated with the PotcarSingle
            potcar_functionals (list): List of potcar functionals associated with
                the PotcarSingle
        """
        if mode == "data":
            check_modes = ["header"]
        elif mode == "file":
            check_modes = ["header", "data"]
        else:
            raise ValueError(f"Bad {mode=}. Choose 'data' or 'file'.")

        identity: dict[str, list] = {"potcar_functionals": [], "potcar_symbols": []}
        for func in self.functional_dir:
            for ref_psp in self.potcar_summary_stats[func].get(self.TITEL.replace(" ", ""), []):
                if self.VRHFIN.replace(" ", "") != ref_psp["VRHFIN"]:
                    continue

                key_match = all(
                    set(ref_psp["keywords"][key]) == set(self._summary_stats["keywords"][key])  # type: ignore[index]
                    for key in check_modes
                )

                data_diff = [
                    abs(ref_psp["stats"][key][stat] - self._summary_stats["stats"][key][stat])  # type: ignore[index]
                    for stat in ["MEAN", "ABSMEAN", "VAR", "MIN", "MAX"]
                    for key in check_modes
                ]

                data_match = all(np.array(data_diff) < data_tol)

                if key_match and data_match:
                    identity["potcar_functionals"].append(func)
                    identity["potcar_symbols"].append(ref_psp["symbol"])

        for key in identity:
            if len(identity[key]) == 0:
                # the two keys are set simultaneously, either key being zero indicates no match
                return [], []
            identity[key] = list(set(identity[key]))

        return identity["potcar_functionals"], identity["potcar_symbols"]

    def identify_potcar_hash_based(self, mode: Literal["data", "file"] = "data"):
        """
        Identify the symbol and compatible functionals associated with this PotcarSingle.

        This method checks the md5 hash of either the POTCAR metadadata (PotcarSingle.md5_header_hash)
        or the entire POTCAR file (PotcarSingle.md5_computed_file_hash) against a database
        of hashes for POTCARs distributed with VASP 5.4.4.

        Args:
            mode ('data' | 'file'): 'data' mode checks the hash of the POTCAR metadata in self.keywords,
                while 'file' mode checks the hash of the entire POTCAR file.

        Returns:
            symbol (list): List of symbols associated with the PotcarSingle
            potcar_functionals (list): List of potcar functionals associated with
                the PotcarSingle
        """
        # Dict to translate the sets in the .json file to the keys used in DictSet
        mapping_dict = {
            "potUSPP_GGA": {
                "pymatgen_key": "PW91_US",
                "vasp_description": "Ultrasoft pseudo potentials"
                "for LDA and PW91 (dated 2002-08-20 and 2002-04-08,"
                "respectively). These files are outdated, not"
                "supported and only distributed as is.",
            },
            "potUSPP_LDA": {
                "pymatgen_key": "LDA_US",
                "vasp_description": "Ultrasoft pseudo potentials"
                "for LDA and PW91 (dated 2002-08-20 and 2002-04-08,"
                "respectively). These files are outdated, not"
                "supported and only distributed as is.",
            },
            "potpaw_GGA": {
                "pymatgen_key": "PW91",
                "vasp_description": "The LDA, PW91 and PBE PAW datasets"
                "(snapshot: 05-05-2010, 19-09-2006 and 06-05-2010,"
                "respectively). These files are outdated, not"
                "supported and only distributed as is.",
            },
            "potpaw_LDA": {
                "pymatgen_key": "Perdew-Zunger81",
                "vasp_description": "The LDA, PW91 and PBE PAW datasets"
                "(snapshot: 05-05-2010, 19-09-2006 and 06-05-2010,"
                "respectively). These files are outdated, not"
                "supported and only distributed as is.",
            },
            "potpaw_LDA.52": {
                "pymatgen_key": "LDA_52",
                "vasp_description": "LDA PAW datasets version 52,"
                "including the early GW variety (snapshot 19-04-2012)."
                "When read by VASP these files yield identical results"
                "as the files distributed in 2012 ('unvie' release).",
            },
            "potpaw_LDA.54": {
                "pymatgen_key": "LDA_54",
                "vasp_description": "LDA PAW datasets version 54,"
                "including the GW variety (original release 2015-09-04)."
                "When read by VASP these files yield identical results as"
                "the files distributed before.",
            },
            "potpaw_PBE": {
                "pymatgen_key": "PBE",
                "vasp_description": "The LDA, PW91 and PBE PAW datasets"
                "(snapshot: 05-05-2010, 19-09-2006 and 06-05-2010,"
                "respectively). These files are outdated, not"
                "supported and only distributed as is.",
            },
            "potpaw_PBE.52": {
                "pymatgen_key": "PBE_52",
                "vasp_description": "PBE PAW datasets version 52,"
                "including early GW variety (snapshot 19-04-2012)."
                "When read by VASP these files yield identical"
                "results as the files distributed in 2012.",
            },
            "potpaw_PBE.54": {
                "pymatgen_key": "PBE_54",
                "vasp_description": "PBE PAW datasets version 54,"
                "including the GW variety (original release 2015-09-04)."
                "When read by VASP these files yield identical results as"
                "the files distributed before.",
            },
            "unvie_potpaw.52": {
                "pymatgen_key": "unvie_LDA_52",
                "vasp_description": "files released previously"
                "for vasp.5.2 (2012-04) and vasp.5.4 (2015-09-04) by univie.",
            },
            "unvie_potpaw.54": {
                "pymatgen_key": "unvie_LDA_54",
                "vasp_description": "files released previously"
                "for vasp.5.2 (2012-04) and vasp.5.4 (2015-09-04) by univie.",
            },
            "unvie_potpaw_PBE.52": {
                "pymatgen_key": "unvie_PBE_52",
                "vasp_description": "files released previously"
                "for vasp.5.2 (2012-04) and vasp.5.4 (2015-09-04) by univie.",
            },
            "unvie_potpaw_PBE.54": {
                "pymatgen_key": "unvie_PBE_52",
                "vasp_description": "files released previously"
                "for vasp.5.2 (2012-04) and vasp.5.4 (2015-09-04) by univie.",
            },
        }

        if mode == "data":
            hash_db = PYMATGEN_POTCAR_HASHES
            potcar_hash = self.md5_header_hash
        elif mode == "file":
            hash_db = VASP_POTCAR_HASHES
            potcar_hash = self.md5_computed_file_hash
        else:
            raise ValueError(f"Bad {mode=}. Choose 'data' or 'file'.")

        identity = hash_db.get(potcar_hash)

        if identity:
            # convert the potcar_functionals from the .json dict into the functional
            # keys that pymatgen uses
            potcar_functionals = [*{mapping_dict[i]["pymatgen_key"] for i in identity["potcar_functionals"]}]

            return potcar_functionals, identity["potcar_symbols"]
        return [], []

    @property
    def hash_sha256_from_file(self) -> str | None:
        """SHA256 hash of the POTCAR file as read from the file. None if no SHA256 hash is found."""
        if sha256 := getattr(self, "SHA256", None):
            return sha256.split()[0]
        return None

    @property
    def sha256_computed_file_hash(self) -> str:
        """Computes a SHA256 hash of the PotcarSingle EXCLUDING lines starting with 'SHA256' and 'COPYR'."""
        # we have to remove lines with the hash itself and the copyright
        # notice to get the correct hash.
        potcar_list = self.data.split("\n")
        potcar_to_hash = [line for line in potcar_list if not line.strip().startswith(("SHA256", "COPYR"))]
        potcar_to_hash_str = "\n".join(potcar_to_hash)
        return sha256(potcar_to_hash_str.encode("utf-8")).hexdigest()

    @property
    def md5_computed_file_hash(self) -> str:
        """md5 hash of the entire PotcarSingle."""
        # usedforsecurity=False needed in FIPS mode (Federal Information Processing Standards)
        # https://github.com/materialsproject/pymatgen/issues/2804
        md5 = hashlib.new("md5", usedforsecurity=False)  # hashlib.md5(usedforsecurity=False) is py39+
        md5.update(self.data.encode("utf-8"))
        return md5.hexdigest()

    @property
    def md5_header_hash(self) -> str:
        """Computes a md5 hash of the metadata defining the PotcarSingle."""
        hash_str = ""
        for k, v in self.keywords.items():
            # for newer POTCARS we have to exclude 'SHA256' and 'COPYR lines
            # since they were not used in the initial hashing
            if k in ("nentries", "Orbitals", "SHA256", "COPYR"):
                continue
            hash_str += f"{k}"
            if isinstance(v, bool):
                hash_str += f"{v}"
            elif isinstance(v, float):
                hash_str += f"{v:.3f}"
            elif isinstance(v, int):
                hash_str += f"{v}"
            elif isinstance(v, (tuple, list)):
                for item in v:
                    if isinstance(item, float):
                        hash_str += f"{item:.3f}"
                    elif isinstance(item, (Orbital, OrbitalDescription)):
                        for item_v in item:
                            if isinstance(item_v, (int, str)):
                                hash_str += f"{item_v}"
                            elif isinstance(item_v, float):
                                hash_str += f"{item_v:.3f}"
                            else:
                                hash_str += f"{item_v}" if item_v else ""
            else:
                hash_str += v.replace(" ", "")

        self.hash_str = hash_str
        # usedforsecurity=False needed in FIPS mode (Federal Information Processing Standards)
        # https://github.com/materialsproject/pymatgen/issues/2804
        md5 = hashlib.new("md5", usedforsecurity=False)  # hashlib.md5(usedforsecurity=False) is py39+
        md5.update(hash_str.lower().encode("utf-8"))
        return md5.hexdigest()

    @property
    def is_valid(self) -> bool:
        """
        Check that POTCAR matches reference metadata.
        Parsed metadata is stored in self._summary_stats as a human-readable dict,
            self._summary_stats = {
                "keywords": {
                    "header": list[str],
                    "data": list[str],
                },
                "stats": {
                    "header": dict[float],
                    "data": dict[float],
                },
            }

        Rationale:
        Each POTCAR is structured as
            Header (self.keywords)
            Data (actual pseudopotential values in data blocks)

        For the Data block of POTCAR, there are unformatted data blocks
        of unknown length and contents/data type, e.g., you might see
            <float> <bool>
            <Data Keyword>
            <int> <int> <float>
            <float> ... <float>
            <Data Keyword>
            <float> ... <float>
        but this is impossible to process algorithmically without a full POTCAR schema.
        Note also that POTCARs can contain **different** data keywords

        All keywords found in the header, essentially self.keywords, and the data block
        (<Data Keyword> above) are stored in self._summary_stats["keywords"]

        To avoid issues of copyright, statistics (mean, mean of abs vals, variance, max, min)
        for the numeric values in the header and data sections of POTCAR are stored
        in self._summary_stats["stats"]

        tol is then used to match statistical values within a tolerance
        """

        possible_potcar_matches = []
        """
        Some POTCARs have an LEXCH (functional used to generate the POTCAR)
        with the expected functional, e.g., the C_d POTCAR for PBE is actually an
        LDA pseudopotential.

        Thus we have to look for matches in all POTCAR dirs, not just the ones with
        consistent values of LEXCH
        """
        for func in self.functional_dir:
            for titel_no_spc in self.potcar_summary_stats[func]:
                if self.TITEL.replace(" ", "") == titel_no_spc:
                    for potcar_subvariant in self.potcar_summary_stats[func][titel_no_spc]:
                        if self.VRHFIN.replace(" ", "") == potcar_subvariant["VRHFIN"]:
                            possible_potcar_matches.append(
                                {
                                    "POTCAR_FUNCTIONAL": func,
                                    "TITEL": titel_no_spc,
                                    **potcar_subvariant,
                                }
                            )

        def parse_fortran_style_str(input_str: str) -> Any:
            """Parse any input string as bool, int, float, or failing that, str.
            Used to parse FORTRAN-generated POTCAR files where it's unknown
            a priori what type of data will be encountered.
            """
            input_str = input_str.strip()

            if input_str.lower() in ("t", "f", "true", "false"):
                return input_str[0].lower() == "t"

            if input_str.upper() == input_str.lower() and input_str[0].isnumeric():
                if "." in input_str:
                    """
                    NB: fortran style floats always include a decimal point.
                        While you can set, e.g., x = 1E4, you cannot print/write x without
                        a decimal point:
                            `write(6,*) x`          -->   `10000.0000` in stdout
                            `write(6,'(E10.0)') x`  -->   segfault
                        The (E10.0) means write an exponential-format number with 10
                            characters before the decimal, and 0 characters after
                    """
                    return float(input_str)
                return int(input_str)
            try:
                return float(input_str)
            except ValueError:
                return input_str

        psp_keys, psp_vals = [], []
        potcar_body = self.data.split("END of PSCTR-controll parameters\n")[1]
        for row in re.split(r"\n+|;", potcar_body):  # FORTRAN allows ; to delimit multiple lines merged into 1 line
            tmp_str = ""
            for raw_val in row.split():
                parsed_val = parse_fortran_style_str(raw_val)
                if isinstance(parsed_val, str):
                    tmp_str += parsed_val.strip()
                elif isinstance(parsed_val, (float, int)):
                    psp_vals.append(parsed_val)
            if len(tmp_str) > 0:
                psp_keys.append(tmp_str.lower())

        keyword_vals = []
        for kwd in self.keywords:
            val = self.keywords[kwd]
            if isinstance(val, bool):
                # has to come first since bools are also ints
                keyword_vals.append(1.0 if val else 0.0)
            elif isinstance(val, (float, int)):
                keyword_vals.append(val)
            elif hasattr(val, "__len__"):
                keyword_vals += [num for num in val if (isinstance(num, (float, int)))]

        def data_stats(data_list: Sequence) -> dict:
            """Used for hash-less and therefore less brittle POTCAR validity checking."""
            arr = np.array(data_list)
            return {
                "MEAN": np.mean(arr),
                "ABSMEAN": np.mean(np.abs(arr)),
                "VAR": np.mean(arr**2),
                "MIN": arr.min(),
                "MAX": arr.max(),
            }

        # NB: to add future summary stats in a way that's consistent with PMG,
        # it's easiest to save the summary stats as an attr of PotcarSingle
        self._summary_stats = {  # for this PotcarSingle instance
            "keywords": {
                "header": [kwd.lower() for kwd in self.keywords],
                "data": psp_keys,
            },
            "stats": {
                "header": data_stats(keyword_vals),
                "data": data_stats(psp_vals),
            },
        }

        data_match_tol = 1e-6
        for ref_psp in possible_potcar_matches:
            key_match = all(
                set(ref_psp["keywords"][key]) == set(self._summary_stats["keywords"][key])  # type: ignore
                for key in ["header", "data"]
            )

            data_diff = [
                abs(ref_psp["stats"][key][stat] - self._summary_stats["stats"][key][stat])  # type: ignore
                for stat in ["MEAN", "ABSMEAN", "VAR", "MIN", "MAX"]
                for key in ["header", "data"]
            ]
            data_match = all(np.array(data_diff) < data_match_tol)

            if key_match and data_match:
                return True

        return False

    def __getattr__(self, attr: str) -> Any:
        """Delegates attributes to keywords. For example, you can use potcarsingle.enmax to get the ENMAX of the POTCAR.

        For float type properties, they are converted to the correct float. By
        default, all energies in eV and all length scales are in Angstroms.
        """
        try:
            return self.keywords[attr.upper()]
        except Exception:
            raise AttributeError(attr)

    def __repr__(self) -> str:
        cls_name = type(self).__name__
        symbol, functional = self.symbol, self.functional
        TITEL, VRHFIN = self.keywords["TITEL"], self.keywords["VRHFIN"]
        TITEL, VRHFIN, n_valence_elec = (self.keywords.get(key) for key in ("TITEL", "VRHFIN", "ZVAL"))
        return f"{cls_name}({symbol=}, {functional=}, {TITEL=}, {VRHFIN=}, {n_valence_elec=:.0f})"


def _gen_potcar_summary_stats(
    append: bool = False, vasp_psp_dir: str | None = None, summary_stats_filename: str = POTCAR_STATS_PATH
):
    """
    This function solely intended to be used for PMG development to regenerate the
    potcar-summary-stats.json.bz2 file used to validate POTCARs

    THIS FUNCTION IS DESTRUCTIVE. It will completely overwrite your potcar-summary-stats.json.bz2.

    Args:
        append (bool): Change whether data is appended to the existing potcar-summary-stats.json.bz2,
            or if a completely new file is generated. Defaults to False.
        PMG_VASP_PSP_DIR (str): Change where this function searches for POTCARs
            defaults to the PMG_VASP_PSP_DIR environment variable if not set. Defaults to None.
        summary_stats_filename (str): Name of the output summary stats file. Defaults to
            '<pymatgen_install_dir>/io/vasp/potcar-summary-stats.json.bz2'.
    """
    func_dir_exist: dict[str, str] = {}
    vasp_psp_dir = vasp_psp_dir or SETTINGS.get("PMG_VASP_PSP_DIR")
    for func in PotcarSingle.functional_dir:
        cpsp_dir = f"{vasp_psp_dir}/{PotcarSingle.functional_dir[func]}"
        if os.path.isdir(cpsp_dir):
            func_dir_exist[func] = PotcarSingle.functional_dir[func]
        else:
            warnings.warn(f"missing {PotcarSingle.functional_dir[func]} POTCAR directory")

    # use append = True if a new POTCAR library is released to add new summary stats
    # without completely regenerating the dict of summary stats
    # use append = False to completely regenerate the summary stats dict
    new_summary_stats = loadfn(summary_stats_filename) if append else {}

    for func in func_dir_exist:
        new_summary_stats.setdefault(func, {})  # initialize dict if key missing

        potcar_list = [
            *glob(f"{vasp_psp_dir}/{func_dir_exist[func]}/POTCAR*"),
            *glob(f"{vasp_psp_dir}/{func_dir_exist[func]}/*/POTCAR*"),
        ]
        for potcar in potcar_list:
            psp = PotcarSingle.from_file(potcar)
            titel_key = psp.TITEL.replace(" ", "")

            # some POTCARs have the same TITEL, but are named differently
            # e.g., there is an "original" PBE POTCAR.Fe_pv and a POTCAR.Fe_pv_new
            # which share a TITEL but differ in their contents
            if titel_key not in new_summary_stats[func]:
                new_summary_stats[func][titel_key] = []
            new_summary_stats[func][titel_key].append(
                {
                    "LEXCH": psp.LEXCH,
                    "VRHFIN": psp.VRHFIN.replace(" ", ""),
                    "symbol": psp.symbol,
                    "ZVAL": psp.ZVAL,
                    **psp._summary_stats,
                }
            )

    dumpfn(new_summary_stats, summary_stats_filename)


class Potcar(list, MSONable):
    """
    Object for reading and writing POTCAR files for calculations. Consists of a
    list of PotcarSingle.
    """

    FUNCTIONAL_CHOICES = tuple(PotcarSingle.functional_dir)

    def __init__(
        self,
        symbols: Sequence[str] | None = None,
        functional: str | None = None,
        sym_potcar_map: dict[str, str] | None = None,
    ):
        """
        Args:
            symbols (list[str]): Element symbols for POTCAR. This should correspond
                to the symbols used by VASP. E.g., "Mg", "Fe_pv", etc.
            functional (str): Functional used. To know what functional options
                there are, use Potcar.FUNCTIONAL_CHOICES. Note that VASP has
                different versions of the same functional. By default, the old
                PBE functional is used. If you want the newer ones, use PBE_52 or
                PBE_54. Note that if you intend to compare your results with the
                Materials Project, you should use the default setting. You can also
                override the default by setting PMG_DEFAULT_FUNCTIONAL in your
                .pmgrc.yaml.
            sym_potcar_map (dict): Allows a user to specify a specific element
                symbol to raw POTCAR mapping.
        """
        if functional is None:
            functional = SETTINGS.get("PMG_DEFAULT_FUNCTIONAL", "PBE")
        super().__init__()
        self.functional = functional
        if symbols is not None:
            self.set_symbols(symbols, functional, sym_potcar_map)

    def __iter__(self) -> Iterator[PotcarSingle]:  # __iter__ only needed to supply type hint
        # so for psingle in Potcar() is correctly inferred as PotcarSingle
        return super().__iter__()

    def as_dict(self):
        """MSONable dict representation"""
        return {
            "functional": self.functional,
            "symbols": self.symbols,
            "@module": type(self).__module__,
            "@class": type(self).__name__,
        }

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation

        Returns:
            Potcar
        """
        return Potcar(symbols=d["symbols"], functional=d["functional"])

    @classmethod
    def from_file(cls, filename: str):
        """
        Reads Potcar from file.

        :param filename: Filename

        Returns:
            Potcar
        """
        with zopen(filename, "rt") as f:
            fdata = f.read()
        potcar = cls()

        functionals = []
        for psingle_str in fdata.split("End of Dataset"):
            if p_strip := psingle_str.strip():
                psingle = PotcarSingle(p_strip + "\nEnd of Dataset\n")
                potcar.append(psingle)
                functionals.append(psingle.functional)
        if len(set(functionals)) != 1:
            raise ValueError("File contains incompatible functionals!")
        potcar.functional = functionals[0]
        return potcar

    def __str__(self) -> str:
        return "\n".join(str(potcar).strip("\n") for potcar in self) + "\n"

    def write_file(self, filename: str) -> None:
        """
        Write Potcar to a file.

        Args:
            filename (str): filename to write to.
        """
        with zopen(filename, "wt") as f:
            f.write(str(self))

    @property
    def symbols(self):
        """Get the atomic symbols of all the atoms in the POTCAR file."""
        return [psingle.symbol for psingle in self]

    @symbols.setter
    def symbols(self, symbols):
        self.set_symbols(symbols, functional=self.functional)

    @property
    def spec(self):
        """Get the atomic symbols and hash of all the atoms in the POTCAR file."""
        return [{"symbol": psingle.symbol, "hash": psingle.md5_computed_file_hash} for psingle in self]

    def set_symbols(
        self, symbols: Sequence[str], functional: str | None = None, sym_potcar_map: dict[str, str] | None = None
    ):
        """
        Initialize the POTCAR from a set of symbols. Currently, the POTCARs can
        be fetched from a location specified in .pmgrc.yaml. Use pmg config
        to add this setting.

        Args:
            symbols (list[str]): A list of element symbols
            functional (str): The functional to use. If None, the setting
                PMG_DEFAULT_FUNCTIONAL in .pmgrc.yaml is used, or if this is
                not set, it will default to PBE.
            sym_potcar_map (dict): A map of symbol:raw POTCAR string. If
                sym_potcar_map is specified, POTCARs will be generated from
                the given map data rather than the config file location.
        """
        del self[:]
        if sym_potcar_map:
            self.extend(PotcarSingle(sym_potcar_map[el]) for el in symbols)
        else:
            self.extend(PotcarSingle.from_symbol_and_functional(el, functional) for el in symbols)


class VaspInput(dict, MSONable):
    """Class to contain a set of vasp input objects corresponding to a run."""

    def __init__(self, incar, kpoints, poscar, potcar, optional_files=None, **kwargs):
        """
        Initializes a VaspInput object with the given input files.

        Args:
            incar (Incar): The Incar object.
            kpoints (Kpoints): The Kpoints object.
            poscar (Poscar): The Poscar object.
            potcar (Potcar): The Potcar object.
            optional_files (dict): Other input files supplied as a dict of {filename: object}.
                The object should follow standard pymatgen conventions in implementing a
                as_dict() and from_dict method.
            **kwargs: Additional keyword arguments to be stored in the VaspInput object.
        """
        super().__init__(**kwargs)
        self.update({"INCAR": incar, "KPOINTS": kpoints, "POSCAR": poscar, "POTCAR": potcar})
        if optional_files is not None:
            self.update(optional_files)

    def __str__(self):
        output = []
        for k, v in self.items():
            output.extend((k, str(v), ""))
        return "\n".join(output)

    def as_dict(self):
        """MSONable dict."""
        dct = {k: v.as_dict() for k, v in self.items()}
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        return dct

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation.

        Returns:
            VaspInput
        """
        dec = MontyDecoder()
        sub_d = {"optional_files": {}}
        for k, v in d.items():
            if k in ["INCAR", "POSCAR", "POTCAR", "KPOINTS"]:
                sub_d[k.lower()] = dec.process_decoded(v)
            elif k not in ["@module", "@class"]:
                sub_d["optional_files"][k] = dec.process_decoded(v)
        return cls(**sub_d)

    def write_input(self, output_dir=".", make_dir_if_not_present=True):
        """
        Write VASP input to a directory.

        Args:
            output_dir (str): Directory to write to. Defaults to current
                directory (".").
            make_dir_if_not_present (bool): Create the directory if not
                present. Defaults to True.
        """
        if make_dir_if_not_present:
            os.makedirs(output_dir, exist_ok=True)
        for k, v in self.items():
            if v is not None:
                with zopen(os.path.join(output_dir, k), "wt") as f:
                    f.write(str(v))

    @classmethod
    def from_directory(cls, input_dir, optional_files=None):
        """
        Read in a set of VASP input from a directory. Note that only the
        standard INCAR, POSCAR, POTCAR and KPOINTS files are read unless
        optional_filenames is specified.

        Args:
            input_dir (str): Directory to read VASP input from.
            optional_files (dict): Optional files to read in as well as a
                dict of {filename: Object type}. Object type must have a
                static method from_file.
        """
        sub_d = {}
        for fname, ftype in [
            ("INCAR", Incar),
            ("KPOINTS", Kpoints),
            ("POSCAR", Poscar),
            ("POTCAR", Potcar),
        ]:
            try:
                full_zpath = zpath(os.path.join(input_dir, fname))
                sub_d[fname.lower()] = ftype.from_file(full_zpath)
            except FileNotFoundError:  # handle the case where there is no KPOINTS file
                sub_d[fname.lower()] = None

        sub_d["optional_files"] = {}
        if optional_files is not None:
            for fname, ftype in optional_files.items():
                sub_d["optional_files"][fname] = ftype.from_file(os.path.join(input_dir, fname))
        return cls(**sub_d)

    def run_vasp(
        self,
        run_dir: PathLike = ".",
        vasp_cmd: list | None = None,
        output_file: PathLike = "vasp.out",
        err_file: PathLike = "vasp.err",
    ):
        """
        Write input files and run VASP.

        :param run_dir: Where to write input files and do the run.
        :param vasp_cmd: Args to be supplied to run VASP. Otherwise, the
            PMG_VASP_EXE in .pmgrc.yaml is used.
        :param output_file: File to write output.
        :param err_file: File to write err.
        """
        self.write_input(output_dir=run_dir)
        vasp_cmd = vasp_cmd or SETTINGS.get("PMG_VASP_EXE")  # type: ignore[assignment]
        if not vasp_cmd:
            raise ValueError("No VASP executable specified!")
        vasp_cmd = [os.path.expanduser(os.path.expandvars(t)) for t in vasp_cmd]
        if not vasp_cmd:
            raise RuntimeError("You need to supply vasp_cmd or set the PMG_VASP_EXE in .pmgrc.yaml to run VASP.")
        with cd(run_dir), open(output_file, "w") as stdout_file, open(err_file, "w", buffering=1) as stderr_file:
            subprocess.check_call(vasp_cmd, stdout=stdout_file, stderr=stderr_file)
