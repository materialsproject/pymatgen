"""
Classes for reading/manipulating/writing VASP input files.
All major VASP input files.
"""

from __future__ import annotations

import codecs
import hashlib
import itertools
import json
import math
import os
import re
import subprocess
import warnings
from collections import Counter, UserDict
from enum import Enum, unique
from glob import glob
from hashlib import sha256
from pathlib import Path
from shutil import copyfileobj
from typing import TYPE_CHECKING, NamedTuple, cast
from zipfile import ZipFile

import numpy as np
import scipy.constants as const
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
from pymatgen.util.typing import Kpoint, Tuple3Floats, Tuple3Ints, Vector3D

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any, ClassVar, Literal

    from numpy.typing import ArrayLike
    from typing_extensions import Self

    from pymatgen.symmetry.bandstructure import HighSymmKpath
    from pymatgen.util.typing import PathLike


__author__ = "Shyue Ping Ong, Geoffroy Hautier, Rickard Armiento, Vincent L Chevrier, Stephen Dacek"
__copyright__ = "Copyright 2011, The Materials Project"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


class Poscar(MSONable):
    """Represent the data in a POSCAR or CONTCAR file.

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
    ) -> None:
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
        if not structure.is_ordered:
            raise ValueError("Disordered structure with partial occupancies cannot be converted into POSCAR!")

        site_properties: dict[str, Any] = {}

        if selective_dynamics is not None:
            selective_dynamics = np.array(selective_dynamics)
            if not selective_dynamics.all():
                site_properties["selective_dynamics"] = selective_dynamics

        if velocities:
            velocities = np.array(velocities)
            if velocities.any():
                site_properties["velocities"] = velocities

        if predictor_corrector:
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

        if lattice_velocities and np.any(lattice_velocities):
            self.structure.properties["lattice_velocities"] = np.asarray(lattice_velocities)

        self.temperature = -1.0

    def __setattr__(self, name: str, value: Any) -> None:
        if name in {"selective_dynamics", "velocities"} and value is not None and len(value) > 0:
            value = np.array(value)
            dim = value.shape
            if dim[1] != 3 or dim[0] != len(self.structure):
                raise ValueError(f"{name} array must be same length as the structure.")
            value = value.tolist()

        super().__setattr__(name, value)

    def __repr__(self) -> str:
        return self.get_str()

    def __str__(self) -> str:
        """String representation of Poscar file."""
        return self.get_str()

    @property
    def velocities(self) -> ArrayLike | None:
        """Velocities in Poscar."""
        return self.structure.site_properties.get("velocities")

    @velocities.setter
    def velocities(self, velocities: ArrayLike | None) -> None:
        self.structure.add_site_property("velocities", velocities)

    @property
    def selective_dynamics(self) -> ArrayLike | None:
        """Selective dynamics in Poscar."""
        return self.structure.site_properties.get("selective_dynamics")

    @selective_dynamics.setter
    def selective_dynamics(self, selective_dynamics: ArrayLike | None) -> None:
        self.structure.add_site_property("selective_dynamics", selective_dynamics)

    @property
    def predictor_corrector(self) -> ArrayLike | None:
        """Predictor corrector in Poscar."""
        return self.structure.site_properties.get("predictor_corrector")

    @predictor_corrector.setter
    def predictor_corrector(self, predictor_corrector: ArrayLike | None) -> None:
        self.structure.add_site_property("predictor_corrector", predictor_corrector)

    @property
    def predictor_corrector_preamble(self) -> str | None:
        """Predictor corrector preamble in Poscar."""
        return self.structure.properties.get("predictor_corrector_preamble")

    @predictor_corrector_preamble.setter
    def predictor_corrector_preamble(self, predictor_corrector_preamble: str | None) -> None:
        self.structure.properties["predictor_corrector"] = predictor_corrector_preamble

    @property
    def lattice_velocities(self) -> ArrayLike | None:
        """Lattice velocities in Poscar (including the current lattice vectors)."""
        return self.structure.properties.get("lattice_velocities")

    @lattice_velocities.setter
    def lattice_velocities(self, lattice_velocities: ArrayLike | None) -> None:
        self.structure.properties["lattice_velocities"] = np.asarray(lattice_velocities)

    @property
    def site_symbols(self) -> list[str]:
        """Sequence of symbols associated with the Poscar. Similar to 6th line in VASP 5+ POSCAR."""
        syms: list[str] = [site.specie.symbol for site in self.structure]
        return [a[0] for a in itertools.groupby(syms)]

    @property
    def natoms(self) -> list[int]:
        """Sequence of number of sites of each type associated with the Poscar.
        Similar to 7th line in VASP 5+ POSCAR or the 6th line in VASP 4 POSCAR.
        """
        syms: list[str] = [site.specie.symbol for site in self.structure]
        return [len(tuple(a[1])) for a in itertools.groupby(syms)]

    @classmethod
    def from_file(
        cls,
        filename: PathLike,
        check_for_potcar: bool = True,
        read_velocities: bool = True,
        **kwargs: dict[str, Any],
    ) -> Self:
        """
        Read POSCAR from a file.

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
        Poscar, where a warning would be issued. For example, H, He, Li, ....
        This will ensure at least a unique element is assigned to each site and
        any analysis that does not require specific elemental properties should work.

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
            warnings.warn(
                "check_for_POTCAR is deprecated. Use check_for_potcar instead.",
                DeprecationWarning,
            )
            check_for_potcar = cast(bool, kwargs.pop("check_for_POTCAR"))

        dirname: str = os.path.dirname(os.path.abspath(filename))
        names: list[str] | None = None
        if (
            check_for_potcar
            and SETTINGS.get("PMG_POTCAR_CHECKS") is not False
            and (potcars := glob(f"{dirname}/*POTCAR*"))
        ):
            try:
                potcar = Potcar.from_file(min(potcars))
                names = [sym.split("_")[0] for sym in potcar.symbols]
                map(get_el_sp, names)  # ensure valid names
            except Exception:
                names = None

        with zopen(filename, mode="rt") as file:
            return cls.from_str(file.read(), names, read_velocities=read_velocities)

    @classmethod
    def from_str(
        cls,
        data: str,
        default_names: list[str] | None = None,
        read_velocities: bool = True,
    ) -> Self:
        """
        Read POSCAR from a string.

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
        chunks: list[str] = re.split(r"\n\s*\n", data.rstrip(), flags=re.MULTILINE)
        try:
            if chunks[0] == "":
                chunks.pop(0)
                chunks[0] = "\n" + chunks[0]
        except IndexError:
            raise ValueError("Empty POSCAR")

        # Parse positions
        lines: list[str] = list(clean_lines(chunks[0].split("\n"), remove_empty_lines=False))
        comment: str = lines[0]
        scale: float = float(lines[1])
        lattice: np.ndarray = np.array([[float(i) for i in line.split()] for line in lines[2:5]])
        if scale < 0:
            # In VASP, a negative scale factor is treated as a volume. We need
            # to translate this to a proper lattice vector scaling.
            vol: float = abs(np.linalg.det(lattice))
            lattice *= (-scale / vol) ** (1 / 3)
        else:
            lattice *= scale

        vasp5_symbols: bool = False
        atomic_symbols: list[str] = []

        try:
            n_atoms: list[int] = [int(i) for i in lines[5].split()]
            ipos: int = 6

        except ValueError:
            vasp5_symbols = True

            # In VASP 6.x.x, part of the POTCAR hash is written to POSCAR-style strings
            # In VASP 6.4.2 and up, the POTCAR symbol is also written, ex.:
            # ```MgSi
            # 1.0
            # -0.000011    4.138704    0.000002
            # -2.981238    2.069353    3.675251
            # 2.942054    2.069351    4.865237
            # Mg_pv/f474ac0d  Si/79d9987ad87```
            # whereas older VASP 5.x.x POSCAR strings would just have `Mg Si` on the last line

            symbols: list[str] = [symbol.split("/")[0].split("_")[0] for symbol in lines[5].split()]

            # Atoms and number of atoms in POSCAR written with VASP appear on
            # multiple lines when atoms of the same type are not grouped together
            # and more than 20 groups are then defined ...
            # Example :
            # Cr16 Fe35 Ni2
            #    1.00000000000000
            #      8.5415010000000002   -0.0077670000000000   -0.0007960000000000
            #     -0.0077730000000000    8.5224019999999996    0.0105580000000000
            #     -0.0007970000000000    0.0105720000000000    8.5356889999999996
            #    Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Cr   Fe   Ni   Fe   Cr   Fe   Cr
            #    Fe   Ni   Fe   Cr   Fe
            #      1   1   2   4   2   1   1   1     2     1     1     1     4     1     1     1     5     3     6     1
            #      2   1   3   2   5
            # Direct
            #   ...
            n_lines_symbols: int = 1
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

            for idx, n_atom in enumerate(n_atoms):
                atomic_symbols.extend([symbols[idx]] * n_atom)

            ipos = 5 + 2 * n_lines_symbols

        pos_type: str = lines[ipos].split()[0]

        has_selective_dynamics: bool = False
        # Selective dynamics
        if pos_type[0] in "sS":
            has_selective_dynamics = True
            ipos += 1
            pos_type = lines[ipos].split()[0]

        cart: bool = pos_type[0] in "cCkK"
        n_sites: int = sum(n_atoms)

        # If default_names is specified (usually coming from a POTCAR), use
        # them. This is in line with VASP's parsing order that the POTCAR
        # specified is the default used.
        if default_names is not None:
            try:
                atomic_symbols = []
                for idx, n_atom in enumerate(n_atoms):
                    atomic_symbols.extend([default_names[idx]] * n_atom)
                vasp5_symbols = True
            except IndexError:
                pass

        if not vasp5_symbols:
            ind: Literal[3, 6] = 6 if has_selective_dynamics else 3
            try:
                # Check if names are appended at the end of the coordinates
                atomic_symbols = [line.split()[ind] for line in lines[ipos + 1 : ipos + 1 + n_sites]]
                # Ensure symbols are valid elements
                if not all(Element.is_valid_symbol(sym) for sym in atomic_symbols):
                    raise ValueError("Non-valid symbols detected.")
                vasp5_symbols = True

            except (ValueError, IndexError):
                # Defaulting to false names
                atomic_symbols = []
                for idx, n_atom in enumerate(n_atoms, start=1):
                    symbol = Element.from_Z(idx).symbol
                    atomic_symbols.extend([symbol] * n_atom)
                warnings.warn(
                    f"Elements in POSCAR cannot be determined. Defaulting to false names {atomic_symbols}.",
                    BadPoscarWarning,
                )

        # Read the atomic coordinates
        coords: list[list[float]] = []
        selective_dynamics: list[list[bool]] | None = [] if has_selective_dynamics else None
        for idx in range(n_sites):
            tokens: list[str] = lines[ipos + 1 + idx].split()
            crd_scale: float = scale if cart else 1
            coords.append([float(j) * crd_scale for j in tokens[:3]])
            if selective_dynamics is not None:
                # Warn when values contain suspicious entries
                if any(value not in {"T", "F"} for value in tokens[3:6]):
                    warnings.warn(
                        "Selective dynamics values must be either 'T' or 'F'.",
                        BadPoscarWarning,
                    )

                # Warn when elements contains Fluorine (F) (#3539)
                if atomic_symbols[idx] == "F" and len(tokens[3:]) >= 4 and "F" in tokens[3:7]:
                    warnings.warn(
                        (
                            "Selective dynamics toggled with Fluorine element detected. "
                            "Make sure the 4th-6th entry each position line is selective dynamics info."
                        ),
                        BadPoscarWarning,
                    )

                selective_dynamics.append([value == "T" for value in tokens[3:6]])

        # Warn when ALL degrees of freedom relaxed (#3539)
        if selective_dynamics is not None and all(all(i is True for i in in_list) for in_list in selective_dynamics):
            warnings.warn(
                "Ignoring selective dynamics tag, as no ionic degrees of freedom were fixed.",
                BadPoscarWarning,
            )

        struct = Structure(
            lattice,
            atomic_symbols,
            coords,
            to_unit_cell=False,
            validate_proximity=False,
            coords_are_cartesian=cart,
        )

        lattice_velocities: list[list[float]] = []
        velocities: list[list[float]] = []
        predictor_corrector: list = []
        predictor_corrector_preamble: str = ""

        if read_velocities:
            # Parse the lattice velocities and current lattice, if present.
            # The header line should contain "Lattice velocities and vectors"
            # There is no space between the coordinates and this section, so
            # it appears in the lines of the first chunk
            if len(lines) > ipos + n_sites + 1 and lines[ipos + n_sites + 1].lower().startswith("l"):
                for line in lines[ipos + n_sites + 3 : ipos + n_sites + 9]:
                    lattice_velocities.append([float(tok) for tok in line.split()])

            # Parse velocities if any
            if len(chunks) > 1:
                for line in chunks[1].strip().split("\n"):
                    velocities.append([float(tok) for tok in line.split()])

            # Parse the predictor-corrector data
            if len(chunks) > 2:
                lines = chunks[2].strip().split("\n")
                # There are 3 sets of 3xN Predictor corrector parameters
                # So can't be stored as a single set of "site_property"

                # First line in chunk is a key in CONTCAR
                # Second line is POTIM
                # Third line is the thermostat parameters
                predictor_corrector_preamble = f"{lines[0]}\n{lines[1]}\n{lines[2]}"
                # Rest is three sets of parameters, each set contains
                # x, y, z predictor-corrector parameters for every atom in order
                lines = lines[3:]
                for st in range(n_sites):
                    d1 = [float(tok) for tok in lines[st].split()]
                    d2 = [float(tok) for tok in lines[st + n_sites].split()]
                    d3 = [float(tok) for tok in lines[st + 2 * n_sites].split()]
                    predictor_corrector.append([d1, d2, d3])

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

    def get_str(
        self,
        direct: bool = True,
        vasp4_compatible: bool = False,
        significant_figures: int = 16,
    ) -> str:
        """Return a string to be written as a POSCAR file. By default, site
        symbols are written, which is compatible for VASP >= 5.

        Args:
            direct (bool): Whether coordinates are output in direct or
                Cartesian. Defaults to True.
            vasp4_compatible (bool): Set to True to omit site symbols on 6th
                line to maintain backward VASP 4.x compatibility. Defaults
                to False.
            significant_figures (int): Number of significant digits to
                output all quantities. Defaults to 16. Note that positions are
                output in fixed point, while velocities are output in
                scientific format.

        Returns:
            str: representation of POSCAR.
        """
        # This corrects for VASP really annoying bug of crashing on lattices
        # which have triple product < 0. We will just invert the lattice
        # vectors.
        lattice = self.structure.lattice
        if np.linalg.det(lattice.matrix) < 0:
            lattice = Lattice(-lattice.matrix)

        # Add comment and lattice
        format_str: str = f"{{:{significant_figures + 5}.{significant_figures}f}}"
        lines: list[str] = [self.comment, "1.0"]
        for vec in lattice.matrix:
            lines.append(" ".join(format_str.format(c) for c in vec))

        # Add element symbols
        if self.true_names and not vasp4_compatible:
            lines.append(" ".join(self.site_symbols))
        lines.append(" ".join(map(str, self.natoms)))

        if self.selective_dynamics:
            lines.append("Selective dynamics")

        lines.append("direct" if direct else "cartesian")

        # Add ion positions and selective dynamics
        for idx, site in enumerate(self.structure):
            coords: ArrayLike = site.frac_coords if direct else site.coords
            line: str = " ".join(format_str.format(c) for c in coords)
            if self.selective_dynamics is not None:
                sd: list[str] = ["T" if j else "F" for j in self.selective_dynamics[idx]]
                line += f" {sd[0]} {sd[1]} {sd[2]}"
            line += f" {site.species_string}"
            lines.append(line)

        if self.lattice_velocities is not None:
            try:
                lines.extend(["Lattice velocities and vectors", "  1"])
                for velo in self.lattice_velocities:
                    # VASP is strict about the format when reading this quantity
                    lines.append(" ".join(f" {val: .7E}" for val in velo))
            except Exception:
                warnings.warn("Lattice velocities are missing or corrupted.", BadPoscarWarning)

        if self.velocities:
            try:
                lines.append("")
                for velo in self.velocities:
                    lines.append(" ".join(format_str.format(val) for val in velo))
            except Exception:
                warnings.warn("Velocities are missing or corrupted.", BadPoscarWarning)

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
                    "Preamble information missing or corrupt. Writing Poscar with no predictor corrector data.",
                    BadPoscarWarning,
                )

        return "\n".join(lines) + "\n"

    get_string = get_str

    def write_file(self, filename: PathLike, **kwargs) -> None:
        """Write POSCAR to a file. The supported kwargs are the same as those for
        the Poscar.get_str method and are passed through directly.
        """
        with zopen(filename, mode="wt") as file:
            file.write(self.get_str(**kwargs))

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
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            Poscar
        """
        return cls(
            Structure.from_dict(dct["structure"]),
            comment=dct["comment"],
            selective_dynamics=dct["selective_dynamics"],
            true_names=dct["true_names"],
            velocities=dct.get("velocities"),
            predictor_corrector=dct.get("predictor_corrector"),
        )

    def set_temperature(self, temperature: float) -> None:
        """
        Initialize the velocities based on Maxwell-Boltzmann distribution.
        Removes linear, but not angular drift (same as VASP).

        Scale the energies to the exact temperature (microcanonical ensemble)
        Velocities are given in A/fs. This is the VASP default when
        direct/cartesian is not specified (even when positions are given in
        direct coordinates).

        Overwrite imported velocities, if any.

        Args:
            temperature (float): Temperature in Kelvin.
        """
        # mean 0 variance 1
        velocities = np.random.default_rng().standard_normal((len(self.structure), 3))

        # In AMU, (N, 1) array
        atomic_masses = np.array([site.specie.atomic_mass.to("kg") for site in self.structure])
        dof = 3 * len(self.structure) - 3

        # Remove linear drift (net momentum)
        velocities -= np.mean(atomic_masses[:, None] * velocities, axis=0) / np.mean(atomic_masses)

        # Scale velocities due to atomic masses
        # mean 0 std proportional to sqrt(1/m)
        velocities /= atomic_masses[:, None] ** (1 / 2)

        # Scale velocities to get correct temperature
        energy = np.sum(1 / 2 * atomic_masses * np.sum(velocities**2, axis=1))
        scale = (temperature * dof / (2 * energy / const.k)) ** (1 / 2)

        velocities *= scale * 1e-5  # these are in A/fs

        self.temperature = temperature
        self.structure.site_properties.pop("selective_dynamics", None)
        self.structure.site_properties.pop("predictor_corrector", None)

        # Set as list[list] to be consistent with the other
        # initializations
        self.structure.add_site_property("velocities", velocities.tolist())


class BadPoscarWarning(UserWarning):
    """Warning class for bad POSCAR entries."""


class Incar(UserDict, MSONable):
    """
    A case-insensitive dictionary to read/write INCAR files with additional helper functions.

    - Keys are stored in uppercase to allow case-insensitive access (set, get, del, update, setdefault).
    - String values are capitalized by default, except for keys specified
        in the `lower_str_keys` of the `proc_val` method.
    """

    def __init__(self, params: dict[str, Any] | None = None) -> None:
        """
        Clean up params and create an Incar object.

        Args:
            params (dict): INCAR parameters as a dictionary.

        Warnings:
            BadIncarWarning: If there are duplicate in keys (case insensitive).
        """
        params = params or {}

        # Check for case-insensitive duplicate keys
        key_counter = Counter(key.strip().upper() for key in params)
        if duplicates := [key for key, count in key_counter.items() if count > 1]:
            warnings.warn(
                f"Duplicate keys found (case-insensitive): {duplicates}",
                BadIncarWarning,
                stacklevel=2,
            )

        # If INCAR contains vector-like MAGMOMS given as a list
        # of floats, convert to a list of lists
        if (params.get("MAGMOM") and isinstance(params["MAGMOM"][0], int | float)) and (
            params.get("LSORBIT") or params.get("LNONCOLLINEAR")
        ):
            val: list[list] = []
            for idx in range(len(params["MAGMOM"]) // 3):
                val.append(params["MAGMOM"][idx * 3 : (idx + 1) * 3])
            params["MAGMOM"] = val

        super().__init__(params)

    def __setitem__(self, key: str, val: Any) -> None:
        """
        Add parameter-val pair to Incar.
        - Clean the parameter and val by stripping leading
            and trailing white spaces.
        - Cast keys to upper case.
        """
        key = key.strip().upper()
        # Cast float/int to str such that proc_val would clean up their types
        val = self.proc_val(key, str(val)) if isinstance(val, str | float | int) else val
        super().__setitem__(key, val)

    def __getitem__(self, key: str) -> Any:
        """
        Get value using a case-insensitive key.
        """
        return super().__getitem__(key.strip().upper())

    def __delitem__(self, key: str) -> None:
        super().__delitem__(key.strip().upper())

    def __contains__(self, key: str) -> bool:
        return super().__contains__(key.upper().strip())

    def __str__(self) -> str:
        return self.get_str(sort_keys=True, pretty=False)

    def __add__(self, other: Self) -> Self:
        """
        Add all the values of another INCAR object to this object.
        Facilitate the use of "standard" INCARs.
        """
        params: dict[str, Any] = dict(self.items())
        for key, val in other.items():
            if key in self and val != self[key]:
                raise ValueError(f"INCARs have conflicting values for {key}: {self[key]} != {val}")
            params[key] = val
        return type(self)(params)

    def get(self, key: str, default: Any = None) -> Any:
        """
        Get a value for a case-insensitive key, return default if not found.
        """
        return super().get(key.strip().upper(), default)

    def as_dict(self) -> dict:
        """MSONable dict."""
        dct = dict(self)
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        return dct

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        """
        Args:
            dct (dict): Serialized Incar.

        Returns:
            Incar
        """
        if dct.get("MAGMOM") and isinstance(dct["MAGMOM"][0], dict):
            dct["MAGMOM"] = [Magmom.from_dict(m) for m in dct["MAGMOM"]]
        return cls({k: v for k, v in dct.items() if k not in ("@module", "@class")})

    def copy(self) -> Self:
        return type(self)(self)

    def get_str(self, sort_keys: bool = False, pretty: bool = False) -> str:
        """Get a string representation of the INCAR. Differ from the
        __str__ method to provide options for pretty printing.

        Args:
            sort_keys (bool): Whether to sort the INCAR parameters
                alphabetically. Defaults to False.
            pretty (bool): Whether to pretty align output.
                Defaults to False.
        """
        keys: list[str] = sorted(self) if sort_keys else list(self)
        lines = []
        for key in keys:
            if key == "MAGMOM" and isinstance(self[key], list):
                value = []

                if isinstance(self[key][0], list | Magmom) and (self.get("LSORBIT") or self.get("LNONCOLLINEAR")):
                    value.append(" ".join(str(i) for j in self[key] for i in j))
                elif self.get("LSORBIT") or self.get("LNONCOLLINEAR"):
                    for _key, group in itertools.groupby(self[key]):
                        value.append(f"3*{len(tuple(group))}*{_key}")
                else:
                    # float() to ensure backwards compatibility between
                    # float magmoms and Magmom objects
                    for _key, group in itertools.groupby(self[key], key=float):
                        value.append(f"{len(tuple(group))}*{_key}")

                lines.append([key, " ".join(value)])
            elif isinstance(self[key], list):
                lines.append([key, " ".join(map(str, self[key]))])
            else:
                lines.append([key, self[key]])

        if pretty:
            return str(tabulate([[line[0], "=", line[1]] for line in lines], tablefmt="plain"))
        return str_delimited(lines, None, " = ") + "\n"

    def write_file(self, filename: PathLike) -> None:
        """Write Incar to a file.

        Args:
            filename (str): filename to write to.
        """
        with zopen(filename, mode="wt") as file:
            file.write(str(self))

    @classmethod
    def from_file(cls, filename: PathLike) -> Self:
        """Read an Incar object from a file.

        Args:
            filename (str): Filename for file

        Returns:
            Incar object
        """
        with zopen(filename, mode="rt") as file:
            return cls.from_str(file.read())

    @classmethod
    def from_str(cls, string: str) -> Self:
        """Read an Incar object from a string.

        Args:
            string (str): Incar string

        Returns:
            Incar object
        """
        params: dict[str, Any] = {}
        for line in clean_lines(string.splitlines()):
            for sline in line.split(";"):
                if match := re.match(r"(\w+)\s*=\s*(.*)", sline.strip()):
                    key: str = match[1].strip()
                    val: str = match[2].strip()
                    params[key] = cls.proc_val(key, val)
        return cls(params)

    @staticmethod
    def proc_val(key: str, val: str) -> list | bool | float | int | str:
        """Helper method to convert INCAR parameters to proper types
        like ints, floats, lists, etc.

        Args:
            key (str): INCAR parameter key.
            val (str): Value of INCAR parameter.
        """
        list_keys = (
            "LDAUU",
            "LDAUL",
            "LDAUJ",
            "MAGMOM",
            "DIPOL",
            "LANGEVIN_GAMMA",
            "QUAD_EFG",
            "EINT",
        )
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
        float_keys = (
            "EDIFF",
            "SIGMA",
            "TIME",
            "ENCUTFOCK",
            "HFSCREEN",
            "POTIM",
            "EDIFFG",
            "AGGAC",
            "PARAM1",
            "PARAM2",
            "ENCUT",
        )
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
            "NSIM",
            "NKRED",
            "NUPDOWN",
            "ISPIND",
            "LDAUTYPE",
            "IVDW",
        )
        lower_str_keys = ("ML_MODE",)

        def smart_int_or_float(num_str: str) -> float:
            """Determine whether a string represents an integer or a float."""
            if "." in num_str or "e" in num_str.lower():
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
                if match := re.match(r"^\.?([T|F|t|f])[A-Za-z]*\.?", val):
                    return match[1].lower() == "t"

                raise ValueError(f"{key} should be a boolean type!")

            if key in float_keys:
                return float(re.search(r"^-?\d*\.?\d*[e|E]?-?\d*", val)[0])  # type: ignore[index]

            if key in int_keys:
                return int(re.match(r"^-?[0-9]+", val)[0])  # type: ignore[index]

            if key in lower_str_keys:
                return val.strip().lower()

        except ValueError:
            pass

        # Not in known keys. We will try a hierarchy of conversions.
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

    def diff(self, other: Self) -> dict[str, dict[str, Any]]:
        """
        Diff function for Incar. Compare two Incars and indicate which
        parameters are the same and which are not. Useful for checking whether
        two runs were done using the same parameters.

        Args:
            other (Incar): The other Incar object to compare to.

        Returns:
            dict[str, dict]: of the following format:
                {"Same" : parameters_that_are_the_same, "Different": parameters_that_are_different}
                Note that the parameters are return as full dictionaries of values. E.g. {"ISIF":3}
        """
        similar_params = {}
        different_params = {}
        for k1, v1 in self.items():
            if k1 not in other:
                different_params[k1] = {"INCAR1": v1, "INCAR2": None}
            elif v1 != other[k1]:
                different_params[k1] = {"INCAR1": v1, "INCAR2": other[k1]}
            else:
                similar_params[k1] = v1

        for k2, v2 in other.items():
            if k2 not in similar_params and k2 not in different_params and k2 not in self:
                different_params[k2] = {"INCAR1": None, "INCAR2": v2}

        return {"Same": similar_params, "Different": different_params}

    def check_params(self) -> None:
        """Check INCAR for invalid tags or values.
        If a tag doesn't exist, calculation will still run, however VASP
        will ignore the tag and set it as default without letting you know.
        """
        # Load INCAR tag/value check reference file
        with open(os.path.join(MODULE_DIR, "incar_parameters.json"), encoding="utf-8") as json_file:
            incar_params = json.loads(json_file.read())

        for tag, val in self.items():
            # Check if the tag exists
            if tag not in incar_params:
                warnings.warn(
                    f"Cannot find {tag} in the list of INCAR tags",
                    BadIncarWarning,
                    stacklevel=2,
                )
                continue

            # Check value type
            param_type: str = incar_params[tag].get("type")
            allowed_values: list[Any] = incar_params[tag].get("values")

            if param_type is not None and not isinstance(val, eval(param_type)):  # noqa: S307
                warnings.warn(f"{tag}: {val} is not a {param_type}", BadIncarWarning, stacklevel=2)

            # Only check value when it's not None,
            # meaning there is recording for corresponding value
            if allowed_values is not None:
                # Note: param_type could be a Union type, e.g. "str | bool"
                if "str" in param_type:
                    allowed_values = [item.capitalize() if isinstance(item, str) else item for item in allowed_values]

                if val not in allowed_values:
                    warnings.warn(
                        f"{tag}: Cannot find {val} in the list of values",
                        BadIncarWarning,
                        stacklevel=2,
                    )


class BadIncarWarning(UserWarning):
    """Warning class for bad INCAR parameters."""


@unique
class KpointsSupportedModes(Enum):
    """Enum type of all supported modes for Kpoint generation."""

    Automatic = 0
    Gamma = 1
    Monkhorst = 2
    Line_mode = 3
    Cartesian = 4
    Reciprocal = 5

    def __str__(self) -> str:
        return str(self.name)

    @classmethod
    def from_str(cls, mode: str) -> Self:
        """
        Args:
            mode: String.

        Returns:
            Kpoints_supported_modes
        """
        initial = mode.lower()[0]
        for key in cls:
            if key.name.lower()[0] == initial:
                return key
        raise ValueError(f"Invalid Kpoint {mode=}")


class Kpoints(MSONable):
    """KPOINTS reader/writer."""

    supported_modes: ClassVar[type[KpointsSupportedModes]] = KpointsSupportedModes

    def __init__(
        self,
        comment: str = "Default gamma",
        num_kpts: int = 0,
        style: KpointsSupportedModes = supported_modes.Gamma,
        kpts: Sequence[Kpoint] = ((1, 1, 1),),
        kpts_shift: Vector3D = (0, 0, 0),
        kpts_weights: list[float] | None = None,
        coord_type: Literal["Reciprocal", "Cartesian"] | None = None,
        labels: list[str] | None = None,
        tet_number: int = 0,
        tet_weight: float = 0,
        tet_connections: list[tuple] | None = None,
    ) -> None:
        """
        Highly flexible constructor for Kpoints object. The flexibility comes
        at the cost of usability and in general, it is recommended that you use
        the default constructor only if you know exactly what you are doing and
        requires the flexibility. For most usage cases, the three automatic
        schemes can be constructed far more easily using the convenience static
        constructors (automatic, gamma_automatic, monkhorst_automatic) and it
        is recommended that you use those.

        The default behavior of the constructor is for a Gamma-centered,
        1x1x1 KPOINTS with no shift.

        Args:
            comment (str): String comment for Kpoints. Defaults to "Default gamma".
            num_kpts: Following VASP method of defining the KPOINTS file, this
                parameter is the number of kpoints specified. If set to 0
                (or negative), VASP automatically generates the KPOINTS.
            style: Style for generating KPOINTS. Use one of the
                Kpoints.supported_modes enum types.
            kpts (2D array): Array of kpoints. Even when only a single
                specification is required, e.g. in the automatic scheme,
                the kpts should still be specified as a 2D array. e.g.
                [(20,),] or [(2, 2, 2),].
            kpts_shift (3x1 array): Shift for kpoints.
            kpts_weights (list[float]): Optional weights for explicit kpoints.
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
        """
        if num_kpts > 0 and not labels and not kpts_weights:
            raise ValueError("For explicit or line-mode kpoints, either the labels or kpts_weights must be specified.")

        self.comment = comment
        self.num_kpts = num_kpts
        self.kpts = kpts  # type: ignore[assignment]
        self.style = style
        self.coord_type = coord_type
        self.kpts_weights = kpts_weights
        self.kpts_shift = tuple(kpts_shift)
        self.labels = labels
        self.tet_number = tet_number
        self.tet_weight = tet_weight
        self.tet_connections = tet_connections

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.as_dict() == other.as_dict()

    def __repr__(self) -> str:
        lines = [self.comment, str(self.num_kpts), self.style.name]

        style = self.style.name.lower()[0]
        if style == "l":
            lines.append(str(self.coord_type))

        for idx, kpt in enumerate(self.kpts):
            lines.append(" ".join(map(str, kpt)))  # type: ignore[arg-type]
            if style == "l":
                if self.labels is not None:
                    lines[-1] += f" ! {self.labels[idx]}"
                if idx % 2 == 1:
                    lines[-1] += "\n"
            elif self.num_kpts > 0 and self.kpts_weights is not None:
                if self.labels is not None:
                    lines[-1] += f" {int(self.kpts_weights[idx])} {self.labels[idx]}"
                else:
                    lines[-1] += f" {int(self.kpts_weights[idx])}"

        # Tetrahedron parameters if the number of tetrahedrons > 0
        if style not in "lagm" and self.tet_number > 0:
            lines.extend(("Tetrahedron", f"{self.tet_number} {self.tet_weight:f}"))
            if self.tet_connections is not None:
                for sym_weight, vertices in self.tet_connections:
                    a, b, c, d = vertices
                    lines.append(f"{sym_weight} {a} {b} {c} {d}")

        # Shifts for automatic kpoints types if not zero
        if self.num_kpts <= 0 and tuple(self.kpts_shift) != (0, 0, 0):
            lines.append(" ".join(map(str, self.kpts_shift)))
        return "\n".join(lines) + "\n"

    @property
    def kpts(self) -> list[Kpoint]:
        """A sequence of Kpoints, where each Kpoint is a tuple of 3 or 1."""
        if all(isinstance(kpt, list | tuple | np.ndarray) and len(kpt) in {1, 3} for kpt in self._kpts):
            return list(map(tuple, self._kpts))  # type: ignore[arg-type]

        if all(isinstance(point, int | float) for point in self._kpts) and len(self._kpts) == 3:
            return [cast(Kpoint, tuple(self._kpts))]

        raise ValueError(f"Invalid Kpoint {self._kpts}.")

    @kpts.setter
    def kpts(self, kpts: Sequence[float | int] | Sequence[Sequence[float | int]]) -> None:
        """
        Args:
            kpts: Sequence[float | int] | Sequence[Sequence[float | int]].
        """
        self._kpts = kpts

    @property
    def style(self) -> KpointsSupportedModes:
        """Style for kpoint generation. One of Kpoints_supported_modes enum."""
        return self._style

    @style.setter
    def style(self, style: str | KpointsSupportedModes) -> None:
        """Set the style for the Kpoints. One of Kpoints_supported_modes enum.

        Args:
            style (str | KpointsSupportedModes): Style
        """
        if isinstance(style, str):
            style = type(self).supported_modes.from_str(style)

        if (
            style
            in {
                type(self).supported_modes.Automatic,
                type(self).supported_modes.Gamma,
                type(self).supported_modes.Monkhorst,
            }
            and len(self.kpts) > 1
        ):
            raise ValueError(
                "For fully automatic or automatic gamma or monk "
                "kpoints, only a single line for the number of "
                "divisions is allowed."
            )

        self._style = style

    @classmethod
    def automatic(cls, subdivisions: int) -> Self:
        """
        Constructor for a fully automatic Kpoint grid, with
        Gamma-centered grids and the number of subdivisions
        along each reciprocal lattice vector determined by the scheme in the
        VASP manual.

        Args:
            subdivisions (int): Number of subdivisions along
                each reciprocal lattice vector.

        Returns:
            Kpoints
        """
        warnings.warn("Please use INCAR KSPACING tag.", DeprecationWarning, stacklevel=2)

        return cls(
            "Fully automatic kpoint scheme",
            0,
            style=cls.supported_modes.Automatic,
            kpts=[
                (subdivisions,),
            ],
        )

    @classmethod
    def gamma_automatic(cls, kpts: Tuple3Ints = (1, 1, 1), shift: Vector3D = (0, 0, 0)) -> Self:
        """
        Construct an automatic Gamma-centered Kpoint grid.

        Args:
            kpts: Subdivisions N_1, N_2 and N_3 along reciprocal lattice
                vectors. Defaults to (1, 1, 1)
            shift: Shift to be applied to the kpoints. Defaults to (0, 0, 0).

        Returns:
            Kpoints
        """
        return cls(
            "Automatic kpoint scheme",
            0,
            cls.supported_modes.Gamma,
            kpts=[kpts],
            kpts_shift=shift,
        )

    @classmethod
    def monkhorst_automatic(cls, kpts: Tuple3Ints = (2, 2, 2), shift: Vector3D = (0, 0, 0)) -> Self:
        """
        Construct an automatic Monkhorst-Pack Kpoint grid.

        Args:
            kpts: Subdivisions N_1, N_2, N_3 along reciprocal lattice
                vectors. Defaults to (2, 2, 2)
            shift: Shift to be applied to the kpoints. Defaults to (0, 0, 0).

        Returns:
            Kpoints
        """
        return cls(
            "Automatic kpoint scheme",
            0,
            cls.supported_modes.Monkhorst,
            kpts=[kpts],
            kpts_shift=shift,
        )

    @classmethod
    def automatic_density(cls, structure: Structure, kppa: float, force_gamma: bool = False) -> Self:
        """Get an automatic Kpoints object based on a structure and a kpoint
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
        lengths: Vector3D = lattice.abc
        ngrid = kppa / len(structure)
        mult: float = (ngrid * lengths[0] * lengths[1] * lengths[2]) ** (1 / 3)

        num_div: Tuple3Ints = cast(Tuple3Ints, [math.floor(max(mult / length, 1)) for length in lengths])

        is_hexagonal: bool = lattice.is_hexagonal()
        is_face_centered: bool = structure.get_space_group_info()[0][0] == "F"
        has_odd: bool = any(idx % 2 == 1 for idx in num_div)
        if has_odd or is_hexagonal or is_face_centered or force_gamma:
            style = cls.supported_modes.Gamma
        else:
            style = cls.supported_modes.Monkhorst

        return cls(
            comment,
            0,
            style,
            [num_div],
            (0, 0, 0),
        )

    @classmethod
    def automatic_gamma_density(cls, structure: Structure, kppa: float) -> Self:
        """Get an automatic Kpoints object based on a structure and a kpoint
        density. Uses Gamma centered meshes always. For GW.

        Algorithm:
            Uses a simple approach scaling the number of divisions along each
            reciprocal lattice vector proportional to its length.

        Args:
            structure: Input structure
            kppa: Grid density
        """
        lattice = structure.lattice
        a, b, c = lattice.abc
        n_grid = kppa / len(structure)

        multip = (n_grid * a * b * c) ** (1 / 3)
        n_div: list[int] = cast(list[int], [round(multip / length) for length in lattice.abc])

        # Ensure that all num_div[i] > 0
        n_div = [idx if idx > 0 else 1 for idx in n_div]

        # VASP documentation recommends to use even grids for n <= 8 and odd grids for n > 8.
        n_div = [idx + idx % 2 if idx <= 8 else idx - idx % 2 + 1 for idx in n_div]

        style = cls.supported_modes.Gamma

        comment = f"pymatgen with grid density = {kppa:.0f} / number of atoms"

        n_kpts = 0
        return cls(
            comment,
            n_kpts,
            style,
            [cast(Tuple3Ints, tuple(n_div))],
            (0, 0, 0),
        )

    @classmethod
    def automatic_density_by_vol(cls, structure: Structure, kppvol: int, force_gamma: bool = False) -> Self:
        """Get an automatic Kpoints object based on a structure and a kpoint
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
        return cls.automatic_density(structure, kppa, force_gamma=force_gamma)

    @classmethod
    def automatic_density_by_lengths(
        cls,
        structure: Structure,
        length_densities: Sequence[float],
        force_gamma: bool = False,
    ) -> Self:
        """Get an automatic Kpoints object based on a structure and a k-point
        density normalized by lattice constants.

        Algorithm:
            For a given dimension, the number of k-points is chosen as
            length_density = # of kpoints * lattice constant, e.g. [50.0, 50.0, 1.0] would
            have k-points of 50/a x 50/b x 1/c.

        Args:
            structure (Structure): Input structure
            length_densities (list[float]): Defines the density of k-points in each
            dimension, e.g. [50.0, 50.0, 1.0].
            force_gamma (bool): Force a gamma centered mesh

        Returns:
            Kpoints
        """
        if len(length_densities) != 3:
            raise ValueError(f"The dimensions of length_densities must be 3, not {len(length_densities)}")

        comment: str = f"k-point density of {length_densities}/[a, b, c]"

        lattice = structure.lattice

        abc = lattice.abc
        num_div: Tuple3Ints = tuple(math.ceil(ld / abc[idx]) for idx, ld in enumerate(length_densities))

        is_hexagonal: bool = lattice.is_hexagonal()
        is_face_centered: bool = structure.get_space_group_info()[0][0] == "F"
        has_odd: bool = any(idx % 2 == 1 for idx in num_div)
        if has_odd or is_hexagonal or is_face_centered or force_gamma:
            style = cls.supported_modes.Gamma
        else:
            style = cls.supported_modes.Monkhorst

        return cls(
            comment,
            0,
            style,
            [
                num_div,
            ],
            (0, 0, 0),
        )

    @classmethod
    def automatic_linemode(cls, divisions: int, ibz: HighSymmKpath) -> Self:
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
                kpoints += [ibz.kpath["kpoints"][path[i]]] * 2
                labels += [path[i]] * 2

            kpoints.append(ibz.kpath["kpoints"][path[-1]])
            labels.append(path[-1])

        return cls(
            "Line_mode KPOINTS file",
            style=cls.supported_modes.Line_mode,
            coord_type="Reciprocal",
            kpts=kpoints,
            labels=labels,
            num_kpts=divisions,
        )

    def copy(self) -> Self:
        """Make a copy of the Kpoints object."""
        return self.from_dict(self.as_dict())

    @classmethod
    def from_file(cls, filename: PathLike) -> Self:
        """Read a Kpoints object from a KPOINTS file.

        Args:
            filename (PathLike): File to read.

        Returns:
            Kpoints object
        """
        with zopen(filename, mode="rt") as file:
            return cls.from_str(file.read())

    @classmethod
    def from_str(cls, string: str) -> Self:
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

        # Automatic Gamma-centered or Monkhorst-Pack KPOINTS, with optional shift
        if style in {"g", "m"}:
            _kpt: list[int] = [int(i) for i in lines[3].split()]
            if len(_kpt) != 3:
                raise ValueError("Invalid Kpoint length.")
            kpt: Tuple3Ints = cast(Tuple3Ints, tuple(_kpt))

            kpts_shift: Vector3D = (0, 0, 0)
            if len(lines) > 4 and coord_pattern.match(lines[4]):
                try:
                    _kpts_shift: list[float] = [float(i) for i in lines[4].split()]
                except ValueError:
                    _kpts_shift = [0, 0, 0]

                if len(_kpts_shift) != 3:
                    raise ValueError("Invalid kpoint shift length.")

                kpts_shift = cast(Vector3D, tuple(_kpts_shift))

            return cls.gamma_automatic(kpt, kpts_shift) if style == "g" else cls.monkhorst_automatic(kpt, kpts_shift)

        # Automatic kpoints with basis
        if num_kpts <= 0:
            _style = cls.supported_modes.Cartesian if style in "ck" else cls.supported_modes.Reciprocal
            _kpts_shift = [float(i) for i in lines[6].split()]
            kpts_shift = cast(Vector3D, tuple(_kpts_shift)) if len(_kpts_shift) == 3 else (0, 0, 0)

            kpts: list[Kpoint] = [
                cast(Kpoint, tuple(float(j) for j in lines[line_idx].split())) for line_idx in range(3, 6)
            ]

            return cls(
                comment=comment,
                num_kpts=num_kpts,
                style=_style,
                kpts=kpts,
                kpts_shift=kpts_shift,
            )

        # Line-mode KPOINTS, usually used with band structures
        if style == "l":
            coord_type: Literal["Cartesian", "Reciprocal"] = (
                "Cartesian" if lines[3].lower()[0] in "ck" else "Reciprocal"
            )
            _style = cls.supported_modes.Line_mode
            _kpts: list[Tuple3Floats] = []
            labels = []
            patt = re.compile(r"([e0-9.\-]+)\s+([e0-9.\-]+)\s+([e0-9.\-]+)\s*!*\s*(.*)")
            for idx in range(4, len(lines)):
                line = lines[idx]
                if match := patt.match(line):
                    _kpts.append((float(match[1]), float(match[2]), float(match[3])))
                    labels.append(match[4].strip())

            return cls(
                comment=comment,
                num_kpts=num_kpts,
                style=_style,
                kpts=_kpts,
                coord_type=coord_type,
                labels=labels,
            )

        # Assume explicit KPOINTS if all else fails.
        _style = cls.supported_modes.Cartesian if style in "ck" else cls.supported_modes.Reciprocal
        kpts = []
        kpts_weights = []
        labels = []
        tet_number = 0
        tet_weight: float = 0
        tet_connections = None

        for idx in range(3, 3 + num_kpts):
            tokens = lines[idx].split()
            kpts.append(cast(Tuple3Floats, tuple(float(i) for i in tokens[:3])))
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
                for idx in range(5 + num_kpts, 5 + num_kpts + tet_number):
                    tokens = lines[idx].split()
                    tet_connections.append((int(tokens[0]), [int(tokens[j]) for j in range(1, 5)]))
        except IndexError:
            pass

        return cls(
            comment=comment,
            num_kpts=num_kpts,
            style=cls.supported_modes[str(_style)],
            kpts=kpts,
            kpts_weights=kpts_weights,
            tet_number=tet_number,
            tet_weight=tet_weight,
            tet_connections=tet_connections,
            labels=labels,
        )

    def write_file(self, filename: PathLike) -> None:
        """Write Kpoints to a file.

        Args:
            filename (PathLike): Filename to write to.
        """
        with zopen(filename, mode="wt") as file:
            file.write(str(self))

    def as_dict(self) -> dict:
        """MSONable dict."""
        dct = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
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

        optional_paras = ("genvec1", "genvec2", "genvec3", "shift")
        for para in optional_paras:
            if para in self.__dict__:
                dct[para] = self.__dict__[para]
        return dct

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            Kpoints
        """
        comment = dct.get("comment", "")
        generation_style = cast(KpointsSupportedModes, dct.get("generation_style"))
        kpts = dct.get("kpoints", ((1, 1, 1),))
        kpts_shift = dct.get("usershift", (0, 0, 0))
        num_kpts = dct.get("nkpoints", 0)
        return cls(
            comment=comment,
            kpts=kpts,
            style=generation_style,
            kpts_shift=kpts_shift,
            num_kpts=num_kpts,
            kpts_weights=dct.get("kpts_weights"),
            coord_type=dct.get("coord_type"),
            labels=dct.get("labels"),
            tet_number=dct.get("tet_number", 0),
            tet_weight=dct.get("tet_weight", 0),
            tet_connections=dct.get("tet_connections"),
        )


def _parse_bool(string: str) -> bool:
    if match := re.match(r"^\.?([TFtf])[A-Za-z]*\.?", string):
        return match[1] in {"T", "t"}
    raise ValueError(f"{string} should be a boolean type!")


def _parse_float(string: str) -> float:
    if match := re.search(r"^-?\d*\.?\d*[eE]?-?\d*", string):
        return float(match[0])
    raise ValueError(f"{string} should be a float type!")


def _parse_int(string: str) -> int:
    if match := re.match(r"^-?[0-9]+", string):
        return int(match[0])
    raise ValueError(f"{string} should be an int type!")


def _parse_list(string: str) -> list[float]:
    """Parse a list of floats from a string."""
    return [float(y) for y in re.split(r"\s+", string.strip()) if not y.isalpha()]


class Orbital(NamedTuple):
    n: int
    l: int  # noqa: E741
    j: float
    E: float
    occ: float


class OrbitalDescription(NamedTuple):
    l: int  # noqa: E741
    E: float
    Type: int
    Rcut: float
    Type2: int | None
    Rcut2: float | None


# Hashes computed from the full POTCAR file contents by pymatgen (not 1st-party VASP hashes)
PYMATGEN_POTCAR_HASHES = loadfn(f"{MODULE_DIR}/vasp_potcar_pymatgen_hashes.json")
# Written to some newer POTCARs by VASP
VASP_POTCAR_HASHES = loadfn(f"{MODULE_DIR}/vasp_potcar_file_hashes.json")
POTCAR_STATS_PATH: str = os.path.join(MODULE_DIR, "potcar-summary-stats.json.bz2")


class PmgVaspPspDirError(ValueError):
    """Error thrown when PMG_VASP_PSP_DIR is not configured, but POTCAR is requested."""


class PotcarSingle:
    """
    Object for a **single** POTCAR. The builder assumes the POTCAR contains
    the complete untouched string "data" and a dict of keywords.

    Attributes:
        data (str): POTCAR data as a string.
        keywords (dict): Keywords parsed from the POTCAR as a dict. All keywords are also
            accessible as attributes in themselves. e.g. potcar.enmax, potcar.encut, etc.

    MD5 hashes of the entire POTCAR file and the actual data are validated
    against a database of known good hashes. Appropriate warnings or errors
    are raised if validation fails.
    """

    # Note: there are multiple releases of the {LDA,PBE} {52,54} POTCARs
    #     the original (UNIVIE) releases include no SHA256 hashes nor COPYR fields
    #     in the PSCTR/header field.
    # We indicate the older release in `functional_dir` as PBE_52, PBE_54, LDA_52, LDA_54.
    # The newer release is indicated as PBE_52_W_HASH, etc.
    functional_dir: ClassVar[dict[str, str]] = {
        "PBE": "POT_GGA_PAW_PBE",
        "PBE_52": "POT_GGA_PAW_PBE_52",
        "PBE_52_W_HASH": "POTPAW_PBE_52",
        "PBE_54": "POT_GGA_PAW_PBE_54",
        "PBE_54_W_HASH": "POTPAW_PBE_54",
        "PBE_64": "POT_PAW_PBE_64",
        "LDA": "POT_LDA_PAW",
        "LDA_52": "POT_LDA_PAW_52",
        "LDA_52_W_HASH": "POTPAW_LDA_52",
        "LDA_54": "POT_LDA_PAW_54",
        "LDA_54_W_HASH": "POTPAW_LDA_54",
        "LDA_64": "POT_LDA_PAW_64",
        "PW91": "POT_GGA_PAW_PW91",
        "LDA_US": "POT_LDA_US",
        "PW91_US": "POT_GGA_US_PW91",
        "Perdew_Zunger81": "POT_LDA_PAW",
    }

    functional_tags: ClassVar[dict[str, dict[Literal["name", "class"], str]]] = {
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

    parse_functions: ClassVar[dict[str, Any]] = {
        "LULTRA": _parse_bool,
        "LUNSCR": _parse_bool,
        "LCOR": _parse_bool,
        "LPAW": _parse_bool,
        "EATOM": _parse_float,
        "RPACOR": _parse_float,
        "POMASS": _parse_float,
        "ZVAL": _parse_float,
        "RCORE": _parse_float,
        "RWIGS": _parse_float,
        "ENMAX": _parse_float,
        "ENMIN": _parse_float,
        "EMMIN": _parse_float,
        "EAUG": _parse_float,
        "DEXC": _parse_float,
        "RMAX": _parse_float,
        "RAUG": _parse_float,
        "RDEP": _parse_float,
        "RDEPT": _parse_float,
        "QCUT": _parse_float,
        "QGAM": _parse_float,
        "RCLOC": _parse_float,
        "IUNSCR": _parse_int,
        "ICORE": _parse_int,
        "NDATA": _parse_int,
        "VRHFIN": str.strip,
        "LEXCH": str.strip,
        "TITEL": str.strip,
        "STEP": _parse_list,
        "RRKJ": _parse_list,
        "GGA": _parse_list,
        "SHA256": str.strip,
        "COPYR": str.strip,
    }

    # Used for POTCAR validation
    _potcar_summary_stats = loadfn(POTCAR_STATS_PATH)

    def __init__(self, data: str, symbol: str | None = None) -> None:
        """
        Args:
            data (str): Complete, single and raw POTCAR file as a string.
            symbol (str): POTCAR symbol corresponding to the filename suffix
                e.g. "Tm_3" for POTCAR.TM_3".
                If not given, pymatgen will attempt to extract the symbol
                from the file itself, but is not always reliable!
        """
        self.data = data

        # VASP parses header in vasprun.xml and this differs from the TITEL
        self.header = data.split("\n")[0].strip()

        match = re.search(r"(?s)(parameters from PSCTR are:.*?END of PSCTR-controll parameters)", data)
        search_lines = match[1] if match else ""

        keywords = {}
        for key, val in re.findall(r"(\S+)\s*=\s*(.*?)(?=;|$)", search_lines, flags=re.MULTILINE):
            try:
                keywords[key] = self.parse_functions[key](val)  # type: ignore[operator]
            except KeyError:
                warnings.warn(f"Ignoring unknown variable type {key}")

        PSCTR: dict[str, Any] = {}

        array_search = re.compile(r"(-*[0-9.]+)")
        orbitals: list[Orbital] = []
        descriptions: list[OrbitalDescription] = []
        if atomic_config_match := re.search(r"(?s)Atomic configuration(.*?)Description", search_lines):
            lines = atomic_config_match[1].splitlines()
            match = re.search(r"([0-9]+)", lines[1])
            num_entries = int(match[1]) if match else 0
            PSCTR["nentries"] = num_entries
            for line in lines[3:]:
                if orbit := array_search.findall(line):
                    orbitals.append(
                        Orbital(
                            int(orbit[0]),
                            int(orbit[1]),
                            float(orbit[2]),
                            float(orbit[3]),
                            float(orbit[4]),
                        )
                    )
            PSCTR["Orbitals"] = tuple(orbitals)

        if description_string := re.search(
            r"(?s)Description\s*\n(.*?)Error from kinetic energy argument \(eV\)",
            search_lines,
        ):
            for line in description_string[1].splitlines():
                if description := array_search.findall(line):
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
            for line in rrkj_kinetic_energy_string[1].splitlines():
                if "=" not in line:
                    rrkj_array += _parse_list(line.strip("\n"))
            if rrkj_array:
                PSCTR["RRKJ"] = tuple(rrkj_array)

        self.keywords = dict(sorted((PSCTR | keywords).items()))

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

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.data == other.data and self.keywords == other.keywords

    def __getattr__(self, attr: str) -> Any:
        """Delegates attributes to keywords. For example, you can use potcarsingle.enmax to get the ENMAX of the POTCAR.

        For float type properties, they are converted to the correct float. By
        default, all energies in eV and all length scales are in Angstroms.
        """
        try:
            return self.keywords[attr.upper()]
        except Exception as exc:
            raise AttributeError(attr) from exc

    def __str__(self) -> str:
        return f"{self.data}\n"

    def __repr__(self) -> str:
        cls_name = type(self).__name__
        symbol, functional = self.symbol, self.functional
        TITEL, VRHFIN, n_valence_elec = (self.keywords.get(key) for key in ("TITEL", "VRHFIN", "ZVAL"))
        return f"{cls_name}({symbol=}, {functional=}, {TITEL=}, {VRHFIN=}, {n_valence_elec=:.0f})"

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

    @property
    def element(self) -> str:
        """Attempt to return the atomic symbol based on the VRHFIN keyword."""
        element = self.keywords["VRHFIN"].split(":")[0].strip()
        try:
            return Element(element).symbol

        except ValueError:
            # VASP incorrectly gives the element symbol for Xe as "X"
            # Some potentials, e.g. Zr_sv, gives the symbol as r.
            if element == "X":
                return "Xe"
            return Element(self.symbol.split("_")[0]).symbol

    @property
    def atomic_no(self) -> int:
        """Attempt to return the atomic number based on the VRHFIN keyword."""
        return Element(self.element).Z

    @property
    def nelectrons(self) -> float:
        """Number of electrons."""
        return self.zval

    @property
    def symbol(self) -> str:
        """The POTCAR symbol, e.g. W_pv."""
        return self._symbol

    @property
    def potential_type(self) -> Literal["NC", "PAW", "US"]:
        """Type of PSP: NC (Norm-conserving), US (Ultra-soft), PAW (Projector augmented wave)."""
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
    def functional_class(self) -> str | None:
        """Functional class associated with PotcarSingle."""
        return self.functional_tags.get(self.LEXCH.lower(), {}).get("class")

    @property
    def hash_sha256_from_file(self) -> str | None:
        """SHA256 hash of the POTCAR file as read from the file. None if no SHA256 hash is found."""
        if _sha256 := getattr(self, "SHA256", None):
            return _sha256.split()[0]
        return None

    @property
    def sha256_computed_file_hash(self) -> str:
        """Compute a SHA256 hash of the PotcarSingle EXCLUDING lines starting with 'SHA256' and 'COPYR'."""
        # We have to remove lines with the hash itself and the copyright
        # notice to get the correct hash.
        potcar_list = self.data.split("\n")
        potcar_to_hash = [line for line in potcar_list if not line.strip().startswith(("SHA256", "COPYR"))]
        potcar_to_hash_str = "\n".join(potcar_to_hash)
        return sha256(potcar_to_hash_str.encode("utf-8")).hexdigest()

    @property
    def md5_computed_file_hash(self) -> str:
        """MD5 hash of the entire PotcarSingle."""
        # usedforsecurity=False needed in FIPS mode (Federal Information Processing Standards)
        # https://github.com/materialsproject/pymatgen/issues/2804
        md5 = hashlib.md5(usedforsecurity=False)
        md5.update(self.data.encode("utf-8"))
        return md5.hexdigest()

    @property
    def md5_header_hash(self) -> str:
        """MD5 hash of the metadata defining the PotcarSingle."""
        hash_str = ""
        for k, v in self.keywords.items():
            # For newer POTCARS we have to exclude 'SHA256' and 'COPYR lines
            # since they were not used in the initial hashing
            if k in {"nentries", "Orbitals", "SHA256", "COPYR"}:
                continue
            hash_str += f"{k}"
            if isinstance(v, bool | int):
                hash_str += f"{v}"
            elif isinstance(v, float):
                hash_str += f"{v:.3f}"
            elif isinstance(v, tuple | list):
                for item in v:
                    if isinstance(item, float):
                        hash_str += f"{item:.3f}"
                    elif isinstance(item, Orbital | OrbitalDescription):
                        for item_v in item:
                            if isinstance(item_v, int | str):
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
        md5 = hashlib.md5(usedforsecurity=False)
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
            }.

        Rationale:
            Each POTCAR is structured as
                Header (self.keywords)
                Data (actual pseudopotential values in data blocks)

            For the Data block of POTCAR, there are unformatted data blocks
            of unknown length and contents/data type, e.g. you might see
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
        # Some POTCARs have an LEXCH (functional used to generate the POTCAR)
        # with the expected functional, e.g. the C_d POTCAR for PBE is actually an
        # LDA pseudopotential.

        # Thus we have to look for matches in all POTCAR dirs, not just the ones with
        # consistent values of LEXCH
        for func in self.functional_dir:
            for titel_no_spc in self._potcar_summary_stats[func]:
                if self.TITEL.replace(" ", "") == titel_no_spc:
                    for potcar_subvariant in self._potcar_summary_stats[func][titel_no_spc]:
                        if self.VRHFIN.replace(" ", "") == potcar_subvariant["VRHFIN"]:
                            possible_match = {
                                "POTCAR_FUNCTIONAL": func,
                                "TITEL": titel_no_spc,
                                **potcar_subvariant,
                            }
                            possible_potcar_matches.append(possible_match)

        def parse_fortran_style_str(input_str: str) -> str | bool | float | int:
            """Parse any input string as bool, int, float, or failing that, str.
            Used to parse FORTRAN-generated POTCAR files where it's unknown
            a priori what type of data will be encountered.
            """
            input_str = input_str.strip()

            if input_str.lower() in {"t", "f", "true", "false"}:
                return input_str[0].lower() == "t"

            if input_str.upper() == input_str.lower() and input_str[0].isnumeric():
                # NB: fortran style floats always include a decimal point.
                #     While you can set, e.g. x = 1E4, you cannot print/write x without
                #     a decimal point:
                #         `write(6,*) x`          -->   `10000.0000` in stdout
                #         `write(6,'(E10.0)') x`  -->   segfault
                #     The (E10.0) means write an exponential-format number with 10
                #         characters before the decimal, and 0 characters after
                return float(input_str) if "." in input_str else int(input_str)

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
                elif isinstance(parsed_val, float | int):
                    psp_vals.append(parsed_val)
            if len(tmp_str) > 0:
                psp_keys.append(tmp_str.lower())

        keyword_vals = []
        for kwd in self.keywords:
            val = self.keywords[kwd]
            if isinstance(val, bool):
                # has to come first since bools are also ints
                keyword_vals.append(1.0 if val else 0.0)
            elif isinstance(val, float | int):
                keyword_vals.append(val)
            elif hasattr(val, "__len__"):
                keyword_vals += [num for num in val if isinstance(num, float | int)]

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
        self._summary_stats: dict[str, dict] = {  # for this PotcarSingle instance
            "keywords": {
                "header": [kwd.lower() for kwd in self.keywords],
                "data": psp_keys,
            },
            "stats": {
                "header": data_stats(keyword_vals),
                "data": data_stats(psp_vals),
            },
        }

        data_match_tol: float = 1e-6
        for ref_psp in possible_potcar_matches:
            key_match = all(
                set(ref_psp["keywords"][key]) == set(self._summary_stats["keywords"][key]) for key in ["header", "data"]
            )

            data_diff = [
                abs(ref_psp["stats"][key][stat] - self._summary_stats["stats"][key][stat])
                for stat in ["MEAN", "ABSMEAN", "VAR", "MIN", "MAX"]
                for key in ["header", "data"]
            ]
            data_match = all(np.array(data_diff) < data_match_tol)

            if key_match and data_match:
                return True

        return False

    def write_file(self, filename: str) -> None:
        """Write PotcarSingle to a file.

        Args:
            filename (str): Filename to write to.
        """
        with zopen(filename, mode="wt") as file:
            file.write(str(self))

    def copy(self) -> Self:
        """Return a copy of the PotcarSingle.

        Returns:
            PotcarSingle
        """
        return type(self)(self.data, symbol=self.symbol)

    @classmethod
    def from_file(cls, filename: PathLike) -> Self:
        """Read PotcarSingle from file.

        Args:
            filename: Filename.

        Returns:
            PotcarSingle
        """
        match = re.search(r"(?<=POTCAR\.)(.*)(?=.gz)", str(filename))
        symbol = match[0] if match else ""

        try:
            with zopen(filename, mode="rt") as file:
                return cls(file.read(), symbol=symbol or None)

        except UnicodeDecodeError:
            warnings.warn("POTCAR contains invalid unicode errors. We will attempt to read it by ignoring errors.")

            with codecs.open(str(filename), "r", encoding="utf-8", errors="ignore") as file:
                return cls(file.read(), symbol=symbol or None)

    @classmethod
    def from_symbol_and_functional(
        cls,
        symbol: str,
        functional: str | None = None,
    ) -> Self:
        """Make a PotcarSingle from a symbol and functional.

        Args:
            symbol (str): Symbol, e.g. Li_sv
            functional (str): Functional, e.g. PBE

        Returns:
            PotcarSingle
        """
        functional = functional or SETTINGS.get("PMG_DEFAULT_FUNCTIONAL", "PBE")
        if functional is None:
            raise ValueError("Cannot get functional.")

        functional_subdir = SETTINGS.get("PMG_VASP_PSP_SUB_DIRS", {}).get(functional, cls.functional_dir[functional])
        PMG_VASP_PSP_DIR = SETTINGS.get("PMG_VASP_PSP_DIR")
        if PMG_VASP_PSP_DIR is None:
            raise PmgVaspPspDirError("Set PMG_VASP_PSP_DIR=<directory-path> in .pmgrc.yaml (needed to find POTCARs)")
        if not os.path.isdir(PMG_VASP_PSP_DIR):
            raise FileNotFoundError(f"{PMG_VASP_PSP_DIR=} does not exist.")

        paths_to_try: list[str] = [
            os.path.join(PMG_VASP_PSP_DIR, functional_subdir, f"POTCAR.{symbol}"),
            os.path.join(PMG_VASP_PSP_DIR, functional_subdir, symbol, "POTCAR"),
        ]
        path = paths_to_try[0]
        for path in paths_to_try:
            path = os.path.expanduser(path)
            path = zpath(path)
            if os.path.isfile(path):
                return cls.from_file(path)

        raise FileNotFoundError(
            f"You do not have the right POTCAR with {functional=} and {symbol=}\n"
            f"in your {PMG_VASP_PSP_DIR=}.\nPaths tried:\n- " + "\n- ".join(paths_to_try)
        )

    def verify_potcar(self) -> tuple[bool, bool]:
        """
        Attempt to verify the integrity of the POTCAR data.

        This method checks the whole file (removing only the SHA256
        metadata) against the SHA256 hash in the header if this is found.
        If no SHA256 hash is found in the file, the file hash (md5 hash of the
        whole file) is checked against all POTCAR file hashes known to pymatgen.

        Returns:
            tuple[bool, bool]: has_sha256 and passed_hash_check.
        """
        if self.hash_sha256_from_file:
            has_sha256 = True
            hash_is_valid = self.hash_sha256_from_file == self.sha256_computed_file_hash

        else:
            has_sha256 = False
            # If no sha256 hash is found in the POTCAR file, compare the whole
            # file with known potcar file hashes.
            md5_file_hash = self.md5_computed_file_hash
            hash_is_valid = md5_file_hash in VASP_POTCAR_HASHES

        return has_sha256, hash_is_valid

    def identify_potcar(
        self,
        mode: Literal["data", "file"] = "data",
        data_tol: float = 1e-6,
    ) -> tuple[list[str], list[str]]:
        """
        Identify the symbol and compatible functionals associated with this PotcarSingle.

        This method checks the summary statistics of either the POTCAR metadadata
        (PotcarSingle._summary_stats[key]["header"] for key in ("keywords", "stats") )
        or the entire POTCAR file (PotcarSingle._summary_stats) against a database
        of hashes for POTCARs distributed with VASP 5.4.4.

        Args:
            mode ("data" or "file"): "data" mode checks the POTCAR header keywords
                and stats only while "file" mode checks the entire summary stats.
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
            for ref_psp in self._potcar_summary_stats[func].get(self.TITEL.replace(" ", ""), []):
                if self.VRHFIN.replace(" ", "") != ref_psp["VRHFIN"]:
                    continue

                key_match = all(
                    set(ref_psp["keywords"][key]) == set(self._summary_stats["keywords"][key]) for key in check_modes
                )

                data_diff = [
                    abs(ref_psp["stats"][key][stat] - self._summary_stats["stats"][key][stat])
                    for stat in ["MEAN", "ABSMEAN", "VAR", "MIN", "MAX"]
                    for key in check_modes
                ]

                data_match = all(np.array(data_diff) < data_tol)

                if key_match and data_match:
                    identity["potcar_functionals"].append(func)
                    identity["potcar_symbols"].append(ref_psp["symbol"])

        for key, values in identity.items():
            if len(values) == 0:
                # The two keys are set simultaneously, either key being zero indicates no match
                return [], []
            identity[key] = list(set(values))

        return identity["potcar_functionals"], identity["potcar_symbols"]

    def identify_potcar_hash_based(
        self,
        mode: Literal["data", "file"] = "data",
    ) -> tuple[list[str], list[str]]:
        """
        Identify the symbol and compatible functionals associated with this PotcarSingle.

        This method checks the MD5 hash of either the POTCAR metadadata (PotcarSingle.md5_header_hash)
        or the entire POTCAR file (PotcarSingle.md5_computed_file_hash) against a database
        of hashes for POTCARs distributed with VASP 5.4.4.

        Args:
            mode ("data" | "file"): "data" mode checks the hash of the POTCAR metadata in self.keywords,
                while "file" mode checks the hash of the entire POTCAR file.

        Returns:
            symbol (list): List of symbols associated with the PotcarSingle
            potcar_functionals (list): List of potcar functionals associated with
                the PotcarSingle
        """
        # Dict to translate the sets in the .json file to the keys used in VaspInputSet
        mapping_dict: dict[str, dict[str, str]] = {
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

        if identity := hash_db.get(potcar_hash):
            # Convert the potcar_functionals from the .json dict into the functional
            # keys that pymatgen uses
            potcar_functionals = [*{mapping_dict[i]["pymatgen_key"] for i in identity["potcar_functionals"]}]

            return potcar_functionals, identity["potcar_symbols"]
        return [], []


def _gen_potcar_summary_stats(
    append: bool = False,
    vasp_psp_dir: str | None = None,
    summary_stats_filename: str | None = POTCAR_STATS_PATH,
) -> dict:
    """
    Regenerate the reference data in potcar-summary-stats.json.bz2 used to validate POTCARs
    by comparing header values and several statistics of copyrighted POTCAR data without
    having to record the POTCAR data itself.

    THIS FUNCTION IS DESTRUCTIVE. It will completely overwrite potcar-summary-stats.json.bz2.

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
    for func, func_dir in PotcarSingle.functional_dir.items():
        cpsp_dir = f"{vasp_psp_dir}/{func_dir}"
        if os.path.isdir(cpsp_dir):
            func_dir_exist[func] = func_dir
        else:
            warnings.warn(f"missing {func_dir} POTCAR directory")

    # Use append = True if a new POTCAR library is released to add new summary stats
    # without completely regenerating the dict of summary stats
    # Use append = False to completely regenerate the summary stats dict
    new_summary_stats = loadfn(summary_stats_filename) if append else {}

    for func, func_dir in func_dir_exist.items():
        new_summary_stats.setdefault(func, {})  # initialize dict if key missing

        potcar_list = glob(f"{vasp_psp_dir}/{func_dir}/POTCAR*") + glob(f"{vasp_psp_dir}/{func_dir}/*/POTCAR*")
        for potcar in potcar_list:
            psp = PotcarSingle.from_file(potcar)
            titel_key = psp.TITEL.replace(" ", "")

            # Some POTCARs have the same TITEL, but are named differently
            # e.g. there is an "original" PBE POTCAR.Fe_pv and a POTCAR.Fe_pv_new
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

    if summary_stats_filename:
        dumpfn(new_summary_stats, summary_stats_filename)

    return new_summary_stats


class Potcar(list, MSONable):
    """Read and write POTCAR files for calculations. Consists of a list of PotcarSingle."""

    FUNCTIONAL_CHOICES: ClassVar[tuple] = tuple(PotcarSingle.functional_dir)

    def __init__(
        self,
        symbols: Sequence[str] | None = None,
        functional: str | None = None,
        sym_potcar_map: dict[str, str] | None = None,
    ) -> None:
        """
        Args:
            symbols (list[str]): Element symbols for POTCAR. This should correspond
                to the symbols used by VASP. e.g. "Mg", "Fe_pv", etc.
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

    def __str__(self) -> str:
        return "\n".join(str(potcar).strip("\n") for potcar in self) + "\n"

    @property
    def symbols(self) -> list[str]:
        """The atomic symbols of all the atoms in the POTCAR file."""
        return [psingle.symbol for psingle in self]

    @symbols.setter
    def symbols(self, symbols: Sequence[str]) -> None:
        self.set_symbols(symbols, functional=self.functional)

    @property
    def spec(self) -> list[dict]:
        """The atomic symbols and hash of all the atoms in the POTCAR file."""
        return [{"symbol": psingle.symbol, "hash": psingle.md5_computed_file_hash} for psingle in self]

    def as_dict(self) -> dict:
        """MSONable dict representation."""
        return {
            "functional": self.functional,
            "symbols": self.symbols,
            "@module": type(self).__module__,
            "@class": type(self).__name__,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            Potcar
        """
        return Potcar(symbols=dct["symbols"], functional=dct["functional"])

    @classmethod
    def from_file(cls, filename: PathLike) -> Self:
        """
        Reads Potcar from file.

        Args:
            filename: Filename

        Returns:
            Potcar
        """
        with zopen(filename, mode="rt") as file:
            fdata = file.read()

        potcar = cls()
        functionals: list[str | None] = []
        for psingle_str in fdata.split("End of Dataset"):
            if p_strip := psingle_str.strip():
                psingle = PotcarSingle(f"{p_strip}\nEnd of Dataset\n")
                potcar.append(psingle)
                functionals.append(psingle.functional)

        if len(set(functionals)) != 1:
            raise ValueError("File contains incompatible functionals!")

        potcar.functional = functionals[0]
        return potcar

    def write_file(self, filename: PathLike) -> None:
        """Write Potcar to a file.

        Args:
            filename (PathLike): filename to write to.
        """
        with zopen(filename, mode="wt") as file:
            file.write(str(self))

    def set_symbols(
        self,
        symbols: Sequence[str],
        functional: str | None = None,
        sym_potcar_map: dict[str, str] | None = None,
    ) -> None:
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


class UnknownPotcarWarning(UserWarning):
    """Warning raised when POTCAR hashes do not pass validation."""


class VaspInput(dict, MSONable):
    """Contain a set of VASP input objects corresponding to a run."""

    def __init__(
        self,
        incar: dict | Incar,
        kpoints: Kpoints | None,
        poscar: Poscar,
        potcar: Potcar | str | None,
        potcar_spec: bool = False,
        optional_files: dict[PathLike, object] | None = None,
        **kwargs,
    ) -> None:
        """
        Initialize a VaspInput object with the given input files.

        Args:
            incar (Incar): The Incar object.
            kpoints (Kpoints): The Kpoints object.
            poscar (Poscar): The Poscar object.
            potcar (Potcar or str): The Potcar object.
            potcar_spec (bool = False) : used to share POTCAR info without license issues.
                True --> POTCAR is a list of symbols, write POTCAR.spec
                False --> POTCAR is a VASP POTCAR, write POTCAR
            optional_files (dict): Other input files supplied as a dict of {filename: object}.
                The object should follow standard pymatgen conventions in implementing a
                as_dict() and from_dict method.
            **kwargs: Additional keyword arguments to be stored in the VaspInput object.
        """
        super().__init__(**kwargs)
        self._potcar_filename = "POTCAR" + (".spec" if potcar_spec else "")
        self.update(
            {
                "INCAR": Incar(incar),
                "KPOINTS": kpoints,
                "POSCAR": poscar,
                self._potcar_filename: potcar,
            }
        )
        if optional_files is not None:
            self.update(optional_files)

    def __str__(self) -> str:
        output: list = []
        for key, val in self.items():
            output.extend((key, str(val), ""))
        return "\n".join(output)

    def as_dict(self) -> dict:
        """MSONable dict."""
        dct = {key: val.as_dict() for key, val in self.items()}
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        return dct

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            VaspInput
        """
        sub_dct: dict[str, dict] = {"optional_files": {}}
        for key, val in dct.items():
            if key in ("INCAR", "POSCAR", "POTCAR", "KPOINTS"):
                sub_dct[key.lower()] = MontyDecoder().process_decoded(val)
            elif key not in ["@module", "@class"]:
                sub_dct["optional_files"][key] = MontyDecoder().process_decoded(val)
        return cls(**sub_dct)  # type: ignore[arg-type]

    def write_input(
        self,
        output_dir: PathLike = ".",
        make_dir_if_not_present: bool = True,
        cif_name: str | None = None,
        zip_name: str | None = None,
        files_to_transfer: dict | None = None,
    ) -> None:
        """Write VASP inputs to a directory.

        Args:
            output_dir (PathLike): Directory to write to.
                Defaults to current directory (".").
            make_dir_if_not_present (bool): Create the directory if not
                present. Defaults to True.
            cif_name (str or None): If a str, the name of the CIF file
                to write the POSCAR to (the POSCAR will also be written).
            zip_name (str or None): If a str, the name of the zip to
                archive the VASP input set to.
            files_to_transfer (dict) : A dictionary of
                    { < input filename >: < output filepath >}.
                This allows the transfer of < input filename > files from
                a previous calculation to < output filepath >.
        """
        if not os.path.isdir(output_dir) and make_dir_if_not_present:
            os.makedirs(output_dir)

        for key, value in self.items():
            if value is not None:
                with zopen(os.path.join(output_dir, key), mode="wt") as file:
                    file.write(str(value))

        if cif_name:
            self["POSCAR"].structure.to(filename=cif_name)

        if zip_name:
            files_to_zip = list(self) + ([cif_name] if cif_name else [])
            with ZipFile(os.path.join(output_dir, zip_name), mode="w") as zip_file:
                for file in files_to_zip:
                    try:
                        zip_file.write(os.path.join(output_dir, file), arcname=file)
                    except FileNotFoundError:
                        pass

                    try:
                        os.remove(os.path.join(output_dir, file))
                    except (FileNotFoundError, PermissionError, IsADirectoryError):
                        pass

        files_to_transfer = files_to_transfer or {}
        for key, val in files_to_transfer.items():
            with (
                zopen(val, "rb") as fin,
                zopen(str(Path(output_dir) / key), "wb") as fout,
            ):
                copyfileobj(fin, fout)

    @classmethod
    def from_directory(
        cls,
        input_dir: PathLike,
        optional_files: dict | None = None,
    ) -> Self:
        """
        Read in a set of VASP inputs from a directory. Note that only the
        standard INCAR, POSCAR, POTCAR and KPOINTS files are read unless
        optional_filenames is specified.

        Args:
            input_dir (PathLike): Directory to read VASP input from.
            optional_files (dict): Optional files to read in as a
                dict of {filename: Object type}. Objects must have
                from_file method.
        """
        sub_dct: dict[str, Any] = {}
        for fname, ftype in (
            ("INCAR", Incar),
            ("KPOINTS", Kpoints),
            ("POSCAR", Poscar),
            ("POTCAR", Potcar),
        ):
            try:
                full_zpath = zpath(os.path.join(input_dir, fname))
                sub_dct[fname.lower()] = ftype.from_file(full_zpath)  # type: ignore[attr-defined]
            except FileNotFoundError:  # handle the case where there is no KPOINTS file
                sub_dct[fname.lower()] = None  # type: ignore[assignment]

        sub_dct["optional_files"] = {  # type: ignore[assignment]
            fname: ftype.from_file(os.path.join(input_dir, fname)) for fname, ftype in (optional_files or {}).items()
        }

        return cls(**sub_dct)

    def copy(self, deep: bool = True) -> Self:
        """Deep copy of VaspInput."""
        if deep:
            return self.from_dict(self.as_dict())
        return type(self)(**{key.lower(): val for key, val in self.items()})

    def run_vasp(
        self,
        run_dir: PathLike = ".",
        vasp_cmd: list | None = None,
        output_file: PathLike = "vasp.out",
        err_file: PathLike = "vasp.err",
    ) -> None:
        """Write input files and run VASP.

        Args:
            run_dir: Where to write input files and do the run.
            vasp_cmd: Args to be supplied to run VASP. Otherwise, the
                PMG_VASP_EXE in .pmgrc.yaml is used.
            output_file: File to write output.
            err_file: File to write err.
        """
        self.write_input(output_dir=run_dir)
        vasp_cmd = vasp_cmd or SETTINGS.get("PMG_VASP_EXE")  # type: ignore[assignment]
        if not vasp_cmd:
            raise ValueError("No VASP executable specified!")
        vasp_cmd = [os.path.expanduser(os.path.expandvars(t)) for t in vasp_cmd]
        if not vasp_cmd:
            raise RuntimeError("You need to supply vasp_cmd or set the PMG_VASP_EXE in .pmgrc.yaml to run VASP.")

        with (
            cd(run_dir),
            open(output_file, mode="w", encoding="utf-8") as stdout_file,
            open(err_file, mode="w", encoding="utf-8", buffering=1) as stderr_file,
        ):
            subprocess.check_call(vasp_cmd, stdout=stdout_file, stderr=stderr_file)

    @property
    def incar(self) -> Incar:
        """INCAR object."""
        return self["INCAR"]

    @property
    def kpoints(self) -> Kpoints | None:
        """KPOINTS object."""
        return self["KPOINTS"]

    @property
    def poscar(self) -> Poscar:
        """POSCAR object."""
        return self["POSCAR"]

    @property
    def potcar(self) -> Potcar | str | None:
        """POTCAR or POTCAR.spec object."""
        return self[self._potcar_filename]
