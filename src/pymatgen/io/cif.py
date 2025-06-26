"""Classes for Cif input and output from Structures."""

from __future__ import annotations

import math
import os
import re
import textwrap
import warnings
from collections import defaultdict, deque
from functools import partial
from inspect import getfullargspec
from io import StringIO
from itertools import groupby
from pathlib import Path
from typing import TYPE_CHECKING, Literal, cast

import numpy as np
from monty.dev import deprecated
from monty.io import zopen
from monty.serialization import loadfn

from pymatgen.core import Composition, DummySpecies, Element, Lattice, PeriodicSite, Species, Structure, get_el_sp
from pymatgen.core.operations import MagSymmOp, SymmOp
from pymatgen.electronic_structure.core import Magmom
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SpacegroupOperations
from pymatgen.symmetry.groups import SYMM_DATA, SpaceGroup
from pymatgen.symmetry.maggroups import MagneticSpaceGroup
from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.util.coord import find_in_coord_list_pbc, in_coord_list_pbc

if TYPE_CHECKING:
    from typing import Any

    from numpy.typing import NDArray
    from typing_extensions import Self

    from pymatgen.core import IStructure
    from pymatgen.util.typing import PathLike

__author__ = "Shyue Ping Ong, Will Richards, Matthew Horton"


class CifBlock:
    """
    Object for storing CIF data. All data is stored in a single dictionary.
    Data inside loops are stored in lists in the data dictionary, and
    information on which keys are grouped together are stored in the loops
    attribute.
    """

    max_len = 70  # not quite 80 so we can deal with semicolons and things

    def __init__(
        self,
        data: dict,
        loops: list[list[str]],
        header: str,
    ) -> None:
        """
        Args:
            data: dict of data to go into the CIF. Values should be convertible to string,
                or lists of these if the key is in a loop
            loops: list of lists of keys, grouped by which loop they should appear in
            header: name of the block (appears after the data_ on the first line).
        """
        self.loops = loops
        self.data = data
        # AJ (@computron) says: CIF Block names can't be more than 75 characters or you get an Exception
        self.header = header[:74]

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.loops == other.loops and self.data == other.data and self.header == other.header

    def __getitem__(self, key: str) -> Any:
        return self.data[key]

    def __str__(self) -> str:
        """Get the CIF string for the data block."""
        out = [f"data_{self.header}"]
        keys = list(self.data)
        written = []
        for key in keys:
            if key in written:
                continue
            for loop in self.loops:
                # search for a corresponding loop
                if key in loop:
                    out.append(self._loop_to_str(loop))
                    written.extend(loop)
                    break
            if key not in written:
                # key didn't belong to a loop
                val = self._format_field(self.data[key])
                if len(key) + len(val) + 3 < self.max_len:
                    out.append(f"{key}   {val}")
                else:
                    out.extend([key, val])
        return "\n".join(out)

    def _loop_to_str(self, loop: list[str]) -> str:
        """Convert a _loop block to string."""
        out = "loop_"
        for line in loop:
            out += "\n " + line

        for fields in zip(*(self.data[k] for k in loop), strict=True):
            line = "\n"
            for val in map(self._format_field, fields):
                if val[0] == ";":
                    out += f"{line}\n{val}"
                    line = "\n"
                elif len(line) + len(val) + 2 < self.max_len:
                    line += f"  {val}"
                else:
                    out += line
                    line = "\n  " + val
            out += line
        return out

    def _format_field(self, val: str) -> str:
        """Format field."""
        val = str(val).strip()

        if not val:
            return '""'

        # Wrap line if max length exceeded
        if len(val) > self.max_len:
            return f";\n{textwrap.fill(val, self.max_len)}\n;"

        # Add quotes if necessary
        if (" " in val or val[0] == "_") and (val[0] != "'" or val[-1] != "'") and (val[0] != '"' or val[-1] != '"'):
            quote = '"' if "'" in val else "'"
            val = quote + val + quote
        return val

    @classmethod
    def _process_string(cls, string: str) -> deque:
        """Process string to remove comments, empty lines and non-ASCII.
        Then break it into a stream of tokens.
        """
        # Remove comments
        string = re.sub(r"(\s|^)#.*$", "", string, flags=re.MULTILINE)

        # Remove empty lines
        string = re.sub(r"^\s*\n", "", string, flags=re.MULTILINE)

        # Remove non-ASCII
        string = string.encode("ascii", "ignore").decode("ascii")

        # Since line breaks in .cif files are mostly meaningless,
        # break up into a stream of tokens to parse, rejoin multiline
        # strings (between semicolons)
        deq: deque = deque()
        multiline: bool = False
        lines: list[str] = []

        # This regex splits on spaces, except when in quotes. Starting quotes must not be
        # preceded by non-whitespace (these get eaten by the first expression). Ending
        # quotes must not be followed by non-whitespace.
        pattern = re.compile(r"""([^'"\s][\S]*)|'(.*?)'(?!\S)|"(.*?)"(?!\S)""")

        for line in string.splitlines():
            if multiline:
                if line.startswith(";"):
                    multiline = False
                    deq.append(("", "", "", " ".join(lines)))
                    lines = []
                    line = line[1:].strip()
                else:
                    lines.append(line)
                    continue

            if line.startswith(";"):
                multiline = True
                lines.append(line[1:].strip())
            else:
                for string in pattern.findall(line):
                    # Location of the data in string depends on whether it was quoted in the input
                    deq.append(tuple(string))
        return deq

    @classmethod
    def from_str(cls, string: str) -> Self:
        """Read CifBlock from string.

        Args:
            string: String representation.

        Returns:
            CifBlock
        """
        deq: deque = cls._process_string(string)
        header: str = deq.popleft()[0][5:]
        data: dict = {}
        loops: list[list[str]] = []

        while deq:
            _str = deq.popleft()
            # CIF keys aren't in quotes, so show up as _str[0]
            if _str[0] == "_eof":
                break

            if _str[0].startswith("_"):
                try:
                    data[_str[0]] = "".join(deq.popleft())
                except IndexError:
                    data[_str[0]] = ""

            elif _str[0].startswith("loop_"):
                columns: list[str] = []
                items: list[str] = []
                while deq:
                    _str = deq[0]
                    if _str[0].startswith("loop_") or not _str[0].startswith("_"):
                        break
                    columns.append("".join(deq.popleft()))
                    data[columns[-1]] = []

                while deq:
                    _str = deq[0]
                    if _str[0].startswith(("loop_", "_")):
                        break
                    items.append("".join(deq.popleft()))

                n = len(items) // len(columns)
                if len(items) % n != 0:
                    raise ValueError(f"{len(items)=} is not a multiple of {n=}")
                loops.append(columns)
                for k, v in zip(columns * n, items, strict=True):
                    data[k].append(v.strip())

            elif issue := "".join(_str).strip():
                warnings.warn(f"Possible issue in CIF file at line: {issue}", stacklevel=2)

        return cls(data, loops, header)


class CifFile:
    """Read and parse CifBlocks from a .cif file or string."""

    def __init__(
        self,
        data: dict[str, CifBlock],
        orig_string: str | None = None,
        comment: str | None = None,
    ) -> None:
        """
        Args:
            data (dict): Of CifBlock objects.
            orig_string (str): The original CIF.
            comment (str): Comment.
        """
        self.data = data
        self.orig_string = orig_string
        self.comment: str = comment or "# generated using pymatgen"

    def __str__(self) -> str:
        out = "\n".join(map(str, self.data.values()))
        return f"{self.comment}\n{out}\n"

    @classmethod
    def from_str(cls, string: str) -> Self:
        """Read CifFile from a string.

        Args:
            string: String representation.

        Returns:
            CifFile
        """
        dct = {}

        for block_str in re.split(r"^\s*data_", f"x\n{string}", flags=re.MULTILINE | re.DOTALL)[1:]:
            # Skip over Cif block that contains powder diffraction data.
            # Some elements in this block were missing from CIF files in
            # Springer materials/Pauling file DBs.
            # This block does not contain any structure information anyway, and
            # CifParser was also not parsing it.
            if "powder_pattern" in re.split(r"\n", block_str, maxsplit=1)[0]:
                continue
            block = CifBlock.from_str(f"data_{block_str}")
            # TODO (@janosh, 2023-10-11) multiple CIF blocks with equal header will overwrite each other,
            # latest taking precedence. maybe something to fix and test e.g. in test_cif_writer_write_file
            dct[block.header] = block

        return cls(dct, string)

    @classmethod
    def from_file(cls, filename: PathLike) -> Self:
        """
        Read CifFile from a filename.

        Args:
            filename: Filename

        Returns:
            CifFile
        """
        with zopen(filename, mode="rt", errors="replace", encoding="utf-8") as file:
            return cls.from_str(file.read())  # type:ignore[arg-type]


class CifParser:
    """
    CIF file parser. Attempt to fix CIFs that are out-of-spec, but will issue warnings
    if corrections applied. These are also stored in the CifParser's warnings attribute.
    CIF file parser. Attempt to fix CIFs that are out-of-spec, but will issue warnings
    if corrections applied. These are also stored in the CifParser's errors attribute.
    """

    def __init__(
        self,
        filename: PathLike | StringIO,
        occupancy_tolerance: float = 1.0,
        site_tolerance: float = 1e-4,
        frac_tolerance: float = 1e-4,
        check_cif: bool = True,
        comp_tol: float = 0.01,
    ) -> None:
        """
        Args:
            filename (PathLike): CIF file, gzipped or bzipped CIF files are fine too.
            occupancy_tolerance (float): If total occupancy of a site is between
                1 and occupancy_tolerance, it will be scaled down to 1.
            site_tolerance (float): This tolerance is used to determine if two sites are at the same position,
                in which case they will be combined to a single disordered site. Defaults to 1e-4.
            frac_tolerance (float): This tolerance is used to determine is a coordinate should be rounded to an ideal
                value. e.g. 0.6667 is rounded to 2/3. This is desired if symmetry operations are going to be applied.
                However, for very large CIF files, this may need to be set to 0.
            check_cif (bool): Whether to check that stoichiometry reported in CIF matches
                that of resulting Structure, and whether elements are missing. Defaults to True.
            comp_tol (float): Tolerance for how closely stoichiometries of CIF file and pymatgen should match.
                Defaults to 0.01. Context: Experimental CIF files often don't report hydrogens positions due to being
                hard-to-locate with X-rays. pymatgen warns if the stoichiometry of the CIF file and the Structure
                don't match to within comp_tol.
        """

        def is_magcif() -> bool:
            """Check if a file is a magCIF file (heuristic)."""
            # Doesn't seem to be a canonical way to test if file is magCIF or
            # not, so instead check for magnetic symmetry datanames
            prefixes = [
                "_space_group_magn",
                "_atom_site_moment",
                "_space_group_symop_magn",
            ]
            for data in self._cif.data.values():
                for key in data.data:
                    for prefix in prefixes:
                        if prefix in key:
                            return True
            return False

        def is_magcif_incommensurate() -> bool:
            """
            Check if a file contains an incommensurate magnetic
            structure (heuristic).
            """
            # Doesn't seem to be a canonical way to test if magCIF file
            # describes incommensurate structure or not, so instead check
            # for common datanames
            if not self.feature_flags["magcif"]:
                return False
            prefixes = ["_cell_modulation_dimension", "_cell_wave_vector"]
            for data in self._cif.data.values():
                for key in data.data:
                    for prefix in prefixes:
                        if prefix in key:
                            return True
            return False

        # Take tolerances
        self._occupancy_tolerance = occupancy_tolerance
        self._site_tolerance = site_tolerance
        self._frac_tolerance = frac_tolerance

        # Read CIF file
        if isinstance(filename, str | Path):
            self._cif = CifFile.from_file(filename)
        elif isinstance(filename, StringIO):
            self._cif = CifFile.from_str(filename.read())
        else:
            raise TypeError("Unsupported file format.")

        # Options related to checking CIFs for missing elements
        # or incorrect stoichiometries
        self.check_cif = check_cif
        self.comp_tol = comp_tol

        # Store features from non-core CIF dictionaries, e.g. magCIF
        self.feature_flags: dict[Literal["magcif", "magcif_incommensurate"], bool] = {}
        self.feature_flags["magcif"] = is_magcif()
        self.feature_flags["magcif_incommensurate"] = is_magcif_incommensurate()

        # Store warnings during parsing
        self.warnings: list[str] = []

        # Pass individual CifBlocks to _sanitize_data
        for key in self._cif.data:
            self._cif.data[key] = self._sanitize_data(self._cif.data[key])

    @classmethod
    def from_str(cls, cif_string: str, **kwargs) -> Self:
        """Create a CifParser from a string.

        Args:
            cif_string (str): String representation of a CIF.

        Returns:
            CifParser
        """
        return cls(StringIO(cif_string), **kwargs)

    def _sanitize_data(self, data: CifBlock) -> CifBlock:
        """Some CIF files do not conform to spec. This method corrects
        known issues, particular in regards to Springer materials/
        Pauling files.

        This method is here so that CifParser can assume its
        input conforms to spec, simplifying its implementation.

        Args:
            data: CifBlock

        Returns:
            data CifBlock
        """
        # Handle formats of data as found in CIF files extracted
        # from the Springer Materials/Pauling File databases,
        # and that are different from standard ICSD formats.

        # Check for implicit hydrogens, warn if any present.
        if "_atom_site_attached_hydrogens" in data.data:
            attached_hydrogens = [str2float(x) for x in data.data["_atom_site_attached_hydrogens"] if str2float(x) != 0]
            if len(attached_hydrogens) > 0:
                self.warnings.append(
                    "Structure has implicit hydrogens defined, parsed structure unlikely to be "
                    "suitable for use in calculations unless hydrogens added."
                )

        # Check if "_atom_site_type_symbol" exists,
        # as some CIFs do not contain this key.
        if "_atom_site_type_symbol" in data.data:
            # Keep track of which data row needs to be removed.
            # Example of a row: Nb,Zr '0.8Nb + 0.2Zr' .2a .m-3m 0 0 0 1 14
            # 'rhombic dodecahedron, Nb<sub>14</sub>'
            # Without this, the above row in a structure would be parsed
            # as an ordered site with only Nb (since
            # CifParser would try to parse the first two characters of the
            # label "Nb,Zr") and occupancy=1.
            # However, this site is meant to be a disordered site with 0.8 of
            # Nb and 0.2 of Zr.
            idxs_to_remove: list[int] = []

            new_atom_site_label: list[str] = []
            new_atom_site_type_symbol: list[str] = []
            new_atom_site_occupancy: list[str] = []
            new_fract_x: list[str] = []
            new_fract_y: list[str] = []
            new_fract_z: list[str] = []

            for idx, el_row in enumerate(data["_atom_site_label"]):
                # CIF files from the Springer Materials/Pauling File have
                # switched the label and symbol. Thus, in the
                # above shown example row, '0.8Nb + 0.2Zr' is the symbol.
                # Below, we split the strings on ' + ' to
                # check if the length (or number of elements) in the label and
                # symbol are equal.
                if len(data["_atom_site_type_symbol"][idx].split(" + ")) > len(el_row.split(" + ")):
                    # Dictionary to hold elements and occupancies
                    els_occu: dict[str, float] = {}

                    # Parse symbol to get element names and occupancy
                    symbol_str: str = data["_atom_site_type_symbol"][idx]
                    symbol_str_lst: list[str] = symbol_str.split(" + ")

                    for _idx, symbol in enumerate(symbol_str_lst):
                        # Remove any bracketed items in the string
                        symbol_str_lst[_idx] = re.sub(r"\([0-9]*\)", "", symbol.strip())

                        # Extract element name and occupancy from the string
                        els_occu[str(re.findall(r"\D+", symbol_str_lst[_idx].strip())[1]).replace("<sup>", "")] = float(
                            "0" + re.findall(r"\.?\d+", symbol_str_lst[_idx].strip())[1]
                        )

                    for et, occu in els_occu.items():
                        # New atom site labels have "fix" appended
                        new_atom_site_label.append(f"{et}_fix{len(new_atom_site_label)}")
                        new_atom_site_type_symbol.append(et)
                        new_atom_site_occupancy.append(str(occu))

                        new_fract_x.append(str(str2float(data["_atom_site_fract_x"][idx])))
                        new_fract_y.append(str(str2float(data["_atom_site_fract_y"][idx])))
                        new_fract_z.append(str(str2float(data["_atom_site_fract_z"][idx])))

                    idxs_to_remove.append(idx)

            # Remove the original row by iterating over all keys in the CIF
            # data looking for lists, which indicates
            # multiple data items, one for each row, and remove items from the
            # list that corresponds to the removed row,
            # so that it's not processed by the rest of this function (which
            # would result in an error).
            for original_key in data.data:
                if isinstance(data.data[original_key], list):
                    for idx in sorted(idxs_to_remove, reverse=True):
                        del data.data[original_key][idx]

            if idxs_to_remove:
                self.warnings.append("Pauling file corrections applied.")

                data.data["_atom_site_label"] += new_atom_site_label
                data.data["_atom_site_type_symbol"] += new_atom_site_type_symbol
                data.data["_atom_site_occupancy"] += new_atom_site_occupancy
                data.data["_atom_site_fract_x"] += new_fract_x
                data.data["_atom_site_fract_y"] += new_fract_y
                data.data["_atom_site_fract_z"] += new_fract_z

        # This fixes inconsistencies in naming of several magCIF tags
        # as a result of magCIF being in widespread use prior to
        # specification being finalized (on advice of Branton Campbell).
        if self.feature_flags["magcif"]:
            # CIF-1 style has all underscores, interim standard
            # had period before magn instead of before the final
            # component (e.g. xyz).
            # We want to standardize keys, to simplify parsing.
            correct_keys = (
                "_space_group_symop_magn_operation.xyz",
                "_space_group_symop_magn_centering.xyz",
                "_space_group_magn.name_BNS",
                "_space_group_magn.number_BNS",
                "_atom_site_moment_crystalaxis_x",
                "_atom_site_moment_crystalaxis_y",
                "_atom_site_moment_crystalaxis_z",
                "_atom_site_moment_label",
            )

            # Cannot mutate dict during enumeration, so store changes
            changes_to_make = {}

            for original_key in data.data:
                for correct_key in correct_keys:
                    # convert to all underscore
                    trial_key = "_".join(correct_key.split("."))
                    test_key = "_".join(original_key.split("."))
                    if trial_key == test_key:
                        changes_to_make[correct_key] = original_key

            # Apply changes
            for correct_key, original_key in changes_to_make.items():
                data.data[correct_key] = data.data[original_key]

            # Map interim_keys to final_keys
            renamed_keys = {
                "_magnetic_space_group.transform_to_standard_Pp_abc": "_space_group_magn.transform_BNS_Pp_abc"
            }
            changes_to_make = {}

            for interim_key, final_key in renamed_keys.items():
                if data.data.get(interim_key):
                    changes_to_make[final_key] = interim_key

            if len(changes_to_make) > 0:
                self.warnings.append("Keys changed to match new magCIF specification.")

            for final_key, interim_key in changes_to_make.items():
                data.data[final_key] = data.data[interim_key]

        # Check for finite precision coordinates (e.g. 0.6667 instead of 0.6666666...),
        # which can cause issues when applying symmetry operations.
        important_fracs = (1 / 3, 2 / 3)
        fracs_to_change = {}
        for label in ("_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"):
            if label in data.data:
                for idx, frac in enumerate(data.data[label]):
                    try:
                        frac = str2float(frac)
                    except Exception:
                        # Coordinate might not be defined, e.g. '?'
                        continue

                    for comparison_frac in important_fracs:
                        if math.isclose(frac / comparison_frac, 1, abs_tol=self._frac_tolerance, rel_tol=0):
                            fracs_to_change[label, idx] = str(comparison_frac)

        if fracs_to_change:
            self.warnings.append(
                f"{len(fracs_to_change)} fractional coordinates rounded to ideal values to avoid issues with "
                "finite precision."
            )
            for (label, idx), val in fracs_to_change.items():
                data.data[label][idx] = val

        return data

    def _unique_coords(
        self,
        coords: list[tuple[float, float, float]],
        magmoms: list[Magmom] | None = None,
        lattice: Lattice | None = None,
        labels: dict[tuple[float, float, float], str] | None = None,
    ) -> tuple[list[NDArray], list[Magmom], list[str]]:
        """Generate unique coordinates using coordinates and symmetry
        positions, and their corresponding magnetic moments if supplied.
        """
        coords_out: list[NDArray] = []
        labels_out: list[str] = []
        labels = labels or {}

        if magmoms:
            if len(magmoms) != len(coords):
                raise ValueError("Length of magmoms and coords don't match.")

            magmoms_out: list[Magmom] = []
            for tmp_coord, tmp_magmom in zip(coords, magmoms, strict=True):
                for op in self.symmetry_operations:
                    coord = op.operate(tmp_coord)
                    coord = np.array([i - math.floor(i) for i in coord])
                    if isinstance(op, MagSymmOp):
                        # Up to this point, magmoms have been defined relative
                        # to crystal axis. Now convert to Cartesian and into
                        # a Magmom object.
                        if lattice is None:
                            raise ValueError("Lattice cannot be None.")
                        magmom = Magmom.from_moment_relative_to_crystal_axes(
                            op.operate_magmom(tmp_magmom), lattice=lattice
                        )
                    else:
                        magmom = Magmom(tmp_magmom)

                    if not in_coord_list_pbc(coords_out, coord, atol=self._site_tolerance):
                        coords_out.append(coord)
                        magmoms_out.append(magmom)
                        labels_out.append(labels.get(tmp_coord, "no_label"))

            return coords_out, magmoms_out, labels_out

        for tmp_coord in coords:
            for op in self.symmetry_operations:
                coord = op.operate(tmp_coord)
                coord = np.array([i - math.floor(i) for i in coord])
                if not in_coord_list_pbc(coords_out, coord, atol=self._site_tolerance):
                    coords_out.append(coord)
                    labels_out.append(labels.get(tmp_coord, "no_label"))

        dummy_magmoms = [Magmom(0)] * len(coords_out)
        return coords_out, dummy_magmoms, labels_out

    def get_lattice(
        self,
        data: CifBlock,
        length_strings=("a", "b", "c"),
        angle_strings=("alpha", "beta", "gamma"),
        lattice_type=None,
    ) -> Lattice | None:
        """Generate the lattice from the provided lattice parameters.
        In the absence of all six lattice parameters, the crystal system
        and necessary parameters are parsed.
        """
        try:
            return self.get_lattice_no_exception(
                data=data,
                angle_strings=angle_strings,
                lattice_type=lattice_type,
                length_strings=length_strings,
            )

        except KeyError:
            # Missing Key search for cell setting
            for lattice_label in (
                "_symmetry_cell_setting",
                "_space_group_crystal_system",
            ):
                if data.data.get(lattice_label):
                    lattice_type = data.data.get(lattice_label, "").lower()
                    try:
                        required_args = getfullargspec(getattr(Lattice, lattice_type)).args

                        lengths = (length for length in length_strings if length in required_args)
                        angles = (a for a in angle_strings if a in required_args)
                        return self.get_lattice(data, lengths, angles, lattice_type=lattice_type)
                    except AttributeError as exc:
                        self.warnings.append(str(exc))
                        warnings.warn(str(exc), stacklevel=2)

                else:
                    return None
        return None

    @staticmethod
    def get_lattice_no_exception(
        data: CifBlock,
        length_strings: tuple[str, str, str] = ("a", "b", "c"),
        angle_strings: tuple[str, str, str] = ("alpha", "beta", "gamma"),
        lattice_type: str | None = None,
    ) -> Lattice:
        """Convert a CifBlock to a pymatgen Lattice.

        Args:
            data: a dictionary of the CIF file
            length_strings: The strings that are used to identify the length parameters in the CIF file.
            angle_strings: The strings that are used to identify the angles in the CIF file.
            lattice_type (str): The type of lattice.

        Returns:
            Lattice object
        """
        lengths = [str2float(data[f"_cell_length_{i}"]) for i in length_strings]
        angles = [str2float(data[f"_cell_angle_{i}"]) for i in angle_strings]
        if not lattice_type:
            return Lattice.from_parameters(*lengths, *angles)
        return getattr(Lattice, lattice_type)(*(lengths + angles))

    def get_symops(self, data: CifBlock) -> list[SymmOp]:
        """
        Get the symmetry operations, in order to generate symmetry
        equivalent positions. If no symops are present, the space
        group symbol is parsed, and symops are generated.
        """
        sym_ops = []
        for symmetry_label in (
            "_symmetry_equiv_pos_as_xyz",
            "_symmetry_equiv_pos_as_xyz_",
            "_space_group_symop_operation_xyz",
            "_space_group_symop_operation_xyz_",
        ):
            if data.data.get(symmetry_label):
                xyz = data.data.get(symmetry_label)
                if xyz is None:
                    raise RuntimeError("Cannot get symmetry_label.")

                if isinstance(xyz, str):
                    msg = "A 1-line symmetry op P1 CIF is detected!"
                    warnings.warn(msg, stacklevel=2)
                    self.warnings.append(msg)
                    xyz = [xyz]
                try:
                    sym_ops = [SymmOp.from_xyz_str(s) for s in xyz]
                    break
                except ValueError:
                    continue

        sub_space_group = partial(re.sub, r"[\s_]", "")

        space_groups = {sub_space_group(key): key for key in SYMM_DATA["space_group_encoding"]}

        if not sym_ops:
            # Try to parse symbol
            for symmetry_label in (
                "_symmetry_space_group_name_H-M",
                "_symmetry_space_group_name_H_M",
                "_symmetry_space_group_name_H-M_",
                "_symmetry_space_group_name_H_M_",
                "_space_group_name_Hall",
                "_space_group_name_Hall_",
                "_space_group_name_H-M_alt",
                "_space_group_name_H-M_alt_",
                "_symmetry_space_group_name_hall",
                "_symmetry_space_group_name_hall_",
                "_symmetry_space_group_name_h-m",
                "_symmetry_space_group_name_h-m_",
            ):
                sg = data.data.get(symmetry_label)
                msg_template = "No _symmetry_equiv_pos_as_xyz type key found. Spacegroup from {} used."

                if sg:
                    sg = sub_space_group(sg)
                    try:
                        if spg := space_groups.get(sg):
                            sym_ops = list(SpaceGroup(spg).symmetry_ops)
                            msg = msg_template.format(symmetry_label)
                            warnings.warn(msg, stacklevel=2)
                            self.warnings.append(msg)
                            break
                    except ValueError:
                        pass

                    try:
                        cod_data = loadfn(
                            os.path.join(
                                os.path.dirname(os.path.dirname(__file__)),
                                "symmetry",
                                "symm_ops.json",
                            )
                        )
                        for _data in cod_data:
                            if sg == re.sub(r"\s+", "", _data["hermann_mauguin"]):
                                xyz = _data["symops"]
                                sym_ops = [SymmOp.from_xyz_str(s) for s in xyz]
                                msg = msg_template.format(symmetry_label)
                                warnings.warn(msg, stacklevel=2)
                                self.warnings.append(msg)
                                break
                    except Exception:
                        continue

                    if sym_ops:
                        break

        if not sym_ops:
            # Try to parse International number
            for symmetry_label in (
                "_space_group_IT_number",
                "_space_group_IT_number_",
                "_symmetry_Int_Tables_number",
                "_symmetry_Int_Tables_number_",
            ):
                if data.data.get(symmetry_label):
                    try:
                        integer = int(str2float(data.data.get(symmetry_label, "")))
                        sym_ops = list(SpaceGroup.from_int_number(integer).symmetry_ops)
                        break
                    except ValueError:
                        continue

        if not sym_ops:
            msg = "No _symmetry_equiv_pos_as_xyz type key found. Defaulting to P1."
            warnings.warn(msg, stacklevel=2)
            self.warnings.append(msg)
            sym_ops = [SymmOp.from_xyz_str("x, y, z")]

        return sym_ops

    def get_magsymops(self, data: CifBlock) -> list[MagSymmOp]:
        """Equivalent to get_symops except for magnetic symmetry groups.
        Separate function since additional operation for time reversal symmetry
        (which changes magnetic moments on sites) needs to be returned.
        """
        # Get BNS label and number for magnetic space group
        bns_name = data.data.get("_space_group_magn.name_BNS", "")
        bns_num = data.data.get("_space_group_magn.number_BNS", "")

        mag_symm_ops = []
        # Check if magCIF file explicitly contains magnetic symmetry operations
        if xyzt := data.data.get("_space_group_symop_magn_operation.xyz"):
            if isinstance(xyzt, str):
                xyzt = [xyzt]
            mag_symm_ops = [MagSymmOp.from_xyzt_str(s) for s in xyzt]

            if data.data.get("_space_group_symop_magn_centering.xyz"):
                xyzt = data.data.get("_space_group_symop_magn_centering.xyz")
                if isinstance(xyzt, str):
                    xyzt = [xyzt]
                centering_symops = [MagSymmOp.from_xyzt_str(s) for s in xyzt]

                all_ops = []
                for op in mag_symm_ops:
                    for centering_op in centering_symops:
                        new_translation = [
                            i - np.floor(i) for i in op.translation_vector + centering_op.translation_vector
                        ]
                        new_time_reversal = op.time_reversal * centering_op.time_reversal

                        all_ops.append(
                            MagSymmOp.from_rotation_and_translation_and_time_reversal(
                                rotation_matrix=op.rotation_matrix,
                                translation_vec=new_translation,
                                time_reversal=cast("Literal[-1, 1]", new_time_reversal),
                            )
                        )
                mag_symm_ops = all_ops

        # Else check if it specifies a magnetic space group
        elif bns_name or bns_num:
            label = bns_name or list(map(int, (bns_num.split("."))))

            if data.data.get("_space_group_magn.transform_BNS_Pp_abc") != "a,b,c;0,0,0":
                jonas_faithful = data.data.get("_space_group_magn.transform_BNS_Pp_abc")
                mag_sg = MagneticSpaceGroup(label, jonas_faithful)

            elif data.data.get("_space_group_magn.transform_BNS_Pp"):
                raise NotImplementedError("Incomplete specification to implement.")
            else:
                mag_sg = MagneticSpaceGroup(label)

            mag_symm_ops = mag_sg.symmetry_ops

        if not mag_symm_ops:
            msg = "No magnetic symmetry detected, using primitive symmetry."
            warnings.warn(msg, stacklevel=2)
            self.warnings.append(msg)
            mag_symm_ops = [MagSymmOp.from_xyzt_str("x, y, z, 1")]

        return mag_symm_ops

    @staticmethod
    def _parse_oxi_states(data: CifBlock) -> dict[str, float] | None:
        """Parse oxidation states from data."""
        try:
            oxi_states = {
                data["_atom_type_symbol"][i]: str2float(data["_atom_type_oxidation_number"][i])
                for i in range(len(data["_atom_type_symbol"]))
            }
            # Attempt to strip oxidation state from _atom_type_symbol
            # in case the label does not contain an oxidation state
            for idx, symbol in enumerate(data["_atom_type_symbol"]):
                oxi_states[re.sub(r"\d?[\+,\-]?$", "", symbol)] = str2float(data["_atom_type_oxidation_number"][idx])

        except (ValueError, KeyError):
            oxi_states = None
        return oxi_states

    @staticmethod
    def _parse_magmoms(data: CifBlock) -> dict[str, NDArray]:
        """Parse atomic magnetic moments from data."""
        try:
            magmoms = {
                data["_atom_site_moment_label"][idx]: np.array(
                    [
                        str2float(data["_atom_site_moment_crystalaxis_x"][idx]),
                        str2float(data["_atom_site_moment_crystalaxis_y"][idx]),
                        str2float(data["_atom_site_moment_crystalaxis_z"][idx]),
                    ]
                )
                for idx in range(len(data["_atom_site_moment_label"]))
            }
        except (ValueError, KeyError):
            return {}
        return magmoms

    def _parse_symbol(self, sym: str) -> str | None:
        """Parse a string with a symbol to extract a string representing an element.

        Args:
            sym (str): A symbol to be parsed.

        Returns:
            A string for the parsed symbol. None if no parsing was possible.
        """
        # Common representations for elements/water in CIF files
        # TODO: fix inconsistent handling of water
        special_syms = {
            "Hw": "H",
            "Ow": "O",
            "Wat": "O",
            "wat": "O",
            "OH": "",
            "OH2": "",
            "NO3": "N",
        }

        parsed_sym = None
        # Try with special symbols, otherwise check the first two letters,
        # then the first letter alone. If everything fails try extracting the
        # first letter.
        m_sp = re.match("|".join(special_syms), sym)
        if m_sp:
            parsed_sym = special_syms[m_sp.group()]
        elif Element.is_valid_symbol(sym[:2].title()):
            parsed_sym = sym[:2].title()
        elif Element.is_valid_symbol(sym[0].upper()):
            parsed_sym = sym[0].upper()
        elif match := re.match(r"w?[A-Z][a-z]*", sym):
            parsed_sym = match.group()

        if parsed_sym is not None and (m_sp or not re.match(rf"{parsed_sym}\d*", sym)):
            msg = f"{sym} parsed as {parsed_sym}"
            warnings.warn(msg, stacklevel=2)
            self.warnings.append(msg)

        return parsed_sym

    def _get_structure(
        self,
        data: CifBlock,
        primitive: bool,
        symmetrized: bool,
        check_occu: bool = False,
        min_thickness: float = 0.01,
    ) -> Structure | None:
        """Generate structure from part of the CIF.

        Args:
            data (CifBlock): The data block to parse.
            primitive (bool): Whether to return primitive unit cells.
            symmetrized (bool): Whether to return SymmetrizedStructure.
            check_occu (bool): Whether to check site for unphysical occupancy > 1.
            min_thickness (float): Minimum thickness in Angstrom to consider structure as valid.
                This is added to guard against unphysical small/thin structure,
                which could result in infinite loop for searching near neighbours.

        Returns:
            Structure or None if not found.
        """

        def get_num_implicit_hydrogens(symbol: str) -> int:
            """Get number of implicit hydrogens."""
            num_h = {"Wat": 2, "wat": 2, "O-H": 1}
            return num_h.get(symbol[:3], 0)

        def get_matching_coord(
            coord_to_species: dict[tuple[float, float, float], Composition],
            coord: tuple[float, float, float],
        ) -> tuple[float, float, float] | Literal[False]:
            """Find site by coordinate."""
            coords: list[tuple[float, float, float]] = list(coord_to_species)
            for op in self.symmetry_operations:
                frac_coord = op.operate(coord)
                indices: NDArray = find_in_coord_list_pbc(coords, frac_coord, atol=self._site_tolerance)
                if len(indices) > 0:
                    return coords[indices[0]]
            return False

        lattice = self.get_lattice(data)

        # Check minimal lattice thickness
        if lattice is not None:
            thickness = [lattice.d_hkl((1, 0, 0)), lattice.d_hkl((0, 1, 0)), lattice.d_hkl((0, 0, 1))]
            if any(t < min_thickness for t in thickness):
                raise ValueError(f"{thickness=} Å below threshold, double check your structure.")

        # If magCIF, get magnetic symmetry moments and magmoms
        # else standard CIF, and use empty magmom dict
        if self.feature_flags["magcif_incommensurate"]:
            raise NotImplementedError("Incommensurate structures not currently supported.")

        if self.feature_flags["magcif"]:
            if lattice is None:
                raise ValueError("Magmoms given in terms of crystal axes in magCIF spec.")
            self.symmetry_operations = self.get_magsymops(data)
            magmoms = self._parse_magmoms(data)

        else:
            self.symmetry_operations = self.get_symops(data)  # type:ignore[assignment]
            magmoms = {}

        oxi_states = self._parse_oxi_states(data)

        coord_to_species: dict[tuple[float, float, float], Composition] = {}
        coord_to_magmoms: dict[tuple[float, float, float], NDArray] = {}
        labels: dict[tuple[float, float, float], str] = {}

        for idx, label in enumerate(data["_atom_site_label"]):
            # If site type symbol exists, use it. Otherwise use the label
            try:
                symbol = self._parse_symbol(data["_atom_site_type_symbol"][idx])
                num_h = get_num_implicit_hydrogens(data["_atom_site_type_symbol"][idx])
            except KeyError:
                symbol = self._parse_symbol(label)
                num_h = get_num_implicit_hydrogens(label)

            if not symbol:
                continue

            # Get oxidation state
            if oxi_states is not None:
                o_s = oxi_states.get(symbol, 0)
                # Use _atom_site_type_symbol if possible for oxidation state
                if "_atom_site_type_symbol" in data.data:  # type: ignore[attr-defined]
                    oxi_symbol = data["_atom_site_type_symbol"][idx]
                    o_s = oxi_states.get(oxi_symbol, o_s)
                try:
                    el = Species(symbol, o_s)
                except Exception:
                    el = DummySpecies(symbol, o_s)
            else:
                el = get_el_sp(symbol)  # type: ignore[assignment]

            # Get occupancy
            try:
                occu: float = str2float(data["_atom_site_occupancy"][idx])
            except (KeyError, ValueError):
                occu = 1

            # If don't check_occu or the occupancy is greater than 0, create comp_dict
            if not check_occu or occu > 0:
                # Create site coordinate
                coord: tuple[float, float, float] = (
                    str2float(data["_atom_site_fract_x"][idx]),
                    str2float(data["_atom_site_fract_y"][idx]),
                    str2float(data["_atom_site_fract_z"][idx]),
                )

                # Create Composition
                comp_dict: dict[Species | str, float] = {el: max(occu, 1e-8)}

                if num_h > 0:
                    comp_dict["H"] = num_h
                    self.warnings.append(
                        "Structure has implicit hydrogens defined, parsed structure unlikely to be "
                        "suitable for use in calculations unless hydrogens added."
                    )
                comp = Composition(comp_dict)

                # Find matching site by coordinate
                match: tuple[float, float, float] | Literal[False] = get_matching_coord(coord_to_species, coord)
                if not match:
                    coord_to_species[coord] = comp
                    coord_to_magmoms[coord] = magmoms.get(label, np.array([0, 0, 0]))
                    labels[coord] = label

                else:
                    coord_to_species[match] += comp
                    # Disordered magnetic currently not supported
                    coord_to_magmoms[match] = None  # type:ignore[assignment]
                    labels[match] = label

        # Check occupancy
        _sum_occupancies: list[float] = [
            sum(comp.values())
            for comp in coord_to_species.values()
            if set(comp.elements) != {Element("O"), Element("H")}
        ]

        if any(occu > 1 for occu in _sum_occupancies):
            msg = (
                f"Some occupancies ({list(filter(lambda x: x > 1, _sum_occupancies))}) sum to > 1! If they are within "
                "the occupancy_tolerance, they will be rescaled. "
                f"The current occupancy_tolerance is set to: {self._occupancy_tolerance}"
            )
            warnings.warn(msg, stacklevel=2)
            self.warnings.append(msg)

        # Collect info for building Structure
        all_species: list[Composition] = []
        all_species_noedit: list[Composition] = []
        all_coords: list[tuple[float, float, float]] = []
        all_magmoms: list[Magmom] = []
        all_hydrogens: list[float] = []
        equivalent_indices: list[int] = []
        all_labels: list[str] = []

        # Check if a magCIF file is disordered
        if self.feature_flags["magcif"]:
            for val in coord_to_magmoms.values():
                if val is None:
                    # Proposed solution to this is to instead store magnetic
                    # moments as Species 'spin' property, instead of site
                    # property, but this introduces ambiguities for end user
                    # (such as unintended use of `spin` and Species will have
                    # fictitious oxidation state).
                    raise NotImplementedError("Disordered magnetic structures not currently supported.")

        if coord_to_species.items():
            for idx, (comp, group) in enumerate(
                groupby(
                    sorted(coord_to_species.items(), key=lambda x: x[1]),
                    key=lambda x: x[1],
                )
            ):
                tmp_coords: list[tuple[float, float, float]] = [site[0] for site in group]
                tmp_magmom: list[Magmom] = [coord_to_magmoms[tmp_coord] for tmp_coord in tmp_coords]

                if self.feature_flags["magcif"]:
                    coords, _magmoms, new_labels = self._unique_coords(
                        tmp_coords,
                        magmoms=tmp_magmom,
                        labels=labels,
                        lattice=lattice,
                    )
                else:
                    coords, _magmoms, new_labels = self._unique_coords(tmp_coords, labels=labels)

                if set(comp.elements) == {Element("O"), Element("H")}:
                    # O with implicit hydrogens
                    im_h = comp["H"]
                    species = Composition({"O": comp["O"]})
                else:
                    im_h = 0
                    species = comp

                # The following might be a more natural representation of equivalent indices,
                # but is not in the format expect by SymmetrizedStructure:
                #   equivalent_indices.append(list(range(len(all_coords), len(coords)+len(all_coords))))
                # The above gives a list like:
                #   [[0, 1, 2, 3], [4, 5, 6, 7, 8, 9, 10, 11]] where the
                # integers are site indices, whereas the version used below will give a version like:
                #   [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
                # which is a list in the same order as the sites, but where if a site has the same integer
                # it is equivalent.
                equivalent_indices += len(coords) * [idx]

                all_hydrogens.extend(len(coords) * [im_h])
                all_coords.extend(coords)  # type:ignore[arg-type]
                all_species.extend(len(coords) * [species])
                all_magmoms.extend(_magmoms)
                all_labels.extend(new_labels)

            # Scale occupancies if necessary
            all_species_noedit = all_species.copy()  # save copy before scaling in case of check_occu=False, used below
            for idx, species in enumerate(all_species):
                total_occu = sum(species.values())
                if check_occu and total_occu > self._occupancy_tolerance:
                    raise ValueError(f"Occupancy {total_occu} exceeded tolerance.")

                if total_occu > 1:
                    all_species[idx] = species / total_occu

        if all_species and len(all_species) == len(all_coords) and len(all_species) == len(all_magmoms):
            site_properties: dict[str, list] = {}
            if any(all_hydrogens):
                if len(all_hydrogens) != len(all_coords):
                    raise ValueError("lengths of all_hydrogens and all_coords mismatch")
                site_properties["implicit_hydrogens"] = all_hydrogens

            if self.feature_flags["magcif"]:
                site_properties["magmom"] = all_magmoms

            if not site_properties:
                site_properties = {}

            if any(all_labels):
                if len(all_labels) != len(all_species):
                    raise ValueError("lengths of all_labels and all_species mismatch")
            else:
                all_labels = None  # type: ignore[assignment]

            struct: Structure = Structure(
                lattice,  # type:ignore[arg-type]
                all_species,
                all_coords,
                site_properties=site_properties,
                labels=all_labels,
            )

            if symmetrized:
                # Wyckoff labels not currently parsed, note that not all CIFs will contain Wyckoff labels
                # TODO: extract Wyckoff labels (or other CIF attributes) and include as site_properties
                wyckoffs = ["Not Parsed"] * len(struct)

                # space groups names are likewise not parsed (again, not all CIFs will contain this information)
                # What is stored are the lists of symmetry operations used to generate the structure
                # TODO: ensure space group labels are stored if present
                sg = SpacegroupOperations("Not Parsed", -1, self.symmetry_operations)
                struct = SymmetrizedStructure(struct, sg, equivalent_indices, wyckoffs)

            if not check_occu:
                if lattice is None:
                    raise RuntimeError("Cannot generate Structure with lattice as None.")

                for idx in range(len(struct)):
                    struct[idx] = PeriodicSite(
                        all_species_noedit[idx],
                        all_coords[idx],
                        lattice,
                        properties=site_properties,
                        label=all_labels[idx],
                        skip_checks=True,
                    )

            if symmetrized or not check_occu:
                return struct

            struct = struct.get_sorted_structure()

            if primitive and self.feature_flags["magcif"]:
                struct = struct.get_primitive_structure(use_site_props=True)
            elif primitive:
                struct = struct.get_primitive_structure()
                struct = struct.get_reduced_structure()

            if self.check_cif:
                cif_failure_reason = self.check(struct)
                if cif_failure_reason is not None:
                    warnings.warn(cif_failure_reason, stacklevel=2)

            return struct
        return None

    def parse_structures(
        self,
        primitive: bool | None = None,
        symmetrized: bool = False,
        check_occu: bool = True,
        on_error: Literal["ignore", "warn", "raise"] = "warn",
    ) -> list[Structure]:
        """Return list of structures in CIF file.

        Args:
            primitive (bool): Whether to return primitive unit cells.
                Defaults to False. With magnetic CIF files, will return primitive
                magnetic cell which may be larger than nuclear primitive cell.
            symmetrized (bool): Whether to return a SymmetrizedStructure which will
                include the equivalent indices and symmetry operations used to
                create the Structure as provided by the CIF (if explicit symmetry
                operations are included in the CIF) or generated from information
                in the CIF (if only space group labels are provided). Note that
                currently Wyckoff labels and space group labels or numbers are
                not included in the generated SymmetrizedStructure, these will be
                notated as "Not Parsed" or -1 respectively.
            check_occu (bool): Whether to check site for unphysical occupancy > 1.
                Useful for experimental results in which occupancy was allowed to
                refine to unphysical values. Warning: unphysical occupancies are
                incompatible with many pymatgen features. Defaults to True.
            on_error ("ignore" | "warn" | "raise"): What to do in case of KeyError
                or ValueError while parsing CIF file. Defaults to "warn".

        Returns:
            list[Structure]: All structures in CIF file.
        """
        if primitive is None:
            primitive = False
            warnings.warn(
                "The default value of primitive was changed from True to False in "
                "https://github.com/materialsproject/pymatgen/pull/3419. CifParser now returns the cell "
                "in the CIF file as is. If you want the primitive cell, please set primitive=True explicitly.",
                stacklevel=2,
            )

        if primitive and symmetrized:
            raise ValueError(
                "Using both 'primitive' and 'symmetrized' arguments is not currently supported "
                "since unexpected behavior might result."
            )

        structures: list[Structure] = []
        for idx, data in enumerate(self._cif.data.values()):
            try:
                if struct := self._get_structure(data, primitive, symmetrized, check_occu=check_occu):
                    structures.append(struct)

            except (KeyError, ValueError) as exc:
                msg = f"No structure parsed for section {idx + 1} in CIF.\n{exc}"
                if on_error == "raise":
                    raise ValueError(msg) from exc
                if on_error == "warn":
                    warnings.warn(msg, stacklevel=2)
                self.warnings.append(msg)

        if self.warnings and on_error == "warn":
            warnings.warn("Issues encountered while parsing CIF: " + "\n".join(self.warnings), stacklevel=2)

        if not structures:
            raise ValueError("Invalid CIF file with no structures!")
        return structures

    @deprecated(
        parse_structures,
        message="The only difference is that primitive defaults to False in the new parse_structures method."
        "So parse_structures(primitive=True) is equivalent to the old behavior of get_structures().",
    )
    def get_structures(self, *args, **kwargs) -> list[Structure]:
        """
        Deprecated, use parse_structures instead. Only difference between
        these two methods is the default primitive=False in parse_structures.
        So parse_structures(primitive=True) is equivalent to the default
        behaviour of get_structures().
        """
        # Extract primitive if passed as arg
        if len(args) > 0:
            kwargs["primitive"] = args[0]
            args = args[1:]
        kwargs.setdefault("primitive", True)
        return self.parse_structures(*args, **kwargs)

    def get_bibtex_string(self) -> str:
        """Get BibTeX reference from CIF file.

        TODO:
            - parse '_publ_section_references' when it exists?
            - CIF specification supports multiple citations.

        Returns:
            BibTeX string.
        """
        from bibtexparser.bibdatabase import BibDatabase
        from bibtexparser.bwriter import BibTexWriter

        bibtex_keys: dict[str, tuple[str, ...]] = {
            "author": ("_publ_author_name", "_citation_author_name"),
            "title": ("_publ_section_title", "_citation_title"),
            "journal": (
                "_journal_name_full",
                "_journal_name_abbrev",
                "_citation_journal_full",
                "_citation_journal_abbrev",
            ),
            "volume": ("_journal_volume", "_citation_journal_volume"),
            "year": ("_journal_year", "_citation_year"),
            "number": ("_journal_number", "_citation_number"),
            "page_first": ("_journal_page_first", "_citation_page_first"),
            "page_last": ("_journal_page_last", "_citation_page_last"),
            "doi": ("_journal_DOI", "_citation_DOI"),
        }

        db = BibDatabase()
        db.entries = []

        for idx, data in enumerate(self._cif.data.values()):
            # Convert to lower-case keys, some CIF files inconsistent
            _data = {k.lower(): v for k, v in data.data.items()}
            entry = {"ENTRYTYPE": "article", "ID": f"cifref{idx}"}

            for field, tags in bibtex_keys.items():
                for tag in tags:
                    if tag in _data:
                        value = _data[tag]
                        entry[field] = value[0] if isinstance(value, list) else value
                        break

            # Convert to bibtex author format ("and" delimited)
            if "author" in entry:
                # Separate out semicolon authors
                if isinstance(entry["author"], str) and ";" in entry["author"]:
                    entry["author"] = entry["author"].split(";")  # type:ignore[assignment]
                if isinstance(entry["author"], list):
                    entry["author"] = " and ".join(entry["author"])

            # Convert to bibtex page range format, use empty string if not specified
            if "page_first" in entry or "page_last" in entry:
                entry["pages"] = f"{entry.get('page_first', '')}--{entry.get('page_last', '')}"
                entry.pop("page_first", None)  # and remove page_first, page_list if present
                entry.pop("page_last", None)

            db.entries.append(entry)

        # NOTE: the following is added to make output consistent with
        # previous pybtex implementation
        writer = BibTexWriter()
        writer.indent = "    "
        writer.display_order = ("author", "title", "journal", "volume", "year", "pages")

        # Replace curly brackets with double quotes (skip the first and last one)
        return re.sub(r"(^\s*\w+\s*=\s*)\{([^{}]*)\}", r'\1"\2"', writer.write(db), flags=re.MULTILINE)

    def as_dict(self) -> dict:
        """MSONable dict."""
        dct: dict = {}
        for key, val in self._cif.data.items():
            dct[key] = {}
            for sub_key, sub_val in val.data.items():
                dct[key][sub_key] = sub_val
        return dct

    @property
    def has_errors(self) -> bool:
        """Whether there are errors/warnings detected in CIF parsing."""
        return len(self.warnings) > 0

    def check(self, structure: Structure) -> str | None:
        """Check whether a Structure created from CIF passes sanity checks.

        Checks:
            - Composition from CIF is valid
            - CIF composition contains only valid elements
            - CIF and structure contain the same elements (often hydrogens
                are omitted from CIFs, as their positions cannot be determined from
                X-ray diffraction, needs more difficult neutron diffraction)
            -  CIF and structure have same relative stoichiometry. Thus
                if CIF reports stoichiometry LiFeO, and the structure has
                composition (LiFeO)4, this check passes.

        Args:
            structure (Structure) : Structure created from CIF.

        Returns:
            str | None: If any check fails, return a human-readable str for the
                reason (e.g., which elements are missing). None if all checks pass.
        """
        cif_as_dict = self.as_dict()
        head_key = next(iter(cif_as_dict))

        cif_formula = None
        for key in ("_chemical_formula_sum", "_chemical_formula_structural"):
            if cif_as_dict[head_key].get(key):
                cif_formula = cif_as_dict[head_key][key]
                break

        # In case of missing CIF formula keys, get non-stoichiometric formula from
        # unique sites and skip relative stoichiometry check (added in gh-3628)
        check_stoichiometry = True
        if cif_formula is None and cif_as_dict[head_key].get("_atom_site_type_symbol"):
            check_stoichiometry = False
            cif_formula = " ".join(cif_as_dict[head_key]["_atom_site_type_symbol"])

        try:
            cif_composition = Composition(cif_formula)
        except Exception as exc:
            return f"Cannot determine chemical composition from CIF! {exc}"

        try:
            orig_comp = cif_composition.remove_charges().as_dict()
            struct_comp = structure.composition.remove_charges().as_dict()
        except Exception as exc:
            return str(exc)

        orig_comp_elts = {str(elt) for elt in orig_comp}
        struct_comp_elts = {str(elt) for elt in struct_comp}
        failure_reason: str | None = None

        # Hard failure: missing elements
        if orig_comp_elts != struct_comp_elts:
            missing = set(orig_comp_elts).difference(set(struct_comp_elts))
            addendum = "from PMG structure composition"
            if not missing:
                addendum = "from CIF-reported composition"
                missing = set(struct_comp_elts).difference(set(orig_comp_elts))
            missing_str = ", ".join([str(x) for x in missing])
            failure_reason = f"Missing elements {missing_str} {addendum}"

        elif any(struct_comp[elt] - orig_comp[elt] != 0 for elt in orig_comp):
            # Check that CIF/PMG stoichiometry has same relative ratios of elements
            if check_stoichiometry:
                ratios = {elt: struct_comp[elt] / orig_comp[elt] for elt in orig_comp_elts}

                same_stoich = all(
                    abs(ratios[elt_a] - ratios[elt_b]) < self.comp_tol
                    for elt_a in orig_comp_elts
                    for elt_b in orig_comp_elts
                )

                if not same_stoich:
                    failure_reason = f"Incorrect stoichiometry:\n  CIF={orig_comp}\n  PMG={struct_comp}\n  {ratios=}"
            else:
                self.warnings += ["Skipping relative stoichiometry check because CIF does not contain formula keys."]

        return failure_reason


def str2float(text: str) -> float:
    """Remove uncertainty brackets from strings and return the float."""
    try:
        # Note that the ending ) is sometimes missing. That is why the code has
        # been modified to treat it as optional. Same logic applies to lists.
        return float(re.sub(r"\(.+\)*", "", text))

    except TypeError:
        if isinstance(text, list) and len(text) == 1:
            return float(re.sub(r"\(.+\)*", "", text[0]))

    except ValueError:
        if text.strip() == ".":
            return 0
        raise
    raise ValueError(f"{text!s} cannot be converted to float")


class CifWriter:
    """A wrapper around CifFile to write CIF files from pymatgen Structure."""

    def __init__(
        self,
        struct: Structure | IStructure,
        symprec: float | None = None,
        write_magmoms: bool = False,
        significant_figures: int = 8,
        angle_tolerance: float = 5,
        refine_struct: bool = True,
        write_site_properties: bool = False,
    ) -> None:
        """
        Args:
            struct (Structure): structure to write.
            symprec (float): If not none, finds the symmetry of the structure
                and writes the CIF with symmetry information. Passes symprec
                to the SpacegroupAnalyzer. See also refine_struct.
            write_magmoms (bool): If True, will write magCIF file. Incompatible
                with symprec
            significant_figures (int): Specifies precision for formatting of floats.
                Defaults to 8.
            angle_tolerance (float): Angle tolerance for symmetry finding. Passes
                angle_tolerance to the SpacegroupAnalyzer. Used only if symprec
                is not None.
            refine_struct: Used only if symprec is not None. If True, get_refined_structure
                is invoked to convert input structure from primitive to conventional.
            write_site_properties (bool): Whether to write the Structure.site_properties
                to the CIF as _atom_site_{property name}. Defaults to False.
        """
        if write_magmoms and symprec is not None:
            warnings.warn(
                "Magnetic symmetry cannot currently be detected by pymatgen, disabling symmetry detection.",
                stacklevel=2,
            )
            symprec = None

        blocks: dict[str, Any] = {}
        spacegroup: tuple[str, int] = ("P 1", 1)
        if symprec is not None:
            spg_analyzer = SpacegroupAnalyzer(struct, symprec, angle_tolerance=angle_tolerance)
            spacegroup = (
                spg_analyzer.get_space_group_symbol(),
                spg_analyzer.get_space_group_number(),
            )

            if refine_struct:
                # Need the refined structure when using symprec. This converts
                # primitive to conventional structures, the standard for CIF.
                struct = spg_analyzer.get_refined_structure()

        lattice = struct.lattice
        comp = struct.composition
        no_oxi_comp = comp.element_composition
        format_str: str = f"{{:.{significant_figures}f}}"
        blocks["_symmetry_space_group_name_H-M"] = spacegroup[0]
        for cell_attr in ("a", "b", "c"):
            blocks[f"_cell_length_{cell_attr}"] = format_str.format(getattr(lattice, cell_attr))
        for cell_attr in ("alpha", "beta", "gamma"):
            blocks[f"_cell_angle_{cell_attr}"] = format_str.format(getattr(lattice, cell_attr))
        blocks["_symmetry_Int_Tables_number"] = spacegroup[1]
        blocks["_chemical_formula_structural"] = no_oxi_comp.reduced_formula
        blocks["_chemical_formula_sum"] = no_oxi_comp.formula
        blocks["_cell_volume"] = format_str.format(lattice.volume)

        _, fu = no_oxi_comp.get_reduced_composition_and_factor()
        blocks["_cell_formula_units_Z"] = str(int(fu))

        if symprec is None:
            blocks["_symmetry_equiv_pos_site_id"] = ["1"]
            blocks["_symmetry_equiv_pos_as_xyz"] = ["x, y, z"]

        else:
            spg_analyzer = SpacegroupAnalyzer(struct, symprec)
            symm_ops: list[SymmOp] = []
            for op in spg_analyzer.get_symmetry_operations():
                v = op.translation_vector
                symm_ops.append(SymmOp.from_rotation_and_translation(op.rotation_matrix, v))

            ops = [op.as_xyz_str() for op in symm_ops]
            blocks["_symmetry_equiv_pos_site_id"] = [f"{i}" for i in range(1, len(ops) + 1)]
            blocks["_symmetry_equiv_pos_as_xyz"] = ops

        loops: list[list[str]] = [
            ["_symmetry_equiv_pos_site_id", "_symmetry_equiv_pos_as_xyz"],
        ]

        try:
            symbol_to_oxi_num = {str(el): float(el.oxi_state or 0) for el in sorted(comp.elements)}
            blocks["_atom_type_symbol"] = list(symbol_to_oxi_num)
            blocks["_atom_type_oxidation_number"] = symbol_to_oxi_num.values()
            loops.append(["_atom_type_symbol", "_atom_type_oxidation_number"])
        except (TypeError, AttributeError):
            symbol_to_oxi_num = {el.symbol: 0 for el in sorted(comp.elements)}

        atom_site_type_symbol = []
        atom_site_symmetry_multiplicity = []
        atom_site_fract_x = []
        atom_site_fract_y = []
        atom_site_fract_z = []
        atom_site_label = []
        atom_site_occupancy = []
        atom_site_moment_label = []
        atom_site_moment_crystalaxis_x = []
        atom_site_moment_crystalaxis_y = []
        atom_site_moment_crystalaxis_z = []
        atom_site_properties: dict[str, list] = defaultdict(list)
        count = 0
        if symprec is None:
            for site in struct:
                for sp, occu in sorted(site.species.items()):
                    atom_site_type_symbol.append(str(sp))
                    atom_site_symmetry_multiplicity.append("1")
                    atom_site_fract_x.append(format_str.format(site.a))
                    atom_site_fract_y.append(format_str.format(site.b))
                    atom_site_fract_z.append(format_str.format(site.c))
                    atom_site_occupancy.append(str(occu))
                    site_label = f"{sp.symbol}{count}"

                    if "magmom" in site.properties:
                        mag = site.properties["magmom"]
                    elif getattr(sp, "spin", None) is not None:
                        mag = sp.spin
                    else:
                        # Use site label if available for regular sites
                        site_label = site.label if site.label != site.species_string else site_label
                        mag = 0

                    atom_site_label.append(site_label)

                    magmom = Magmom(mag)
                    if write_magmoms and abs(magmom) > 0:
                        moment = Magmom.get_moment_relative_to_crystal_axes(magmom, lattice)
                        atom_site_moment_label.append(f"{sp.symbol}{count}")
                        atom_site_moment_crystalaxis_x.append(format_str.format(moment[0]))
                        atom_site_moment_crystalaxis_y.append(format_str.format(moment[1]))
                        atom_site_moment_crystalaxis_z.append(format_str.format(moment[2]))

                    if write_site_properties:
                        for key, val in site.properties.items():
                            atom_site_properties[key].append(format_str.format(val))

                    count += 1

        else:
            # The following just presents a deterministic ordering
            unique_sites = [
                (
                    min(sites, key=lambda site: tuple(abs(x) for x in site.frac_coords)),
                    len(sites),
                )
                for sites in spg_analyzer.get_symmetrized_structure().equivalent_sites  # type: ignore[reportPossiblyUnboundVariable]
            ]
            for site, mult in sorted(
                unique_sites,
                key=lambda t: (
                    t[0].species.average_electroneg,
                    -t[1],
                    t[0].a,
                    t[0].b,
                    t[0].c,
                ),
            ):
                for sp, occu in site.species.items():
                    atom_site_type_symbol.append(str(sp))
                    atom_site_symmetry_multiplicity.append(f"{mult}")
                    atom_site_fract_x.append(format_str.format(site.a))
                    atom_site_fract_y.append(format_str.format(site.b))
                    atom_site_fract_z.append(format_str.format(site.c))
                    site_label = site.label if site.label != site.species_string else f"{sp.symbol}{count}"
                    atom_site_label.append(site_label)
                    atom_site_occupancy.append(str(occu))
                    count += 1

        if len(set(atom_site_label)) != len(atom_site_label):
            warnings.warn(
                "Site labels are not unique, which is not compliant with the CIF spec "
                "(https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_label.html):"
                f"`{atom_site_label}`.",
                stacklevel=2,
            )

        blocks["_atom_site_type_symbol"] = atom_site_type_symbol
        blocks["_atom_site_label"] = atom_site_label
        blocks["_atom_site_symmetry_multiplicity"] = atom_site_symmetry_multiplicity
        blocks["_atom_site_fract_x"] = atom_site_fract_x
        blocks["_atom_site_fract_y"] = atom_site_fract_y
        blocks["_atom_site_fract_z"] = atom_site_fract_z
        blocks["_atom_site_occupancy"] = atom_site_occupancy
        loop_labels = [
            "_atom_site_type_symbol",
            "_atom_site_label",
            "_atom_site_symmetry_multiplicity",
            "_atom_site_fract_x",
            "_atom_site_fract_y",
            "_atom_site_fract_z",
            "_atom_site_occupancy",
        ]
        if write_site_properties:
            for key, vals in atom_site_properties.items():
                blocks[f"_atom_site_{key}"] = vals
                loop_labels += [f"_atom_site_{key}"]
        loops.append(loop_labels)

        if write_magmoms:
            blocks["_atom_site_moment_label"] = atom_site_moment_label
            blocks["_atom_site_moment_crystalaxis_x"] = atom_site_moment_crystalaxis_x
            blocks["_atom_site_moment_crystalaxis_y"] = atom_site_moment_crystalaxis_y
            blocks["_atom_site_moment_crystalaxis_z"] = atom_site_moment_crystalaxis_z
            loops.append(
                [
                    "_atom_site_moment_label",
                    "_atom_site_moment_crystalaxis_x",
                    "_atom_site_moment_crystalaxis_y",
                    "_atom_site_moment_crystalaxis_z",
                ]
            )
        dct = {comp.reduced_formula: CifBlock(blocks, loops, comp.reduced_formula)}
        self._cf = CifFile(dct)

    def __str__(self) -> str:
        """The CIF as a string."""
        return str(self._cf)

    @property
    def cif_file(self) -> CifFile:
        """CifFile associated with the CifWriter."""
        return self._cf

    def write_file(
        self,
        filename: PathLike,
        mode: Literal["wt", "at"] = "wt",
    ) -> None:
        """Write the CIF file."""
        with zopen(filename, mode=mode, encoding="utf-8") as file:
            file.write(str(self))  # type:ignore[arg-type]
