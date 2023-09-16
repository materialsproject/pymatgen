"""Wrapper classes for Cif input and output from Structures."""

from __future__ import annotations

import math
import os
import re
import textwrap
import warnings
from collections import deque
from functools import partial
from inspect import getfullargspec as getargspec
from io import StringIO
from itertools import groupby
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
from monty.io import zopen
from monty.serialization import loadfn
from monty.string import remove_non_ascii

from pymatgen.core.composition import Composition
from pymatgen.core.lattice import Lattice
from pymatgen.core.operations import MagSymmOp, SymmOp
from pymatgen.core.periodic_table import DummySpecies, Element, Species, get_el_sp
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Magmom
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SpacegroupOperations
from pymatgen.symmetry.groups import SYMM_DATA, SpaceGroup
from pymatgen.symmetry.maggroups import MagneticSpaceGroup
from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.util.coord import find_in_coord_list_pbc, in_coord_list_pbc

if TYPE_CHECKING:
    from pymatgen.core.trajectory import Vector3D

__author__ = "Shyue Ping Ong, Will Richards, Matthew Horton"

sub_spgrp = partial(re.sub, r"[\s_]", "")

space_groups = {sub_spgrp(key): key for key in SYMM_DATA["space_group_encoding"]}  # type: ignore

space_groups.update({sub_spgrp(key): key for key in SYMM_DATA["space_group_encoding"]})  # type: ignore


class CifBlock:
    """
    Object for storing cif data. All data is stored in a single dictionary.
    Data inside loops are stored in lists in the data dictionary, and
    information on which keys are grouped together are stored in the loops
    attribute.
    """

    maxlen = 70  # not quite 80 so we can deal with semicolons and things

    def __init__(self, data, loops, header):
        """
        Args:
            data: dict of data to go into the cif. Values should be convertible to string,
                or lists of these if the key is in a loop
            loops: list of lists of keys, grouped by which loop they should appear in
            header: name of the block (appears after the data_ on the first line).
        """
        self.loops = loops
        self.data = data
        # AJ says: CIF Block names cannot be more than 75 characters or you
        # get an Exception
        self.header = header[:74]

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, CifBlock):
            return NotImplemented
        return self.loops == other.loops and self.data == other.data and self.header == other.header

    def __getitem__(self, key):
        return self.data[key]

    def __str__(self):
        """Returns the cif string for the data block."""
        out = [f"data_{self.header}"]
        keys = list(self.data)
        written = []
        for key in keys:
            if key in written:
                continue
            for loop in self.loops:
                # search for a corresponding loop
                if key in loop:
                    out.append(self._loop_to_string(loop))
                    written.extend(loop)
                    break
            if key not in written:
                # k didn't belong to a loop
                v = self._format_field(self.data[key])
                if len(key) + len(v) + 3 < self.maxlen:
                    out.append(f"{key}   {v}")
                else:
                    out.extend([key, v])
        return "\n".join(out)

    def _loop_to_string(self, loop):
        out = "loop_"
        for line in loop:
            out += "\n " + line
        for fields in zip(*(self.data[k] for k in loop)):
            line = "\n"
            for val in map(self._format_field, fields):
                if val[0] == ";":
                    out += line + "\n" + val
                    line = "\n"
                elif len(line) + len(val) + 2 < self.maxlen:
                    line += "  " + val
                else:
                    out += line
                    line = "\n  " + val
            out += line
        return out

    def _format_field(self, v):
        v = str(v).strip()
        if len(v) > self.maxlen:
            return ";\n" + textwrap.fill(v, self.maxlen) + "\n;"
        # add quotes if necessary
        if v == "":
            return '""'
        if (" " in v or v[0] == "_") and not (v[0] == "'" and v[-1] == "'") and not (v[0] == '"' and v[-1] == '"'):
            q = '"' if "'" in v else "'"
            v = q + v + q
        return v

    @classmethod
    def _process_string(cls, string):
        # remove comments
        string = re.sub(r"(\s|^)#.*$", "", string, flags=re.MULTILINE)
        # remove empty lines
        string = re.sub(r"^\s*\n", "", string, flags=re.MULTILINE)
        # remove non_ascii
        string = remove_non_ascii(string)
        # since line breaks in .cif files are mostly meaningless,
        # break up into a stream of tokens to parse, rejoining multiline
        # strings (between semicolons)
        q = deque()
        multiline = False
        ml = []
        # this regex splits on spaces, except when in quotes.
        # starting quotes must not be preceded by non-whitespace
        # (these get eaten by the first expression)
        # ending quotes must not be followed by non-whitespace
        p = re.compile(r"""([^'"\s][\S]*)|'(.*?)'(?!\S)|"(.*?)"(?!\S)""")
        for line in string.splitlines():
            if multiline:
                if line.startswith(";"):
                    multiline = False
                    q.append(("", "", "", " ".join(ml)))
                    ml = []
                    line = line[1:].strip()
                else:
                    ml.append(line)
                    continue
            if line.startswith(";"):
                multiline = True
                ml.append(line[1:].strip())
            else:
                for s in p.findall(line):
                    # s is tuple. location of the data in the tuple
                    # depends on whether it was quoted in the input
                    q.append(tuple(s))
        return q

    @classmethod
    @np.deprecate(message="Use from_str instead")
    def from_string(cls, *args, **kwargs):
        return cls.from_str(*args, **kwargs)

    @classmethod
    def from_str(cls, string):
        """
        Reads CifBlock from string.

        :param string: String representation.

        Returns:
            CifBlock
        """
        q = cls._process_string(string)
        header = q.popleft()[0][5:]
        data = {}
        loops = []
        while q:
            s = q.popleft()
            # cif keys aren't in quotes, so show up in s[0]
            if s[0] == "_eof":
                break
            if s[0].startswith("_"):
                try:
                    data[s[0]] = "".join(q.popleft())
                except IndexError:
                    data[s[0]] = ""
            elif s[0].startswith("loop_"):
                columns = []
                items = []
                while q:
                    s = q[0]
                    if s[0].startswith("loop_") or not s[0].startswith("_"):
                        break
                    columns.append("".join(q.popleft()))
                    data[columns[-1]] = []
                while q:
                    s = q[0]
                    if s[0].startswith(("loop_", "_")):
                        break
                    items.append("".join(q.popleft()))
                n = len(items) // len(columns)
                assert len(items) % n == 0
                loops.append(columns)
                for k, v in zip(columns * n, items):
                    data[k].append(v.strip())
            elif issue := "".join(s).strip():
                warnings.warn(f"Possible issue in CIF file at line: {issue}")
        return cls(data, loops, header)


class CifFile:
    """Reads and parses CifBlocks from a .cif file or string."""

    def __init__(self, data, orig_string=None, comment=None):
        """
        Args:
            data (dict): Of CifBlock objects.
            orig_string (str): The original cif string.
            comment (str): Comment string.
        """
        self.data = data
        self.orig_string = orig_string
        self.comment = comment or "# generated using pymatgen"

    def __str__(self):
        out = "\n".join(map(str, self.data.values()))
        return f"{self.comment}\n{out}\n"

    @classmethod
    @np.deprecate(message="Use from_str instead")
    def from_string(cls, *args, **kwargs):
        return cls.from_str(*args, **kwargs)

    @classmethod
    def from_str(cls, string):
        """
        Reads CifFile from a string.

        :param string: String representation.

        Returns:
            CifFile
        """
        dct = {}
        for x in re.split(r"^\s*data_", f"x\n{string}", flags=re.MULTILINE | re.DOTALL)[1:]:
            # Skip over Cif block that contains powder diffraction data.
            # Some elements in this block were missing from CIF files in
            # Springer materials/Pauling file DBs.
            # This block anyway does not contain any structure information, and
            # CifParser was also not parsing it.
            if "powder_pattern" in re.split(r"\n", x, maxsplit=1)[0]:
                continue
            c = CifBlock.from_str("data_" + x)
            dct[c.header] = c
        return cls(dct, string)

    @classmethod
    def from_file(cls, filename):
        """
        Reads CifFile from a filename.

        :param filename: Filename

        Returns:
            CifFile
        """
        with zopen(str(filename), "rt", errors="replace") as f:
            return cls.from_str(f.read())


class CifParser:
    """
    Parses a CIF file. Attempts to fix CIFs that are out-of-spec, but will
    issue warnings if corrections applied. These are also stored in the
    CifParser's errors attribute.
    """

    def __init__(
        self,
        filename: str | StringIO,
        occupancy_tolerance: float = 1.0,
        site_tolerance: float = 1e-4,
        frac_tolerance: float = 1e-4,
    ) -> None:
        """
        Args:
            filename (str): CIF filename, gzipped or bzipped CIF files are fine too.
            occupancy_tolerance (float): If total occupancy of a site is between 1 and occupancy_tolerance, the
                occupancies will be scaled down to 1.
            site_tolerance (float): This tolerance is used to determine if two sites are sitting in the same position,
                in which case they will be combined to a single disordered site. Defaults to 1e-4.
            frac_tolerance (float): This tolerance is used to determine is a coordinate should be rounded to an ideal
                value. E.g., 0.6667 is rounded to 2/3. This is desired if symmetry operations are going to be applied.
                However, for very large CIF files, this may need to be set to 0.
        """
        self._occupancy_tolerance = occupancy_tolerance
        self._site_tolerance = site_tolerance
        self._frac_tolerance = frac_tolerance
        if isinstance(filename, (str, Path)):
            self._cif = CifFile.from_file(filename)
        else:
            self._cif = CifFile.from_str(filename.read())

        # store if CIF contains features from non-core CIF dictionaries
        # e.g. magCIF
        self.feature_flags = {}
        self.warnings: list[str] = []

        def is_magcif() -> bool:
            """Checks to see if file appears to be a magCIF file (heuristic)."""
            # Doesn't seem to be a canonical way to test if file is magCIF or
            # not, so instead check for magnetic symmetry datanames
            prefixes = ["_space_group_magn", "_atom_site_moment", "_space_group_symop_magn"]
            for d in self._cif.data.values():
                for k in d.data:
                    for prefix in prefixes:
                        if prefix in k:
                            return True
            return False

        self.feature_flags["magcif"] = is_magcif()

        def is_magcif_incommensurate() -> bool:
            """
            Checks to see if file contains an incommensurate magnetic
            structure (heuristic).
            """
            # Doesn't seem to be a canonical way to test if magCIF file
            # describes incommensurate structure or not, so instead check
            # for common datanames
            if not self.feature_flags["magcif"]:
                return False
            prefixes = ["_cell_modulation_dimension", "_cell_wave_vector"]
            for d in self._cif.data.values():
                for k in d.data:
                    for prefix in prefixes:
                        if prefix in k:
                            return True
            return False

        self.feature_flags["magcif_incommensurate"] = is_magcif_incommensurate()

        for k in self._cif.data:
            # pass individual CifBlocks to _sanitize_data
            self._cif.data[k] = self._sanitize_data(self._cif.data[k])

    @classmethod
    @np.deprecate(message="Use from_str instead")
    def from_string(cls, *args, **kwargs):
        return cls.from_str(*args, **kwargs)

    @staticmethod
    def from_str(cif_string: str, **kwargs) -> CifParser:
        """
        Creates a CifParser from a string.

        Args:
            cif_string (str): String representation of a CIF.
            **kwargs: Passthrough of all kwargs supported by CifParser.

        Returns:
            CifParser
        """
        stream = StringIO(cif_string)
        return CifParser(stream, **kwargs)

    def _sanitize_data(self, data):
        """
        Some CIF files do not conform to spec. This function corrects
        known issues, particular in regards to Springer materials/
        Pauling files.

        This function is here so that CifParser can assume its
        input conforms to spec, simplifying its implementation.
        :param data: CifBlock

        Returns:
            data CifBlock
        """
        """
        This part of the code deals with handling formats of data as found in
        CIF files extracted from the Springer Materials/Pauling File
        databases, and that are different from standard ICSD formats.
        """
        # check for implicit hydrogens, warn if any present
        if "_atom_site_attached_hydrogens" in data.data:
            attached_hydrogens = [str2float(x) for x in data.data["_atom_site_attached_hydrogens"] if str2float(x) != 0]
            if len(attached_hydrogens) > 0:
                self.warnings.append(
                    "Structure has implicit hydrogens defined, parsed structure unlikely to be "
                    "suitable for use in calculations unless hydrogens added."
                )

        # Check to see if "_atom_site_type_symbol" exists, as some test CIFs do
        # not contain this key.
        if "_atom_site_type_symbol" in data.data:
            # Keep a track of which data row needs to be removed.
            # Example of a row: Nb,Zr '0.8Nb + 0.2Zr' .2a .m-3m 0 0 0 1 14
            # 'rhombic dodecahedron, Nb<sub>14</sub>'
            # Without this code, the above row in a structure would be parsed
            # as an ordered site with only Nb (since
            # CifParser would try to parse the first two characters of the
            # label "Nb,Zr") and occupancy=1.
            # However, this site is meant to be a disordered site with 0.8 of
            # Nb and 0.2 of Zr.
            idxs_to_remove = []

            new_atom_site_label = []
            new_atom_site_type_symbol = []
            new_atom_site_occupancy = []
            new_fract_x = []
            new_fract_y = []
            new_fract_z = []

            for idx, el_row in enumerate(data["_atom_site_label"]):
                # CIF files from the Springer Materials/Pauling File have
                # switched the label and symbol. Thus, in the
                # above shown example row, '0.8Nb + 0.2Zr' is the symbol.
                # Below, we split the strings on ' + ' to
                # check if the length (or number of elements) in the label and
                # symbol are equal.
                if len(data["_atom_site_type_symbol"][idx].split(" + ")) > len(el_row.split(" + ")):
                    # Dictionary to hold extracted elements and occupancies
                    els_occu = {}

                    # parse symbol to get element names and occupancy and store
                    # in "els_occu"
                    symbol_str = data["_atom_site_type_symbol"][idx]
                    symbol_str_lst = symbol_str.split(" + ")
                    for elocc_idx, sym in enumerate(symbol_str_lst):
                        # Remove any bracketed items in the string
                        symbol_str_lst[elocc_idx] = re.sub(r"\([0-9]*\)", "", sym.strip())

                        # Extract element name and its occupancy from the
                        # string, and store it as a
                        # key-value pair in "els_occ".
                        els_occu[
                            str(re.findall(r"\D+", symbol_str_lst[elocc_idx].strip())[1]).replace("<sup>", "")
                        ] = float("0" + re.findall(r"\.?\d+", symbol_str_lst[elocc_idx].strip())[1])

                    x = str2float(data["_atom_site_fract_x"][idx])
                    y = str2float(data["_atom_site_fract_y"][idx])
                    z = str2float(data["_atom_site_fract_z"][idx])

                    for et, occu in els_occu.items():
                        # new atom site labels have 'fix' appended
                        new_atom_site_label.append(f"{et}_fix{len(new_atom_site_label)}")
                        new_atom_site_type_symbol.append(et)
                        new_atom_site_occupancy.append(str(occu))
                        new_fract_x.append(str(x))
                        new_fract_y.append(str(y))
                        new_fract_z.append(str(z))

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

            if len(idxs_to_remove) > 0:
                self.warnings.append("Pauling file corrections applied.")

                data.data["_atom_site_label"] += new_atom_site_label
                data.data["_atom_site_type_symbol"] += new_atom_site_type_symbol
                data.data["_atom_site_occupancy"] += new_atom_site_occupancy
                data.data["_atom_site_fract_x"] += new_fract_x
                data.data["_atom_site_fract_y"] += new_fract_y
                data.data["_atom_site_fract_z"] += new_fract_z
        """
        This fixes inconsistencies in naming of several magCIF tags
        as a result of magCIF being in widespread use prior to
        specification being finalized (on advice of Branton Campbell).
        """
        if self.feature_flags["magcif"]:
            # CIF-1 style has all underscores, interim standard
            # had period before magn instead of before the final
            # component (e.g. xyz)
            # we want to standardize on a specific key, to simplify
            # parsing code
            correct_keys = [
                "_space_group_symop_magn_operation.xyz",
                "_space_group_symop_magn_centering.xyz",
                "_space_group_magn.name_BNS",
                "_space_group_magn.number_BNS",
                "_atom_site_moment_crystalaxis_x",
                "_atom_site_moment_crystalaxis_y",
                "_atom_site_moment_crystalaxis_z",
                "_atom_site_moment_label",
            ]

            # cannot mutate dict during enumeration, so store changes we want to make
            changes_to_make = {}

            for original_key in data.data:
                for correct_key in correct_keys:
                    # convert to all underscore
                    trial_key = "_".join(correct_key.split("."))
                    test_key = "_".join(original_key.split("."))
                    if trial_key == test_key:
                        changes_to_make[correct_key] = original_key

            # make changes
            for correct_key, original_key in changes_to_make.items():
                data.data[correct_key] = data.data[original_key]

            # renamed_keys maps interim_keys to final_keys
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

        # check for finite precision frac coordinates (e.g. 0.6667 instead of 0.6666666...7)
        # this can sometimes cause serious issues when applying symmetry operations
        important_fracs = (1 / 3, 2 / 3)
        fracs_to_change = {}
        for label in ("_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"):
            if label in data.data:
                for idx, frac in enumerate(data.data[label]):
                    try:
                        frac = str2float(frac)
                    except Exception:
                        # coordinate might not be defined e.g. '?'
                        continue
                    for comparison_frac in important_fracs:
                        if abs(1 - frac / comparison_frac) < self._frac_tolerance:
                            fracs_to_change[(label, idx)] = str(comparison_frac)
        if fracs_to_change:
            self.warnings.append(
                "Some fractional coordinates rounded to ideal values to avoid issues with finite precision."
            )
            for (label, idx), val in fracs_to_change.items():
                data.data[label][idx] = val

        return data

    def _unique_coords(
        self,
        coords: list[Vector3D],
        magmoms: list[Magmom] | None = None,
        lattice: Lattice | None = None,
        labels: dict[Vector3D, str] | None = None,
    ):
        """
        Generate unique coordinates using coord and symmetry positions
        and also their corresponding magnetic moments, if supplied.
        """
        coords_out: list[np.ndarray] = []
        labels_out = []
        labels = labels or {}

        if magmoms:
            magmoms_out = []
            if len(magmoms) != len(coords):
                raise ValueError
            for tmp_coord, tmp_magmom in zip(coords, magmoms):
                for op in self.symmetry_operations:
                    coord = op.operate(tmp_coord)
                    coord = np.array([i - math.floor(i) for i in coord])
                    if isinstance(op, MagSymmOp):
                        # Up to this point, magmoms have been defined relative
                        # to crystal axis. Now convert to Cartesian and into
                        # a Magmom object.
                        magmom = Magmom.from_moment_relative_to_crystal_axes(
                            op.operate_magmom(tmp_magmom), lattice=lattice
                        )
                    else:
                        magmom = Magmom(tmp_magmom)
                    if not in_coord_list_pbc(coords_out, coord, atol=self._site_tolerance):
                        coords_out.append(coord)
                        magmoms_out.append(magmom)
                        labels_out.append(labels.get(tmp_coord))
            return coords_out, magmoms_out, labels_out

        for tmp_coord in coords:
            for op in self.symmetry_operations:
                coord = op.operate(tmp_coord)
                coord = np.array([i - math.floor(i) for i in coord])
                if not in_coord_list_pbc(coords_out, coord, atol=self._site_tolerance):
                    coords_out.append(coord)
                    labels_out.append(labels.get(tmp_coord))

        dummy_magmoms = [Magmom(0)] * len(coords_out)
        return coords_out, dummy_magmoms, labels_out

    def get_lattice(
        self,
        data,
        length_strings=("a", "b", "c"),
        angle_strings=("alpha", "beta", "gamma"),
        lattice_type=None,
    ):
        """
        Generate the lattice from the provided lattice parameters. In
        the absence of all six lattice parameters, the crystal system
        and necessary parameters are parsed.
        """
        try:
            return self.get_lattice_no_exception(
                data=data, angle_strings=angle_strings, lattice_type=lattice_type, length_strings=length_strings
            )

        except KeyError:
            # Missing Key search for cell setting
            for lattice_label in ["_symmetry_cell_setting", "_space_group_crystal_system"]:
                if data.data.get(lattice_label):
                    lattice_type = data.data.get(lattice_label).lower()
                    try:
                        required_args = getargspec(getattr(Lattice, lattice_type)).args

                        lengths = (length for length in length_strings if length in required_args)
                        angles = (a for a in angle_strings if a in required_args)
                        return self.get_lattice(data, lengths, angles, lattice_type=lattice_type)
                    except AttributeError as exc:
                        self.warnings.append(str(exc))
                        warnings.warn(exc)

                else:
                    return None
        return None

    @staticmethod
    def get_lattice_no_exception(
        data, length_strings=("a", "b", "c"), angle_strings=("alpha", "beta", "gamma"), lattice_type=None
    ):
        """
        Take a dictionary of CIF data and returns a pymatgen Lattice object.

        Args:
            data: a dictionary of the CIF file
            length_strings: The strings that are used to identify the length parameters in the CIF file.
            angle_strings: The strings that are used to identify the angles in the CIF file.
            lattice_type: The type of lattice.  This is a string, and can be any of the following:

        Returns:
            Lattice object
        """
        lengths = [str2float(data["_cell_length_" + i]) for i in length_strings]
        angles = [str2float(data["_cell_angle_" + i]) for i in angle_strings]
        if not lattice_type:
            return Lattice.from_parameters(*lengths, *angles)
        return getattr(Lattice, lattice_type)(*(lengths + angles))

    def get_symops(self, data):
        """
        In order to generate symmetry equivalent positions, the symmetry
        operations are parsed. If the symops are not present, the space
        group symbol is parsed, and symops are generated.
        """
        symops = []
        for symmetry_label in [
            "_symmetry_equiv_pos_as_xyz",
            "_symmetry_equiv_pos_as_xyz_",
            "_space_group_symop_operation_xyz",
            "_space_group_symop_operation_xyz_",
        ]:
            if data.data.get(symmetry_label):
                xyz = data.data.get(symmetry_label)
                if isinstance(xyz, str):
                    msg = "A 1-line symmetry op P1 CIF is detected!"
                    warnings.warn(msg)
                    self.warnings.append(msg)
                    xyz = [xyz]
                try:
                    symops = [SymmOp.from_xyz_str(s) for s in xyz]
                    break
                except ValueError:
                    continue
        if not symops:
            # Try to parse symbol
            for symmetry_label in [
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
            ]:
                sg = data.data.get(symmetry_label)
                msg_template = "No _symmetry_equiv_pos_as_xyz type key found. Spacegroup from {} used."

                if sg:
                    sg = sub_spgrp(sg)
                    try:
                        spg = space_groups.get(sg)
                        if spg:
                            symops = SpaceGroup(spg).symmetry_ops
                            msg = msg_template.format(symmetry_label)
                            warnings.warn(msg)
                            self.warnings.append(msg)
                            break
                    except ValueError:
                        # Ignore any errors
                        pass

                    try:
                        cod_data = loadfn(
                            os.path.join(os.path.dirname(os.path.dirname(__file__)), "symmetry", "symm_ops.json")
                        )
                        for d in cod_data:
                            if sg == re.sub(r"\s+", "", d["hermann_mauguin"]):
                                xyz = d["symops"]
                                symops = [SymmOp.from_xyz_str(s) for s in xyz]
                                msg = msg_template.format(symmetry_label)
                                warnings.warn(msg)
                                self.warnings.append(msg)
                                break
                    except Exception:
                        continue

                    if symops:
                        break
        if not symops:
            # Try to parse International number
            for symmetry_label in [
                "_space_group_IT_number",
                "_space_group_IT_number_",
                "_symmetry_Int_Tables_number",
                "_symmetry_Int_Tables_number_",
            ]:
                if data.data.get(symmetry_label):
                    try:
                        i = int(str2float(data.data.get(symmetry_label)))
                        symops = SpaceGroup.from_int_number(i).symmetry_ops
                        break
                    except ValueError:
                        continue

        if not symops:
            msg = "No _symmetry_equiv_pos_as_xyz type key found. Defaulting to P1."
            warnings.warn(msg)
            self.warnings.append(msg)
            symops = [SymmOp.from_xyz_str(s) for s in ["x", "y", "z"]]

        return symops

    def get_magsymops(self, data):
        """
        Equivalent to get_symops except for magnetic symmetry groups.
        Separate function since additional operation for time reversal symmetry
        (which changes magnetic moments on sites) needs to be returned.
        """
        mag_symm_ops = []
        bns_name = data.data.get("_space_group_magn.name_BNS")  # get BNS label for MagneticSpaceGroup()
        bns_num = data.data.get("_space_group_magn.number_BNS")  # get BNS number for MagneticSpaceGroup()

        # check to see if magCIF file explicitly contains magnetic symmetry operations
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
                                time_reversal=new_time_reversal,
                            )
                        )
                mag_symm_ops = all_ops

        # else check to see if it specifies a magnetic space group
        elif bns_name or bns_num:
            label = bns_name if bns_name else list(map(int, (bns_num.split("."))))

            if data.data.get("_space_group_magn.transform_BNS_Pp_abc") != "a,b,c;0,0,0":
                jf = data.data.get("_space_group_magn.transform_BNS_Pp_abc")
                msg = MagneticSpaceGroup(label, jf)

            elif data.data.get("_space_group_magn.transform_BNS_Pp"):
                return NotImplementedError("Incomplete specification to implement.")
            else:
                msg = MagneticSpaceGroup(label)

            mag_symm_ops = msg.symmetry_ops

        if not mag_symm_ops:
            msg = "No magnetic symmetry detected, using primitive symmetry."
            warnings.warn(msg)
            self.warnings.append(msg)
            mag_symm_ops = [MagSymmOp.from_xyzt_str("x, y, z, 1")]

        return mag_symm_ops

    @staticmethod
    def parse_oxi_states(data):
        """Parse oxidation states from data dictionary."""
        try:
            oxi_states = {
                data["_atom_type_symbol"][i]: str2float(data["_atom_type_oxidation_number"][i])
                for i in range(len(data["_atom_type_symbol"]))
            }
            # attempt to strip oxidation state from _atom_type_symbol
            # in case the label does not contain an oxidation state
            for i, symbol in enumerate(data["_atom_type_symbol"]):
                oxi_states[re.sub(r"\d?[\+,\-]?$", "", symbol)] = str2float(data["_atom_type_oxidation_number"][i])

        except (ValueError, KeyError):
            oxi_states = None
        return oxi_states

    @staticmethod
    def parse_magmoms(data, lattice=None):
        """Parse atomic magnetic moments from data dictionary."""
        if lattice is None:
            raise Exception("Magmoms given in terms of crystal axes in magCIF spec.")
        try:
            magmoms = {
                data["_atom_site_moment_label"][i]: np.array(
                    [
                        str2float(data["_atom_site_moment_crystalaxis_x"][i]),
                        str2float(data["_atom_site_moment_crystalaxis_y"][i]),
                        str2float(data["_atom_site_moment_crystalaxis_z"][i]),
                    ]
                )
                for i in range(len(data["_atom_site_moment_label"]))
            }
        except (ValueError, KeyError):
            return None
        return magmoms

    def _parse_symbol(self, sym):
        """
        Parse a string with a symbol to extract a string representing an element.

        Args:
            sym (str): A symbol to be parsed.

        Returns:
            A string with the parsed symbol. None if no parsing was possible.
        """
        # Common representations for elements/water in cif files
        # TODO: fix inconsistent handling of water
        special = {
            "Hw": "H",
            "Ow": "O",
            "Wat": "O",
            "wat": "O",
            "OH": "",
            "OH2": "",
            "NO3": "N",
        }

        parsed_sym = None
        # try with special symbols, otherwise check the first two letters,
        # then the first letter alone. If everything fails try extracting the
        # first letters.
        m_sp = re.match("|".join(special), sym)
        if m_sp:
            parsed_sym = special[m_sp.group()]
        elif Element.is_valid_symbol(sym[:2].title()):
            parsed_sym = sym[:2].title()
        elif Element.is_valid_symbol(sym[0].upper()):
            parsed_sym = sym[0].upper()
        else:
            m = re.match(r"w?[A-Z][a-z]*", sym)
            if m:
                parsed_sym = m.group()

        if parsed_sym is not None and (m_sp or not re.match(rf"{parsed_sym}\d*", sym)):
            msg = f"{sym} parsed as {parsed_sym}"
            warnings.warn(msg)
            self.warnings.append(msg)

        return parsed_sym

    def _get_structure(
        self, data: dict[str, Any], primitive: bool, symmetrized: bool, check_occu: bool = False
    ) -> Structure | None:
        """Generate structure from part of the cif."""

        def get_num_implicit_hydrogens(sym):
            num_h = {"Wat": 2, "wat": 2, "O-H": 1}
            return num_h.get(sym[:3], 0)

        lattice = self.get_lattice(data)

        # if magCIF, get magnetic symmetry moments and magmoms
        # else standard CIF, and use empty magmom dict
        if self.feature_flags["magcif_incommensurate"]:
            raise NotImplementedError("Incommensurate structures not currently supported.")
        if self.feature_flags["magcif"]:
            self.symmetry_operations = self.get_magsymops(data)
            magmoms = self.parse_magmoms(data, lattice=lattice)
        else:
            self.symmetry_operations = self.get_symops(data)
            magmoms = {}

        oxi_states = self.parse_oxi_states(data)

        coord_to_species = {}  # type: ignore
        coord_to_magmoms = {}
        labels = {}

        def get_matching_coord(coord):
            keys = list(coord_to_species)
            coords = np.array(keys)
            for op in self.symmetry_operations:
                frac_coord = op.operate(coord)
                indices = find_in_coord_list_pbc(coords, frac_coord, atol=self._site_tolerance)
                if len(indices) > 0:
                    return keys[indices[0]]
            return False

        for idx, label in enumerate(data["_atom_site_label"]):
            try:
                # If site type symbol exists, use it. Otherwise, we use the label.
                symbol = self._parse_symbol(data["_atom_site_type_symbol"][idx])
                num_h = get_num_implicit_hydrogens(data["_atom_site_type_symbol"][idx])
            except KeyError:
                symbol = self._parse_symbol(label)
                num_h = get_num_implicit_hydrogens(label)
            if not symbol:
                continue

            if oxi_states is not None:
                o_s = oxi_states.get(symbol, 0)
                # use _atom_site_type_symbol if possible for oxidation state
                if "_atom_site_type_symbol" in data.data:  # type: ignore[attr-defined]
                    oxi_symbol = data["_atom_site_type_symbol"][idx]
                    o_s = oxi_states.get(oxi_symbol, o_s)
                try:
                    el = Species(symbol, o_s)
                except Exception:
                    el = DummySpecies(symbol, o_s)
            else:
                el = get_el_sp(symbol)  # type: ignore

            x = str2float(data["_atom_site_fract_x"][idx])
            y = str2float(data["_atom_site_fract_y"][idx])
            z = str2float(data["_atom_site_fract_z"][idx])
            magmom = magmoms.get(label, np.array([0, 0, 0]))

            try:
                occu = str2float(data["_atom_site_occupancy"][idx])
            except (KeyError, ValueError):
                occu = 1
            # If check_occu is True or the occupancy is greater than 0, create comp_d
            if not check_occu or occu > 0:
                coord = (x, y, z)
                match = get_matching_coord(coord)
                comp_dict = {el: max(occu, 1e-8)}

                if num_h > 0:
                    comp_dict["H"] = num_h  # type: ignore
                    self.warnings.append(
                        "Structure has implicit hydrogens defined, parsed structure unlikely to be "
                        "suitable for use in calculations unless hydrogens added."
                    )
                comp = Composition(comp_dict)

                if not match:
                    coord_to_species[coord] = comp
                    coord_to_magmoms[coord] = magmom
                    labels[coord] = label
                else:
                    coord_to_species[match] += comp
                    # disordered magnetic not currently supported
                    coord_to_magmoms[match] = None
                    labels[match] = label
        sum_occu = [
            sum(c.values()) for c in coord_to_species.values() if set(c.elements) != {Element("O"), Element("H")}
        ]
        if any(occu > 1 for occu in sum_occu):
            msg = (
                f"Some occupancies ({sum_occu}) sum to > 1! If they are within "
                "the occupancy_tolerance, they will be rescaled. "
                f"The current occupancy_tolerance is set to: {self._occupancy_tolerance}"
            )
            warnings.warn(msg)
            self.warnings.append(msg)

        all_species = []
        all_coords = []
        all_magmoms = []
        all_hydrogens = []
        equivalent_indices = []
        all_labels = []

        # check to see if magCIF file is disordered
        if self.feature_flags["magcif"]:
            for v in coord_to_magmoms.values():
                if v is None:
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
                tmp_coords = [site[0] for site in group]
                tmp_magmom = [coord_to_magmoms[tmp_coord] for tmp_coord in tmp_coords]

                if self.feature_flags["magcif"]:
                    coords, magmoms, new_labels = self._unique_coords(
                        tmp_coords, magmoms=tmp_magmom, labels=labels, lattice=lattice
                    )
                else:
                    coords, magmoms, new_labels = self._unique_coords(tmp_coords, labels=labels)

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
                all_coords.extend(coords)
                all_species.extend(len(coords) * [species])
                all_magmoms.extend(magmoms)
                all_labels.extend(new_labels)

            # rescale occupancies if necessary
            all_species_noedit = all_species[:]  # save copy before scaling in case of check_occu=False, used below
            for idx, species in enumerate(all_species):
                total_occu = sum(species.values())
                if 1 < total_occu <= self._occupancy_tolerance:
                    all_species[idx] = species / total_occu

        if all_species and len(all_species) == len(all_coords) and len(all_species) == len(all_magmoms):
            site_properties = {}
            if any(all_hydrogens):
                assert len(all_hydrogens) == len(all_coords)
                site_properties["implicit_hydrogens"] = all_hydrogens

            if self.feature_flags["magcif"]:
                site_properties["magmom"] = all_magmoms

            if len(site_properties) == 0:
                site_properties = None  # type: ignore

            if any(all_labels):
                assert len(all_labels) == len(all_species)
            else:
                all_labels = None  # type: ignore

            struct = Structure(lattice, all_species, all_coords, site_properties=site_properties, labels=all_labels)

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
                for idx in range(len(struct)):
                    struct[idx] = PeriodicSite(
                        all_species_noedit[idx], all_coords[idx], lattice, properties=site_properties, skip_checks=True
                    )

            if symmetrized or not check_occu:
                return struct

            struct = struct.get_sorted_structure()

            if primitive and self.feature_flags["magcif"]:
                struct = struct.get_primitive_structure(use_site_props=True)
            elif primitive:
                struct = struct.get_primitive_structure()
                struct = struct.get_reduced_structure()

            return struct
        return None

    def get_structures(
        self,
        primitive: bool = True,
        symmetrized: bool = False,
        check_occu: bool = True,
        on_error: Literal["ignore", "warn", "raise"] = "warn",
    ) -> list[Structure]:
        """Return list of structures in CIF file.

        Args:
            primitive (bool): Set to False to return conventional unit cells.
                Defaults to True. With magnetic CIF files, will return primitive
                magnetic cell which may be larger than nuclear primitive cell.
            symmetrized (bool): If True, return a SymmetrizedStructure which will
                include the equivalent indices and symmetry operations used to
                create the Structure as provided by the CIF (if explicit symmetry
                operations are included in the CIF) or generated from information
                in the CIF (if only space group labels are provided). Note that
                currently Wyckoff labels and space group labels or numbers are
                not included in the generated SymmetrizedStructure, these will be
                notated as "Not Parsed" or -1 respectively.
            check_occu (bool): If False, site occupancy will not be checked, allowing unphysical
                occupancy != 1. Useful for experimental results in which occupancy was allowed
                to refine to unphysical values. Warning: unphysical site occupancies are incompatible
                with many pymatgen features. Defaults to True.
            on_error ('ignore' | 'warn' | 'raise'): What to do in case of KeyError or ValueError
                while parsing CIF file. Defaults to 'warn'.

        Returns:
            list[Structure]: All structures in CIF file.
        """
        if not check_occu:  # added in https://github.com/materialsproject/pymatgen/pull/2836
            warnings.warn("Structures with unphysical site occupancies are not compatible with many pymatgen features.")
        if primitive and symmetrized:
            raise ValueError(
                "Using both 'primitive' and 'symmetrized' arguments is not currently supported "
                "since unexpected behavior might result."
            )

        structures = []
        for idx, dct in enumerate(self._cif.data.values()):
            try:
                struct = self._get_structure(dct, primitive, symmetrized, check_occu=check_occu)
                if struct:
                    structures.append(struct)
            except (KeyError, ValueError) as exc:
                # A user reported a problem with cif files produced by Avogadro
                # in which the atomic coordinates are in Cartesian coords.
                msg = f"No structure parsed for section {idx + 1} in CIF.\n{exc}"
                if on_error == "raise":
                    raise ValueError(msg) from exc
                if on_error == "warn":
                    warnings.warn(msg)
                self.warnings.append(msg)
                # continue silently if on_error == "ignore"

        # if on_error == "raise" we don't get to here so no need to check
        if self.warnings and on_error == "warn":
            warnings.warn("Issues encountered while parsing CIF: " + "\n".join(self.warnings))

        if len(structures) == 0:
            raise ValueError("Invalid CIF file with no structures!")
        return structures

    def get_bibtex_string(self):
        """
        Get BibTeX reference from CIF file.
        :param data:

        Returns:
            BibTeX string.
        """
        try:
            from pybtex.database import BibliographyData, Entry
        except ImportError:
            raise RuntimeError("Bibliographic data extraction requires pybtex.")

        bibtex_keys = {
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

        entries = {}

        # TODO: parse '_publ_section_references' when it exists?
        # TODO: CIF specification supports multiple citations.

        for idx, data in enumerate(self._cif.data.values()):
            # convert to lower-case keys, some cif files inconsistent
            data = {k.lower(): v for k, v in data.data.items()}

            bibtex_entry = {}

            for field, tags in bibtex_keys.items():
                for tag in tags:
                    if tag in data:
                        if isinstance(data[tag], list):
                            bibtex_entry[field] = data[tag][0]
                        else:
                            bibtex_entry[field] = data[tag]

            # convert to bibtex author format ('and' delimited)
            if "author" in bibtex_entry:
                # separate out semicolon authors
                if isinstance(bibtex_entry["author"], str) and ";" in bibtex_entry["author"]:
                    bibtex_entry["author"] = bibtex_entry["author"].split(";")

                if isinstance(bibtex_entry["author"], list):
                    bibtex_entry["author"] = " and ".join(bibtex_entry["author"])

            # convert to bibtex page range format, use empty string if not specified
            if ("page_first" in bibtex_entry) or ("page_last" in bibtex_entry):
                bibtex_entry["pages"] = bibtex_entry.get("page_first", "") + "--" + bibtex_entry.get("page_last", "")
                bibtex_entry.pop("page_first", None)  # and remove page_first, page_list if present
                bibtex_entry.pop("page_last", None)

            # cite keys are given as cif-reference-idx in order they are found
            entries[f"cifref{idx}"] = Entry("article", list(bibtex_entry.items()))

        return BibliographyData(entries).to_string(bib_format="bibtex")

    def as_dict(self):
        """MSONable dict"""
        dct = {}
        for k, v in self._cif.data.items():
            dct[k] = {}
            for k2, v2 in v.data.items():
                dct[k][k2] = v2
        return dct

    @property
    def has_errors(self):
        """Whether there are errors/warnings detected in CIF parsing."""
        return len(self.warnings) > 0


class CifWriter:
    """A wrapper around CifFile to write CIF files from pymatgen structures."""

    def __init__(
        self,
        struct,
        symprec=None,
        write_magmoms=False,
        significant_figures=8,
        angle_tolerance=5.0,
        refine_struct=True,
    ):
        """
        Args:
            struct (Structure): structure to write
            symprec (float): If not none, finds the symmetry of the structure
                and writes the cif with symmetry information. Passes symprec
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
        """
        if write_magmoms and symprec:
            warnings.warn("Magnetic symmetry cannot currently be detected by pymatgen,disabling symmetry detection.")
            symprec = None

        format_str = f"{{:.{significant_figures}f}}"

        block = {}
        loops = []
        spacegroup = ("P 1", 1)
        if symprec is not None:
            sf = SpacegroupAnalyzer(struct, symprec, angle_tolerance=angle_tolerance)
            spacegroup = (sf.get_space_group_symbol(), sf.get_space_group_number())

            if refine_struct:
                # Needs the refined structure when using symprec. This converts
                # primitive to conventional structures, the standard for CIF.
                struct = sf.get_refined_structure()

        latt = struct.lattice
        comp = struct.composition
        no_oxi_comp = comp.element_composition
        block["_symmetry_space_group_name_H-M"] = spacegroup[0]
        for cell_attr in ["a", "b", "c"]:
            block["_cell_length_" + cell_attr] = format_str.format(getattr(latt, cell_attr))
        for cell_attr in ["alpha", "beta", "gamma"]:
            block["_cell_angle_" + cell_attr] = format_str.format(getattr(latt, cell_attr))
        block["_symmetry_Int_Tables_number"] = spacegroup[1]
        block["_chemical_formula_structural"] = no_oxi_comp.reduced_formula
        block["_chemical_formula_sum"] = no_oxi_comp.formula
        block["_cell_volume"] = format_str.format(latt.volume)

        reduced_comp, fu = no_oxi_comp.get_reduced_composition_and_factor()
        block["_cell_formula_units_Z"] = str(int(fu))

        if symprec is None:
            block["_symmetry_equiv_pos_site_id"] = ["1"]
            block["_symmetry_equiv_pos_as_xyz"] = ["x, y, z"]
        else:
            sf = SpacegroupAnalyzer(struct, symprec)

            symmops = []
            for op in sf.get_symmetry_operations():
                v = op.translation_vector
                symmops.append(SymmOp.from_rotation_and_translation(op.rotation_matrix, v))

            ops = [op.as_xyz_string() for op in symmops]
            block["_symmetry_equiv_pos_site_id"] = [f"{i}" for i in range(1, len(ops) + 1)]
            block["_symmetry_equiv_pos_as_xyz"] = ops

        loops.append(["_symmetry_equiv_pos_site_id", "_symmetry_equiv_pos_as_xyz"])

        try:
            symbol_to_oxinum = {str(el): float(el.oxi_state) for el in sorted(comp.elements)}
            block["_atom_type_symbol"] = list(symbol_to_oxinum)
            block["_atom_type_oxidation_number"] = symbol_to_oxinum.values()
            loops.append(["_atom_type_symbol", "_atom_type_oxidation_number"])
        except (TypeError, AttributeError):
            symbol_to_oxinum = {el.symbol: 0 for el in sorted(comp.elements)}

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
                        moment = Magmom.get_moment_relative_to_crystal_axes(magmom, latt)
                        atom_site_moment_label.append(f"{sp.symbol}{count}")
                        atom_site_moment_crystalaxis_x.append(format_str.format(moment[0]))
                        atom_site_moment_crystalaxis_y.append(format_str.format(moment[1]))
                        atom_site_moment_crystalaxis_z.append(format_str.format(moment[2]))

                    count += 1
        else:
            # The following just presents a deterministic ordering.
            unique_sites = [
                (
                    sorted(sites, key=lambda s: tuple(abs(x) for x in s.frac_coords))[0],
                    len(sites),
                )
                for sites in sf.get_symmetrized_structure().equivalent_sites
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

        block["_atom_site_type_symbol"] = atom_site_type_symbol
        block["_atom_site_label"] = atom_site_label
        block["_atom_site_symmetry_multiplicity"] = atom_site_symmetry_multiplicity
        block["_atom_site_fract_x"] = atom_site_fract_x
        block["_atom_site_fract_y"] = atom_site_fract_y
        block["_atom_site_fract_z"] = atom_site_fract_z
        block["_atom_site_occupancy"] = atom_site_occupancy
        loops.append(
            [
                "_atom_site_type_symbol",
                "_atom_site_label",
                "_atom_site_symmetry_multiplicity",
                "_atom_site_fract_x",
                "_atom_site_fract_y",
                "_atom_site_fract_z",
                "_atom_site_occupancy",
            ]
        )
        if write_magmoms:
            block["_atom_site_moment_label"] = atom_site_moment_label
            block["_atom_site_moment_crystalaxis_x"] = atom_site_moment_crystalaxis_x
            block["_atom_site_moment_crystalaxis_y"] = atom_site_moment_crystalaxis_y
            block["_atom_site_moment_crystalaxis_z"] = atom_site_moment_crystalaxis_z
            loops.append(
                [
                    "_atom_site_moment_label",
                    "_atom_site_moment_crystalaxis_x",
                    "_atom_site_moment_crystalaxis_y",
                    "_atom_site_moment_crystalaxis_z",
                ]
            )
        dct = {}
        dct[comp.reduced_formula] = CifBlock(block, loops, comp.reduced_formula)
        self._cf = CifFile(dct)

    @property
    def ciffile(self):
        """Returns: CifFile associated with the CifWriter."""
        return self._cf

    def __str__(self):
        """Returns the cif as a string."""
        return str(self._cf)

    def write_file(self, filename):
        """Write the cif file."""
        with zopen(filename, "wt") as f:
            f.write(str(self))


def str2float(text):
    """Remove uncertainty brackets from strings and return the float."""
    try:
        # Note that the ending ) is sometimes missing. That is why the code has
        # been modified to treat it as optional. Same logic applies to lists.
        return float(re.sub(r"\(.+\)*", "", text))
    except TypeError:
        if isinstance(text, list) and len(text) == 1:
            return float(re.sub(r"\(.+\)*", "", text[0]))
    except ValueError as exc:
        if text.strip() == ".":
            return 0
        raise exc
    raise ValueError(f"{text} cannot be converted to float")
