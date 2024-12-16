from __future__ import annotations

import re
from dataclasses import dataclass, field
from difflib import get_close_matches
from enum import EnumMeta
from typing import TYPE_CHECKING, Any, NamedTuple
from warnings import warn

import numpy as np
from monty.io import zopen

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.castep.constants import CellKeyword, OTFPseudopotentialLibrary, ParamKeyword
from pymatgen.io.core import InputFile

if TYPE_CHECKING:
    from enum import Enum
    from os import PathLike

"""
This module provides an interface between pymatgen and the CASTEP (http://www.castep.org)
electronic structure code for the purposes of generating input files or parsing output files, pymatgen
does not offer any capability to run CASTEP directly. To run or use CASTEP, an appropriate license has
to be obtained from the CASTEP developers.

"First principles methods using CASTEP", Zeitschrift fuer Kristallographie 220(5-6) pp. 567-570 (2005)
S. J. Clark, M. D. Segall, C. J. Pickard, P. J. Hasnip, M. J. Probert, K. Refson, M. C. Payne
"""

# Thank you to Adam Jackson for some of the initial implementation of this module
# which has been borrowed from sumo (https://github.com/SMTG-UCL/sumo) with permission.

_bohr_to_angstrom = 0.5291772
_ry_to_ev = 13.605693009

# From 2002 CODATA values
to_angstrom = {"ang": 1, "bohr": 0.5291772108, "m": 1e10, "cm": 1e8, "nm": 10}


class Tag(NamedTuple):
    value: Any
    comment: str | None


class Block(NamedTuple):
    values: list[Any]
    comments: list[str] | None


def _get_enum(key: str | Enum, enum: EnumMeta) -> EnumMeta:
    """
    Utility function to create a specified enum from a string.
    """
    if isinstance(key, EnumMeta):
        return key
    if isinstance(key, str):
        try:
            return enum[key.upper()]
        except KeyError:
            error_message = f"{key} is not a valid keyword."
            possible_matches = get_close_matches(key.upper(), enum.__members__.keys(), n=3)
            if possible_matches:
                error_message += f" Perhaps you mean: {', '.join(possible_matches)}?"
            raise KeyError(error_message)
    raise ValueError(f"{key} is provided incorrectly, it is not an {enum} or a string and cannot be parsed.")


@dataclass
class Cell(InputFile):
    """
    An interface for CASTEP's .cell file which defines atomic
    positions, k-point grids etc. The CASTEP documentation is
    the canonical resource for what the .cell file can contain.

    A .cell file contains "tags" (single-line values) and "blocks"
    (multi-line values). Valid names for tags and blocks are
    given by the CellKeyword Enum in pymatgen.io.castep.constants

    Tags are stored as a dict, where each key is a tag and each value is a
    Tag NamedTuple with the attributes 'value' and 'comment'.

    Blocks are also stored as a dict, where each key is the block name
    (e.g. 'bs_kpoint_list') and the value is a Block NamedTuple with attributes
    'values' (a list of lists) and 'comments' (a commensurate list of str).
    """

    blocks: dict[CellKeyword, Block] = field(default_factory=dict)
    tags: dict[CellKeyword, Tag] = field(default_factory=dict)

    def set_block(self, name: str | CellKeyword, block: list[list[str]] | Block):
        """
        Set a block in the Cell.

        Args:
            name: A valid block name. Can be given as a string (any case) or
                a CellKeyword. If the block name is not known, an exception
                will be raised.
            block: A Block, a named tuple contained values and comments. Can
                be given as a list of strings, in which case the comment will
                be left empty.
        """
        name = _get_enum(name, CellKeyword)
        if isinstance(block, list):
            block = Block(values=block, comments=None)
        self.blocks[name] = block

    def set_tag(self, name: str | CellKeyword, tag: str | Tag):
        """
        Set a tag in the Cell.

        Args:
            name: A valid block name. Can be given as a string (any case) or
                a CellKeyword. If the block name is not known, an exception
                will be raised.
            tag: A Tag, a named tuple contained value and comments. Can
                be given as a list of strings, in which case the comment will
                be left empty.
        """
        name = _get_enum(name, CellKeyword)
        if isinstance(tag, str):
            tag = Tag(value=tag, comment=None)
        self.tags[name] = tag

    def set_pseudopotential(self, library: str | OTFPseudopotentialLibrary):
        """
        Convenience method to add one of the default pseudopotential libraries for all elements.

        Args:
            library: one of the standard CASTEP pseudpotential names, as specified by
            the OTFPseudopotentialLibrary Enum, e.g. "C19"
        """
        if isinstance(library, str):
            library = OTFPseudopotentialLibrary[library]
        self.set_block(CellKeyword.SPECIES_POT, [[library.name]])

    def get_str(self):
        """
        Returns: The string representation Cell, i.e. the contents of the .cell file
        """

        lines = []

        for tag, content in self.tags.items():
            if content.comment in (None, ""):
                lines.append("{: <24}: {}".format(tag.name.lower(), " ".join(content.value)))
            else:
                lines.append("{: <24}: {: <16} ! {}".format(tag.name.lower(), " ".join(content.value), content.comment))
        for block, content in self.blocks.items():
            lines.append(f"\n%block {block.name.lower()}")
            comments = [""] * len(content.values) if content.comments is None else content.comments

            for row, comment in zip(content.values, comments):
                line = " ".join(map(str, row))
                if comment != "":
                    line = f"{line: <30} ! {comment}"
                lines.append(line)

            lines.append(f"%endblock {block.name.lower()}")

        return "\n".join(lines)

    def write_file(self, filename: str | PathLike) -> None:
        """
        Write the .cell file

        Args:
            filename: filename for the .cell file
        """
        with zopen(filename, "wt") as f:
            f.write(str(self))

    @staticmethod
    def _data_comment_from_line(line: str, in_block: bool = False):
        """Parse a line from a CASTEP input file
        Args:
            line (str): line from CASTEP .cell or .param file
            in_block (bool): Include all tokens before comment in
                data and leave tag as None.
        Returns:
            (tag, data, comment) where tag is a str, data is a list
            of str and comment is a str
        """
        comment_split_line = re.split("[#;!]", line)

        line_content = comment_split_line[0]
        comment = "" if len(comment_split_line) == 1 else line[len(line_content) :]

        if in_block:
            tag, data = None, line.split()
        else:
            tag, *data = line_content.split()
            if len(data) == 0:
                data = [""]

            # Clean up delimiter: 'TAG DAT' 'TAG= DAT' 'TAG : DAT'
            # are all possible so check both symbols in both places
            data = [item for item in data if item not in ":="]
            if tag[-1] in ":=":
                tag = tag[:-1]

        return tag, data, comment

    @classmethod
    def from_str(cls, contents: str) -> Cell:
        """
        Load from the string representation of a .cell file.

        Args:
            contents: .cell file string

        Returns: a Cell
        """

        lines = [line.strip() for line in contents.splitlines()]

        # Remove lines which are entirely commented or empty
        def _is_not_empty_line(line):
            if len(line) == 0 or line[0] in "#;!":
                return False
            return not (len(line) > 6 and line[:7] == "COMMENT")

        lines = list(filter(_is_not_empty_line, lines))

        tags, blocks = {}, {}
        current_block_values, current_block_comments = [], []
        in_block = False
        current_block_label = None

        for line in lines:
            if len(line.split()) == 0:
                continue
            if line[:6].lower() == "%block":
                if in_block:
                    raise OSError("Cell file contains nested blocks. This possibility was not anticipated.")
                current_block_label = CellKeyword[line.split()[1].upper()]
                in_block = True
                continue
            if line[:9].lower() == "%endblock":
                if not in_block:
                    raise OSError(f"Cannot cope with line {line}: not currently in a block.")
                if CellKeyword[line.split()[1].upper()] != current_block_label:
                    raise OSError(
                        f"Endblock {line.split()[1]} does not match current block "
                        f"{current_block_label}: cannot interpret cell file."
                    )
                blocks[current_block_label] = Block(current_block_values, current_block_comments)
                current_block_values, current_block_comments = [], []
                in_block = False

            elif in_block:
                _, data, comment = cls._data_comment_from_line(line, in_block=True)
                current_block_values.append(data)
                current_block_comments.append(comment)

            else:
                tag, data, comment = cls._data_comment_from_line(line)
                tags[CellKeyword[tag.upper()]] = Tag(data, comment)

        return cls(tags=tags, blocks=blocks)

    def set_structure(self, structure: Structure):
        """
        Set a structure to your .cell file. If a structure is already specified,
        it will be replaced.

        Args:
            structure: Structure object defining your crystal structure
        """

        for key in (
            CellKeyword.LATTICE_ABC,
            CellKeyword.LATTICE_CART,
            CellKeyword.POSITIONS_FRAC,
            CellKeyword.POSITIONS_ABS,
        ):
            if key in self.blocks:
                del self.blocks[key]
                warn("Structure already specified in .cell file, replacing with new structure.")

        self.set_block(CellKeyword.LATTICE_CART, [["ang"]] + [list(map(str, row)) for row in structure.lattice.matrix])
        self.set_block(
            CellKeyword.POSITIONS_FRAC,
            [
                [site.species_string, *map(str, site.frac_coords), "spin={}".format(site.properties["magmom"])]
                if "magmom" in site.properties
                else [site.species_string, *map(str, site.frac_coords)]
                for site in structure.sites
            ],
        )

        # clear existing cached structure
        self._structure = None
        # and cache the new one
        self._structure = self.structure

    @classmethod
    def from_structure(
        cls, structure, blocks: dict[CellKeyword, Block] | None = None, tags: dict[CellKeyword, Tag] | None = None
    ):
        """
        Initialize a .cell file with a provided input structure.

        Args:
            structure: pymatgen Structure
            blocks: any additional blocks, keyed by block keyword
            tags: any additional tags, keyed by tag keyword

        Returns: a Cell object
        """

        tags = tags or {}
        blocks = blocks or {}

        cell = cls(blocks=blocks, tags=tags)
        cell.set_structure(structure)

        return cell

    @property
    def structure(self) -> Structure:
        """
        Returns: the Structure object defined by this .cell file
        """

        if hasattr(self, "_structure") and self._structure:
            return self._structure

        # Get lattice vectors
        if CellKeyword.LATTICE_CART in self.blocks:
            lattice_cart = self.blocks[CellKeyword.LATTICE_CART].values  # noqa: PD011
            if lattice_cart[0][0].lower() in ("ang", "nm", "cm", "m", "bohr", "a0"):
                unit = lattice_cart[0][0].lower()
                vectors = lattice_cart[1:]
            else:
                unit = "ang"
                vectors = lattice_cart
            vectors = np.array([list(map(float, row)) for row in vectors])
            if vectors.shape != (3, 3):
                raise ValueError("lattice_cart should contain a 3x3 matrix")

            vectors *= to_angstrom[unit]
            lattice = Lattice(vectors)
        elif CellKeyword.LATTICE_ABC in self.blocks:
            lattice_abc = self.blocks[CellKeyword.LATTICE_ABC].values  # noqa: PD011
            if lattice_abc[0][0].lower() in ("ang", "nm", "cm", "m", "bohr", "a0"):
                unit = lattice_abc[0][0].lower()
                lengths_and_angles = lattice_abc[1:]
            else:
                unit = "ang"
                lengths_and_angles = lattice_abc[1:]
            if len(lengths_and_angles) != 2:
                raise ValueError("lattice_abc should have two rows")
            lengths_and_angles = [list(map(float, row)) for row in lengths_and_angles]
            lengths_and_angles[0] = [x * to_angstrom[unit] for x in lengths_and_angles[0]]
            lattice = Lattice.from_lengths_and_angles(*lengths_and_angles)
        else:
            raise ValueError("Couldn't find a lattice in cell file!")

        if CellKeyword.POSITIONS_FRAC in self.blocks:
            elements_coords = []
            magmoms = {}
            for idx, row in enumerate(self.blocks[CellKeyword.POSITIONS_FRAC].values):
                elements_coords.append((row[0], list(map(float, row[1:4]))))
                if "spin=" in row:
                    # assumes e.g. "spin=+2" with no space around equals sign
                    # TODO: should just handle any numeral after SPIN (strip letters)
                    magmoms[idx] = float(row.split("spin=")[1])
            elements, coords = zip(*elements_coords)
            structure = Structure(lattice, elements, coords, coords_are_cartesian=False)
        elif CellKeyword.POSITIONS_ABS in self.blocks:
            positions_abs = self.blocks[CellKeyword.POSITIONS_ABS].values  # noqa: PD011
            if positions_abs[0][0].lower() in ("ang", "nm", "cm", "m", "bohr", "a0"):
                unit = positions_abs[0][0].lower()
                positions_abs = positions_abs[1:]
            else:
                unit = "ang"
                # TODO: make sure unit variable is used
            elements_coords = []
            magmoms = {}
            for idx, row in enumerate(positions_abs):
                elements_coords.append(row[0], list(map(float, row[1:4])))
                if "spin=" in row:
                    # assumes e.g. "spin=+2" with no space around equals sign
                    magmoms[idx] = float(row.split("spin=")[1])
            elements, coords = zip(*elements_coords)
            structure = Structure(lattice, elements, coords, coords_are_cartesian=True)
        else:
            raise ValueError("Couldn't find any atomic positions in cell file!")

        if magmoms:
            structure.add_site_property("magmom", [magmoms.get(idx, None) for idx in range(len(structure))])

        # store to avoid re-build if this property is called again
        self._structure = structure

        return structure


@dataclass
class Param(InputFile):
    """
    An interface for CASTEP's .param file which defines
    calculation parameters. The CASTEP documentation is
    the canonical resource for what the .param file can contain.

    A .param file contains "tags" (single-line values).
    Valid names for tags and blocks are given by the ParamKeyword
    Enum in pymatgen.io.castep.constants

    Tags are stored as a dict, where each key is a tag and each value is a
    Tag NamedTuple with the attributes 'value' and 'comment'.
    """

    tags: dict[ParamKeyword, Tag] = field(default_factory=dict)

    def set_tag(self, name: str | ParamKeyword, tag: str | Tag):
        """
        Set a tag in the Param.

        Args:
            name: A valid block name. Can be given as a string (any case) or
                a ParamKeyword. If the block name is not known, an exception
                will be raised.
            tag: A Tag, a named tuple contained value and comments. Can
                be given as a list of strings, in which case the comment will
                be left empty.
        """
        name = _get_enum(name, ParamKeyword)
        if isinstance(tag, str):
            tag = Tag(value=tag, comment=None)
        self.tags[name] = tag

    @classmethod
    def from_tags(cls, tags: dict[str | ParamKeyword, Any]):
        """
        Convenience method to create a Param file using a dictionary
        of tags and their associated values.

        Args:
            tags: dictionary of parameters, keyed by parameter name (case insensitive)
                or keyed by ParamKeyword Enum,
                e.g. {"task": "singlepoint", "basis_precision": "fine"}

        Returns: a Param object
        """
        param = cls()
        for name, value in tags.items():
            param.set_tag(name, value)
        return param

    def get_str(self) -> str:
        """
        Returns: The string representation Param, i.e. the contents of the .param file
        """

        lines = []
        for keyword, tag in self.tags.items():
            if tag.comment:
                lines.append(f"{keyword:20} = {tag.value} ! {tag.comment}")
            else:
                lines.append(f"{keyword:20} = {tag.value}")

        return "\n".join(lines)

    @classmethod
    def from_str(cls, contents: str) -> Param:
        tags = {}
        lines = contents.splitlines()

        for line in lines:
            line = line.strip()
            if not line.startswith("!"):  # comment line
                keyword, line = line.split("=")
                value = line.split("!")
                if len(value) == 2:
                    value, comment = value
                    comment = comment.strip()
                else:
                    value = value[0]
                    comment = None
                keyword = ParamKeyword(keyword.strip().upper())
                tags[keyword] = Tag(value=value.strip(), comment=comment)

        return cls(tags)
