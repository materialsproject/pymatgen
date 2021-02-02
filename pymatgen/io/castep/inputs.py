import collections
import re
from dataclasses import dataclass, field
from typing import Dict, Any, Union, List
from warnings import warn

import numpy as np
from monty.io import zopen
from monty.json import MSONable

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.castep.constants import CellKeyword, ParamKeyword

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

Tag = collections.namedtuple("Tag", "value comment")
Block = collections.namedtuple("Block", "values comments")


@dataclass
class Cell(MSONable):
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

    blocks: Dict[CellKeyword, Block] = field(default_factory=dict)
    tags: Dict[CellKeyword, Tag] = field(default_factory=dict)

    def add_block(self, name: Union[str, CellKeyword], block: Union[List[List[str]], Block]):
        """
        Add a block to the Cell.

        Args:
            name: A valid block name. Can be given as a string (any case) or
                a CellKeyword. If the block name is not known, an exception
                will be raised.
            block: A Block, a named tuple contained values and comments. Can
                be given as a list of strings, in which case the comment will
                be left empty.
        """
        if isinstance(name, str):
            name = CellKeyword[name.upper()]
        if isinstance(block, list):
            block = Block(values=block, comments=None)
        self.blocks[name] = block

    def add_tag(self, name: Union[str, CellKeyword], tag: Union[str, Tag]):
        """
        Add a tag to the Cell.

        Args:
            name: A valid block name. Can be given as a string (any case) or
                a CellKeyword. If the block name is not known, an exception
                will be raised.
            tag: A Tag, a named tuple contained value and comments. Can
                be given as a list of strings, in which case the comment will
                be left empty.
        """
        if isinstance(name, str):
            name = CellKeyword[name.upper()]
        if isinstance(tag, str):
            tag = Block(value=tag, comments=None)
        self.tags[name] = tag

    def __str__(self):
        """
        Returns: The string representation Cell, i.e. the contents of the .cell file
        """

        lines = []

        for tag, content in self.tags.items():
            if content.comment in (None, ""):
                lines.append("{0: <24}: {1}".format(tag.name.lower(), " ".join(content.value)))
            else:
                lines.append(
                    "{0: <24}: {1: <16} ! {2}".format(tag.name.lower(), " ".join(content.value), content.comment)
                )
        for block, content in self.blocks.items():
            lines.append("\n%block {}".format(block.name.lower()))
            if content.comments is None:
                comments = [""] * len(content.values)
            else:
                comments = content.comments

            for row, comment in zip(content.values, comments):
                line = " ".join(map(str, row))
                if comment != "":
                    line = "{0: <30} ! {1}".format(line, comment)
                lines.append(line)

            lines.append("%endblock {}".format(block.name.lower()))

        return "\n".join(lines)

    def write_file(self, filename):
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
        if len(comment_split_line) == 1:
            comment = ""
        else:
            comment = line[len(line_content) :]

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
    def from_file(cls, filename):
        """
        Load a .cell file.

        Args:
            filename: filename of .cell file

        Returns: a Cell
        """

        with zopen(filename, "rt") as f:
            lines = [line.strip() for line in f]

        # Remove lines which are entirely commented or empty
        def _is_not_empty_line(line):
            if len(line) == 0:
                return False
            elif line[0] in "#;!":
                return False
            elif len(line) > 6 and line[:7] == "COMMENT":
                return False
            else:
                return True

        lines = list(filter(_is_not_empty_line, lines))

        tags, blocks = {}, {}
        current_block_values, current_block_comments = [], []
        in_block = False
        current_block_label = None

        for line in lines:
            if len(line.split()) == 0:
                continue
            elif line[:6].lower() == "%block":
                if in_block:
                    raise IOError("Cell file contains nested blocks. This possibility was not anticipated.")
                else:
                    current_block_label = CellKeyword[line.split()[1].upper()]
                    in_block = True
                    continue
            elif line[:9].lower() == "%endblock":
                if not in_block:
                    raise IOError("Cannot cope with line {}: not currently in a block.".format(line))
                if CellKeyword[line.split()[1].upper()] != current_block_label:
                    raise IOError(
                        "Endblock {} does not match current block "
                        "{}: cannot interpret cell file.".format(line.split()[1], current_block_label)
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

    def add_structure(self, structure: Structure):
        """
        Add a structure to your .cell file. If a structure is already specified,
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

        self.add_block(CellKeyword.LATTICE_CART, [["ang"]] + [list(map(str, row)) for row in structure.lattice.matrix])
        self.add_block(
            CellKeyword.POSITIONS_FRAC,
            [
                [site.species_string, *map(str, site.frac_coords), "spin={}".format(site.properties["magmom"])]
                if "magmom" in site.properties
                else [site.species_string, *map(str, site.frac_coords)]
                for site in structure.sites
            ],
        )
        self._structure = self.structure

    @classmethod
    def from_structure(cls, structure, blocks: Dict[CellKeyword, Block] = None, tags: Dict[CellKeyword, Tag] = None):
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
        cell.add_structure(structure)

        return cell

    @property
    def structure(self) -> Structure:
        """
        Returns: the Structure object defined by this .cell file
        """

        if hasattr(self, "_structure") and getattr(self, "_structure"):
            return self._structure

        # Get lattice vectors
        if CellKeyword.LATTICE_CART in self.blocks:
            lattice_cart = self.blocks[CellKeyword.LATTICE_CART].values
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
            lattice_abc = self.blocks[CellKeyword.LATTICE_ABC].values
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
                    magmoms[idx] = float(row.split("spin=")[1])
            elements, coords = zip(*elements_coords)
            structure = Structure(lattice, elements, coords, coords_are_cartesian=False)
        elif CellKeyword.POSITIONS_ABS in self.blocks:
            # TODO: handle magmom!
            positions_abs = self.blocks[CellKeyword.POSITIONS_ABS].values
            if positions_abs[0][0].lower() in ("ang", "nm", "cm", "m", "bohr", "a0"):
                unit = positions_abs[0][0].lower()
                positions_abs = positions_abs[1:]
            else:
                unit = "ang"
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
            structure.add_site_property(
                "magmom", [magmoms[idx] if idx in magmoms else None for idx in range(len(structure))]
            )

        # store to avoid re-build if this property is called again
        self._structure = structure

        return structure


@dataclass
class Param(MSONable):
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

    tags: Dict[ParamKeyword, Tag] = field(default_factory=dict)