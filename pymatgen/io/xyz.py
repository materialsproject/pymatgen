"""Module implementing an XYZ file object class."""

from __future__ import annotations

import re
from io import StringIO
from typing import TYPE_CHECKING, cast

import pandas as pd
from monty.io import zopen

from pymatgen.core import Molecule, Structure
from pymatgen.core.structure import SiteCollection

if TYPE_CHECKING:
    from collections.abc import Sequence
    from pathlib import Path

    from typing_extensions import Self


class XYZ:
    """
    Basic class for importing and exporting Molecules or Structures in XYZ
    format.

    Note:
        Exporting periodic structures in the XYZ format will lose information
        about the periodicity. Essentially, only Cartesian coordinates are
        written in this format and no information is retained about the
        lattice.
    """

    def __init__(self, mol: Molecule | Structure | Sequence[Molecule | Structure], coord_precision: int = 6) -> None:
        """
        Args:
            mol (Molecule | Structure): Input molecule or structure or list thereof.
            coord_precision: Precision to be used for coordinates.
        """
        self._mols = cast(list[SiteCollection], [mol] if isinstance(mol, SiteCollection) else mol)
        self.precision = coord_precision

    @property
    def molecule(self) -> Molecule:
        """
        Returns molecule associated with this XYZ. In case of multi-frame
        XYZ, returns the last frame.
        """
        return self._mols[-1]  # type: ignore[return-value]

    @property
    def all_molecules(self) -> list[Molecule]:
        """Returns all the frames of molecule associated with this XYZ."""
        return self._mols  # type: ignore[return-value]

    @staticmethod
    def _from_frame_str(contents) -> Molecule:
        """Convert a single frame XYZ string to a molecule."""
        lines = contents.split("\n")
        n_sites = int(lines[0])
        coords = []
        sp = []
        coord_pattern = re.compile(r"(\w+)\s+([0-9\-\+\.*^eEdD]+)\s+([0-9\-\+\.*^eEdD]+)\s+([0-9\-\+\.*^eEdD]+)")
        for idx in range(2, 2 + n_sites):
            if match := coord_pattern.search(lines[idx]):
                sp.append(match.group(1))  # this is 1-indexed
                # this is 0-indexed
                # in case of 0.0D+00 or 0.00d+01 old double precision writing
                # replace d or D by e for ten power exponent,
                # and some files use *^ convention in place of e
                xyz = [val.lower().replace("d", "e").replace("*^", "e") for val in match.groups()[1:4]]
                coords.append([float(val) for val in xyz])
        return Molecule(sp, coords)

    @classmethod
    def from_str(cls, contents: str) -> Self:
        """
        Creates XYZ object from a string.

        Args:
            contents: String representing an XYZ file.

        Returns:
            XYZ object
        """
        if contents[-1] != "\n":
            contents += "\n"
        white_space = r"[ \t\r\f\v]"
        n_atoms_line = white_space + r"*\d+" + white_space + r"*\n"
        comment_line = r"[^\n]*\n"
        coord_lines = r"(\s*\w+\s+[0-9\-\+\.*^eEdD]+\s+[0-9\-\+\.*^eEdD]+\s+[0-9\-\+\.*^eEdD]+.*\n)+"
        frame_pattern_text = n_atoms_line + comment_line + coord_lines
        pat = re.compile(frame_pattern_text, re.MULTILINE)
        mols = []
        for xyz_match in pat.finditer(contents):
            xyz_text = xyz_match.group(0)
            mols.append(XYZ._from_frame_str(xyz_text))
        return cls(mols)

    @classmethod
    def from_file(cls, filename: str | Path) -> Self:
        """
        Creates XYZ object from a file.

        Args:
            filename: XYZ filename

        Returns:
            XYZ object
        """
        with zopen(filename, mode="rt") as file:
            return cls.from_str(file.read())

    def as_dataframe(self):
        """
        Generates a coordinates data frame with columns: atom, x, y, and z
        In case of multiple frame XYZ, returns the last frame.

        Returns:
            pandas.DataFrame
        """
        lines = str(self)
        sio = StringIO(lines)
        df_xyz = pd.read_csv(
            sio, header=None, skiprows=(0, 1), comment="#", delim_whitespace=True, names=("atom", "x", "y", "z")
        )
        df_xyz.index += 1
        return df_xyz

    def _frame_str(self, frame_mol):
        output = [str(len(frame_mol)), frame_mol.formula]
        prec = self.precision
        fmt = f"{{}} {{:.{prec}f}} {{:.{prec}f}} {{:.{prec}f}}"
        for site in frame_mol:
            output.append(fmt.format(site.specie, site.x, site.y, site.z))
        return "\n".join(output)

    def __str__(self):
        return "\n".join(self._frame_str(mol) for mol in self._mols)

    def write_file(self, filename: str) -> None:
        """
        Writes XYZ to file.

        Args:
            filename (str): File name of output file.
        """
        with zopen(filename, mode="wt") as file:
            file.write(str(self))
