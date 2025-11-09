"""This module provides input and output from the CSSR file format."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from monty.io import zopen

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure

if TYPE_CHECKING:
    from pathlib import Path

    from typing_extensions import Self

    from pymatgen.core.structure import IStructure

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jan 24, 2012"


class Cssr:
    """
    Basic object for working with Cssr file. Right now, only conversion from
    a Structure to a Cssr file is supported.
    """

    def __init__(self, structure: Structure | IStructure):
        """
        Args:
            structure (Structure | IStructure): A structure to create the Cssr object.
        """
        if not structure.is_ordered:
            raise ValueError("Cssr file can only be constructed from ordered structure")
        self.structure = structure

    def __str__(self):
        a, b, c = self.structure.lattice.abc
        alpha, beta, gamma = self.structure.lattice.angles
        output = [
            f"{a:.4f} {b:.4f} {c:.4f}",
            f"{alpha:.2f} {beta:.2f} {gamma:.2f} SPGR =  1 P 1    OPT = 1",
            f"{len(self.structure)} 0",
            f"0 {self.structure.formula}",
        ]
        for idx, site in enumerate(self.structure, start=1):
            output.append(f"{idx} {site.specie} {site.a:.4f} {site.b:.4f} {site.c:.4f}")
        return "\n".join(output)

    def write_file(self, filename):
        """Write out a CSSR file.

        Args:
            filename (str): Filename to write to.
        """
        with zopen(filename, mode="wt", encoding="utf-8") as file:
            file.write(str(self) + "\n")

    @classmethod
    def from_str(cls, string: str) -> Self:
        """
        Reads a string representation to a Cssr object.

        Args:
            string (str): A string representation of a CSSR.

        Returns:
            Cssr object.
        """
        lines = string.split("\n")
        tokens = lines[0].split()
        lengths = [float(tok) for tok in tokens]
        tokens = lines[1].split()
        angles = [float(tok) for tok in tokens[:3]]
        lattice = Lattice.from_parameters(*lengths, *angles)
        sp, coords = [], []
        for line in lines[4:]:
            if match := re.match(
                r"\d+\s+(\w+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)",
                line.strip(),
            ):
                sp.append(match[1])
                coords.append([float(match[i]) for i in range(2, 5)])
        return cls(Structure(lattice, sp, coords))

    @classmethod
    def from_file(cls, filename: str | Path) -> Self:
        """
        Reads a CSSR file to a Cssr object.

        Args:
            filename (str): Filename to read from.

        Returns:
            Cssr object.
        """
        with zopen(filename, mode="rt", encoding="utf-8") as file:
            return cls.from_str(file.read())  # type:ignore[arg-type]
