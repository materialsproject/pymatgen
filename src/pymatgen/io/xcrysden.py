"""Support for reading XCrysDen files."""

from __future__ import annotations

from typing import TYPE_CHECKING

from pymatgen.core import Element, IStructure, Structure

if TYPE_CHECKING:
    from typing_extensions import Self

    from pymatgen.core.structure import IStructure

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"


class XSF:
    """Parse XCrysden files."""

    def __init__(self, structure: Structure | IStructure):
        """
        Args:
            structure (Structure): Structure object.
        """
        self.structure = structure

    def to_str(self, atom_symbol: bool = True) -> str:
        """Get a string with the structure in XSF format
        See http://www.xcrysden.org/doc/XSF.html.

        Args:
            atom_symbol (bool): Uses atom symbol instead of atomic number. Defaults to True.
        """
        lines: list[str] = []

        lines.extend(("CRYSTAL", "# Primitive lattice vectors in Angstrom", "PRIMVEC"))
        cell = self.structure.lattice.matrix
        for i in range(3):
            lines.append(f" {cell[i][0]:.14f} {cell[i][1]:.14f} {cell[i][2]:.14f}")

        cart_coords = self.structure.cart_coords
        lines.extend(
            (
                "# Cartesian coordinates in Angstrom.",
                "PRIMCOORD",
                f" {len(cart_coords)} 1",
            )
        )

        for site, coord in zip(self.structure, cart_coords, strict=True):
            sp = site.specie.symbol if atom_symbol else f"{site.specie.Z}"
            x, y, z = coord
            lines.append(f"{sp} {x:20.14f} {y:20.14f} {z:20.14f}")
            if "vect" in site.properties:
                vx, vy, vz = site.properties["vect"]
                lines[-1] += f" {vx:20.14f} {vy:20.14f} {vz:20.14f}"

        return "\n".join(lines)

    @classmethod
    def from_str(cls, input_string: str, cls_=None) -> Self:
        """
        Initialize a `Structure` object from a string with data in XSF format.

        Args:
            input_string: String with the structure in XSF format.
                See http://www.xcrysden.org/doc/XSF.html
            cls_: Structure class to be created. default: pymatgen structure

        Example file:
            CRYSTAL                                        see (1)
            these are primitive lattice vectors (in Angstroms)
            PRIMVEC
            0.0000000    2.7100000    2.7100000         see (2)
            2.7100000    0.0000000    2.7100000
            2.7100000    2.7100000    0.0000000

            these are conventional lattice vectors (in Angstroms)
            CONVVEC
            5.4200000    0.0000000    0.0000000         see (3)
            0.0000000    5.4200000    0.0000000
            0.0000000    0.0000000    5.4200000

            these are atomic coordinates in a primitive unit cell  (in Angstroms)
            PRIMCOORD
            2 1                                            see (4)
            16      0.0000000     0.0000000     0.0000000  see (5)
            30      1.3550000    -1.3550000    -1.3550000
        """
        lattice, coords, species = [], [], []
        lines = input_string.splitlines()

        for idx, line in enumerate(lines, start=1):
            if "PRIMVEC" in line:
                for j in range(idx, idx + 3):
                    lattice.append([float(c) for c in lines[j].split()])

            if "PRIMCOORD" in line:
                n_sites = int(lines[idx].split()[0])

                for j in range(idx + 1, idx + 1 + n_sites):
                    tokens = lines[j].split()
                    Z = Element(tokens[0]).Z if tokens[0].isalpha() else int(tokens[0])
                    species.append(Z)
                    coords.append([float(j) for j in tokens[1:4]])
                break
        else:
            raise ValueError("Invalid XSF data")

        if cls_ is None:
            cls_ = Structure

        return cls(cls_(lattice, species, coords, coords_are_cartesian=True))
