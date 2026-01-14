"""This module implements symmetry-related structure forms."""

from __future__ import annotations

from typing import TYPE_CHECKING, cast

import numpy as np
from monty.json import MontyEncoder
from tabulate import tabulate

from pymatgen.core.structure import FileFormats, PeriodicSite, Structure

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

    from numpy.typing import ArrayLike
    from typing_extensions import Self

    from pymatgen.core import IStructure
    from pymatgen.symmetry.analyzer import SpacegroupOperations
    from pymatgen.util.typing import PathLike


class SymmetrizedStructure(Structure):
    """This class represents a symmetrized structure, i.e. a structure
    where the space group and symmetry operations are defined. This class is
    typically not called but instead is typically obtained by calling
    pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_symmetrized_structure.

    Attributes:
        equivalent_indices (list[list[int]]): Indices of the sites in the structure that are
            considered equivalent based on the symmetry operations of the space group.
    """

    def __init__(
        self,
        structure: Structure | IStructure,
        spacegroup: SpacegroupOperations,
        equivalent_positions: ArrayLike,
        wyckoff_letters: Sequence[str],
    ) -> None:
        """
        Args:
            structure (Structure): Original structure
            spacegroup (SpacegroupOperations): An input SpacegroupOperations from SpacegroupAnalyzer.
            equivalent_positions (list[int]): Equivalent positions from SpacegroupAnalyzer.
            wyckoff_letters (list[str]): Wyckoff letters.
        """
        self.spacegroup = spacegroup
        uniq, inverse = np.unique(equivalent_positions, return_inverse=True)
        inverse = np.atleast_1d(inverse)  # Ensures `inverse` is 1D
        self.site_labels = equivalent_positions

        super().__init__(
            structure.lattice,
            [site.species for site in structure],
            structure.frac_coords,
            site_properties=structure.site_properties,
            properties=structure.properties,
            labels=structure.labels,
        )

        equivalent_indices: list[list[int]] = [[] for _ in range(len(uniq))]
        equivalent_sites: list[list[PeriodicSite]] = [[] for _ in range(len(uniq))]
        wyckoff_symbols: list[list[str]] = [[] for _ in range(len(uniq))]
        for idx, inv_ in enumerate(inverse):
            equivalent_indices[inv_].append(idx)
            equivalent_sites[inv_].append(self.sites[idx])
            wyckoff_symbols[inv_].append(wyckoff_letters[idx])
        self.equivalent_indices = equivalent_indices
        self.equivalent_sites = equivalent_sites
        self.wyckoff_letters = wyckoff_letters
        self.wyckoff_symbols = [f"{len(symb)}{symb[0]}" for symb in wyckoff_symbols]

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        outs = [
            "SymmetrizedStructure",
            f"Full Formula ({self.formula})",
            f"Reduced Formula: {self.reduced_formula}",
            f"Spacegroup: {self.spacegroup.int_symbol} ({self.spacegroup.int_number})",
            f"abc   : {' '.join(f'{val:>10.6f}' for val in self.lattice.abc)}",
            f"angles: {' '.join(f'{val:>10.6f}' for val in self.lattice.angles)}",
        ]

        if self._charge:
            outs.append(f"Overall Charge: {self._charge:+}")
        outs.append(f"Sites ({len(self)})")
        data = []
        props = self.site_properties
        keys = sorted(props)
        for idx, sites in enumerate(self.equivalent_sites):
            site = sites[0]
            row = [str(idx), site.species_string]
            row.extend([f"{j:>10.6f}" for j in site.frac_coords])
            row.append(self.wyckoff_symbols[idx])
            row += [props[key][idx] for key in keys]
            data.append(row)
        outs.append(tabulate(data, headers=["#", "SP", "a", "b", "c", "Wyckoff", *keys]))
        return "\n".join(outs)

    def copy(self) -> Self:
        """Make a copy of the SymmetrizedStructure."""
        return type(self)(
            self,
            spacegroup=self.spacegroup,
            equivalent_positions=self.site_labels,
            wyckoff_letters=self.wyckoff_letters,
        )

    def find_equivalent_sites(self, site: PeriodicSite) -> list[PeriodicSite]:
        """Find all symmetrically equivalent sites for a particular site.

        Args:
            site (PeriodicSite): A site in the structure.

        Raises:
            ValueError: if site is not in the structure.

        Returns:
            list[PeriodicSite]: all symmetrically equivalent sites.
        """
        for sites in self.equivalent_sites:
            if site in sites:
                return sites

        raise ValueError("Site not in structure")

    def as_dict(self) -> dict[str, Any]:
        """MSONable dict."""
        structure = Structure.from_sites(self.sites)
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "structure": structure.as_dict(),
            "spacegroup": self.spacegroup,
            "equivalent_positions": self.site_labels,
            "wyckoff_letters": self.wyckoff_letters,
        }

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            SymmetrizedStructure
        """
        return cls(
            Structure.from_dict(dct["structure"]),
            spacegroup=dct["spacegroup"],
            equivalent_positions=dct["equivalent_positions"],
            wyckoff_letters=dct["wyckoff_letters"],
        )

    def to(self, filename: PathLike = "", fmt: FileFormats = "", **kwargs) -> str:
        """Use `MontyEncoder` as default JSON encoder."""
        filename, fmt = str(filename), cast("FileFormats", fmt.lower())
        if fmt == "json":
            kwargs.setdefault("cls", MontyEncoder)

        return super().to(filename, fmt, **kwargs)
