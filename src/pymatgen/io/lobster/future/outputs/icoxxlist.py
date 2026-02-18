from __future__ import annotations

import re
from collections import defaultdict
from itertools import islice
from typing import TYPE_CHECKING

import numpy as np

from pymatgen.electronic_structure.core import Spin
from pymatgen.io.lobster.future.core import LobsterInteractionsHolder
from pymatgen.io.lobster.future.outputs.coxxcar import COXXCAR
from pymatgen.io.lobster.future.utils import parse_orbital_from_text
from pymatgen.io.lobster.future.versioning import version_processor

if TYPE_CHECKING:
    from typing import ClassVar, Literal

    from numpy.typing import NDArray

    from pymatgen.io.lobster.future.types import LobsterInteractionData


class ICOXXLIST(LobsterInteractionsHolder):
    """Reader for ICOXX data files (ICOHPLIST, ICOOPLIST, ICOBILIST).

    Parses interaction data from ICOXXLIST files, including spin-resolved values.

    Attributes:
        interactions (list[LobsterInteractionData]): List of parsed interactions.
        spins (list[Spin]): List of spins present in the file.
        data (NDArray[np.floating]): Array of ICOXX values for each interaction and spin.
        icoxxlist_type (str): Type of ICOXXLIST file ("COHP", "COOP", "COBI").
        is_lcfo (bool): Whether the file is in LCFO format (if applicable).
    """

    interactions_regex: ClassVar[str] = (
        r"(?i)\s*(\d+)\s+(\S+\s+\S+)\s+(\d+\.\d+)\s+(\-?\d+\s+\-?\d+\s+\-?\d+)?\s+(\-?\d+\.\d+)(?:\s+(\-?\d+\.\d+))?"
    )
    icoxxlist_type: ClassVar[str]

    @version_processor(max_version="5.1")
    def parse_file_legacy(self) -> None:
        """Parse ICOXXLIST file using legacy format (versions ≤5.1).

        Extracts interaction data, including spin-resolved values, and populates
        the `interactions`, `spins`, and `data` attributes.

        Raises:
            ValueError: If the file contains invalid spin values or cannot parse
                interaction lines.
        """
        self.interactions = []
        self.spins = []

        interaction_counter = 0
        for line in self.iterate_lines():
            if not line:
                continue

            if line.startswith(f"{self.icoxxlist_type.upper()}#"):
                if spin_regex := re.search(r"(?i)for spin\s+(\d)", line):
                    spin_regex = int(spin_regex.group(1))
                else:
                    continue

                if spin_regex == 1:
                    self.spins.append(Spin.up)
                elif spin_regex == 2:
                    self.spins.append(Spin.down)
                    interaction_counter = 0
                else:
                    raise ValueError(f"Invalid spin value {spin_regex} in line: {line}")
            else:
                if matches := re.search(self.interactions_regex, line):
                    matches = matches.groups()
                else:
                    raise ValueError(f"Could not parse interaction line: {line}")

                first_center, second_center = matches[1].split()

                first_orbital = parse_orbital_from_text(first_center)
                second_orbital = parse_orbital_from_text(second_center)

                index = int(matches[0])
                centers = [first_center, second_center]
                length = float(matches[2])

                cells = (
                    [[], []]
                    if matches[3] is None
                    else [[0, 0, 0], [int(i) for i in matches[3].split()]]
                )

                bond_tmp: LobsterInteractionData = {
                    "index": index,
                    "centers": centers,
                    "cells": cells,
                    "orbitals": [first_orbital, second_orbital],
                    "length": length,
                    "icoxx": {
                        self.spins[-1]: float(matches[4]),
                    },
                }

                if self.spins[-1] == Spin.up:
                    self.interactions.append(bond_tmp)
                elif self.spins[-1] == Spin.down:
                    interaction = self.interactions[interaction_counter]
                    if "icoxx" in interaction:
                        interaction["icoxx"][Spin.down] = float(matches[4])
                    else:
                        raise ValueError(
                            f"Down spin ICOXX value found without corresponding up spin value in line: {line}"
                        )

                interaction_counter += 1

        self.data = np.full(
            (len(self.interactions), len(self.spins)), np.nan, dtype=np.float64
        )

        for i, interaction in enumerate(self.interactions):
            if "icoxx" in interaction:
                icoxx = interaction["icoxx"]
            else:
                raise ValueError(f"No ICOXX data found for interaction: {interaction}")

            if Spin.up in icoxx:
                self.data[i, 0] = icoxx[Spin.up]
            if Spin.down in icoxx:
                self.data[i, 1] = icoxx[Spin.down]

    @version_processor(min_version="5.1")
    def parse_file(self) -> None:
        """Parse ICOXXLIST file using modern format (versions ≥5.1).

        Extracts interaction data, including spin-resolved values, and populates
        the `interactions`, `spins`, and `data` attributes.

        Raises:
            ValueError: If the file contains invalid spin values or cannot parse
                interaction lines.
        """
        self.interactions = []
        self.spins = []

        for line in islice(self.iterate_lines(), 1, None):
            if not line:
                continue

            if spin_regex := re.findall(r"(?i)for spin\s+(\d)", line):
                self.spins.append(Spin.up)

                if len(spin_regex) == 2:
                    self.spins.append(Spin.down)
            else:
                if matches := re.search(self.interactions_regex, line):
                    matches = matches.groups()
                else:
                    raise ValueError(f"Could not parse interaction line: {line}")

                first_center, second_center = matches[1].split()

                first_orbital = parse_orbital_from_text(first_center)
                second_orbital = parse_orbital_from_text(second_center)

                index = int(matches[0])
                centers = [
                    first_center.replace(f"_{first_orbital}", ""),
                    second_center.replace(f"_{second_orbital}", ""),
                ]
                length = float(matches[2])

                cells = (
                    [[], []]
                    if matches[3] is None
                    else [[0, 0, 0], [int(i) for i in matches[3].split()]]
                )

                bond_tmp: LobsterInteractionData = {
                    "index": index,
                    "centers": centers,
                    "cells": cells,
                    "orbitals": [first_orbital, second_orbital],
                    "length": length,
                    "icoxx": {
                        Spin.up: float(matches[4]),
                    },
                }

                if len(self.spins) == 2:
                    bond_tmp["icoxx"][Spin.down] = float(matches[5])

                self.interactions.append(bond_tmp)

        self.data = np.full(
            (len(self.interactions), len(self.spins)), np.nan, dtype=np.float64
        )

        for i, interaction in enumerate(self.interactions):
            if "icoxx" in interaction:
                icoxx = interaction["icoxx"]
            else:
                raise ValueError(f"No ICOXX data found for interaction: {interaction}")

            if Spin.up in icoxx:
                self.data[i, 0] = icoxx[Spin.up]
            if Spin.down in icoxx:
                self.data[i, 1] = icoxx[Spin.down]

    def get_data_by_properties(
        self: LobsterInteractionsHolder,
        indices: list[int] | None = None,
        centers: list[str] | None = None,
        cells: list[list[int]] | None = None,
        orbitals: list[str] | None = None,
        length: tuple[float, float] | None = None,
        spins: list[Literal[1, -1]] | None = None,
    ) -> NDArray[np.floating]:
        """Get the data for bonds matching specified properties.

        Args:
            indices (list[int] | None): Indices of bonds to retrieve.
            centers (list[str] | None): Centers of bonds to retrieve.
            cells (list[list[int]] | None): Cells of bonds to retrieve.
            orbitals (list[str] | None): Orbitals of bonds to retrieve.
            length (tuple[float, float] | None): Length range to filter.
            spins (list[Spin] | None): Spins to retrieve.

        Returns:
            NDArray[np.floating]: Array of data for specified bonds.
        """
        interaction_indices = self.get_interaction_indices_by_properties(
            indices, centers, cells, orbitals, length
        )

        spins = spins or self.spins
        spin_indices = [0 if spin == Spin.up else 1 for spin in spins]

        return self.data[np.ix_(interaction_indices, spin_indices)]

    def process_data_into_interactions(self) -> None:
        """Populate each interaction dict with 'coxx' and 'icoxx' views.

        Assigns numpy views into `self.data` for each spin channel.
        """
        spin_indices = {spin: i for i, spin in enumerate(self.spins)}

        for i, interaction in enumerate(self.interactions):
            interaction["icoxx"] = {}
            for spin, index in spin_indices.items():
                interaction["icoxx"][spin] = float(self.data[i, index])


class ICOHPLIST(ICOXXLIST):
    """Reader for ICOHPLIST.lobster files.

    Attributes:
        icoxxlist_type (str): Type of ICOXXLIST file ("COHP").
        is_lcfo (bool): Whether the file is in LCFO format.
    """

    icoxxlist_type: ClassVar[str] = "COHP"
    is_lcfo: ClassVar[bool] = False

    @classmethod
    def get_default_filename(cls) -> str:
        """Return the default filename for `ICOHPLIST`."""
        return "ICOHPLIST.LCFO.lobster" if cls.is_lcfo else "ICOHPLIST.lobster"


class ICOOPLIST(ICOXXLIST):
    """Reader for ICOOPLIST.lobster files.

    Attributes:
        icoxxlist_type (str): Type of ICOXXLIST file ("COOP").
    """

    icoxxlist_type: ClassVar[str] = "COOP"

    @classmethod
    def get_default_filename(cls) -> str:
        """Return the default filename for `ICOOPLIST`."""
        return "ICOOPLIST.lobster"


class ICOBILIST(ICOXXLIST):
    """Reader for ICOBILIST.lobster files.

    Attributes:
        icoxxlist_type (str): Type of ICOXXLIST file ("COBI").
        is_lcfo (bool): Whether the file is in LCFO format.
    """

    icoxxlist_type: ClassVar[str] = "COBI"
    is_lcfo: ClassVar[bool] = False

    @classmethod
    def get_default_filename(cls) -> str:
        """Return the default filename for ICOBILIST."""
        return "ICOBILIST.LCFO.lobster" if cls.is_lcfo else "ICOBILIST.lobster"


class ICOHPLIST_LCFO(ICOHPLIST):
    """Reader for ICOHPLIST.LCFO.lobster files.

    Attributes:
        is_lcfo (bool): Always True for LCFO format.
    """

    is_lcfo: ClassVar[bool] = True


class ICOBILIST_LCFO(ICOBILIST):
    """Reader for ICOBILIST.LCFO.lobster files.

    Attributes:
        is_lcfo (bool): Always True for LCFO format.
    """

    is_lcfo: ClassVar[bool] = True


class NcICOBILIST(LobsterInteractionsHolder):
    """Reader for NcICOBILIST.lobster files.

    Parses non-conventional ICOBI interaction data.

    Attributes:
        interactions (list): List of parsed interactions.
        spins (list[Spin]): List of spins present in the file.
        data (NDArray[np.floating]): Array of ICOXX values for each interaction and spin.
    """

    interactions_regex: ClassVar[str] = COXXCAR.interactions_regex

    @version_processor(min_version="5.1")
    def parse_file(self) -> None:
        """Parse the `NcICOBILIST` file.

        Extracts interaction data, including spin-resolved values, and populates
        the `interactions`, `spins`, and `data` attributes.

        Raises:
            ValueError: If the file contains invalid spin values or cannot parse
                interaction lines.
        """
        self.interactions = []
        self.spins = []

        interaction_counter = 0

        for line in self.iterate_lines():
            if not line:
                continue

            if line.startswith(("COBI#", "for spin")):
                interaction_counter = 0
                if match := re.search(r"(?i)for spin\s+(\d)", line):
                    spin_regex = int(match.group(1))
                else:
                    continue

                if spin_regex == 1:
                    self.spins.append(Spin.up)
                elif spin_regex == 2:
                    self.spins.append(Spin.down)
                else:
                    raise ValueError(f"Invalid spin value {spin_regex} in line: {line}")
            else:
                line = re.split(r"\s+(?![^\[]*\])", line)

                bond_tmp = defaultdict(list)
                length = None

                index = int(line[0])

                nc_icobi_value = float(line[-2])

                for center in line[-1].split("->"):
                    if match := re.search(self.interactions_regex, center):
                        match = match.groups()
                    else:
                        raise ValueError(
                            f"Could not parse interaction center line: {line}"
                        )

                    bond_tmp["centers"].append(match[0])

                    if match[1]:
                        bond_tmp["cells"].append([int(x) for x in match[1].split()])
                    else:
                        bond_tmp["cells"].append([])

                    bond_tmp["orbitals"].append(match[2])

                    if match[3]:
                        length = float(match[3])

                current_spin = self.spins[-1]
                if current_spin == Spin.up:
                    self.interactions.append(
                        {
                            "index": index,
                            "centers": bond_tmp["centers"],
                            "cells": bond_tmp["cells"],
                            "orbitals": bond_tmp["orbitals"],
                            "length": length,
                            "icoxx": {
                                Spin.up: nc_icobi_value,
                            },
                        }
                    )
                elif current_spin == Spin.down:
                    interaction = self.interactions[interaction_counter]
                    if "icoxx" in interaction:
                        interaction["icoxx"][Spin.down] = nc_icobi_value
                else:
                    raise ValueError(
                        f"Invalid spin value {current_spin} in line: {line}"
                    )

                interaction_counter += 1

        self.data = np.full(
            (len(self.interactions), len(self.spins)), np.nan, dtype=np.float64
        )

        for i, interaction in enumerate(self.interactions):
            if "icoxx" in interaction:
                icoxx = interaction["icoxx"]
            else:
                raise ValueError(f"No ICOXX data found for interaction: {interaction}")

            if Spin.up in icoxx:
                self.data[i, 0] = icoxx[Spin.up]
            if Spin.down in icoxx:
                self.data[i, 1] = icoxx[Spin.down]

    get_data_by_properties = ICOXXLIST.get_data_by_properties

    process_data_into_interactions = ICOXXLIST.process_data_into_interactions

    @classmethod
    def get_default_filename(cls) -> str:
        """Return the default filename for NcICOBILIST."""
        return "NcICOBILIST.lobster"
