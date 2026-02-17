from __future__ import annotations

from abc import abstractmethod
from collections import Counter
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, cast

import numpy as np
from monty.io import zopen
from monty.json import MontyDecoder, MSONable
from typing_extensions import Self

from pymatgen.io.lobster.future.constants import LOBSTER_VERSION
from pymatgen.io.lobster.future.utils import convert_spin_keys, restore_spin_keys

if TYPE_CHECKING:
    from collections.abc import Callable, Generator

    from pymatgen.io.lobster.future.types import LobsterInteraction, LobsterInteractionData, Spin
    from pymatgen.util.typing import PathLike


class LobsterFile(MSONable):
    """
    Representation of a LOBSTER file.

    This class provides a framework for parsing and processing LOBSTER output files.
    It supports version-specific processing through a registry of version processors.

    Attributes:
        filename (PathLike): Name or path of the file. Defaults to the class default.
        lobster_version (str): Version string parsed from the file or the default LOBSTER_VERSION.
        version_processors (ClassVar[dict[tuple[str, str | None], Callable]]): Registry of version processors.
        spins (list[Spin] | None): List of Spin objects if the file contains spin-polarized data, else None.
    """

    version_processors: ClassVar[dict[tuple[str, str | None], Callable]]

    spins: list[Spin] | None = None

    def __init__(self, filename: PathLike | None = None, process_immediately: bool = True):
        """
        Initialize a LobsterFile instance.

        Args:
            filename (PathLike | None): Path to the file. If None, uses the default filename.
            process_immediately (bool): Whether to process the file immediately upon initialization. Defaults to True.
        """
        self.filename = Path(filename or self.get_default_filename()).expanduser().resolve()
        self.lobster_version = self.get_file_version() or LOBSTER_VERSION

        if process_immediately:
            self.process()

    def __init_subclass__(cls, **kwargs: Any) -> None:
        """
        Automatically registers version processors for subclasses.
        This method scans the subclass for methods decorated with @version_processor
        and adds them to the version_processors registry.
        """
        super().__init_subclass__(**kwargs)

        version_processors = {}
        for base in cls.__bases__:
            if hasattr(base, "version_processors"):
                version_processors.update(base.version_processors)

        for value in cls.__dict__.values():
            if hasattr(value, "version_info"):
                version_processors[value.version_info] = value

        cls.version_processors = version_processors

    def process(self) -> None:
        """
        Processes the file using the appropriate version processor.

        Selects the best matching version processor from the `version_processors` registry and invokes it.

        Raises:
            ValueError: If no processor matches the file version.
            RuntimeError: If the selected processor raises an exception.
        """
        eligible_methods = []

        if not self.filename.exists():
            raise FileNotFoundError(f"The file {self.filename} does not exist.")

        for (
            minimum_version,
            maximum_version,
        ), processor in self.version_processors.items():
            if LobsterFile.check_version(self.lobster_version, minimum_version, maximum_version):
                eligible_methods.append((minimum_version, processor))

        if not eligible_methods:
            raise ValueError(f"No processor found for version {self.lobster_version}")

        best_method = max(eligible_methods, key=lambda x: x[0])[-1]

        try:
            best_method(self)
        except Exception as e:
            processor_name = getattr(best_method, "__name__", str(best_method))

            raise RuntimeError(
                f"Error occurred during file processing with {processor_name} (version {self.lobster_version}): {e}"
            ) from e

    @staticmethod
    def check_version(actual: str, minimum: str, maximum: str | None) -> bool:
        """
        Checks whether a version string falls within a min/max inclusive range.

        Args:
            actual (str): Version string to check (e.g., "1.2.3").
            minimum (str): Minimum acceptable version string (exclusive).
            maximum (str | None): Maximum acceptable version string (exclusive) or None for no upper bound.

        Returns:
            bool: True if `actual` is > `minimum` and < `maximum` (if provided), otherwise False.
        """
        actual_parts = tuple(map(int, actual.split(".")))
        minimum_parts = tuple(map(int, minimum.split(".")))

        if actual_parts < minimum_parts:
            return False

        if maximum is not None:
            maximum_parts = tuple(map(int, maximum.split(".")))

            if actual_parts > maximum_parts:
                return False

        return True

    @abstractmethod
    def get_file_version(self) -> str | None:
        """
        Retrieves the file version. Override in subclasses to extract version from file content if possible.

        Returns:
            str | None: Version string (e.g., "1.2.3") if found, else None.
        """
        ...

    @classmethod
    @abstractmethod
    def get_default_filename(cls) -> str:
        """
        Returns the default filename for this LobsterFile subclass.

        Returns:
            str: The default filename.
        """
        ...

    @cached_property
    def lines(self) -> list[str]:
        """
        Returns all lines from the file as a list of strings.

        Returns:
            list[str]: Lines from the file.

        Raises:
            ValueError: If the file is empty.
        """
        with zopen(self.filename, mode="rt", encoding="utf-8") as file:
            lines = file.read().splitlines()

            if len(lines) == 0:
                raise ValueError(f"The file {self.filename} is empty.")

            return cast("list[str]", lines)

    def iterate_lines(self) -> Generator[str, None, None]:
        """
        Iterates over lines in the file, yielding each line as a string.

        Yields:
            str: Each line in the file, stripped of whitespace.
        """
        with zopen(self.filename, mode="rt", encoding="utf-8") as file:
            for line in file:
                yield cast("str", line.strip())

    def as_dict(self) -> dict[str, Any]:
        """
        Serializes the LobsterFile object to a dictionary.
        Spin keys in dictionaries are converted to strings for JSON compatibility.

        Returns:
            dict[str, Any]: Dictionary with keys "@module", "@class", "@version", and all attributes of the object.
        """
        dictionary = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "@version": None,
        }

        for k, v in vars(self).items():
            dictionary[k] = convert_spin_keys(v)

        return dictionary

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> Self:
        """
        Deserializes a LobsterFile object from a dictionary. Spin keys in dictionaries are restored from strings.

        Args:
            d (dict[str, Any]): Dictionary produced by as_dict or similar serialization.

        Returns:
            LobsterFile: Instance of LobsterFile with attributes populated from the dictionary.
        """
        instance = cls.__new__(cls)

        decoded_dictionary = {
            k: restore_spin_keys(MontyDecoder().process_decoded(v)) for k, v in d.items() if not k.startswith("@")
        }

        for k, v in decoded_dictionary.items():
            setattr(instance, k, v)

        instance.filename = Path(instance.filename)

        return instance

    @property
    def has_spin(self) -> bool:
        """
        Indicates whether the file could contain spin-polarized data.

        Returns:
            bool: True if this file type supports spin, False otherwise.
        """
        return self.spins is not None

    @property
    def is_spin_polarized(self) -> bool:
        """
        Indicates whether the file contains spin-polarized data.

        Returns:
            bool: True if multiple spins are present, False otherwise.
        """
        return self.has_spin and len(self.spins) > 1


class LobsterInteractionsHolder:
    """
    Container for LOBSTER interaction data. This class holds interaction metadata. It provides methods for filtering and
    retrieving interaction data based on various criteria.

    Attributes:
        interactions (list[LobsterInteractionData]): List of interaction metadata dicts.
    """

    def __init__(self):
        """
        Initializes the LobsterInteractionsHolder object. Attributes are expected to be set by the caller or by
        deserialization.
        """
        self.interactions: list[LobsterInteractionData]

    def get_interaction_indices_by_properties(
        self,
        indices: list[int] | None = None,
        centers: list[str] | None = None,
        cells: list[list[int]] | None = None,
        orbitals: list[str] | None = None,
        length: tuple[float, float] | None = None,
    ) -> list[int]:
        """
        Returns indices of interactions that match provided property filters.

        Args:
            indices (list[int] | None): Optional sequence of interaction indices to match.
            centers (list[str] | None): Optional sequence of center name substrings; interaction must contain each
            center substring the required number of times.
            cells (list[list[int]] | None): Optional sequence of cell vectors to match against interaction cells.
            orbitals (list[str] | None): Optional sequence of orbital name substrings; interaction must contain each
            orbital substring the required number of times.
            length (tuple[float, float] | None): Optional (min, max) tuple to filter interactions by length inclusive.

        Returns:
            list[int]: Sorted list of interaction indices that match all supplied filters. If no filters are supplied,
            returns an empty list.
        """
        matching_sets = []

        if indices is not None:
            matching_sets.append({i for i, bond in enumerate(self.interactions) if bond["index"] in indices})

        if centers is not None:
            center_counts = Counter(centers)
            matching_centers = set()

            for i, bond in enumerate(self.interactions):
                bond_centers = bond.get("centers", [])

                if all(
                    sum(center_suffix in b for b in bond_centers) >= required_count
                    for center_suffix, required_count in center_counts.items()
                ):
                    matching_centers.add(i)

            matching_sets.append(matching_centers)

        if cells is not None:
            matching_sets.append(
                {
                    i
                    for i, bond in enumerate(self.interactions)
                    if bond.get("cells")
                    and all(np.all(np.equal(bond.get("cells"), cell), axis=1).any() for cell in cells if cell)
                }
            )

        if orbitals is not None:

            if not orbitals:
                matching_sets.append(
                    {i for i, bond in enumerate(self.interactions) 
                    if all(o is None for o in bond.get("orbitals", []))}
                )
            else:
                orbital_counts = Counter(orbitals)
                matching_orbitals = set()

                for i, bond in enumerate(self.interactions):
                    bond_orbitals = bond.get("orbitals", [])

                    if all(
                        sum(orbital_suffix in b for b in bond_orbitals if b) >= required_count
                        for orbital_suffix, required_count in orbital_counts.items()
                    ):
                        matching_orbitals.add(i)

                matching_sets.append(matching_orbitals)

        if length is not None:
            matching_sets.append(
                {
                    i
                    for i, bond in enumerate(self.interactions)
                    if (this_length := bond.get("length")) is not None and length[0] <= this_length <= length[1]
                }
            )

        return sorted(set.intersection(*matching_sets)) if matching_sets else []

    @staticmethod
    def get_label_from_interaction(
        interaction: LobsterInteraction,
        include_centers: bool = True,
        include_orbitals: bool = True,
        include_cells: bool = False,
        include_length: bool = False,
    ) -> str:
        """
        Generates a label string for a given interaction.

        Args:
            interaction (LobsterInteraction): Interaction metadata dictionary.
        Returns:
            str: Formatted label string representing the interaction.
        """
        parts = []

        for center, orbital, cell in zip(
            interaction["centers"],
            interaction["orbitals"],
            interaction["cells"],
            strict=True,
        ):
            tmp = ""
            if include_centers:
                tmp += center

            if include_cells and cell is not None:
                tmp += f"[{' '.join(map(str, cell))}]"

            if include_orbitals and orbital:
                tmp += f"[{orbital}]"

            parts.append(tmp)

        if not parts:
            raise ValueError(f"Could not generate label from interaction {interaction}")

        if include_length and interaction["length"] is not None:
            parts[-1] += f"({interaction['length']:.3f})"

        return "->".join(parts)
