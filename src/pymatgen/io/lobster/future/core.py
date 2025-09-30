from __future__ import annotations

from abc import abstractmethod
from collections import Counter
from functools import cached_property
from typing import TYPE_CHECKING, Any, ClassVar, TypeVar, cast

import numpy as np
from monty.io import zopen
from monty.json import MontyDecoder, MSONable

from pymatgen.io.lobster.future.constants import LOBSTER_VERSION

if TYPE_CHECKING:
    from collections.abc import Callable, Generator
    from typing import Literal

    from numpy.typing import NDArray

    from pymatgen.electronic_structure.core import Spin
    from pymatgen.io.lobster.future.types import LobsterInteractionData
    from pymatgen.util.typing import PathLike

Self = TypeVar("Self", bound="LobsterFile")


class LobsterFile(MSONable):
    """
    Representation of a LOBSTER file.

    This class provides a framework for parsing and processing LOBSTER output files.
    It supports version-specific processing through a registry of version processors.

    Attributes:
        filename (PathLike): Name or path of the file. Defaults to the class default.
        lobster_version (str): Version string parsed from the file or the default LOBSTER_VERSION.
        version_processors (ClassVar[dict[tuple[str, str | None], Callable]]): Registry of version processors.
    """

    version_processors: ClassVar[dict[tuple[str, str | None], Callable]]

    def __init__(
        self, filename: PathLike | None = None, process_immediately: bool = True
    ):
        """
        Initialize a LobsterFile instance.

        Args:
            filename (PathLike | None): Path to the file. If None, uses the default filename.
            process_immediately (bool): Whether to process the file immediately upon initialization. Defaults to True.
        """
        self.filename = filename or self.get_default_filename()
        self.lobster_version = self.get_file_version() or LOBSTER_VERSION

        if process_immediately:
            self.process()

    def __init_subclass__(cls, **kwargs):
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

        for (
            minimum_version,
            maximum_version,
        ), processor in self.version_processors.items():
            if LobsterFile.check_version(
                self.lobster_version, minimum_version, maximum_version
            ):
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
            minimum (str): Minimum acceptable version string (inclusive).
            maximum (str | None): Maximum acceptable version string (inclusive) or None for no upper bound.

        Returns:
            bool: True if `actual` is >= `minimum` and <= `maximum` (if provided), otherwise False.
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
        Retrieves the file version.

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
        Serializes the LobsterFile to a JSON-serializable dictionary.

        Returns:
            dict[str, Any]: Dictionary with keys "@module", "@class", "@version", "kwargs", and "attributes".
        """
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "@version": None,
            "kwargs": {"filename": str(self.filename)},
            "attributes": {
                "lobster_version": self.lobster_version,
            },
        }

    @classmethod
    def from_dict(cls: type[Self], d: dict[str, Any]) -> Self:  # NOQA: PYI019
        """
        Deserializes a LobsterFile from a dictionary.

        Args:
            d (dict[str, Any]): Dictionary produced by as_dict or similar serialization.

        Returns:
            LobsterFile: Instance of LobsterFile with attributes populated from the dictionary.
        """
        decoded_dictionary = {
            k: MontyDecoder().process_decoded(v)
            for k, v in d.items()
            if not k.startswith("@")
        }
        kwargs: dict[str, Any] = decoded_dictionary.get("kwargs", {})
        kwargs.setdefault("process_immediately", False)

        instance = cls(**kwargs)

        attributes: dict[str, Any] = decoded_dictionary.get("attributes", {})
        for attr, value in attributes.items():
            setattr(instance, attr, value)

        return instance


class LobsterDataHolder:
    """
    Container for LOBSTER interaction data and associated arrays.

    This class holds interaction metadata, numeric data arrays, and spin information
    for LOBSTER output files. It provides methods for filtering and retrieving interaction
    data based on various criteria.

    Attributes:
        interactions (list[LobsterInteractionData]): List of interaction metadata dicts.
        data (NDArray[np.floating]): Array containing numeric data for interactions.
        spins (list[Spin]): List of Spin objects corresponding to data spin channels.
    """

    def __init__(self):
        """
        Initializes the LobsterDataHolder.

        Attributes are expected to be set by the caller or by deserialization.
        """
        self.interactions: list[LobsterInteractionData]
        self.data: NDArray[np.floating]
        self.spins: list[Spin]

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
            matching_sets.append(
                {
                    i
                    for i, bond in enumerate(self.interactions)
                    if bond["index"] in indices
                }
            )

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
                    and all(
                        np.all(np.equal(bond.get("cells"), cell), axis=1).any()
                        for cell in cells
                        if cell
                    )
                }
            )

        if orbitals is not None:
            orbital_counts = Counter(orbitals)
            matching_orbitals = set()

            for i, bond in enumerate(self.interactions):
                bond_orbitals = bond.get("orbitals", [])

                if all(
                    sum(orbital_suffix in b for b in bond_orbitals if b)
                    >= required_count
                    for orbital_suffix, required_count in orbital_counts.items()
                ):
                    matching_orbitals.add(i)

            matching_sets.append(matching_orbitals)

        if length is not None:
            matching_sets.append(
                {
                    i
                    for i, bond in enumerate(self.interactions)
                    if (this_length := bond.get("length")) is not None
                    and length[0] <= this_length <= length[1]
                }
            )

        return sorted(set.intersection(*matching_sets)) if matching_sets else []

    def get_header_from_interaction_indices(
        self,
        interaction_indices: int | list[int],
        spins: list[Spin] | None = None,
        data_type: Literal["coxx", "icoxx"] | None = None,
    ) -> list[str]:
        """
        Builds header strings for specified interaction indices.

        Args:
            interaction_indices (int | list[int]): Single index or sequence of interaction indices to include.
            spins (list[Spin] | None): Optional sequence of Spin objects to include in headers; if None,
            all spins are used.
            data_type (Literal["coxx", "icoxx"] | None): Optional data type string or None to include both "coxx" and
            "icoxx".

        Returns:
            list[str]: Header strings constructed for each requested interaction, spin, and data type.
        """
        if isinstance(interaction_indices, int):
            interaction_indices = [interaction_indices]

        if spins is None:
            spins = self.spins

        header = []
        for s in spins:
            for i, bond in enumerate(self.interactions):
                if i in interaction_indices:
                    for d in [data_type] if data_type else ["coxx", "icoxx"]:
                        header.append(f"No.{bond['index']}:")

                        length = f"({bond['length']})" if bond["length"] else ""

                        tmp_string = []

                        for center, cell, orbital in zip(
                            bond["centers"],
                            bond["cells"],
                            bond["orbitals"],
                            strict=True,
                        ):
                            orbital = f"[{orbital}]" if orbital else ""
                            cell = "[{:d} {:d} {:d}]".format(*cell) if cell else ""

                            tmp_string.append(f"{center}{cell}{orbital}")

                        header[-1] += f"{'->'.join(tmp_string)}[{s.name}][{d}]{length}"

        return header
