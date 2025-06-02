"""A representation of FHI-aims output (based on ASE output parser)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from monty.json import MontyDecoder, MSONable

from pymatgen.io.aims.parsers import (
    read_aims_header_info,
    read_aims_header_info_from_content,
    read_aims_output,
    read_aims_output_from_content,
)

if TYPE_CHECKING:
    from collections.abc import Sequence
    from pathlib import Path
    from typing import Any

    from typing_extensions import Self

    from pymatgen.core import Molecule, Structure

__author__ = "Andrey Sobolev and Thomas A. R. Purcell"
__version__ = "1.0"
__email__ = "andrey.n.sobolev@gmail.com and purcellt@arizona.edu"
__date__ = "November 2023"


class AimsOutput(MSONable):
    """The main output file for FHI-aims."""

    def __init__(
        self,
        results: Molecule | Structure | Sequence[Molecule | Structure],
        metadata: dict[str, Any],
        structure_summary: dict[str, Any],
    ) -> None:
        """
        Args:
            results (Molecule or Structure or Sequence[Molecule or Structure]):  A list
                of all images in an output file
            metadata (Dict[str, Any]): The metadata of the executable used to perform
                the calculation
            structure_summary (Dict[str, Any]): The summary of the starting
                atomic structure.
        """
        self._results = results
        self._metadata = metadata
        self._structure_summary = structure_summary

    def as_dict(self) -> dict[str, Any]:
        """Create a dict representation of the outputs for MSONable."""
        dct: dict[str, Any] = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
        }

        dct["results"] = self._results
        dct["metadata"] = self._metadata
        dct["structure_summary"] = self._structure_summary
        return dct

    @classmethod
    def from_outfile(cls, outfile: str | Path) -> Self:
        """Construct an AimsOutput from an output file.

        Args:
            outfile: str | Path: The aims.out file to parse

        Returns:
            The AimsOutput object for the output file
        """
        metadata, structure_summary = read_aims_header_info(outfile)
        results = read_aims_output(outfile, index=slice(0, None))

        return cls(results, metadata, structure_summary)

    @classmethod
    def from_str(cls, content: str) -> Self:
        """Construct an AimsOutput from an output file.

        Args:
            content (str): The content of the aims.out file

        Returns:
            The AimsOutput for the output file content
        """
        metadata, structure_summary = read_aims_header_info_from_content(content)
        results = read_aims_output_from_content(content, index=slice(0, None))

        return cls(results, metadata, structure_summary)

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        """Construct an AimsOutput from a dictionary.

        Args:
            dct (dict[str, Any]): The dictionary used to create AimsOutput

        Returns:
            AimsOutput
        """
        decoded = {k: MontyDecoder().process_decoded(v) for k, v in dct.items() if not k.startswith("@")}
        for struct in decoded["results"]:
            struct.properties = {k: MontyDecoder().process_decoded(v) for k, v in struct.properties.items()}

        return cls(
            decoded["results"],
            decoded["metadata"],
            decoded["structure_summary"],
        )

    def get_results_for_image(self, image_ind: int) -> Structure | Molecule:
        """Get the results dictionary for a particular image or slice of images.

        Args:
            image_ind (int): The index of the image to get the results for

        Returns:
            The results of the image with index images_ind
        """
        return self._results[image_ind]

    @property
    def structure_summary(self) -> dict[str, Any]:
        """The summary of the material/molecule that the calculations represent."""
        return self._structure_summary

    @property
    def metadata(self) -> dict[str, Any]:
        """The system metadata."""
        return self._metadata

    @property
    def n_images(self) -> int:
        """The number of images in results."""
        return len(self._results)

    @property
    def initial_structure(self) -> Structure | Molecule:
        """The initial structure for the calculations."""
        return self._structure_summary["initial_structure"]

    @property
    def final_structure(self) -> Structure | Molecule:
        """The final structure for the calculation."""
        return self._results[-1]

    @property
    def structures(self) -> Sequence[Structure | Molecule]:
        """All images in the output file."""
        return self._results

    @property
    def fermi_energy(self) -> float:
        """The Fermi energy for the final structure in the calculation."""
        return self.get_results_for_image(-1).properties["fermi_energy"]

    @property
    def vbm(self) -> float:
        """The HOMO level for the final structure in the calculation."""
        return self.get_results_for_image(-1).properties["vbm"]

    @property
    def cbm(self) -> float:
        """The LUMO level for the final structure in the calculation."""
        return self.get_results_for_image(-1).properties["cbm"]

    @property
    def band_gap(self) -> float:
        """The band gap for the final structure in the calculation."""
        return self.get_results_for_image(-1).properties["gap"]

    @property
    def direct_band_gap(self) -> float:
        """The direct band gap for the final structure in the calculation."""
        return self.get_results_for_image(-1).properties["direct_gap"]

    @property
    def final_energy(self) -> float:
        """The total energy for the final structure in the calculation."""
        return self.get_results_for_image(-1).properties["energy"]

    @property
    def completed(self) -> bool:
        """Did the calculation complete."""
        return len(self._results) > 0

    @property
    def aims_version(self) -> str:
        """The version of FHI-aims used for the calculation."""
        return self._metadata["version_number"]

    @property
    def forces(self) -> Sequence[tuple[float, float, float]] | None:
        """The forces for the final image of the calculation."""
        force_array = self.get_results_for_image(-1).site_properties.get("force", None)
        if isinstance(force_array, np.ndarray):
            return force_array.tolist()

        return force_array

    @property
    def stress(
        self,
    ) -> np.ndarray:
        """The stress for the final image of the calculation."""
        return self.get_results_for_image(-1).properties.get("stress", None)

    @property
    def stresses(
        self,
    ) -> list[list[float]] | None:
        """The atomic virial stresses for the final image of the calculation."""
        stresses_array = self.get_results_for_image(-1).site_properties.get("atomic_virial_stress", None)
        if isinstance(stresses_array, np.ndarray):
            return stresses_array.tolist()
        return stresses_array  # type:ignore[return-value]

    @property
    def all_forces(self) -> list[list[tuple[float, float, float]]]:
        """The forces for all images in the calculation."""
        all_forces_array = [res.site_properties.get("force", None) for res in self._results]
        return [af.tolist() if isinstance(af, np.ndarray) else af for af in all_forces_array]
