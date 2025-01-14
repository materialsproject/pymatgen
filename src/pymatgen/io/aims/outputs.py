"""A representation of FHI-aims output (based on ASE output parser)."""

from __future__ import annotations

import gzip
import warnings
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from monty.json import MontyDecoder, MSONable
from pyfhiaims.outputs.stdout import AimsParseError, AimsStdout

from pymatgen.core import Lattice, Structure

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

    from typing_extensions import Self

    from pymatgen.core import Molecule
    from pymatgen.util.typing import Matrix3D, Vector3D

__author__ = "Andrey Sobolev and Thomas A. R. Purcell"
__version__ = "1.0"
__email__ = "andrey.n.sobolev@gmail.com and purcellt@arizona.edu"
__date__ = "November 2023"


AIMS_OUTPUT_KEY_MAP = {
    "free_energy": "energy",  # TARP These are the force consistent energies
}


def remap_outputs(results: dict[str, Any]) -> dict[str, Any]:
    """Remap FHIAimsOutput keys to AimsOutput keys"""
    to_ret = results.copy()
    for key, val in AIMS_OUTPUT_KEY_MAP.items():
        to_ret[val] = to_ret.pop(key, None)

    return to_ret


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

        Raises:
            AimsParseError if file does not exist
        """
        aims_out = None
        for path in [Path(outfile), Path(f"{outfile}.gz")]:
            if not path.exists():
                continue
            if path.suffix == ".gz":
                with gzip.open(f"{outfile}.gz", mode="rt") as file:
                    aims_out = AimsStdout(file)
            else:
                with open(outfile) as file:
                    aims_out = AimsStdout(file)

        if aims_out is None:
            raise AimsParseError(f"File {outfile} not found.")

        metadata = aims_out.metadata_summary
        structure_summary = aims_out.header_summary

        structure_summary["initial_structure"] = structure_summary.pop("initial_geometry").structure
        for site in structure_summary["initial_structure"]:
            if abs(site.properties.get("magmom", 0.0)) < 1e-10:
                site.properties.pop("magmom", None)

        lattice = structure_summary.pop("initial_lattice", None)
        if lattice is not None:
            lattice = Lattice(lattice)
        structure_summary["initial_lattice"] = lattice

        results = []
        for image in aims_out:
            image_results = remap_outputs(image._results)
            structure = image.geometry.structure
            site_prop_keys = {
                "forces": "force",
                "stresses": "atomic_virial_stress",
                "hirshfeld_charges": "hirshfeld_charge",
                "hirshfeld_volumes": "hirshfeld_volume",
                "hirshfeld_atomic_dipoles": "hirshfeld_atomic_dipole",
                "mulliken_charges": "charge",
                "mulliken_spins": "magmom",
            }
            properties = {prop: image_results[prop] for prop in image_results if prop not in site_prop_keys}
            site_properties = {}
            for prop, site_key in site_prop_keys.items():
                if prop in image_results:
                    site_properties[site_key] = image_results[prop]

            if ((magmom := site_properties.get("magmom")) is not None) and np.abs(
                np.sum(magmom) - properties["magmom"]
            ) < 1e-3:
                warnings.warn(
                    "Total magnetic moment and sum of Mulliken spins are not consistent",
                    stacklevel=2,
                )
            if isinstance(structure, Structure):
                structure._properties.update(properties)
            else:
                structure.properties.update(properties)
            for st, site in enumerate(structure.sites):
                site.properties = {key: val[st] for key, val in site_properties.items()}

            results.append(structure)

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
    def forces(self) -> Sequence[Vector3D] | None:
        """The forces for the final image of the calculation."""
        force_array = self.get_results_for_image(-1).site_properties.get("force", None)
        if isinstance(force_array, np.ndarray):
            return force_array.tolist()

        return force_array

    @property
    def stress(self) -> Matrix3D:
        """The stress for the final image of the calculation."""
        return self.get_results_for_image(-1).properties.get("stress", None)

    @property
    def stresses(self) -> Sequence[Matrix3D] | None:
        """The atomic virial stresses for the final image of the calculation."""
        stresses_array = self.get_results_for_image(-1).site_properties.get("atomic_virial_stress", None)
        if isinstance(stresses_array, np.ndarray):
            return stresses_array.tolist()
        return stresses_array

    @property
    def all_forces(self) -> list[list[Vector3D]]:
        """The forces for all images in the calculation."""
        all_forces_array = [res.site_properties.get("force", None) for res in self._results]
        return [af.tolist() if isinstance(af, np.ndarray) else af for af in all_forces_array]
