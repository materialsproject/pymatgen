"""AIMS output parser, taken from ASE with modifications."""

from __future__ import annotations

import gzip
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, cast

import numpy as np

from pymatgen.core import Lattice, Molecule, Structure
from pymatgen.core.tensors import Tensor
from pymatgen.util.typing import Tuple3Floats

if TYPE_CHECKING:
    from collections.abc import Generator, Sequence
    from io import TextIOWrapper
    from typing import Any

    from pymatgen.util.typing import Matrix3D, Vector3D

__author__ = "Thomas A. R. Purcell and Andrey Sobolev"
__version__ = "1.0"
__email__ = "purcellt@arizona.edu and andrey.n.sobolev@gmail.com"
__date__ = "November 2023"

# TARP: Originally an object, but type hinting needs this to be an int
LINE_NOT_FOUND = -1000
EV_PER_A3_TO_KBAR = 1.60217653e-19 * 1e22


class ParseError(Exception):
    """Parse error during reading of a file."""


class AimsParseError(Exception):
    """Exception raised if an error occurs when parsing an Aims output file."""

    def __init__(self, message: str) -> None:
        """Initialize the error with the message, message."""
        self.message = message
        super().__init__(self.message)


# Read aims.out files
SCALAR_PROPERTY_TO_LINE_KEY = {
    "free_energy": ["| Electronic free energy"],
    "number_of_iterations": ["| Number of self-consistency cycles"],
    "magnetic_moment": ["N_up - N_down"],
    "n_atoms": ["| Number of atoms"],
    "n_bands": [
        "Number of Kohn-Sham states",
        "Reducing total number of  Kohn-Sham states",
        "Reducing total number of Kohn-Sham states",
    ],
    "n_electrons": ["The structure contains"],
    "n_kpts": ["| Number of k-points"],
    "n_spins": ["| Number of spin channels"],
    "electronic_temp": ["Occupation type:"],
    "fermi_energy": ["| Chemical potential (Fermi level)"],
}


@dataclass
class AimsOutChunk:
    """Base class for AimsOutChunks.

    Attributes:
        lines (list[str]): The list of all lines in the chunk
    """

    lines: list[str] = field(default_factory=list)

    def reverse_search_for(self, keys: list[str], line_start: int = 0) -> int:
        """Find the last time one of the keys appears in self.lines.

        Args:
            keys (list[str]): The key strings to search for in self.lines
            line_start (int): The lowest index to search for in self.lines

        Returns:
            The last time one of the keys appears in self.lines
        """
        for idx, line in enumerate(self.lines[line_start:][::-1]):
            if any(key in line for key in keys):
                return len(self.lines) - idx - 1

        return LINE_NOT_FOUND

    def search_for_all(self, key: str, line_start: int = 0, line_end: int = -1) -> list[int]:
        """Find the all times the key appears in self.lines.

        Args:
            key (str): The key string to search for in self.lines
            line_start (int): The first line to start the search from
            line_end (int): The last line to end the search at

        Returns:
            All times the key appears in the lines
        """
        line_index = []
        for ll, line in enumerate(self.lines[line_start:line_end]):
            if key in line:
                line_index.append(ll + line_start)
        return line_index

    def parse_scalar(self, property: str) -> float | None:  # noqa: A002
        """Parse a scalar property from the chunk.

        Args:
            property (str): The property key to parse

        Returns:
            The scalar value of the property or None if not found
        """
        line_start = self.reverse_search_for(SCALAR_PROPERTY_TO_LINE_KEY[property])

        if line_start == LINE_NOT_FOUND:
            return None

        line = self.lines[line_start]
        return float(line.split(":")[-1].strip().split()[0])


@dataclass
class AimsOutHeaderChunk(AimsOutChunk):
    """The header of the aims.out file containing general information."""

    lines: list[str] = field(default_factory=list)
    _cache: dict[str, Any] = field(default_factory=dict)

    @property
    def commit_hash(self) -> str:
        """The commit hash for the FHI-aims version."""
        line_start = self.reverse_search_for(["Commit number"])
        if line_start == LINE_NOT_FOUND:
            raise AimsParseError("This file does not appear to be an aims-output file")

        return self.lines[line_start].split(":")[1].strip()

    @property
    def aims_uuid(self) -> str:
        """The aims-uuid for the calculation."""
        line_start = self.reverse_search_for(["aims_uuid"])
        if line_start == LINE_NOT_FOUND:
            raise AimsParseError("This file does not appear to be an aims-output file")

        return self.lines[line_start].split(":")[1].strip()

    @property
    def version_number(self) -> str:
        """The commit hash for the FHI-aims version."""
        line_start = self.reverse_search_for(["FHI-aims version"])
        if line_start == LINE_NOT_FOUND:
            raise AimsParseError("This file does not appear to be an aims-output file")

        return self.lines[line_start].split(":")[1].strip()

    @property
    def fortran_compiler(self) -> str | None:
        """The fortran compiler used to make FHI-aims."""
        line_start = self.reverse_search_for(["Fortran compiler      :"])
        if line_start == LINE_NOT_FOUND:
            raise AimsParseError("This file does not appear to be an aims-output file")

        return self.lines[line_start].split(":")[1].split("/")[-1].strip()

    @property
    def c_compiler(self) -> str | None:
        """The C compiler used to make FHI-aims."""
        line_start = self.reverse_search_for(["C compiler            :"])
        if line_start == LINE_NOT_FOUND:
            return None

        return self.lines[line_start].split(":")[1].split("/")[-1].strip()

    @property
    def fortran_compiler_flags(self) -> str | None:
        """The fortran compiler flags used to make FHI-aims."""
        line_start = self.reverse_search_for(["Fortran compiler flags"])
        if line_start == LINE_NOT_FOUND:
            raise AimsParseError("This file does not appear to be an aims-output file")

        return self.lines[line_start].split(":")[1].strip()

    @property
    def c_compiler_flags(self) -> str | None:
        """The C compiler flags used to make FHI-aims."""
        line_start = self.reverse_search_for(["C compiler flags"])
        if line_start == LINE_NOT_FOUND:
            return None

        return self.lines[line_start].split(":")[1].strip()

    @property
    def build_type(self) -> list[str]:
        """The optional build flags passed to cmake."""
        line_end = self.reverse_search_for(["Linking against:"])
        line_inds = self.search_for_all("Using", line_end=line_end)

        return [" ".join(self.lines[ind].split()[1:]).strip() for ind in line_inds]

    @property
    def linked_against(self) -> list[str]:
        """All libraries used to link the FHI-aims executable."""
        line_start = self.reverse_search_for(["Linking against:"])
        if line_start == LINE_NOT_FOUND:
            return []

        linked_libs = [self.lines[line_start].split(":")[1].strip()]
        line_start += 1
        while "lib" in self.lines[line_start]:
            linked_libs.append(self.lines[line_start].strip())
            line_start += 1

        return linked_libs

    @property
    def initial_lattice(self) -> Lattice | None:
        """The initial lattice vectors from the aims.out file."""
        line_start = self.reverse_search_for(["| Unit cell:"])
        if line_start == LINE_NOT_FOUND:
            return None

        return Lattice(
            np.array(
                [[float(inp) for inp in line.split()[-3:]] for line in self.lines[line_start + 1 : line_start + 4]]
            )
        )

    @property
    def initial_structure(self) -> Structure | Molecule:
        """The initial structure.

        Using the FHI-aims output file recreate the initial structure for
        the calculation.
        """
        lattice = self.initial_lattice

        line_start = self.reverse_search_for(["Atomic structure:"])
        if line_start == LINE_NOT_FOUND:
            raise AimsParseError("No information about the structure in the chunk.")

        line_start += 2

        coords = np.zeros((self.n_atoms, 3))
        species = [""] * self.n_atoms
        for ll, line in enumerate(self.lines[line_start : line_start + self.n_atoms]):
            inp = line.split()
            coords[ll, :] = [float(pos) for pos in inp[4:7]]
            species[ll] = inp[3]

        site_properties = {"charge": self.initial_charges}
        if self.initial_magnetic_moments is not None:
            site_properties["magmoms"] = self.initial_magnetic_moments

        if lattice:
            return Structure(
                lattice,
                species,
                coords,
                np.sum(self.initial_charges),
                coords_are_cartesian=True,
                site_properties=site_properties,
            )

        return Molecule(
            species,
            coords,
            np.sum(self.initial_charges),
            site_properties=site_properties,
        )

    @property
    def initial_charges(self) -> Sequence[float]:
        """The initial charges for the structure."""
        if "initial_charges" not in self._cache:
            self._parse_initial_charges_and_moments()
        return self._cache["initial_charges"]

    @property
    def initial_magnetic_moments(self) -> Sequence[float]:
        """The initial magnetic Moments."""
        if "initial_magnetic_moments" not in self._cache:
            self._parse_initial_charges_and_moments()
        return self._cache["initial_magnetic_moments"]

    def _parse_initial_charges_and_moments(self) -> None:
        """Parse the initial charges and magnetic moments from a file."""
        charges = np.zeros(self.n_atoms)
        magmoms = None
        line_start = self.reverse_search_for(["Initial charges", "Initial moments and charges"])
        if line_start != LINE_NOT_FOUND:
            line_start += 2
            magmoms = np.zeros(self.n_atoms)
            for ll, line in enumerate(self.lines[line_start : line_start + self.n_atoms]):
                inp = line.split()
                if len(inp) == 4:
                    charges[ll] = float(inp[2])
                    magmoms = None
                else:
                    charges[ll] = float(inp[3])
                    magmoms[ll] = float(inp[2])

        self._cache["initial_charges"] = charges
        self._cache["initial_magnetic_moments"] = magmoms

    @property
    def is_md(self) -> bool:
        """Is the output for a molecular dynamics calculation?"""
        return self.reverse_search_for(["Complete information for previous time-step:"]) != LINE_NOT_FOUND

    @property
    def is_relaxation(self) -> bool:
        """Is the output for a relaxation?"""
        return self.reverse_search_for(["Geometry relaxation:"]) != LINE_NOT_FOUND

    def _parse_k_points(self) -> None:
        """Parse the list of k-points used in the calculation."""
        n_kpts = self.parse_scalar("n_kpts")
        if n_kpts is None:
            self._cache |= {"k_points": None, "k_point_weights": None}
            return
        n_kpts = int(n_kpts)

        line_start = self.reverse_search_for(["| K-points in task"])
        line_end = self.reverse_search_for(["| k-point:"])
        if LINE_NOT_FOUND in {line_start, line_end} or (line_end - line_start != n_kpts):
            self._cache |= {"k_points": None, "k_point_weights": None}
            return

        k_points = np.zeros((n_kpts, 3))
        k_point_weights = np.zeros(n_kpts)
        for kk, line in enumerate(self.lines[line_start + 1 : line_end + 1]):
            k_points[kk] = [float(inp) for inp in line.split()[4:7]]
            k_point_weights[kk] = float(line.split()[-1])

        self._cache |= {"k_points": k_points, "k_point_weights": k_point_weights}

    @property
    def n_atoms(self) -> int:
        """The number of atoms for the material."""
        n_atoms = self.parse_scalar("n_atoms")
        if n_atoms is None:
            raise AimsParseError("No information about the number of atoms in the header.")
        return int(n_atoms)

    @property
    def n_bands(self) -> int | None:
        """The number of Kohn-Sham states for the chunk."""
        line_start = self.reverse_search_for(SCALAR_PROPERTY_TO_LINE_KEY["n_bands"])

        if line_start == LINE_NOT_FOUND:
            raise AimsParseError("No information about the number of Kohn-Sham states in the header.")

        line = self.lines[line_start]
        if "| Number of Kohn-Sham states" in line:
            return int(line.split(":")[-1].strip().split()[0])

        return int(line.split()[-1].strip()[:-1])

    @property
    def n_electrons(self) -> int | None:
        """The number of electrons for the chunk."""
        line_start = self.reverse_search_for(SCALAR_PROPERTY_TO_LINE_KEY["n_electrons"])

        if line_start == LINE_NOT_FOUND:
            raise AimsParseError("No information about the number of electrons in the header.")

        line = self.lines[line_start]
        return int(float(line.split()[-2]))

    @property
    def n_k_points(self) -> int | None:
        """The number of k_ppoints for the calculation."""
        n_kpts = self.parse_scalar("n_kpts")
        if n_kpts is None:
            return None

        return int(n_kpts)

    @property
    def n_spins(self) -> int | None:
        """The number of spin channels for the chunk."""
        n_spins = self.parse_scalar("n_spins")
        if n_spins is None:
            raise AimsParseError("No information about the number of spin channels in the header.")
        return int(n_spins)

    @property
    def electronic_temperature(self) -> float:
        """The electronic temperature for the chunk."""
        line_start = self.reverse_search_for(SCALAR_PROPERTY_TO_LINE_KEY["electronic_temp"])
        # TARP: Default FHI-aims value
        if line_start == LINE_NOT_FOUND:
            return 0.00

        line = self.lines[line_start]
        return float(line.split("=")[-1].strip().split()[0])

    @property
    def k_points(self) -> Sequence[Vector3D]:
        """All k-points listed in the calculation."""
        if "k_points" not in self._cache:
            self._parse_k_points()

        return self._cache["k_points"]

    @property
    def k_point_weights(self) -> Sequence[float]:
        """The k-point weights for the calculation."""
        if "k_point_weights" not in self._cache:
            self._parse_k_points()

        return self._cache["k_point_weights"]

    @property
    def header_summary(self) -> dict[str, Any]:
        """Dictionary summarizing the information inside the header."""
        return {
            "initial_structure": self.initial_structure,
            "initial_lattice": self.initial_lattice,
            "is_relaxation": self.is_relaxation,
            "is_md": self.is_md,
            "n_atoms": self.n_atoms,
            "n_bands": self.n_bands,
            "n_electrons": self.n_electrons,
            "n_spins": self.n_spins,
            "electronic_temperature": self.electronic_temperature,
            "n_k_points": self.n_k_points,
            "k_points": self.k_points,
            "k_point_weights": self.k_point_weights,
        }

    @property
    def metadata_summary(self) -> dict[str, list[str] | str | None]:
        """Dictionary containing all metadata for FHI-aims build."""
        return {
            "commit_hash": self.commit_hash,
            "aims_uuid": self.aims_uuid,
            "version_number": self.version_number,
            "fortran_compiler": self.fortran_compiler,
            "c_compiler": self.c_compiler,
            "fortran_compiler_flags": self.fortran_compiler_flags,
            "c_compiler_flags": self.c_compiler_flags,
            "build_type": self.build_type,
            "linked_against": self.linked_against,
        }


class AimsOutCalcChunk(AimsOutChunk):
    """A part of the aims.out file corresponding to a single structure."""

    def __init__(self, lines: list[str], header: AimsOutHeaderChunk) -> None:
        """Construct the AimsOutCalcChunk.

        Args:
            lines (list[str]): The lines used for the structure
            header (.AimsOutHeaderChunk):  A summary of the relevant information from
                the aims.out header
        """
        super().__init__(lines)
        self._header = header.header_summary
        self._cache: dict[str, Any] = {}

    def _parse_structure(self) -> Structure | Molecule:
        """Parse a structure object from the file.

        For the given section of the aims output file generate the
        calculated structure.

        Returns:
            The structure or molecule for the calculation
        """
        species, coords, velocities, lattice = self._parse_lattice_atom_pos()

        site_properties: dict[str, Sequence[Any]] = {}
        if len(velocities) > 0:
            site_properties["velocity"] = np.array(velocities)

        results = self.results
        site_prop_keys = {
            "forces": "force",
            "stresses": "atomic_virial_stress",
            "hirshfeld_charges": "hirshfeld_charge",
            "hirshfeld_volumes": "hirshfeld_volume",
            "hirshfeld_atomic_dipoles": "hirshfeld_atomic_dipole",
            "mulliken_charges": "charge",
            "mulliken_spins": "magmom",
        }
        properties = {prop: results[prop] for prop in results if prop not in site_prop_keys}
        for prop, site_key in site_prop_keys.items():
            if prop in results:
                site_properties[site_key] = results[prop]

        if ((magmom := site_properties.get("magmom")) is not None) and np.abs(
            np.sum(magmom) - properties["magmom"]
        ) < 1e-3:
            warnings.warn(
                "Total magnetic moment and sum of Mulliken spins are not consistent",
                UserWarning,
                stacklevel=1,
            )

        if lattice is not None:
            return Structure(
                lattice,
                species,
                coords,
                site_properties=site_properties,
                properties=properties,
                coords_are_cartesian=True,
            )

        return Molecule(
            species,
            coords,
            site_properties=site_properties,
            properties=properties,
        )

    def _parse_lattice_atom_pos(
        self,
    ) -> tuple[list[str], list[Vector3D], list[Vector3D], Lattice | None]:
        """Parse the lattice and atomic positions of the structure.

        Returns:
            list[str]: The species symbols for the atoms in the structure
            list[Vector3D]: The Cartesian coordinates of the atoms
            list[Vector3D]: The velocities of the atoms
            Lattice or None: The Lattice for the structure
        """
        lattice_vectors = []
        velocities: list[Vector3D] = []
        species: list[str] = []
        coords: list[Vector3D] = []

        start_keys = [
            "Atomic structure (and velocities) as used in the preceding time step",
            "Updated atomic structure",
            "Atomic structure that was used in the preceding time step of the wrapper",
        ]
        line_start = self.reverse_search_for(start_keys)
        if line_start == LINE_NOT_FOUND:
            species = [sp.symbol for sp in self.initial_structure.species]
            coords = self.initial_structure.cart_coords.tolist()
            velocities = list(self.initial_structure.site_properties.get("velocity", []))
            lattice = self.initial_lattice

            return species, coords, velocities, lattice

        line_start += 1

        line_end = self.reverse_search_for(
            ['Writing the current geometry to file "geometry.in.next_step"'],
            line_start,
        )
        if line_end == LINE_NOT_FOUND:
            line_end = len(self.lines)

        for line in self.lines[line_start:line_end]:
            if "lattice_vector   " in line:
                lattice_vectors.append([float(inp) for inp in line.split()[1:]])
            elif "atom   " in line:
                line_split = line.split()
                species.append(line_split[4])
                coords.append(cast(Tuple3Floats, tuple(float(inp) for inp in line_split[1:4])))
            elif "velocity   " in line:
                velocities.append(cast(Tuple3Floats, tuple(float(inp) for inp in line.split()[1:4])))

        lattice = Lattice(lattice_vectors) if len(lattice_vectors) == 3 else None
        return species, coords, velocities, lattice

    @property
    def species(self) -> list[str]:
        """The list of atomic symbols for all atoms in the structure."""
        if "species" not in self._cache:
            (
                self._cache["species"],
                self._cache["coords"],
                self._cache["velocities"],
                self._cache["lattice"],
            ) = self._parse_lattice_atom_pos()
        return self._cache["species"]

    @property
    def coords(self) -> list[Vector3D]:
        """The cartesian coordinates of the atoms."""
        if "coords" not in self._cache:
            (
                self._cache["species"],
                self._cache["coords"],
                self._cache["velocities"],
                self._cache["lattice"],
            ) = self._parse_lattice_atom_pos()
        return self._cache["coords"]

    @property
    def velocities(self) -> list[Vector3D]:
        """The velocities of the atoms."""
        if "velocities" not in self._cache:
            (
                self._cache["species"],
                self._cache["coords"],
                self._cache["velocities"],
                self._cache["lattice"],
            ) = self._parse_lattice_atom_pos()
        return self._cache["velocities"]

    @property
    def lattice(self) -> Lattice:
        """The Lattice object for the structure."""
        if "lattice" not in self._cache:
            (
                self._cache["species"],
                self._cache["coords"],
                self._cache["velocities"],
                self._cache["lattice"],
            ) = self._parse_lattice_atom_pos()
        return self._cache["lattice"]

    @property
    def forces(self) -> np.ndarray | None:
        """The forces from the aims.out file."""
        line_start = self.reverse_search_for(["Total atomic forces"])
        if line_start == LINE_NOT_FOUND:
            return None

        line_start += 1

        return np.array(
            [[float(inp) for inp in line.split()[-3:]] for line in self.lines[line_start : line_start + self.n_atoms]]
        )

    @property
    def stresses(self) -> np.ndarray | None:
        """The stresses from the aims.out file and convert to kBar."""
        line_start = self.reverse_search_for(["Per atom stress (eV) used for heat flux calculation"])
        if line_start == LINE_NOT_FOUND:
            return None
        line_start += 3
        stresses = []
        for line in self.lines[line_start : line_start + self.n_atoms]:
            xx, yy, zz, xy, xz, yz = (float(d) for d in line.split()[2:8])
            stresses.append(Tensor.from_voigt([xx, yy, zz, yz, xz, xy]))

        return np.array(stresses) * EV_PER_A3_TO_KBAR

    @property
    def stress(self) -> Matrix3D | None:
        """The stress from the aims.out file and convert to kBar."""
        line_start = self.reverse_search_for(
            ["Analytical stress tensor - Symmetrized", "Numerical stress tensor"]
        )  # Offset to relevant lines
        if line_start == LINE_NOT_FOUND:
            return None

        stress = [[float(inp) for inp in line.split()[2:5]] for line in self.lines[line_start + 5 : line_start + 8]]
        return np.array(stress) * EV_PER_A3_TO_KBAR

    @property
    def is_metallic(self) -> bool:
        """Is the system is metallic."""
        line_start = self.reverse_search_for(
            ["material is metallic within the approximate finite broadening function (occupation_type)"]
        )
        return line_start != LINE_NOT_FOUND

    @property
    def energy(self) -> float:
        """The energy from the aims.out file."""
        if self.initial_lattice is not None and self.is_metallic:
            line_ind = self.reverse_search_for(["Total energy corrected"])
        else:
            line_ind = self.reverse_search_for(["Total energy uncorrected"])
        if line_ind == LINE_NOT_FOUND:
            raise AimsParseError("No energy is associated with the structure.")

        return float(self.lines[line_ind].split()[5])

    @property
    def dipole(self) -> Vector3D | None:
        """The electric dipole moment from the aims.out file."""
        line_start = self.reverse_search_for(["Total dipole moment [eAng]"])
        if line_start == LINE_NOT_FOUND:
            return None

        line = self.lines[line_start]
        return np.array([float(inp) for inp in line.split()[6:9]])

    @property
    def dielectric_tensor(self) -> Matrix3D | None:
        """The dielectric tensor from the aims.out file."""
        line_start = self.reverse_search_for(["PARSE DFPT_dielectric_tensor"])
        if line_start == LINE_NOT_FOUND:
            return None

        # we should find the tensor in the next three lines:
        lines = self.lines[line_start + 1 : line_start + 4]

        # make ndarray and return
        return np.array([np.fromstring(line, sep=" ") for line in lines])

    @property
    def polarization(self) -> Vector3D | None:
        """The polarization vector from the aims.out file."""
        line_start = self.reverse_search_for(["| Cartesian Polarization"])
        if line_start == LINE_NOT_FOUND:
            return None
        line = self.lines[line_start]
        return np.array([float(s) for s in line.split()[-3:]])

    def _parse_homo_lumo(self) -> dict[str, float]:
        """Parse the HOMO/LUMO values and get band gap if periodic."""
        line_start = self.reverse_search_for(["Highest occupied state (VBM)"])
        homo = float(self.lines[line_start].split(" at ")[1].split("eV")[0].strip())

        line_start = self.reverse_search_for(["Lowest unoccupied state (CBM)"])
        lumo = float(self.lines[line_start].split(" at ")[1].split("eV")[0].strip())

        line_start = self.reverse_search_for(["verall HOMO-LUMO gap"])
        homo_lumo_gap = float(self.lines[line_start].split(":")[1].split("eV")[0].strip())

        line_start = self.reverse_search_for(["Smallest direct gap"])
        if line_start == LINE_NOT_FOUND:
            return {
                "vbm": homo,
                "cbm": lumo,
                "gap": homo_lumo_gap,
                "direct_gap": homo_lumo_gap,
            }

        direct_gap = float(self.lines[line_start].split(":")[1].split("eV")[0].strip())
        return {
            "vbm": homo,
            "cbm": lumo,
            "gap": homo_lumo_gap,
            "direct_gap": direct_gap,
        }

    def _parse_hirshfeld(
        self,
    ) -> None:
        """Parse the Hirshfled charges volumes, and dipole moments."""
        line_start = self.reverse_search_for(["Performing Hirshfeld analysis of fragment charges and moments."])
        if line_start == LINE_NOT_FOUND:
            self._cache |= {
                "hirshfeld_charges": None,
                "hirshfeld_volumes": None,
                "hirshfeld_atomic_dipoles": None,
                "hirshfeld_dipole": None,
            }
            return

        line_inds = self.search_for_all("Hirshfeld charge", line_start, -1)
        hirshfeld_charges = np.array([float(self.lines[ind].split(":")[1]) for ind in line_inds])

        line_inds = self.search_for_all("Hirshfeld volume", line_start, -1)
        hirshfeld_volumes = np.array([float(self.lines[ind].split(":")[1]) for ind in line_inds])

        line_inds = self.search_for_all("Hirshfeld dipole vector", line_start, -1)
        hirshfeld_atomic_dipoles = np.array(
            [[float(inp) for inp in self.lines[ind].split(":")[1].split()] for ind in line_inds]
        )

        if self.lattice is None:
            hirshfeld_dipole = np.sum(
                hirshfeld_charges.reshape((-1, 1)) * self.coords,
                axis=1,
            )
        else:
            hirshfeld_dipole = None

        self._cache |= {
            "hirshfeld_charges": hirshfeld_charges,
            "hirshfeld_volumes": hirshfeld_volumes,
            "hirshfeld_atomic_dipoles": hirshfeld_atomic_dipoles,
            "hirshfeld_dipole": hirshfeld_dipole,
        }

    def _parse_mulliken(
        self,
    ) -> None:
        """Parse the Mulliken charges and spins."""
        line_start = self.reverse_search_for(["Performing Mulliken charge analysis"])
        if line_start == LINE_NOT_FOUND:
            self._cache.update(mulliken_charges=None, mulliken_spins=None)
            return

        line_start = self.reverse_search_for(["Summary of the per-atom charge analysis"])
        mulliken_charges = np.array(
            [float(self.lines[ind].split()[3]) for ind in range(line_start + 3, line_start + 3 + self.n_atoms)]
        )

        line_start = self.reverse_search_for(["Summary of the per-atom spin analysis"])
        if line_start == LINE_NOT_FOUND:
            mulliken_spins = None
        else:
            mulliken_spins = np.array(
                [float(self.lines[ind].split()[2]) for ind in range(line_start + 3, line_start + 3 + self.n_atoms)]
            )

        self._cache.update(
            {
                "mulliken_charges": mulliken_charges,
                "mulliken_spins": mulliken_spins,
            }
        )

    @property
    def structure(self) -> Structure | Molecule:
        """The pytmagen SiteCollection of the chunk."""
        if "structure" not in self._cache:
            self._cache["structure"] = self._parse_structure()
        return self._cache["structure"]

    @property
    def results(self) -> dict[str, Any]:
        """Convert an AimsOutChunk to a Results Dictionary."""
        results = {
            "energy": self.energy,
            "free_energy": self.free_energy,
            "forces": self.forces,
            "stress": self.stress,
            "stresses": self.stresses,
            "magmom": self.magmom,
            "dipole": self.dipole,
            "fermi_energy": self.E_f,
            "n_iter": self.n_iter,
            "mulliken_charges": self.mulliken_charges,
            "mulliken_spins": self.mulliken_spins,
            "hirshfeld_charges": self.hirshfeld_charges,
            "hirshfeld_dipole": self.hirshfeld_dipole,
            "hirshfeld_volumes": self.hirshfeld_volumes,
            "hirshfeld_atomic_dipoles": self.hirshfeld_atomic_dipoles,
            "dielectric_tensor": self.dielectric_tensor,
            "polarization": self.polarization,
            "vbm": self.vbm,
            "cbm": self.cbm,
            "gap": self.gap,
            "direct_gap": self.direct_gap,
        }

        return {key: value for key, value in results.items() if value is not None}

    # Properties from the aims.out header
    @property
    def initial_structure(self) -> Structure | Molecule:
        """The initial structure for the calculation."""
        return self._header["initial_structure"]

    @property
    def initial_lattice(self) -> Lattice | None:
        """The initial Lattice of the structure."""
        return self._header["initial_lattice"]

    @property
    def n_atoms(self) -> int:
        """The number of atoms in the structure."""
        return self._header["n_atoms"]

    @property
    def n_bands(self) -> int:
        """The number of Kohn-Sham states for the chunk."""
        return self._header["n_bands"]

    @property
    def n_electrons(self) -> int:
        """The number of electrons for the chunk."""
        return self._header["n_electrons"]

    @property
    def n_spins(self) -> int:
        """The number of spin channels for the chunk."""
        return self._header["n_spins"]

    @property
    def electronic_temperature(self) -> float:
        """The electronic temperature for the chunk."""
        return self._header["electronic_temperature"]

    @property
    def n_k_points(self) -> int:
        """The number of k_ppoints for the calculation."""
        return self._header["n_k_points"]

    @property
    def k_points(self) -> Sequence[Vector3D]:
        """All k-points listed in the calculation."""
        return self._header["k_points"]

    @property
    def k_point_weights(self) -> Sequence[float]:
        """The k-point weights for the calculation."""
        return self._header["k_point_weights"]

    @property
    def free_energy(self) -> float | None:
        """The free energy of the calculation."""
        return self.parse_scalar("free_energy")

    @property
    def n_iter(self) -> int | None:
        """The number of steps needed to converge the SCF cycle for the chunk."""
        val = self.parse_scalar("number_of_iterations")
        if val is not None:
            return int(val)
        return None

    @property
    def magmom(self) -> float | None:
        """The magnetic moment of the structure."""
        return self.parse_scalar("magnetic_moment")

    @property
    def E_f(self) -> float | None:
        """The Fermi energy."""
        return self.parse_scalar("fermi_energy")

    @property
    def converged(self) -> bool:
        """True if the calculation is converged."""
        return (len(self.lines) > 0) and ("Have a nice day." in self.lines[-5:])

    @property
    def mulliken_charges(self) -> Sequence[float] | None:
        """The Mulliken charges of the system"""
        if "mulliken_charges" not in self._cache:
            self._parse_mulliken()
        return self._cache["mulliken_charges"]

    @property
    def mulliken_spins(self) -> Sequence[float] | None:
        """The Mulliken spins of the system"""
        if "mulliken_spins" not in self._cache:
            self._parse_mulliken()
        return self._cache["mulliken_spins"]

    @property
    def hirshfeld_charges(self) -> Sequence[float] | None:
        """The Hirshfeld charges of the system."""
        if "hirshfeld_charges" not in self._cache:
            self._parse_hirshfeld()
        return self._cache["hirshfeld_charges"]

    @property
    def hirshfeld_atomic_dipoles(self) -> Sequence[Vector3D] | None:
        """The Hirshfeld atomic dipoles of the system."""
        if "hirshfeld_atomic_dipoles" not in self._cache:
            self._parse_hirshfeld()
        return self._cache["hirshfeld_atomic_dipoles"]

    @property
    def hirshfeld_volumes(self) -> Sequence[float] | None:
        """The Hirshfeld atomic dipoles of the system."""
        if "hirshfeld_volumes" not in self._cache:
            self._parse_hirshfeld()
        return self._cache["hirshfeld_volumes"]

    @property
    def hirshfeld_dipole(self) -> None | Vector3D:
        """The Hirshfeld dipole of the system."""
        if "hirshfeld_dipole" not in self._cache:
            self._parse_hirshfeld()

        return self._cache["hirshfeld_dipole"]

    @property
    def vbm(self) -> float:
        """The valance band maximum."""
        return self._parse_homo_lumo()["vbm"]

    @property
    def cbm(self) -> float:
        """The conduction band minimnum."""
        return self._parse_homo_lumo()["cbm"]

    @property
    def gap(self) -> float:
        """The band gap."""
        return self._parse_homo_lumo()["gap"]

    @property
    def direct_gap(self) -> float:
        """The direct bandgap."""
        return self._parse_homo_lumo()["direct_gap"]


def get_lines(content: str | TextIOWrapper) -> list[str]:
    """Get a list of lines from a str or file of content.

    Args:
        content: the content of the file to parse

    Returns:
        The list of lines
    """
    if isinstance(content, str):
        return [line.strip() for line in content.split("\n")]
    return [line.strip() for line in content.readlines()]


def get_header_chunk(content: str | TextIOWrapper) -> AimsOutHeaderChunk:
    """Get the header chunk for an output.

    Args:
        content (str or TextIOWrapper): the content to parse

    Returns:
        The AimsHeaderChunk of the file
    """
    lines = get_lines(content)
    header = []
    stopped = False
    # Stop the header once the first SCF cycle begins
    for line in lines:
        header.append(line)
        if (
            "Convergence:    q app. |  density  | eigen (eV) | Etot (eV)" in line
            or "Begin self-consistency iteration #" in line
        ):
            stopped = True
            break

    if not stopped:
        raise ParseError("No SCF steps present, calculation failed at setup.")

    return AimsOutHeaderChunk(header)


def get_aims_out_chunks(content: str | TextIOWrapper, header_chunk: AimsOutHeaderChunk) -> Generator:
    """Yield unprocessed chunks (header, lines) for each AimsOutChunk image.

    Args:
        content (str or TextIOWrapper): the content to parse
        header_chunk (AimsOutHeaderChunk): The AimsOutHeader for the calculation

    Yields:
        The next AimsOutChunk
    """
    lines = get_lines(content)[len(header_chunk.lines) :]
    if len(lines) == 0:
        return

    # If the calculation is relaxation the updated structural information
    # occurs before the re-initialization
    if header_chunk.is_relaxation:
        chunk_end_line = "Geometry optimization: Attempting to predict improved coordinates."
    else:
        chunk_end_line = "Begin self-consistency loop: Re-initialization"

    # If SCF is not converged then do not treat the next chunk_end_line as a
    # new chunk until after the SCF is re-initialized
    ignore_chunk_end_line = False
    line_iter = iter(lines)
    while True:
        try:
            line = next(line_iter).strip()  # Raises StopIteration on empty file
        except StopIteration:
            break

        chunk_lines = []
        while chunk_end_line not in line or ignore_chunk_end_line:
            chunk_lines.append(line)
            # If SCF cycle not converged or numerical stresses are requested,
            # don't end chunk on next Re-initialization
            patterns = [
                ("Self-consistency cycle not yet converged - restarting mixer to attempt better convergence."),
                (
                    "Components of the stress tensor (for mathematical "
                    "background see comments in numerical_stress.f90)."
                ),
                "Calculation of numerical stress completed",
            ]
            if any(pattern in line for pattern in patterns):
                ignore_chunk_end_line = True
            elif "Begin self-consistency loop: Re-initialization" in line:
                ignore_chunk_end_line = False

            try:
                line = next(line_iter).strip()
            except StopIteration:
                break
        yield AimsOutCalcChunk(chunk_lines, header_chunk)


def check_convergence(chunks: list[AimsOutCalcChunk], non_convergence_ok: bool = False) -> bool:
    """Check if the aims output file is for a converged calculation.

    Args:
        chunks(list[.AimsOutCalcChunk]): The list of chunks for the aims calculations
        non_convergence_ok(bool): True if it is okay for the calculation to not be converged
        chunks: list[AimsOutCalcChunk]:
        non_convergence_ok: bool:  (Default value = False)

    Returns:
        True if the calculation is converged
    """
    if not non_convergence_ok and not chunks[-1].converged:
        raise ParseError("The calculation did not complete successfully")
    return True


def read_aims_header_info_from_content(
    content: str,
) -> tuple[dict[str, list[str] | None | str], dict[str, Any]]:
    """Read the FHI-aims header information.

    Args:
      content (str): The content of the output file to check

    Returns:
        The metadata for the header of the aims calculation
    """
    header_chunk = get_header_chunk(content)
    return header_chunk.metadata_summary, header_chunk.header_summary


def read_aims_header_info(
    filename: str | Path,
) -> tuple[dict[str, None | list[str] | str], dict[str, Any]]:
    """Read the FHI-aims header information.

    Args:
      filename(str or Path): The file to read

    Returns:
        The metadata for the header of the aims calculation
    """
    content = None
    for path in [Path(filename), Path(f"{filename}.gz")]:
        if not path.exists():
            continue
        if path.suffix == ".gz":
            with gzip.open(filename, mode="rt") as file:
                content = file.read()
        else:
            with open(filename) as file:
                content = file.read()

    if content is None:
        raise FileNotFoundError(f"The requested output file {filename} does not exist.")

    return read_aims_header_info_from_content(content)


def read_aims_output_from_content(
    content: str, index: int | slice = -1, non_convergence_ok: bool = False
) -> Structure | Molecule | Sequence[Structure | Molecule]:
    """Read and aims output file from the content of a file.

    Args:
      content (str): The content of the file to read
      index: int | slice:  (Default value = -1)
      non_convergence_ok: bool:  (Default value = False)

    Returns:
        The list of images to get
    """
    header_chunk = get_header_chunk(content)
    chunks = list(get_aims_out_chunks(content, header_chunk))
    if header_chunk.is_relaxation and any("Final atomic structure:" in line for line in chunks[-1].lines):
        chunks[-2].lines += chunks[-1].lines
        chunks = chunks[:-1]

    check_convergence(chunks, non_convergence_ok)
    # Relaxations have an additional footer chunk due to how it is split
    images = [chunk.structure for chunk in chunks]
    return images[index]


def read_aims_output(
    filename: str | Path,
    index: int | slice = -1,
    non_convergence_ok: bool = False,
) -> Structure | Molecule | Sequence[Structure | Molecule]:
    """Import FHI-aims output files with all data available.

    Includes all structures for relaxations and MD runs with FHI-aims

    Args:
        filename(str or Path): The file to read
        index(int or slice): The index of the images to read
        non_convergence_ok(bool): True if the calculations do not have to be converged

    Returns:
        The list of images to get
    """
    content = None
    for path in [Path(filename), Path(f"{filename}.gz")]:
        if not path.exists():
            continue
        if path.suffix == ".gz":
            with gzip.open(path, mode="rt") as file:
                content = file.read()
        else:
            with open(path) as file:
                content = file.read()

    if content is None:
        raise FileNotFoundError(f"The requested output file {filename} does not exist.")

    return read_aims_output_from_content(content, index, non_convergence_ok)
