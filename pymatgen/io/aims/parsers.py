"""AIMS output parser, taken from ASE with modifications."""
from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from pymatgen.core import Lattice, Molecule, Structure

if TYPE_CHECKING:
    from collections.abc import Sequence

    from emmet.core.math import Vector3D

    from pymatgen.core.structure import SiteCollection

# TARP: Originally an object, but type hinting needs this to be an int
LINE_NOT_FOUND = -1000
EV_PER_A3_TO_KBAR = 1.60217653e-19 * 1e22


class ParseError(Exception):
    """Parse error during reading of a file"""


class AimsParseError(Exception):
    """Exception raised if an error occurs when parsing an Aims output file."""

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


# Read aims.out files
scalar_property_to_line_key = {
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

    Parameters
    ----------
    lines: list[str]
        The list of all lines in the aims.out file
    """

    lines: list[str] = field(default_factory=list)

    def reverse_search_for(self, keys: list[str], line_start: int = 0) -> int:
        """Find the last time one of the keys appears in self.lines.

        Parameters
        ----------
        keys: list[str]
            The key strings to search for in self.lines
        line_start: int
            The lowest index to search for in self.lines

        Returns
        -------
        int
            The last time one of the keys appears in self.lines
        """
        for ll, line in enumerate(self.lines[line_start:][::-1]):
            if any(key in line for key in keys):
                return len(self.lines) - ll - 1

        return LINE_NOT_FOUND

    def search_for_all(self, key: str, line_start: int = 0, line_end: int = -1) -> list[int]:
        """Find the all times the key appears in self.lines.

        Parameters
        ----------
        key: str
            The key string to search for in self.lines
        line_start: int
            The first line to start the search from
        line_end: int
            The last line to end the search at

        Returns
        -------
        list[ints]
            All times the key appears in the lines
        """
        line_index = []
        for ll, line in enumerate(self.lines[line_start:line_end]):
            if key in line:
                line_index.append(ll + line_start)
        return line_index

    def parse_scalar(self, property: str) -> float | None:
        """Parse a scalar property from the chunk.

        Parameters
        ----------
        property: str
            The property key to parse

        Returns
        -------
        float
            The scalar value of the property
        """
        line_start = self.reverse_search_for(scalar_property_to_line_key[property])

        if line_start == LINE_NOT_FOUND:
            return None

        line = self.lines[line_start]
        return float(line.split(":")[-1].strip().split()[0])


@dataclass
class AimsOutHeaderChunk(AimsOutChunk):
    """The header of the aims.out file containing general information.

    Parameters
    ----------
    lines: list[str]
        The list of all lines in the aims.out file
    _cache: dict[str, Any]
        The cache storing previously calculated results
    """

    lines: list[str] = field(default_factory=list)
    _cache: dict[str, Any] = field(default_factory=dict)

    @property
    def commit_hash(self):
        """Get the commit hash for the FHI-aims version."""
        line_start = self.reverse_search_for(["Commit number"])
        if line_start == LINE_NOT_FOUND:
            raise AimsParseError("This file does not appear to be an aims-output file")

        return self.lines[line_start].split(":")[1].strip()

    @property
    def aims_uuid(self):
        """Get the aims-uuid for the calculation."""
        line_start = self.reverse_search_for(["aims_uuid"])
        if line_start == LINE_NOT_FOUND:
            raise AimsParseError("This file does not appear to be an aims-output file")

        return self.lines[line_start].split(":")[1].strip()

    @property
    def version_number(self):
        """Get the commit hash for the FHI-aims version."""
        line_start = self.reverse_search_for(["FHI-aims version"])
        if line_start == LINE_NOT_FOUND:
            raise AimsParseError("This file does not appear to be an aims-output file")

        return self.lines[line_start].split(":")[1].strip()

    @property
    def fortran_compiler(self):
        """Get the fortran compiler used to make FHI-aims."""
        line_start = self.reverse_search_for(["Fortran compiler      :"])
        if line_start == LINE_NOT_FOUND:
            raise AimsParseError("This file does not appear to be an aims-output file")

        return self.lines[line_start].split(":")[1].split("/")[-1].strip()

    @property
    def c_compiler(self):
        """Get the C compiler used to make FHI-aims."""
        line_start = self.reverse_search_for(["C compiler            :"])
        if line_start == LINE_NOT_FOUND:
            return None

        return self.lines[line_start].split(":")[1].split("/")[-1].strip()

    @property
    def fortran_compiler_flags(self):
        """Get the fortran compiler flags used to make FHI-aims."""
        line_start = self.reverse_search_for(["Fortran compiler flags"])
        if line_start == LINE_NOT_FOUND:
            raise AimsParseError("This file does not appear to be an aims-output file")

        return self.lines[line_start].split(":")[1].strip()

    @property
    def c_compiler_flags(self):
        """Get the C compiler flags used to make FHI-aims."""
        line_start = self.reverse_search_for(["C compiler flags"])
        if line_start == LINE_NOT_FOUND:
            return None

        return self.lines[line_start].split(":")[1].strip()

    @property
    def build_type(self):
        """Get the optional build flags passed to cmake."""
        line_end = self.reverse_search_for(["Linking against:"])
        line_inds = self.search_for_all("Using", line_end=line_end)

        return [" ".join(self.lines[ind].split()[1:]).strip() for ind in line_inds]

    @property
    def linked_against(self):
        """Get all libraries used to link the FHI-aims executable."""
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
    def initial_lattice(self):
        """Parse the initial lattice vectors from the aims.out file."""
        line_start = self.reverse_search_for(["| Unit cell:"])
        if line_start == LINE_NOT_FOUND:
            return None

        return Lattice(
            np.array(
                [[float(inp) for inp in line.split()[-3:]] for line in self.lines[line_start + 1 : line_start + 4]]
            )
        )

    @property
    def initial_structure(self):
        """Create an SiteCollection object for the initial structure.

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
            return Structure(lattice, species, coords, np.sum(self.initial_charges), site_properties=site_properties)

        return Molecule(species, coords, np.sum(self.initial_charges), site_properties=site_properties)

    @property
    def initial_charges(self):
        if "initial_charges" not in self._cache:
            self._parse_initial_charges_and_moments()
        return self._cache["initial_charges"]

    @property
    def initial_magnetic_moments(self):
        if "initial_magnetic_moments" not in self._cache:
            self._parse_initial_charges_and_moments()
        return self._cache["initial_magnetic_moments"]

    def _parse_initial_charges_and_moments(self):
        """Parse the initial charges and magnetic moments

        Once parsed store them in the cache
        """
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
    def is_md(self):
        """Determine if calculation is a molecular dynamics calculation."""
        return self.reverse_search_for(["Complete information for previous time-step:"]) != LINE_NOT_FOUND

    @property
    def is_relaxation(self):
        """Determine if the calculation is a geometry optimization or not."""
        return self.reverse_search_for(["Geometry relaxation:"]) != LINE_NOT_FOUND

    def _parse_k_points(self):
        """Get the list of k-points used in the calculation."""
        n_kpts = self.parse_scalar("n_kpts")
        if n_kpts is None:
            return {
                "k_points": None,
                "k_point_weights": None,
            }
        n_kpts = int(n_kpts)

        line_start = self.reverse_search_for(["| K-points in task"])
        line_end = self.reverse_search_for(["| k-point:"])
        if (line_start == LINE_NOT_FOUND) or (line_end == LINE_NOT_FOUND) or (line_end - line_start != n_kpts):
            return {
                "k_points": None,
                "k_point_weights": None,
            }

        k_points = np.zeros((n_kpts, 3))
        k_point_weights = np.zeros(n_kpts)
        for kk, line in enumerate(self.lines[line_start + 1 : line_end + 1]):
            k_points[kk] = [float(inp) for inp in line.split()[4:7]]
            k_point_weights[kk] = float(line.split()[-1])

        self._cache.update(
            {
                "k_points": k_points,
                "k_point_weights": k_point_weights,
            }
        )
        return None

    @property
    def n_atoms(self) -> int:
        """The number of atoms for the material."""
        n_atoms = self.parse_scalar("n_atoms")
        if n_atoms is None:
            raise AimsParseError("No information about the number of atoms in the header.")
        return int(n_atoms)

    @property
    def n_bands(self):
        """The number of Kohn-Sham states for the chunk."""
        line_start = self.reverse_search_for(scalar_property_to_line_key["n_bands"])

        if line_start == LINE_NOT_FOUND:
            raise AimsParseError("No information about the number of Kohn-Sham states in the header.")

        line = self.lines[line_start]
        if "| Number of Kohn-Sham states" in line:
            return int(line.split(":")[-1].strip().split()[0])

        return int(line.split()[-1].strip()[:-1])

    @property
    def n_electrons(self):
        """The number of electrons for the chunk."""
        line_start = self.reverse_search_for(scalar_property_to_line_key["n_electrons"])

        if line_start == LINE_NOT_FOUND:
            raise AimsParseError("No information about the number of electrons in the header.")

        line = self.lines[line_start]
        return int(float(line.split()[-2]))

    @property
    def n_k_points(self):
        """The number of k_ppoints for the calculation."""
        n_kpts = self.parse_scalar("n_kpts")
        if n_kpts is None:
            return None

        return int(n_kpts)

    @property
    def n_spins(self):
        """The number of spin channels for the chunk."""
        n_spins = self.parse_scalar("n_spins")
        if n_spins is None:
            raise AimsParseError("No information about the number of spin channels in the header.")
        return int(n_spins)

    @property
    def electronic_temperature(self):
        """The electronic temperature for the chunk."""
        line_start = self.reverse_search_for(scalar_property_to_line_key["electronic_temp"])
        if line_start == LINE_NOT_FOUND:
            return 0.10

        line = self.lines[line_start]
        return float(line.split("=")[-1].strip().split()[0])

    @property
    def k_points(self):
        """All k-points listed in the calculation."""
        if "k_points" not in self._cache:
            self._parse_k_points()

        return self._cache["k_points"]

    @property
    def k_point_weights(self):
        """The k-point weights for the calculation."""
        if "k_point_weights" not in self._cache:
            self._parse_k_points()

        return self._cache["k_point_weights"]

    @property
    def header_summary(self):
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
    def metadata_summary(self) -> dict[str, str]:
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

    def __init__(self, lines, header):
        """Construct the AimsOutCalcChunk.

        Parameters
        ----------
        lines: list[str]
            The lines used for the structure
        header: .AimsOutHeaderChunk
            A summary of the relevant information from the aims.out header
        """
        super().__init__(lines)
        self._header = header.header_summary
        self._cache = {}

    def _parse_structure(self):
        """Parse a structure object from the file.

        For the given section of the aims output file generate the
        calculated structure.
        """
        lattice = self.lattice
        velocities = self.velocities
        species = self.species
        coords = self.coords

        site_properties = dict()
        if len(velocities) > 0:
            site_properties["velocity"] = np.array(velocities)

        results = self.results
        site_prop_keys = {
            "forces": "force",
            "stresses": "atomic_virial_stress",
            "hirshfeld_charges": "hirshfeld_charge",
            "hirshfeld_volumes": "hirshfeld_volume",
            "hirshfeld_atomic_dipoles": "hirshfeld_atomic_dipole",
        }
        properties = {prop: results[prop] for prop in results if prop not in site_prop_keys}
        for prop, site_key in site_prop_keys.items():
            if prop in results:
                site_properties[site_key] = results[prop]

        if len(site_properties.keys()) == 0:
            site_properties = None

        if lattice is not None:
            return Structure(lattice, species, coords, site_properties=site_properties, properties=properties)
        return Molecule(species, coords, site_properties=site_properties, properties=properties)

    @property
    def _parse_lattice_atom_pos(self) -> tuple[list[str], list[Vector3D], list[Vector3D], Lattice | None]:
        """Get the lattice of the structure

        Returns
        -------
        species: list[str]
            The species of all atoms
        coords: list[Vector3D]
            The cartesian coordinates of the atoms
        velocities: list[Vector3D]
            The velocities of all atoms
        lattice: Lattice or None
            The lattice of the system
        """
        lattice_vectors = []
        velocities = []
        species = []
        coords = []

        start_keys = [
            "Atomic structure (and velocities) as used in the preceding time step",
            "Updated atomic structure",
            "Atomic structure that was used in the preceding time step of the wrapper",
        ]
        line_start = self.reverse_search_for(start_keys)
        if line_start == LINE_NOT_FOUND:
            return self.initial_structure

        line_start += 1

        line_end = self.reverse_search_for(['Writing the current geometry to file "geometry.in.next_step"'], line_start)
        if line_end == LINE_NOT_FOUND:
            line_end = len(self.lines)

        for line in self.lines[line_start:line_end]:
            if "lattice_vector   " in line:
                lattice_vectors.append([float(inp) for inp in line.split()[1:]])
            elif "atom   " in line:
                line_split = line.split()
                species.append(line_split[4])
                coords.append([float(inp) for inp in line_split[1:4]])
            elif "velocity   " in line:
                velocities.append([float(inp) for inp in line.split()[1:]])

        return species, coords, velocities, Lattice(lattice_vectors)

    @property
    def species(self):
        if "species" not in self._cache:
            self._cache.update(self._parse_lattice_atom_pos())
        return self._cache["species"]

    @property
    def coords(self):
        if "coords" not in self._cache:
            self._cache.update(self._parse_lattice_atom_pos())
        return self._cache["coords"]

    @property
    def velocities(self):
        if "velocities" not in self._cache:
            self._cache.update(self._parse_lattice_atom_pos())
        return self._cache["velocities"]

    @property
    def lattice(self):
        if "lattice" not in self._cache:
            self._cache.update(self._parse_lattice_atom_pos())
        return self._cache["lattice"]

    @property
    def forces(self):
        """Parse the forces from the aims.out file."""
        line_start = self.reverse_search_for(["Total atomic forces"])
        if line_start == LINE_NOT_FOUND:
            return None

        line_start += 1

        return np.array(
            [[float(inp) for inp in line.split()[-3:]] for line in self.lines[line_start : line_start + self.n_atoms]]
        )

    @property
    def stresses(self):
        """Parse the stresses from the aims.out file and convert to kbar."""
        line_start = self.reverse_search_for(["Per atom stress (eV) used for heat flux calculation"])
        if line_start == LINE_NOT_FOUND:
            return None
        line_start += 3
        stresses = []
        for line in self.lines[line_start : line_start + self.n_atoms]:
            xx, yy, zz, xy, xz, yz = (float(d) for d in line.split()[2:8])
            stresses.append([xx, yy, zz, yz, xz, xy])

        return np.array(stresses) * EV_PER_A3_TO_KBAR

    @property
    def stress(self):
        """Parse the stress from the aims.out file and convert to kbar."""
        from ase.stress import full_3x3_to_voigt_6_stress

        line_start = self.reverse_search_for(
            [
                "Analytical stress tensor - Symmetrized",
                "Numerical stress tensor",
            ]
        )  # Offset to relevant lines
        if line_start == LINE_NOT_FOUND:
            return None

        stress = [[float(inp) for inp in line.split()[2:5]] for line in self.lines[line_start + 5 : line_start + 8]]
        return full_3x3_to_voigt_6_stress(stress) * EV_PER_A3_TO_KBAR

    @property
    def is_metallic(self):
        """Checks if the system is metallic."""
        line_start = self.reverse_search_for(
            ["material is metallic within the approximate finite broadening function (occupation_type)"]
        )
        return line_start != LINE_NOT_FOUND

    @property
    def energy(self):
        """Parse the energy from the aims.out file."""
        if self.initial_lattice is not None and self.is_metallic:
            line_ind = self.reverse_search_for(["Total energy corrected"])
        else:
            line_ind = self.reverse_search_for(["Total energy uncorrected"])
        if line_ind == LINE_NOT_FOUND:
            raise AimsParseError("No energy is associated with the structure.")

        return float(self.lines[line_ind].split()[5])

    @property
    def dipole(self):
        """Parse the electric dipole moment from the aims.out file."""
        line_start = self.reverse_search_for(["Total dipole moment [eAng]"])
        if line_start == LINE_NOT_FOUND:
            return None

        line = self.lines[line_start]
        return np.array([float(inp) for inp in line.split()[6:9]])

    @property
    def dielectric_tensor(self):
        """Parse the dielectric tensor from the aims.out file."""
        line_start = self.reverse_search_for(["PARSE DFPT_dielectric_tensor"])
        if line_start == LINE_NOT_FOUND:
            return None

        # we should find the tensor in the next three lines:
        lines = self.lines[line_start + 1 : line_start + 4]

        # make ndarray and return
        return np.array([np.fromstring(line, sep=" ") for line in lines])

    @property
    def polarization(self):
        """Parse the polarization vector from the aims.out file."""
        line_start = self.reverse_search_for(["| Cartesian Polarization"])
        if line_start == LINE_NOT_FOUND:
            return None
        line = self.lines[line_start]
        return np.array([float(s) for s in line.split()[-3:]])

    def _parse_homo_lumo(self):
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

    def _parse_hirshfeld(self):
        """Parse the Hirshfled charges volumes, and dipole moments."""
        line_start = self.reverse_search_for(["Performing Hirshfeld analysis of fragment charges and moments."])
        if line_start == LINE_NOT_FOUND:
            return {
                "charges": None,
                "volumes": None,
                "atomic_dipoles": None,
                "dipole": None,
            }

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
        return {
            "charges": hirshfeld_charges,
            "volumes": hirshfeld_volumes,
            "atomic_dipoles": hirshfeld_atomic_dipoles,
            "dipole": hirshfeld_dipole,
        }

    @property
    def structure(self) -> Structure | Molecule:
        """The pytmagen SiteCollection of the output file."""
        if "structure" not in self._cache:
            self._cache["structure"] = self._parse_structure()
        return self._cache["structure"]

    @property
    def results(self):
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
    def initial_structure(self):
        """Return the initial structure of the calculation."""
        return self._header["initial_structure"]

    @property
    def initial_lattice(self):
        """Return the initial lattice vectors for the structure."""
        return self._header["initial_lattice"]

    @property
    def n_atoms(self):
        """Return the number of atoms for the material."""
        return self._header["n_atoms"]

    @property
    def n_bands(self):
        """Return the number of Kohn-Sham states for the chunk."""
        return self._header["n_bands"]

    @property
    def n_electrons(self):
        """Return the number of electrons for the chunk."""
        return self._header["n_electrons"]

    @property
    def n_spins(self):
        """Return the number of spin channels for the chunk."""
        return self._header["n_spins"]

    @property
    def electronic_temperature(self):
        """Return the electronic temperature for the chunk."""
        return self._header["electronic_temperature"]

    @property
    def n_k_points(self):
        """Return the number of electrons for the chunk."""
        return self._header["n_k_points"]

    @property
    def k_points(self):
        """Return the number of spin channels for the chunk."""
        return self._header["k_points"]

    @property
    def k_point_weights(self):
        """Return tk_point_weights electronic temperature for the chunk."""
        return self._header["k_point_weights"]

    @property
    def free_energy(self):
        """Return the free energy for the chunk."""
        return self.parse_scalar("free_energy")

    @property
    def n_iter(self):
        """Return the number of SCF iterations.

        The number of steps needed to converge the SCF cycle for the chunk.
        """
        return self.parse_scalar("number_of_iterations")

    @property
    def magmom(self):
        """Return the magnetic moment for the chunk."""
        return self.parse_scalar("magnetic_moment")

    @property
    def E_f(self):
        """Return he Fermi energy for the chunk."""
        return self.parse_scalar("fermi_energy")

    @property
    def converged(self):
        """Return True if the chunk is a fully converged final structure."""
        return (len(self.lines) > 0) and ("Have a nice day." in self.lines[-5:])

    @property
    def hirshfeld_charges(self):
        """Return the Hirshfeld charges for the chunk."""
        return self._parse_hirshfeld()["charges"]

    @property
    def hirshfeld_atomic_dipoles(self):
        """Return the Hirshfeld atomic dipole moments for the chunk."""
        return self._parse_hirshfeld()["atomic_dipoles"]

    @property
    def hirshfeld_volumes(self):
        """Return the Hirshfeld volume for the chunk."""
        return self._parse_hirshfeld()["volumes"]

    @property
    def hirshfeld_dipole(self):
        """Return the Hirshfeld systematic dipole moment for the chunk."""
        if self.lattice is not None:
            return self._parse_hirshfeld()["dipole"]

        return None

    @property
    def vbm(self):
        """Return the HOMO (VBM) of the calculation."""
        return self._parse_homo_lumo()["vbm"]

    @property
    def cbm(self):
        """Return the LUMO (CBM) of the calculation."""
        return self._parse_homo_lumo()["cbm"]

    @property
    def gap(self):
        """Return the HOMO-LUMO gap (band gap) of the calculation."""
        return self._parse_homo_lumo()["gap"]

    @property
    def direct_gap(self):
        """Return the direct band gap of the calculation."""
        return self._parse_homo_lumo()["direct_gap"]


def get_header_chunk(fd):
    """Return the header information from the aims.out file."""
    header = []
    line = ""

    # Stop the header once the first SCF cycle begins
    while (
        "Convergence:    q app. |  density  | eigen (eV) | Etot (eV)" not in line
        and "Begin self-consistency iteration #" not in line
    ):
        try:
            line = next(fd).strip()  # Raises StopIteration on empty file
        except StopIteration:
            raise ParseError("No SCF steps present, calculation failed at setup.") from None

        header.append(line)
    return AimsOutHeaderChunk(header)


def get_aims_out_chunks(fd, header_chunk):
    """Yield unprocessed chunks (header, lines) for each AimsOutChunk image."""
    try:
        line = next(fd).strip()  # Raises StopIteration on empty file
    except StopIteration:
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
    while True:
        try:
            line = next(fd).strip()  # Raises StopIteration on empty file
        except StopIteration:
            break

        lines = []
        while chunk_end_line not in line or ignore_chunk_end_line:
            lines.append(line)
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
                line = next(fd).strip()
            except StopIteration:
                break
        yield AimsOutCalcChunk(lines, header_chunk)


def check_convergence(chunks: list[AimsOutCalcChunk], non_convergence_ok: bool = False) -> bool:
    """Check if the aims output file is for a converged calculation.

    Parameters
    ----------
    chunks: list[.AimsOutCalcChunk]
        The list of chunks for the aims calculations
    non_convergence_ok: bool
        True if it is okay for the calculation to not be converged

    Returns
    -------
    bool
        True if the calculation is converged
    """
    if not non_convergence_ok and not chunks[-1].converged:
        raise ParseError("The calculation did not complete successfully")
    return True


def read_aims_header_info(
    filename: str | Path,
) -> tuple[dict[str, str], dict[str, Any]]:
    """Read the FHI-aims header information.

    Parameters
    ----------
    filename: str or Path
        The file to read

    Returns
    -------
    The calculation metadata and the system summary
    """
    header_chunk = None
    for path in [Path(filename), Path(f"{filename}.gz")]:
        if not path.exists():
            continue
        if path.suffix == ".gz":
            with gzip.open(filename, "rt") as fd:
                header_chunk = get_header_chunk(fd)
        else:
            with open(filename) as fd:
                header_chunk = get_header_chunk(fd)

    if header_chunk is None:
        raise FileNotFoundError(f"The requested output file {filename} does not exist.")

    system_summary = header_chunk.header_summary
    metadata = header_chunk.metadata_summary
    return metadata, system_summary


def read_aims_output(
    filename: str | Path,
    index: int | slice = -1,
    non_convergence_ok: bool = False,
) -> SiteCollection | Sequence[SiteCollection]:
    """Import FHI-aims output files with all data available.

    Includes all structures for relaxations and MD runs with FHI-aims

    Parameters
    ----------
    filename: str or Path
        The file to read
    index: int or slice
        The index of the images to read
    non_convergence_ok: bool
        True if the calculations do not have to be converged

    Returns
    -------
    The selected Structure objects
    """
    chunks = None
    for path in [Path(filename), Path(f"{filename}.gz")]:
        if not path.exists():
            continue
        if path.suffix == ".gz":
            with gzip.open(path, "rt") as fd:
                header_chunk = get_header_chunk(fd)
                chunks = list(get_aims_out_chunks(fd, header_chunk))
        else:
            with open(path) as fd:
                header_chunk = get_header_chunk(fd)
                chunks = list(get_aims_out_chunks(fd, header_chunk))

    if chunks is None:
        raise FileNotFoundError(f"The requested output file {filename} does not exist.")

    check_convergence(chunks, non_convergence_ok)

    # Relaxations have an additional footer chunk due to how it is split
    if header_chunk.is_relaxation:
        images = [chunk.structure for chunk in chunks[:-1]]
    else:
        images = [chunk.structure for chunk in chunks]
    return images[index]
