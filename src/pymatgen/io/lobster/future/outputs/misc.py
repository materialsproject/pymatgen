from __future__ import annotations

import itertools
import re
from itertools import islice
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from pymatgen.electronic_structure.core import Orbital, Spin
from pymatgen.io.lobster.future.constants import LOBSTER_ORBITALS
from pymatgen.io.lobster.future.core import LobsterFile
from pymatgen.io.lobster.future.utils import make_json_compatible, parse_orbital_from_text
from pymatgen.io.lobster.future.versioning import version_processor
from pymatgen.io.vasp.outputs import VolumetricData

if TYPE_CHECKING:
    from typing import Any, Literal

    from numpy import floating

    from pymatgen.core.structure import Structure
    from pymatgen.io.lobster.future.types import LobsterMatrixData
    from pymatgen.util.typing import PathLike


class Wavefunction(LobsterFile):
    """Parser for wave function files from LOBSTER.

    Reads wave function files and creates VolumetricData objects.

    Attributes:
        grid (tuple[int, int, int]): Grid for the wave function [Nx+1, Ny+1, Nz+1].
        points (list[tuple[float, float, float]]): Points in real space.
        reals (list[float]): Real parts of the wave function.
        imaginaries (list[float]): Imaginary parts of the wave function.
        distances (list[float]): Distances to the first point in the wave function file.
        structure (Structure): Structure object associated with the calculation.
    """

    def __init__(
        self, filename: PathLike, structure: Structure, process_immediately: bool = True
    ) -> None:
        """Initialize the Wavefunction parser.

        Args:
            filename (PathLike): The wavecar file from LOBSTER.
            structure (Structure): The Structure object.
            process_immediately (bool): Whether to parse the file immediately. Defaults to True.
        """
        super().__init__(filename, process_immediately=process_immediately)

        self.structure = structure

    @version_processor()
    def parse_file(
        self,
    ) -> None:
        """Parse wave function file.

        Reads the wave function file and extracts grid, points, real and imaginary parts,
        and distances.

        Raises:
            ValueError: If the number of real or imaginary parts does not match the expected grid size.
        """
        lines = self.lines

        self.points = []
        self.distances = []
        self.reals = []
        self.imaginaries = []

        line_parts = lines[0].split()

        self.grid: tuple[int, int, int] = (
            int(line_parts[7]),
            int(line_parts[8]),
            int(line_parts[9]),
        )

        for line in lines[1:]:
            line_parts = line.split()
            if len(line_parts) >= 6:
                self.points.append(
                    (float(line_parts[0]), float(line_parts[1]), float(line_parts[2]))
                )
                self.distances.append(float(line_parts[3]))
                self.reals.append(float(line_parts[4]))
                self.imaginaries.append(float(line_parts[5]))

        if (
            len(self.reals) != self.grid[0] * self.grid[1] * self.grid[2]
            or len(self.imaginaries) != self.grid[0] * self.grid[1] * self.grid[2]
        ):
            raise ValueError("Something went wrong while reading the file")

    def set_volumetric_data(
        self, grid: tuple[int, int, int], structure: Structure
    ) -> None:
        """Create VolumetricData instances for real, imaginary, and density parts.

        Args:
            grid (tuple[int, int, int]): Grid on which wavefunction was calculated.
            structure (Structure): Structure object.

        Raises:
            ValueError: If the wavefunction file does not contain all relevant points.
        """
        Nx = grid[0] - 1
        Ny = grid[1] - 1
        Nz = grid[2] - 1
        a = structure.lattice.matrix[0]
        b = structure.lattice.matrix[1]
        c = structure.lattice.matrix[2]
        new_x = []
        new_y = []
        new_z = []
        new_real = []
        new_imaginary = []
        new_density = []

        for runner, (x, y, z) in enumerate(
            itertools.product(range(Nx + 1), range(Ny + 1), range(Nz + 1))
        ):
            x_here = x / float(Nx) * a[0] + y / float(Ny) * b[0] + z / float(Nz) * c[0]
            y_here = x / float(Nx) * a[1] + y / float(Ny) * b[1] + z / float(Nz) * c[1]
            z_here = x / float(Nx) * a[2] + y / float(Ny) * b[2] + z / float(Nz) * c[2]

            if x != Nx and y != Ny and z != Nz:
                if (
                    not np.isclose(self.points[runner][0], x_here, 1e-3)
                    and not np.isclose(self.points[runner][1], y_here, 1e-3)
                    and not np.isclose(self.points[runner][2], z_here, 1e-3)
                ):
                    raise ValueError(
                        "The provided wavefunction from Lobster does not contain all relevant"
                        " points. "
                        "Please use a line similar to: printLCAORealSpaceWavefunction kpoint 1 "
                        "coordinates 0.0 0.0 0.0 coordinates 1.0 1.0 1.0 box bandlist 1 "
                    )

                new_x.append(x_here)
                new_y.append(y_here)
                new_z.append(z_here)

                new_real.append(self.reals[runner])
                new_imaginary.append(self.imaginaries[runner])
                new_density.append(
                    self.reals[runner] ** 2 + self.imaginaries[runner] ** 2
                )

        self.final_real = np.reshape(new_real, [Nx, Ny, Nz])
        self.final_imaginary = np.reshape(new_imaginary, [Nx, Ny, Nz])
        self.final_density = np.reshape(new_density, [Nx, Ny, Nz])

        self.volumetricdata_real = VolumetricData(structure, {"total": self.final_real})
        self.volumetricdata_imaginary = VolumetricData(
            structure, {"total": self.final_imaginary}
        )
        self.volumetricdata_density = VolumetricData(
            structure, {"total": self.final_density}
        )

    def get_volumetricdata_real(self) -> VolumetricData:
        """Get VolumetricData object for the real part of the wave function.

        Returns:
            VolumetricData: Real part volumetric data.
        """
        if not hasattr(self, "volumetricdata_real"):
            self.set_volumetric_data(self.grid, self.structure)
        return self.volumetricdata_real

    def get_volumetricdata_imaginary(self) -> VolumetricData:
        """Get VolumetricData object for the imaginary part of the wave function.

        Returns:
            VolumetricData: Imaginary part volumetric data.
        """
        if not hasattr(self, "volumetricdata_imaginary"):
            self.set_volumetric_data(self.grid, self.structure)
        return self.volumetricdata_imaginary

    def get_volumetricdata_density(self) -> VolumetricData:
        """Get VolumetricData object for the density part of the wave function.

        Returns:
            VolumetricData: Density volumetric data.
        """
        if not hasattr(self, "volumetricdata_density"):
            self.set_volumetric_data(self.grid, self.structure)
        return self.volumetricdata_density

    def write_file(
        self,
        filename: PathLike = "WAVECAR.vasp",
        part: Literal["real", "imaginary", "density"] = "real",
    ) -> None:
        """Save the wave function in a file readable by VESTA.

        Args:
            filename (PathLike): Output file name. Defaults to "WAVECAR.vasp".
            part (Literal["real", "imaginary", "density"]): Which part to save. Defaults to "real".

        Raises:
            ValueError: If the specified part is not "real", "imaginary", or "density".
        """
        if not (
            hasattr(self, "volumetricdata_real")
            and hasattr(self, "volumetricdata_imaginary")
            and hasattr(self, "volumetricdata_density")
        ):
            self.set_volumetric_data(self.grid, self.structure)

        if part == "real":
            self.volumetricdata_real.write_file(filename)
        elif part == "imaginary":
            self.volumetricdata_imaginary.write_file(filename)
        elif part == "density":
            self.volumetricdata_density.write_file(filename)
        else:
            raise ValueError('part can be only "real" or "imaginary" or "density"')

    def as_dict(self) -> dict[str, Any]:
        """Convert the `Wavefunction` object to a dictionary for serialization.

        Returns:
            dict[str, Any]: Dictionary representation of the object.
        """
        dictionary = super().as_dict()

        dictionary["kwargs"]["structure"] = self.structure.as_dict()

        dictionary["attributes"]["grid"] = self.grid
        dictionary["attributes"]["points"] = self.points
        dictionary["attributes"]["reals"] = self.reals
        dictionary["attributes"]["imaginaries"] = self.imaginaries
        dictionary["attributes"]["distances"] = self.distances

        return dictionary


class MadelungEnergies(LobsterFile):
    """Parser for MadelungEnergies.lobster files.

    Attributes:
        madelung_energies_mulliken (float): Madelung energy (Mulliken).
        madelung_energies_loewdin (float): Madelung energy (Loewdin).
        ewald_splitting (float): Ewald splitting parameter.
    """

    @version_processor()
    def parse_file(self) -> None:
        """Parse MadelungEnergies.lobster file.

        Extracts the Ewald splitting parameter and Madelung energies.

        Returns:
            None
        """
        line = self.lines[5]

        line_parts = line.split()

        self.ewald_splitting = float(line_parts[0])
        self.madelung_energies_mulliken = float(line_parts[1])
        self.madelung_energies_loewdin = float(line_parts[2])

    @classmethod
    def get_default_filename(cls) -> str:
        """Get the default filename for MadelungEnergies.

        Returns:
            str: Default filename.
        """
        return "MadelungEnergies.lobster"

    def as_dict(self) -> dict[str, Any]:
        """Convert the `MadelungEnergies` object to a dictionary for serialization.

        Returns:
            dict[str, Any]: Dictionary representation of the object.
        """
        dictionary = super().as_dict()

        dictionary["attributes"][
            "madelung_energies_mulliken"
        ] = self.madelung_energies_mulliken
        dictionary["attributes"][
            "madelung_energies_loewdin"
        ] = self.madelung_energies_loewdin
        dictionary["attributes"]["ewald_splitting"] = self.ewald_splitting

        return dictionary


class SitePotentials(LobsterFile):
    """Parser for SitePotentials.lobster files.

    Attributes:
        centers (list[str]): Atom centers.
        site_potentials_mulliken (list[float]): Mulliken site potentials.
        site_potentials_loewdin (list[float]): Loewdin site potentials.
        madelung_energies_mulliken (float): Madelung energy (Mulliken).
        madelung_energies_loewdin (float): Madelung energy (Loewdin).
        ewald_splitting (float): Ewald splitting parameter.
    """

    @version_processor()
    def parse_file(self) -> None:
        """Parse SitePotentials.lobster file.

        Extracts site potentials, Madelung energies, and Ewald splitting parameter.

        Returns:
            None
        """
        self.centers = []
        self.site_potentials_mulliken = []
        self.site_potentials_loewdin = []

        for line in self.iterate_lines():
            if ewald_splitting := re.search(r"splitting parameter\s+(\S+)", line):
                self.ewald_splitting = float(ewald_splitting.group(1))

            if madelung_energies := re.search(
                r"Madelung Energy \(eV\)\s*(\S+)\s+(\S+)", line
            ):
                self.madelung_energies_mulliken = float(madelung_energies.group(1))
                self.madelung_energies_loewdin = float(madelung_energies.group(2))

            if data := re.search(r"(\d+)\s+([a-zA-Z]{1,2})\s+(\S+)\s+(\S+)", line):
                data = data.groups()
                self.centers.append(data[1] + data[0])
                self.site_potentials_mulliken.append(float(data[2]))
                self.site_potentials_loewdin.append(float(data[3]))

    @classmethod
    def get_default_filename(cls) -> str:
        """Get the default filename for SitePotentials.

        Returns:
            str: Default filename.
        """
        return "SitePotentials.lobster"

    def as_dict(self) -> dict[str, Any]:
        """Convert the `SitePotentials` object to a dictionary for serialization.

        Returns:
            dict[str, Any]: Dictionary representation of the object.
        """
        dictionary = super().as_dict()

        dictionary["attributes"]["centers"] = self.centers
        dictionary["attributes"][
            "site_potentials_mulliken"
        ] = self.site_potentials_mulliken
        dictionary["attributes"][
            "site_potentials_loewdin"
        ] = self.site_potentials_loewdin
        dictionary["attributes"][
            "madelung_energies_mulliken"
        ] = self.madelung_energies_mulliken
        dictionary["attributes"][
            "madelung_energies_loewdin"
        ] = self.madelung_energies_loewdin
        dictionary["attributes"]["ewald_splitting"] = self.ewald_splitting

        return dictionary


def get_orb_from_str(orbs: list[str]) -> tuple[str, list[tuple[int, Orbital]]]:
    """Get Orbitals from string representations.

    Args:
        orbs (list[str]): Orbitals, e.g. ["2p_x", "3s"].

    Returns:
        tuple[str, list[tuple[int, Orbital]]]: Orbital label and list of orbitals.
    """
    orbitals = [(int(orb[0]), Orbital(LOBSTER_ORBITALS.index(orb[1:]))) for orb in orbs]

    orb_label = ""
    for iorb, orbital in enumerate(orbitals):
        if iorb == 0:
            orb_label += f"{orbital[0]}{orbital[1].name}"
        else:
            orb_label += f"-{orbital[0]}{orbital[1].name}"

    return orb_label, orbitals


class LobsterMatrices(LobsterFile):
    """Parser for LOBSTER matrix files.

    Attributes:
        matrix_type (str): Type of matrix (hamilton, coefficient, transfer, overlap).
        centers (list[str]): Atom centers.
        orbitals (list[str]): Orbitals.
        matrices (LobsterMatrixData): Matrix data for each k-point and spin.
        efermi (float): Fermi energy (for Hamilton matrices).
    """

    matrix_types: ClassVar[set[str]] = {
        "hamilton",
        "coefficient",
        "transfer",
        "overlap",
    }

    def __init__(
        self,
        filename: PathLike | None = None,
        matrix_type: str | None = None,
        efermi: float | None = None,
        process_immediately: bool = True,
    ) -> None:
        """Initialize LOBSTER matrices parser.

        Args:
            filename: Path to the matrix file
            matrix_type: Type of matrix. If None, inferred from filename
            efermi: Fermi level in eV (required for Hamilton matrices)
            process_immediately: Whether to parse the file immediately
        """
        super().__init__(filename=filename, process_immediately=False)

        self.efermi = efermi

        self.matrix_type = matrix_type or self.get_matrix_type()

        self.centers: list[str] = []
        self.orbitals: list[str] = []
        self.matrices: LobsterMatrixData = {}

        if self.matrix_type == "hamilton" and self.efermi is None:
            raise ValueError("Fermi energy (eV) required for Hamilton matrices")

        if process_immediately:
            self.parse_file()

    def get_matrix_type(self) -> str:
        """Infer matrix type from filename.

        Returns:
            str: Matrix type.
        """
        name = str(self.filename).lower()

        for matrix_type in self.matrix_types:
            if matrix_type in name:
                return matrix_type

        raise ValueError(f"Cannot infer matrix type from filename: {self.filename}")

    @version_processor()
    def parse_file(self) -> None:
        """Parse matrix data and set instance attributes.

        Returns:
            None
        """
        header_regex_pattern = (
            r"kpoint\s+(\d+)"
            if self.matrix_type == "overlap"
            else r"(\d+)\s+kpoint\s+(\d+)"
        )

        current_kpoint, current_spin = None, None
        multiplier = 1

        lines_generator = self.iterate_lines()
        for line in lines_generator:
            if header_match := re.search(header_regex_pattern, line):
                header_match = header_match.groups()
                if self.matrix_type != "overlap":
                    current_spin = Spin.up if header_match[0] == "1" else Spin.down

                current_kpoint = int(header_match[-1])
            elif "real parts" in line.lower():
                multiplier = 1
            elif "imag parts" in line.lower():
                multiplier = 1j
            elif line.startswith("basisfunction"):
                num_parts = (
                    len(re.findall(r"band\s+\d+", line))
                    if "band" in line
                    else len(line.split()[1:])
                )

                if current_kpoint not in self.matrices:
                    if current_kpoint is None:
                        raise ValueError(
                            "Could not read any k-point before matrix data."
                        )

                    self.matrices[current_kpoint] = {
                        current_spin: np.zeros((num_parts, num_parts), dtype=complex)
                    }
                elif current_spin not in self.matrices[current_kpoint]:
                    self.matrices[current_kpoint][current_spin] = np.zeros(
                        (num_parts, num_parts), dtype=complex
                    )

                values = []
                for _ in range(num_parts):
                    line_split = next(lines_generator).split()

                    values.append([float(val) * multiplier for val in line_split[1:]])

                    if (
                        len(self.centers) != num_parts
                        and len(self.orbitals) != num_parts
                    ):
                        self.centers.append(line_split[0].split("_")[0].title())
                        orbital = parse_orbital_from_text(line_split[0])

                        if orbital is None:
                            raise ValueError(
                                f"Could not read orbital format: {line_split[0]} when parsing header line: {line}"
                            )

                        self.orbitals.append(orbital)

                self.matrices[current_kpoint][current_spin] += np.array(
                    values, dtype=complex
                )

    def get_onsite_values(
        self, center: str | None = None, orbital: str | None = None
    ) -> dict | float | floating:
        """Get onsite values for specific centers/orbitals.

        Args:
            center (str | None): Specific center or None for all.
            orbital (str | None): Specific orbital or None for all.

        Returns:
            dict | float | floating: Dict of values or single value if both specified.
        """
        results = {}

        energy_shift = self.efermi if self.matrix_type == "hamilton" else 0

        for i, (c, o) in enumerate(zip(self.centers, self.orbitals, strict=True)):
            if (center is None or c == center) and (orbital is None or o == orbital):
                values = [
                    m[i, i].real - energy_shift
                    for kpoint in self.matrices.values()
                    for m in kpoint.values()
                ]
                avg_value = np.mean(values)

                if center and orbital:
                    return avg_value

                results[f"{c}_{o}"] = avg_value

        return results

    @classmethod
    def get_default_filename(cls) -> str:
        """Get the default filename for the LobsterMatrices class.

        Returns:
            str: Default filename.
        """
        return "hamiltonMatrices.lobster"

    def as_dict(self) -> dict[str, Any]:
        """Convert the `LobsterMatrices` object to a dictionary for serialization.

        Returns:
            dict[str, Any]: Dictionary representation of the object.
        """
        dictionary = super().as_dict()

        dictionary["kwargs"]["efermi"] = self.efermi

        dictionary["attributes"]["matrix_type"] = self.matrix_type
        dictionary["attributes"]["centers"] = self.centers
        dictionary["attributes"]["orbitals"] = self.orbitals
        dictionary["attributes"]["matrices"] = make_json_compatible(self.matrices)

        return dictionary


class POLARIZATION(LobsterFile):
    """Parser for POLARIZATION.lobster file.

    Attributes:
        rel_mulliken_pol_vector (dict[str, float]): Relative Mulliken polarization vector.
        rel_loewdin_pol_vector (dict[str, float]): Relative Loewdin polarization vector.
    """

    @version_processor()
    def parse_file(self) -> None:
        """Parse POLARIZATION.lobster file.

        Returns:
            None
        """
        self.rel_mulliken_pol_vector = {}
        self.rel_loewdin_pol_vector = {}

        for line in islice(self.iterate_lines(), 3, None):
            cleanlines = [idx for idx in line.split(" ") if idx != ""]
            if cleanlines and len(cleanlines) == 3:
                self.rel_mulliken_pol_vector[cleanlines[0]] = float(cleanlines[1])
                self.rel_loewdin_pol_vector[cleanlines[0]] = float(cleanlines[2])
            if cleanlines and len(cleanlines) == 4:
                self.rel_mulliken_pol_vector[cleanlines[0].replace(":", "")] = (
                    cleanlines[1].replace("\u03bc", "u")
                )
                self.rel_loewdin_pol_vector[cleanlines[2].replace(":", "")] = (
                    cleanlines[3].replace("\u03bc", "u")
                )

    @classmethod
    def get_default_filename(cls) -> str:
        """Get the default filename for the Polarization class.

        Returns:
            str: Default filename.
        """
        return "POLARIZATION.lobster"

    def as_dict(self) -> dict[str, Any]:
        """Convert the `POLARIZATION` object to a dictionary for serialization.

        Returns:
            dict[str, Any]: Dictionary representation of the object.
        """
        dictionary = super().as_dict()

        dictionary["attributes"][
            "rel_mulliken_pol_vector"
        ] = self.rel_mulliken_pol_vector
        dictionary["attributes"]["rel_loewdin_pol_vector"] = self.rel_loewdin_pol_vector

        return dictionary


class BWDF(LobsterFile):
    """Parser for BWDF.lobster/BWDFCOHP.lobster files.

    Attributes:
        centers (NDArray): Bond length centers for the distribution.
        bwdf (dict[Spin, NDArray]): Bond weighted distribution function.
        bin_width (float): Bin width used for computing the distribution by LOBSTER.
    """

    is_cohp: ClassVar[bool] = False

    def __init__(
        self,
        filename: PathLike | None = None,
        process_immediately: bool = True,
    ) -> None:
        """
        Args:
            filename (PathLike): The BWDF file from LOBSTER, typically "BWDF.lobster"
                or "BWDFCOHP.lobster".
        """
        self.bwdf = {}
        self.centers = np.array([])
        self.data = np.array([[]])

        super().__init__(filename=filename, process_immediately=process_immediately)

    @version_processor()
    def parse_file(self) -> None:
        """Parse BWDF.lobster/BWDFCOHP.lobster file.

        Returns:
            None
        """
        self.bwdf = {}
        self.data = np.genfromtxt(self.lines[1:], dtype=float)

        self.process_data_into_bwdf_centers()

    def process_data_into_bwdf_centers(self) -> None:
        """Process data into bwdf and centers.

        Returns:
            None
        """
        self.centers = self.data[:, 0]
        self.bwdf[Spin.up] = self.data[:, 1]

        if self.data.shape[1] > 2:
            self.bwdf[Spin.down] = self.data[:, 2]

    @classmethod
    def get_default_filename(cls) -> str:
        """Get the default filename for the BWDF class.

        Returns:
            str: Default filename.
        """
        return "BWDF.lobster"

    def as_dict(self) -> dict[str, Any]:
        """Convert the `BWDF` object to a dictionary for serialization.

        Returns:
            dict[str, Any]: Dictionary representation of the object.
        """
        dictionary = super().as_dict()

        dictionary["attributes"]["data"] = self.data

        return dictionary

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> BWDF:
        """Create a `BWDF` object from a dictionary.

        Args:
            d (dict): Dictionary representation of a `BWDF` object.

        Returns:
            BWDF: The created object.
        """
        instance = super().from_dict(d)
        instance.process_data_into_bwdf_centers()

        return instance


class BWDFCOHP(BWDF):
    """Parser for BWDFCOHP.lobster files.

    Returns:
        None
    """

    @classmethod
    def get_default_filename(cls) -> str:
        """Get the default filename for the BWDFCOHP class.

        Returns:
            str: Default filename.
        """
        return "BWDFCOHP.lobster"
