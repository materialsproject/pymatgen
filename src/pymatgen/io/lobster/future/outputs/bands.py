from __future__ import annotations

import re
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, Self, cast

import numpy as np
from monty.json import MSONable

from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.lobster.future.constants import LOBSTER_VERSION
from pymatgen.io.lobster.future.core import LobsterFile
from pymatgen.io.lobster.future.utils import natural_sort, parse_orbital_from_text
from pymatgen.io.lobster.future.versioning import version_processor
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.outputs import Vasprun

if TYPE_CHECKING:
    from pymatgen.core.structure import IStructure
    from pymatgen.io.lobster.future.types import LobsterBandOverlaps, LobsterFatband
    from pymatgen.util.typing import PathLike


class BandOverlaps(LobsterFile):
    """Parser for bandOverlaps.lobster files.

    Parses band overlap information produced by LOBSTER and stores it in a structured dictionary keyed by spin.
    See the :class:`~pymatgen.io.lobster.future.types.LobsterBandOverlaps` type for details.

    Attributes:
        band_overlaps (dict[Spin, dict]):
            "k_points", "max_deviations", and "matrices" holding the corresponding data.
            - "k_points" (list[list[float]]): List of k-point coordinates.
            - "max_deviations" (list[float]): List of maximal deviations for each k-point.
            - "matrices" (list[np.ndarray]): List of overlap matrices for each k-point.

            each holding data for each spin channel.
    """

    @version_processor(max_version="3.2")
    def parse_file_v3_2_legacy(self) -> None:
        """Parse bandOverlaps.lobster file for LOBSTER versions ≤3.2.

        Uses legacy spin numbering [0, 1] for parsing.
        """
        self.parse_file(spin_numbers=[0, 1])

    @version_processor(min_version="4.0")
    def parse_file_v4_0(self) -> None:
        """Parse bandOverlaps.lobster file for LOBSTER versions ≥4.0.

        Uses updated spin numbering [1, 2] for parsing.
        """
        self.parse_file(spin_numbers=[1, 2])

    def parse_file(self, spin_numbers: list[int]) -> None:
        """Read all lines of the file and populate `self.band_overlaps`.

        Args:
            spin_numbers (list[int]): Two integers indicating the spin numbering used
                in the file (e.g., [0, 1] for legacy or [1, 2] for newer versions).

        Raises:
            ValueError: If no data is found for a key in the bandOverlaps file.
        """
        n_kpoints = {Spin.up: 0, Spin.down: 0}
        matrix_size = None

        current_spin = Spin.up
        for line in self.iterate_lines():
            if f"Overlap Matrix (abs) of the orthonormalized projected bands for spin {spin_numbers[0]}" in line:
                current_spin = Spin.up
            elif f"Overlap Matrix (abs) of the orthonormalized projected bands for spin {spin_numbers[1]}" in line:
                current_spin = Spin.down
            elif "maxDeviation" in line:
                n_kpoints[current_spin] += 1
            elif matrix_size is None:
                try:
                    float(line.split()[0])
                    matrix_size = len(line.split())
                except ValueError:
                    continue

        if matrix_size is None:
            raise ValueError("No data found for band overlaps in the file.")

        self.band_overlaps: LobsterBandOverlaps = {
            "k_points": {},
            "max_deviations": {},
            "matrices": {},
        }

        self.spins = [Spin.up]
        if n_kpoints[Spin.down] > 0:
            self.spins.append(Spin.down)

        for spin in self.spins:
            n = n_kpoints[spin]
            self.band_overlaps["k_points"][spin] = np.empty((n, 3), dtype=np.float64)
            self.band_overlaps["max_deviations"][spin] = np.empty(n, dtype=np.float64)
            self.band_overlaps["matrices"][spin] = np.empty((n, matrix_size, matrix_size), dtype=np.float64)

        current_spin = Spin.up
        kpoint_idx = {Spin.up: 0, Spin.down: 0}
        matrix_row = 0
        for line in self.iterate_lines():
            if f"Overlap Matrix (abs) of the orthonormalized projected bands for spin {spin_numbers[0]}" in line:
                current_spin = Spin.up
            elif f"Overlap Matrix (abs) of the orthonormalized projected bands for spin {spin_numbers[1]}" in line:
                current_spin = Spin.down
            elif "k-point" in line:
                self.band_overlaps["k_points"][current_spin][kpoint_idx[current_spin]] = [
                    float(el) for el in line.strip().split()[-3:]
                ]
                kpoint_idx[current_spin] += 1
            elif "maxDeviation" in line:
                maxdev = line.split(" ")[-1]
                self.band_overlaps["max_deviations"][current_spin][kpoint_idx[current_spin] - 1] = float(maxdev)

                matrix_row = 0
            elif line.strip():
                try:
                    parts = [float(el) for el in re.split(r"\s+", line.strip())]
                except ValueError:
                    raise ValueError(f"Incomplete or non-numeric data found in bandOverlaps file at line: {line}")

                if len(parts) == matrix_size:
                    self.band_overlaps["matrices"][current_spin][kpoint_idx[current_spin] - 1, matrix_row] = parts

                matrix_row += 1

    def has_good_quality_max_deviation(self, limit_max_deviation: float = 0.1) -> bool:
        """Check if the maxDeviation values are within a given limit.

        Args:
            limit_max_deviation (float): Upper limit for acceptable max_deviation.

        Returns:
            bool: True if all recorded max_deviation values are <= limit_max_deviation.
        """
        return all(
            deviation <= limit_max_deviation for deviation in self.band_overlaps["max_deviations"].get(Spin.up, [])
        ) and all(
            deviation <= limit_max_deviation for deviation in self.band_overlaps["max_deviations"].get(Spin.down, [])
        )

    def has_good_quality_check_occupied_bands(
        self,
        number_occ_bands_spin_up: int,
        number_occ_bands_spin_down: int | None = None,
        spin_polarized: bool = False,
        limit_deviation: float = 0.1,
    ) -> bool:
        """Check if deviation from the ideal overlap for occupied bands is acceptable.

        Args:
            number_occ_bands_spin_up (int): Number of occupied bands for spin up.
            number_occ_bands_spin_down (int | None): Number of occupied bands for spin down.
                Required if spin_polarized is True.
            spin_polarized (bool): Whether the calculation is spin-polarized.
            limit_deviation (float): Acceptable absolute tolerance for deviations.

        Raises:
            ValueError: If `number_occ_bands_spin_down` is not specified for spin-polarized calculations.

        Returns:
            bool: True if all occupied-band submatrices are close to identity within the tolerance.
        """
        if spin_polarized and number_occ_bands_spin_down is None:
            raise ValueError("number_occ_bands_spin_down has to be specified")

        for spin in (Spin.up, Spin.down) if spin_polarized else (Spin.up,):
            if spin is Spin.up:
                num_occ_bands = number_occ_bands_spin_up
            else:
                if number_occ_bands_spin_down is None:
                    raise ValueError("number_occ_bands_spin_down has to be specified")
                num_occ_bands = number_occ_bands_spin_down

            for overlap_matrix in self.band_overlaps["matrices"][spin]:
                sub_array = np.asarray(overlap_matrix)[:num_occ_bands, :num_occ_bands]

                if not np.allclose(sub_array, np.identity(num_occ_bands), atol=limit_deviation, rtol=0):
                    return False

        return True

    @classmethod
    def get_default_filename(cls) -> str:
        """Get the default filename for the `BandOverlaps` class.

        Returns:
            str: Default filename.
        """
        return "bandOverlaps.lobster"

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> Self:
        """Reconstruct a BandOverlaps instance from a dictionary.

        Args:
            d (dict[str, Any]): Dictionary representation of a BandOverlaps instance.

        Returns:
            BandOverlaps: Reconstructed instance.
        """
        instance = super().from_dict(d)

        instance.band_overlaps["k_points"] = {
            spin: np.asarray(k_points, dtype=np.float64)
            for spin, k_points in instance.band_overlaps["k_points"].items()
        }
        instance.band_overlaps["max_deviations"] = {
            spin: np.asarray(deviations, dtype=np.float64)
            for spin, deviations in instance.band_overlaps["max_deviations"].items()
        }
        instance.band_overlaps["matrices"] = {
            spin: np.asarray(matrices, dtype=np.float64)
            for spin, matrices in instance.band_overlaps["matrices"].items()
        }

        return instance


class Fatbands(MSONable):
    """Reader for multiple FATBAND_*.lobster files in a directory.

    Collects FATBAND files, reads VASP outputs for the Fermi level and kpoints, and aggregates per-file parsed data.

    Attributes:
        efermi (float): Fermi level read from vasprun.xml.
        spins (list[Spin]): Spins present in the FATBAND files.
        kpoints (Kpoints): Parsed KPOINTS used for the lobster FatBand calculations.
        filenames (list[Path]): Sorted list of matched filenames.
        structure (Structure): Structure object used for projections.
        reciprocal_lattice (Lattice): Reciprocal lattice of the structure.
        lobster_version (str): LOBSTER version string used for parsing.
        fatbands (list[dict]): Aggregated parsed fatband data after process() is called.
    """

    def __init__(
        self,
        directory: PathLike = ".",
        structure: IStructure | None = None,
        kpoints_file: PathLike = "KPOINTS",
        vasprun_file: PathLike = "vasprun.xml",
        process_immediately: bool = True,
        lobster_version: str | None = None,
    ) -> None:
        """Initialize the Fatbands reader.

        Args:
            directory (PathLike): Path to directory containing FATBAND files.
            structure (IStructure | None): Structure object. If None, POSCAR.lobster is read from directory.
            kpoints_file (PathLike): Name of the KPOINTS file to be read in the directory.
            vasprun_file (PathLike): Name of the vasprun.xml file to be read in the directory.
            process_immediately (bool): If True, process FATBAND files immediately after initialization.
            lobster_version (str | None): Optional LOBSTER version string. If None, the default LOBSTER_VERSION is used.

        Raises:
            FileNotFoundError: If required files are missing in the directory.
            ValueError: If no FATBAND files are found or KPOINTS file lacks weights.
        """
        self.directory = Path(directory)

        self.filenames = sorted(
            [str(i) for i in self.directory.glob(Fatband.get_default_filename())],
            key=natural_sort,
        )

        if len(self.filenames) == 0:
            raise ValueError("No FATBAND files found in the provided directory")

        self.efermi = Vasprun(
            filename=self.directory / vasprun_file,
            parse_eigen=False,
            parse_potcar_file=False,
        ).efermi
        self.spins = [Spin.up]
        full_kpoints = Kpoints.from_file(self.directory / kpoints_file)

        if full_kpoints.kpts_weights is not None:
            filtered_data = [
                (k, w, n)
                for k, w, n in zip(
                    full_kpoints.kpts,
                    full_kpoints.kpts_weights,
                    (full_kpoints.labels or [None] * len(full_kpoints.kpts)),
                    strict=True,
                )
                if w == 0
            ]

            new_kpts, new_weights, new_labels = zip(*filtered_data, strict=True) if filtered_data else ([], [], [])

            coord_type = full_kpoints.coord_type

            if coord_type is None:
                pass
            elif coord_type not in {"Reciprocal", "Cartesian"}:
                raise ValueError("KPOINTS coord_type must be 'Reciprocal' or 'Cartesian' for `Fatbands` parsing.")

            coord_type = cast("Literal['Reciprocal', 'Cartesian'] | None", coord_type)

            self.kpoints = Kpoints(
                comment=full_kpoints.comment,
                num_kpts=len(new_kpts),
                style=full_kpoints.style,
                kpts=list(new_kpts),
                kpts_weights=list(new_weights),
                labels=list(new_labels) if any(new_labels) else None,
                coord_type=coord_type,
            )
        else:
            raise ValueError("KPOINTS file must contain weights for `Fatbands` parsing.")

        if structure is None:
            try:
                self.structure = Structure.from_file(Path(directory, "POSCAR.lobster"))
            except FileNotFoundError:
                raise FileNotFoundError("No POSCAR.lobster file found in directory, structure has to be given")
        else:
            self.structure = structure

        self.reciprocal_lattice = self.structure.lattice.reciprocal_lattice

        self.fatbands: list[LobsterFatband] = []

        self.lobster_version = lobster_version or LOBSTER_VERSION

        if process_immediately:
            self.process()

    def process(self) -> None:
        """Parse all FATBAND files and aggregate fatband data.

        Raises:
            ValueError: If the number of kpoints does not match or if there is a mix of spin-polarized and
            non-spin-polarized files.
        """
        is_spin_polarized = None

        for filename in self.filenames:
            fatband = Fatband(
                filename=filename,
                process_immediately=False,
            )
            fatband.lobster_version = self.lobster_version

            fatband.process()
            fatband_data = fatband.fatband

            for spin in fatband.spins:
                if len(fatband_data["projections"][spin]) != self.kpoints.num_kpts:
                    raise ValueError(
                        f"Number of kpoints ({self.kpoints.num_kpts}) does not "
                        f"match number of kpoints for {filename} "
                        f"({len(fatband_data['projections'][spin])})"
                    )

            if is_spin_polarized is None:
                is_spin_polarized = len(fatband_data["projections"]) > 1
            elif is_spin_polarized != (len(fatband_data["projections"]) > 1):
                raise ValueError("Mix of spin polarized and non-spin polarized FATBAND files")

            self.fatbands.append(fatband_data)

        if is_spin_polarized:
            self.spins.append(Spin.down)

    as_dict = LobsterFile.as_dict

    has_spin = LobsterFile.has_spin

    is_spin_polarized = LobsterFile.is_spin_polarized


class Fatband(LobsterFile):
    """Parser for a single FATBAND_*.lobster file.

    Parses a single FATBAND file and stores:
        center (str): Central atom/species label parsed from filename.
        orbital (str): Orbital descriptor parsed from filename.
        nbands (int): Number of bands in the FATBAND file.
        fatband (LobsterFatband): Parsed fatband data dictionary. Please see
        :class:`~pymatgen.io.lobster.future.types.LobsterFatband` for details.

    The parsed data is available in the fatband attribute after parse_file().
    """

    def __init__(
        self,
        filename: PathLike,
        process_immediately: bool = True,
        lobster_version: str | None = None,
    ) -> None:
        """Initialize a Fatband parser.

        Args:
            filename (PathLike): Path to the FATBAND file to parse.
            process_immediately (bool): If True, parse the file during initialization.
            lobster_version (str | None): LOBSTER version string to use for parsing. If None, attempts to detect
            from file or falls back to default.

        Raises:
            ValueError: If the orbital name cannot be parsed from the filename.
        """
        self.center = Path(filename).name.split("_")[1].title()

        if orbital := parse_orbital_from_text(Path(filename).stem):
            self.orbital = orbital
        else:
            raise ValueError(
                f"Could not parse orbital from filename {filename}. "
                "Ensure it follows the FATBAND_<center>_<orbital>.lobster pattern."
            )

        super().__init__(
            filename=filename,
            process_immediately=process_immediately,
            lobster_version=lobster_version,
        )

    @version_processor()
    def parse_file(self) -> None:
        """Parse the FATBAND file and populate the fatband attribute."""
        fatband: dict[str, dict[Spin, list[Any]]] = {
            "energies": {Spin.up: []},
            "projections": {Spin.up: []},
        }
        self.spins = [Spin.up]

        current_spin = Spin.up
        for idx, line in enumerate(self.iterate_lines()):
            if idx == 0:
                self.nbands = int(line.split()[-1])
                continue

            if line.startswith("#"):
                current_spin = Spin.up

                fatband["energies"][current_spin].append([])
                fatband["projections"][current_spin].append([])

                continue

            if len(fatband["projections"][current_spin][-1]) == self.nbands:
                current_spin = Spin.down

                if current_spin not in self.spins:
                    self.spins.append(current_spin)
                    fatband["energies"][current_spin] = [[]]
                    fatband["projections"][current_spin] = [[]]
                else:
                    fatband["energies"][current_spin].append([])
                    fatband["projections"][current_spin].append([])

            if data := re.findall(r"[+-]?(?:[0-9]*[.])?[0-9]+", line):
                fatband["energies"][current_spin][-1].append(float(data[1]))
                fatband["projections"][current_spin][-1].append(float(data[-1]))

        self.fatband: LobsterFatband = {
            "center": self.center,
            "orbital": self.orbital,
            "energies": fatband["energies"],
            "projections": fatband["projections"],
        }

        self.convert_to_numpy_arrays()

    def convert_to_numpy_arrays(self) -> None:
        """Convert lists in band_overlaps to numpy arrays."""
        for spin in self.spins:
            self.fatband["energies"][spin] = np.asarray(self.fatband["energies"][spin], dtype=np.float64)
            self.fatband["projections"][spin] = np.asarray(self.fatband["projections"][spin], dtype=np.float64)

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> Self:
        """Reconstruct a Fatband instance from a dictionary.

        Args:
            d (dict[str, Any]): Dictionary representation of a Fatband instance.

        Returns:
            Fatband: Reconstructed instance.
        """
        instance = super().from_dict(d)
        instance.convert_to_numpy_arrays()

        return instance

    @classmethod
    def get_default_filename(cls) -> str:
        """Get the default filename for the Fatband class.

        Returns:
            str: Default filename pattern.
        """
        return "FATBAND_*.lobster"
