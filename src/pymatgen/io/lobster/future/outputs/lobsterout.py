from __future__ import annotations

import re
from typing import TYPE_CHECKING
from warnings import warn

from pymatgen.io.lobster.future.constants import LOBSTER_VERSION
from pymatgen.io.lobster.future.core import LobsterFile
from pymatgen.io.lobster.future.versioning import version_processor

if TYPE_CHECKING:
    from typing import Literal


class LobsterOut(LobsterFile):
    """Parser for `lobsterout` file from LOBSTER.

    This class reads the `lobsterout` file and extracts information about
    basis functions, spillings, warnings, timing, and file presence. It supports
    parsing for different LOBSTER versions and provides attributes to access
    the parsed data.

    Attributes:
        basis_functions (list[str]): Basis functions used in the LOBSTER run.
        basis_type (list[str]): Basis types used in the LOBSTER run.
        charge_spilling (list[float]): Charge spilling for each spin channel.
        dft_program (str): DFT program used for the calculation.
        elements (list[str]): Elements present in the calculation.
        has_charge (bool): Whether `CHARGE.lobster` is present.
        has_cohpcar (bool): Whether `COHPCAR.lobster` and `ICOHPLIST.lobster` are present.
        has_madelung (bool): Whether `SitePotentials.lobster` and `MadelungEnergies.lobster` are present.
        has_coopcar (bool): Whether `COOPCAR.lobster` and `ICOOPLIST.lobster` are present.
        has_cobicar (bool): Whether `COBICAR.lobster` and `ICOBILIST.lobster` are present.
        has_doscar (bool): Whether `DOSCAR.lobster` is present.
        has_doscar_lso (bool): Whether `DOSCAR.LSO.lobster` is present.
        has_projection (bool): Whether `projectionData.lobster` is present.
        has_bandoverlaps (bool): Whether `bandOverlaps.lobster` is present.
        has_density_of_energies (bool): Whether `DensityOfEnergy.lobster` is present.
        has_fatbands (bool): Whether fatband calculation was performed.
        has_grosspopulation (bool): Whether `GROSSPOP.lobster` is present.
        has_polarization (bool): Whether `POLARIZATION.lobster` is present.
        info_lines (list[str]): Additional information on the run.
        info_orthonormalization (list[str]): Information on orthonormalization.
        is_restart_from_projection (bool): Whether calculation was restarted from a projection file.
        lobster_version (str): The LOBSTER version.
        number_of_spins (int): Number of spins.
        number_of_threads (int): Number of threads used.
        timing (dict[str, float]): Timing information.
        total_spilling (list[float]): Total spilling for each spin channel.
        warning_lines (list[str]): All warning messages.
    """

    @version_processor(min_version="5.1")
    def _process_v5_1(self) -> None:
        """Process `lobsterout` (version ≥ 5.1).

        This method extracts file presence flags and other attributes
        specific to LOBSTER version 5.1 and later.

        Note:
            This method is automatically invoked for files with version ≥5.1.
        """
        lines = self.lines

        self.has_cohpcar = (
            "writing COOPCAR.lobster..." in lines
            and "SKIPPING writing COOPCAR.lobster..." not in lines
        )
        self.has_coopcar = (
            "writing COHPCAR.lobster..." in lines
            and "SKIPPING writing COHPCAR.lobster..." not in lines
        )
        self.has_cobicar = (
            "writing COBICAR.lobster..." in lines
            and "SKIPPING writing COBICAR.lobster..." not in lines
        )

        self._process_common()

    def _process_common(self) -> None:
        """Process common parts of `lobsterout` for all versions.

        This method extracts general information such as timing, warnings,
        basis functions, and file presence flags that are applicable across
        all LOBSTER versions.
        """
        lines = self.lines

        self.is_restart_from_projection = (
            "loading projection from projectionData.lobster..." in lines
        )

        self.has_error = "ERROR:" in lines

        if self.has_error:
            self.error_lines = [line for line in lines if line.startswith("ERROR:")]
            raise RuntimeError(
                f"LOBSTER calculation ended with errors:\n{self.error_lines}"
            )

        self.number_of_threads = self._get_threads(lines)
        self.dft_program = self._get_dft_program(lines)

        self.number_of_spins = self._get_number_of_spins(lines)
        self.charge_spilling, self.total_spilling = self._get_spillings(
            lines=lines, number_of_spins=self.number_of_spins
        )

        self.elements, self.basis_type, self.basis_functions = (
            self._get_elements_basistype_basisfunctions(lines)
        )

        wall_time, user_time, sys_time = self._get_timing(lines)
        self.timing = {
            "wall_time": wall_time,
            "user_time": user_time,
            "sys_time": sys_time,
        }

        self.warning_lines = self._get_all_warning_lines(lines)

        self.info_orthonormalization = self._get_warning_orthonormalization(lines)

        self.info_lines = self._get_all_info_lines(lines)

        self.has_doscar = (
            "writing DOSCAR.lobster..." in lines
            and "SKIPPING writing DOSCAR.lobster..." not in lines
        )
        self.has_doscar_lso = (
            "writing DOSCAR.LSO.lobster..." in lines
            and "SKIPPING writing DOSCAR.LSO.lobster..." not in lines
        )

        self.has_cobicar_lcfo = "writing COBICAR.LCFO.lobster..." in lines
        self.has_cohpcar_lcfo = "writing COHPCAR.LCFO.lobster..." in lines
        self.has_coopcar_lcfo = "writing COOPCAR.LCFO.lobster..." in lines
        self.has_doscar_lcfo = "writing DOSCAR.LCFO.lobster..." in lines
        self.has_polarization = (
            "writing polarization to POLARIZATION.lobster..." in lines
        )
        self.has_charge = "SKIPPING writing CHARGE.lobster..." not in lines
        self.has_projection = "saving projection to projectionData.lobster..." in lines
        self.has_bandoverlaps = (
            "WARNING: I dumped the band overlap matrices to the file bandOverlaps.lobster."
            in lines
        )
        self.has_fatbands = self._has_fatband(lines)
        self.has_grosspopulation = (
            "writing CHARGE.lobster and GROSSPOP.lobster..." in lines
        )
        self.has_density_of_energies = "writing DensityOfEnergy.lobster..." in lines
        self.has_madelung = (
            "writing SitePotentials.lobster and MadelungEnergies.lobster..." in lines
            and "skipping writing SitePotentials.lobster and MadelungEnergies.lobster..."
            not in lines
        )
        self.has_mofecar = "Writing MOFECAR.lobster and IMOFELIST.lobster..." in lines

    @version_processor(max_version="5.0")
    def _process_legacy(self) -> None:
        """Process `lobsterout` for legacy versions (≤5.0).

        This method extracts file presence flags and other attributes
        specific to LOBSTER versions ≤5.0.

        Note:
            This method is automatically invoked for files with version ≤5.0.
        """
        lines = self.lines

        self.has_cohpcar = (
            "writing COOPCAR.lobster and ICOOPLIST.lobster..." in lines
            and "SKIPPING writing COOPCAR.lobster and ICOOPLIST.lobster..." not in lines
        )
        self.has_coopcar = (
            "writing COHPCAR.lobster and ICOHPLIST.lobster..." in lines
            and "SKIPPING writing COHPCAR.lobster and ICOHPLIST.lobster..." not in lines
        )
        self.has_cobicar = (
            "writing COBICAR.lobster and ICOBILIST.lobster..." in lines
            and "SKIPPING writing COBICAR.lobster and ICOBILIST.lobster..." not in lines
        )

        self._process_common()

    def process(self) -> None:
        """Parse the `lobsterout` file and populate attributes.

        This method determines the LOBSTER version and invokes the appropriate
        version-specific processing method.

        Raises:
            RuntimeError: If the LOBSTER version cannot be determined.
        """
        self.lobster_version = self.get_lobster_version(self.lines)

        super().process()

    @staticmethod
    def get_lobster_version(lines: list[str]) -> str:
        """Get the LOBSTER version from the `lobsterout` lines.

        Args:
            lines (list[str]): Lines of the `lobsterout` file.

        Returns:
            str: The LOBSTER version.

        Raises:
            RuntimeError: If the version line is not found.
        """
        for line in lines:
            if version := re.search(r"(?i)LOBSTER\s(?:v(\d+\.\d+\.\d+))", line):
                return version.group(1)

        warn(
            f"Could not find LOBSTER version in lobsterout. Defaulting to v{LOBSTER_VERSION}",
            stacklevel=2,
        )

        return LOBSTER_VERSION

    @staticmethod
    def _has_fatband(lines: list[str]) -> bool:
        """Check whether the calculation includes fatband data.

        Args:
            lines (list[str]): Lines of the `lobsterout` file.

        Returns:
            bool: True if fatband data is present, False otherwise.
        """
        for line in lines:
            line_parts = line.split()
            if len(line_parts) > 1 and line_parts[1] == "FatBand":
                return True
        return False

    @staticmethod
    def _get_dft_program(lines: list[str]) -> str | None:
        """Get the DFT program used for the calculation.

        Args:
            lines (list[str]): Lines of the `lobsterout` file.

        Returns:
            str | None: The name of the DFT program, or None if not found.
        """
        for line in lines:
            line_parts = line.split()
            if len(line_parts) > 4 and line_parts[3] == "program...":
                return line_parts[4]
        return None

    @staticmethod
    def _get_number_of_spins(lines: list[str]) -> Literal[1, 2]:
        """Get the number of spin channels.

        Args:
            lines (list[str]): Lines of the `lobsterout` file.

        Returns:
            Literal[1, 2]: 1 if single spin channel, 2 if two spin channels.
        """
        return 2 if "spillings for spin channel 2" in lines else 1

    @staticmethod
    def _get_threads(lines: list[str]) -> int:
        """Get the number of CPU threads used.

        Args:
            lines (list[str]): Lines of the `lobsterout` file.

        Returns:
            int: Number of threads.

        Raises:
            ValueError: If the number of threads cannot be determined.
        """
        for line in lines:
            line_parts = line.split()
            if len(line_parts) > 11 and line_parts[11] in {"threads", "thread"}:
                return int(line_parts[10])
        raise ValueError("Threads not found.")

    @staticmethod
    def _get_spillings(
        lines: list[str],
        number_of_spins: Literal[1, 2],
    ) -> tuple[list[float], list[float]]:
        """Get charge spillings and total spillings.

        Args:
            lines (list[str]): Lines of the `lobsterout` file.
            number_of_spins (Literal[1, 2]): Number of spin channels.

        Returns:
            tuple[list[float], list[float]]: Charge spillings and total spillings
            for each spin channel.
        """
        charge_spillings = []
        total_spillings = []
        for line in lines:
            line_parts = line.split()
            if len(line_parts) > 2 and line_parts[2] == "spilling:":
                if line_parts[1] == "charge":
                    charge_spillings.append(
                        float(line_parts[3].replace("%", "")) / 100.0
                    )
                elif line_parts[1] == "total":
                    total_spillings.append(
                        float(line_parts[3].replace("%", "")) / 100.0
                    )

            if (
                len(charge_spillings) == number_of_spins
                and len(total_spillings) == number_of_spins
            ):
                break

        return charge_spillings, total_spillings

    @staticmethod
    def _get_elements_basistype_basisfunctions(
        lines: list[str],
    ) -> tuple[list[str], list[str], list[list[str]]]:
        """Get elements, basis types, and basis functions.

        Args:
            lines (list[str]): Lines of the `lobsterout` file.

        Returns:
            tuple[list[str], list[str], list[list[str]]]: Elements, basis types,
            and basis functions used in the calculation.
        """
        begin = False
        end = False
        elements: list[str] = []
        basistypes: list[str] = []
        basisfunctions: list[list[str]] = []
        for line in lines:
            if begin and not end:
                line_parts = line.split()
                if line_parts[0] not in {
                    "INFO:",
                    "WARNING:",
                    "setting",
                    "calculating",
                    "post-processing",
                    "saving",
                    "spillings",
                    "writing",
                }:
                    elements.append(line_parts[0])
                    basistypes.append(line_parts[1].replace("(", "").replace(")", ""))
                    # Last sign is ''
                    basisfunctions.append(line_parts[2:])
                else:
                    end = True

            if "setting up local basis functions..." in line:
                begin = True
        return elements, basistypes, basisfunctions

    @staticmethod
    def _get_timing(
        lines: list[str],
    ) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
        """Get wall time, user time, and system time.

        Args:
            lines (list[str]): Lines of the `lobsterout` file.

        Returns:
            tuple[dict[str, str], dict[str, str], dict[str, str]]: Dictionaries
            containing timing information for wall, user, and system times.
        """
        begin = False
        user_times, wall_times, sys_times = [], [], []

        for line in lines:
            line_parts = line.split()
            if "finished" in line_parts:
                begin = True
            if begin:
                if "wall" in line_parts:
                    wall_times = line_parts[2:10]
                if "user" in line_parts:
                    user_times = line_parts[:8]
                if "sys" in line_parts:
                    sys_times = line_parts[:8]

        wall_time_dict = {
            "h": wall_times[0],
            "min": wall_times[2],
            "s": wall_times[4],
            "ms": wall_times[6],
        }
        user_time_dict = {
            "h": user_times[0],
            "min": user_times[2],
            "s": user_times[4],
            "ms": user_times[6],
        }
        sys_time_dict = {
            "h": sys_times[0],
            "min": sys_times[2],
            "s": sys_times[4],
            "ms": sys_times[6],
        }

        return wall_time_dict, user_time_dict, sys_time_dict

    @staticmethod
    def _get_warning_orthonormalization(lines: list[str]) -> list[str]:
        """Get orthonormalization warnings.

        Args:
            lines (list[str]): Lines of the `lobsterout` file.

        Returns:
            list[str]: List of orthonormalization warnings.
        """
        orthowarnings = []
        for line in lines:
            line_parts = line.split()
            if "orthonormalized" in line_parts:
                orthowarnings.append(" ".join(line_parts[1:]))
        return orthowarnings

    @staticmethod
    def _get_all_warning_lines(lines: list[str]) -> list[str]:
        """Get all warning lines.

        Args:
            lines (list[str]): Lines of the `lobsterout` file.

        Returns:
            list[str]: List of warning messages.
        """
        warnings_ = []
        for line in lines:
            line_parts = line.split()
            if len(line_parts) > 0 and line_parts[0] == "WARNING:":
                warnings_.append(" ".join(line_parts[1:]))
        return warnings_

    @staticmethod
    def _get_all_info_lines(lines: list[str]) -> list[str]:
        """Get all informational lines.

        Args:
            lines (list[str]): Lines of the `lobsterout` file.

        Returns:
            list[str]: List of informational messages.
        """
        infos = []
        for line in lines:
            line_parts = line.split()
            if len(line_parts) > 0 and line_parts[0] == "INFO:":
                infos.append(" ".join(line_parts[1:]))
        return infos

    @classmethod
    def get_default_filename(cls) -> str:
        """Get the default filename for `LobsterOut`.

        Returns:
            str: The default filename.
        """
        return "lobsterout"
