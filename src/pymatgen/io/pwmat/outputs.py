from __future__ import annotations

import copy
import linecache
from io import StringIO
from typing import TYPE_CHECKING

import numpy as np
from monty.io import zopen
from monty.json import MSONable

from pymatgen.io.pwmat.inputs import ACstrExtractor, AtomConfig, LineLocator

if TYPE_CHECKING:
    from pymatgen.core import Structure
    from pymatgen.util.typing import PathLike, Tuple3Ints


__author__ = "Hanyu Liu"
__email__ = "domainofbuaa@gmail.com"
__date__ = "2024-1-16"


class Movement(MSONable):
    """Parser for data in MOVEMENT which records trajectory during MD."""

    def __init__(
        self,
        filename: PathLike,
        ionic_step_skip: int | None = None,
        ionic_step_offset: int | None = None,
    ):
        """Initialization function.

        Args:
            filename (PathLike): The path of MOVEMENT
            ionic_step_skip (int | None, optional): If ionic_step_skip is a number > 1,
                only every ionic_step_skip ionic steps will be read for
                structure and energies. This is very useful if you are parsing
                very large MOVEMENT files. Defaults to None.
            ionic_step_offset (int | None, optional): Used together with ionic_step_skip.
                If set, the first ionic step read will be offset by the amount of
                ionic_step_offset. Defaults to None.
        """
        self.filename: PathLike = filename
        self.ionic_step_skip: int | None = ionic_step_skip
        self.ionic_step_offset: int | None = ionic_step_offset
        self.split_mark: str = "--------------------------------------"

        self.chunk_sizes, self.chunk_starts = self._get_chunk_info()
        self.n_ionic_steps: int = len(self.chunk_sizes)

        self.ionic_steps: list[dict] = self._parse_sefv()
        if self.ionic_step_offset and self.ionic_step_skip:
            self.ionic_steps = self.ionic_steps[self.ionic_step_offset :: self.ionic_step_skip]

    def _get_chunk_info(self) -> tuple[list[int], list[int]]:
        """Split MOVEMENT into many chunks, so that program process it chunk by chunk.

        Returns:
            tuple[list[int], list[int]]:
                chunk_sizes (list[int]): The number of lines occupied by structural
                    information in each step.
                chunk_starts (list[int]): The starting line number for structural
                    information in each step.
        """
        chunk_sizes: list[int] = []
        row_idxs: list[int] = LineLocator.locate_all_lines(self.filename, self.split_mark)
        chunk_sizes.append(row_idxs[0])
        for ii in range(1, len(row_idxs)):
            chunk_sizes.append(row_idxs[ii] - row_idxs[ii - 1])
        chunk_sizes_bak: list[int] = copy.deepcopy(chunk_sizes)
        chunk_sizes_bak.insert(0, 0)
        chunk_starts: list[int] = np.cumsum(chunk_sizes_bak).tolist()
        chunk_starts.pop(-1)
        return chunk_sizes, chunk_starts

    @property
    def atom_configs(self) -> list[Structure]:
        """AtomConfig for structures contained in MOVEMENT file.

        Returns:
            list[Structure]: List of Structure objects for the structure at each ionic step.
        """
        return [step["atom_config"] for step in self.ionic_steps]

    @property
    def e_tots(self) -> np.ndarray:
        """Total energies of each ionic step structures contained in MOVEMENT.

        Returns:
            np.ndarray: Total energy of of each ionic step structure,
                with shape of (n_ionic_steps,).
        """
        return np.array([step["e_tot"] for step in self.ionic_steps])

    @property
    def atom_forces(self) -> np.ndarray:
        """Forces on atoms in each structures contained in MOVEMENT.

        Returns:
            np.ndarray: The forces on atoms of each ionic step structure,
                with shape of (n_ionic_steps, n_atoms, 3).
        """
        return np.array([step["atom_forces"] for step in self.ionic_steps])

    @property
    def e_atoms(self) -> np.ndarray:
        """Individual energies of atoms in each ionic step structures
        contained in MOVEMENT.

        Returns:
            np.ndarray: The individual energy of atoms in each ionic step structure,
                with shape of (n_ionic_steps, n_atoms).
        """
        return np.array([step["eatoms"] for step in self.ionic_steps if ("eatoms" in step)])

    @property
    def virials(self) -> np.ndarray:
        """Virial tensor of each ionic step structure contained in MOVEMENT.

        Returns:
            np.ndarray: The virial tensor of each ionic step structure,
                with shape of (n_ionic_steps, 3, 3)
        """
        return np.array([step["virial"] for step in self.ionic_steps if ("virial" in step)])

    def _parse_sefv(self) -> list[dict]:
        """Parse the MOVEMENT file, return information ionic step structure containing
        structures, energies, forces on atoms and virial tensor.

        Returns:
            list[dict]: Structure containing structures, energies, forces on atoms
                and virial tensor. The corresponding keys are 'atom_config', 'e_tot',
                'atom_forces' and 'virial'.
        """
        ionic_steps: list[dict] = []
        with zopen(self.filename, "rt") as mvt:
            tmp_step: dict = {}
            for ii in range(self.n_ionic_steps):
                tmp_chunk: str = ""
                for _ in range(self.chunk_sizes[ii]):
                    tmp_chunk += mvt.readline()
                tmp_step["atom_config"] = AtomConfig.from_str(tmp_chunk)
                tmp_step["e_tot"] = ACstrExtractor(tmp_chunk).get_e_tot()[0]
                tmp_step["atom_forces"] = ACstrExtractor(tmp_chunk).get_atom_forces().reshape(-1, 3)
                e_atoms: np.ndarray | None = ACstrExtractor(tmp_chunk).get_atom_forces()
                if e_atoms is not None:
                    tmp_step["atom_energies"] = ACstrExtractor(tmp_chunk).get_atom_energies()
                else:
                    print(f"Ionic step #{ii} : Energy deposition is turn down.")
                virial: np.ndarray | None = ACstrExtractor(tmp_chunk).get_virial()
                if virial is not None:
                    tmp_step["virial"] = virial.reshape(3, 3)
                else:
                    print(f"Ionic step #{ii} : No virial information.")
                ionic_steps.append(tmp_step)
        return ionic_steps


class OutFermi(MSONable):
    """Extract fermi energy (eV) from OUT.FERMI."""

    def __init__(self, filename: PathLike):
        """Initialization function.

        Args:
            filename (PathLike): The absolute path of OUT.FERMI file.
        """
        self.filename: PathLike = filename
        with zopen(self.filename, "rt") as file:
            self._e_fermi: float = np.round(float(file.readline().split()[-2].strip()), 3)

    @property
    def e_fermi(self) -> float:
        """The fermi energy level.

        Returns:
            float: Fermi energy level.
        """
        return self._e_fermi


class Report(MSONable):
    """Extract information of spin, kpoints, bands, eigenvalues from REPORT file."""

    def __init__(self, filename: PathLike):
        """Initialization function.

        Args:
            filename (PathLike): The absolute path of REPORT file.
        """
        self.filename: PathLike = filename
        self._spin, self._num_kpts, self._num_bands = self._parse_band()
        self._eigenvalues = self._parse_eigen()
        self._kpts, self._kpts_weight, self._hsps = self._parse_kpt()

    def _parse_band(self) -> Tuple3Ints:
        """Parse REPORT file to obtain spin switches, the number of kpoints
        and the number of bands.

        Returns:
            tuple[int, int, int]: containing:
                spin (int): Whether turn on spin or not
                    1: turn down the spin
                    2: turn on the spin.
                num_kpts (int): The number of kpoints.
                num_bands (int): The number of bands.
        """
        content: str = "SPIN"
        row_idx: int = LineLocator.locate_all_lines(file_path=self.filename, content=content)[0]
        spin = int(linecache.getline(str(self.filename), row_idx).split()[-1].strip())

        content = "NUM_KPT"
        row_idx = LineLocator.locate_all_lines(file_path=self.filename, content=content)[0]
        num_kpts = int(linecache.getline(str(self.filename), row_idx).split()[-1].strip())

        content = "NUM_BAND"
        row_idx = LineLocator.locate_all_lines(file_path=self.filename, content=content)[0]
        num_bands = int(linecache.getline(str(self.filename), row_idx).split()[-1].strip())
        return spin, num_kpts, num_bands

    def _parse_eigen(self) -> np.ndarray:
        """Parse REPORT file to obtain information about eigenvalues.

        Returns:
            np.ndarray: Eigenvalues with shape of (1 or 2, n_kpoints, n_bands).
                The first index represents spin, the second index represents
                kpoints, the third index represents band.
        """
        num_rows: int = int(np.ceil(self._num_bands / 5))
        content: str = "eigen energies, in eV"
        rows_lst: list[int] = LineLocator.locate_all_lines(file_path=self.filename, content=content)
        rows_array: np.ndarray = np.array(rows_lst).reshape(self._spin, -1)
        eigenvalues: np.ndarray = np.zeros((self._spin, self._num_kpts, self._num_bands))

        for ii in range(self._spin):
            for jj in range(self._num_kpts):
                tmp_eigenvalues_str = ""
                for kk in range(num_rows):
                    tmp_eigenvalues_str += linecache.getline(str(self.filename), rows_array[ii][jj] + kk + 1)
                tmp_eigenvalues_array = np.array([float(eigen_value) for eigen_value in tmp_eigenvalues_str.split()])
                for kk in range(self._num_bands):
                    eigenvalues[ii][jj][kk] = tmp_eigenvalues_array[kk]
        return eigenvalues

    def _parse_kpt(self) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
        """Parse REPORT file to obtain information about kpoints.

        Returns:
            3-tuple containing:
                kpts (np.ndarray): The fractional coordinates of kpoints.
                kpts_weight (np.ndarray): The weight of kpoints.
                hsps (dict[str, np.ndarray]): The name and coordinates of high symmetric points.
        """
        num_rows: int = int(self._num_kpts)
        content: str = "total number of K-point:"
        row_idx: int = LineLocator.locate_all_lines(self.filename, content)[0]
        kpts: np.ndarray = np.zeros((self._num_kpts, 3))
        kpts_weight: np.ndarray = np.zeros(self._num_kpts)
        hsps: dict[str, np.array] = {}
        for ii in range(num_rows):
            #  0.00000     0.00000    0.00000     0.03704           G
            tmp_row_lst: list[str] = linecache.getline(str(self.filename), row_idx + ii + 1).split()
            for jj in range(3):
                kpts[ii][jj] = float(tmp_row_lst[jj].strip())
            kpts_weight[ii] = float(tmp_row_lst[3].strip())

            if len(tmp_row_lst) == 5:
                hsps |= {
                    tmp_row_lst[4]: np.array(
                        [
                            float(tmp_row_lst[0]),
                            float(tmp_row_lst[1]),
                            float(tmp_row_lst[2]),
                        ]
                    )
                }
        return kpts, kpts_weight, hsps

    @property
    def spin(self) -> int:
        """The spin switches.

        Returns:
            int: Spin switches. 1 represents turn on spin, 2 represents turn down spin.
        """
        return self._spin

    @property
    def n_kpoints(self) -> int:
        """The number of k-points."""
        return self._num_kpts

    @property
    def n_bands(self) -> int:
        """The number of bands."""
        return self._num_bands

    @property
    def eigenvalues(self) -> np.ndarray:
        """The eigenvalues.

        Returns:
            np.ndarray: The first index represents spin, the second index
                represents kpoint, the third index represents band.
        """
        return self._eigenvalues

    @property
    def kpoints(self) -> np.ndarray:
        """The fractional coordinates of kpoints."""
        return self._kpts

    @property
    def kpoints_weight(self) -> np.ndarray:
        """The weight of kpoints."""
        return self._kpts_weight

    @property
    def hsps(self) -> dict[str, np.ndarray]:
        """The labels and fractional coordinates of high symmetry points as dict[str, np.ndarray].
        Empty dict when task is not line-mode kpath.
        """
        return self._hsps


class DosSpin(MSONable):
    """Extract information of DOS from DOS_SPIN file:
    - DOS.totalspin, DOS.totalspin_projected
    - DOS.spinup, DOS.spinup_projected
    - DOS.spindown, DOS.spindown_projected.
    """

    def __init__(self, filename: PathLike):
        self.filename: PathLike = filename
        self._labels, self._dos = self._parse()

    def _parse(self):
        """Parse the DOS_SPIN file to get name and values of partial DOS.

        Returns:
            labels (list[str]): The label of DOS, e.g. Total, Cr-3S, ...
            dos (np.array): Value of density of state.
        """
        labels: list[str] = []
        labels = linecache.getline(str(self.filename), 1).split()[1:]
        dos_str: str = ""
        with zopen(self.filename, mode="rt") as file:
            file.readline()
            dos_str = file.read()
        dos: np.ndarray = np.loadtxt(StringIO(dos_str))
        return labels, dos

    @property
    def labels(self) -> list[str]:
        """The name of the partial density of states."""
        return self._labels

    @property
    def dos(self) -> np.ndarray:
        """Value of density of state."""
        return self._dos

    def get_partial_dos(self, part: str) -> np.ndarray:
        """Get partial dos for give element or orbital.

        Args:
            part (str): The name of partial dos.
                e.g. 'Energy', 'Total', 'Cr-3S', 'Cr-3P',
                    'Cr-4S', 'Cr-3D', 'I-4D', 'I-5S', 'I-5P', 'Cr-3S', 'Cr-3Pz',
                    'Cr-3Px', 'Cr-3Py', 'Cr-4S', 'Cr-3Dz2','Cr-3Dxz', 'Cr-3Dyz',
                    'Cr-3D(x^2-y^2)', 'Cr-3Dxy', 'I-4Dz2', 'I-4Dxz', 'I-4Dyz',
                    'I-4D(x^2-y^2)', 'I-4Dxy', 'I-5S', 'I-5Pz', 'I-5Px', 'I-5Py'

        Returns:
            partial_dos: np.array
        """
        part_upper: str = part.upper()
        labels_upper: list[str] = [tmp_label.upper() for tmp_label in self._labels]
        idx_dos = labels_upper.index(part_upper)
        return self._dos[:, idx_dos]
