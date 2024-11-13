"""This module provides classes to define a phonon band structure."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import numpy as np
from monty.json import MSONable

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.bandstructure import Kpoint

if TYPE_CHECKING:
    from collections.abc import Sequence
    from os import PathLike
    from typing import Any

    from numpy.typing import ArrayLike
    from typing_extensions import Self

    from pymatgen.util.typing import Tuple3Ints


def get_reasonable_repetitions(n_atoms: int) -> Tuple3Ints:
    """Choose the number of repetitions in a supercell
    according to the number of atoms in the system.
    """
    if n_atoms < 4:
        return 3, 3, 3
    if 4 <= n_atoms < 15:
        return 2, 2, 2
    if 15 <= n_atoms < 50:
        return 2, 2, 1

    return 1, 1, 1


def eigenvectors_from_displacements(disp: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """Calculate the eigenvectors from the atomic displacements."""
    return np.einsum("nax,a->nax", disp, masses**0.5)


def estimate_band_connection(prev_eigvecs, eigvecs, prev_band_order) -> list[int]:
    """A function to order the phonon eigenvectors taken from phonopy."""
    metric = np.abs(np.dot(prev_eigvecs.conjugate().T, eigvecs))
    connection_order = []
    for overlaps in metric:
        max_val = max_idx = 0
        for idx in reversed(range(len(metric))):
            val = overlaps[idx]
            if idx in connection_order:
                continue
            if val > max_val:
                max_val = val
                max_idx = idx
        connection_order.append(max_idx)

    return [connection_order[x] for x in prev_band_order]


class PhononBandStructure(MSONable):
    """This is the most generic phonon band structure data possible
    it's defined by a list of qpoints + frequencies for each of them.
    Additional information may be given for frequencies at Gamma, where
    non-analytical contribution may be taken into account.
    """

    def __init__(
        self,
        qpoints: Sequence[Kpoint],
        frequencies: ArrayLike,
        lattice: Lattice,
        nac_frequencies: Sequence[Sequence] | None = None,
        eigendisplacements: ArrayLike = None,
        nac_eigendisplacements: Sequence[Sequence] | None = None,
        labels_dict: dict | None = None,
        coords_are_cartesian: bool = False,
        structure: Structure | None = None,
    ) -> None:
        """
        Args:
            qpoints: list of qpoint as numpy arrays, in frac_coords of the
                given lattice by default
            frequencies: list of phonon frequencies in THz as a numpy array with shape
                (3*len(structure), len(qpoints)). The First index of the array
                refers to the band and the second to the index of the qpoint.
            lattice: The reciprocal lattice as a pymatgen Lattice object.
                Pymatgen uses the physics convention of reciprocal lattice vectors
                WITH a 2*pi coefficient.
            nac_frequencies: Frequencies with non-analytical contributions at Gamma in THz.
                A list of tuples. The first element of each tuple should be a list
                defining the direction (not necessarily a versor, will be normalized
                internally). The second element containing the 3*len(structure)
                phonon frequencies with non-analytical correction for that direction.
            eigendisplacements: the phonon eigendisplacements associated to the
                frequencies in Cartesian coordinates. A numpy array of complex
                numbers with shape (3*len(structure), len(qpoints), len(structure), 3).
                The first index of the array refers to the band, the second to the index
                of the qpoint, the third to the atom in the structure and the fourth
                to the Cartesian coordinates.
            nac_eigendisplacements: the phonon eigendisplacements associated to the
                non-analytical frequencies in nac_frequencies in Cartesian coordinates.
                A list of tuples. The first element of each tuple should be a list
                defining the direction. The second element containing a numpy array of
                complex numbers with shape (3*len(structure), len(structure), 3).
            labels_dict: (dict[str, Kpoint]): this links a qpoint (in frac coords or
                Cartesian coordinates depending on the coords) to a label.
            coords_are_cartesian (bool): Whether the qpoint coordinates are Cartesian. Defaults to False.
            structure: The crystal structure (as a pymatgen Structure object)
                associated with the band structure. This is needed to calculate element/orbital
                projections of the band structure.
        """
        self.lattice_rec = lattice
        self.qpoints: list[Kpoint] = []
        self.labels_dict = {}
        self.structure = structure
        if eigendisplacements is None:
            eigendisplacements = np.array([])
        self.eigendisplacements = eigendisplacements
        if labels_dict is None:
            labels_dict = {}

        for q_pt in qpoints:
            label = None  # check below if this qpoint has an assigned label
            for key in labels_dict:
                if np.linalg.norm(q_pt - np.array(labels_dict[key])) < 0.0001:
                    label = key
                    self.labels_dict[label] = Kpoint(
                        q_pt,
                        lattice,
                        label=label,
                        coords_are_cartesian=coords_are_cartesian,
                    )
            self.qpoints += [
                Kpoint(
                    q_pt,
                    lattice,
                    label=label,
                    coords_are_cartesian=coords_are_cartesian,
                )
            ]
        self.bands = np.asarray(frequencies)
        self.nb_bands = len(self.bands)
        self.nb_qpoints = len(self.qpoints)

        # normalize directions for nac_frequencies and nac_eigendisplacements
        self.nac_frequencies: list[tuple[list[float], np.ndarray]] = []
        self.nac_eigendisplacements: list[tuple[list[float], np.ndarray]] = []
        if nac_frequencies is not None:
            for freq in nac_frequencies:
                self.nac_frequencies.append(([idx / np.linalg.norm(freq[0]) for idx in freq[0]], freq[1]))
        if nac_eigendisplacements is not None:
            for freq in nac_eigendisplacements:
                self.nac_eigendisplacements.append(([idx / np.linalg.norm(freq[0]) for idx in freq[0]], freq[1]))

    def get_gamma_point(self) -> Kpoint | None:
        """Get the Gamma q-point as a Kpoint object (or None if not found)."""
        for q_point in self.qpoints:
            if np.allclose(q_point.frac_coords, (0, 0, 0)):
                return q_point

        return None

    def min_freq(self) -> tuple[Kpoint, float]:
        """Get the q-point where the minimum frequency is reached and its value."""
        idx = np.unravel_index(np.argmin(self.bands), self.bands.shape)

        return self.qpoints[idx[1]], self.bands[idx]

    def max_freq(self) -> tuple[Kpoint, float]:
        """Get the q-point where the maximum frequency is reached and its value."""
        idx = np.unravel_index(np.argmax(self.bands), self.bands.shape)

        return self.qpoints[idx[1]], self.bands[idx]

    def width(self, with_imaginary: bool = False) -> float:
        """Get the difference between the maximum and minimum frequencies anywhere in the
        band structure, not necessarily at identical same q-points. If with_imaginary is False,
        only positive frequencies are considered.
        """
        if with_imaginary:
            return np.max(self.bands) - np.min(self.bands)
        mask_pos = self.bands >= 0
        return self.bands[mask_pos].max() - self.bands[mask_pos].min()

    def has_imaginary_freq(self, tol: float = 0.01) -> bool:
        """True if imaginary frequencies are present anywhere in the band structure. Always True if
        has_imaginary_gamma_freq is True.

        Args:
            tol: Tolerance for determining if a frequency is imaginary. Defaults to 0.01.
        """
        return bool(self.min_freq()[1] + tol < 0)

    def has_imaginary_gamma_freq(self, tol: float = 0.01) -> bool:
        """Check if there are imaginary modes at the gamma point and all close points.

        Args:
            tol: Tolerance for determining if a frequency is imaginary. Defaults to 0.01.
        """
        # Calculate the radial distance from the gamma point for each q-point
        close_points = [q_pt for q_pt in self.qpoints if np.linalg.norm(q_pt.frac_coords) < tol]

        # check for negative frequencies at all q-points close to the gamma point
        for qpoint in close_points:
            idx = self.qpoints.index(qpoint)
            if any(freq < -tol for freq in self.bands[:, idx]):
                return True

        return False

    @property
    def has_nac(self) -> bool:
        """True if nac_frequencies are present (i.e. the band structure has been
        calculated taking into account Born-charge-derived non-analytical corrections at Gamma).
        """
        return len(self.nac_frequencies) > 0

    @property
    def has_eigendisplacements(self) -> bool:
        """True if eigendisplacements are present."""
        return len(self.eigendisplacements) > 0

    def get_nac_frequencies_along_dir(self, direction: Sequence) -> np.ndarray | None:
        """Get the nac_frequencies for the given direction (not necessarily a versor).
        None if the direction is not present or nac_frequencies has not been calculated.

        Args:
            direction: the direction as a list of 3 elements

        Returns:
            the frequencies as a numpy array o(3*len(structure), len(qpoints)).
            None if not found.
        """
        versor = [idx / np.linalg.norm(direction) for idx in direction]
        for dist, freq in self.nac_frequencies:
            if np.allclose(versor, dist):
                return freq

        return None

    def get_nac_eigendisplacements_along_dir(self, direction) -> np.ndarray | None:
        """Get the nac_eigendisplacements for the given direction (not necessarily a versor).
        None if the direction is not present or nac_eigendisplacements has not been calculated.

        Args:
            direction: the direction as a list of 3 elements

        Returns:
            the eigendisplacements as a numpy array of complex numbers with shape
            (3*len(structure), len(structure), 3). None if not found.
        """
        versor = [idx / np.linalg.norm(direction) for idx in direction]
        for dist, eigen_disp in self.nac_eigendisplacements:
            if np.allclose(versor, dist):
                return eigen_disp

        return None

    def asr_breaking(self, tol_eigendisplacements: float = 1e-5) -> np.ndarray | None:
        """Get the breaking of the acoustic sum rule for the three acoustic modes,
        if Gamma is present. None otherwise.
        If eigendisplacements are available they are used to determine the acoustic
        modes: selects the bands corresponding to the eigendisplacements that
        represent to a translation within tol_eigendisplacements. If these are not
        identified or eigendisplacements are missing the first 3 modes will be used
        (indices [:3]).
        """
        for idx in range(self.nb_qpoints):
            if np.allclose(self.qpoints[idx].frac_coords, (0, 0, 0)):
                if self.has_eigendisplacements:
                    acoustic_modes_index = []
                    for j in range(self.nb_bands):
                        eig = self.eigendisplacements[j][idx]
                        if np.max(np.abs(eig[1:] - eig[:1])) < tol_eigendisplacements:
                            acoustic_modes_index.append(j)
                    # if acoustic modes are not correctly identified return use
                    # the first three modes
                    if len(acoustic_modes_index) != 3:
                        acoustic_modes_index = [0, 1, 2]
                    return self.bands[acoustic_modes_index, idx]

                return self.bands[:3, idx]

        return None

    def as_dict(self) -> dict[str, Any]:
        """MSONable dict."""
        dct: dict[str, Any] = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "lattice_rec": self.lattice_rec.as_dict(),
            # qpoints are not Kpoint objects dicts but are frac coords. This makes
            # the dict smaller and avoids the repetition of the lattice
            "qpoints": [q_pt.as_dict()["fcoords"] for q_pt in self.qpoints],
        }
        dct["bands"] = self.bands.tolist()
        dct["labels_dict"] = {}
        for kpoint_letter, kpoint_object in self.labels_dict.items():
            dct["labels_dict"][kpoint_letter] = kpoint_object.as_dict()["fcoords"]

        # split the eigendisplacements to real and imaginary part for serialization
        dct["eigendisplacements"] = {
            "real": np.real(self.eigendisplacements).tolist(),
            "imag": np.imag(self.eigendisplacements).tolist(),
        }
        dct["nac_eigendisplacements"] = [
            (direction, {"real": np.real(e).tolist(), "imag": np.imag(e).tolist()})
            for direction, e in self.nac_eigendisplacements
        ]
        dct["nac_frequencies"] = [(direction, f.tolist()) for direction, f in self.nac_frequencies]

        if self.structure:
            dct["structure"] = self.structure.as_dict()

        return dct

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        """
        Args:
            dct (dict): Dict representation of PhononBandStructure.

        Returns:
            PhononBandStructure
        """
        lattice_rec = Lattice(dct["lattice_rec"]["matrix"])
        eigendisplacements = (
            np.array(dct["eigendisplacements"]["real"]) + np.array(dct["eigendisplacements"]["imag"]) * 1j
        )
        nac_eigendisplacements = [
            (direction, np.array(e["real"]) + np.array(e["imag"]) * 1j)
            for direction, e in dct["nac_eigendisplacements"]
        ]
        nac_frequencies = [(direction, np.array(f)) for direction, f in dct["nac_frequencies"]]
        structure = Structure.from_dict(dct["structure"]) if "structure" in dct else None
        return cls(
            dct["qpoints"],
            np.array(dct["bands"]),
            lattice_rec,
            nac_frequencies,
            eigendisplacements,
            nac_eigendisplacements,
            dct["labels_dict"],
            structure=structure,
        )


class PhononBandStructureSymmLine(PhononBandStructure):
    r"""Store phonon band structures along selected (symmetry) lines in the Brillouin zone.
    We call the different symmetry lines (ex: \\Gamma to Z) "branches".
    """

    def __init__(
        self,
        qpoints: Sequence[Kpoint],
        frequencies: ArrayLike,
        lattice: Lattice,
        has_nac: bool = False,
        eigendisplacements: ArrayLike = None,
        labels_dict: dict | None = None,
        coords_are_cartesian: bool = False,
        structure: Structure | None = None,
    ) -> None:
        """
        Args:
            qpoints: list of qpoints as numpy arrays, in frac_coords of the
                given lattice by default
            frequencies: list of phonon frequencies in eV as a numpy array with shape
                (3*len(structure), len(qpoints))
            lattice: The reciprocal lattice as a pymatgen Lattice object.
                Pymatgen uses the physics convention of reciprocal lattice vectors
                WITH a 2*pi coefficient
            has_nac: specify if the band structure has been produced taking into account
                non-analytical corrections at Gamma. If True frequencies at Gamma from
                different directions will be stored in naf. Default False.
            eigendisplacements: the phonon eigendisplacements associated to the
                frequencies in Cartesian coordinates. A numpy array of complex
                numbers with shape (3*len(structure), len(qpoints), len(structure), 3).
                he First index of the array refers to the band, the second to the index
                of the qpoint, the third to the atom in the structure and the fourth
                to the Cartesian coordinates.
            labels_dict: (dict) of {} this links a qpoint (in frac coords or
                Cartesian coordinates depending on the coords) to a label.
            coords_are_cartesian: Whether the qpoint coordinates are cartesian.
            structure: The crystal structure (as a pymatgen Structure object)
                associated with the band structure. This is needed if we
                provide projections to the band structure.
        """
        super().__init__(
            qpoints=qpoints,
            frequencies=frequencies,
            lattice=lattice,
            nac_frequencies=None,
            eigendisplacements=eigendisplacements,
            nac_eigendisplacements=None,
            labels_dict=labels_dict,
            coords_are_cartesian=coords_are_cartesian,
            structure=structure,
        )
        self._reuse_init(eigendisplacements, frequencies, has_nac, qpoints)

    def __repr__(self) -> str:
        bands, labels = self.bands.shape, list(self.labels_dict)
        return f"{type(self).__name__}({bands=}, {labels=})"

    def _reuse_init(
        self,
        eigendisplacements: ArrayLike,
        frequencies: ArrayLike,
        has_nac: bool,
        qpoints: Sequence[Kpoint],
    ) -> None:
        self.distance = []
        self.branches = []
        one_group: list = []
        branches_tmp = []
        # get labels and distance for each qpoint
        previous_qpoint = self.qpoints[0]
        previous_distance = 0.0
        previous_label = self.qpoints[0].label
        for idx in range(self.nb_qpoints):
            label = self.qpoints[idx].label
            if label is not None and previous_label is not None:
                self.distance += [previous_distance]
            else:
                self.distance += [
                    np.linalg.norm(self.qpoints[idx].cart_coords - previous_qpoint.cart_coords) + previous_distance
                ]
            previous_qpoint = self.qpoints[idx]
            previous_distance = self.distance[idx]
            if label and previous_label:
                if len(one_group) != 0:
                    branches_tmp += [one_group]
                one_group = []
            previous_label = label
            one_group += [idx]
        if len(one_group) != 0:
            branches_tmp += [one_group]
        for branch in branches_tmp:
            self.branches += [
                {
                    "start_index": branch[0],
                    "end_index": branch[-1],
                    "name": f"{self.qpoints[branch[0]].label}-{self.qpoints[branch[-1]].label}",
                }
            ]
        # extract the frequencies with non-analytical contribution at gamma
        if has_nac:
            naf = []
            nac_eigendisplacements = []
            for idx in range(self.nb_qpoints):
                # get directions with nac irrespectively of the label_dict. NB: with labels
                # the gamma point is expected to appear twice consecutively.
                if np.allclose(qpoints[idx], (0, 0, 0)):
                    if idx > 0 and not np.allclose(qpoints[idx - 1], (0, 0, 0)):
                        q_dir = self.qpoints[idx - 1]
                        direction = q_dir.frac_coords / np.linalg.norm(q_dir.frac_coords)
                        naf.append((direction, frequencies[:, idx]))
                        if self.has_eigendisplacements:
                            nac_eigendisplacements.append((direction, eigendisplacements[:, idx]))
                    if idx < len(qpoints) - 1 and not np.allclose(qpoints[idx + 1], (0, 0, 0)):
                        q_dir = self.qpoints[idx + 1]
                        direction = q_dir.frac_coords / np.linalg.norm(q_dir.frac_coords)
                        naf.append((direction, frequencies[:, idx]))
                        if self.has_eigendisplacements:
                            nac_eigendisplacements.append((direction, eigendisplacements[:, idx]))

            self.nac_frequencies = np.array(naf, dtype=object)
            self.nac_eigendisplacements = np.array(nac_eigendisplacements, dtype=object)

    def get_equivalent_qpoints(self, index: int) -> list[int]:
        """Get the list of qpoint indices equivalent (meaning they are the
        same frac coords) to the given one.

        Args:
            index (int): the qpoint index

        Returns:
            list[int]: equivalent indices

        TODO: now it uses the label we might want to use coordinates instead
        (in case there was a mislabel)
        """
        # if the qpoint has no label it can't have a repetition along the band
        # structure line object

        if self.qpoints[index].label is None:
            return [index]

        list_index_qpoints = []
        for idx in range(self.nb_qpoints):
            if self.qpoints[idx].label == self.qpoints[index].label:
                list_index_qpoints.append(idx)

        return list_index_qpoints

    def get_branch(self, index: int) -> list[dict[str, str | int]]:
        r"""Get in what branch(es) is the qpoint. There can be several branches.

        Args:
            index (int): the qpoint index

        Returns:
            list[dict[str, str | int]]: [{"name","start_index","end_index","index"}]
                indicating all branches in which the qpoint is. It takes into account
                the fact that one qpoint (e.g., \\Gamma) can be in several branches
        """
        lst = []
        for pt_idx in self.get_equivalent_qpoints(index):
            for branch in self.branches:
                start_idx, end_idx = branch["start_index"], branch["end_index"]
                if start_idx <= pt_idx <= end_idx:
                    lst.append(
                        {
                            "name": branch["name"],
                            "start_index": start_idx,
                            "end_index": end_idx,
                            "index": pt_idx,
                        }
                    )
        return lst

    def write_phononwebsite(self, filename: str | PathLike) -> None:
        """Write a JSON file for the phononwebsite:
        https://henriquemiranda.github.io/phononwebsite.
        """
        with open(filename, mode="w") as file:
            json.dump(self.as_phononwebsite(), file)

    def as_phononwebsite(self) -> dict:
        """Return a dictionary with the phononwebsite format:
        https://henriquemiranda.github.io/phononwebsite.
        """
        if self.structure is None:
            raise RuntimeError("Structure is required for as_phononwebsite")
        dct = {}

        # define the lattice
        dct["lattice"] = self.structure.lattice._matrix.tolist()

        # define atoms
        atom_pos_car = []
        atom_pos_red = []
        atom_types = []
        for site in self.structure:
            atom_pos_car.append(site.coords.tolist())
            atom_pos_red.append(site.frac_coords.tolist())
            atom_types.append(site.species_string)

        # default for now
        dct["repetitions"] = get_reasonable_repetitions(len(atom_pos_car))

        dct["natoms"] = len(atom_pos_car)
        dct["atom_pos_car"] = atom_pos_car
        dct["atom_pos_red"] = atom_pos_red
        dct["atom_types"] = atom_types
        dct["atom_numbers"] = self.structure.atomic_numbers
        dct["formula"] = self.structure.formula
        dct["name"] = self.structure.formula

        # get qpoints
        qpoints = []
        for q_pt in self.qpoints:
            qpoints.append(list(q_pt.frac_coords))
        dct["qpoints"] = qpoints

        # get labels
        hsq_dict = {}
        for nq, q_pt in enumerate(self.qpoints):
            if q_pt.label is not None:
                hsq_dict[nq] = q_pt.label

        # get distances
        dist = nq_start = 0
        distances = [dist]
        line_breaks = []
        for nq in range(1, len(qpoints)):
            q1 = np.array(qpoints[nq])
            q2 = np.array(qpoints[nq - 1])
            # detect jumps
            if (nq in hsq_dict) and (nq - 1 in hsq_dict):
                if hsq_dict[nq] != hsq_dict[nq - 1]:
                    hsq_dict[nq - 1] += "|" + hsq_dict[nq]
                del hsq_dict[nq]
                line_breaks.append((nq_start, nq))
                nq_start = nq
            else:
                dist += np.linalg.norm(q1 - q2)
            distances.append(dist)
        line_breaks.append((nq_start, len(qpoints)))
        dct["distances"] = distances
        dct["line_breaks"] = line_breaks
        dct["highsym_qpts"] = list(hsq_dict.items())

        # eigenvalues
        thz2cm1 = 33.35641
        bands = self.bands.copy() * thz2cm1
        dct["eigenvalues"] = bands.T.tolist()

        # eigenvectors
        eigen_vecs = self.eigendisplacements.copy()
        eigen_vecs /= np.linalg.norm(eigen_vecs[0, 0])
        eigen_vecs = eigen_vecs.swapaxes(0, 1)
        eigen_vecs = np.array([eigen_vecs.real, eigen_vecs.imag])
        eigen_vecs = np.rollaxis(eigen_vecs, 0, 5)
        dct["vectors"] = eigen_vecs.tolist()

        return dct

    def band_reorder(self) -> None:
        """Re-order the eigenvalues according to the similarity of the eigenvectors."""
        eigen_displacements = self.eigendisplacements
        eig = self.bands

        n_phonons, n_qpoints = self.bands.shape
        order = np.zeros([n_qpoints, n_phonons], dtype=np.int64)
        order[0] = np.array(range(n_phonons))

        # Get the atomic masses
        if self.structure is None:
            raise RuntimeError("Structure is required for band_reorder")
        atomic_masses = [site.specie.atomic_mass for site in self.structure]

        # Get order
        for nq in range(1, n_qpoints):
            old_eig_vecs = eigenvectors_from_displacements(eigen_displacements[:, nq - 1], atomic_masses)
            new_eig_vecs = eigenvectors_from_displacements(eigen_displacements[:, nq], atomic_masses)
            order[nq] = estimate_band_connection(
                old_eig_vecs.reshape([n_phonons, n_phonons]).T,
                new_eig_vecs.reshape([n_phonons, n_phonons]).T,
                order[nq - 1],
            )

        # reorder
        for nq in range(1, n_qpoints):
            eivq = eigen_displacements[:, nq]
            eigq = eig[:, nq]
            eigen_displacements[:, nq] = eivq[order[nq]]
            eig[:, nq] = eigq[order[nq]]

    def as_dict(self) -> dict:
        """Get MSONable dict."""
        dct = super().as_dict()
        # remove nac_frequencies and nac_eigendisplacements as they are reconstructed
        # in the __init__ when the dict is deserialized
        nac_frequencies = dct.pop("nac_frequencies")
        dct.pop("nac_eigendisplacements")
        dct["has_nac"] = len(nac_frequencies) > 0
        return dct

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            PhononBandStructureSymmLine
        """
        lattice_rec = Lattice(dct["lattice_rec"]["matrix"])
        eigendisplacements = (
            np.array(dct["eigendisplacements"]["real"]) + np.array(dct["eigendisplacements"]["imag"]) * 1j
        )
        return cls(
            dct["qpoints"],
            np.array(dct["bands"]),
            lattice_rec,
            dct["has_nac"],
            eigendisplacements,
            dct["labels_dict"],
            structure=(Structure.from_dict(dct["structure"]) if "structure" in dct else None),
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, PhononBandStructureSymmLine):
            return NotImplemented
        return (
            self.bands.shape == other.bands.shape
            and np.allclose(self.bands, other.bands)
            and self.lattice_rec == other.lattice_rec
            # and self.qpoints == other.qpoints
            and self.labels_dict == other.labels_dict
            and self.structure == other.structure
        )
