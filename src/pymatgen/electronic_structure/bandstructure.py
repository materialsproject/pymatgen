"""This module provides classes to define things related to band structures."""

from __future__ import annotations

import itertools
import math
import re
import warnings
from collections import defaultdict
from typing import TYPE_CHECKING, overload

import numpy as np
from monty.json import MSONable

from pymatgen.core import Element, Lattice, Structure, get_el_sp
from pymatgen.electronic_structure.core import Orbital, Spin
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord import pbc_diff

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any

    from numpy.typing import ArrayLike, NDArray
    from typing_extensions import Self

    from pymatgen.util.typing import SpeciesLike

__author__ = "Geoffroy Hautier, Shyue Ping Ong, Michael Kocher"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "March 14, 2012"


class Kpoint(MSONable):
    """A kpoint defined with a lattice and frac or Cartesian coordinates,
    similar to the Site object in pymatgen.core.structure.
    """

    def __init__(
        self,
        coords: ArrayLike,
        lattice: Lattice,
        to_unit_cell: bool = False,
        coords_are_cartesian: bool = False,
        label: str | None = None,
    ) -> None:
        """
        Args:
            coords (NDArray): Coordinate of the Kpoint.
            lattice (Lattice): The reciprocal lattice of the kpoint.
            to_unit_cell (bool): Translate fractional coordinate to the basic unit
                cell, i.e., all fractional coordinates satisfy 0 <= a < 1.
                Defaults to False.
            coords_are_cartesian (bool): Whether the coordinates given are
                in Cartesian (True) or fractional coordinates (by default fractional).
            label (str): The label of the Kpoint if any (None by default).
        """
        self._lattice = lattice
        self._frac_coords = lattice.get_fractional_coords(coords) if coords_are_cartesian else np.asarray(coords)
        self._label = label

        if to_unit_cell:
            for idx, fc in enumerate(self._frac_coords):
                self._frac_coords[idx] -= math.floor(fc)

        self._cart_coords = lattice.get_cartesian_coords(self._frac_coords)

    def __str__(self) -> str:
        """String with fractional, Cartesian coordinates and label."""
        return f"{self.frac_coords} {self.cart_coords} {self.label}"

    def __eq__(self, other: object) -> bool:
        """Whether two Kpoints are equal."""
        if not isinstance(other, type(self)):
            return NotImplemented

        return (
            np.allclose(self.frac_coords, other.frac_coords)
            and self.lattice == other.lattice
            and self.label == other.label
        )

    @property
    def lattice(self) -> Lattice:
        """The lattice associated with the kpoint, as a Lattice object."""
        return self._lattice

    @property
    def label(self) -> str | None:
        """The label associated with the kpoint."""
        return self._label

    @label.setter
    def label(self, label: str | None) -> None:
        """Set the label of the kpoint."""
        self._label = label

    @property
    def frac_coords(self) -> NDArray:
        """The fractional coordinates of the kpoint as a NumPy array."""
        return np.copy(self._frac_coords)

    @property
    def cart_coords(self) -> NDArray:
        """The Cartesian coordinates of the kpoint as a NumPy array."""
        return np.copy(self._cart_coords)

    @property
    def a(self) -> float:
        """Fractional a coordinate of the kpoint."""
        return self._frac_coords[0]

    @property
    def b(self) -> float:
        """Fractional b coordinate of the kpoint."""
        return self._frac_coords[1]

    @property
    def c(self) -> float:
        """Fractional c coordinate of the kpoint."""
        return self._frac_coords[2]

    def as_dict(self) -> dict[str, Any]:
        """JSON-serializable dict representation of the kpoint."""
        return {
            "lattice": self.lattice.as_dict(),
            "fcoords": self.frac_coords.tolist(),
            "ccoords": self.cart_coords.tolist(),
            "label": self.label,
            "@module": type(self).__module__,
            "@class": type(self).__name__,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Create from a dict.

        Args:
            dct (dict): A dict with all data for a kpoint object.

        Returns:
            Kpoint
        """
        lattice = Lattice.from_dict(dct["lattice"])
        return cls(
            coords=dct["fcoords"],
            lattice=lattice,
            coords_are_cartesian=False,
            label=dct["label"],
        )


class BandStructure:
    """Generic band structure data, defined by a list of Kpoints
        and corresponding energies for each of them.

    Attributes:
        kpoints (list[Kpoint]): Kpoints in the band structure.
        lattice_rec (Lattice): The reciprocal lattice of the band structure.
        efermi (float): The Fermi level.
        is_spin_polarized (bool): Whether the band structure is spin-polarized.
        bands (dict[Spin, NDArray]): The energy eigenvalues. Note that the use of an
            array is necessary for computational and memory efficiency due to the large
            amount of numerical data. The indices of the array are (band_index, kpoint_index).
        nb_bands (int): The number of bands in the band structure.
        structure (Structure): The structure.
        projections (dict[Spin, NDArray]): The projections. Note that the use of an
            array is necessary for computational and memory efficiency due to the large
            amount of numerical data. The indices of the array are (band_index, kpoint_index,
            orbital_index, ion_index).
    """

    def __init__(
        self,
        kpoints: ArrayLike,
        eigenvals: Mapping[Spin, ArrayLike],
        lattice: Lattice,
        efermi: float,
        labels_dict: Mapping[str, Kpoint] | None = None,
        coords_are_cartesian: bool = False,
        structure: Structure | None = None,
        projections: Mapping[Spin, NDArray] | None = None,
    ) -> None:
        """
        Args:
            kpoints (NDArray): Kpoint as NumPy array, in frac_coords of the
                given lattice by default.
            eigenvals (dict): Energies for spin up and spin down as
                {Spin.up:[][], Spin.down:[][]}, the first index of the array
                [][] refers to the band and the second to the index of the kpoint.
                The kpoints are ordered according to the kpoints array.
                If the band structure is not spin polarized, we
                only store one data set under Spin.up.
            lattice (Lattice): The reciprocal lattice. Pymatgen uses the physics
                convention of reciprocal lattice vectors with a 2*pi coefficient.
            efermi (float): The Fermi level.
            labels_dict (dict[str, Kpoint]): Dict mapping label to Kpoint.
            coords_are_cartesian (bool): Whether coordinates are cartesian.
            structure (Structure): The crystal structure
                associated with the band structure. This is needed if we
                provide projections to the band structure.
            projections (dict[Spin, NDArray]): Orbital projections. The
                indices of the array are (band_index, kpoint_index, orbital_index,
                ion_index). If the band structure is not spin polarized, we only
                store one data set under Spin.up.
        """
        self.efermi = efermi
        self.lattice_rec = lattice
        self.kpoints = []
        self.labels_dict = {}
        self.structure = structure
        self.projections = projections or {}
        self.projections = {k: np.array(v) for k, v in self.projections.items()}

        if labels_dict is None:
            labels_dict = {}

        if self.projections and self.structure is None:
            raise RuntimeError("if projections are provided a structure object is also required")

        for kpt in kpoints:
            # Check if this Kpoint has a label
            label = None
            for c in labels_dict:
                if np.linalg.norm(kpt - np.array(labels_dict[c])) < 0.0001:
                    label = c
                    self.labels_dict[label] = Kpoint(
                        kpt,
                        lattice,
                        label=label,
                        coords_are_cartesian=coords_are_cartesian,
                    )
            self.kpoints.append(Kpoint(kpt, lattice, label=label, coords_are_cartesian=coords_are_cartesian))
        self.bands = {spin: np.array(v) for spin, v in eigenvals.items()}
        self.nb_bands = len(self.bands[Spin.up])
        self.is_spin_polarized = len(self.bands) == 2

    def get_projection_on_elements(self) -> dict[Spin, list[list[dict[str, float]]]]:
        """Get projections on elements.

        Returns:
            dict[Spin, NDArray]: Dict in {Spin.up:[][{Element: [values]}],
                Spin.down: [][{Element: [values]}]} format.
                If there is no projections in the band structure, return {}.
        """
        if self.structure is None:
            raise ValueError("structure is None.")
        result: dict[Spin, list[list[dict[str, float]]]] = {}
        for spin, val in self.projections.items():
            result[spin] = [[defaultdict(float) for _ in range(len(self.kpoints))] for _ in range(self.nb_bands)]
            for i, j, k in itertools.product(
                range(self.nb_bands),
                range(len(self.kpoints)),
                range(len(self.structure)),
            ):
                result[spin][i][j][str(self.structure[k].specie)] += np.sum(val[i, j, :, k])
        return result

    def get_projections_on_elements_and_orbitals(self, el_orb_spec: dict[str, list[str]]):
        """Get projections on elements and specific orbitals.

        Args:
            el_orb_spec (dict[str, list[str]]): Elements and orbitals to project onto.
                Format is {Element: [orbitals]}, e.g. {"Cu": ["d", "s"]}.

        Returns:
            dict[str, list[str]: Projections on elements in the
                {Spin.up: [][{Element: {orb: values}}],
                Spin.down: [][{Element: {orb: values}}]} format.
                If there is no projections in the band structure, return {}.
        """
        if self.structure is None:
            raise ValueError("Structure is required for this method")
        result: dict[Spin, list] = {}
        species_orb_spec = {get_el_sp(el): orbs for el, orbs in el_orb_spec.items()}
        for spin, v in self.projections.items():
            result[spin] = [
                [{str(e): defaultdict(float) for e in species_orb_spec} for _ in range(len(self.kpoints))]
                for _ in range(self.nb_bands)
            ]

            for i, j, k in itertools.product(
                range(self.nb_bands),
                range(len(self.kpoints)),
                range(len(self.structure)),
            ):
                sp = self.structure[k].specie
                for orb_i in range(len(v[i][j])):
                    o = Orbital(orb_i).name[0]
                    if sp in species_orb_spec and o in species_orb_spec[sp]:
                        result[spin][i][j][str(sp)][o] += v[i][j][orb_i][k]
        return result

    def is_metal(self, efermi_tol: float = 1e-4) -> bool:
        """Check if the band structure indicates a metal,
        by looking at if the fermi level crosses a band.

        Returns:
            bool: True if is metal.
        """
        for vals in self.bands.values():
            for idx in range(self.nb_bands):
                if np.any(vals[idx, :] - self.efermi < -efermi_tol) and np.any(vals[idx, :] - self.efermi > efermi_tol):
                    return True
        return False

    def get_vbm(self) -> dict[str, Any]:
        """Get data about the valence band maximum (VBM).

        Returns:
            dict with keys "band_index", "kpoint_index", "kpoint", "energy":
                - "band_index" (dict): A dict with spin keys pointing to a list of the
                indices of the band containing the VBM (please note that you
                can have several bands sharing the VBM) {Spin.up:[],
                Spin.down:[]}.
                - "kpoint_index": The list of indices in self.kpoints for the
                kpoint VBM. Please note that there can be several
                kpoint_indices relating to the same kpoint (e.g., Gamma can
                occur at different spots in the band structure line plot).
                - "kpoint" (Kpoint): The kpoint.
                - "energy" (float): The energy of the VBM.
                - "projections": The projections along sites and orbitals of the
                VBM if any projection data is available (else it is an empty
                dictionary). The format is similar to the projections field in
                BandStructure: {spin:{'Orbital': [proj]}} where the array
                [proj] is ordered according to the sites in structure.
        """
        if self.is_metal():
            return {
                "band_index": [],
                "kpoint_index": [],
                "kpoint": [],
                "energy": None,
                "projections": {},
            }

        max_tmp = -float("inf")
        index = kpoint_vbm = None
        for value in self.bands.values():
            for idx, j in zip(*np.where(value < self.efermi), strict=True):
                if value[idx, j] > max_tmp:
                    max_tmp = float(value[idx, j])
                    index = j
                    kpoint_vbm = self.kpoints[j]

        list_ind_kpts = []
        if kpoint_vbm is not None and kpoint_vbm.label is not None:
            for idx, kpt in enumerate(self.kpoints):
                if kpt.label == kpoint_vbm.label:
                    list_ind_kpts.append(idx)
        else:
            list_ind_kpts.append(index)

        # Get all other bands sharing the VBM
        list_ind_band = defaultdict(list)
        for spin in self.bands:
            for idx in range(self.nb_bands):
                if math.isclose(self.bands[spin][idx][index], max_tmp, abs_tol=1e-3, rel_tol=0):
                    list_ind_band[spin].append(idx)
        proj = {}
        for spin, value in self.projections.items():
            if len(list_ind_band[spin]) == 0:
                continue
            proj[spin] = value[list_ind_band[spin][0]][list_ind_kpts[0]]

        return {
            "band_index": list_ind_band,
            "kpoint_index": list_ind_kpts,
            "kpoint": kpoint_vbm,
            "energy": max_tmp,
            "projections": proj,
        }

    def get_cbm(self) -> dict[str, Any]:
        """Get data about the conduction band minimum (CBM).

        Returns:
            dict with keys "band_index", "kpoint_index", "kpoint", "energy":
                - "band_index" (dict): A dict with spin keys pointing to a list of the
                indices of the band containing the CBM (please note that you
                can have several bands sharing the CBM) {Spin.up:[], Spin.down:[]}.
                - "kpoint_index": The list of indices in self.kpoints for the
                kpoint CBM. Please note that there can be several
                kpoint_indices relating to the same kpoint (e.g., Gamma can
                occur at different spots in the band structure line plot).
                - "kpoint" (Kpoint): The kpoint.
                - "energy" (float): The energy of the CBM.
                - "projections": The projections along sites and orbitals of the
                CBM if any projection data is available (else it is an empty
                dictionary). The format is similar to the projections field in
                BandStructure: {spin:{'Orbital': [proj]}} where the array
                [proj] is ordered according to the sites in structure.
        """
        if self.is_metal():
            return {
                "band_index": [],
                "kpoint_index": [],
                "kpoint": [],
                "energy": None,
                "projections": {},
            }

        max_tmp = float("inf")
        index = kpoint_cbm = None
        for value in self.bands.values():
            for idx, j in zip(*np.where(value >= self.efermi), strict=True):
                if value[idx, j] < max_tmp:
                    max_tmp = float(value[idx, j])
                    index = j
                    kpoint_cbm = self.kpoints[j]

        list_index_kpoints = []
        if kpoint_cbm is not None and kpoint_cbm.label is not None:
            for idx, kpt in enumerate(self.kpoints):
                if kpt.label == kpoint_cbm.label:
                    list_index_kpoints.append(idx)
        else:
            list_index_kpoints.append(index)

        # Get all other bands sharing the CBM
        list_index_band = defaultdict(list)
        for spin in self.bands:
            for idx in range(self.nb_bands):
                if math.isclose(self.bands[spin][idx][index], max_tmp, abs_tol=1e-3, rel_tol=0):
                    list_index_band[spin].append(idx)
        proj = {}
        for spin, value in self.projections.items():
            if len(list_index_band[spin]) == 0:
                continue
            proj[spin] = value[list_index_band[spin][0]][list_index_kpoints[0]]

        return {
            "band_index": list_index_band,
            "kpoint_index": list_index_kpoints,
            "kpoint": kpoint_cbm,
            "energy": max_tmp,
            "projections": proj,
        }

    def get_band_gap(self) -> dict[str, Any]:
        r"""Get band gap.

        Returns:
            dict with keys "energy", "direct", "transition":
                "energy" (float): Band gap energy.
                "direct" (bool): Whether the gap is direct.
                "transition" (str): Kpoint labels of the transition (e.g., "\\Gamma-X").
        """
        if self.is_metal():
            return {"energy": 0.0, "direct": False, "transition": None}

        cbm = self.get_cbm()
        vbm = self.get_vbm()
        result = {
            "direct": False,
            "transition": None,
            "energy": cbm["energy"] - vbm["energy"],
        }

        if (cbm["kpoint"].label is not None and cbm["kpoint"].label == vbm["kpoint"].label) or np.linalg.norm(
            cbm["kpoint"].cart_coords - vbm["kpoint"].cart_coords
        ) < 0.01:
            result["direct"] = True

        result["transition"] = "-".join(
            [
                (str(c.label) if c.label is not None else f"({','.join(f'{c.frac_coords[i]:.3f}' for i in range(3))})")
                for c in [vbm["kpoint"], cbm["kpoint"]]
            ]
        )

        return result

    def get_direct_band_gap_dict(self) -> dict[Spin, dict[str, Any]]:
        """Get information about the direct band gap.

        Returns:
            dict[Spin, dict[str, Any]]: The band gaps indexed by spin
                along with their band indices and kpoint index.
        """
        if self.is_metal():
            raise ValueError("get_direct_band_gap_dict should only be used with non-metals")

        direct_gap_dict = {}
        for spin, v in self.bands.items():
            above = v[np.all(v > self.efermi, axis=1)]
            min_above = np.min(above, axis=0)
            below = v[np.all(v < self.efermi, axis=1)]
            max_below = np.max(below, axis=0)
            diff = min_above - max_below
            kpoint_index = np.argmin(diff)
            band_indices = [
                np.argmax(below[:, kpoint_index]),
                np.argmin(above[:, kpoint_index]) + len(below),
            ]
            direct_gap_dict[spin] = {
                "value": diff[kpoint_index],
                "kpoint_index": kpoint_index,
                "band_indices": band_indices,
            }
        return direct_gap_dict

    def get_direct_band_gap(self) -> float:
        """Get the direct band gap.

        Returns:
            float: The direct band gap value.
        """
        if self.is_metal():
            return 0.0

        dg = self.get_direct_band_gap_dict()
        return min(v["value"] for v in dg.values())

    def get_sym_eq_kpoints(
        self,
        kpoint: NDArray,
        cartesian: bool = False,
        tol: float = 1e-2,
    ) -> NDArray | None:
        """Get unique symmetrically equivalent Kpoints.

        Args:
            kpoint (1x3 array): Coordinate of the Kpoint.
            cartesian (bool): Whether kpoint is in Cartesian or fractional coordinates.
            tol (float): Tolerance below which coordinates are considered equal.

        Returns:
            (1x3 NDArray) | None: None if structure is not available.
        """
        if not self.structure:
            return None

        sg = SpacegroupAnalyzer(self.structure)
        symm_ops = sg.get_point_group_operations(cartesian=cartesian)
        points = np.dot(kpoint, [m.rotation_matrix for m in symm_ops])
        rm_list = []
        # Identify and remove duplicates from equivalent k-points
        for i in range(len(points) - 1):
            for j in range(i + 1, len(points)):
                if np.allclose(pbc_diff(points[i], points[j]), [0, 0, 0], tol):
                    rm_list.append(i)
                    break
        return np.delete(points, rm_list, axis=0)

    def get_kpoint_degeneracy(
        self,
        kpoint: NDArray,
        cartesian: bool = False,
        tol: float = 1e-2,
    ) -> int | None:
        """Get degeneracy of a given kpoint based on structure symmetry.

        Args:
            kpoint (1x3 NDArray): Coordinate of the k-point.
            cartesian (bool): Whether kpoint is in Cartesian or fractional coordinates.
            tol (float): Tolerance below which coordinates are considered equal.

        Returns:
            int | None: Degeneracy, or None if structure is not available.
        """
        all_kpts = self.get_sym_eq_kpoints(kpoint, cartesian, tol=tol)
        return len(all_kpts) if all_kpts is not None else None

    def as_dict(self) -> dict[str, Any]:
        """JSON-serializable dict representation of BandStructure."""
        dct: dict[str, Any] = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "lattice_rec": self.lattice_rec.as_dict(),
            "efermi": self.efermi,
            "kpoints": [],
        }
        # kpoints are not kpoint objects dicts but are frac coords (this makes
        # the dict smaller and avoids the repetition of the lattice).
        for k in self.kpoints:
            dct["kpoints"].append(k.as_dict()["fcoords"])

        dct["bands"] = {str(int(spin)): self.bands[spin].tolist() for spin in self.bands}
        dct["is_metal"] = self.is_metal()
        vbm = self.get_vbm()
        dct["vbm"] = {
            "energy": vbm["energy"],
            "kpoint_index": vbm["kpoint_index"],
            "band_index": {str(int(spin)): vbm["band_index"][spin] for spin in vbm["band_index"]},
            "projections": {str(spin): v.tolist() for spin, v in vbm["projections"].items()},
        }
        cbm = self.get_cbm()
        dct["cbm"] = {
            "energy": cbm["energy"],
            "kpoint_index": cbm["kpoint_index"],
            "band_index": {str(int(spin)): cbm["band_index"][spin] for spin in cbm["band_index"]},
            "projections": {str(spin): v.tolist() for spin, v in cbm["projections"].items()},
        }
        dct["band_gap"] = self.get_band_gap()
        dct["labels_dict"] = {}
        dct["is_spin_polarized"] = self.is_spin_polarized

        # MongoDB does not accept keys starting with "$", add a space to fix this.
        for c, label in self.labels_dict.items():
            mongo_key = f" {c}" if c.startswith("$") else c
            dct["labels_dict"][mongo_key] = label.as_dict()["fcoords"]
        dct["projections"] = {}
        if len(self.projections) != 0 and self.structure is not None:
            dct["structure"] = self.structure.as_dict()
            dct["projections"] = {str(int(spin)): np.array(v).tolist() for spin, v in self.projections.items()}
        return dct

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        """Create from a dict.

        Args:
            dct: A dict with all data for a BandStructure.

        Returns:
            BandStructure
        """
        # Strip the label to recover initial string
        # (see trick used in as_dict to handle "$"" chars)
        labels_dict = {k.strip(): v for k, v in dct["labels_dict"].items()}

        if isinstance(next(iter(dct["bands"].values())), dict):
            eigenvals = {Spin(int(k)): np.array(dct["bands"][k]["data"]) for k in dct["bands"]}
        else:
            eigenvals = {Spin(int(k)): dct["bands"][k] for k in dct["bands"]}

        structure = None
        if "structure" in dct:
            structure = Structure.from_dict(dct["structure"])

        projections = {}
        try:
            if dct.get("projections"):
                if isinstance(dct["projections"]["1"][0][0], dict):
                    raise ValueError("Old band structure dict format detected!")
                projections = {Spin(int(spin)): np.array(v) for spin, v in dct["projections"].items()}

            return cls(
                dct["kpoints"],
                eigenvals,
                Lattice(dct["lattice_rec"]["matrix"]),
                dct["efermi"],
                labels_dict,
                structure=structure,
                projections=projections,
            )

        except Exception:
            warnings.warn(
                "Trying from_dict failed. Now we are trying the old "
                "format. Please convert your BS dicts to the new "
                "format. The old format will be retired in pymatgen "
                "5.0.",
                stacklevel=2,
            )
            return cls.from_old_dict(dct)

    @classmethod
    def from_old_dict(cls, dct: dict[str, Any]) -> Self:
        """
        Args:
            dct (dict): A dict with all data for a BandStructure object.

        Returns:
            BandStructure
        """
        # Strip the label to recover initial string
        # (see trick used in as_dict to handle "$" chars)
        labels_dict = {k.strip(): v for k, v in dct["labels_dict"].items()}
        projections: dict = {}
        structure = None
        if dct.get("projections"):
            structure = Structure.from_dict(dct["structure"])
            projections = {}
            for spin in dct["projections"]:
                dd = []
                for ii in range(len(dct["projections"][spin])):
                    ddd = []
                    for jj in range(len(dct["projections"][spin][ii])):
                        dddd = []
                        for kk in range(len(dct["projections"][spin][ii][jj])):
                            ddddd = []
                            orb = Orbital(kk).name
                            for ll in range(len(dct["projections"][spin][ii][jj][orb])):
                                ddddd.append(dct["projections"][spin][ii][jj][orb][ll])
                            dddd.append(np.array(ddddd))
                        ddd.append(np.array(dddd))
                    dd.append(np.array(ddd))
                projections[Spin(int(spin))] = np.array(dd)

        return cls(
            dct["kpoints"],
            {Spin(int(k)): dct["bands"][k] for k in dct["bands"]},
            Lattice(dct["lattice_rec"]["matrix"]),
            dct["efermi"],
            labels_dict,
            structure=structure,
            projections=projections,
        )


class BandStructureSymmLine(BandStructure, MSONable):
    r"""Store band structures along selected (symmetry) lines in the Brillouin zone.
    We call the different symmetry lines (ex: \\Gamma to Z) "branches".
    """

    def __init__(
        self,
        kpoints: ArrayLike,
        eigenvals: Mapping[Spin, ArrayLike],
        lattice: Lattice,
        efermi: float,
        labels_dict: Mapping[str, Kpoint],
        coords_are_cartesian: bool = False,
        structure: Structure | None = None,
        projections: Mapping[Spin, NDArray] | None = None,
    ) -> None:
        """
        Args:
            kpoints (NDArray): Array of kpoint, in frac_coords of the
                given lattice by default
            eigenvals (dict[Spin, list]): Energies for spin up and spin down
                {Spin.up:[][],Spin.down:[][]}, the first index of the array
                [][] refers to the band and the second to the index of the
                kpoint. The kpoints are ordered according to the order of the
                kpoints array. If the band structure is not spin polarized, we
                only store one data set under Spin.up.
            lattice (Lattice): The reciprocal lattice. Pymatgen uses the physics
                convention of reciprocal lattice vectors with a 2*pi coefficient.
            efermi (float): The Fermi level.
            labels_dict (dict[str, Kpoint]): Dict mapping label to Kpoint.
            coords_are_cartesian (bool): Whether coordinates are cartesian.
            structure (Structure): The crystal structure associated with the
                band structure. This is needed if we provide projections to
                the band structure.
            projections (dict[Spin, NDArray]): Orbital projections as {spin: array}.
                The indices of the array are [band_index, kpoint_index, orbital_index,
                ion_index].If the band structure is not spin polarized, we only
                store one data set under Spin.up.
        """
        super().__init__(
            kpoints,
            eigenvals,
            lattice,
            efermi,
            labels_dict,
            coords_are_cartesian,
            structure,
            projections,
        )
        self.distance = []
        self.branches = []
        one_group: list = []
        branches_tmp = []
        # Get labels and distance for each kpoint
        previous_kpoint = self.kpoints[0]
        previous_distance = 0.0

        previous_label = self.kpoints[0].label
        for i, kpt in enumerate(self.kpoints):
            label = kpt.label
            if label is not None and previous_label is not None:
                self.distance.append(previous_distance)
            else:
                self.distance.append(
                    np.linalg.norm(kpt.cart_coords - previous_kpoint.cart_coords)  # type: ignore[arg-type]
                    + previous_distance
                )
            previous_kpoint = kpt
            previous_distance = self.distance[i]
            if label and previous_label:
                if len(one_group) != 0:
                    branches_tmp.append(one_group)
                one_group = []
            previous_label = label
            one_group.append(i)

        if len(one_group) != 0:
            branches_tmp.append(one_group)
        for branch in branches_tmp:
            self.branches.append(
                {
                    "start_index": branch[0],
                    "end_index": branch[-1],
                    "name": f"{self.kpoints[branch[0]].label}-{self.kpoints[branch[-1]].label}",
                }
            )

        self.is_spin_polarized = False
        if len(self.bands) == 2:
            self.is_spin_polarized = True

    def get_equivalent_kpoints(self, index: int) -> list[int]:
        """Get kpoint indices equivalent (having the same coords) to the given one.

        Args:
            index (int): The kpoint index

        Returns:
            list[int]: Equivalent indices.

        TODO: now it uses the label, we might want to use coordinates
        instead in case there was a mislabel.
        """
        # If the kpoint has no label it can't have a repetition
        # along the BandStructureSymmLine object
        if self.kpoints[index].label is None:
            return [index]

        list_index_kpoints = []
        for idx, kpt in enumerate(self.kpoints):
            if kpt.label == self.kpoints[index].label:
                list_index_kpoints.append(idx)

        return list_index_kpoints

    def get_branch(self, index: int) -> list[dict[str, Any]]:
        """Get what branch(es) is the kpoint. It takes into account the
            fact that one kpoint (e.g., Gamma) can be in several branches.

        Args:
            index (int): The kpoint index.

        Returns:
            A list of dicts [{"name", "start_index", "end_index", "index"}]
                indicating all branches in which the k_point is.
        """
        to_return = []
        for idx in self.get_equivalent_kpoints(index):
            for branch in self.branches:
                if branch["start_index"] <= idx <= branch["end_index"]:
                    to_return.append(
                        {
                            "name": branch["name"],
                            "start_index": branch["start_index"],
                            "end_index": branch["end_index"],
                            "index": idx,
                        }
                    )
        return to_return

    def apply_scissor(self, new_band_gap: float) -> Self:
        """Apply a scissor operator (shift of the CBM) to fit the given band gap.
        If it's a metal, we look for the band crossing the Fermi level
        and shift this one up. This will not work all the time for metals!

        Args:
            new_band_gap (float): The band gap the scissor band structure need to have.

        Returns:
            BandStructureSymmLine: With the applied scissor shift.
        """
        if self.is_metal():
            # Move then the highest index band crossing the Fermi level find this band...
            max_index = -1000
            # spin_index = None
            for idx in range(self.nb_bands):
                below = above = False
                for j in range(len(self.kpoints)):
                    if self.bands[Spin.up][idx][j] < self.efermi:
                        below = True
                    if self.bands[Spin.up][idx][j] > self.efermi:
                        above = True
                if above and below and idx > max_index:
                    max_index = idx
                    # spin_index = Spin.up
                if self.is_spin_polarized:
                    below = above = False
                    for j in range(len(self.kpoints)):
                        if self.bands[Spin.down][idx][j] < self.efermi:
                            below = True
                        if self.bands[Spin.down][idx][j] > self.efermi:
                            above = True
                    if above and below and idx > max_index:
                        max_index = idx
                        # spin_index = Spin.down
            old_dict = self.as_dict()
            shift = new_band_gap
            for spin in old_dict["bands"]:
                for k in range(len(old_dict["bands"][spin])):
                    for v in range(len(old_dict["bands"][spin][k])):
                        if k >= max_index:
                            old_dict["bands"][spin][k][v] += shift

        else:
            shift = new_band_gap - self.get_band_gap()["energy"]
            old_dict = self.as_dict()
            for spin in old_dict["bands"]:
                for k in range(len(old_dict["bands"][spin])):
                    for v in range(len(old_dict["bands"][spin][k])):
                        if old_dict["bands"][spin][k][v] >= old_dict["cbm"]["energy"]:
                            old_dict["bands"][spin][k][v] += shift
            old_dict["efermi"] += shift

        return self.from_dict(old_dict)

    def as_dict(self) -> dict[str, Any]:
        """JSON-serializable dict representation of BandStructureSymmLine."""
        dct = super().as_dict()
        dct["branches"] = self.branches
        return dct


class LobsterBandStructureSymmLine(BandStructureSymmLine):
    """LOBSTER subclass of BandStructure with customized functions."""

    def as_dict(self) -> dict[str, Any]:
        """JSON-serializable dict representation of BandStructureSymmLine."""
        dct: dict[str, Any] = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "lattice_rec": self.lattice_rec.as_dict(),
            "efermi": self.efermi,
            "kpoints": [],
        }
        # kpoints are not kpoint objects dicts but are frac coords (this makes
        # the dict smaller and avoids the repetition of the lattice
        for k in self.kpoints:
            dct["kpoints"].append(k.as_dict()["fcoords"])
        dct["branches"] = self.branches
        dct["bands"] = {str(int(spin)): self.bands[spin].tolist() for spin in self.bands}
        dct["is_metal"] = self.is_metal()
        vbm = self.get_vbm()
        dct["vbm"] = {
            "energy": vbm["energy"],
            "kpoint_index": [int(x) for x in vbm["kpoint_index"]],
            "band_index": {str(int(spin)): vbm["band_index"][spin] for spin in vbm["band_index"]},
            "projections": {str(spin): v for spin, v in vbm["projections"].items()},
        }
        cbm = self.get_cbm()
        dct["cbm"] = {
            "energy": cbm["energy"],
            "kpoint_index": [int(x) for x in cbm["kpoint_index"]],
            "band_index": {str(int(spin)): cbm["band_index"][spin] for spin in cbm["band_index"]},
            "projections": {str(spin): v for spin, v in cbm["projections"].items()},
        }
        dct["band_gap"] = self.get_band_gap()
        dct["labels_dict"] = {}
        dct["is_spin_polarized"] = self.is_spin_polarized

        # MongoDB does not accept keys starting with "$", add a space to fix this.
        for c, label in self.labels_dict.items():
            mongo_key = f" {c}" if c.startswith("$") else c
            dct["labels_dict"][mongo_key] = label.as_dict()["fcoords"]
        if len(self.projections) != 0 and self.structure is not None:
            dct["structure"] = self.structure.as_dict()
            dct["projections"] = {str(int(spin)): np.array(v).tolist() for spin, v in self.projections.items()}
        return dct

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        """
        Args:
            dct (dict): All data for a LobsterBandStructureSymmLine object.

        Returns:
            A LobsterBandStructureSymmLine object.
        """
        try:
            # Strip the label to recover initial string
            # (see trick used in as_dict to handle "$" chars)
            labels_dict = {k.strip(): v for k, v in dct["labels_dict"].items()}
            projections = {}
            structure = None
            if dct.get("projections"):
                if isinstance(dct["projections"]["1"][0][0], dict):
                    raise ValueError("Old band structure dict format detected!")
                structure = Structure.from_dict(dct["structure"])
                projections = {Spin(int(spin)): np.array(v) for spin, v in dct["projections"].items()}

            return cls(
                dct["kpoints"],
                {Spin(int(k)): dct["bands"][k] for k in dct["bands"]},
                Lattice(dct["lattice_rec"]["matrix"]),
                dct["efermi"],
                labels_dict,
                structure=structure,
                projections=projections,
            )
        except Exception:
            warnings.warn(
                "Trying from_dict failed. Now we are trying the old "
                "format. Please convert your BS dicts to the new "
                "format. The old format will be retired in pymatgen "
                "5.0.",
                stacklevel=2,
            )
            return cls.from_old_dict(dct)

    @classmethod
    def from_old_dict(cls, dct: dict[str, Any]) -> Self:
        """
        Args:
            dct (dict): All data for a LobsterBandStructureSymmLine object.

        Returns:
            A LobsterBandStructureSymmLine object
        """
        # Strip the label to recover initial string
        # (see trick used in as_dict to handle "$" chars)
        labels_dict = {k.strip(): v for k, v in dct["labels_dict"].items()}
        projections: dict = {}
        structure = None
        if "projections" in dct and len(dct["projections"]) != 0:
            structure = Structure.from_dict(dct["structure"])
            projections = {}
            for spin in dct["projections"]:
                dd = []
                for i in range(len(dct["projections"][spin])):
                    ddd = []
                    for j in range(len(dct["projections"][spin][i])):
                        ddd.append(dct["projections"][spin][i][j])
                    dd.append(np.array(ddd))
                projections[Spin(int(spin))] = np.array(dd)

        return cls(
            dct["kpoints"],
            {Spin(int(k)): dct["bands"][k] for k in dct["bands"]},
            Lattice(dct["lattice_rec"]["matrix"]),
            dct["efermi"],
            labels_dict,
            structure=structure,
            projections=projections,
        )

    def get_projection_on_elements(self) -> dict[Spin, list]:
        """Get projections on elements. It sums over all available orbitals
        for each element.

        Returns:
            dict[Spin, list]: dict in the {Spin.up:[][{Element:values}],
                Spin.down:[][{Element:values}]} format.
                If there is no projections in the band structure, return {}.
        """
        result: dict[Spin, list] = {}
        for spin, v in self.projections.items():
            result[spin] = [[defaultdict(float) for _ in range(len(self.kpoints))] for _ in range(self.nb_bands)]
            for i, j in itertools.product(range(self.nb_bands), range(len(self.kpoints))):
                for key, item in v[i][j].items():
                    for item2 in item.values():
                        specie = str(Element(re.split(r"[0-9]+", key)[0]))
                        result[spin][i][j][specie] += item2
        return result

    def get_projections_on_elements_and_orbitals(
        self,
        el_orb_spec: dict[SpeciesLike, list],
    ) -> dict[Spin, list]:
        """Get projections on elements and specific orbitals.

        Args:
            el_orb_spec (dict): Elements and Orbitals for which we want
                to project on. It is given as {Element: [orbitals]},
                e.g. {"Si": ["3s", "3p"]} or {"Si": ["3s", "3p_x", "3p_y", "3p_z']}
                depending on input files.

        Returns:
            A dictionary of projections on elements in the
            {Spin.up:[][{Element:{orb:values}}],
            Spin.down:[][{Element:{orb:values}}]} format
            if there is no projections in the band structure returns an empty
            dict.
        """
        result: dict[Spin, list] = {}
        el_orb_spec = {get_el_sp(el): orbs for el, orbs in el_orb_spec.items()}
        for spin, v in self.projections.items():
            result[spin] = [
                [{str(e): defaultdict(float) for e in el_orb_spec} for _ in range(len(self.kpoints))]
                for _ in range(self.nb_bands)
            ]

            for i, j in itertools.product(range(self.nb_bands), range(len(self.kpoints))):
                for key, item in v[i][j].items():
                    for key2, item2 in item.items():
                        specie = str(Element(re.split(r"[0-9]+", key)[0]))
                        if get_el_sp(str(specie)) in el_orb_spec and key2 in el_orb_spec[get_el_sp(str(specie))]:
                            result[spin][i][j][specie][key2] += item2
        return result


@overload
def get_reconstructed_band_structure(
    list_bs: list[BandStructure],
    efermi: float | None = None,
) -> BandStructure:
    pass


@overload
def get_reconstructed_band_structure(
    list_bs: list[BandStructureSymmLine],
    efermi: float | None = None,
) -> BandStructureSymmLine:
    pass


def get_reconstructed_band_structure(
    list_bs: list[BandStructure] | list[BandStructureSymmLine],
    efermi: float | None = None,
) -> BandStructure | BandStructureSymmLine:
    """Merge multiple BandStructure(SymmLine) objects to a single one.

    This is typically useful when you split non self-consistent band
    structure runs to several independent jobs and want to merge the results.

    Args:
        list_bs (list): BandStructure or BandStructureSymmLine objects.
        efermi (float): The Fermi level of the reconstructed band structure.
            If None, an average of all the Fermi levels in each
            object in the list_bs is used.

    Returns:
        A BandStructure or BandStructureSymmLine object (depending on
        the type of the objects in list_bs).
    """
    if efermi is None:
        efermi = sum(b.efermi for b in list_bs) / len(list_bs)

    rec_lattice = list_bs[0].lattice_rec
    nb_bands = min(list_bs[i].nb_bands for i in range(len(list_bs)))

    kpoints = np.concatenate([[kpt.frac_coords for kpt in bs.kpoints] for bs in list_bs])
    dicts = [bs.labels_dict for bs in list_bs]
    labels_dict = {key: val.frac_coords for dct in dicts for key, val in dct.items()}

    eigenvals = {Spin.up: np.concatenate([bs.bands[Spin.up][:nb_bands] for bs in list_bs], axis=1)}
    if list_bs[0].is_spin_polarized:
        eigenvals[Spin.down] = np.concatenate([bs.bands[Spin.down][:nb_bands] for bs in list_bs], axis=1)

    projections = {}
    if len(list_bs[0].projections) != 0:
        projs = [bs.projections[Spin.up][:nb_bands] for bs in list_bs]
        projections[Spin.up] = np.concatenate(projs, axis=1)

        if list_bs[0].is_spin_polarized:
            projs = [bs.projections[Spin.down][:nb_bands] for bs in list_bs]
            projections[Spin.down] = np.concatenate(projs, axis=1)

    if isinstance(list_bs[0], BandStructureSymmLine):
        return BandStructureSymmLine(
            kpoints,
            eigenvals,
            rec_lattice,
            efermi,
            labels_dict,  # type:ignore[arg-type]
            structure=list_bs[0].structure,
            projections=projections,
        )
    return BandStructure(
        kpoints,
        eigenvals,
        rec_lattice,
        efermi,
        labels_dict,  # type:ignore[arg-type]
        structure=list_bs[0].structure,
        projections=projections,
    )
