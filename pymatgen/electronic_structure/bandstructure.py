# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes to define everything related to band structures.
"""

from __future__ import annotations

import collections
import itertools
import math
import re
import warnings

import numpy as np
from monty.json import MSONable

from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element, get_el_sp
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Orbital, Spin
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord import pbc_diff

__author__ = "Geoffroy Hautier, Shyue Ping Ong, Michael Kocher"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "March 14, 2012"


class Kpoint(MSONable):
    """
    Class to store kpoint objects. A kpoint is defined with a lattice and frac
    or Cartesian coordinates syntax similar than the site object in
    pymatgen.core.structure.
    """

    def __init__(
        self,
        coords,
        lattice,
        to_unit_cell=False,
        coords_are_cartesian=False,
        label=None,
    ):
        """
        Args:
            coords: coordinate of the kpoint as a numpy array
            lattice: A pymatgen.core.lattice.Lattice lattice object representing
                the reciprocal lattice of the kpoint
            to_unit_cell: Translates fractional coordinate to the basic unit
                cell, i.e., all fractional coordinates satisfy 0 <= a < 1.
                Defaults to False.
            coords_are_cartesian: Boolean indicating if the coordinates given are
                in Cartesian or fractional coordinates (by default fractional)
            label: the label of the kpoint if any (None by default)
        """
        self._lattice = lattice
        self._fcoords = lattice.get_fractional_coords(coords) if coords_are_cartesian else coords
        self._label = label

        if to_unit_cell:
            for i, fc in enumerate(self._fcoords):
                self._fcoords[i] -= math.floor(fc)

        self._ccoords = lattice.get_cartesian_coords(self._fcoords)

    @property
    def lattice(self):
        """
        The lattice associated with the kpoint. It's a
        pymatgen.core.lattice.Lattice object
        """
        return self._lattice

    @property
    def label(self):
        """
        The label associated with the kpoint
        """
        return self._label

    @property
    def frac_coords(self):
        """
        The fractional coordinates of the kpoint as a numpy array
        """
        return np.copy(self._fcoords)

    @property
    def cart_coords(self):
        """
        The Cartesian coordinates of the kpoint as a numpy array
        """
        return np.copy(self._ccoords)

    @property
    def a(self):
        """
        Fractional a coordinate of the kpoint
        """
        return self._fcoords[0]

    @property
    def b(self):
        """
        Fractional b coordinate of the kpoint
        """
        return self._fcoords[1]

    @property
    def c(self):
        """
        Fractional c coordinate of the kpoint
        """
        return self._fcoords[2]

    def __str__(self):
        """
        Returns a string with fractional, Cartesian coordinates and label
        """
        return f"{self.frac_coords} {self.cart_coords} {self.label}"

    def as_dict(self):
        """
        JSON-serializable dict representation of a kpoint
        """
        return {
            "lattice": self.lattice.as_dict(),
            "fcoords": self.frac_coords.tolist(),
            "ccoords": self.cart_coords.tolist(),
            "label": self.label,
            "@module": type(self).__module__,
            "@class": type(self).__name__,
        }

    @classmethod
    def from_dict(cls, d):
        """
        Create from dict.

        Args:
            A dict with all data for a kpoint object.

        Returns:
            A Kpoint object
        """
        return cls(
            coords=d["fcoords"],
            lattice=Lattice.from_dict(d["lattice"]),
            coords_are_cartesian=False,
            label=d["label"],
        )


class BandStructure:
    """
    This is the most generic band structure data possible
    it's defined by a list of kpoints + energies for each of them

    .. attribute:: kpoints:
        the list of kpoints (as Kpoint objects) in the band structure

    .. attribute:: lattice_rec

        the reciprocal lattice of the band structure.

    .. attribute:: efermi

        the fermi energy

    .. attribute::  is_spin_polarized

        True if the band structure is spin-polarized, False otherwise

    .. attribute:: bands

        The energy eigenvalues as a {spin: ndarray}. Note that the use of an
        ndarray is necessary for computational as well as memory efficiency
        due to the large amount of numerical data. The indices of the ndarray
        are [band_index, kpoint_index].

    .. attribute:: nb_bands

        returns the number of bands in the band structure

    .. attribute:: structure

        returns the structure

    .. attribute:: projections

        The projections as a {spin: ndarray}. Note that the use of an
        ndarray is necessary for computational as well as memory efficiency
        due to the large amount of numerical data. The indices of the ndarray
        are [band_index, kpoint_index, orbital_index, ion_index].
    """

    def __init__(
        self,
        kpoints: np.ndarray,
        eigenvals: dict[Spin, np.ndarray],
        lattice: Lattice,
        efermi: float,
        labels_dict=None,
        coords_are_cartesian: bool = False,
        structure: Structure | None = None,
        projections: dict[Spin, np.ndarray] | None = None,
    ) -> None:
        """
        Args:
            kpoints: list of kpoint as numpy arrays, in frac_coords of the
                given lattice by default
            eigenvals: dict of energies for spin up and spin down
                {Spin.up:[][],Spin.down:[][]}, the first index of the array
                [][] refers to the band and the second to the index of the
                kpoint. The kpoints are ordered according to the order of the
                kpoints array. If the band structure is not spin polarized, we
                only store one data set under Spin.up
            lattice: The reciprocal lattice as a pymatgen Lattice object.
                Pymatgen uses the physics convention of reciprocal lattice vectors
                WITH a 2*pi coefficient
            efermi (float): fermi energy
            labels_dict: (dict) of {} this links a kpoint (in frac coords or
                Cartesian coordinates depending on the coords) to a label.
            coords_are_cartesian: Whether coordinates are cartesian.
            structure: The crystal structure (as a pymatgen Structure object)
                associated with the band structure. This is needed if we
                provide projections to the band structure
            projections: dict of orbital projections as {spin: ndarray}. The
                indices of the ndarrayare [band_index, kpoint_index, orbital_index,
                ion_index].If the band structure is not spin polarized, we only
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

        if len(self.projections) != 0 and self.structure is None:
            raise Exception("if projections are provided a structure object needs also to be given")

        for k in kpoints:
            # let see if this kpoint has been assigned a label
            label = None
            for c in labels_dict:
                if np.linalg.norm(k - np.array(labels_dict[c])) < 0.0001:
                    label = c
                    self.labels_dict[label] = Kpoint(
                        k,
                        lattice,
                        label=label,
                        coords_are_cartesian=coords_are_cartesian,
                    )
            self.kpoints.append(Kpoint(k, lattice, label=label, coords_are_cartesian=coords_are_cartesian))
        self.bands = {spin: np.array(v) for spin, v in eigenvals.items()}
        self.nb_bands = len(eigenvals[Spin.up])
        self.is_spin_polarized = len(self.bands) == 2

    def get_projection_on_elements(self):
        """
        Method returning a dictionary of projections on elements.

        Returns:
            a dictionary in the {Spin.up:[][{Element:values}],
            Spin.down:[][{Element:values}]} format
            if there is no projections in the band structure
            returns an empty dict
        """
        result = {}
        structure = self.structure
        for spin, v in self.projections.items():
            result[spin] = [
                [collections.defaultdict(float) for i in range(len(self.kpoints))] for j in range(self.nb_bands)
            ]
            for i, j, k in itertools.product(
                range(self.nb_bands),
                range(len(self.kpoints)),
                range(structure.num_sites),
            ):
                result[spin][i][j][str(structure[k].specie)] += np.sum(v[i, j, :, k])
        return result

    def get_projections_on_elements_and_orbitals(self, el_orb_spec):
        """
        Method returning a dictionary of projections on elements and specific
        orbitals

        Args:
            el_orb_spec: A dictionary of Elements and Orbitals for which we want
                to have projections on. It is given as: {Element:[orbitals]},
                e.g., {'Cu':['d','s']}

        Returns:
            A dictionary of projections on elements in the
            {Spin.up:[][{Element:{orb:values}}],
            Spin.down:[][{Element:{orb:values}}]} format
            if there is no projections in the band structure returns an empty
            dict.
        """
        result = {}
        structure = self.structure
        el_orb_spec = {get_el_sp(el): orbs for el, orbs in el_orb_spec.items()}
        for spin, v in self.projections.items():
            result[spin] = [
                [{str(e): collections.defaultdict(float) for e in el_orb_spec} for i in range(len(self.kpoints))]
                for j in range(self.nb_bands)
            ]

            for i, j, k in itertools.product(
                range(self.nb_bands),
                range(len(self.kpoints)),
                range(structure.num_sites),
            ):
                sp = structure[k].specie
                for orb_i in range(len(v[i][j])):
                    o = Orbital(orb_i).name[0]
                    if sp in el_orb_spec:
                        if o in el_orb_spec[sp]:
                            result[spin][i][j][str(sp)][o] += v[i][j][orb_i][k]
        return result

    def is_metal(self, efermi_tol=1e-4) -> bool:
        """
        Check if the band structure indicates a metal by looking if the fermi
        level crosses a band.

        Returns:
            True if a metal, False if not
        """
        for values in self.bands.values():
            for i in range(self.nb_bands):
                if np.any(values[i, :] - self.efermi < -efermi_tol) and np.any(values[i, :] - self.efermi > efermi_tol):
                    return True
        return False

    def get_vbm(self):
        """
        Returns data about the VBM.

        Returns:
            dict as {"band_index","kpoint_index","kpoint","energy"}
            - "band_index": A dict with spin keys pointing to a list of the
            indices of the band containing the VBM (please note that you
            can have several bands sharing the VBM) {Spin.up:[],
            Spin.down:[]}
            - "kpoint_index": The list of indices in self.kpoints for the
            kpoint VBM. Please note that there can be several
            kpoint_indices relating to the same kpoint (e.g., Gamma can
            occur at different spots in the band structure line plot)
            - "kpoint": The kpoint (as a kpoint object)
            - "energy": The energy of the VBM
            - "projections": The projections along sites and orbitals of the
            VBM if any projection data is available (else it is an empty
            dictionary). The format is similar to the projections field in
            BandStructure: {spin:{'Orbital': [proj]}} where the array
            [proj] is ordered according to the sites in structure
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
        index = None
        kpointvbm = None
        for value in self.bands.values():
            for i, j in zip(*np.where(value < self.efermi)):
                if value[i, j] > max_tmp:
                    max_tmp = float(value[i, j])
                    index = j
                    kpointvbm = self.kpoints[j]

        list_ind_kpts = []
        if kpointvbm.label is not None:
            for i, kpt in enumerate(self.kpoints):
                if kpt.label == kpointvbm.label:
                    list_ind_kpts.append(i)
        else:
            list_ind_kpts.append(index)
        # get all other bands sharing the vbm
        list_ind_band = collections.defaultdict(list)
        for spin in self.bands:
            for i in range(self.nb_bands):
                if math.fabs(self.bands[spin][i][index] - max_tmp) < 0.001:
                    list_ind_band[spin].append(i)
        proj = {}
        for spin, value in self.projections.items():
            if len(list_ind_band[spin]) == 0:
                continue
            proj[spin] = value[list_ind_band[spin][0]][list_ind_kpts[0]]
        return {
            "band_index": list_ind_band,
            "kpoint_index": list_ind_kpts,
            "kpoint": kpointvbm,
            "energy": max_tmp,
            "projections": proj,
        }

    def get_cbm(self):
        """
        Returns data about the CBM.

        Returns:
            {"band_index","kpoint_index","kpoint","energy"}
            - "band_index": A dict with spin keys pointing to a list of the
            indices of the band containing the CBM (please note that you
            can have several bands sharing the CBM) {Spin.up:[],
            Spin.down:[]}
            - "kpoint_index": The list of indices in self.kpoints for the
            kpoint CBM. Please note that there can be several
            kpoint_indices relating to the same kpoint (e.g., Gamma can
            occur at different spots in the band structure line plot)
            - "kpoint": The kpoint (as a kpoint object)
            - "energy": The energy of the CBM
            - "projections": The projections along sites and orbitals of the
            CBM if any projection data is available (else it is an empty
            dictionary). The format is similar to the projections field in
            BandStructure: {spin:{'Orbital': [proj]}} where the array
            [proj] is ordered according to the sites in structure
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

        index = None
        kpointcbm = None
        for value in self.bands.values():
            for i, j in zip(*np.where(value >= self.efermi)):
                if value[i, j] < max_tmp:
                    max_tmp = float(value[i, j])
                    index = j
                    kpointcbm = self.kpoints[j]

        list_index_kpoints = []
        if kpointcbm.label is not None:
            for i, kpt in enumerate(self.kpoints):
                if kpt.label == kpointcbm.label:
                    list_index_kpoints.append(i)
        else:
            list_index_kpoints.append(index)

        # get all other bands sharing the cbm
        list_index_band = collections.defaultdict(list)
        for spin in self.bands:
            for i in range(self.nb_bands):
                if math.fabs(self.bands[spin][i][index] - max_tmp) < 0.001:
                    list_index_band[spin].append(i)
        proj = {}
        for spin, value in self.projections.items():
            if len(list_index_band[spin]) == 0:
                continue
            proj[spin] = value[list_index_band[spin][0]][list_index_kpoints[0]]

        return {
            "band_index": list_index_band,
            "kpoint_index": list_index_kpoints,
            "kpoint": kpointcbm,
            "energy": max_tmp,
            "projections": proj,
        }

    def get_band_gap(self):
        r"""
        Returns band gap data.

        Returns:
            A dict {"energy","direct","transition"}:
            "energy": band gap energy
            "direct": A boolean telling if the gap is direct or not
            "transition": kpoint labels of the transition (e.g., "\\Gamma-X")
        """
        if self.is_metal():
            return {"energy": 0.0, "direct": False, "transition": None}
        cbm = self.get_cbm()
        vbm = self.get_vbm()
        result = dict(direct=False, energy=0.0, transition=None)

        result["energy"] = cbm["energy"] - vbm["energy"]

        if (cbm["kpoint"].label is not None and cbm["kpoint"].label == vbm["kpoint"].label) or np.linalg.norm(
            cbm["kpoint"].cart_coords - vbm["kpoint"].cart_coords
        ) < 0.01:
            result["direct"] = True

        result["transition"] = "-".join(
            [
                str(c.label)
                if c.label is not None
                else "(" + ",".join(f"{c.frac_coords[i]:.3f}" for i in range(3)) + ")"
                for c in [vbm["kpoint"], cbm["kpoint"]]
            ]
        )

        return result

    def get_direct_band_gap_dict(self):
        """
        Returns a dictionary of information about the direct
        band gap

        Returns:
            a dictionary of the band gaps indexed by spin
            along with their band indices and k-point index
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

    def get_direct_band_gap(self):
        """
        Returns the direct band gap.

        Returns:
             the value of the direct band gap
        """
        if self.is_metal():
            return 0.0
        dg = self.get_direct_band_gap_dict()
        return min(v["value"] for v in dg.values())

    def get_sym_eq_kpoints(self, kpoint, cartesian=False, tol: float = 1e-2):
        """
        Returns a list of unique symmetrically equivalent k-points.

        Args:
            kpoint (1x3 array): coordinate of the k-point
            cartesian (bool): kpoint is in Cartesian or fractional coordinates
            tol (float): tolerance below which coordinates are considered equal

        Returns:
            ([1x3 array] or None): if structure is not available returns None
        """
        if not self.structure:
            return None
        sg = SpacegroupAnalyzer(self.structure)
        symmops = sg.get_point_group_operations(cartesian=cartesian)
        points = np.dot(kpoint, [m.rotation_matrix for m in symmops])
        rm_list = []
        # identify and remove duplicates from the list of equivalent k-points:
        for i in range(len(points) - 1):
            for j in range(i + 1, len(points)):
                if np.allclose(pbc_diff(points[i], points[j]), [0, 0, 0], tol):
                    rm_list.append(i)
                    break
        return np.delete(points, rm_list, axis=0)

    def get_kpoint_degeneracy(self, kpoint, cartesian=False, tol: float = 1e-2):
        """
        Returns degeneracy of a given k-point based on structure symmetry
        Args:
            kpoint (1x3 array): coordinate of the k-point
            cartesian (bool): kpoint is in Cartesian or fractional coordinates
            tol (float): tolerance below which coordinates are considered equal

        Returns:
            (int or None): degeneracy or None if structure is not available
        """
        all_kpts = self.get_sym_eq_kpoints(kpoint, cartesian, tol=tol)
        if all_kpts is not None:
            return len(all_kpts)
        return None

    def as_dict(self):
        """
        JSON-serializable dict representation of BandStructure.
        """
        d = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "lattice_rec": self.lattice_rec.as_dict(),
            "efermi": self.efermi,
            "kpoints": [],
        }
        # kpoints are not kpoint objects dicts but are frac coords (this makes
        # the dict smaller and avoids the repetition of the lattice
        for k in self.kpoints:
            d["kpoints"].append(k.as_dict()["fcoords"])

        d["bands"] = {str(int(spin)): self.bands[spin].tolist() for spin in self.bands}
        d["is_metal"] = self.is_metal()
        vbm = self.get_vbm()
        d["vbm"] = {
            "energy": vbm["energy"],
            "kpoint_index": vbm["kpoint_index"],
            "band_index": {str(int(spin)): vbm["band_index"][spin] for spin in vbm["band_index"]},
            "projections": {str(spin): v.tolist() for spin, v in vbm["projections"].items()},
        }
        cbm = self.get_cbm()
        d["cbm"] = {
            "energy": cbm["energy"],
            "kpoint_index": cbm["kpoint_index"],
            "band_index": {str(int(spin)): cbm["band_index"][spin] for spin in cbm["band_index"]},
            "projections": {str(spin): v.tolist() for spin, v in cbm["projections"].items()},
        }
        d["band_gap"] = self.get_band_gap()
        d["labels_dict"] = {}
        d["is_spin_polarized"] = self.is_spin_polarized

        # MongoDB does not accept keys starting with $. Add a blank space to fix the problem
        for c, label in self.labels_dict.items():
            mongo_key = c if not c.startswith("$") else " " + c
            d["labels_dict"][mongo_key] = label.as_dict()["fcoords"]
        d["projections"] = {}
        if len(self.projections) != 0:
            d["structure"] = self.structure.as_dict()
            d["projections"] = {str(int(spin)): np.array(v).tolist() for spin, v in self.projections.items()}
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Create from dict.

        Args:
            A dict with all data for a band structure object.

        Returns:
            A BandStructure object
        """
        # Strip the label to recover initial string
        # (see trick used in as_dict to handle $ chars)
        labels_dict = {k.strip(): v for k, v in d["labels_dict"].items()}
        projections = {}
        structure = None
        if isinstance(list(d["bands"].values())[0], dict):
            eigenvals = {Spin(int(k)): np.array(d["bands"][k]["data"]) for k in d["bands"]}
        else:
            eigenvals = {Spin(int(k)): d["bands"][k] for k in d["bands"]}

        if "structure" in d:
            structure = Structure.from_dict(d["structure"])

        try:
            if d.get("projections"):
                if isinstance(d["projections"]["1"][0][0], dict):
                    raise ValueError("Old band structure dict format detected!")
                projections = {Spin(int(spin)): np.array(v) for spin, v in d["projections"].items()}

            return cls(
                d["kpoints"],
                eigenvals,
                Lattice(d["lattice_rec"]["matrix"]),
                d["efermi"],
                labels_dict,
                structure=structure,
                projections=projections,
            )

        except Exception:
            warnings.warn(
                "Trying from_dict failed. Now we are trying the old "
                "format. Please convert your BS dicts to the new "
                "format. The old format will be retired in pymatgen "
                "5.0."
            )
            return cls.from_old_dict(d)

    @classmethod
    def from_old_dict(cls, d):
        """
        Args:
            d (dict): A dict with all data for a band structure symm line
                object.
        Returns:
            A BandStructureSymmLine object
        """
        # Strip the label to recover initial string (see trick used in as_dict to handle $ chars)
        labels_dict = {k.strip(): v for k, v in d["labels_dict"].items()}
        projections = {}
        structure = None
        if "projections" in d and len(d["projections"]) != 0:
            structure = Structure.from_dict(d["structure"])
            projections = {}
            for spin in d["projections"]:
                dd = []
                for i in range(len(d["projections"][spin])):
                    ddd = []
                    for j in range(len(d["projections"][spin][i])):
                        dddd = []
                        for k in range(len(d["projections"][spin][i][j])):
                            ddddd = []
                            orb = Orbital(k).name
                            for l in range(len(d["projections"][spin][i][j][orb])):
                                ddddd.append(d["projections"][spin][i][j][orb][l])
                            dddd.append(np.array(ddddd))
                        ddd.append(np.array(dddd))
                    dd.append(np.array(ddd))
                projections[Spin(int(spin))] = np.array(dd)

        return BandStructure(
            d["kpoints"],
            {Spin(int(k)): d["bands"][k] for k in d["bands"]},
            Lattice(d["lattice_rec"]["matrix"]),
            d["efermi"],
            labels_dict,
            structure=structure,
            projections=projections,
        )


class BandStructureSymmLine(BandStructure, MSONable):
    r"""
    This object stores band structures along selected (symmetry) lines in the
    Brillouin zone. We call the different symmetry lines (ex: \\Gamma to Z)
    "branches".
    """

    def __init__(
        self,
        kpoints,
        eigenvals,
        lattice,
        efermi,
        labels_dict,
        coords_are_cartesian=False,
        structure=None,
        projections=None,
    ):
        """
        Args:
            kpoints: list of kpoint as numpy arrays, in frac_coords of the
                given lattice by default
            eigenvals: dict of energies for spin up and spin down
                {Spin.up:[][],Spin.down:[][]}, the first index of the array
                [][] refers to the band and the second to the index of the
                kpoint. The kpoints are ordered according to the order of the
                kpoints array. If the band structure is not spin polarized, we
                only store one data set under Spin.up.
            lattice: The reciprocal lattice.
                Pymatgen uses the physics convention of reciprocal lattice vectors
                WITH a 2*pi coefficient
            efermi: fermi energy
            label_dict: (dict) of {} this link a kpoint (in frac coords or
                Cartesian coordinates depending on the coords).
            coords_are_cartesian: Whether coordinates are cartesian.
            structure: The crystal structure (as a pymatgen Structure object)
                associated with the band structure. This is needed if we
                provide projections to the band structure.
            projections: dict of orbital projections as {spin: ndarray}. The
                indices of the ndarrayare [band_index, kpoint_index, orbital_index,
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
        one_group = []
        branches_tmp = []
        # get labels and distance for each kpoint
        previous_kpoint = self.kpoints[0]
        previous_distance = 0.0

        previous_label = self.kpoints[0].label
        for i, kpt in enumerate(self.kpoints):
            label = kpt.label
            if label is not None and previous_label is not None:
                self.distance.append(previous_distance)
            else:
                self.distance.append(np.linalg.norm(kpt.cart_coords - previous_kpoint.cart_coords) + previous_distance)
            previous_kpoint = kpt
            previous_distance = self.distance[i]
            if label:
                if previous_label:
                    if len(one_group) != 0:
                        branches_tmp.append(one_group)
                    one_group = []
            previous_label = label
            one_group.append(i)

        if len(one_group) != 0:
            branches_tmp.append(one_group)
        for b in branches_tmp:
            self.branches.append(
                {
                    "start_index": b[0],
                    "end_index": b[-1],
                    "name": str(self.kpoints[b[0]].label) + "-" + str(self.kpoints[b[-1]].label),
                }
            )

        self.is_spin_polarized = False
        if len(self.bands) == 2:
            self.is_spin_polarized = True

    def get_equivalent_kpoints(self, index):
        """
        Returns the list of kpoint indices equivalent (meaning they are the
        same frac coords) to the given one.

        Args:
            index: the kpoint index

        Returns:
            a list of equivalent indices

        TODO: now it uses the label we might want to use coordinates instead
        (in case there was a mislabel)
        """
        # if the kpoint has no label it can"t have a repetition along the band
        # structure line object

        if self.kpoints[index].label is None:
            return [index]

        list_index_kpoints = []
        for i, kpt in enumerate(self.kpoints):
            if kpt.label == self.kpoints[index].label:
                list_index_kpoints.append(i)

        return list_index_kpoints

    def get_branch(self, index):
        r"""
        Returns in what branch(es) is the kpoint. There can be several
        branches.

        Args:
            index: the kpoint index

        Returns:
            A list of dictionaries [{"name","start_index","end_index","index"}]
            indicating all branches in which the k_point is. It takes into
            account the fact that one kpoint (e.g., \\Gamma) can be in several
            branches
        """
        to_return = []
        for i in self.get_equivalent_kpoints(index):
            for b in self.branches:
                if b["start_index"] <= i <= b["end_index"]:
                    to_return.append(
                        {
                            "name": b["name"],
                            "start_index": b["start_index"],
                            "end_index": b["end_index"],
                            "index": i,
                        }
                    )
        return to_return

    def apply_scissor(self, new_band_gap):
        """
        Apply a scissor operator (shift of the CBM) to fit the given band gap.
        If it's a metal. We look for the band crossing the fermi level
        and shift this one up. This will not work all the time for metals!

        Args:
            new_band_gap: the band gap the scissor band structure need to have.

        Returns:
            a BandStructureSymmLine object with the applied scissor shift
        """
        if self.is_metal():
            # moves then the highest index band crossing the fermi level
            # find this band...
            max_index = -1000
            # spin_index = None
            for i in range(self.nb_bands):
                below = False
                above = False
                for j in range(len(self.kpoints)):
                    if self.bands[Spin.up][i][j] < self.efermi:
                        below = True
                    if self.bands[Spin.up][i][j] > self.efermi:
                        above = True
                if above and below:
                    if i > max_index:
                        max_index = i
                        # spin_index = Spin.up
                if self.is_spin_polarized:
                    below = False
                    above = False
                    for j in range(len(self.kpoints)):
                        if self.bands[Spin.down][i][j] < self.efermi:
                            below = True
                        if self.bands[Spin.down][i][j] > self.efermi:
                            above = True
                    if above and below:
                        if i > max_index:
                            max_index = i
                            # spin_index = Spin.down
            old_dict = self.as_dict()
            shift = new_band_gap
            for spin in old_dict["bands"]:
                for k in range(len(old_dict["bands"][spin])):
                    for v in range(len(old_dict["bands"][spin][k])):
                        if k >= max_index:
                            old_dict["bands"][spin][k][v] = old_dict["bands"][spin][k][v] + shift
        else:

            shift = new_band_gap - self.get_band_gap()["energy"]
            old_dict = self.as_dict()
            for spin in old_dict["bands"]:
                for k in range(len(old_dict["bands"][spin])):
                    for v in range(len(old_dict["bands"][spin][k])):
                        if old_dict["bands"][spin][k][v] >= old_dict["cbm"]["energy"]:
                            old_dict["bands"][spin][k][v] = old_dict["bands"][spin][k][v] + shift
            old_dict["efermi"] = old_dict["efermi"] + shift
        return self.from_dict(old_dict)

    def as_dict(self):
        """
        JSON-serializable dict representation of BandStructureSymmLine.
        """
        d = super().as_dict()
        d["branches"] = self.branches
        return d


class LobsterBandStructureSymmLine(BandStructureSymmLine):
    """
    Lobster subclass of BandStructure with customized functions.
    """

    def as_dict(self):
        """
        JSON-serializable dict representation of BandStructureSymmLine.
        """
        d = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "lattice_rec": self.lattice_rec.as_dict(),
            "efermi": self.efermi,
            "kpoints": [],
        }
        # kpoints are not kpoint objects dicts but are frac coords (this makes
        # the dict smaller and avoids the repetition of the lattice
        for k in self.kpoints:
            d["kpoints"].append(k.as_dict()["fcoords"])
        d["branches"] = self.branches
        d["bands"] = {str(int(spin)): self.bands[spin].tolist() for spin in self.bands}
        d["is_metal"] = self.is_metal()
        vbm = self.get_vbm()
        d["vbm"] = {
            "energy": vbm["energy"],
            "kpoint_index": [int(x) for x in vbm["kpoint_index"]],
            "band_index": {str(int(spin)): vbm["band_index"][spin] for spin in vbm["band_index"]},
            "projections": {str(spin): v for spin, v in vbm["projections"].items()},
        }
        cbm = self.get_cbm()
        d["cbm"] = {
            "energy": cbm["energy"],
            "kpoint_index": [int(x) for x in cbm["kpoint_index"]],
            "band_index": {str(int(spin)): cbm["band_index"][spin] for spin in cbm["band_index"]},
            "projections": {str(spin): v for spin, v in cbm["projections"].items()},
        }
        d["band_gap"] = self.get_band_gap()
        d["labels_dict"] = {}
        d["is_spin_polarized"] = self.is_spin_polarized
        # MongoDB does not accept keys starting with $. Add a blank space to fix the problem
        for c, label in self.labels_dict.items():
            mongo_key = c if not c.startswith("$") else " " + c
            d["labels_dict"][mongo_key] = label.as_dict()["fcoords"]
        if len(self.projections) != 0:
            d["structure"] = self.structure.as_dict()
            d["projections"] = {str(int(spin)): np.array(v).tolist() for spin, v in self.projections.items()}
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Args:
            d (dict): A dict with all data for a band structure symm line
                object.

        Returns:
            A BandStructureSymmLine object
        """
        try:
            # Strip the label to recover initial string (see trick used in as_dict to handle $ chars)
            labels_dict = {k.strip(): v for k, v in d["labels_dict"].items()}
            projections = {}
            structure = None
            if d.get("projections"):
                if isinstance(d["projections"]["1"][0][0], dict):
                    raise ValueError("Old band structure dict format detected!")
                structure = Structure.from_dict(d["structure"])
                projections = {Spin(int(spin)): np.array(v) for spin, v in d["projections"].items()}

            return LobsterBandStructureSymmLine(
                d["kpoints"],
                {Spin(int(k)): d["bands"][k] for k in d["bands"]},
                Lattice(d["lattice_rec"]["matrix"]),
                d["efermi"],
                labels_dict,
                structure=structure,
                projections=projections,
            )
        except Exception:
            warnings.warn(
                "Trying from_dict failed. Now we are trying the old "
                "format. Please convert your BS dicts to the new "
                "format. The old format will be retired in pymatgen "
                "5.0."
            )
            return LobsterBandStructureSymmLine.from_old_dict(d)

    @classmethod
    def from_old_dict(cls, d):
        """
        Args:
            d (dict): A dict with all data for a band structure symm line
                object.
        Returns:
            A BandStructureSymmLine object
        """
        # Strip the label to recover initial string (see trick used in as_dict to handle $ chars)
        labels_dict = {k.strip(): v for k, v in d["labels_dict"].items()}
        projections = {}
        structure = None
        if "projections" in d and len(d["projections"]) != 0:
            structure = Structure.from_dict(d["structure"])
            projections = {}
            for spin in d["projections"]:
                dd = []
                for i in range(len(d["projections"][spin])):
                    ddd = []
                    for j in range(len(d["projections"][spin][i])):
                        ddd.append(d["projections"][spin][i][j])
                    dd.append(np.array(ddd))
                projections[Spin(int(spin))] = np.array(dd)

        return LobsterBandStructureSymmLine(
            d["kpoints"],
            {Spin(int(k)): d["bands"][k] for k in d["bands"]},
            Lattice(d["lattice_rec"]["matrix"]),
            d["efermi"],
            labels_dict,
            structure=structure,
            projections=projections,
        )

    def get_projection_on_elements(self):
        """
        Method returning a dictionary of projections on elements.
        It sums over all available orbitals for each element.

        Returns:
            a dictionary in the {Spin.up:[][{Element:values}],
            Spin.down:[][{Element:values}]} format
            if there is no projections in the band structure
            returns an empty dict
        """
        result = {}
        for spin, v in self.projections.items():
            result[spin] = [
                [collections.defaultdict(float) for i in range(len(self.kpoints))] for j in range(self.nb_bands)
            ]
            for i, j in itertools.product(range(self.nb_bands), range(len(self.kpoints))):
                for key, item in v[i][j].items():
                    for item2 in item.values():
                        specie = str(Element(re.split(r"[0-9]+", key)[0]))
                        result[spin][i][j][specie] += item2
        return result

    def get_projections_on_elements_and_orbitals(self, el_orb_spec):
        """
        Method returning a dictionary of projections on elements and specific
        orbitals

        Args:
            el_orb_spec: A dictionary of Elements and Orbitals for which we want
                to have projections on. It is given as: {Element:[orbitals]},
                e.g., {'Si':['3s','3p']} or {'Si':['3s','3p_x', '3p_y', '3p_z']} depending on input files

        Returns:
            A dictionary of projections on elements in the
            {Spin.up:[][{Element:{orb:values}}],
            Spin.down:[][{Element:{orb:values}}]} format
            if there is no projections in the band structure returns an empty
            dict.
        """
        result = {}
        el_orb_spec = {get_el_sp(el): orbs for el, orbs in el_orb_spec.items()}
        for spin, v in self.projections.items():
            result[spin] = [
                [{str(e): collections.defaultdict(float) for e in el_orb_spec} for i in range(len(self.kpoints))]
                for j in range(self.nb_bands)
            ]

            for i, j in itertools.product(range(self.nb_bands), range(len(self.kpoints))):
                for key, item in v[i][j].items():
                    for key2, item2 in item.items():
                        specie = str(Element(re.split(r"[0-9]+", key)[0]))
                        if get_el_sp(str(specie)) in el_orb_spec:
                            if key2 in el_orb_spec[get_el_sp(str(specie))]:
                                result[spin][i][j][specie][key2] += item2
        return result


def get_reconstructed_band_structure(list_bs, efermi=None):
    """
    This method takes a list of band structures and reconstructs
    one band structure object from all of them.

    This is typically very useful when you split non self consistent
    band structure runs in several independent jobs and want to merge back
    the results

    Args:
        list_bs: A list of BandStructure or BandStructureSymmLine objects.
        efermi: The Fermi energy of the reconstructed band structure. If
            None is assigned an average of all the Fermi energy in each
            object in the list_bs is used.

    Returns:
        A BandStructure or BandStructureSymmLine object (depending on
        the type of the list_bs objects)
    """
    if efermi is None:
        efermi = sum(b.efermi for b in list_bs) / len(list_bs)

    kpoints = []
    labels_dict = {}
    rec_lattice = list_bs[0].lattice_rec
    nb_bands = min(list_bs[i].nb_bands for i in range(len(list_bs)))

    kpoints = np.concatenate([[k.frac_coords for k in bs.kpoints] for bs in list_bs])
    dicts = [bs.labels_dict for bs in list_bs]
    labels_dict = {k: v.frac_coords for d in dicts for k, v in d.items()}

    eigenvals = {}
    eigenvals[Spin.up] = np.concatenate([bs.bands[Spin.up][:nb_bands] for bs in list_bs], axis=1)

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
            labels_dict,
            structure=list_bs[0].structure,
            projections=projections,
        )
    return BandStructure(
        kpoints,
        eigenvals,
        rec_lattice,
        efermi,
        labels_dict,
        structure=list_bs[0].structure,
        projections=projections,
    )
