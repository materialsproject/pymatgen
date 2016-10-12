# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import numpy as np
import math
import itertools
import collections
import warnings

from pymatgen.core.periodic_table import get_el_sp
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.electronic_structure.core import Spin, Orbital
from monty.json import MSONable

"""
This module provides classes to define everything related to band structures.
"""

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
    or cartesian coordinates syntax similar than the site object in
    pymatgen.core.structure.

    Args:
        coords: coordinate of the kpoint as a numpy array
        lattice: A pymatgen.core.lattice.Lattice lattice object representing
            the reciprocal lattice of the kpoint
        to_unit_cell: Translates fractional coordinate to the basic unit
            cell, i.e., all fractional coordinates satisfy 0 <= a < 1.
            Defaults to False.
        coords_are_cartesian: Boolean indicating if the coordinates given are
            in cartesian or fractional coordinates (by default fractional)
        label: the label of the kpoint if any (None by default)
    """

    def __init__(self, coords, lattice, to_unit_cell=False,
                 coords_are_cartesian=False, label=None):
        self._lattice = lattice
        self._fcoords = lattice.get_fractional_coords(coords) \
            if coords_are_cartesian else coords
        self._label = label

        if to_unit_cell:
            for i in range(len(self._fcoords)):
                self._fcoords[i] -= math.floor(self._fcoords[i])

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
        The cartesian coordinates of the kpoint as a numpy array
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
        Returns a string with fractional, cartesian coordinates and label
        """
        return "{} {} {}".format(self.frac_coords, self.cart_coords,
                                 self.label)

    def as_dict(self):
        """
        Json-serializable dict representation of a kpoint
        """
        return {"lattice": self.lattice.as_dict(),
                "fcoords": list(self.frac_coords),
                "ccoords": list(self.cart_coords), "label": self.label,
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}


class BandStructure(object):
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
        efermi: fermi energy
        labels_dict: (dict) of {} this links a kpoint (in frac coords or
            cartesian coordinates depending on the coords) to a label.
        coords_are_cartesian: Whether coordinates are cartesian.
        structure: The crystal structure (as a pymatgen Structure object)
            associated with the band structure. This is needed if we
            provide projections to the band structure
        projections: dict of orbital projections for spin up and spin down
            {Spin.up:[][{Orbital:[]}],Spin.down:[][{Orbital:[]}]. The
            format follows the one from eigenvals: The first index of the
            array refers to the band and the second to the index of the
            kpoint. The kpoints are ordered according to the order of the
            kpoints array. For each band and kpoint, we associate a
            dictionary indicating projections on orbitals and on different
            sites the keys of the dictionary are Orbital objects and the
            values are the projections on each site ordered as in the
            structure object. If the band structure is not spin polarized,
            we only store one data set under Spin.up.
    """

    def __init__(self, kpoints, eigenvals, lattice, efermi, labels_dict=None,
                 coords_are_cartesian=False, structure=None, projections=None):
        self.efermi = efermi
        self.lattice_rec = lattice
        self.kpoints = []
        self.labels_dict = {}
        self.structure = structure
        self.projections = projections or {}

        if labels_dict is None:
            labels_dict = {}

        if len(self.projections) != 0 and self.structure is None:
            raise Exception("if projections are provided a structure object"
                            " needs also to be given")

        for k in kpoints:
            # let see if this kpoint has been assigned a label
            label = None
            for c in labels_dict:
                if np.linalg.norm(k - np.array(labels_dict[c])) < 0.0001:
                    label = c
                    self.labels_dict[label] = Kpoint(
                        k, lattice, label=label,
                        coords_are_cartesian=coords_are_cartesian)
            self.kpoints.append(
                Kpoint(k, lattice, label=label,
                       coords_are_cartesian=coords_are_cartesian))
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
            result[spin] = [[collections.defaultdict(float)
                             for i in range(len(self.kpoints))]
                            for j in range(self.nb_bands)]
            for i, j, k in itertools.product(range(self.nb_bands),
                                             range(len(self.kpoints)),
                                             range(structure.num_sites)):
                result[spin][i][j][str(structure[k].specie)] += np.sum(
                    v[i, j, :, k])
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
            result[spin] = [[{str(e): collections.defaultdict(float)
                              for e in el_orb_spec}
                             for i in range(len(self.kpoints))]
                            for j in range(self.nb_bands)]

            for i, j, k in itertools.product(
                    range(self.nb_bands), range(len(self.kpoints)),
                    range(structure.num_sites)):
                sp = structure[k].specie
                for orb_i in range(len(v[i][j])):
                    o = Orbital(orb_i).name[0]
                    if sp in el_orb_spec:
                        if o in el_orb_spec[sp]:
                            result[spin][i][j][str(sp)][o] += v[i][j][
                                orb_i][k]
        return result

    def is_metal(self):
        """
        Check if the band structure indicates a metal by looking if the fermi
        level crosses a band.

        Returns:
            True if a metal, False if not
        """
        for spin, values in self.bands.items():
            for i in range(self.nb_bands):
                if np.any(values[i, :] < self.efermi) and \
                        np.any(values[i, :] > self.efermi):
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
            kpoint vbm. Please note that there can be several
            kpoint_indices relating to the same kpoint (e.g., Gamma can
            occur at different spots in the band structure line plot)
            - "kpoint": The kpoint (as a kpoint object)
            - "energy": The energy of the VBM
            - "projections": The projections along sites and orbitals of the
            VBM if any projection data is available (else it is an empty
            dictionnary). The format is similar to the projections field in
            BandStructure: {spin:{'Orbital': [proj]}} where the array
            [proj] is ordered according to the sites in structure
    """
        if self.is_metal():
            return {"band_index": [], "kpoint_index": [],
                    "kpoint": [], "energy": None, "projections": {}}
        max_tmp = -float("inf")
        index = None
        kpointvbm = None
        for spin, v in self.bands.items():
            for i, j in zip(*np.where(v < self.efermi)):
                if v[i, j] > max_tmp:
                    max_tmp = float(v[i, j])
                    index = j
                    kpointvbm = self.kpoints[j]

        list_ind_kpts = []
        if kpointvbm.label is not None:
            for i in range(len(self.kpoints)):
                if self.kpoints[i].label == kpointvbm.label:
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
        for spin, v in self.projections.items():
            if len(list_ind_band[spin]) == 0:
                continue
            proj[spin] = v[list_ind_band[spin][0]][list_ind_kpts[0]]
        return {'band_index': list_ind_band,
                'kpoint_index': list_ind_kpts,
                'kpoint': kpointvbm, 'energy': max_tmp,
                'projections': proj}

    def get_cbm(self):
        """
        Returns data about the CBM.

        Returns:
            {"band_index","kpoint_index","kpoint","energy"}
            - "band_index": A dict with spin keys pointing to a list of the
            indices of the band containing the VBM (please note that you
            can have several bands sharing the VBM) {Spin.up:[],
            Spin.down:[]}
            - "kpoint_index": The list of indices in self.kpoints for the
            kpoint vbm. Please note that there can be several
            kpoint_indices relating to the same kpoint (e.g., Gamma can
            occur at different spots in the band structure line plot)
            - "kpoint": The kpoint (as a kpoint object)
            - "energy": The energy of the VBM
            - "projections": The projections along sites and orbitals of the
            VBM if any projection data is available (else it is an empty
            dictionnary). The format is similar to the projections field in
            BandStructure: {spin:{'Orbital': [proj]}} where the array
            [proj] is ordered according to the sites in structure
        """
        if self.is_metal():
            return {"band_index": [], "kpoint_index": [],
                    "kpoint": [], "energy": None, "projections": {}}
        max_tmp = float("inf")

        index = None
        kpointcbm = None
        for spin, v in self.bands.items():
            for i, j in zip(*np.where(v > self.efermi)):
                if v[i, j] < max_tmp:
                    max_tmp = float(v[i, j])
                    index = j
                    kpointcbm = self.kpoints[j]

        list_index_kpoints = []
        if kpointcbm.label is not None:
            for i in range(len(self.kpoints)):
                if self.kpoints[i].label == kpointcbm.label:
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
        for spin, v in self.projections.items():
            if len(list_index_band[spin]) == 0:
                continue
            proj[spin] = v[list_index_band[spin][0]][list_index_kpoints[0]]

        return {'band_index': list_index_band,
                'kpoint_index': list_index_kpoints,
                'kpoint': kpointcbm, 'energy': max_tmp,
                'projections': proj}

    def get_band_gap(self):
        """
        Returns band gap data.

        Returns:
            A dict {"energy","direct","transition"}:
            "energy": band gap energy
            "direct": A boolean telling if the gap is direct or not
            "transition": kpoint labels of the transition (e.g., "\Gamma-X")
        """
        if self.is_metal():
            return {"energy": 0.0, "direct": False, "transition": None}
        cbm = self.get_cbm()
        vbm = self.get_vbm()
        result = dict(direct=False, energy=0.0, transition=None)

        result["energy"] = cbm["energy"] - vbm["energy"]

        if (cbm["kpoint"].label is not None and cbm["kpoint"].label == vbm[
            "kpoint"].label) \
                or np.linalg.norm(cbm["kpoint"].cart_coords
                                          - vbm["kpoint"].cart_coords) < 0.01:
            result["direct"] = True

        result["transition"] = "-".join(
            [str(c.label) if c.label is not None else
             str("(") + ",".join(["{0:.3f}".format(c.frac_coords[i])
                                  for i in range(3)])
             + str(")") for c in [vbm["kpoint"], cbm["kpoint"]]])

        return result

    def get_direct_band_gap(self):
        """
        Returns the direct band gap.

        Returns:
             the value of the direct band gap
        """
        if self.is_metal():
            return 0.0
        cb = collections.defaultdict(list)
        vb = collections.defaultdict(list)
        for spin, v in self.bands.items():
            for i in range(len(v)):
                for j in range(len(self.kpoints)):
                    if self.bands[Spin.up][i][j] > self.efermi:
                        cb[spin].append(v[i][j])
                        vb[spin].append(v[i - 1][j])
        diff = []
        for i in range(len(self.kpoints)):
            diff.append(min(
                [cb[Spin.up][i], cb.get(Spin.down, cb[Spin.up])[i]])
                        - max(
                [vb[Spin.up][i], vb.get(Spin.down, vb[Spin.up])[i]]))
        return min(diff)

    def as_dict(self):
        """
        Json-serializable dict representation of BandStructureSymmLine.
        """
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "lattice_rec": self.lattice_rec.as_dict(), "efermi": self.efermi,
             "kpoints": []}
        # kpoints are not kpoint objects dicts but are frac coords (this makes
        # the dict smaller and avoids the repetition of the lattice
        for k in self.kpoints:
            d["kpoints"].append(k.as_dict()["fcoords"])
        d["bands"] = {str(int(spin)): self.bands[spin]
                      for spin in self.bands}
        d["is_metal"] = self.is_metal()
        vbm = self.get_vbm()
        d["vbm"] = {"energy": vbm["energy"],
                    "kpoint_index": vbm["kpoint_index"],
                    "band_index": {str(int(spin)): vbm["band_index"][spin]
                                   for spin in vbm["band_index"]},
                    'projections': {str(spin): {str(orb):
                                                    vbm['projections'][spin][
                                                        orb]
                                                for orb in
                                                vbm['projections'][spin]}
                                    for spin in vbm['projections']}}
        cbm = self.get_cbm()
        d['cbm'] = {'energy': cbm['energy'],
                    'kpoint_index': cbm['kpoint_index'],
                    'band_index': {str(int(spin)): cbm['band_index'][spin]
                                   for spin in cbm['band_index']},
                    'projections': {str(spin): {str(orb):
                                                    cbm['projections'][spin][
                                                        orb]
                                                for orb in
                                                cbm['projections'][spin]}
                                    for spin in cbm['projections']}}
        d['band_gap'] = self.get_band_gap()
        d['labels_dict'] = {}
        d['is_spin_polarized'] = self.is_spin_polarized
        for c in self.labels_dict:
            d['labels_dict'][c] = self.labels_dict[c].as_dict()['fcoords']
        d['projections'] = {}
        if len(self.projections) != 0:
            d['structure'] = self.structure.as_dict()
            d['projections'] = {str(int(spin)): v.tolist()
                                for spin, v in self.projections.items()}
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
        labels_dict = d['labels_dict']
        projections = {}
        structure = None
        if 'structure' in d:
            structure = Structure.from_dict(d['structure'])
        if d.get('projections'):
            projections = {Spin(int(spin)): np.array(v)
                           for spin, v in d["projections"].items()}

        return BandStructure(
            d['kpoints'], {Spin(int(k)): d['bands'][k]
                           for k in d['bands']},
            Lattice(d['lattice_rec']['matrix']), d['efermi'],
            labels_dict, structure=structure, projections=projections)


class BandStructureSymmLine(BandStructure, MSONable):
    """
    This object stores band structures along selected (symmetry) lines in the
    Brillouin zone. We call the different symmetry lines (ex: \Gamma to Z)
    "branches".

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
            cartesian coordinates depending on the coords).
        coords_are_cartesian: Whether coordinates are cartesian.
        structure: The crystal structure (as a pymatgen Structure object)
            associated with the band structure. This is needed if we
            provide projections to the band structure.
        projections: dict of orbital projections for spin up and spin down
            {Spin.up:[][{Orbital:[]}],Spin.down:[][{Orbital:[]}]. The
            format follows the one from eigenvals: the first index of the
            array refers to the band and the second to the index of the
            kpoint. The kpoints are ordered according to the order of the
            kpoints array. For each band and kpoint, we associate a
            dictionary indicating projections on orbitals and on different
            sites the keys of the dictionary are Orbital objects and the
            values are the projections on each site ordered as in the
            structure object. If the band structure is not spin polarized,
            we only store one data set under Spin.up.
    """

    def __init__(self, kpoints, eigenvals, lattice, efermi, labels_dict,
                 coords_are_cartesian=False, structure=None,
                 projections=None):
        super(BandStructureSymmLine, self).__init__(
            kpoints, eigenvals, lattice, efermi, labels_dict,
            coords_are_cartesian, structure, projections)
        self.distance = []
        self.branches = []
        one_group = []
        branches_tmp = []
        # get labels and distance for each kpoint
        previous_kpoint = self.kpoints[0]
        previous_distance = 0.0

        previous_label = self.kpoints[0].label
        for i in range(len(self.kpoints)):
            label = self.kpoints[i].label
            if label is not None and previous_label is not None:
                self.distance.append(previous_distance)
            else:
                self.distance.append(
                    np.linalg.norm(self.kpoints[i].cart_coords -
                                   previous_kpoint.cart_coords) +
                    previous_distance)
            previous_kpoint = self.kpoints[i]
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
                {"start_index": b[0], "end_index": b[-1],
                 "name": str(self.kpoints[b[0]].label) + "-" +
                         str(self.kpoints[b[-1]].label)})

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
        for i in range(len(self.kpoints)):
            if self.kpoints[i].label == self.kpoints[index].label:
                list_index_kpoints.append(i)

        return list_index_kpoints

    def get_branch(self, index):
        """
        Returns in what branch(es) is the kpoint. There can be several
        branches.

        Args:
            index: the kpoint index

        Returns:
            A list of dictionaries [{"name","start_index","end_index","index"}]
            indicating all branches in which the k_point is. It takes into
            account the fact that one kpoint (e.g., \Gamma) can be in several
            branches
        """
        to_return = []
        for i in self.get_equivalent_kpoints(index):
            for b in self.branches:
                if b["start_index"] <= i <= b["end_index"]:
                    to_return.append({"name": b["name"],
                                      "start_index": b["start_index"],
                                      "end_index": b["end_index"],
                                      "index": i})
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
            for spin in old_dict['bands']:
                for k in range(len(old_dict['bands'][spin])):
                    for v in range(len(old_dict['bands'][spin][k])):
                        if k >= max_index:
                            old_dict['bands'][spin][k][v] = \
                                old_dict['bands'][spin][k][v] + shift
        else:

            shift = new_band_gap - self.get_band_gap()['energy']
            old_dict = self.as_dict()
            for spin in old_dict['bands']:
                for k in range(len(old_dict['bands'][spin])):
                    for v in range(len(old_dict['bands'][spin][k])):
                        if old_dict['bands'][spin][k][v] >= \
                                old_dict['cbm']['energy']:
                            old_dict['bands'][spin][k][v] = \
                                old_dict['bands'][spin][k][v] + shift
            old_dict['efermi'] = old_dict['efermi'] + shift
            return BandStructureSymmLine.from_dict(old_dict)

    def as_dict(self):
        """
        Json-serializable dict representation of BandStructureSymmLine.
        """

        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "lattice_rec": self.lattice_rec.as_dict(), "efermi": self.efermi,
             "kpoints": []}
        # kpoints are not kpoint objects dicts but are frac coords (this makes
        # the dict smaller and avoids the repetition of the lattice
        for k in self.kpoints:
            d["kpoints"].append(k.as_dict()["fcoords"])
        d["branches"] = self.branches
        d["bands"] = {str(int(spin)): self.bands[spin].tolist()
                      for spin in self.bands}
        d["is_metal"] = self.is_metal()
        vbm = self.get_vbm()
        d["vbm"] = {"energy": vbm["energy"],
                    "kpoint_index": vbm["kpoint_index"],
                    "band_index": {str(int(spin)): vbm["band_index"][spin]
                                   for spin in vbm["band_index"]},
                    'projections': {str(spin): v.tolist() for spin, v in vbm[
                        'projections'].items()}}
        cbm = self.get_cbm()
        d['cbm'] = {'energy': cbm['energy'],
                    'kpoint_index': cbm['kpoint_index'],
                    'band_index': {str(int(spin)): cbm['band_index'][spin]
                                   for spin in cbm['band_index']},
                    'projections': {str(spin): v.tolist() for spin, v in cbm[
                        'projections'].items()}}
        d['band_gap'] = self.get_band_gap()
        d['labels_dict'] = {}
        d['is_spin_polarized'] = self.is_spin_polarized
        # MongoDB does not accept keys starting with $. Add a blanck space to fix the problem
        for c in self.labels_dict:
            mongo_key = c if not c.startswith("$") else " " + c
            d['labels_dict'][mongo_key] = self.labels_dict[c].as_dict()[
                'fcoords']
        if len(self.projections) != 0:
            d['structure'] = self.structure.as_dict()
            d['projections'] = {str(int(spin)): v.tolist()
                                for spin, v in self.projections.items()}
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
            labels_dict = {k.strip(): v for k, v in d['labels_dict'].items()}
            projections = {}
            structure = None
            if d.get('projections'):
                structure = Structure.from_dict(d['structure'])
                projections = {Spin(int(spin)): np.array(v)
                               for spin, v in d["projections"].items()}

            return BandStructureSymmLine(
                d['kpoints'], {Spin(int(k)): d['bands'][k]
                               for k in d['bands']},
                Lattice(d['lattice_rec']['matrix']), d['efermi'],
                labels_dict, structure=structure, projections=projections)
        except:
            warnings.warn("Trying from_dict failed. Now we are trying the old "
                          "format. Please convert your BS dicts to the new "
                          "format. The old format will be retired in pymatgen "
                          "5.0.")
            return BandStructureSymmLine.from_old_dict(d)

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
        labels_dict = {k.strip(): v for k, v in d['labels_dict'].items()}
        projections = {}
        structure = None
        if 'projections' in d and len(d['projections']) != 0:
            structure = Structure.from_dict(d['structure'])
            projections = {}
            for spin in d['projections']:
                dd = []
                for i in range(len(d['projections'][spin])):
                    ddd = []
                    for j in range(len(d['projections'][spin][i])):
                        dddd = []
                        for k in range(len(d['projections'][spin][i][j])):
                            ddddd = []
                            orb = Orbital(k).name
                            for l in range(len(d['projections'][spin][i][j][
                                                   orb])):
                                ddddd.append(d['projections'][spin][i][j][
                                                orb][l])
                            dddd.append(ddddd)
                        ddd.append(dddd)
                    dd.append(ddd)
                projections[spin] = dd

        return BandStructureSymmLine(
            d['kpoints'], {Spin(int(k)): d['bands'][k]
                           for k in d['bands']},
            Lattice(d['lattice_rec']['matrix']), d['efermi'],
            labels_dict, structure=structure, projections=projections)


def get_reconstructed_band_structure(list_bs, efermi=None):
    """
    This method takes a list of band structures
    and reconstruct one band structure object from all of them

    this is typically very useful when you split non self consistent
    band structure runs in several independent jobs and want to merge back
    the results

    Args:
        list_bs: A list of BandStructure
        efermi: The fermi energy of the reconstructed band structure. If
            None is assigned an average of all the fermi energy in each
            object in the list_bs is used.

    Returns:
        A BandStructure or BandStructureSymmLine object (depending on
        the type of the list_bs objects)
    """
    if efermi is None:
        efermi = sum([b.efermi for b in list_bs]) / len(list_bs)

    kpoints = []
    labels_dict = {}
    rec_lattice = list_bs[0]._lattice_rec
    nb_bands = min([list_bs[i]._nb_bands for i in range(len(list_bs))])

    for bs in list_bs:
        for k in bs._kpoints:
            kpoints.append(k.frac_coords)
        for k, v in bs._labels_dict.items():
            labels_dict[k] = v.frac_coords
    eigenvals = {Spin.up: [list_bs[0]._bands[Spin.up][i]
                           for i in range(nb_bands)]}
    for i in range(nb_bands):
        for bs in list_bs[1:]:
            for e in bs._bands[Spin.up][i]:
                eigenvals[Spin.up][i].append(e)
    if list_bs[0].is_spin_polarized:
        eigenvals[Spin.down] = [list_bs[0]._bands[Spin.down][i]
                                for i in range(nb_bands)]
        for i in range(nb_bands):
            for bs in list_bs[1:]:
                for e in bs._bands[Spin.down][i]:
                    eigenvals[Spin.down][i].append(e)
    projections = {}
    if len(list_bs[0]._projections) != 0:
        projections = {Spin.up: [list_bs[0]._projections[Spin.up][i]
                                 for i in range(nb_bands)]}
        for i in range(nb_bands):
            for bs in list_bs[1:]:
                projections[Spin.up][i].extend(bs._projections[Spin.up][i])
        if list_bs[0].is_spin_polarized:
            projections[Spin.down] = [list_bs[0]._projections[Spin.down][i]
                                      for i in range(nb_bands)]
            for i in range(nb_bands):
                for bs in list_bs[1:]:
                    projections[Spin.down][i].extend(
                        bs._projections[Spin.down][i])

    if isinstance(list_bs[0], BandStructureSymmLine):
        return BandStructureSymmLine(kpoints, eigenvals, rec_lattice,
                                     efermi, labels_dict,
                                     structure=list_bs[0]._structure,
                                     projections=projections)
    else:
        return BandStructure(kpoints, eigenvals, rec_lattice, efermi,
                             labels_dict, structure=list_bs[0]._structure,
                             projections=projections)
