# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.electronic_structure.bandstructure import Kpoint
from monty.json import MSONable, MontyDecoder, MontyEncoder

"""
This module provides classes to define a phonon band structure.
"""


class PhononBandStructure(MSONable):
    """
    This is the most generic phonon band structure data possible
    it's defined by a list of qpoints + frequencies for each of them.
    Additional information may be given for frequencies at Gamma, where
    non-analytical contribution may be taken into account.
    """

    def __init__(self, qpoints, frequencies, lattice, nac_frequencies=None, eigendisplacements=None,
                 nac_eigendisplacements=None, labels_dict=None, coords_are_cartesian=False,
                 structure=None):
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
                frequencies in cartesian coordinates. A numpy array of complex
                numbers with shape (3*len(structure), len(qpoints), len(structure), 3).
                he First index of the array refers to the band, the second to the index
                of the qpoint, the third to the atom in the structure and the fourth
                to the cartesian coordinates.
            nac_eigendisplacements: the phonon eigendisplacements associated to the
                non-analytical frequencies in nac_frequencies in cartesian coordinates.
                A list of tuples. The first element of each tuple should be a list
                defining the direction. The second element containing a numpy array of
                complex numbers with shape (3*len(structure), len(structure), 3).
            labels_dict: (dict) of {} this links a qpoint (in frac coords or
                cartesian coordinates depending on the coords) to a label.
            coords_are_cartesian: Whether the qpoint coordinates are cartesian.
            structure: The crystal structure (as a pymatgen Structure object)
                associated with the band structure. This is needed if we
                provide projections to the band structure
        """
        self.lattice_rec = lattice
        self.qpoints = []
        self.labels_dict = {}
        self.structure = structure
        if eigendisplacements is None:
            eigendisplacements = np.array([])
        self.eigendisplacements = eigendisplacements
        if labels_dict is None:
            labels_dict = {}

        for q in qpoints:
            # let see if this qpoint has been assigned a label
            label = None
            for c in labels_dict:
                if np.linalg.norm(q - np.array(labels_dict[c])) < 0.0001:
                    label = c
                    self.labels_dict[label] = Kpoint(
                        q, lattice, label=label,
                        coords_are_cartesian=coords_are_cartesian)
            self.qpoints.append(
                Kpoint(q, lattice, label=label,
                       coords_are_cartesian=coords_are_cartesian))
        self.bands = frequencies
        self.nb_bands = len(self.bands)
        self.nb_qpoints = len(self.qpoints)

        # normalize directions for nac_frequencies and nac_eigendisplacements
        self.nac_frequencies = []
        self.nac_eigendisplacements = []
        if nac_frequencies is not None:
            for t in nac_frequencies:
                self.nac_frequencies.append(([i / np.linalg.norm(t[0]) for i in t[0]], t[1]))
        if nac_eigendisplacements is not None:
            for t in nac_eigendisplacements:
                self.nac_eigendisplacements.append(([i/np.linalg.norm(t[0]) for i in t[0]], t[1]))

    def min_freq(self):
        """
        Returns the point where the minimum frequency is reached and its value
        """
        i = np.unravel_index(np.argmin(self.bands), self.bands.shape)

        return self.qpoints[i[0]], self.bands[i]

    def has_imaginary_freq(self, tol=1e-5):
        """
        True if imaginary frequencies are present in the BS.
        """

        return self.min_freq()[1] + tol < 0

    @property
    def has_nac(self):
        """
        True if nac_frequencies are present.
        """
        return len(self.nac_frequencies) > 0

    @property
    def has_eigendisplacements(self):
        """
        True if eigendisplacements are present.
        """
        return len(self.eigendisplacements) > 0

    def get_nac_frequencies_along_dir(self, direction):
        """
        Returns the nac_frequencies for the given direction (not necessarily a versor).
        None if the direction is not present or nac_frequencies has not been calculated.

        Args:
            direction: the direction as a list of 3 elements
        Returns:
            the frequencies as a numpy array o(3*len(structure), len(qpoints)).
            None if not found.
        """
        versor = [i / np.linalg.norm(direction) for i in direction]
        for d, f in self.nac_frequencies:
            if np.allclose(versor, d):
                return f

        return None

    def get_nac_eigendisplacements_along_dir(self, direction):
        """
        Returns the nac_eigendisplacements for the given direction (not necessarily a versor).
        None if the direction is not present or nac_eigendisplacements has not been calculated.

        Args:
            direction: the direction as a list of 3 elements
        Returns:
            the eigendisplacements as a numpy array of complex numbers with shape
            (3*len(structure), len(structure), 3). None if not found.
        """
        versor = [i / np.linalg.norm(direction) for i in direction]
        for d, e in self.nac_eigendisplacements:
            if np.allclose(versor, d):
                return e

        return None

    def asr_breaking(self, tol_eigendisplacements=1e-5):
        """
        Returns the breaking of the acoustic sum rule for the three acoustic modes,
        if Gamma is present. None otherwise.
        If eigendisplacements are available they are used to determine the acoustic
        modes: selects the bands corresponding  to the eigendisplacements that
        represent to a translation within tol_eigendisplacements. If these are not
        identified or eigendisplacements are missing the first 3 modes will be used
        (indices [0:3]).
        """

        for i in range(self.nb_qpoints):
            if np.allclose(self.qpoints[i].frac_coords, (0, 0, 0)):

                if self.has_eigendisplacements:
                    acoustic_modes_index = []
                    for j in range(self.nb_bands):
                        eig = self.eigendisplacements[j][i]
                        if np.max(np.abs(eig[1:] - eig[:1])) < tol_eigendisplacements:
                            acoustic_modes_index.append(j)
                    # if acoustic modes are not correctly identified return use
                    # the first three modes
                    if len(acoustic_modes_index) != 3:
                        acoustic_modes_index = [0, 1, 2]
                    return self.bands[acoustic_modes_index, i]
                else:
                    return self.bands[:3, i]

        return None

    def as_dict(self):
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "lattice_rec": self.lattice_rec.as_dict(),
             "qpoints": []}
        #qpoints are not Kpoint objects dicts but are frac coords.Tthis makes
        #the dict smaller and avoids the repetition of the lattice
        for q in self.qpoints:
            d["qpoints"].append(q.as_dict()["fcoords"])
        d["bands"] = self.bands.tolist()
        d['labels_dict'] = {}
        for c in self.labels_dict:
            d['labels_dict'][c] = self.labels_dict[c].as_dict()['fcoords']
        # split the eigendisplacements to real and imaginary part for serialization
        d['eigendisplacements'] = dict(real=np.real(self.eigendisplacements).tolist(),
                                       imag=np.imag(self.eigendisplacements).tolist())
        d['nac_eigendisplacements'] = [(direction, dict(real=np.real(e).tolist(), imag=np.imag(e).tolist()))
                                       for direction, e in self.nac_eigendisplacements]
        d['nac_frequencies'] = [(direction, f.tolist()) for direction, f in self.nac_frequencies]

        if self.structure:
            d['structure'] = self.structure.as_dict()

        return d

    @classmethod
    def from_dict(cls, d):
        lattice_rec = Lattice(d['lattice_rec']['matrix'])
        eigendisplacements = np.array(d['eigendisplacements']['real']) + np.array(d['eigendisplacements']['imag'])*1j
        nac_eigendisplacements = [(direction, np.array(e['real']) + np.array(e['imag'])*1j)
                                  for direction, e in d['nac_eigendisplacements']]
        nac_frequencies = [(direction, np.array(f)) for direction, f in d['nac_frequencies']]
        structure = Structure.from_dict(d['structure']) if 'structure' in d else None
        return cls(d['qpoints'], np.array(d['bands']), lattice_rec, nac_frequencies, eigendisplacements,
                   nac_eigendisplacements, d['labels_dict'], structure=structure)


class PhononBandStructureSymmLine(PhononBandStructure):
    """
    This object stores phonon band structures along selected (symmetry) lines in the
    Brillouin zone. We call the different symmetry lines (ex: \Gamma to Z)
    "branches".
    """

    def __init__(self, qpoints, frequencies, lattice, has_nac=False, eigendisplacements=None,
                 labels_dict=None, coords_are_cartesian=False, structure=None):
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
                non-analytical corrections at Gamma. If True frequenciens at Gamma from
                diffent directions will be stored in naf. Default False.
            eigendisplacements: the phonon eigendisplacements associated to the
                frequencies in cartesian coordinates. A numpy array of complex
                numbers with shape (3*len(structure), len(qpoints), len(structure), 3).
                he First index of the array refers to the band, the second to the index
                of the qpoint, the third to the atom in the structure and the fourth
                to the cartesian coordinates.
            labels_dict: (dict) of {} this links a qpoint (in frac coords or
                cartesian coordinates depending on the coords) to a label.
            coords_are_cartesian: Whether the qpoint coordinates are cartesian.
            structure: The crystal structure (as a pymatgen Structure object)
                associated with the band structure. This is needed if we
                provide projections to the band structure
        """

        super(PhononBandStructureSymmLine, self).__init__(
            qpoints, frequencies, lattice, None, eigendisplacements,
            None, labels_dict, coords_are_cartesian, structure)

        self.distance = []
        self.branches = []
        one_group = []
        branches_tmp = []

        #get labels and distance for each qpoint
        previous_qpoint = self.qpoints[0]
        previous_distance = 0.0
        previous_label = self.qpoints[0].label
        for i in range(self.nb_qpoints):
            label = self.qpoints[i].label
            if label is not None and previous_label is not None:
                self.distance.append(previous_distance)
            else:
                self.distance.append(
                    np.linalg.norm(self.qpoints[i].cart_coords -
                                   previous_qpoint.cart_coords) +
                    previous_distance)
            previous_qpoint = self.qpoints[i]
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
                "name": str(self.qpoints[b[0]].label) + "-" +
                        str(self.qpoints[b[-1]].label)})

        # extract the frequencies with non-analytical contribution at gamma
        if has_nac:
            naf = []
            nac_eigendisplacements = []
            for i in range(self.nb_qpoints):
                # get directions with nac irrespectively of the label_dict. NB: with labels
                # the gamma point is expected to appear twice consecutively.
                if np.allclose(qpoints[i], (0, 0, 0)):
                    if i>0 and not np.allclose(qpoints[i-1], (0, 0, 0)):
                        q_dir = self.qpoints[i - 1]
                        direction = [q_dir.frac_coords / np.linalg.norm(q_dir.frac_coords)]
                        naf.append((direction, frequencies[:, i]))
                        nac_eigendisplacements.append((direction, eigendisplacements[:, i]))
                    if i<len(frequencies)-1 and not np.allclose(qpoints[i+1], (0, 0, 0)):
                        q_dir = self.qpoints[i + 1]
                        direction = [q_dir.frac_coords / np.linalg.norm(q_dir.frac_coords)]
                        naf.append((direction, frequencies[:, i]))
                        nac_eigendisplacements.append((direction, eigendisplacements[:, i]))

            self.nac_frequencies = np.array(naf)
            self.nac_eigendisplacements = np.array(nac_eigendisplacements)

    def get_equivalent_qpoints(self, index):
        """
        Returns the list of qpoint indices equivalent (meaning they are the
        same frac coords) to the given one.

        Args:
            index: the qpoint index

        Returns:
            a list of equivalent indices

        TODO: now it uses the label we might want to use coordinates instead
        (in case there was a mislabel)
        """
        #if the qpoint has no label it can"t have a repetition along the band
        #structure line object

        if self.qpoints[index].label is None:
            return [index]

        list_index_qpoints = []
        for i in range(self.nb_qpoints):
            if self.qpoints[i].label == self.qpoints[index].label:
                list_index_qpoints.append(i)

        return list_index_qpoints

    def get_branch(self, index):
        """
        Returns in what branch(es) is the qpoint. There can be several
        branches.

        Args:
            index: the qpoint index

        Returns:
            A list of dictionaries [{"name","start_index","end_index","index"}]
            indicating all branches in which the qpoint is. It takes into
            account the fact that one qpoint (e.g., \Gamma) can be in several
            branches
        """
        to_return = []
        for i in self.get_equivalent_qpoints(index):
            for b in self.branches:
                    if b["start_index"] <= i <= b["end_index"]:
                        to_return.append({"name": b["name"],
                                          "start_index": b["start_index"],
                                          "end_index": b["end_index"],
                                          "index": i})
        return to_return

    def as_dict(self):
        d = super(PhononBandStructureSymmLine, self).as_dict()
        # remove nac_frequencies and nac_eigendisplacements as they are reconstructed
        # in the __init__ when the dict is deserialized
        nac_frequencies = d.pop('nac_frequencies')
        d.pop('nac_eigendisplacements')
        d['has_nac'] = len(nac_frequencies) > 0
        return d

    @classmethod
    def from_dict(cls, d):
        lattice_rec = Lattice(d['lattice_rec']['matrix'])
        eigendisplacements = np.array(d['eigendisplacements']['real']) + np.array(d['eigendisplacements']['imag'])*1j
        structure = Structure.from_dict(d['structure']) if 'structure' in d else None
        return cls(d['qpoints'], np.array(d['bands']), lattice_rec, d['has_nac'], eigendisplacements,
                   d['labels_dict'], structure=structure)
