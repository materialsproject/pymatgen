# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import numpy as np
import math
import itertools
import collections
import json

from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.electronic_structure.bandstructure import Kpoint
from monty.json import MSONable, MontyDecoder, MontyEncoder

"""
This module provides classes to define a phonon band structure.
"""

class PhBandStructure(MSONable):
    """
    This is the most generic phonon band structure data possible
    it's defined by a list of kpoints + frequencies for each of them.
    Additional information may be given for frequencies at Gamma, where
    non-analytical contribution may be taken into account.
    """
    #TODO define  units for frequency
    def __init__(self, kpoints, frequencies, lattice, nac_frequencies=None, eigendisplacements=None,
                 nac_eigendisplacements=None, labels_dict=None, coords_are_cartesian=False,
                 structure=None):
        """
        Args:
            kpoints: list of kpoint as numpy arrays, in frac_coords of the
                given lattice by default
            frequencies: list of phonon frequencies in eV as a numpy array with shape
                (len(kpoints), 3*len(structure))
            lattice: The reciprocal lattice as a pymatgen Lattice object.
                Pymatgen uses the physics convention of reciprocal lattice vectors
                WITH a 2*pi coefficient.
            nac_frequencies: Frequencies with non-analytical contributions at Gamma.
                A list of tuples. The first element of each tuple should be a list
                defining the direction (not necessarily a versor, will be normalized
                internally). The second element containing the 3*len(structure)
                phonon frequecies with non-analytical correction for that direction.
            eigendisplacements: the phonon eigendisplacements associated to the
                frequencies in cartesian coordinates. A numpy array of complex
                numbers with shape (len(kpoints), 3*len(structure), 3).
            nac_eigendisplacements: the phonon eigendisplacements associated to the
                non-analytical frequencies in nac_frequencies in cartesian coordinates.
                A list of tuples. The first element of each tuple should be a list
                defining the direction. The second element containing a numpy array of
                complex numbers with shape (3*len(structure), len(structure), 3).
            labels_dict: (dict) of {} this links a kpoint (in frac coords or
                cartesian coordinates depending on the coords) to a label.
            coords_are_cartesian: Whether the k-point coordinates are cartesian.
            structure: The crystal structure (as a pymatgen Structure object)
                associated with the band structure. This is needed if we
                provide projections to the band structure
        """
        self.lattice_rec = lattice
        self.kpoints = []
        self.labels_dict = {}
        self.structure = structure
        if eigendisplacements is None:
            eigendisplacements = []
        self.eigendisplacements = eigendisplacements
        if labels_dict is None:
            labels_dict = {}

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
        self.bands = frequencies
        self.nb_bands = len(self.bands)
        self.nb_kpoints = len(self.kpoints)

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
        i = np.argmin(self.bands)

        return self.kpoints[i[0]], self.bands[i]

    def has_imaginary_freq(self, tol=1e-5):
        """
        True imaginary frequencies are present in the BS.
        """

        return self.min_freq()[1] + tol < 0

    def asr_breaking(self):
        """
        Returns the breaking of the acoustic sum rule for the three acoustic modes,
        if Gamma is present. None otherwise.
        if eigendisplacements are available they are used to determine the acoustic modes,
        otherwise the first 3 modes will be used.
        """
        for i in range(self.nb_kpoints):
            if np.allclose(self.kpoints[i].frac_coords, (0, 0, 0)):

                if self.eigendisplacements:
                    acoustic_modes_index = []
                    for j in range(self.nb_bands):
                        eig = self.eigendisplacements[i][j]
                        if max(np.abs(eig[1:] - eig[:1])) < 1e-8:
                            acoustic_modes_index.append(j)
                    return self.bands[acoustic_modes_index, i]
                else:
                    return self.bands[:3, i]

        return None

    def as_dict(self):
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "lattice_rec": self.lattice_rec.as_dict(),
             "kpoints": []}
        #kpoints are not kpoint objects dicts but are frac coords (this makes
        #the dict smaller and avoids the repetition of the lattice
        for k in self.kpoints:
            d["kpoints"].append(k.as_dict()["fcoords"])
        d["bands"] = self.bands
        d['labels_dict'] = {}
        for c in self.labels_dict:
            d['labels_dict'][c] = self.labels_dict[c].frac_coords
        # split the eigendisplacements to real and imaginary part for serialization
        d['eigendisplacements'] = dict(real=np.real(self.eigendisplacements),
                                       imag=np.imag(self.eigendisplacements))
        d['nac_eigendisplacements'] = [(direction, dict(real=np.real(e), imag=np.imag(e)))
                                       for direction, e in self.nac_eigendisplacements]
        d['nac_frequencies'] = self.nac_frequencies

        if self.structure:
            d['structure'] = self.structure.as_dict()

        # converts numpy arrays
        return json.loads(json.dumps(d, cls=MontyEncoder))

    @classmethod
    def from_dict(cls, d):
        dec = MontyDecoder()
        d = dec.process_decoded(d)
        lattice_rec = d["lattice_rec"]
        eigendisplacements = d['eigendisplacements']['real'] + d['eigendisplacements']['imag']*1j
        nac_eigendisplacements = [(direction, e['real']+ e['imag']*1j)
                                  for direction, e in d['nac_eigendisplacements']]
        structure = d['structure'] if 'structure' in d else None
        return cls(d['kpoints'], d['bands'], lattice_rec, d['nac_frequencies'], eigendisplacements,
                   nac_eigendisplacements, d['labels_dict'], structure=structure)


class PhBandStructureSymmLine(PhBandStructure):
    """
    This object stores phonon band structures along selected (symmetry) lines in the
    Brillouin zone. We call the different symmetry lines (ex: \Gamma to Z)
    "branches".
    """

    def __init__(self, kpoints, frequencies, lattice, has_nac=False, eigendisplacements=None,
                 labels_dict=None, coords_are_cartesian=False, structure=None):
        """
        Args:
            kpoints: list of kpoint as numpy arrays, in frac_coords of the
                given lattice by default
            frequencies: list of phonon frequencies in eV as a numpy array with shape
                (len(kpoints), 3*len(structure))
            lattice: The reciprocal lattice as a pymatgen Lattice object.
                Pymatgen uses the physics convention of reciprocal lattice vectors
                WITH a 2*pi coefficient
            has_nac: specify if the band structure has been produced taking into account
                non-analytical corrections at Gamma. If True frequenciens at Gamma from
                diffent directions will be stored in naf. Default False.
            eigendisplacements: the phonon eigendisplacements associated to the
                frequencies in cartesian coordinates. A numpy array of complex
                numbers with shape (len(kpoints), 3*len(structure), 3).
            labels_dict: (dict) of {} this links a kpoint (in frac coords or
                cartesian coordinates depending on the coords) to a label.
            coords_are_cartesian: Whether the k-point coordinates are cartesian.
            structure: The crystal structure (as a pymatgen Structure object)
                associated with the band structure. This is needed if we
                provide projections to the band structure
        """

        super(PhBandStructureSymmLine, self).__init__(
            kpoints, frequencies, lattice, None, eigendisplacements,
            None, labels_dict, coords_are_cartesian, structure)

        self.distance = []
        self.branches = []
        one_group = []
        branches_tmp = []

        #get labels and distance for each kpoint
        previous_kpoint = self.kpoints[0]
        previous_distance = 0.0
        previous_label = self.kpoints[0].label
        for i in range(self.nb_kpoints):
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

        # extract the frequencies with non-analytical contribution at gamma
        if has_nac:
            naf = []
            nac_eigendisplacements = []
            for i in range(self.nb_kpoints):
                # get directions with nac irrespectively of the label_dict. nb: with labels
                # the gamma point is expected to appear twice consecutively.
                if np.allclose(kpoints[i], (0, 0, 0)):
                    if i>0 and not np.allclose(kpoints[i-1], (0, 0, 0)):
                        k_dir = self.kpoints[i-1]
                        direction = [k_dir.frac_coords / np.linalg.norm(k_dir.frac_coords)]
                        naf.append((direction, frequencies[:, i]))
                        nac_eigendisplacements.append((direction, eigendisplacements[i]))
                    if i<len(frequencies)-1 and not np.allclose(kpoints[i+1], (0, 0, 0)):
                        k_dir = self.kpoints[i+1]
                        direction = [k_dir.frac_coords / np.linalg.norm(k_dir.frac_coords)]
                        naf.append((direction, frequencies[:, i]))
                        nac_eigendisplacements.append((direction, eigendisplacements[i]))

            self.nac_frequencies = np.array(naf)
            self.nac_eigendisplacements = np.array(nac_eigendisplacements)

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
        #if the kpoint has no label it can"t have a repetition along the band
        #structure line object

        if self.kpoints[index].label is None:
            return [index]

        list_index_kpoints = []
        for i in range(self.nb_kpoints):
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

    def as_dict(self):
        d = super(PhBandStructureSymmLine, self).as_dict()
        # remove nac_frequencies and nac_eigendisplacements as they are reconstructed
        # in the __init__ when the dict is deserialized
        nac_frequencies = d.pop('nac_frequencies')
        d.pop('nac_eigendisplacements')
        d['has_nac'] = len(nac_frequencies) > 0
        return d

    @classmethod
    def from_dict(cls, d):
        dec = MontyDecoder()
        d = dec.process_decoded(d)
        lattice_rec = d["lattice_rec"]
        eigendisplacements = d['eigendisplacements']['real'] + d['eigendisplacements']['imag']*1j
        structure = d['structure'] if 'structure' in d else None
        return cls(d['kpoints'], d['bands'], lattice_rec, d['has_nac'], eigendisplacements,
                   d['labels_dict'], structure=structure)

