# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import numpy as np

import six

from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.core.structure import Structure
from pymatgen.util.coord_utils import get_linear_interpolated_value
from monty.json import MSONable

"""
This module defines classes to represent the phonon density of states, etc.
"""

class PhDos(MSONable):
    """
    Basic DOS object. All other DOS objects are extended versions of this
    object.

    Args:
        frequencies: A sequences of frequencies
        densities A list representing the density of states.

    .. attribute: frequencies

        The sequence of frequencies

    .. attribute: densities

        A list of phonon densities

    """

    def __init__(self, frequencies, densities):
        self.frequencies = np.array(frequencies)
        self.densities = np.array(densities)

    def get_smeared_densities(self, sigma):
        """
        Returns the densities, but with a Gaussian smearing of
        std dev sigma applied.

        Args:
            sigma: Std dev of Gaussian smearing function.

        Returns:
            Gaussian-smeared densities.
        """

        from scipy.ndimage.filters import gaussian_filter1d
        diff = [self.frequencies[i + 1] - self.frequencies[i]
                for i in range(len(self.frequencies) - 1)]
        avgdiff = sum(diff) / len(diff)

        smeared_dens = gaussian_filter1d(self.densities, sigma / avgdiff)
        return smeared_dens

    def __add__(self, other):
        """
        Adds two DOS together. Checks that frequency scales are the same.
        Otherwise, a ValueError is thrown.

        Args:
            other: Another DOS object.

        Returns:
            Sum of the two DOSs.
        """
        if not all(np.equal(self.frequencies, other.frequencies)):
            raise ValueError("Frequencies of both DOS are not compatible!")
        densities = self.frequencies + other.frequencies
        return Dos(self.frequencies, densities)

    def __radd__(self, other):
        """
        Reflected addition of two DOS objects

        Args:
            other: Another DOS object.

        Returns:
            Sum of the two DOSs.
        """

        return self.__add__(other)

    def get_interpolated_value(self, frequency):
        """
        Returns interpolated density for a particular frequency.

        Args:
            frequency: frequency to return the density for.
        """
        return  get_linear_interpolated_value(self.frequencies,
                                              self.densities, frequency)

    def __str__(self):
        """
        Returns a string which can be easily plotted (using gnuplot).
        """
        stringarray = ["#{:30s} {:30s}".format("Frequency", "Density")]
        for i, frequency in enumerate(self.frequencies):
            stringarray.append("{:.5f} {:.5f}"
                               .format(frequency, self.densities[i]))
        return "\n".join(stringarray)

    @classmethod
    def from_dict(cls, d):
        """
        Returns PhDos object from dict representation of PhDos.
        """
        return cls(d["freuqencies"],d["densities"])

    def as_dict(self):
        """
        Json-serializable dict representation of PhDos.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "frequency": list(self.frequencies),
                "densities": self.densities}


class CompletePhDos(PhDos):
    """
    This wrapper class defines a total dos, and also provides a list of PDos.

    Args:
        structure: Structure associated with this particular DOS.
        total_dos: total Dos for structure
        pdoss: The pdoss are supplied as an {Site: Densities}

    .. attribute:: structure

        Structure associated with the CompleteDos.

    .. attribute:: pdos

        Dict of partial densities of the form {Site:Densities}
    """
    def __init__(self, structure, total_dos, pdoss):
        super(CompletePhDos, self).__init__(
            frequencies=total_dos.frequencies, densities=total_dos.densities)
        self.pdos = pdoss
        self.structure = structure

    def get_element_dos(self):
        """
        Get element projected Dos.

        Returns:
            dict of {Element: Dos}
        """

        el_dos = {}
        for site, atom_dos in self.pdos.items():
            el = site.specie
            if el not in el_dos:
                el_dos[el] = atom_dos
            else:
                el_dos[el] += add_densities(el_dos[el], pdos)
        return {el: PhDos(self.frequencies, densities)
                for el, densities in el_dos.items()}

    @classmethod
    def from_dict(cls, d):
        """
        Returns CompleteDos object from dict representation.
        """
        tdos = Dos.from_dict(d)
        struct = Structure.from_dict(d["structure"])
        pdoss = {}
        for at, pdos in zip(struct, d["pdos"]):
            pdoss[at] = pdoss
            pdoss[at] = orb_dos

        return cls(struct, tdos, pdoss)

    def as_dict(self):
        """
        Json-serializable dict representation of CompletePhDos.
        """
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "structure": self.structure.as_dict(),
             "frequencies": list(self.frequencies),
             "densities": list(self.densities),
             "pdos": []}
        if len(self.pdos) > 0:
            l = []
            for at in self.structure:
                l.append(self.pdos[at])
            d["pdos"].append(l)
        return d

    def __str__(self):
        return "Complete DOS for " + str(self.structure)
