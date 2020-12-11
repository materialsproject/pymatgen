# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Entries are containers for calculated information, which is used in
many analyses. This module contains entry related tools and implements
the base Entry class, which is the basic entity that can be used to
store calculated information. Other Entry classes such as ComputedEntry
and PDEntry inherit from this class.
"""

import copy
import hashlib
from abc import ABCMeta, abstractmethod
from typing import Optional

from monty.json import MSONable

from pymatgen.core.composition import Composition

__author__ = "Shyue Ping Ong, Anubhav Jain, Ayush Gupta"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Mar 03, 2020"


class Entry(MSONable, metaclass=ABCMeta):
    """
    A lightweight object containing the energy associated with
    a specific chemical composition. This base class is not
    intended to be instantiated directly. Note that classes
    which inherit from Entry must define a .energy property.

    """

    def __init__(
        self,
        composition: Composition,
        energy: float,
        per_atom: bool = False,
    ):
        """
        Initializes an Entry.

        Args:
            composition (Composition): Composition of the entry. For
                flexibility, this can take the form of all the typical input
                taken by a Composition, including a {symbol: amt} dict,
                a string formula, and others.
            energy (float): Energy of the entry.
            per_atom (bool): Whether the energy given is per atom.
        """
        self._composition = Composition(composition)
        if per_atom:
            self._energy = energy * self.composition.num_atoms
        else:
            self._energy = energy

    @property
    def is_element(self) -> bool:
        """
        :return: Whether composition of entry is an element.
        """
        return self._composition.is_element

    @property
    def composition(self) -> Composition:
        """
        :return: the composition of the entry.
        """
        return self._composition

    @property
    @abstractmethod
    def energy(self) -> float:
        """
        :return: the energy of the entry.
        """

    @property
    def energy_per_atom(self) -> float:
        """
        :return: the energy per atom of the entry.
        """
        return self.energy / self.composition.num_atoms

    def __str__(self):
        return self.__repr__()

    def normalize(
        self, mode: str = "formula_unit", inplace: bool = True
    ) -> Optional["Entry"]:
        """
        Normalize the entry's composition and energy.

        Args:
            mode: "formula_unit" is the default, which normalizes to
                composition.reduced_formula. The other option is "atom", which
                normalizes such that the composition amounts sum to 1.
            inplace: "True" is the default which normalises the current Entry object.
                Setting inplace to "False" returns a normalized copy of the Entry object.
        """
        if inplace:
            factor = self._normalization_factor(mode)
            self._composition /= factor
            self._energy /= factor
            return None

        entry = copy.deepcopy(self)
        entry.normalize(mode, inplace=True)
        return entry

    def _normalization_factor(self, mode: str = "formula_unit") -> float:
        # NOTE here we use composition rather than _composition in order to ensure
        # that we have the expected behaviour downstream in cases where composition
        # is overwritten (GrandPotPDEntry, TransformedPDEntry)
        if mode == "atom":
            factor = self.composition.num_atoms
        elif mode == "formula_unit":
            factor = self.composition.get_reduced_composition_and_factor()[1]
        else:
            raise ValueError(
                "`{}` is not an allowed option for normalization".format(mode)
            )

        return factor

    def as_dict(self) -> dict:
        """
        :return: MSONable dict.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "energy": self._energy,
                "composition": self.composition.as_dict()}

    def __eq__(self, other):
        # NOTE Scaled duplicates i.e. physically equivalent materials
        # are not equal unless normalized separately
        if self is other:
            return True

        if isinstance(other, self.__class__):
            # NOTE this is not performant and should be overwritten or
            # preceeded if faster rigorous checks are available.
            return self.as_dict() == other.as_dict()
        return False

    def __hash__(self):
        data_md5 = hashlib.md5((
            f"{self.__class__.__name__}"
            f"{self._composition.reduced_formula}"
            f"{self._energy}").encode('utf-8')
        ).hexdigest()
        return int(data_md5, 16)