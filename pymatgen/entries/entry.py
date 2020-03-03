# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements the base Entry class, which is the basic entity
that can be used to store calculated information. Other Entry classes 
such as ComputedEntry and PDEntry inherit from this class.
"""
from pymatgen.core.composition import Composition
from monty.json import MSONable
from abc import ABCMeta, abstractmethod


__author__ = "Shyue Ping Ong, Anubhav Jain, Ayush Gupta"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Mar 02, 2020"


class Entry(MSONable, metaclass=ABCMeta):
    """
    A lightweight Entry object containing various types 
    of data associated with a specific chemical composition.

    """

    def __init__(self,
                 composition: Composition,
                 energy: float):
        """
        Initializes an Entry.

        Args:
            composition (Composition): Composition of the entry. For
                flexibility, this can take the form of all the typical input
                taken by a Composition, including a {symbol: amt} dict,
                a string formula, and others.
            energy (float): Energy of the entry.
            correction (float): A correction to be applied to the energy.
                This is used to modify the energy for certain analyses.
                Defaults to 0.0.
        """
        self._energy = energy
        self.composition = Composition(composition)

    @property
    def is_element(self) -> bool:
        """
        :return: Whether composition of entry is an element.
        """
        return self.composition.is_element

    @property
    @abstractmethod
    def energy(self) -> float:
        """
        :return: the energy of the entry. Each subclass must define its own energy property.
        """

    @property
    def energy_per_atom(self) -> float:
        """
        :return: the *corrected* energy per atom of the entry.
        """
        return self.energy / self.composition.num_atoms

    def __str__(self):
        return self.__repr__()

    def normalize(self, mode: str = "formula_unit") -> None:
        """
        Normalize the entry's composition, energy and any corrections.

        Args:
            mode: "formula_unit" is the default, which normalizes to
                composition.reduced_formula. The other option is "atom", which
                normalizes such that the composition amounts sum to 1.
        """
        factor = self.normalization_factor(mode)
        self.composition /= factor
        self._energy /= factor

    def normalization_factor(self, mode: str = "formula_unit") -> float:
        if mode == "atom":
            factor = self.composition.num_atoms
        else:
            comp, factor = self.composition.get_reduced_composition_and_factor()
        return factor

    def as_dict(self) -> dict:
        """
        :return: MSONable dict.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "energy": self._energy,
                "composition": self.composition.as_dict()}
