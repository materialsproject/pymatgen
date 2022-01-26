# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Entries are containers for calculated information, which is used in
many analyses. This module contains entry related tools and implements
the base Entry class, which is the basic entity that can be used to
store calculated information. Other Entry classes such as ComputedEntry
and PDEntry inherit from this class.
"""

from abc import ABCMeta, abstractmethod
from typing import Dict, Literal, Union

import numpy as np
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
        composition: Union[Composition, str, Dict[str, float]],
        energy: float,
    ):
        """
        Initializes an Entry.

        Args:
            composition (Composition): Composition of the entry. For
                flexibility, this can take the form of all the typical input
                taken by a Composition, including a {symbol: amt} dict,
                a string formula, and others.
            energy (float): Energy of the entry.
        """
        self._composition = Composition(composition)
        self._energy = energy

    @property
    def is_element(self) -> bool:
        """
        :return: Whether composition of entry is an element.
        """
        # NOTE _composition rather than composition as GrandPDEntry
        # edge case exists if we have a compound where chempots are
        # given for all bar one element type
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

    def __repr__(self):
        return f"{self.__class__.__name__} : {self.composition} with energy = {self.energy:.4f}"

    def __str__(self):
        return self.__repr__()

    def normalize(self, mode: Literal["formula_unit", "atom"] = "formula_unit") -> "Entry":
        """
        Normalize the entry's composition and energy.

        Args:
            mode ("formula_unit" | "atom"): "formula_unit" (the default) normalizes to composition.reduced_formula.
                "atom" normalizes such that the composition amounts sum to 1.
        """

        factor = self._normalization_factor(mode)
        new_composition = self._composition / factor
        new_energy = self._energy / factor

        new_entry_dict = self.as_dict()
        new_entry_dict["composition"] = new_composition.as_dict()
        new_entry_dict["energy"] = new_energy

        return self.from_dict(new_entry_dict)

    def _normalization_factor(self, mode: Literal["formula_unit", "atom"] = "formula_unit") -> float:
        # NOTE here we use composition rather than _composition in order to ensure
        # that we have the expected behavior downstream in cases where composition
        # is overwritten (GrandPotPDEntry, TransformedPDEntry)
        if mode == "atom":
            factor = self.composition.num_atoms
        elif mode == "formula_unit":
            factor = self.composition.get_reduced_composition_and_factor()[1]
        else:
            raise ValueError(f"{mode} is not an allowed option for normalization")

        return factor

    def as_dict(self) -> dict:
        """
        :return: MSONable dict.
        """
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "energy": self._energy,
            "composition": self._composition.as_dict(),
        }

    def __eq__(self, other):
        # NOTE: Scaled duplicates i.e. physically equivalent materials
        # are not equal unless normalized separately.
        if self is other:
            return True

        # Equality is defined based on composition and energy
        # If structures are involved, it is assumed that a {composition, energy} is
        # vanishingly unlikely to be the same if the structures are different

        if not np.allclose(self.energy, other.energy):
            return False

        return self.composition == other.composition

    def __hash__(self):
        # NOTE truncate _energy to 8 dp to ensure same robustness
        # as np.allclose
        return hash(f"{self.__class__.__name__}{self._composition.formula}{self._energy:.8f}")
