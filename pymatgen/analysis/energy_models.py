# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements a EnergyModel abstract class and some basic
implementations. Basically, an EnergyModel is any model that returns an
"energy" for any given structure.
"""

import abc

from monty.json import MSONable

from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

__version__ = "0.1"


class EnergyModel(MSONable, metaclass=abc.ABCMeta):
    """
    Abstract structure filter class.
    """

    @abc.abstractmethod
    def get_energy(self, structure) -> float:
        """
        :param structure: Structure
        :return: Energy value
        """
        return 0.0

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation
        :return: EnergyModel
        """
        return cls(**d["init_args"])


class EwaldElectrostaticModel(EnergyModel):
    """
    Wrapper around EwaldSum to calculate the electrostatic energy.
    """

    def __init__(self, real_space_cut=None, recip_space_cut=None, eta=None, acc_factor=8.0):
        """
        Initializes the model. Args have the same definitions as in
        :class:`pymatgen.analysis.ewald.EwaldSummation`.

        Args:
            real_space_cut (float): Real space cutoff radius dictating how
                many terms are used in the real space sum. Defaults to None,
                which means determine automagically using the formula given
                in gulp 3.1 documentation.
            recip_space_cut (float): Reciprocal space cutoff radius.
                Defaults to None, which means determine automagically using
                the formula given in gulp 3.1 documentation.
            eta (float): Screening parameter. Defaults to None, which means
                determine automatically.
            acc_factor (float): No. of significant figures each sum is
                converged to.
        """
        self.real_space_cut = real_space_cut
        self.recip_space_cut = recip_space_cut
        self.eta = eta
        self.acc_factor = acc_factor

    def get_energy(self, structure):
        """
        :param structure: Structure
        :return: Energy value
        """
        e = EwaldSummation(
            structure,
            real_space_cut=self.real_space_cut,
            recip_space_cut=self.recip_space_cut,
            eta=self.eta,
            acc_factor=self.acc_factor,
        )
        return e.total_energy

    def as_dict(self):
        """
        :return: MSONable dict
        """
        return {
            "version": __version__,
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "init_args": {
                "real_space_cut": self.real_space_cut,
                "recip_space_cut": self.recip_space_cut,
                "eta": self.eta,
                "acc_factor": self.acc_factor,
            },
        }


class SymmetryModel(EnergyModel):
    """
    Sets the energy to the -ve of the spacegroup number. Higher symmetry =>
    lower "energy".

    Args have same meaning as in
    :class:`pymatgen.symmetry.finder.SpacegroupAnalyzer`.
    """

    def __init__(self, symprec=0.1, angle_tolerance=5):
        """
        Args:
            symprec (float): Symmetry tolerance. Defaults to 0.1.
            angle_tolerance (float): Tolerance for angles. Defaults to 5 degrees.
        """
        self.symprec = symprec
        self.angle_tolerance = angle_tolerance

    def get_energy(self, structure):
        """
        :param structure: Structure
        :return: Energy value
        """
        f = SpacegroupAnalyzer(structure, symprec=self.symprec, angle_tolerance=self.angle_tolerance)
        return -f.get_space_group_number()

    def as_dict(self):
        """
        :return: MSONable dict
        """
        return {
            "version": __version__,
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "init_args": {
                "symprec": self.symprec,
                "angle_tolerance": self.angle_tolerance,
            },
        }


class IsingModel(EnergyModel):
    """
    A very simple Ising model, with r^2 decay.
    """

    def __init__(self, j, max_radius):
        """
        Args:
            j (float): The interaction parameter. E = J * spin1 * spin2.
            radius (float): max_radius for the interaction.
        """
        self.j = j
        self.max_radius = max_radius

    def get_energy(self, structure):
        """
        :param structure: Structure
        :return: Energy value
        """
        all_nn = structure.get_all_neighbors(r=self.max_radius)
        energy = 0
        for i, nns in enumerate(all_nn):
            s1 = getattr(structure[i].specie, "spin", 0)
            for nn in nns:
                energy += self.j * s1 * getattr(nn.specie, "spin", 0) / (nn.nn_distance ** 2)
        return energy

    def as_dict(self):
        """
        :return: MSONable dict
        """
        return {
            "version": __version__,
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "init_args": {"j": self.j, "max_radius": self.max_radius},
        }


class NsitesModel(EnergyModel):
    """
    Sets the energy to the number of sites. More sites => higher "energy".
    Used to rank structures from smallest number of sites to largest number
    of sites after enumeration.
    """

    def get_energy(self, structure):
        """
        :param structure: Structure
        :return: Energy value
        """
        return len(structure)

    def as_dict(self):
        """
        :return: MSONable dict
        """
        return {
            "version": __version__,
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "init_args": {},
        }
