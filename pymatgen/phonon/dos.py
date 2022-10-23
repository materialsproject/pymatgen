# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module defines classes to represent the phonon density of states, etc.
"""

import numpy as np
import scipy.constants as const
from monty.functools import lazy_property
from monty.json import MSONable

from pymatgen.core.structure import Structure
from pymatgen.util.coord import get_linear_interpolated_value

BOLTZ_THZ_PER_K = const.value("Boltzmann constant in Hz/K") / const.tera  # Boltzmann constant in THz/K
THZ_TO_J = const.value("hertz-joule relationship") * const.tera


def coth(x):
    """
    Coth function.

    Args:
        x (): value

    Returns:
        coth(x)
    """
    return 1.0 / np.tanh(x)


class PhononDos(MSONable):
    """
    Basic DOS object. All other DOS objects are extended versions of this
    object.
    """

    def __init__(self, frequencies, densities):
        """
        Args:
            frequencies: A sequences of frequencies in THz
            densities: A list representing the density of states.
        """
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

        diff = [self.frequencies[i + 1] - self.frequencies[i] for i in range(len(self.frequencies) - 1)]
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
        densities = self.densities + other.densities
        return PhononDos(self.frequencies, densities)

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
        return get_linear_interpolated_value(self.frequencies, self.densities, frequency)

    def __str__(self):
        """
        Returns a string which can be easily plotted (using gnuplot).
        """
        stringarray = [f"#{'Frequency':30s} {'Density':30s}"]
        for i, frequency in enumerate(self.frequencies):
            stringarray.append(f"{frequency:.5f} {self.densities[i]:.5f}")
        return "\n".join(stringarray)

    @classmethod
    def from_dict(cls, d):
        """
        Returns PhononDos object from dict representation of PhononDos.
        """
        return cls(d["frequencies"], d["densities"])

    def as_dict(self):
        """
        JSON-serializable dict representation of PhononDos.
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "frequencies": list(self.frequencies),
            "densities": list(self.densities),
        }

    @lazy_property
    def ind_zero_freq(self):
        """
        Index of the first point for which the frequencies are equal or greater than zero.
        """
        ind = np.searchsorted(self.frequencies, 0)
        if ind >= len(self.frequencies):
            raise ValueError("No positive frequencies found")
        return ind

    @lazy_property
    def _positive_frequencies(self):
        """
        Numpy array containing the list of positive frequencies
        """
        return self.frequencies[self.ind_zero_freq :]

    @lazy_property
    def _positive_densities(self):
        """
        Numpy array containing the list of densities corresponding to positive frequencies
        """
        return self.densities[self.ind_zero_freq :]

    def cv(self, t, structure=None):
        """
        Constant volume specific heat C_v at temperature T obtained from the integration of the DOS.
        Only positive frequencies will be used.
        Result in J/(K*mol-c). A mol-c is the abbreviation of a mole-cell, that is, the number
        of Avogadro times the atoms in a unit cell. To compare with experimental data the result
        should be divided by the number of unit formulas in the cell. If the structure is provided
        the division is performed internally and the result is in J/(K*mol)

        Args:
            t: a temperature in K
            structure: the structure of the system. If not None it will be used to determine the number of
                formula units
        Returns:
            Constant volume specific heat C_v
        """

        if t == 0:
            return 0

        freqs = self._positive_frequencies
        dens = self._positive_densities

        def csch2(x):
            return 1.0 / (np.sinh(x) ** 2)

        wd2kt = freqs / (2 * BOLTZ_THZ_PER_K * t)
        cv = np.trapz(wd2kt**2 * csch2(wd2kt) * dens, x=freqs)
        cv *= const.Boltzmann * const.Avogadro

        if structure:
            formula_units = structure.composition.num_atoms / structure.composition.reduced_composition.num_atoms
            cv /= formula_units

        return cv

    def entropy(self, t, structure=None):
        """
        Vibrational entropy at temperature T obtained from the integration of the DOS.
        Only positive frequencies will be used.
        Result in J/(K*mol-c). A mol-c is the abbreviation of a mole-cell, that is, the number
        of Avogadro times the atoms in a unit cell. To compare with experimental data the result
        should be divided by the number of unit formulas in the cell. If the structure is provided
        the division is performed internally and the result is in J/(K*mol)

        Args:
            t: a temperature in K
            structure: the structure of the system. If not None it will be used to determine the number of
                formula units
        Returns:
            Vibrational entropy
        """

        if t == 0:
            return 0

        freqs = self._positive_frequencies
        dens = self._positive_densities

        wd2kt = freqs / (2 * BOLTZ_THZ_PER_K * t)
        s = np.trapz((wd2kt * coth(wd2kt) - np.log(2 * np.sinh(wd2kt))) * dens, x=freqs)

        s *= const.Boltzmann * const.Avogadro

        if structure:
            formula_units = structure.composition.num_atoms / structure.composition.reduced_composition.num_atoms
            s /= formula_units

        return s

    def internal_energy(self, t, structure=None):
        """
        Phonon contribution to the internal energy at temperature T obtained from the integration of the DOS.
        Only positive frequencies will be used.
        Result in J/mol-c. A mol-c is the abbreviation of a mole-cell, that is, the number
        of Avogadro times the atoms in a unit cell. To compare with experimental data the result
        should be divided by the number of unit formulas in the cell. If the structure is provided
        the division is performed internally and the result is in J/mol

        Args:
            t: a temperature in K
            structure: the structure of the system. If not None it will be used to determine the number of
                formula units
        Returns:
            Phonon contribution to the internal energy
        """

        if t == 0:
            return self.zero_point_energy(structure=structure)

        freqs = self._positive_frequencies
        dens = self._positive_densities

        wd2kt = freqs / (2 * BOLTZ_THZ_PER_K * t)
        e = np.trapz(freqs * coth(wd2kt) * dens, x=freqs) / 2

        e *= THZ_TO_J * const.Avogadro

        if structure:
            formula_units = structure.composition.num_atoms / structure.composition.reduced_composition.num_atoms
            e /= formula_units

        return e

    def helmholtz_free_energy(self, t, structure=None):
        """
        Phonon contribution to the Helmholtz free energy at temperature T obtained from the integration of the DOS.
        Only positive frequencies will be used.
        Result in J/mol-c. A mol-c is the abbreviation of a mole-cell, that is, the number
        of Avogadro times the atoms in a unit cell. To compare with experimental data the result
        should be divided by the number of unit formulas in the cell. If the structure is provided
        the division is performed internally and the result is in J/mol

        Args:
            t: a temperature in K
            structure: the structure of the system. If not None it will be used to determine the number of
                formula units
        Returns:
            Phonon contribution to the Helmholtz free energy
        """

        if t == 0:
            return self.zero_point_energy(structure=structure)

        freqs = self._positive_frequencies
        dens = self._positive_densities

        wd2kt = freqs / (2 * BOLTZ_THZ_PER_K * t)
        f = np.trapz(np.log(2 * np.sinh(wd2kt)) * dens, x=freqs)

        f *= const.Boltzmann * const.Avogadro * t

        if structure:
            formula_units = structure.composition.num_atoms / structure.composition.reduced_composition.num_atoms
            f /= formula_units

        return f

    def zero_point_energy(self, structure=None):
        """
        Zero point energy of the system. Only positive frequencies will be used.
        Result in J/mol-c. A mol-c is the abbreviation of a mole-cell, that is, the number
        of Avogadro times the atoms in a unit cell. To compare with experimental data the result
        should be divided by the number of unit formulas in the cell. If the structure is provided
        the division is performed internally and the result is in J/mol

        Args:
            structure: the structure of the system. If not None it will be used to determine the number of
                formula units
        Returns:
            Phonon contribution to the internal energy
        """

        freqs = self._positive_frequencies
        dens = self._positive_densities

        zpe = 0.5 * np.trapz(freqs * dens, x=freqs)
        zpe *= THZ_TO_J * const.Avogadro

        if structure:
            formula_units = structure.composition.num_atoms / structure.composition.reduced_composition.num_atoms
            zpe /= formula_units

        return zpe


class CompletePhononDos(PhononDos):
    """
    This wrapper class defines a total dos, and also provides a list of PDos.

    .. attribute:: pdos

        Dict of partial densities of the form {Site:Densities}
    """

    def __init__(self, structure: Structure, total_dos, pdoss):
        """
        Args:
            structure: Structure associated with this particular DOS.
            total_dos: total Dos for structure
            pdoss: The pdoss are supplied as an {Site: Densities}
        """
        super().__init__(frequencies=total_dos.frequencies, densities=total_dos.densities)
        self.pdos = {site: np.array(dens) for site, dens in pdoss.items()}
        self.structure = structure

    def get_site_dos(self, site):
        """
        Get the Dos for a site.

        Args:
            site: Site in Structure associated with CompletePhononDos.

        Returns:
            PhononDos containing summed orbital densities for site.
        """
        return PhononDos(self.frequencies, self.pdos[site])

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
                el_dos[el] = np.array(atom_dos)
            else:
                el_dos[el] += np.array(atom_dos)
        return {el: PhononDos(self.frequencies, densities) for el, densities in el_dos.items()}

    @classmethod
    def from_dict(cls, d):
        """
        Returns CompleteDos object from dict representation.
        """
        tdos = PhononDos.from_dict(d)
        struct = Structure.from_dict(d["structure"])
        pdoss = {}
        for at, pdos in zip(struct, d["pdos"]):
            pdoss[at] = pdos

        return cls(struct, tdos, pdoss)

    def as_dict(self):
        """
        JSON-serializable dict representation of CompletePhononDos.
        """
        d = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "structure": self.structure.as_dict(),
            "frequencies": list(self.frequencies),
            "densities": list(self.densities),
            "pdos": [],
        }
        if len(self.pdos) > 0:
            for at in self.structure:
                d["pdos"].append(list(self.pdos[at]))
        return d

    def __str__(self):
        return "Complete phonon DOS for " + str(self.structure)
