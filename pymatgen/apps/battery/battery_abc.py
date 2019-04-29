# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from collections.abc import Sequence
import abc

from monty.json import MSONable

from scipy.constants import N_A

"""
This module defines the abstract base classes for battery-related classes.
Regardless of the kind of electrode, conversion or insertion, there are many
common definitions and properties, e.g., average voltage, capacity, etc. which
can be defined in a general way. The Abc for battery classes implements some of
these common definitions to allow sharing of common logic between them.
"""


__author__ = "Anubhav Jain, Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Feb 1, 2012"
__status__ = "Beta"


class AbstractVoltagePair:
    """
    An Abstract Base Class for a Voltage Pair.
    """
    __metaclass__ = abc.ABCMeta

    @property
    @abc.abstractmethod
    def voltage(self):
        return self._voltage

    @property
    @abc.abstractmethod
    def mAh(self):
        return self._mAh

    @property
    @abc.abstractmethod
    def mass_charge(self):
        return self._mass_charge

    @property
    @abc.abstractmethod
    def mass_discharge(self):
        return self._mass_discharge

    @property
    @abc.abstractmethod
    def vol_charge(self):
        return self._vol_charge

    @property
    @abc.abstractmethod
    def vol_discharge(self):
        return self._vol_discharge

    @property
    @abc.abstractmethod
    def frac_charge(self):
        return self._frac_charge

    @property
    @abc.abstractmethod
    def frac_discharge(self):
        return self._frac_discharge

    @property
    @abc.abstractmethod
    def working_ion_entry(self):
        return self._working_ion_entry


class AbstractElectrode(Sequence, MSONable):
    """
    An Abstract Base Class representing an Electrode. It is essentially a
    sequence of VoltagePairs. Generally, subclasses only need to implement
    three abstract properties: voltage_pairs, working_ion and
    working_ion_entry.

    The general concept is that all other battery properties such as capacity,
    etc. are derived from voltage pairs.

    One of the major challenges with representing battery materials is keeping
    track of the normalization between different entries. For example, one
    entry might be TiO2 with one unit cell whereas another is LiTi2O4 with two
    unit cells. When computing battery properties, it is needed to always use
    a universal reference state otherwise you have normalization errors (e.g.,
    the energy of LiTi2O4 must be divided by two to be compared with TiO2).

    For properties such as volume, mass, or mAh transferred within the voltage
    pair, a universal convention is necessary. AbstractElectrode can query for
    extrinsic properties of several different AbstractVoltagePairs belonging to
    a single charge/discharge path and be confident that the normalization is
    being carried out properly throughout, even if more AbstractVoltagePairs
    are added later.

    The universal normalization is defined by the reduced structural framework
    of the entries, which is common along the entire charge/discharge path. For
    example, LiTi2O4 has a reduced structural framework of TiO2. Another
    example is Li9V6P16O58 which would have a reduced structural framework of
    V3P8O29. Note that reduced structural frameworks need not be
    charge-balanced or physical, e.g. V3P8O29 is not charge-balanced, they are
    just a tool for normalization.

    Example: for a LiTi2O4 -> TiO2 AbstractVoltagePair, extrinsic quantities
    like mAh or cell volumes are given per TiO2 formula unit.

    Developers implementing a new battery (other than the two general ones
    already implemented) need to implement a VoltagePair and an Electrode.
    """

    __metaclass__ = abc.ABCMeta

    @property
    @abc.abstractmethod
    def voltage_pairs(self):
        """
        Returns all the VoltagePairs
        """
        return

    @property
    @abc.abstractmethod
    def working_ion(self):
        """
        The working ion as an Element object
        """
        return

    @property
    @abc.abstractmethod
    def working_ion_entry(self):
        """
        The working ion as an Entry object
        """
        return

    def __getitem__(self, index):
        return self.voltage_pairs[index]

    def __contains__(self, obj):
        return obj in self.voltage_pairs

    def __iter__(self):
        return self.voltage_pairs.__iter__()

    def __len__(self):
        return len(self.voltage_pairs)

    @property
    def max_delta_volume(self):
        """
        Maximum volume change along insertion
        """
        vols = [v.vol_charge for v in self.voltage_pairs]
        vols.extend([v.vol_discharge for v in self.voltage_pairs])
        return max(vols) / min(vols) - 1

    @property
    def num_steps(self):
        """
        The number of distinct voltage steps in from fully charge to discharge
        based on the stable intermediate states
        """
        return len(self.voltage_pairs)

    @property
    def max_voltage(self):
        """
        Highest voltage along insertion
        """
        return max([p.voltage for p in self.voltage_pairs])

    @property
    def min_voltage(self):
        """
        Lowest voltage along insertion
        """
        return min([p.voltage for p in self.voltage_pairs])

    @property
    def max_voltage_step(self):
        """
        Maximum absolute difference in adjacent voltage steps
        """
        steps = [self.voltage_pairs[i].voltage
                 - self.voltage_pairs[i + 1].voltage
                 for i in range(len(self.voltage_pairs) - 1)]
        return max(steps) if len(steps) > 0 else 0

    @property
    def normalization_mass(self):
        return self.voltage_pairs[-1].mass_discharge

    @property
    def normalization_volume(self):
        return self.voltage_pairs[-1].vol_discharge

    def get_average_voltage(self, min_voltage=None, max_voltage=None):
        """
        Average voltage for path satisfying between a min and max voltage.

        Args:
            min_voltage (float): The minimum allowable voltage for a given
                step.
            max_voltage (float): The maximum allowable voltage allowable for a
                given step.

        Returns:
            Average voltage in V across the insertion path (a subset of the
            path can be chosen by the optional arguments)
        """
        pairs_in_range = self._select_in_voltage_range(min_voltage,
                                                       max_voltage)
        if len(pairs_in_range) == 0:
            return 0
        total_cap_in_range = sum([p.mAh for p in pairs_in_range])
        total_edens_in_range = sum([p.mAh * p.voltage for p in pairs_in_range])
        return total_edens_in_range / total_cap_in_range

    def get_capacity_grav(self, min_voltage=None, max_voltage=None,
                          use_overall_normalization=True):
        """
        Get the gravimetric capacity of the electrode.

        Args:
            min_voltage (float): The minimum allowable voltage for a given
                step.
            max_voltage (float): The maximum allowable voltage allowable for a
                given step.
            use_overall_normalization (booL): If False, normalize by the
                discharged state of only the voltage pairs matching the voltage
                criteria. if True, use default normalization of the full
                electrode path.

        Returns:
            Gravimetric capacity in mAh/g across the insertion path (a subset
            of the path can be chosen by the optional arguments).
        """
        pairs_in_range = self._select_in_voltage_range(min_voltage,
                                                       max_voltage)
        normalization_mass = self.normalization_mass \
            if use_overall_normalization or len(pairs_in_range) == 0 \
            else pairs_in_range[-1].mass_discharge
        return sum([pair.mAh for pair in pairs_in_range]) / normalization_mass

    def get_capacity_vol(self, min_voltage=None, max_voltage=None,
                         use_overall_normalization=True):
        """
        Get the volumetric capacity of the electrode.

        Args:
            min_voltage (float): The minimum allowable voltage for a given
                step.
            max_voltage (float): The maximum allowable voltage allowable for a
                given step.
            use_overall_normalization (booL): If False, normalize by the
                discharged state of only the voltage pairs matching the voltage
                criteria. if True, use default normalization of the full
                electrode path.

        Returns:
            Volumetric capacity in mAh/cc across the insertion path (a subset
            of the path can be chosen by the optional arguments)
        """
        pairs_in_range = self._select_in_voltage_range(min_voltage,
                                                       max_voltage)
        normalization_vol = self.normalization_volume \
            if use_overall_normalization or len(pairs_in_range) == 0 \
            else pairs_in_range[-1].vol_discharge
        return sum([pair.mAh for pair in pairs_in_range]) / normalization_vol \
            * 1e24 / N_A

    def get_specific_energy(self, min_voltage=None, max_voltage=None,
                            use_overall_normalization=True):
        """
        Returns the specific energy of the battery in mAh/g.

        Args:
            min_voltage (float): The minimum allowable voltage for a given
                step.
            max_voltage (float): The maximum allowable voltage allowable for a
                given step.
            use_overall_normalization (booL): If False, normalize by the
                discharged state of only the voltage pairs matching the voltage
                criteria. if True, use default normalization of the full
                electrode path.

        Returns:
            Specific energy in Wh/kg across the insertion path (a subset of
            the path can be chosen by the optional arguments)
        """
        return self.get_capacity_grav(min_voltage, max_voltage,
                                      use_overall_normalization) \
            * self.get_average_voltage(min_voltage, max_voltage)

    def get_energy_density(self, min_voltage=None, max_voltage=None,
                           use_overall_normalization=True):
        """
        Args:
            min_voltage (float): The minimum allowable voltage for a given
                step.
            max_voltage (float): The maximum allowable voltage allowable for a
                given step.
            use_overall_normalization (booL): If False, normalize by the
                discharged state of only the voltage pairs matching the voltage
                criteria. if True, use default normalization of the full
                electrode path.

        Returns:
            Energy density in Wh/L across the insertion path (a subset of the
            path can be chosen by the optional arguments).
        """
        return self.get_capacity_vol(min_voltage, max_voltage,
                                     use_overall_normalization) \
            * self.get_average_voltage(min_voltage, max_voltage)

    def _select_in_voltage_range(self, min_voltage=None, max_voltage=None):
        """
        Selects VoltagePairs within a certain voltage range.

        Args:
            min_voltage (float): The minimum allowable voltage for a given
                step.
            max_voltage (float): The maximum allowable voltage allowable for a
                given step.

        Returns:
            A list of VoltagePair objects
        """
        min_voltage = min_voltage if min_voltage is not None \
            else self.min_voltage
        max_voltage = max_voltage if max_voltage is not None \
            else self.max_voltage
        return list(filter(lambda p: min_voltage <= p.voltage <= max_voltage,
                           self.voltage_pairs))
