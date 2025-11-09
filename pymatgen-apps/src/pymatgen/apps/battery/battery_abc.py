"""This module defines the abstract base classes for battery-related classes.
Regardless of the kind of electrode, conversion or insertion, there are many
common definitions and properties, e.g. average voltage, capacity, etc. which
can be defined in a general way. The Abc for battery classes implements some of
these common definitions to allow sharing of common logic between them.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

from monty.json import MSONable
from scipy.constants import N_A

from pymatgen.core import Composition

if TYPE_CHECKING:
    from pymatgen.entries import Entry
    from pymatgen.util.typing import SpeciesLike

__author__ = "Anubhav Jain, Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Feb 1, 2012"
__status__ = "Beta"


@dataclass
class AbstractVoltagePair(MSONable):
    """An Abstract Base Class for a Voltage Pair.

    Attributes:
        voltage : Voltage of voltage pair.
        mAh: Energy in mAh.
        mass_charge: Mass of charged pair.
        mass_discharge: Mass of discharged pair.
        vol_charge: Vol of charged pair.
        vol_discharge: Vol of discharged pair.
        frac_charge: Frac of working ion in charged pair.
        frac_discharge: Frac of working ion in discharged pair.
        working_ion_entry: Working ion as an entry.
        framework_formula : The compositions of one formula unit of the host material
    """

    voltage: float
    mAh: float
    mass_charge: float
    mass_discharge: float
    vol_charge: float
    vol_discharge: float
    frac_charge: float
    frac_discharge: float
    working_ion_entry: Entry
    framework_formula: str

    def __post_init__(self):
        # ensure the frame work is a reduced composition
        fw = Composition(self.framework_formula)
        self.framework_formula = fw.reduced_formula

    @property
    def working_ion(self) -> SpeciesLike:
        """Working ion as pymatgen Element object."""
        return self.working_ion_entry.elements[0]

    @property
    def framework(self) -> Composition:
        """The composition object representing the framework."""
        return Composition(self.framework_formula)

    @property
    def x_charge(self) -> float:
        """The number of working ions per formula unit of host in the charged state."""
        return self.frac_charge * self.framework.num_atoms / (1 - self.frac_charge)

    @property
    def x_discharge(self) -> float:
        """The number of working ions per formula unit of host in the discharged state."""
        return self.frac_discharge * self.framework.num_atoms / (1 - self.frac_discharge)


@dataclass
class AbstractElectrode(Sequence, MSONable):
    """An Abstract Base Class representing an Electrode. It is essentially a
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

    Attributes:
        voltage_pairs: Objects that represent each voltage step
        working_ion: Representation of the working ion that only contains element type
        working_ion_entry: Representation of the working_ion that contains the energy
        framework_formula: The compositions of one formula unit of the host material
    """

    voltage_pairs: tuple[AbstractVoltagePair, ...]
    working_ion_entry: Entry
    framework_formula: str  # should be made into Composition whenever the as_dict and from dict are fixed

    def __post_init__(self):
        # ensure the frame work is a reduced composition
        self.framework_formula = self.framework.reduced_formula

    def __getitem__(self, index):
        return self.voltage_pairs[index]

    def __contains__(self, obj):
        return obj in self.voltage_pairs

    def __iter__(self):
        return iter(self.voltage_pairs)

    def __len__(self):
        return len(self.voltage_pairs)

    @property
    def working_ion(self):
        """Working ion as pymatgen Element object."""
        return self.working_ion_entry.elements[0]

    @property
    def framework(self):
        """The composition object representing the framework."""
        return Composition(self.framework_formula)

    @property
    def x_charge(self) -> float:
        """The number of working ions per formula unit of host in the charged state."""
        return self.voltage_pairs[0].x_charge

    @property
    def x_discharge(self) -> float:
        """The number of working ions per formula unit of host in the discharged state."""
        return self.voltage_pairs[-1].x_discharge

    @property
    def max_delta_volume(self):
        """Maximum volume change along insertion."""
        vols = [v.vol_charge for v in self.voltage_pairs]
        vols.extend([v.vol_discharge for v in self.voltage_pairs])
        return max(vols) / min(vols) - 1

    @property
    def num_steps(self):
        """The number of distinct voltage steps in from fully charge to discharge
        based on the stable intermediate states.
        """
        return len(self.voltage_pairs)

    @property
    def max_voltage(self):
        """Highest voltage along insertion."""
        return max(p.voltage for p in self.voltage_pairs)

    @property
    def min_voltage(self):
        """Lowest voltage along insertion."""
        return min(p.voltage for p in self.voltage_pairs)

    @property
    def max_voltage_step(self):
        """Maximum absolute difference in adjacent voltage steps."""
        steps = [
            self.voltage_pairs[i].voltage - self.voltage_pairs[i + 1].voltage
            for i in range(len(self.voltage_pairs) - 1)
        ]
        return max(steps) if len(steps) > 0 else 0

    @property
    def normalization_mass(self):
        """The mass used for normalization. This is the mass of the discharged
        electrode of the last voltage pair.
        """
        return self.voltage_pairs[-1].mass_discharge

    @property
    def normalization_volume(self):
        """The mass used for normalization. This is the vol of the discharged
        electrode of the last voltage pair.
        """
        return self.voltage_pairs[-1].vol_discharge

    def get_sub_electrodes(self, adjacent_only=True):
        """If this electrode contains multiple voltage steps, then it is possible
        to use only a subset of the voltage steps to define other electrodes.
        Must be implemented for each electrode object.

        Args:
            adjacent_only: Only return electrodes from compounds that are
                adjacent on the convex hull, i.e. no electrodes returned
                will have multiple voltage steps if this is set true

        Returns:
            A list of Electrode objects
        """
        raise NotImplementedError(
            "The get_sub_electrodes function must be implemented for each concrete electrode "
            f"class {type(self).__name__}"
        )

    def get_average_voltage(self, min_voltage=None, max_voltage=None):
        """Average voltage for path satisfying between a min and max voltage.

        Args:
            min_voltage (float): The minimum allowable voltage for a given
                step.
            max_voltage (float): The maximum allowable voltage allowable for a
                given step.

        Returns:
            Average voltage in V across the insertion path (a subset of the
            path can be chosen by the optional arguments)
        """
        pairs_in_range = self._select_in_voltage_range(min_voltage, max_voltage)
        if len(pairs_in_range) == 0:
            return 0
        total_cap_in_range = sum(p.mAh for p in pairs_in_range)
        total_edens_in_range = sum(p.mAh * p.voltage for p in pairs_in_range)
        return total_edens_in_range / total_cap_in_range

    def get_capacity_grav(self, min_voltage=None, max_voltage=None, use_overall_normalization=True):
        """Get the gravimetric capacity of the electrode.

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
        pairs_in_range = self._select_in_voltage_range(min_voltage, max_voltage)
        normalization_mass = (
            self.normalization_mass
            if use_overall_normalization or len(pairs_in_range) == 0
            else pairs_in_range[-1].mass_discharge
        )
        return sum(pair.mAh for pair in pairs_in_range) / normalization_mass

    def get_capacity_vol(self, min_voltage=None, max_voltage=None, use_overall_normalization=True):
        """Get the volumetric capacity of the electrode.

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
        pairs_in_range = self._select_in_voltage_range(min_voltage, max_voltage)
        normalization_vol = (
            self.normalization_volume
            if use_overall_normalization or len(pairs_in_range) == 0
            else pairs_in_range[-1].vol_discharge
        )
        return sum(pair.mAh for pair in pairs_in_range) / normalization_vol * 1e24 / N_A

    def get_specific_energy(self, min_voltage=None, max_voltage=None, use_overall_normalization=True):
        """Get the specific energy of the battery in mAh/g.

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
        return self.get_capacity_grav(min_voltage, max_voltage, use_overall_normalization) * self.get_average_voltage(
            min_voltage, max_voltage
        )

    def get_energy_density(self, min_voltage=None, max_voltage=None, use_overall_normalization=True):
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
        return self.get_capacity_vol(min_voltage, max_voltage, use_overall_normalization) * self.get_average_voltage(
            min_voltage, max_voltage
        )

    def _select_in_voltage_range(self, min_voltage=None, max_voltage=None):
        """Selects VoltagePairs within a certain voltage range.

        Args:
            min_voltage (float): The minimum allowable voltage for a given
                step.
            max_voltage (float): The maximum allowable voltage allowable for a
                given step.

        Returns:
            A list of VoltagePair objects
        """
        min_voltage = min_voltage if min_voltage is not None else self.min_voltage
        max_voltage = max_voltage if max_voltage is not None else self.max_voltage
        return list(filter(lambda p: min_voltage <= p.voltage <= max_voltage, self.voltage_pairs))

    def get_summary_dict(self, print_subelectrodes=True) -> dict:
        """Generate a summary dict.

        Args:
            print_subelectrodes: Also print data on all the possible
                subelectrodes.

        Returns:
            A summary of this electrode's properties in dict format.
        """
        dct = {
            "average_voltage": self.get_average_voltage(),
            "max_voltage": self.max_voltage,
            "min_voltage": self.min_voltage,
            "max_delta_volume": self.max_delta_volume,
            "max_voltage_step": self.max_voltage_step,
            "capacity_grav": self.get_capacity_grav(),
            "capacity_vol": self.get_capacity_vol(),
            "energy_grav": self.get_specific_energy(),
            "energy_vol": self.get_energy_density(),
            "working_ion": self.working_ion.symbol,
            "nsteps": self.num_steps,
            "fracA_charge": self.voltage_pairs[0].frac_charge,
            "fracA_discharge": self.voltage_pairs[-1].frac_discharge,
            "framework_formula": self.framework_formula,
        }

        if print_subelectrodes:

            def f_dict(c):
                return c.get_summary_dict(print_subelectrodes=False)

            dct["adj_pairs"] = list(map(f_dict, self.get_sub_electrodes(adjacent_only=True)))
            dct["all_pairs"] = list(map(f_dict, self.get_sub_electrodes(adjacent_only=False)))
        return dct
