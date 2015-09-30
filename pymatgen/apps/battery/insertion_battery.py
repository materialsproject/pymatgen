# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module is used for analysis of materials with potential application as
intercalation batteries.
"""


__author__ = "Anubhav Jain, Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Jan 13, 2012"
__status__ = "Beta"

import itertools

from pymatgen.core.composition import Composition
from pymatgen.core.units import Charge, Time
from pymatgen.core.physical_constants import AVOGADROS_CONST
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.apps.battery.battery_abc import AbstractElectrode, \
    AbstractVoltagePair
from pymatgen.core.periodic_table import Element


class InsertionElectrode(AbstractElectrode):
    """
    A set of topotactically related compounds, with different amounts of a
    single element, e.g. TiO2 and LiTiO2, that can be used to define an
    insertion battery electrode.
    """

    def __init__(self, entries, working_ion_entry):
        """
        Create a new InsertionElectrode.

        Args:
            entries: A list of ComputedStructureEntries (or subclasses)
                representing the different topotactic states of the battery,
                e.g. TiO2 and LiTiO2.
            working_ion_entry: A single ComputedEntry or PDEntry
                representing the element that carries charge across the
                battery, e.g. Li.
        """
        self._entries = entries
        self._working_ion = working_ion_entry.composition.elements[0]
        self._working_ion_entry = working_ion_entry

        #Prepare to make phase diagram: determine elements and set their energy
        #to be very high
        elements = set()
        for entry in entries:
            elements.update(entry.composition.elements)

        #Set an artificial energy for each element for convex hull generation
        element_energy = max([entry.energy_per_atom for entry in entries]) + 10

        pdentries = []
        pdentries.extend(entries)
        pdentries.extend([PDEntry(Composition({el:1}), element_energy)
                          for el in elements])

        #Make phase diagram to determine which entries are stable vs. unstable
        pd = PhaseDiagram(pdentries)

        lifrac = lambda e: e.composition.get_atomic_fraction(self._working_ion)

        #stable entries ordered by amount of Li asc
        self._stable_entries = tuple(sorted([e for e in pd.stable_entries
                                             if e in entries], key=lifrac))

        #unstable entries ordered by amount of Li asc
        self._unstable_entries = tuple(sorted([e for e in pd.unstable_entries
                                               if e in entries], key=lifrac))

        #create voltage pairs
        self._vpairs = tuple([InsertionVoltagePair(self._stable_entries[i],
                                                   self._stable_entries[i + 1],
                                                   working_ion_entry)
                              for i in range(len(self._stable_entries) - 1)])

    @property
    def working_ion(self):
        """
        The working ion as an Element object
        """
        return self._working_ion

    @property
    def working_ion_entry(self):
        return self._working_ion_entry

    @property
    def voltage_pairs(self):
        return self._vpairs

    def get_stable_entries(self, charge_to_discharge=True):
        """
        Get the stable entries.

        Args:
            charge_to_discharge: order from most charge to most discharged
                state? Default to True.

        Returns:
            A list of stable entries in the electrode, ordered by amount of the
            working ion.
        """
        list_copy = list(self._stable_entries)
        return list_copy if charge_to_discharge else list_copy.reverse()

    def get_unstable_entries(self, charge_to_discharge=True):
        """
        Returns the unstable entries for the electrode.

        Args:
            charge_to_discharge: Order from most charge to most discharged
                state? Defaults to True.

        Returns:
            A list of unstable entries in the electrode, ordered by amount of
            the working ion.
        """
        list_copy = list(self._unstable_entries)
        return list_copy if charge_to_discharge else list_copy.reverse()

    def get_all_entries(self, charge_to_discharge=True):
        """
        Return all entries input for the electrode.

        Args:
            charge_to_discharge:
                order from most charge to most discharged state? Defaults to
                True.

        Returns:
            A list of all entries in the electrode (both stable and unstable),
            ordered by amount of the working ion.
        """
        all_entries = list(self.get_stable_entries())
        all_entries.extend(self.get_unstable_entries())
        #sort all entries by amount of working ion ASC
        fsrt = lambda e: e.composition.get_atomic_fraction(self.working_ion)
        all_entries = sorted([e for e in all_entries],
                             key=fsrt)
        return all_entries if charge_to_discharge else all_entries.reverse()

    @property
    def fully_charged_entry(self):
        """
        The most charged entry along the topotactic path.
        """
        return self._stable_entries[0]

    @property
    def fully_discharged_entry(self):
        """
        The most discharged entry along the topotactic path.
        """
        return self._stable_entries[-1]

    def get_max_instability(self, min_voltage=None, max_voltage=None):
        """
        The maximum instability along a path for a specific voltage range.

        Args:
            min_voltage: The minimum allowable voltage.
            max_voltage: The maximum allowable voltage.

        Returns:
            Maximum decomposition energy of all compounds along the insertion
            path (a subset of the path can be chosen by the optional arguments)
        """
        data = []
        for pair in self._select_in_voltage_range(min_voltage, max_voltage):
            if pair.decomp_e_charge is not None:
                data.append(pair.decomp_e_charge)
            if pair.decomp_e_discharge is not None:
                data.append(pair.decomp_e_discharge)
        return max(data) if len(data) > 0 else None

    def get_min_instability(self, min_voltage=None, max_voltage=None):
        """
        The minimum instability along a path for a specific voltage range.

        Args:
            min_voltage: The minimum allowable voltage.
            max_voltage: The maximum allowable voltage.

        Returns:
            Minimum decomposition energy of all compounds along the insertion
            path (a subset of the path can be chosen by the optional arguments)
        """
        data = []
        for pair in self._select_in_voltage_range(min_voltage, max_voltage):
            if pair.decomp_e_charge is not None:
                data.append(pair.decomp_e_charge)
            if pair.decomp_e_discharge is not None:
                data.append(pair.decomp_e_discharge)
        return min(data) if len(data) > 0 else None

    def get_max_muO2(self, min_voltage=None, max_voltage=None):
        """
        Maximum critical oxygen chemical potential along path.

        Args:
            min_voltage: The minimum allowable voltage.
            max_voltage: The maximum allowable voltage.

        Returns:
            Maximum critical oxygen chemical of all compounds along the
            insertion path (a subset of the path can be chosen by the optional
            arguments).
        """
        data = []
        for pair in self._select_in_voltage_range(min_voltage, max_voltage):
            if pair.muO2_discharge is not None:
                data.append(pair.pair.muO2_discharge)
            if pair.muO2_charge is not None:
                data.append(pair.muO2_charge)
        return max(data) if len(data) > 0 else None

    def get_min_muO2(self, min_voltage=None, max_voltage=None):
        """
        Minimum critical oxygen chemical potential along path.

        Args:
            min_voltage: The minimum allowable voltage for a given step
            max_voltage: The maximum allowable voltage allowable for a given
                step

        Returns:
            Minimum critical oxygen chemical of all compounds along the
            insertion path (a subset of the path can be chosen by the optional
            arguments).
        """
        data = []
        for pair in self._select_in_voltage_range(min_voltage, max_voltage):
            if pair.pair.muO2_discharge is not None:
                data.append(pair.pair.muO2_discharge)
            if pair.muO2_charge is not None:
                data.append(pair.muO2_charge)
        return min(data) if len(data) > 0 else None

    def get_sub_electrodes(self, adjacent_only=True, include_myself=True):
        """
        If this electrode contains multiple voltage steps, then it is possible
        to use only a subset of the voltage steps to define other electrodes.
        For example, an LiTiO2 electrode might contain three subelectrodes:
        [LiTiO2 --> TiO2, LiTiO2 --> Li0.5TiO2, Li0.5TiO2 --> TiO2]
        This method can be used to return all the subelectrodes with some
        options

        Args:
            adjacent_only: Only return electrodes from compounds that are
                adjacent on the convex hull, i.e. no electrodes returned
                will have multiple voltage steps if this is set True.
            include_myself: Include this identical electrode in the list of
                results.

        Returns:
            A list of InsertionElectrode objects
        """
        battery_list = []
        pair_it = self._vpairs if adjacent_only \
            else itertools.combinations_with_replacement(self._vpairs, 2)

        ion = self._working_ion

        for pair in pair_it:
            entry_charge = pair.entry_charge if adjacent_only \
                else pair[0].entry_charge
            entry_discharge = pair.entry_discharge if adjacent_only \
                else pair[1].entry_discharge

            chg_frac = entry_charge.composition.get_atomic_fraction(ion)
            dischg_frac = entry_discharge.composition.get_atomic_fraction(ion)

            def in_range(entry):
                frac = entry.composition.get_atomic_fraction(ion)
                return chg_frac <= frac <= dischg_frac

            if include_myself or entry_charge != self.fully_charged_entry \
                    or entry_discharge != self.fully_discharged_entry:
                unstable_entries = filter(in_range,
                                          self.get_unstable_entries())
                stable_entries = filter(in_range, self.get_stable_entries())
                all_entries = list(stable_entries)
                all_entries.extend(unstable_entries)
                battery_list.append(self.__class__(all_entries,
                                                   self.working_ion_entry))
        return battery_list

    def as_dict_summary(self, print_subelectrodes=True):
        """
        Generate a summary dict.

        Args:
            print_subelectrodes: Also print data on all the possible
                subelectrodes.

        Returns:
            A summary of this electrode"s properties in dict format.
        """
        chg_comp = self.fully_charged_entry.composition
        dischg_comp = self.fully_discharged_entry.composition
        ion = self.working_ion
        d = {"average_voltage": self.get_average_voltage(),
             "max_voltage": self.max_voltage,
             "min_voltage": self.min_voltage,
             "max_delta_volume": self.max_delta_volume,
             "max_voltage_step": self.max_voltage_step,
             "capacity_grav": self.get_capacity_grav(),
             "capacity_vol": self.get_capacity_vol(),
             "energy_grav": self.get_specific_energy(),
             "energy_vol": self.get_energy_density(),
             "working_ion": self._working_ion.symbol,
             "nsteps": self.num_steps,
             "framework": self._vpairs[0].framework.to_data_dict,
             "formula_charge": chg_comp.reduced_formula,
             "formula_discharge": dischg_comp.reduced_formula,
             "fracA_charge": chg_comp.get_atomic_fraction(ion),
             "fracA_discharge": dischg_comp.get_atomic_fraction(ion),
             "max_instability": self.get_max_instability(),
             "min_instability": self.get_min_instability()}
        if print_subelectrodes:
            f_dict = lambda c: c.as_dict_summary(print_subelectrodes=False)
            d["adj_pairs"] = map(f_dict,
                                 self.get_sub_electrodes(adjacent_only=True))
            d["all_pairs"] = map(f_dict,
                                 self.get_sub_electrodes(adjacent_only=False))
        return d

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        output = []
        chg_form = self.fully_charged_entry.composition.reduced_formula
        dischg_form = self.fully_discharged_entry.composition.reduced_formula
        output.append("InsertionElectrode with endpoints at {} and {}".format(
                      chg_form, dischg_form))
        output.append("Avg. volt. = {} V".format(self.get_average_voltage()))
        output.append("Grav. cap. = {} mAh/g".format(self.get_capacity_grav()))
        output.append("Vol. cap. = {}".format(self.get_capacity_vol()))
        return  "\n".join(output)

    @classmethod
    def from_dict(cls, d):
        from monty.json import MontyDecoder
        dec = MontyDecoder()
        return cls(dec.process_decoded(d["entries"]),
                   dec.process_decoded(d["working_ion_entry"]))

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "entries": [entry.as_dict() for entry in self._entries],
                "working_ion_entry": self.working_ion_entry.as_dict()}


class InsertionVoltagePair(AbstractVoltagePair):
    """
    Defines an Insertion Voltage Pair.

    Args:
        entry1: Entry corresponding to one of the entries in the voltage step.
        entry2: Entry corresponding to the other entry in the voltage step.
        working_ion_entry: A single ComputedEntry or PDEntry representing
            the element that carries charge across the battery, e.g. Li.
    """

    def __init__(self, entry1, entry2, working_ion_entry):
        #initialize some internal variables
        working_element = working_ion_entry.composition.elements[0]

        entry_charge = entry1
        entry_discharge = entry2
        if entry_charge.composition.get_atomic_fraction(working_element) \
                > entry2.composition.get_atomic_fraction(working_element):
            (entry_charge, entry_discharge) = (entry_discharge, entry_charge)

        comp_charge = entry_charge.composition
        comp_discharge = entry_discharge.composition

        ion_sym = working_element.symbol

        frame_charge_comp = Composition({el: comp_charge[el]
                                         for el in comp_charge
                                         if el.symbol != ion_sym})
        frame_discharge_comp = Composition({el: comp_discharge[el]
                                            for el in comp_discharge
                                            if el.symbol != ion_sym})

        #Data validation

        #check that the ion is just a single element
        if not working_ion_entry.composition.is_element:
            raise ValueError("VoltagePair: The working ion specified must be "
                             "an element")

        #check that at least one of the entries contains the working element
        if not comp_charge.get_atomic_fraction(working_element) > 0 and \
                not comp_discharge.get_atomic_fraction(working_element) > 0:
            raise ValueError("VoltagePair: The working ion must be present in "
                             "one of the entries")

        #check that the entries do not contain the same amount of the workin
        #element
        if comp_charge.get_atomic_fraction(working_element) == \
                comp_discharge.get_atomic_fraction(working_element):
            raise ValueError("VoltagePair: The working ion atomic percentage "
                             "cannot be the same in both the entries")

        #check that the frameworks of the entries are equivalent
        if not frame_charge_comp.reduced_formula == \
                frame_discharge_comp.reduced_formula:
            raise ValueError("VoltagePair: the specified entries must have the"
                             " same compositional framework")

        #Initialize normalization factors, charged and discharged entries

        valence_list = Element(ion_sym).oxidation_states
        working_ion_valence = max(valence_list)


        (self.framework,
         norm_charge) = frame_charge_comp.get_reduced_composition_and_factor()
        norm_discharge = \
            frame_discharge_comp.get_reduced_composition_and_factor()[1]

        self._working_ion_entry = working_ion_entry

        #Initialize normalized properties
        self._vol_charge = entry_charge.structure.volume / norm_charge
        self._vol_discharge = entry_discharge.structure.volume / norm_discharge

        comp_charge = entry_charge.composition
        comp_discharge = entry_discharge.composition

        self._mass_charge = comp_charge.weight / norm_charge
        self._mass_discharge = comp_discharge.weight / norm_discharge

        self._num_ions_transferred = \
            (comp_discharge[working_element] / norm_discharge) \
            - (comp_charge[working_element] / norm_charge)

        self._voltage = \
            (((entry_charge.energy / norm_charge) -
             (entry_discharge.energy / norm_discharge)) / \
            self._num_ions_transferred + working_ion_entry.energy_per_atom) / working_ion_valence
        self._mAh = self._num_ions_transferred * Charge(1, "e").to("C") * \
            Time(1, "s").to("h") * AVOGADROS_CONST * 1000 * working_ion_valence

        #Step 4: add (optional) hull and muO2 data
        self.decomp_e_charge = \
            entry_charge.data.get("decomposition_energy", None)
        self.decomp_e_discharge = \
            entry_discharge.data.get("decomposition_energy", None)

        self.muO2_charge = entry_charge.data.get("muO2", None)
        self.muO2_discharge = entry_discharge.data.get("muO2", None)

        self.entry_charge = entry_charge
        self.entry_discharge = entry_discharge
        self.normalization_charge = norm_charge
        self.normalization_discharge = norm_discharge
        self._frac_charge = comp_charge.get_atomic_fraction(working_element)
        self._frac_discharge = \
            comp_discharge.get_atomic_fraction(working_element)

    @property
    def frac_charge(self):
        return self._frac_charge

    @property
    def frac_discharge(self):
        return self._frac_discharge

    @property
    def voltage(self):
        return self._voltage

    @property
    def mAh(self):
        return self._mAh

    @property
    def mass_charge(self):
        return self._mass_charge

    @property
    def mass_discharge(self):
        return self._mass_discharge

    @property
    def vol_charge(self):
        return self._vol_charge

    @property
    def vol_discharge(self):
        return self._vol_discharge

    @property
    def working_ion_entry(self):
        return self._working_ion_entry

    def __repr__(self):
        output = ["Insertion voltage pair with working ion {}"
                  .format(self._working_ion_entry.composition.reduced_formula),
                  "V = {}, mAh = {}".format(self.voltage, self.mAh),
                  "mass_charge = {}, mass_discharge = {}"
                  .format(self.mass_charge, self.mass_discharge),
                  "vol_charge = {}, vol_discharge = {}"
                  .format(self.vol_charge, self.vol_discharge),
                  "frac_charge = {}, frac_discharge = {}"
                  .format(self.frac_charge, self.frac_discharge)]
        return "\n".join(output)

    def __str__(self):
        return self.__repr__()
