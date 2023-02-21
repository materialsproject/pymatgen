# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module is used for analysis of materials with potential application as
intercalation batteries.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Iterable

from scipy.constants import N_A

from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.apps.battery.battery_abc import AbstractElectrode, AbstractVoltagePair
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.core.units import Charge, Time
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry

__author__ = "Anubhav Jain, Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"


@dataclass
class InsertionElectrode(AbstractElectrode):
    """
    A set of topotactically related compounds, with different amounts of a
    single element, e.g. TiO2 and LiTiO2, that can be used to define an
    insertion battery electrode.
    """

    stable_entries: Iterable[ComputedEntry]
    unstable_entries: Iterable[ComputedEntry]

    @classmethod
    def from_entries(
        cls,
        entries: Iterable[ComputedEntry | ComputedStructureEntry],
        working_ion_entry: ComputedEntry | ComputedStructureEntry | PDEntry,
        strip_structures: bool = False,
    ):
        """
        Create a new InsertionElectrode.

        Args:
            entries: A list of ComputedEntries, ComputedStructureEntries, or
                subclasses representing the different topotactic states
                of the battery, e.g. TiO2 and LiTiO2.
            working_ion_entry: A single ComputedEntry or PDEntry
                representing the element that carries charge across the
                battery, e.g. Li.
            strip_structures: Since the electrode document only uses volume we can make the
                electrode object significantly leaner by dropping the structure data.
                If this parameter is set to True, the ComputedStructureEntry will be
                replaced with a ComputedEntry and the volume will be stored in
                ComputedEntry.data['volume']. If entries provided are ComputedEntries,
                must set strip_structures=False.
        """
        if strip_structures:
            ents = []
            for ient in entries:
                dd = ient.as_dict()
                ent = ComputedEntry.from_dict(dd)
                ent.data["volume"] = ient.structure.volume
                ents.append(ent)
            entries = ents

        _working_ion = working_ion_entry.composition.elements[0]
        _working_ion_entry = working_ion_entry

        # Prepare to make phase diagram: determine elements and set their energy
        # to be very high
        elements = set()
        for entry in entries:
            elements.update(entry.composition.elements)

        # Set an artificial high energy for each element for convex hull generation
        element_energy = max(entry.energy_per_atom for entry in entries) + 10

        pdentries: list[ComputedEntry | ComputedStructureEntry | PDEntry] = []
        pdentries.extend(entries)
        pdentries.extend([PDEntry(Composition({el: 1}), element_energy) for el in elements])

        # Make phase diagram to determine which entries are stable vs. unstable.
        # For each working ion concentration, we want one stable entry
        # to use in forming voltage pairs. PhaseDiagram allows for easy comparison
        # of entry energies.
        pd = PhaseDiagram(pdentries)

        def lifrac(e):
            return e.composition.get_atomic_fraction(_working_ion)

        # stable entries ordered by amount of Li asc
        _stable_entries = tuple(sorted((e for e in pd.stable_entries if e in entries), key=lifrac))

        # unstable entries ordered by amount of Li asc
        _unstable_entries = tuple(sorted((e for e in pd.unstable_entries if e in entries), key=lifrac))

        # create voltage pairs
        _vpairs: tuple[AbstractVoltagePair, ...] = tuple(
            InsertionVoltagePair.from_entries(
                _stable_entries[i],
                _stable_entries[i + 1],
                working_ion_entry,
            )
            for i in range(len(_stable_entries) - 1)
        )
        framework = _vpairs[0].framework
        return cls(  # pylint: disable=E1123
            voltage_pairs=_vpairs,
            working_ion_entry=_working_ion_entry,
            stable_entries=_stable_entries,
            unstable_entries=_unstable_entries,
            framework_formula=framework.reduced_formula,
        )

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
        list_copy = list(self.stable_entries)
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
        list_copy = list(self.unstable_entries)
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
        # sort all entries by amount of working ion ASC
        all_entries = sorted(
            all_entries,
            key=lambda e: e.composition.get_atomic_fraction(self.working_ion),
        )
        return all_entries if charge_to_discharge else all_entries.reverse()

    @property
    def fully_charged_entry(self):
        """
        The most charged entry along the topotactic path.
        """
        return self.stable_entries[0]

    @property
    def fully_discharged_entry(self):
        """
        The most discharged entry along the topotactic path.
        """
        return self.stable_entries[-1]

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
            if getattr(pair, "decomp_e_charge", None) is not None:
                data.append(pair.decomp_e_charge)
            if getattr(pair, "decomp_e_discharge", None) is not None:
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
            if getattr(pair, "decomp_e_charge", None) is not None:
                data.append(pair.decomp_e_charge)
            if getattr(pair, "decomp_e_discharge", None) is not None:
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
                data.extend([d["chempot"] for d in pair.muO2_discharge])
            if pair.muO2_charge is not None:
                data.extend([d["chempot"] for d in pair.muO2_discharge])
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
            if pair.muO2_discharge is not None:
                data.extend([d["chempot"] for d in pair.muO2_discharge])
            if pair.muO2_charge is not None:
                data.extend([d["chempot"] for d in pair.muO2_discharge])
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
        pair_it = (
            self.voltage_pairs if adjacent_only else itertools.combinations_with_replacement(self.voltage_pairs, 2)
        )

        ion = self.working_ion

        for pair in pair_it:
            entry_charge = pair.entry_charge if adjacent_only else pair[0].entry_charge
            entry_discharge = pair.entry_discharge if adjacent_only else pair[1].entry_discharge

            def in_range(entry):
                chg_frac = entry_charge.composition.get_atomic_fraction(ion)
                dischg_frac = entry_discharge.composition.get_atomic_fraction(ion)
                frac = entry.composition.get_atomic_fraction(ion)
                return chg_frac <= frac <= dischg_frac

            if (
                include_myself
                or entry_charge != self.fully_charged_entry
                or entry_discharge != self.fully_discharged_entry
            ):
                unstable_entries = filter(in_range, self.get_unstable_entries())
                stable_entries = filter(in_range, self.get_stable_entries())
                all_entries = list(stable_entries)
                all_entries.extend(unstable_entries)
                battery_list.append(type(self).from_entries(all_entries, self.working_ion_entry))
        return battery_list

    def get_summary_dict(self, print_subelectrodes=True) -> dict:
        """
        Generate a summary dict.
        Populates the summary dict with the basic information from the parent method then populates more information.
        Since the parent method calls self.get_summary_dict(print_subelectrodes=True) for the subelectrodes.
        The current method will be called from within super().get_summary_dict.

        Args:
            print_subelectrodes: Also print data on all the possible
                subelectrodes.

        Returns:
            A summary of this electrode's properties in dict format.
        """
        d = super().get_summary_dict(print_subelectrodes=print_subelectrodes)

        chg_comp = self.fully_charged_entry.composition
        dischg_comp = self.fully_discharged_entry.composition

        d.update(
            {
                "id_charge": self.fully_charged_entry.entry_id,
                "formula_charge": chg_comp.reduced_formula,
                "id_discharge": self.fully_discharged_entry.entry_id,
                "formula_discharge": dischg_comp.reduced_formula,
                "max_instability": self.get_max_instability(),
                "min_instability": self.get_min_instability(),
                "material_ids": [itr_ent.entry_id for itr_ent in self.get_all_entries()],
                "stable_material_ids": [itr_ent.entry_id for itr_ent in self.get_stable_entries()],
                "unstable_material_ids": [itr_ent.entry_id for itr_ent in self.get_unstable_entries()],
            }
        )
        if all("decomposition_energy" in itr_ent.data for itr_ent in self.get_all_entries()):
            d.update(
                {
                    "stability_charge": self.fully_charged_entry.data["decomposition_energy"],
                    "stability_discharge": self.fully_discharged_entry.data["decomposition_energy"],
                    "stability_data": {
                        itr_ent.entry_id: itr_ent.data["decomposition_energy"] for itr_ent in self.get_all_entries()
                    },
                }
            )

        if all("muO2" in itr_ent.data for itr_ent in self.get_all_entries()):
            d.update({"muO2_data": {itr_ent.entry_id: itr_ent.data["muO2"] for itr_ent in self.get_all_entries()}})

        return d

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        output = []
        chg_form = self.fully_charged_entry.composition.reduced_formula
        dischg_form = self.fully_discharged_entry.composition.reduced_formula
        output.append(f"InsertionElectrode with endpoints at {chg_form} and {dischg_form}")
        output.append(f"Avg. volt. = {self.get_average_voltage()} V")
        output.append(f"Grav. cap. = {self.get_capacity_grav()} mAh/g")
        output.append(f"Vol. cap. = {self.get_capacity_vol()}")
        return "\n".join(output)

    @classmethod
    def from_dict_legacy(cls, d):
        """
        Args:
            d (dict): Dict representation

        Returns:
            InsertionElectrode
        """
        from monty.json import MontyDecoder

        dec = MontyDecoder()
        return InsertionElectrode(  # pylint: disable=E1120
            dec.process_decoded(d["entries"]),
            dec.process_decoded(d["working_ion_entry"]),
        )

    def as_dict_legacy(self):
        """
        Returns: MSONAble dict
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "entries": [entry.as_dict() for entry in self.get_all_entries()],
            "working_ion_entry": self.working_ion_entry.as_dict(),
        }


@dataclass
class InsertionVoltagePair(AbstractVoltagePair):
    """
    Defines an Insertion Voltage Pair.
    """

    entry_charge: ComputedEntry
    entry_discharge: ComputedEntry

    @classmethod
    def from_entries(cls, entry1, entry2, working_ion_entry):
        """
        Args:
            entry1: Entry corresponding to one of the entries in the voltage step.
            entry2: Entry corresponding to the other entry in the voltage step.
            working_ion_entry: A single ComputedEntry or PDEntry representing
                the element that carries charge across the battery, e.g. Li.
        """
        # initialize some internal variables
        working_element = working_ion_entry.composition.elements[0]

        entry_charge = entry1
        entry_discharge = entry2
        if entry_charge.composition.get_atomic_fraction(working_element) > entry2.composition.get_atomic_fraction(
            working_element
        ):
            (entry_charge, entry_discharge) = (entry_discharge, entry_charge)

        comp_charge = entry_charge.composition
        comp_discharge = entry_discharge.composition

        ion_sym = working_element.symbol

        frame_charge_comp = Composition({el: comp_charge[el] for el in comp_charge if el.symbol != ion_sym})
        frame_discharge_comp = Composition({el: comp_discharge[el] for el in comp_discharge if el.symbol != ion_sym})

        # Data validation

        # check that the ion is just a single element
        if not working_ion_entry.composition.is_element:
            raise ValueError("VoltagePair: The working ion specified must be an element")

        # check that at least one of the entries contains the working element
        if (
            not comp_charge.get_atomic_fraction(working_element) > 0
            and not comp_discharge.get_atomic_fraction(working_element) > 0
        ):
            raise ValueError("VoltagePair: The working ion must be present in one of the entries")

        # check that the entries do not contain the same amount of the working element
        if comp_charge.get_atomic_fraction(working_element) == comp_discharge.get_atomic_fraction(working_element):
            raise ValueError("VoltagePair: The working ion atomic percentage cannot be the same in both the entries")

        # check that the frameworks of the entries are equivalent
        if frame_charge_comp.reduced_formula != frame_discharge_comp.reduced_formula:
            raise ValueError("VoltagePair: the specified entries must have the same compositional framework")

        # Initialize normalization factors, charged and discharged entries

        valence_list = Element(ion_sym).oxidation_states
        working_ion_valence = abs(max(valence_list))

        (
            framework,
            norm_charge,
        ) = frame_charge_comp.get_reduced_composition_and_factor()
        norm_discharge = frame_discharge_comp.get_reduced_composition_and_factor()[1]

        # Initialize normalized properties
        if hasattr(entry_charge, "structure"):
            _vol_charge = entry_charge.structure.volume / norm_charge
        else:
            _vol_charge = entry_charge.data.get("volume")

        if hasattr(entry_discharge, "structure"):
            _vol_discharge = entry_discharge.structure.volume / norm_discharge
        else:
            _vol_discharge = entry_discharge.data.get("volume")

        comp_charge = entry_charge.composition
        comp_discharge = entry_discharge.composition

        _mass_charge = comp_charge.weight / norm_charge
        _mass_discharge = comp_discharge.weight / norm_discharge

        _num_ions_transferred = (comp_discharge[working_element] / norm_discharge) - (
            comp_charge[working_element] / norm_charge
        )

        _voltage = (
            ((entry_charge.energy / norm_charge) - (entry_discharge.energy / norm_discharge)) / _num_ions_transferred
            + working_ion_entry.energy_per_atom
        ) / working_ion_valence
        _mAh = _num_ions_transferred * Charge(1, "e").to("C") * Time(1, "s").to("h") * N_A * 1000 * working_ion_valence

        _frac_charge = comp_charge.get_atomic_fraction(working_element)
        _frac_discharge = comp_discharge.get_atomic_fraction(working_element)

        vpair = InsertionVoltagePair(  # pylint: disable=E1123
            voltage=_voltage,
            mAh=_mAh,
            mass_charge=_mass_charge,
            mass_discharge=_mass_discharge,
            vol_charge=_vol_charge,
            vol_discharge=_vol_discharge,
            frac_charge=_frac_charge,
            frac_discharge=_frac_discharge,
            working_ion_entry=working_ion_entry,
            entry_charge=entry_charge,
            entry_discharge=entry_discharge,
            framework_formula=framework.reduced_formula,
        )

        # Step 4: add (optional) hull and muO2 data
        vpair.decomp_e_charge = entry_charge.data.get("decomposition_energy", None)
        vpair.decomp_e_discharge = entry_discharge.data.get("decomposition_energy", None)

        vpair.muO2_charge = entry_charge.data.get("muO2", None)
        vpair.muO2_discharge = entry_discharge.data.get("muO2", None)

        return vpair

    def __repr__(self):
        output = [
            f"Insertion voltage pair with working ion {self.working_ion_entry.composition.reduced_formula}",
            f"V = {self.voltage}, mAh = {self.mAh}",
            f"mass_charge = {self.mass_charge}, mass_discharge = {self.mass_discharge}",
            f"vol_charge = {self.vol_charge}, vol_discharge = {self.vol_discharge}",
            f"frac_charge = {self.frac_charge}, frac_discharge = {self.frac_discharge}",
        ]
        return "\n".join(output)

    def __str__(self):
        return self.__repr__()
