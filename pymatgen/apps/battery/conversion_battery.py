# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module contains the classes to build a ConversionElectrode.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from scipy.constants import N_A

from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.reaction_calculator import BalancedReaction
from pymatgen.apps.battery.battery_abc import AbstractElectrode, AbstractVoltagePair
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.core.units import Charge, Time
from pymatgen.entries.computed_entries import ComputedEntry


@dataclass
class ConversionElectrode(AbstractElectrode):
    """
    Class representing a ConversionElectrode, since it is dataclass
    this object can be constructed for the attributes.
    However, it is usually easier to construct a ConversionElectrode using one of the classmethods
    constructors provided.

    Attribute:
        voltage_pairs: The voltage pairs making up the Conversion Electrode.
        working_ion_entry: A single ComputedEntry or PDEntry
            representing the element that carries charge across the
            battery, e.g. Li.
        initial_comp_formula: Starting composition for ConversionElectrode represented
            as a string/formula.
    """

    initial_comp_formula: str

    @property
    def initial_comp(self) -> Composition:
        """
        The pymatgen Composition representation of the initial composition
        """
        return Composition(self.initial_comp_formula)

    @classmethod
    def from_composition_and_pd(cls, comp, pd, working_ion_symbol="Li", allow_unstable=False):
        """
        Convenience constructor to make a ConversionElectrode from a
        composition and a phase diagram.

        Args:
            comp:
                Starting composition for ConversionElectrode, e.g.,
                Composition("FeF3")
            pd:
                A PhaseDiagram of the relevant system (e.g., Li-Fe-F)
            working_ion_symbol:
                Element symbol of working ion. Defaults to Li.
            allow_unstable:
                Allow compositions that are unstable
        """
        working_ion = Element(working_ion_symbol)
        entry = None
        working_ion_entry = None
        for e in pd.stable_entries:
            if e.composition.reduced_formula == comp.reduced_formula:
                entry = e
            elif e.is_element and e.composition.reduced_formula == working_ion_symbol:
                working_ion_entry = e

        if not allow_unstable and not entry:
            raise ValueError(f"Not stable compound found at composition {comp}.")

        profile = pd.get_element_profile(working_ion, comp)
        # Need to reverse because voltage goes form most charged to most
        # discharged.
        profile.reverse()
        if len(profile) < 2:
            return None
        working_ion_entry = working_ion_entry
        working_ion = working_ion_entry.composition.elements[0].symbol
        normalization_els = {}
        for el, amt in comp.items():
            if el != Element(working_ion):
                normalization_els[el] = amt
        framework = comp.as_dict()
        if working_ion in framework:
            framework.pop(working_ion)
        framework = Composition(framework)

        vpairs = [
            ConversionVoltagePair.from_steps(
                profile[i],
                profile[i + 1],
                normalization_els,
                framework_formula=framework.reduced_formula,
            )
            for i in range(len(profile) - 1)
        ]

        return ConversionElectrode(  # pylint: disable=E1123
            voltage_pairs=vpairs,
            working_ion_entry=working_ion_entry,
            initial_comp_formula=comp.reduced_formula,
            framework_formula=framework.reduced_formula,
        )

    @classmethod
    def from_composition_and_entries(cls, comp, entries_in_chemsys, working_ion_symbol="Li", allow_unstable=False):
        """
        Convenience constructor to make a ConversionElectrode from a
        composition and all entries in a chemical system.

        Args:
            comp: Starting composition for ConversionElectrode, e.g.,
                Composition("FeF3")
            entries_in_chemsys: Sequence containing all entries in a
               chemical system. E.g., all Li-Fe-F containing entries.
            working_ion_symbol: Element symbol of working ion. Defaults to Li.
            allow_unstable: If True, allow any composition to be used as the
                    starting point of a conversion voltage curve, this is useful
                    for comparing with insertion electrodes
        """
        pd = PhaseDiagram(entries_in_chemsys)
        return ConversionElectrode.from_composition_and_pd(comp, pd, working_ion_symbol, allow_unstable)

    def get_sub_electrodes(self, adjacent_only=True):
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
                will have multiple voltage steps if this is set true

        Returns:
            A list of ConversionElectrode objects
        """

        # voltage_pairs = vpairs, working_ion_entry = working_ion_entry,
        # _initial_comp_formula = comp.reduced_formula, framework_formula = framework.reduced_formula
        if adjacent_only:
            return [
                ConversionElectrode(  # pylint: disable=E1123
                    voltage_pairs=self.voltage_pairs[i : i + 1],
                    working_ion_entry=self.working_ion_entry,
                    initial_comp_formula=self.initial_comp_formula,
                    framework_formula=self.framework_formula,
                )
                for i in range(len(self.voltage_pairs))
            ]
        sub_electrodes = []
        for i in range(len(self.voltage_pairs)):
            for j in range(i, len(self.voltage_pairs)):
                sub_electrodes.append(
                    ConversionElectrode(  # pylint: disable=E1123
                        voltage_pairs=self.voltage_pairs[i : j + 1],
                        working_ion_entry=self.working_ion_entry,
                        initial_comp_formula=self.initial_comp_formula,
                        framework_formula=self.framework_formula,
                    )
                )
        return sub_electrodes

    def is_super_electrode(self, conversion_electrode):
        """
        Checks if a particular conversion electrode is a sub electrode of the
        current electrode. Starting from a more lithiated state may result in
        a subelectrode that is essentially on the same path. For example, a
        ConversionElectrode formed by starting from an FePO4 composition would
        be a super_electrode of a ConversionElectrode formed from an LiFePO4
        composition.
        """
        for pair1 in conversion_electrode:
            rxn1 = pair1.rxn
            all_formulas1 = {
                rxn1.all_comp[i].reduced_formula for i in range(len(rxn1.all_comp)) if abs(rxn1.coeffs[i]) > 1e-5
            }
            for pair2 in self:
                rxn2 = pair2.rxn
                all_formulas2 = {
                    rxn2.all_comp[i].reduced_formula for i in range(len(rxn2.all_comp)) if abs(rxn2.coeffs[i]) > 1e-5
                }
                if all_formulas1 == all_formulas2:
                    break
            else:
                return False
        return True

    def __eq__(self, conversion_electrode):
        """
        Check if two electrodes are exactly the same:
        """
        if len(self) != len(conversion_electrode):
            return False

        for pair1 in conversion_electrode:
            rxn1 = pair1.rxn
            all_formulas1 = {
                rxn1.all_comp[i].reduced_formula for i in range(len(rxn1.all_comp)) if abs(rxn1.coeffs[i]) > 1e-5
            }
            for pair2 in self:
                rxn2 = pair2.rxn
                all_formulas2 = {
                    rxn2.all_comp[i].reduced_formula for i in range(len(rxn2.all_comp)) if abs(rxn2.coeffs[i]) > 1e-5
                }
                if all_formulas1 == all_formulas2:
                    break
            else:
                return False
        return True

    def __hash__(self):
        return 7

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        output = [
            f"Conversion electrode with formula {self.initial_comp.reduced_formula} and nsteps {self.num_steps}",
            f"Avg voltage {self.get_average_voltage()} V, min voltage {self.min_voltage} V, "
            f"max voltage {self.max_voltage} V",
            f"Capacity (grav.) {self.get_capacity_grav()} mAh/g, capacity (vol.) {self.get_capacity_vol()} Ah/l",
            f"Specific energy {self.get_specific_energy()} Wh/kg, energy density {self.get_energy_density()} Wh/l",
        ]
        return "\n".join(output)

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
        d["reactions"] = []
        d["reactant_compositions"] = []
        comps = []
        frac = []
        for pair in self.voltage_pairs:
            rxn = pair.rxn
            frac.append(pair.frac_charge)
            frac.append(pair.frac_discharge)
            d["reactions"].append(str(rxn))
            for i, v in enumerate(rxn.coeffs):
                if abs(v) > 1e-5 and rxn.all_comp[i] not in comps:
                    comps.append(rxn.all_comp[i])
                if abs(v) > 1e-5 and rxn.all_comp[i].reduced_formula != d["working_ion"]:
                    reduced_comp = rxn.all_comp[i].reduced_composition
                    comp_dict = reduced_comp.as_dict()
                    d["reactant_compositions"].append(comp_dict)
        return d


@dataclass
class ConversionVoltagePair(AbstractVoltagePair):
    """
    A VoltagePair representing a Conversion Reaction with a defined voltage.
    Typically not initialized directly but rather used by ConversionElectrode.
    Attributes:
        rxn (BalancedReaction): BalancedReaction for the step
        voltage (float): Voltage for the step
        mAh (float): Capacity of the step
        vol_charge (float): Volume of charged state
        vol_discharge (float): Volume of discharged state
        mass_charge (float): Mass of charged state
        mass_discharge (float): Mass of discharged state
        frac_charge (float): Fraction of working ion in the charged state
        frac_discharge (float): Fraction of working ion in the discharged state
        entries_charge ([ComputedEntry]): Entries in the charged state
        entries_discharge ([ComputedEntry]): Entries in discharged state
        working_ion_entry (ComputedEntry): Entry of the working ion.
    """

    rxn: BalancedReaction
    entries_charge: Iterable[ComputedEntry]
    entries_discharge: Iterable[ComputedEntry]

    @classmethod
    def from_steps(cls, step1, step2, normalization_els, framework_formula=None):
        """
        Creates a ConversionVoltagePair from two steps in the element profile
        from a PD analysis.

        Args:
            step1: Starting step
            step2: Ending step
            normalization_els: Elements to normalize the reaction by. To
                ensure correct capacities.
        """
        working_ion_entry = step1["element_reference"]
        working_ion = working_ion_entry.composition.elements[0].symbol
        working_ion_valence = max(Element(working_ion).oxidation_states)
        voltage = (-step1["chempot"] + working_ion_entry.energy_per_atom) / working_ion_valence
        mAh = (
            (step2["evolution"] - step1["evolution"])
            * Charge(1, "e").to("C")
            * Time(1, "s").to("h")
            * N_A
            * 1000
            * working_ion_valence
        )
        licomp = Composition(working_ion)
        prev_rxn = step1["reaction"]
        reactants = {comp: abs(prev_rxn.get_coeff(comp)) for comp in prev_rxn.products if comp != licomp}

        curr_rxn = step2["reaction"]
        products = {comp: abs(curr_rxn.get_coeff(comp)) for comp in curr_rxn.products if comp != licomp}

        reactants[licomp] = step2["evolution"] - step1["evolution"]

        rxn = BalancedReaction(reactants, products)

        for el, amt in normalization_els.items():
            if rxn.get_el_amount(el) > 1e-6:
                rxn.normalize_to_element(el, amt)
                break

        prev_mass_dischg = (
            sum(prev_rxn.all_comp[i].weight * abs(prev_rxn.coeffs[i]) for i in range(len(prev_rxn.all_comp))) / 2
        )
        vol_charge = sum(
            abs(prev_rxn.get_coeff(e.composition)) * e.structure.volume
            for e in step1["entries"]
            if e.composition.reduced_formula != working_ion
        )
        mass_discharge = (
            sum(curr_rxn.all_comp[i].weight * abs(curr_rxn.coeffs[i]) for i in range(len(curr_rxn.all_comp))) / 2
        )
        mass_charge = prev_mass_dischg
        mass_discharge = mass_discharge
        vol_discharge = sum(
            abs(curr_rxn.get_coeff(e.composition)) * e.structure.volume
            for e in step2["entries"]
            if e.composition.reduced_formula != working_ion
        )

        totalcomp = Composition({})
        for comp in prev_rxn.products:
            if comp.reduced_formula != working_ion:
                totalcomp += comp * abs(prev_rxn.get_coeff(comp))
        frac_charge = totalcomp.get_atomic_fraction(Element(working_ion))

        totalcomp = Composition({})
        for comp in curr_rxn.products:
            if comp.reduced_formula != working_ion:
                totalcomp += comp * abs(curr_rxn.get_coeff(comp))
        frac_discharge = totalcomp.get_atomic_fraction(Element(working_ion))

        rxn = rxn
        entries_charge = step1["entries"]
        entries_discharge = step2["entries"]

        return ConversionVoltagePair(  # pylint: disable=E1123
            rxn=rxn,
            voltage=voltage,
            mAh=mAh,
            vol_charge=vol_charge,
            vol_discharge=vol_discharge,
            mass_charge=mass_charge,
            mass_discharge=mass_discharge,
            frac_charge=frac_charge,
            frac_discharge=frac_discharge,
            entries_charge=entries_charge,
            entries_discharge=entries_discharge,
            working_ion_entry=working_ion_entry,
            framework_formula=framework_formula,
        )

    def __repr__(self):
        output = [
            f"Conversion voltage pair with working ion {self.working_ion_entry.composition.reduced_formula}",
            f"Reaction : {self.rxn}",
            f"V = {self.voltage}, mAh = {self.mAh}",
            f"frac_charge = {self.frac_charge}, frac_discharge = {self.frac_discharge}",
            f"mass_charge = {self.mass_charge}, mass_discharge = {self.mass_discharge}",
            f"vol_charge = {self.vol_charge}, vol_discharge = {self.vol_discharge}",
        ]
        return "\n".join(output)

    def __str__(self):
        return self.__repr__()
