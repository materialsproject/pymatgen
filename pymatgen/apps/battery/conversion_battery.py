# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module contains the classes to build a ConversionElectrode.
"""
from typing import Iterable, Dict
from dataclasses import dataclass

from monty.dev import deprecated
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
        _initial_comp_formula: Starting composition for ConversionElectrode represented
            as a string/formula.
    """

    _initial_comp_formula: str

    @property
    def initial_comp(self) -> Composition:
        """
        The pymatgen Composition representation of the initial composition
        """
        return Composition(self._initial_comp_formula)

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
            raise ValueError("Not stable compound found at composition {}.".format(comp))

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

        return cls(
            voltage_pairs=vpairs,
            working_ion_entry=working_ion_entry,
            _initial_comp_formula=comp.reduced_formula,
            _framework_formula=framework.reduced_formula,
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
        return cls.from_composition_and_pd(comp, pd, working_ion_symbol, allow_unstable)

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
        # _initial_comp_formula = comp.reduced_formula, _framework_formula = framework.reduced_formula
        if adjacent_only:
            return [
                self.__class__(
                    voltage_pairs=self.voltage_pairs[i : i + 1],
                    working_ion_entry=self.working_ion_entry,
                    _initial_comp_formula=self._initial_comp_formula,
                    _framework_formula=self._framework_formula,
                )
                for i in range(len(self.voltage_pairs))
            ]
        sub_electrodes = []
        for i in range(len(self.voltage_pairs)):
            for j in range(i, len(self.voltage_pairs)):
                sub_electrodes.append(
                    self.__class__(
                        voltage_pairs=self.voltage_pairs[i : j + 1],
                        working_ion_entry=self.working_ion_entry,
                        _initial_comp_formula=self._initial_comp_formula,
                        _framework_formula=self._framework_formula,
                    )
                )
        return sub_electrodes

    def is_super_electrode(self, conversion_electrode):
        """
        Checks if a particular conversion electrode is a sub electrode of the
        current electrode. Starting from a more lithiated state may result in
        a subelectrode that is essentially on the same path.  For example, a
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
            "Conversion electrode with formula {} and nsteps {}".format(
                self.initial_comp.reduced_formula, self.num_steps
            ),
            "Avg voltage {} V, min voltage {} V, max voltage {} V".format(
                self.get_average_voltage(), self.min_voltage, self.max_voltage
            ),
            "Capacity (grav.) {} mAh/g, capacity (vol.) {} Ah/l".format(
                self.get_capacity_grav(), self.get_capacity_vol()
            ),
            "Specific energy {} Wh/kg, energy density {} Wh/l".format(
                self.get_specific_energy(), self.get_energy_density()
            ),
        ]
        return "\n".join(output)

    def get_summary_dict(self, print_subelectrodes=True) -> Dict:
        """
        Generate a summary dict.
        Populates the summary dict with the basic information from the parent method then populates more information.
        Since the parent method calls self.get_summary_dict(print_subelectrodes=True) for the subelectrodes.
        The current methode will be called from within super().get_summary_dict.

        Args:
            print_subelectrodes: Also print data on all the possible
                subelectrodes.

        Returns:
            A summary of this electrode"s properties in dict format.
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

    @deprecated(
        replacement=get_summary_dict,
        message="Name and logic changed, will be as_dict_summary will be removed in the futurn.",
    )
    def as_dict_summary(self, print_subelectrodes=True):
        """
        Args:
            print_subelectrodes:
                Also print data on all the possible subelectrodes

        Returns:
            a summary of this electrode"s properties in dictionary format
        """

        d = {}
        framework_comp = Composition(
            {k: v for k, v in self.initial_comp.items() if k.symbol != self.working_ion.symbol}
        )

        d["framework"] = framework_comp.to_data_dict
        d["framework_pretty"] = framework_comp.reduced_formula
        d["average_voltage"] = self.get_average_voltage()
        d["max_voltage"] = self.max_voltage
        d["min_voltage"] = self.min_voltage
        d["max_delta_volume"] = self.max_delta_volume
        d["max_instability"] = 0
        d["max_voltage_step"] = self.max_voltage_step
        d["nsteps"] = self.num_steps
        d["capacity_grav"] = self.get_capacity_grav()
        d["capacity_vol"] = self.get_capacity_vol()
        d["energy_grav"] = self.get_specific_energy()
        d["energy_vol"] = self.get_energy_density()
        d["working_ion"] = self.working_ion.symbol
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
        d["fracA_charge"] = min(frac)
        d["fracA_discharge"] = max(frac)
        d["nsteps"] = self.num_steps
        if print_subelectrodes:

            def f_dict(c):
                return c.get_summary_dict(print_subelectrodes=False)

            d["adj_pairs"] = list(map(f_dict, self.get_sub_electrodes(adjacent_only=True)))
            d["all_pairs"] = list(map(f_dict, self.get_sub_electrodes(adjacent_only=False)))
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
            sum([prev_rxn.all_comp[i].weight * abs(prev_rxn.coeffs[i]) for i in range(len(prev_rxn.all_comp))]) / 2
        )
        vol_charge = sum(
            [
                abs(prev_rxn.get_coeff(e.composition)) * e.structure.volume
                for e in step1["entries"]
                if e.composition.reduced_formula != working_ion
            ]
        )
        mass_discharge = (
            sum([curr_rxn.all_comp[i].weight * abs(curr_rxn.coeffs[i]) for i in range(len(curr_rxn.all_comp))]) / 2
        )
        mass_charge = prev_mass_dischg
        mass_discharge = mass_discharge
        vol_discharge = sum(
            [
                abs(curr_rxn.get_coeff(e.composition)) * e.structure.volume
                for e in step2["entries"]
                if e.composition.reduced_formula != working_ion
            ]
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
        entries_charge = step2["entries"]
        entries_discharge = step1["entries"]

        return cls(
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
            _framework_formula=framework_formula,
        )

    def __repr__(self):
        output = [
            "Conversion voltage pair with working ion {}".format(self.working_ion_entry.composition.reduced_formula),
            "Reaction : {}".format(self.rxn),
            "V = {}, mAh = {}".format(self.voltage, self.mAh),
            "frac_charge = {}, frac_discharge = {}".format(self.frac_charge, self.frac_discharge),
            "mass_charge = {}, mass_discharge = {}".format(self.mass_charge, self.mass_discharge),
            "vol_charge = {}, vol_discharge = {}".format(self.vol_charge, self.vol_discharge),
        ]
        return "\n".join(output)

    def __str__(self):
        return self.__repr__()
