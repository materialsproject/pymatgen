# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from scipy.constants import N_A

from pymatgen.core.periodic_table import Element
from pymatgen.core.units import Charge, Time

from pymatgen.analysis.reaction_calculator import BalancedReaction
from pymatgen.core.composition import Composition
from pymatgen.apps.battery.battery_abc import AbstractElectrode, \
    AbstractVoltagePair
from pymatgen.analysis.phase_diagram import PhaseDiagram
from monty.json import MontyDecoder
"""
This module contains the classes to build a ConversionElectrode.
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Feb 1, 2012"
__status__ = "Beta"


class ConversionElectrode(AbstractElectrode):
    """
    Class representing a ConversionElectrode.
    """

    def __init__(self, voltage_pairs, working_ion_entry, initial_comp):
        """
        General constructor for ConversionElectrode. However, it is usually
        easier to construct a ConversionElectrode using one of the static
        constructors provided.

        Args:
            voltage_pairs: The voltage pairs making up the Conversion
                Electrode.
            working_ion_entry: A single ComputedEntry or PDEntry
                representing the element that carries charge across the
                battery, e.g. Li.
            initial_comp: Starting composition for ConversionElectrode.
        """
        self._composition = initial_comp
        self._working_ion_entry = working_ion_entry
        ion_el = self._working_ion_entry.composition.elements[0]
        self._working_ion = ion_el.symbol
        self._vpairs = voltage_pairs

    @staticmethod
    def from_composition_and_pd(comp, pd, working_ion_symbol="Li", allow_unstable=False):
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
            elif e.is_element and \
                    e.composition.reduced_formula == working_ion_symbol:
                working_ion_entry = e

        if not allow_unstable and not entry:
            raise ValueError("Not stable compound found at composition {}."
                             .format(comp))

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
        vpairs = [ConversionVoltagePair.from_steps(profile[i], profile[i + 1],
                                                   normalization_els)
                  for i in range(len(profile) - 1)]
        return ConversionElectrode(vpairs, working_ion_entry, comp)

    @staticmethod
    def from_composition_and_entries(comp, entries_in_chemsys,
                                     working_ion_symbol="Li", allow_unstable=False):
        """
        Convenience constructor to make a ConversionElectrode from a
        composition and all entries in a chemical system.

        Args:
            comp: Starting composition for ConversionElectrode, e.g.,
                Composition("FeF3")
            entries_in_chemsys: Sequence containing all entries in a
               chemical system. E.g., all Li-Fe-F containing entries.
            working_ion_symbol: Element symbol of working ion. Defaults to Li.
        """
        pd = PhaseDiagram(entries_in_chemsys)
        return ConversionElectrode.from_composition_and_pd(comp, pd,
                                                           working_ion_symbol, allow_unstable)

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

        if adjacent_only:
            return [self.__class__(self._vpairs[i:i + 1],
                                   self._working_ion_entry, self._composition)
                    for i in range(len(self._vpairs))]
        sub_electrodes = []
        for i in range(len(self._vpairs)):
            for j in range(i, len(self._vpairs)):
                sub_electrodes.append(self.__class__(self._vpairs[i:j + 1],
                                                     self._working_ion_entry,
                                                     self._composition))
        return sub_electrodes

    @property
    def composition(self):
        return self._composition

    @property
    def working_ion(self):
        """
        The working ion as an Element object
        """
        return self._working_ion_entry.composition.elements[0]

    @property
    def working_ion_entry(self):
        return self._working_ion_entry

    @property
    def voltage_pairs(self):
        return self._vpairs

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
            all_formulas1 = set([rxn1.all_comp[i].reduced_formula
                                 for i in range(len(rxn1.all_comp))
                                 if abs(rxn1.coeffs[i]) > 1e-5])
            for pair2 in self:
                rxn2 = pair2.rxn
                all_formulas2 = set([rxn2.all_comp[i].reduced_formula
                                     for i in range(len(rxn2.all_comp))
                                     if abs(rxn2.coeffs[i]) > 1e-5])
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
            all_formulas1 = set([rxn1.all_comp[i].reduced_formula
                                 for i in range(len(rxn1.all_comp))
                                 if abs(rxn1.coeffs[i]) > 1e-5])
            for pair2 in self:
                rxn2 = pair2.rxn
                all_formulas2 = set([rxn2.all_comp[i].reduced_formula
                                     for i in range(len(rxn2.all_comp))
                                     if abs(rxn2.coeffs[i]) > 1e-5])
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
        output = ["Conversion electrode with formula {} and nsteps {}"
                  .format(self._composition.reduced_formula, self.num_steps),
                  "Avg voltage {} V, min voltage {} V, max voltage {} V"
                  .format(self.get_average_voltage(), self.min_voltage,
                          self.max_voltage),
                  "Capacity (grav.) {} mAh/g, capacity (vol.) {} Ah/l"
                  .format(self.get_capacity_grav(),
                          self.get_capacity_vol()),
                  "Specific energy {} Wh/kg, energy density {} Wh/l"
                  .format(self.get_specific_energy(),
                          self.get_energy_density())]
        return "\n".join(output)

    @classmethod
    def from_dict(cls, d):
        dec = MontyDecoder()
        return cls(dec.process_decoded(d["voltage_pairs"]),
                   dec.process_decoded(d["working_ion_entry"]),
                   Composition(d["initial_comp"]))

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "voltage_pairs": [v.as_dict() for v in self._vpairs],
                "working_ion_entry": self.working_ion_entry.as_dict(),
                "initial_comp": self._composition.as_dict()}

    def get_summary_dict(self, print_subelectrodes=True):
        """
        Args:
            print_subelectrodes:
                Also print data on all the possible subelectrodes

        Returns:
            a summary of this electrode"s properties in dictionary format
        """

        d = {}
        framework_comp = Composition({k: v
                                      for k, v in self._composition.items()
                                      if k.symbol != self.working_ion.symbol})

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
        for pair in self._vpairs:
            rxn = pair.rxn
            frac.append(pair.frac_charge)
            frac.append(pair.frac_discharge)
            d["reactions"].append(str(rxn))
            for i in range(len(rxn.coeffs)):
                if abs(rxn.coeffs[i]) > 1e-5 and rxn.all_comp[i] not in comps:
                    comps.append(rxn.all_comp[i])
                if abs(rxn.coeffs[i]) > 1e-5 and \
                        rxn.all_comp[i].reduced_formula != d["working_ion"]:
                    reduced_comp = rxn.all_comp[i].reduced_composition
                    comp_dict = reduced_comp.as_dict()
                    d["reactant_compositions"].append(comp_dict)
        d["fracA_charge"] = min(frac)
        d["fracA_discharge"] = max(frac)
        d["nsteps"] = self.num_steps
        if print_subelectrodes:
            f_dict = lambda c: c.get_summary_dict(print_subelectrodes=False)
            d["adj_pairs"] = list(map(f_dict,
                                 self.get_sub_electrodes(adjacent_only=True)))
            d["all_pairs"] = list(map(f_dict,
                                 self.get_sub_electrodes(adjacent_only=False)))
        return d


class ConversionVoltagePair(AbstractVoltagePair):
    """
    A VoltagePair representing a Conversion Reaction with a defined voltage.
    Typically not initialized directly but rather used by ConversionElectrode.

    Args:
        balanced_rxn (BalancedReaction): BalancedReaction for the step
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

    def __init__(self, balanced_rxn, voltage, mAh, vol_charge, vol_discharge,
                 mass_charge, mass_discharge, frac_charge, frac_discharge,
                 entries_charge, entries_discharge, working_ion_entry):
        self._working_ion_entry = working_ion_entry
        working_ion = self._working_ion_entry.composition.elements[0].symbol
        self._voltage = voltage
        self._mAh = mAh
        self._vol_charge = vol_charge
        self._mass_charge = mass_charge
        self._mass_discharge = mass_discharge
        self._vol_discharge = vol_discharge
        self._frac_charge = frac_charge
        self._frac_discharge = frac_discharge

        self._rxn = balanced_rxn
        self._working_ion = working_ion
        self._entries_charge = entries_charge
        self._entries_discharge = entries_discharge

    @staticmethod
    def from_steps(step1, step2, normalization_els):
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
        voltage = (-step1["chempot"] + working_ion_entry.energy_per_atom)/working_ion_valence
        mAh = (step2["evolution"] - step1["evolution"]) \
            * Charge(1, "e").to("C") * Time(1, "s").to("h") * N_A * 1000*working_ion_valence
        licomp = Composition(working_ion)
        prev_rxn = step1["reaction"]
        reactants = {comp: abs(prev_rxn.get_coeff(comp))
                     for comp in prev_rxn.products if comp != licomp}

        curr_rxn = step2["reaction"]
        products = {comp: abs(curr_rxn.get_coeff(comp))
                    for comp in curr_rxn.products if comp != licomp}

        reactants[licomp] = (step2["evolution"] - step1["evolution"])

        rxn = BalancedReaction(reactants, products)

        for el, amt in normalization_els.items():
            if rxn.get_el_amount(el) > 1e-6:
                rxn.normalize_to_element(el, amt)
                break

        prev_mass_dischg = sum([prev_rxn.all_comp[i].weight
                                * abs(prev_rxn.coeffs[i])
                                for i in range(len(prev_rxn.all_comp))]) / 2
        vol_charge = sum([abs(prev_rxn.get_coeff(e.composition))
                          * e.structure.volume
                          for e in step1["entries"]
                          if e.composition.reduced_formula != working_ion])
        mass_discharge = sum([curr_rxn.all_comp[i].weight
                              * abs(curr_rxn.coeffs[i])
                              for i in range(len(curr_rxn.all_comp))]) / 2
        mass_charge = prev_mass_dischg
        mass_discharge = mass_discharge
        vol_discharge = sum([abs(curr_rxn.get_coeff(e.composition))
                             * e.structure.volume
                             for e in step2["entries"]
                             if e.composition.reduced_formula != working_ion])

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

        return ConversionVoltagePair(rxn, voltage, mAh, vol_charge,
                                     vol_discharge, mass_charge,
                                     mass_discharge,
                                     frac_charge, frac_discharge,
                                     entries_charge, entries_discharge,
                                     working_ion_entry)

    @property
    def working_ion(self):
        return self._working_ion

    @property
    def entries_charge(self):
        return self._entries_charge

    @property
    def entries_discharge(self):
        return self._entries_discharge

    @property
    def frac_charge(self):
        return self._frac_charge

    @property
    def frac_discharge(self):
        return self._frac_discharge

    @property
    def rxn(self):
        return self._rxn

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
        output = ["Conversion voltage pair with working ion {}"
                  .format(self._working_ion_entry.composition.reduced_formula),
                  "Reaction : {}".format(self._rxn),
                  "V = {}, mAh = {}".format(self.voltage, self.mAh),
                  "frac_charge = {}, frac_discharge = {}"
                  .format(self.frac_charge, self.frac_discharge),
                  "mass_charge = {}, mass_discharge = {}"
                  .format(self.mass_charge, self.mass_discharge),
                  "vol_charge = {}, vol_discharge = {}"
                  .format(self.vol_charge, self.vol_discharge)]
        return "\n".join(output)

    def __str__(self):
        return self.__repr__()

    @classmethod
    def from_dict(cls, d):
        dec = MontyDecoder()
        working_ion_entry = dec.process_decoded(d["working_ion_entry"])
        balanced_rxn = dec.process_decoded(d["balanced_rxn"])
        entries_charge = dec.process_decoded(d["entries_charge"])
        entries_discharge = dec.process_decoded(d["entries_discharge"])
        return ConversionVoltagePair(balanced_rxn, d["voltage"], d["mAh"],
                                   d["vol_charge"], d["vol_discharge"],
                                   d["mass_charge"], d["mass_discharge"],
                                   d["frac_charge"], d["frac_discharge"],
                                   entries_charge, entries_discharge,
                                   working_ion_entry)

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "working_ion_entry": self._working_ion_entry.as_dict(),
                "voltage": self._voltage, "mAh": self._mAh,
                "vol_charge": self._vol_charge,
                "mass_charge": self._mass_charge,
                "mass_discharge": self._mass_discharge,
                "vol_discharge": self._vol_discharge,
                "frac_charge": self._frac_charge,
                "frac_discharge": self._frac_discharge,
                "balanced_rxn": self._rxn.as_dict(),
                "entries_charge": [e.as_dict() for e in self._entries_charge],
                "entries_discharge": [e.as_dict() for e in
                                      self._entries_discharge]}
