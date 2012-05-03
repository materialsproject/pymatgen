#!/usr/bin/env python

'''
This module is used for analysis of materials with potential application as
intercalation batteries.
'''

from __future__ import division

__author__ = "Anubhav Jain, Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Jan 13, 2012"
__status__ = "Beta"

import itertools

from pymatgen.core.structure import Composition
from pymatgen.core.physical_constants import ELECTRON_TO_AMPERE_HOURS
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.apps.battery.battery_abc import AbstractElectrode, AbstractVoltagePair


class InsertionElectrode(AbstractElectrode):
    '''
    A set of topotactically related compounds, with different amounts of a
    single element, e.g. TiO2 and LiTiO2, that can be used to define an
    insertion battery electrode.
    '''

    def __init__(self, entries, working_ion_entry):
        '''
        Create a new InsertionElectrode.
        
        Args:
            entries:
                A list of ComputedStructureEntries (or subclasses) representing
                the different topotactic states of the battery, e.g. TiO2 and
                LiTiO2.
            working_ion_entry:
                A single ComputedEntry or PDEntry representing the element that
                carries charge across the battery, e.g. Li.
        '''

        self._entries = entries
        self._working_ion = working_ion_entry.composition.elements[0]
        self._working_ion_entry = working_ion_entry

        #Prepare to make phase diagram: determine elements and set their energy to be very high
        elements = set()
        map(elements.update, [entry.composition.elements for entry in entries])

        """
        Set an artifical energy for each element for convex hull generation
        """
        element_energy = max([entry.energy_per_atom for entry in entries]) + 10

        pdentries = []
        pdentries.extend(entries)
        pdentries.extend([PDEntry(Composition({el:1}), element_energy) for el in elements])

        #Make phase diagram to determine which entries are stable vs. unstable
        pd = PhaseDiagram(pdentries)

        lifrac = lambda entry: entry.composition.get_atomic_fraction(self._working_ion)

        #stable entries ordered by amount of Li asc
        self._stable_entries = tuple(sorted([e for e in pd.stable_entries if e in entries], key=lifrac))

        #unstable entries ordered by amount of Li asc
        self._unstable_entries = tuple(sorted([e for e in pd.unstable_entries if e in entries], key=lifrac))

        #create voltage pairs
        self._vpairs = tuple([InsertionVoltagePair(self._stable_entries[i], self._stable_entries[i + 1], working_ion_entry) for i in range(len(self._stable_entries) - 1)])

    @property
    def working_ion(self):
        '''
        The working ion as an Element object
        '''
        return self._working_ion

    @property
    def working_ion_entry(self):
        return self._working_ion_entry

    @property
    def voltage_pairs(self):
        return self._vpairs

    def get_stable_entries(self, charge_to_discharge=True):
        '''   
        Args:
            charge_to_discharge:
                order from most charge to most discharged state? Default to
                True.
        
        Returns:
            A list of stable entries in the electrode, ordered by amount of the
            working ion.
        '''
        list_copy = list(self._stable_entries)
        return list_copy if charge_to_discharge else list_copy.reverse()

    def get_unstable_entries(self, charge_to_discharge=True):
        '''
        Returns the unstable entries for the electrode.
           
        Args:
            charge_to_discharge:
                order from most charge to most discharged state? Defaults to
                True.
        
        Returns:
            A list of unstable entries in the electrode, ordered by amount of
            the working ion.
        '''
        list_copy = list(self._unstable_entries)
        return list_copy if charge_to_discharge else list_copy.reverse()

    def get_all_entries(self, charge_to_discharge=True):
        '''
        Return all entries input for the electrode.
        
        Args:
            charge_to_discharge:
                order from most charge to most discharged state? Defaults to
                True.
        
        Returns:
            A list of all entries in the electrode (both stable and unstable),
            ordered by amount of the working ion.
        '''
        all_entries = list(self.stable_entries())
        all_entries.extend(self.unstable_entries())
        #sort all entries by amount of working ion ASC
        all_entries = sorted([e for e in all_entries], key=lambda entry: entry.composition.get_atomic_fraction(self.working_ion))
        return all_entries if charge_to_discharge else all_entries.reverse()

    @property
    def fully_charged_entry(self):
        '''
        The most charged entry along the topotactic path.
        '''
        return self._stable_entries[0]

    @property
    def fully_discharged_entry(self):
        '''
        The most discharged entry along the topotactic path.
        '''
        return self._stable_entries[-1]

    def get_max_instability(self, min_voltage=None, max_voltage=None):
        '''
        The maximum instability along a path for a specific voltage range.
        
        Args:
            min_voltage:
                The minimum allowable voltage.
            max_voltage:
                The maximum allowable voltage.
                
        Returns:
            Maximum decomposition energy of all compounds along the insertion
            path (a subset of the path can be chosen by the optional arguments)
        '''
        return max([max(pair.decomp_e_discharge, pair.decomp_e_charge) for pair in self._select_in_voltage_range(min_voltage, max_voltage)])

    def get_min_instability(self, min_voltage=None, max_voltage=None):
        '''
        The minimum instability along a path for a specific voltage range.
        
        Args:
            min_voltage:
                The minimum allowable voltage.
            max_voltage:
                The maximum allowable voltage.
                
        Returns:
            Minimum decomposition energy of all compounds along the insertion
            path (a subset of the path can be chosen by the optional arguments)
        '''
        return min([min(pair.decomp_e_discharge, pair.decomp_e_charge) for pair in self._select_in_voltage_range(min_voltage, max_voltage)])

    def get_max_muO2(self, min_voltage=None, max_voltage=None):
        '''
        Maximum critical oxygen chemical potential along path.
        
        Args:
            min_voltage:
                The minimum allowable voltage.
            max_voltage:
                The maximum allowable voltage.
                
        Returns:
            Maximum critical oxygen chemical of all compounds along the
            insertion path (a subset of the path can be chosen by the optional
            arguments).
        '''
        return max([max(pair.muO2_discharge, pair.muO2_charge) for pair in self._select_in_voltage_range(min_voltage, max_voltage)])

    def get_min_muO2(self, min_voltage=None, max_voltage=None):
        '''
        Minimum critical oxygen chemical potential along path.
        
        Args:
            min_voltage:
                the minimum allowable voltage for a given step
            max_voltage:
                the maximum allowable voltage allowable for a given step
                
        Returns:
            Minimum critical oxygen chemical of all compounds along the
            insertion path (a subset of the path can be chosen by the optional
            arguments).
        '''
        return min([min(pair.muO2_discharge, pair.muO2_charge) for pair in self._select_in_voltage_range(min_voltage, max_voltage)])

    def get_sub_electrodes(self, adjacent_only=True, include_myself=True):
        '''
        If this electrode contains multiple voltage steps, then it is possible
        to use only a subset of the voltage steps to define other electrodes.
        For example, an LiTiO2 electrode might contain three subelectrodes:
        [LiTiO2 --> TiO2, LiTiO2 --> Li0.5TiO2, Li0.5TiO2 --> TiO2]
        This method can be used to return all the subelectrodes with some
        options
        
        Args:
            adjacent_only:
                Only return electrodes from compounds that are adjacent on the
                convex hull, i.e. no electrodes returned will have multiple
                voltage steps if this is set True.
            include_myself:
                Include this identical electrode in the list of results.
        
        Returns:
            A list of InsertionElectrode objects
        '''
        battery_list = []
        pair_it = self._vpairs if adjacent_only else itertools.combinations_with_replacement(self._vpairs, 2)

        for pair in pair_it:
            entry_charge = pair.entry_charge if adjacent_only else pair[0].entry_charge
            entry_discharge = pair.entry_discharge if adjacent_only else pair[1].entry_discharge

            def in_range(entry):
                frac = entry.composition.get_atomic_fraction(self.working_ion)
                return frac >= entry_charge.composition.get_atomic_fraction(self.working_ion) and frac <= entry_discharge.composition.get_atomic_fraction(self.working_ion)

            if include_myself or entry_charge != self.fully_charged_entry or entry_discharge != self.fully_discharged_entry:
                unstable_entries = filter(in_range, self.get_unstable_entries())
                stable_entries = filter(in_range, self.get_stable_entries())
                all_entries = list(stable_entries)
                all_entries.extend(unstable_entries)
                battery_list.append(InsertionElectrode(all_entries, self.working_ion_entry))

        return battery_list

    def to_dict_summary(self, print_subelectrodes=True):
        '''
        Arguments:
            print_subelectrodes:
                Also print data on all the possible subelectrodes
        
        Returns:
            A summary of this electrode's properties in dictionary format.
        '''

        d = {}
        d['average_voltage'] = self.get_average_voltage()
        d['max_voltage'] = self.max_voltage
        d['min_voltage'] = self.min_voltage
        d['max_delta_volume'] = self.max_delta_volume
        d['max_voltage_step'] = self.max_voltage_step
        d['capacity_grav'] = self.get_capacity_grav()
        d['capacity_vol'] = self.get_capacity_vol()
        d['energy_grav'] = self.get_specific_energy()
        d['energy_vol'] = self.get_energy_density()
        d['working_ion'] = self._working_ion.symbol
        d['nsteps'] = self.num_steps
        d['framework'] = self._vpairs[0].framework.to_data_dict
        d['formula_charge'] = self.fully_charged_entry.composition.reduced_formula
        d['formula_discharge'] = self.fully_discharged_entry.composition.reduced_formula
        d['fracA_charge'] = self.fully_charged_entry.composition.get_atomic_fraction(self.working_ion)
        d['fracA_discharge'] = self.fully_discharged_entry.composition.get_atomic_fraction(self.working_ion)
        d['max_instability'] = self.get_max_instability()
        d['min_instability'] = self.get_min_instability()
        if print_subelectrodes:
            f_dict = lambda c: c.to_dict_summary(print_subelectrodes=False)
            d['adj_pairs'] = map(f_dict, self.get_sub_electrodes(adjacent_only=True))
            d['all_pairs'] = map(f_dict, self.get_sub_electrodes(adjacent_only=False))
        return d

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return 'InsertionElectrode with endpoints at %s and %s, average voltage %f, capacity (grav.) %f, capacity (vol.) %f' % (self.fully_charged_entry.composition.reduced_formula, self.fully_discharged_entry.composition.reduced_formula, self.get_average_voltage(), self.get_capacity_grav(), self.get_capacity_vol())

    @staticmethod
    def from_dict(d):
        from pymatgen.serializers.json_coders import PMGJSONDecoder
        dec = PMGJSONDecoder()
        return InsertionElectrode(dec.process_decoded(d['entries']),
                                  dec.process_decoded(d['working_ion_entry']))

    @property
    def to_dict(self):
        d = {}
        d['module'] = self.__class__.__module__
        d['class'] = self.__class__.__name__
        d['entries'] = [entry.to_dict for entry in self._entries]
        d['working_ion_entry'] = self.working_ion_entry.to_dict
        return d


class InsertionVoltagePair(AbstractVoltagePair):
    '''
    Defines an Insertion Voltage Pair.
    '''

    def __init__(self, entry1, entry2, working_ion_entry):
        """
        Args:
            entry1:
                Entry corresponding to one of the entries in the voltage step.
            entry2:
                Entry corresponding to the other entry in the voltage step.
            working_ion_entry:
                A single ComputedEntry or PDEntry representing the element that
                carries charge across the battery, e.g. Li.
        """
        #initialize some internal variables
        working_element = working_ion_entry.composition.elements[0]
        frame1_dict = entry1.composition.to_dict
        if working_element.symbol in frame1_dict:
            del frame1_dict[working_element.symbol]
        frame1_comp = Composition.from_dict(frame1_dict)

        frame2_dict = entry2.composition.to_dict
        if working_element.symbol in frame2_dict:
            del frame2_dict[working_element.symbol]
        frame2_comp = Composition.from_dict(frame2_dict)

        #Data validation

        #check that the ion is just a single element
        if not working_ion_entry.composition.is_element:
            raise ValueError('VoltagePair: The working ion specified must be an element')

        #check that at least one of the entries contains the working element
        if not entry1.composition.get_atomic_fraction(working_element) > 0 and not entry2.composition.get_atomic_fraction(working_element) > 0:
            raise ValueError('VoltagePair: The working ion must be present in one of the entries')

        #check that the entries do not contain the same amount of the working element
        if entry1.composition.get_atomic_fraction(working_element) == entry2.composition.get_atomic_fraction(working_element):
            raise ValueError('VoltagePair: The working ion atomic percentage cannot be the same in both the entries')

        #check that the frameworks of the entries are equivalent other than the amount of working element
        if not frame1_comp.reduced_formula == frame2_comp.reduced_formula:
            raise ValueError('VoltagePair: the specified entries must have the same compositional framework')

        #Initialize normalization factors, charged and discharged entries
        entry1_normalization = frame1_comp.get_reduced_formula_and_factor()[1]
        entry2_normalization = frame2_comp.get_reduced_formula_and_factor()[1]

        entry1_is_charged = entry1.composition.get_atomic_fraction(working_element) < entry2.composition.get_atomic_fraction(working_element)

        entry_charge = entry1 if entry1_is_charged else entry2
        entry_discharge = entry2 if entry1_is_charged else entry1
        normalization_charge = entry1_normalization if entry1_is_charged else entry2_normalization
        normalization_discharge = entry2_normalization if entry1_is_charged else entry1_normalization

        self.framework = Composition.from_formula(frame2_comp.get_reduced_formula_and_factor()[0])
        self._working_ion_entry = working_ion_entry

        #Initialize normalized properties
        self._vol_charge = entry_charge.structure.volume / normalization_charge
        self._vol_discharge = entry_discharge.structure.volume / normalization_discharge

        comp_charge = entry_charge.composition
        comp_discharge = entry_discharge.composition

        self._mass_charge = comp_charge.weight / normalization_charge
        self._mass_discharge = comp_discharge.weight / normalization_discharge

        self._num_ions_transferred = (comp_discharge[working_element] / normalization_discharge) - (comp_charge[working_element] / normalization_charge)

        self._voltage = ((entry_charge.energy / normalization_charge) - (entry_discharge.energy / normalization_discharge)) / self._num_ions_transferred + working_ion_entry.energy_per_atom
        self._mAh = self._num_ions_transferred * ELECTRON_TO_AMPERE_HOURS * 1000

        #Step 4: add (optional) hull and muO2 data
        self.decomp_e_charge = entry_charge.data['decomposition_energy'] if 'decomposition_energy' in entry_charge.data else None
        self.decomp_e_discharge = entry_discharge.data['decomposition_energy'] if 'decomposition_energy' in entry_discharge.data else None

        self._muO2_charge = entry_charge.data['muO2'] if 'muO2' in entry_charge.data else None
        self._muO2_discharge = entry_discharge.data['muO2'] if 'muO2' in entry_discharge.data else None

        self.entry_charge = entry_charge
        self.entry_discharge = entry_discharge
        self.normalization_charge = normalization_charge
        self.normalization_discharge = normalization_discharge

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
