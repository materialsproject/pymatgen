#!/usr/bin/env python

'''
Created on Mar 19, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong, Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 19, 2012"


import os
import ConfigParser
from collections import defaultdict

from pymatgen.core.periodic_table import Element
from pymatgen.entries.post_processors_abc import EntryPostProcessor
from pymatgen.io.vaspio_set import MaterialsProjectVaspInputSet, MITVaspInputSet

class MaterialsProjectCompatibility(EntryPostProcessor):
    """
    This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. Note that this should only be used for VASP calculations using the
    MaterialsProject parameters (see pymatgen.io.vaspio_set 
    MaterialsProjectVaspInputSet). Using this compatibility scheme on runs with
    different parameters is not valid.
    """

    def __init__(self, compat_type="Advanced"):
        """
        Arguments:
            compat_type:
                Two options, GGA or Advanced.  GGA means all GGA+U entries are 
                excluded.  Advanced means mixing scheme is implemented to make 
                entries compatible with each other, but entries which are 
                supposed to be done in GGA+U will have the equivalent GGA entries 
                excluded. For example, Fe oxides should have a U value under the 
                Advanced scheme. A GGA Fe oxide run will therefore be excluded 
                under the scheme.
        """
        self.compat_type = compat_type
        module_dir = os.path.dirname(os.path.abspath(__file__))
        self._config = ConfigParser.SafeConfigParser()
        self._config.optionxform = str
        self._config.readfp(open(os.path.join(module_dir, "Compatibility.cfg")))
        u_corrections = dict(self._config.items('AdvancedUCorrections'))
        u_corrections_sulfides = dict(self._config.items('AdvancedUCorrectionsSulfides'))
        cpd_energies = dict(self._config.items('AdvancedCompoundEnergies'))

        self._u_corrections = dict()
        for key, val in u_corrections.items():
            self._u_corrections[Element(key)] = float(val)
        self._u_corrections_sulfides = dict()
        for key, val in u_corrections_sulfides.items():
            self._u_corrections_sulfides[Element(key)] = float(val)
        self._cpd_energies = dict()
        for key, val in cpd_energies.items():
            self._cpd_energies[key] = float(val)

        input_set = MaterialsProjectVaspInputSet()
        self._valid_potcars = set(input_set.potcar_settings.values())
        self._oxide_u = {Element(k): v for k, v in input_set.incar_settings["LDAUU"].items()}

        if compat_type == "GGA":
            self._u_corrections = dict()
            self._u_corrections_sulfides = dict()
            self._oxide_u = defaultdict(int)


    def has_u_element_oxides(self, comp):
        if Element("O") not in comp:
            return False
        for el in comp.elements:
            if el in self._oxide_u:
                return True
        return False

    def process_entry(self, entry):
        """
        Process a single entry with the chosen Compatibility scheme.
        
        Args:
            entry - An ComputedEntry object.
        
        Returns:
            An adjusted entry if entry is compatible, otherwise None is returned.
        """
        if entry.parameters.get('run_type', 'GGA') == "HF":
            return None

        ucorr = self._u_corrections
        cpdenergies = self._cpd_energies
        u_settings = entry.parameters['hubbards']
        u_settings = defaultdict(int) if u_settings == None else u_settings
        comp = entry.composition
        #Check that POTCARs are valid
        rform = comp.reduced_formula
        if rform not in cpdenergies:
            psp_settings = set([sym.split(" ")[1] for sym in entry.parameters['potcar_symbols']])
            if not self._valid_potcars.issuperset(psp_settings):

                #print "invalid psp"
                return None


        if comp.is_element:
            #correct all elements that are wrong, e.g. O2 molecule
            if rform in cpdenergies:
                entry.structureid = -comp.keys()[0].Z
                entry.correction = cpdenergies[rform] * comp.num_atoms - entry.uncorrected_energy
            return entry
        elif self.has_u_element_oxides(comp):
            correction = 0
            for el in comp.elements:
                if el in ucorr:
                    if el.symbol in u_settings and u_settings[el.symbol] == self._oxide_u[el]:
                        correction += float(ucorr[el]) * comp[el]
                    else:
                        return None
            entry.correction = correction
            return entry

        return entry

    def process_entries(self, entries):
        """
        Process a sequence of entries with the chosen Compatibility scheme.
        
        Args:
            entries - A sequence of entries.
        
        Returns:
            An list of adjusted entries.  Entries in the original list which 
            are not compatible are excluded.
        """
        proc_entries = list()
        for entry in entries:
            proc_entry = self.process_entry(entry)
            if proc_entry != None:
                proc_entries.append(proc_entry)
        return proc_entries

    @property
    def corrected_compound_formulas(self):
        return self._cpd_energies.keys()

    def __str__(self):
        return self.compat_type + " Compatibility corrects element states and mixes GGA/GGA+U calcs"



class MITCompatibility(MaterialsProjectCompatibility):
    """
    This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. Note that this should only be used for VASP calculations using the
    MIT parameters (see pymatgen.io.vaspio_set MITVaspInputSet). Using 
    this compatibility scheme on runs with different parameters is not valid.
    """

    def __init__(self, compat_type="Advanced"):
        """
        Arguments:
            compat_type:
                Two options, GGA or Advanced.  GGA means all GGA+U entries are 
                excluded.  Advanced means mixing scheme is implemented to make 
                entries compatible with each other, but entries which are 
                supposed to be done in GGA+U will have the equivalent GGA entries 
                excluded. For example, Fe oxides should have a U value under the 
                Advanced scheme. A GGA Fe oxide run will therefore be excluded 
                under the scheme.
        """
        self.compat_type = compat_type
        module_dir = os.path.dirname(os.path.abspath(__file__))
        self._config = ConfigParser.SafeConfigParser()
        self._config.optionxform = str
        self._config.readfp(open(os.path.join(module_dir, "MITCompatibility.cfg")))
        u_corrections = dict(self._config.items('AdvancedUCorrections'))
        u_corrections_sulfides = dict(self._config.items('AdvancedUCorrectionsSulfides'))
        cpd_energies = dict(self._config.items('AdvancedCompoundEnergies'))

        self._u_corrections = dict()
        for key, val in u_corrections.items():
            self._u_corrections[Element(key)] = float(val)
        self._u_corrections_sulfides = dict()
        for key, val in u_corrections_sulfides.items():
            self._u_corrections_sulfides[Element(key)] = float(val)
        self._cpd_energies = dict()
        for key, val in cpd_energies.items():
            self._cpd_energies[key] = float(val)

        input_set = MITVaspInputSet()
        self._valid_potcars = set(input_set.potcar_settings.values())
        self._oxide_u = {Element(k): v for k, v in input_set.incar_settings["LDAUU"].items()}

        if compat_type == "GGA":
            self._u_corrections = dict()
            self._u_corrections_sulfides = dict()
            self._oxide_u = defaultdict(int)

    def has_u_element_oxides(self, comp):
        if Element("O") not in comp and Element("F") not in comp:
            return False
        for el in comp.elements:
            if el in self._oxide_u:
                return True
        return False
