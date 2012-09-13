#!/usr/bin/env python

"""
This module implements Compatibility corrections for mixing runs of different
functionals.
"""

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

from pymatgen.core.structure import Composition
from pymatgen.entries.post_processors_abc import EntryPostProcessor
from pymatgen.io.vaspio_set import VaspInputSet


class Compatibility(EntryPostProcessor):
    """
    This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. This is a base class from which other specific compatibility
    schemes are implemented.
    """

    def __init__(self, input_set_name, compat_type):
        """
        Args:
            input_set_name:
                The name of the input set to use. Can be either
                MaterialsProject or MITMatgen.
            compat_type:
                Two options, GGA or Advanced.  GGA means all GGA+U entries are
                excluded.  Advanced means mixing scheme is implemented to make
                entries compatible with each other, but entries which are
                supposed to be done in GGA+U will have the equivalent GGA
                entries excluded. For example, Fe oxides should have a U value
                under the Advanced scheme. A GGA Fe oxide run will therefore be
                excluded under the scheme.
        """
        self.compat_type = compat_type
        self.input_set_name = input_set_name
        self.input_set = VaspInputSet(input_set_name)

        module_dir = os.path.dirname(os.path.abspath(__file__))
        self._config = ConfigParser.SafeConfigParser()
        self._config.optionxform = str
        self._config.readfp(open(os.path.join(module_dir,
                                              "Compatibility.cfg")))
        u_corrections = {}
        for el in self.input_set.incar_settings["LDAUU"].keys():
            name = "{}{}UCorrections{}".format(input_set_name, compat_type, el)
            if name in self._config.sections():
                corr = dict(self._config.items(name))
                u_corrections[el] = {k: float(v) for k, v in corr.items()}

        cpd_energies = dict(self._config.items("{}{}CompoundEnergies"
                                               .format(input_set_name,
                                                       compat_type)))

        self._u_corrections = u_corrections
        self._cpd_energies = {k: float(v) for k, v in cpd_energies.items()}

        self._valid_potcars = set(self.input_set.potcar_settings.values())
        self._u_settings = self.input_set.incar_settings["LDAUU"]

        if compat_type == "GGA":
            self._u_corrections = {}
            self._u_settings = {}

    def requires_hubbard(self, comp):
        """
        Check if a particular composition requies U parameters to be set.

        Args:
            comp:
                Composition

        Returns:
            True if hubbard U parameter required. False otherwise.
        """
        comp = Composition(comp)
        elements = sorted([el for el in comp.elements if comp[el] > 0],
                              key=lambda el: el.X)
        most_electroneg = elements[-1].symbol

        usettings = self._u_settings.get(most_electroneg, {})

        return any([usettings.get(el.symbol, 0) for el in comp.elements])

    def process_entry(self, entry):
        """
        Process a single entry with the chosen Compatibility scheme.

        Args:
            entry:
                A ComputedEntry object.

        Returns:
            An adjusted entry if entry is compatible, otherwise None is
            returned.
        """
        if entry.parameters.get("run_type", "GGA") == "HF":
            return None

        ucorr = self._u_corrections
        cpdenergies = self._cpd_energies
        calc_u = entry.parameters["hubbards"]
        calc_u = defaultdict(int) if calc_u == None else calc_u
        comp = entry.composition
        #Check that POTCARs are valid
        rform = comp.reduced_formula
        if rform not in cpdenergies:
            psp_settings = set([sym.split(" ")[1]
                                for sym in entry.parameters["potcar_symbols"]])
            if not self._valid_potcars.issuperset(psp_settings):
                return None

        #correct all compounds that are wrong, e.g. O2 molecule
        if rform in cpdenergies:
            entry.structureid = -comp.keys()[0].Z
            entry.correction = cpdenergies[rform] * comp.num_atoms \
                - entry.uncorrected_energy
        else:
            elements = sorted([el for el in comp.elements if comp[el] > 0],
                              key=lambda el: el.X)
            most_electroneg = elements[-1].symbol
            correction = 0

            ucorr = self._u_corrections.get(most_electroneg, {})
            usettings = self._u_settings.get(most_electroneg, {})

            for el in comp.elements:
                sym = el.symbol
                #Check for bad U values
                if calc_u.get(sym, 0) != usettings.get(sym, 0):
                    return None
                if sym in ucorr:
                    correction += float(ucorr[sym]) * comp[el]

            entry.correction = correction
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
        return "{} {} Compatibility".format(self.input_set_name,
                                            self.compat_type)


class MaterialsProjectCompatibility(Compatibility):
    """
    This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. Note that this should only be used for VASP calculations using the
    MaterialsProject parameters (see pymatgen.io.vaspio_set
    MaterialsProjectVaspInputSet). Using this compatibility scheme on runs with
    different parameters is not valid.
    """

    def __init__(self, compat_type="Advanced"):
        """
        Args:
            compat_type:
                Two options, GGA or Advanced.  GGA means all GGA+U entries are
                excluded.  Advanced means mixing scheme is implemented to make
                entries compatible with each other, but entries which are
                supposed to be done in GGA+U will have the equivalent GGA
                entries excluded. For example, Fe oxides should have a U value
                under the Advanced scheme. A GGA Fe oxide run will therefore be
                excluded under the scheme.
        """
        Compatibility.__init__(self, "MaterialsProject", compat_type)


class MITCompatibility(MaterialsProjectCompatibility):
    """
    This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. Note that this should only be used for VASP calculations using the
    MIT parameters (see pymatgen.io.vaspio_set MITVaspInputSet). Using
    this compatibility scheme on runs with different parameters is not valid.
    """

    def __init__(self, compat_type="Advanced"):
        """
        Args:
            compat_type:
                Two options, GGA or Advanced.  GGA means all GGA+U entries are
                excluded.  Advanced means mixing scheme is implemented to make
                entries compatible with each other, but entries which are
                supposed to be done in GGA+U will have the equivalent GGA
                entries excluded. For example, Fe oxides should have a U value
                under the Advanced scheme. A GGA Fe oxide run will therefore be
                excluded under the scheme.
        """
        Compatibility.__init__(self, "MITMatgen", compat_type)
