#!/usr/bin/env python
"""
This module implements corrections for H+, and ozonides
"""

__author__ = "Sai Jayaraman" 
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Sai Jayaraman"
__email__ = "sjayaram@mit.edu"
__date__ = "Mar 5, 2013"

from pymatgen.entries.post_processors_abc import EntryPostProcessor

import os
import ConfigParser


class AqueousCompatibility(EntryPostProcessor):
    """
    This class implements aqueous phase compound corrections for elements, and H2O.
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

        module_dir = os.path.dirname(os.path.abspath(__file__))
        self._config = ConfigParser.SafeConfigParser()
        self._config.optionxform = str
        self._config.readfp(open(os.path.join(module_dir,
                                              "AqueousCompatibility.cfg")))

        cpd_energies = dict(self._config.items(
            "{}{}CompoundEnergies".format(input_set_name, compat_type)))

        self.cpd_energies = {k: float(v) for k, v in cpd_energies.items()}

    @property
    def corrected_compound_formulas(self):
        return

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
        comp = entry.composition
        rform = comp.reduced_formula
        cpdenergies = self.cpd_energies
        if rform in cpdenergies:
            entry.structureid = -comp.keys()[0].Z
            entry.correction = cpdenergies[rform] * comp.num_atoms \
                - entry.uncorrected_energy
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
        return filter(None, map(self.process_entry, entries))


class MaterialsProjectAqueousCompatibility(AqueousCompatibility):
    """
    This class implements aqueous phase compound corrections for elements, and H2O.
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
        AqueousCompatibility.__init__(self, "MaterialsProject", compat_type)


class MITAqueousCompatibility(MaterialsProjectAqueousCompatibility):
    """
    This class implements aqueous phase compound corrections for elements, and H2O.
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
        AqueousCompatibility.__init__(self, "MITMatgen", compat_type)

