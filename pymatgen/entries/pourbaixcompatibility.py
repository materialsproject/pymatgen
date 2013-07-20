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

import itertools
from pymatgen.entries.post_processors_abc import EntryPostProcessor
from pymatgen.analysis.structure_analyzer import contains_peroxide
import os
import ConfigParser


class PourbaixCompatibility(EntryPostProcessor):
    """
    This class implements corrections for peroxides, superoxides, ozonides, and O-H bonds
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
                                              "PourbaixCompatibility.cfg")))

        self.corr_str = str(input_set_name) + str(compat_type) +\
         str("U") + "Corrections"
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
        n_ion = 1
        cpdenergies = self.cpd_energies
        if rform in cpdenergies:
            entry.structureid = -comp.keys()[0].Z
            entry.correction = cpdenergies[rform] * comp.num_atoms \
                - entry.uncorrected_energy
        else:
            if contains_hydroxide(entry.structure):
                corr_str = "Hydroxide"
                n_ion = min(comp["O"], comp["H"])
            elif contains_ozonide(entry.structure, 1.15):
                return None
            elif contains_superoxide(entry.structure, 1.15):
                corr_str = "Superoxide"
                n_ion = (comp["O"]) / 2
            elif contains_peroxide(entry.structure, 1.15):
                corr_str = "Peroxide"
                n_ion = (comp["O"]) / 2
            else:
                return entry
            corr_str = self.corr_str + corr_str
            corr = self._config.items(corr_str)
            entry.correction += float(corr[0][1]) * float(n_ion)
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


def contains_ozonide(structure, relative_cutoff=1.2):
    """
    Determines if a structure contains ozonide anions.

    Args:
        structure:
            Input structure.
        relative_cutoff:
            The ozonide bond distance is 1.35 Angstrom. Relative_cutoff * 1.35
            stipulates the maximum distance two O atoms must be to each other
            to be considered a ozonide. The angle should be 108-110 degrees
            (can we measure this?)

    Returns:
        Boolean indicating if structure contains a ozonide anion.
    """
    max_dist = relative_cutoff * 1.35
    o_sites = []
    for site in structure:
        syms = [sp.symbol for sp in site.species_and_occu.keys()]
        if "O" in syms:
            o_sites.append(site)

    is_ozonide = False
    for i, j in itertools.combinations(o_sites, 2):
        if i.distance(j) < max_dist:
            is_ozonide = True
            break
    from pymatgen.core.periodic_table import Element
    if is_ozonide:
        if (int(structure.composition["O"] / structure.composition.get_reduced_composition_and_factor()[1]) % 3 == 0):
            sum = 0
            for elt in [el for el in structure.composition.elements if el is not Element("O")]:
                sum += structure.composition[elt.symbol]
            if sum / structure.composition.get_reduced_composition_and_factor()[1] > 1:
                return False
            else:
                return True
    return False


def contains_superoxide(structure, relative_cutoff=1.2):
    """
    Determines if a structure contains superoxide anions

    Args:
        structure:
            Input structure.
        relative_cutoff:
            The superoxide bond distance is 1.35 Angstrom. Relative cutoff * 1.35
            stipulates the maximum distance two O atoms must be to each other
            to be considered a superoxide.
    """
    max_dist = relative_cutoff * 1.35
    o_sites = []
    for site in structure:
        syms = [sp.symbol for sp in site.species_and_occu.keys()]
        if "O" in syms:
            o_sites.append(site)
    is_superoxide = False
    for i, j in itertools.combinations(o_sites, 2):
        if i.distance(j) < max_dist and i.distance(j) < 1.49:
            is_superoxide = True
            break
    if is_superoxide:
        if not(contains_ozonide(structure, relative_cutoff)):
            return True
    return False


def contains_hydroxide(structure, relative_cutoff=1.2):
    """
    Determines if a structure contains hydroxide anions

    Args:
        structure: 
            Input structure
        relative_cutoff: 
            The O-H bond distance is ~0.93 Angstrom. Relative cutoff * 0.93 
            stipulates the maximum distance an O and H must be from each other
            to be considered a hydroxide.
    """
    max_dist = relative_cutoff * 0.93
    o_sites = []
    h_sites = []
    for site in structure:
        syms = [sp.symbol for sp in site.species_and_occu.keys()]
        if "O" in syms:
            o_sites.append(site)
        if "H" in syms:
            h_sites.append(site)
    maybe_oh = False
    for i in o_sites:
        for j in h_sites:
            if i.distance(j) < max_dist:
                maybe_oh = True
    if maybe_oh:
        if len([elt.symbol for elt in structure.composition.keys() if elt.symbol not in ["O", "H"]]) == 0.0 :
            maybe_oh = False

    return maybe_oh


class MaterialsProjectPourbaixCompatibility(PourbaixCompatibility):
    """
    This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. Note that this should only be used for VASP calculations using the
    MaterialsProject parameters (see pymatgen.io.vaspio_set
    MaterialsProjectVaspInputSet). Using this compatibility scheme on runs with
    different parameters is not valid.
    """

    def __init__(self, compat_type="GGA"):
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
        PourbaixCompatibility.__init__(self, "MaterialsProject", compat_type)


class MITPourbaixCompatibility(MaterialsProjectPourbaixCompatibility):
    """
    This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. Note that this should only be used for VASP calculations using the
    MIT parameters (see pymatgen.io.vaspio_set MITVaspInputSet). Using
    this compatibility scheme on runs with different parameters is not valid.
    """

    def __init__(self, compat_type="GGA"):
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
        PourbaixCompatibility.__init__(self, "MITMatgen", compat_type)
