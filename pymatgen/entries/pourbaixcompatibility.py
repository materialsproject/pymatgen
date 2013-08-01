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
from pymatgen.util.coord_utils import pbc_all_distances
from pymatgen.core.periodic_table import Element
import numpy as np

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
            try:
                ox_type = entry.data["oxide_type"]
                if ox_type == "hydroxide":
                    corr_str = "Hydroxide"
                    n_ion = min(comp["O"], comp["H"])
                elif ox_type == "superoxide":
                    corr_str = "Superoxide"
                    n_ion = (comp["O"]) / 2.0
                elif ox_type == "peroxide":
                    corr_str = "Peroxide"
                    n_ion = (comp["O"]) / 2.0
                elif ox_type == "ozonide":
                    return None
                else:
                    return entry
            except:
                ox_type = oxide_type(entry.structure, 1.1)
                if ox_type == "ozonide":
                    return None
                elif ox_type == "superoxide":
                    corr_str = "Superoxide"
                    n_ion = (comp["O"]) / 2.0
                elif ox_type == "peroxide":
                    corr_str = "Peroxide"
                    n_ion = (comp["O"]) / 2.0
                elif ox_type == "hydroxide":
                    corr_str = "Hydroxide"
                    n_ion = min(comp["O"], comp["H"])
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


def oxide_type(structure, relative_cutoff=1.2):
    """
    Determines if an oxide is a peroxide/superoxide/ozonide/normal oxide
    
    Args:
        structure:
            Input structure.
        relative_cutoff:
            Relative_cutoff * act. cutoff stipulates the max. distance two 
            O atoms must be from each other. 
    """

    o_sites_frac_coords = []
    h_sites_frac_coords = []
    lattice = structure.lattice

#    Check if H2O/O2/H2O2/non-oxide is being checked
    el_other_than_OH = [el.symbol for el in structure.composition.elements if el.symbol not in ["H", "O"]]
    if (not el_other_than_OH) | (len(el_other_than_OH) == len(structure.composition.elements)):
        return "None"

    for site in structure:
        syms = [sp. symbol for sp in site.species_and_occu.keys()]
        if "O" in syms:
            o_sites_frac_coords.append(site.frac_coords)
        if "H" in syms:
            h_sites_frac_coords.append(site.frac_coords)

#   We are not correcting hydrides
    if not o_sites_frac_coords:
        return "None"

    if h_sites_frac_coords:
        dist_matrix = pbc_all_distances(lattice, o_sites_frac_coords, h_sites_frac_coords)
        if np.any(dist_matrix < relative_cutoff * 0.93):
            return "hydroxide"
    dist_matrix = pbc_all_distances(lattice, o_sites_frac_coords, o_sites_frac_coords)
    np.fill_diagonal(dist_matrix, 1000)
    is_superoxide = False
    is_peroxide = False
    is_ozonide = False
    if np.any(dist_matrix < relative_cutoff * 1.35):
        is_superoxide = True
    elif np.any(dist_matrix < relative_cutoff * 1.49):
        is_peroxide = True
    if is_superoxide:
        if (int(structure.composition["O"] / structure.composition.get_reduced_composition_and_factor()[1]) % 3 == 0):
            sum = 0
            for elt in [el for el in structure.composition.elements if el is not Element("O")]:
                sum += structure.composition[elt.symbol]
            if sum / structure.composition.get_reduced_composition_and_factor()[1] <= 1:
                is_superoxide = False
                is_ozonide = True
    if is_ozonide: 
        return "ozonide"
    elif is_superoxide:
        return "superoxide"
    elif is_peroxide:
        return "peroxide"
    else:
        return "oxide"


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
