#!/usr/bin/env python

"""
This module implements Compatibility corrections for mixing runs of different
functionals.
"""

from __future__ import division

__author__ = "Shyue Ping Ong, Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 19, 2012"


import os
import ConfigParser
from collections import defaultdict

from pymatgen.core.composition import Composition
from pymatgen.entries.post_processors_abc import EntryPostProcessor
from pymatgen.io.vaspio_set import MITVaspInputSet, MPVaspInputSet


class Compatibility(EntryPostProcessor):
    """
    This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. This is a base class from which other specific compatibility
    schemes are implemented.

    For compatibility to be checked, the entry supplied have two additional
    restrictions in terms of its parameters key:

    1. Entry.parameters must contain a "hubbards" key which is a dict of all
       non-zero Hubbard U values used in the calculation. For example,
       if you ran a Fe2O3 calculation with Materials Project parameters,
       this would look like entry.parameters["hubbards"] = {"Fe": 5.3}
       If the "hubbards" key is missing, a GGA run is assumed.
    2. Entry.parameters must contain a "potcar_symbols" key that is a list of
       all POTCARs used in the run. Again, using the example of an Fe2O3 run
       using Materials Project parameters, this would look like
       entry.parameters["potcar_symbols"] = ['PAW_PBE Fe_pv 06Sep2000',
       'PAW_PBE O 08Apr2002'].

    It should be noted that ComputedEntries assimilated using the
    pymatgen.apps.borg package and obtained via the MaterialsProject REST
    interface using the pymatgen.matproj.rest package will automatically have
    these fields populated.
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
        if input_set_name == "MaterialsProject":
            self.input_set = MPVaspInputSet()
        elif input_set_name == "MITMatgen":
            self.input_set = MITVaspInputSet()
        else:
            raise ValueError("Invalid input set name {}".format(input_set_name))

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

        cpd_energies = dict(self._config.items(
            "{}{}CompoundEnergies".format(input_set_name, compat_type)))

        self.u_corrections = u_corrections
        self.cpd_energies = {k: float(v) for k, v in cpd_energies.items()}

        self.valid_potcars = set(self.input_set.potcar_settings.values())
        self.u_settings = self.input_set.incar_settings["LDAUU"]

        if compat_type == "GGA":
            self.u_corrections = {}
            self.u_settings = {}

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

        usettings = self.u_settings.get(most_electroneg, {})

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

        Raises:
            ValueError if entry do not contain "potcar_symbols" key.
        """
        if entry.parameters.get("run_type", "GGA") == "HF":
            return None

        cpdenergies = self.cpd_energies
        calc_u = entry.parameters.get("hubbards", None)
        calc_u = defaultdict(int) if calc_u is None else calc_u
        comp = entry.composition
        #Check that POTCARs are valid
        rform = comp.reduced_formula
        if rform not in cpdenergies:
            try:
                psp_settings = set([sym.split(" ")[1]
                                    for sym
                                    in entry.parameters["potcar_symbols"]])
            except KeyError:
                raise ValueError("Compatibility can only be checked for "
                                 "entries with a \"potcar_symbols\" in "
                                 "entry.parameters")
            if not self.valid_potcars.issuperset(psp_settings):
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

            ucorr = self.u_corrections.get(most_electroneg, {})
            usettings = self.u_settings.get(most_electroneg, {})

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
        return filter(None, map(self.process_entry, entries))

    @property
    def corrected_compound_formulas(self):
        return self.cpd_energies.keys()

    def __str__(self):
        return "{} {} Compatibility".format(self.input_set_name,
                                            self.compat_type)


class MaterialsProjectCompatibility(Compatibility):
    """
    This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. Note that this should only be used for VASP calculations using the
    MaterialsProject parameters (see pymatgen.io.vaspio_set.MPVaspInputSet).
    Using this compatibility scheme on runs with different parameters is not
    valid.
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
