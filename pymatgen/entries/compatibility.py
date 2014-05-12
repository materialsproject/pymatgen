#!/usr/bin/env python

"""
This module implements Compatibility corrections for mixing runs of different
functionals.
"""

from __future__ import division

__author__ = "Shyue Ping Ong, Anubhav Jain, Sai Jayaraman"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 19, 2012"


import os
import ConfigParser
from collections import defaultdict

from pymatgen.io.vaspio_set import MITVaspInputSet, MPVaspInputSet
from pymatgen.core.periodic_table import Element
from pymatgen.analysis.structure_analyzer import oxide_type

import abc


class CompatibilityError(Exception):
    """
    Exception class for Compatibility. Raised by attempting correction
    on incompatible calculation
    """
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class Correction(object):
    """
    A Correction class is a pre-defined scheme for correction a computed
    entry based on the type and chemistry of the structure and the
    calculation parameters. All Correction classes must implement a
    correct_entry method.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def get_correction(self, entry):
        """
        Returns correction for a single entry.

        Args:
            entry: A ComputedEntry object.

        Returns:
            The energy correction to be applied.

        Raises:
            CompatibilityError if entry is not compatible.
        """
        return

    def correct_entry(self, entry):
        """
        Corrects a single entry.

        Args:
            entry: A ComputedEntry object.

        Returns:
            An processed entry.

        Raises:
            CompatibilityError if entry is not compatible.
        """
        entry.correction += self.get_correction(entry)
        return entry


class PotcarCorrection(Correction):
    """
    Checks that POTCARs are valid within a pre-defined input set. This
    ensures that calculations performed using different InputSets are not
    compared against each other.

    Entry.parameters must contain a "potcar_symbols" key that is a list of
    all POTCARs used in the run. Again, using the example of an Fe2O3 run
    using Materials Project parameters, this would look like
    entry.parameters["potcar_symbols"] = ['PAW_PBE Fe_pv 06Sep2000',
    'PAW_PBE O 08Apr2002'].

    Args:
        name: Name of settings to use. Current valid settings are MP or
            MIT, which is the relevant settings based on the MP or MIT
            VaspInputSets.

    Raises:
        ValueError if entry do not contain "potcar_symbols" key.
        CombatibilityError if wrong potcar symbols
    """
    def __init__(self, name):
        if name == "MP":
            input_set = MPVaspInputSet()
        elif name == "MIT":
            input_set = MITVaspInputSet()
        else:
            raise ValueError("Only MIT and MP POTCAR corrections are "
                             "supported currently.")
        self.valid_potcars = set(input_set.potcar_settings.values())
        self.name = name

    def get_correction(self, entry):
        try:
            psp_settings = set([sym.split(" ")[1]
                                for sym
                                in entry.parameters["potcar_symbols"]])
        except KeyError:
            raise ValueError(
                "PotcarCorrection can only be checked for entries with a "
                "\"potcar_symbols\" in entry.parameters")
        if not self.valid_potcars.issuperset(psp_settings):
            raise CompatibilityError('Incompatible potcar')
        return 0

    def __str__(self):
        return "{} Potcar Correction".format(self.name)


class GasCorrection(Correction):
    """
    Correct gas energies to obtain the right formation energies. Note that
    this depends on calculations being run within the same input set.

    Args:
        name: Name of settings to use. Current valid settings are MP or
            MIT, which is the relevant settings based on the MP or MIT
            VaspInputSets.
        correct_peroxide: Specify whether peroxide/superoxide/ozonide 
            corrections are to be applied or not. 
    """
    def __init__(self, name, correct_peroxide=True):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        config = ConfigParser.SafeConfigParser()
        config.optionxform = str
        config.readfp(open(os.path.join(module_dir, "Compatibility.cfg")))
        cpd_energies = dict(
            config.items("{}CompoundEnergies".format(name)))
        self.cpd_energies = {k: float(v) for k, v in cpd_energies.items()}
        self.oxide_correction = {
            k: float(v) for k, v
            in config.items("{}OxideCorrection".format(name))}
        self.name = name
        self.correct_peroxide = correct_peroxide

    def get_correction(self, entry):
        comp = entry.composition

        rform = entry.composition.reduced_formula
        if rform in self.cpd_energies:
            return self.cpd_energies[rform] * comp.num_atoms \
                - entry.uncorrected_energy

        correction = 0
        #Check for oxide, peroxide, superoxide, and ozonide corrections.
        if self.correct_peroxide:
            if len(comp) >= 2 and Element("O") in comp:
                if "oxide_type" in entry.data:
                    if entry.data["oxide_type"] in self.oxide_correction:
                        ox_corr = self.oxide_correction[
                            entry.data["oxide_type"]]
                        correction += ox_corr * comp["O"]
                    if entry.data["oxide_type"] == "hydroxide":
                        ox_corr = self.oxide_correction["oxide"]
                        correction += ox_corr * comp["O"]
    
                elif hasattr(entry, "structure"):
                    ox_type, nbonds = oxide_type(entry.structure, 1.05,
                                                 return_nbonds=True)
                    if ox_type in self.oxide_correction:
                        correction += self.oxide_correction[ox_type] * \
                            nbonds
                    elif ox_type == "hydroxide":
                        correction += self.oxide_correction["oxide"] * comp["O"]
                else:
                    if rform in UCorrection.common_peroxides:
                        correction += self.oxide_correction["peroxide"] * \
                            comp["O"]
                    elif rform in UCorrection.common_superoxides:
                        correction += self.oxide_correction["superoxide"] * \
                            comp["O"]
                    elif rform in UCorrection.ozonides:
                        correction += self.oxide_correction["ozonide"] * \
                            comp["O"]
                    elif Element("O") in comp.elements and len(comp.elements) > 1:
                        correction += self.oxide_correction['oxide'] * comp["O"]
        else:
            correction += self.oxide_correction['oxide'] * comp["O"]

        return correction

    def __str__(self):
        return "{} Gas Correction".format(self.name)


class AqueousCorrection(Correction):
    """
    This class implements aqueous phase compound corrections for elements
    and H2O.

    Args:
        name: The name of the input set to use. Can be either MP or MIT.
    """
    def __init__(self, name):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        config = ConfigParser.SafeConfigParser()
        config.optionxform = str
        config.readfp(open(os.path.join(module_dir,
                                        "Compatibility.cfg")))
        cpd_energies = dict(
            config.items("{}AqueousCompoundEnergies".format(name)))
        self.cpd_energies = {k: float(v) for k, v in cpd_energies.items()}
        self.name = name

    def get_correction(self, entry):
        comp = entry.composition
        rform = comp.reduced_formula
        cpdenergies = self.cpd_energies
        correction = 0
        if rform in cpdenergies:
            if rform in ["H2", "H2O"]:
                correction = cpdenergies[rform] * comp.num_atoms \
                    - entry.uncorrected_energy - entry.correction
            else:
                correction += cpdenergies[rform] * comp.num_atoms
        if not rform == "H2O":
            correction += 0.5 * 2.46 * min(comp["H"]/2.0, comp["O"])
        return correction

    def __str__(self):
        return "{} Aqueous Correction".format(self.name)


class UCorrection(Correction):
    """
    This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. Entry.parameters must contain a "hubbards" key which is a dict
    of all non-zero Hubbard U values used in the calculation. For example,
    if you ran a Fe2O3 calculation with Materials Project parameters,
    this would look like entry.parameters["hubbards"] = {"Fe": 5.3}
    If the "hubbards" key is missing, a GGA run is assumed.

    It should be noted that ComputedEntries assimilated using the
    pymatgen.apps.borg package and obtained via the MaterialsProject REST
    interface using the pymatgen.matproj.rest package will automatically have
    these fields populated.

    Args:
        name: Name of settings to use. Current valid settings are MP or
            MIT, which is the relevant settings based on the MP or MIT
            VaspInputSets.
        compat_type: Two options, GGA or Advanced.  GGA means all GGA+U
            entries are excluded.  Advanced means mixing scheme is
            implemented to make entries compatible with each other,
            but entries which are supposed to be done in GGA+U will have the
            equivalent GGA entries excluded. For example, Fe oxides should
            have a U value under the Advanced scheme. A GGA Fe oxide run
            will therefore be excluded under the scheme.
    """
    common_peroxides = ["Li2O2", "Na2O2", "K2O2", "Cs2O2", "Rb2O2", "BeO2",
                        "MgO2", "CaO2", "SrO2", "BaO2"]
    common_superoxides = ["LiO2", "NaO2", "KO2", "RbO2", "CsO2"]
    ozonides = ["LiO3", "NaO3", "KO3", "NaO5"]

    def __init__(self, name, compat_type):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        config = ConfigParser.SafeConfigParser()
        config.optionxform = str
        config.readfp(open(os.path.join(module_dir, "Compatibility.cfg")))
        if name == "MP":
            self.input_set = MPVaspInputSet()
        elif name == "MIT":
            self.input_set = MITVaspInputSet()
        else:
            raise ValueError("Invalid input set name {}".format(name))

        u_corrections = {}
        for el in self.input_set.incar_settings["LDAUU"].keys():
            sect_name = "{}{}UCorrections{}".format(name, compat_type, el)
            if sect_name in config.sections():
                corr = dict(config.items(sect_name))
                u_corrections[el] = {k: float(v) for k, v in corr.items()}

        self.u_corrections = u_corrections
        self.u_settings = self.input_set.incar_settings["LDAUU"]

        if compat_type == "GGA":
            self.u_corrections = {}
            self.u_settings = {}

        self.name = name
        self.compat_type = compat_type

    def get_correction(self, entry):
        if entry.parameters.get("run_type", "GGA") == "HF":
            raise CompatibilityError('Invalid run type')

        calc_u = entry.parameters.get("hubbards", None)
        calc_u = defaultdict(int) if calc_u is None else calc_u
        comp = entry.composition

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
                raise CompatibilityError('Invalid U value on {}'.format(sym))
            if sym in ucorr:
                correction += float(ucorr[sym]) * comp[el]

        return correction

    def __str__(self):
        return "{} {} Correction".format(self.name, self.compat_type)


class Compatibility(object):
    """
    The Compatibility class combines a list of corrections to be applied to
    an entry or a set of entries. Note that some of the Corrections have
    interdependencies. For example, PotcarCorrection must always be used
    before any other compatibility. Also, GasCorrection("MP") must be used
    with PotcarCorrection("MP") (similarly with "MIT"). Typically,
    you should use the specific MaterialsProjectCompatibility and
    MITCompatibility subclasses instead.

    Args:
        corrections: List of corrections to apply.
    """
    def __init__(self, corrections):
        self.corrections = corrections

    def process_entry(self, entry):
        """
        Process a single entry with the chosen Corrections.

        Args:
            entry: A ComputedEntry object.

        Returns:
            An adjusted entry if entry is compatible, otherwise None is
            returned.
        """
        try:
            corrections = self.get_corrections_dict(entry)
        except CompatibilityError:
            return None
        entry.correction = sum(corrections.values())
        return entry

    def get_corrections_dict(self, entry):
        """
        Returns the corrections applied to a particular entry.

        Args:
            entry: A ComputedEntry object.

        Returns:
            ({correction_name: value})
        """
        corrections = {}
        for c in self.corrections:
            val = c.get_correction(entry)
            if val != 0:
                corrections[str(c)] = val
        return corrections

    def process_entries(self, entries):
        """
        Process a sequence of entries with the chosen Compatibility scheme.

        Args:
            entries: A sequence of entries.

        Returns:
            An list of adjusted entries.  Entries in the original list which
            are not compatible are excluded.
        """
        return filter(None, map(self.process_entry, entries))


class MaterialsProjectCompatibility(Compatibility):
    """
    This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. Note that this should only be used for VASP calculations using the
    MaterialsProject parameters (see pymatgen.io.vaspio_set.MPVaspInputSet).
    Using this compatibility scheme on runs with different parameters is not
    valid.

    Args:
        compat_type: Two options, GGA or Advanced.  GGA means all GGA+U
            entries are excluded.  Advanced means mixing scheme is
            implemented to make entries compatible with each other,
            but entries which are supposed to be done in GGA+U will have the
            equivalent GGA entries excluded. For example, Fe oxides should
            have a U value under the Advanced scheme. A GGA Fe oxide run
            will therefore be excluded under the scheme.
        correct_peroxide: Specify whether peroxide/superoxide/ozonide 
            corrections are to be applied or not. 
    """

    def __init__(self, compat_type="Advanced", correct_peroxide=True):
        name = "MP"
        Compatibility.__init__(
            self, [PotcarCorrection(name), GasCorrection(name, correct_peroxide=correct_peroxide),
                   UCorrection(name, compat_type)])


class MITCompatibility(Compatibility):
    """
    This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. Note that this should only be used for VASP calculations using the
    MIT parameters (see pymatgen.io.vaspio_set MITVaspInputSet). Using
    this compatibility scheme on runs with different parameters is not valid.

    Args:
        compat_type: Two options, GGA or Advanced.  GGA means all GGA+U
            entries are excluded.  Advanced means mixing scheme is
            implemented to make entries compatible with each other,
            but entries which are supposed to be done in GGA+U will have the
            equivalent GGA entries excluded. For example, Fe oxides should
            have a U value under the Advanced scheme. A GGA Fe oxide run
            will therefore be excluded under the scheme.
        correct_peroxide: Specify whether peroxide/superoxide/ozonide 
            corrections are to be applied or not. 
    """

    def __init__(self, compat_type="Advanced", correct_peroxide=True):
        name = "MIT"
        Compatibility.__init__(
            self, [PotcarCorrection(name), GasCorrection(name, correct_peroxide=correct_peroxide),
                   UCorrection(name, compat_type)])
