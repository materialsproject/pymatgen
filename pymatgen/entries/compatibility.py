# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
This module implements Compatibility corrections for mixing runs of different
functionals.
"""

import abc
import copy
import os
import warnings
from collections import defaultdict
from typing import Optional, Sequence, Union, List, Type

import numpy as np
from monty.design_patterns import cached_class
from monty.json import MSONable
from monty.serialization import loadfn
from uncertainties import ufloat

from pymatgen.analysis.structure_analyzer import oxide_type, sulfide_type
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.entries.entry_tools import EntrySet
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core.periodic_table import Element
from pymatgen.entries.computed_entries import (
    CompositionEnergyAdjustment,
    ComputedEntry,
    ComputedStructureEntry,
    ConstantEnergyAdjustment,
    TemperatureEnergyAdjustment,
)
from pymatgen.io.vasp.sets import MITRelaxSet, MPRelaxSet
from pymatgen.util.sequence import PBar

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
MU_H2O = -2.4583  # Free energy of formation of water, eV/H2O, used by MaterialsProjectAqueousCompatibility

__author__ = "Amanda Wang, Ryan Kingsbury, Shyue Ping Ong, Anubhav Jain, Stephen Dacek, Sai Jayaraman"
__copyright__ = "Copyright 2012-2020, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "April 2020"


class CompatibilityError(Exception):
    """
    Exception class for Compatibility. Raised by attempting correction
    on incompatible calculation
    """

    pass


class Correction(metaclass=abc.ABCMeta):
    """
    A Correction class is a pre-defined scheme for correction a computed
    entry based on the type and chemistry of the structure and the
    calculation parameters. All Correction classes must implement a
    correct_entry method.
    """

    @abc.abstractmethod
    def get_correction(self, entry):
        """
        Returns correction and uncertainty for a single entry.

        Args:
            entry: A ComputedEntry object.

        Returns:
            The energy correction to be applied and the uncertainty of the correction.

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
        new_corr = self.get_correction(entry)
        old_std_dev = entry.correction_uncertainty
        if np.isnan(old_std_dev):
            old_std_dev = 0
        old_corr = ufloat(entry.correction, old_std_dev)
        updated_corr = new_corr + old_corr

        if updated_corr.nominal_value != 0 and updated_corr.std_dev == 0:
            # if there are no error values available for the corrections applied,
            # set correction uncertainty to not a number
            uncertainty = np.nan
        else:
            uncertainty = updated_corr.std_dev

        entry.energy_adjustments.append(ConstantEnergyAdjustment(updated_corr.nominal_value, uncertainty))

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
    """

    def __init__(self, input_set, check_hash=False):
        """
        Args:
            input_set: InputSet object used to generate the runs (used to check
                for correct potcar symbols)

            check_hash (bool): If true, uses the potcar hash to check for valid
                potcars. If false, uses the potcar symbol (Less reliable).
                Defaults to True

        Raises:
            ValueError if entry do not contain "potcar_symbols" key.
            CombatibilityError if wrong potcar symbols
        """
        potcar_settings = input_set.CONFIG["POTCAR"]
        if isinstance(list(potcar_settings.values())[-1], dict):
            if check_hash:
                self.valid_potcars = {k: d["hash"] for k, d in potcar_settings.items()}
            else:
                self.valid_potcars = {k: d["symbol"] for k, d in potcar_settings.items()}
        else:
            if check_hash:
                raise ValueError("Cannot check hashes of potcars, since hashes are not included in the entry.")
            self.valid_potcars = potcar_settings

        self.input_set = input_set
        self.check_hash = check_hash

    def get_correction(self, entry) -> ufloat:
        """
        :param entry: A ComputedEntry/ComputedStructureEntry
        :return: Correction, Uncertainty.
        """
        if self.check_hash:
            if entry.parameters.get("potcar_spec"):
                psp_settings = {d.get("hash") for d in entry.parameters["potcar_spec"] if d}
            else:
                raise ValueError("Cannot check hash without potcar_spec field")
        else:
            if entry.parameters.get("potcar_spec"):
                psp_settings = {d.get("titel").split()[1] for d in entry.parameters["potcar_spec"] if d}
            else:
                psp_settings = {sym.split()[1] for sym in entry.parameters["potcar_symbols"] if sym}

        if {self.valid_potcars.get(str(el)) for el in entry.composition.elements} != psp_settings:
            raise CompatibilityError("Incompatible potcar")
        return ufloat(0.0, 0.0)

    def __str__(self):
        return "{} Potcar Correction".format(self.input_set.__name__)


@cached_class
class GasCorrection(Correction):
    """
    Correct gas energies to obtain the right formation energies. Note that
    this depends on calculations being run within the same input set.
    Used by legacy MaterialsProjectCompatibility and MITCompatibility.
    """

    def __init__(self, config_file):
        """
        Args:
            config_file: Path to the selected compatibility.yaml config file.
        """
        c = loadfn(config_file)
        self.name = c["Name"]
        self.cpd_energies = c["Advanced"]["CompoundEnergies"]

    def get_correction(self, entry) -> ufloat:
        """
        :param entry: A ComputedEntry/ComputedStructureEntry
        :return: Correction.
        """
        comp = entry.composition

        correction = ufloat(0.0, 0.0)

        # set error to 0 because old MPCompatibility doesn't have errors

        # only correct GGA or GGA+U entries
        if entry.parameters.get("run_type", None) not in ["GGA", "GGA+U"]:
            return ufloat(0.0, 0.0)

        rform = entry.composition.reduced_formula
        if rform in self.cpd_energies:
            correction += self.cpd_energies[rform] * comp.num_atoms - entry.uncorrected_energy

        return correction

    def __str__(self):
        return "{} Gas Correction".format(self.name)


@cached_class
class AnionCorrection(Correction):
    """
    Correct anion energies to obtain the right formation energies. Note that
    this depends on calculations being run within the same input set.

    Used by legacy MaterialsProjectCompatibility and MITCompatibility.
    """

    def __init__(self, config_file, correct_peroxide=True):
        """
        Args:
            config_file: Path to the selected compatibility.yaml config file.
            correct_peroxide: Specify whether peroxide/superoxide/ozonide
                corrections are to be applied or not.
        """
        c = loadfn(config_file)
        self.oxide_correction = c["OxideCorrections"]
        self.sulfide_correction = c.get("SulfideCorrections", defaultdict(float))
        self.name = c["Name"]
        self.correct_peroxide = correct_peroxide

    def get_correction(self, entry) -> ufloat:
        """
        :param entry: A ComputedEntry/ComputedStructureEntry
        :return: Correction.
        """
        comp = entry.composition
        if len(comp) == 1:  # Skip element entry
            return ufloat(0.0, 0.0)

        correction = ufloat(0.0, 0.0)

        # only correct GGA or GGA+U entries
        if entry.parameters.get("run_type", None) not in ["GGA", "GGA+U"]:
            return ufloat(0.0, 0.0)

        # Check for sulfide corrections
        if Element("S") in comp:
            sf_type = "sulfide"

            if entry.data.get("sulfide_type"):
                sf_type = entry.data["sulfide_type"]
            elif hasattr(entry, "structure"):
                warnings.warn(sf_type)
                sf_type = sulfide_type(entry.structure)

            # use the same correction for polysulfides and sulfides
            if sf_type == "polysulfide":
                sf_type = "sulfide"

            if sf_type in self.sulfide_correction:
                correction += self.sulfide_correction[sf_type] * comp["S"]

        # Check for oxide, peroxide, superoxide, and ozonide corrections.
        if Element("O") in comp:
            if self.correct_peroxide:
                if entry.data.get("oxide_type"):
                    if entry.data["oxide_type"] in self.oxide_correction:
                        ox_corr = self.oxide_correction[entry.data["oxide_type"]]
                        correction += ox_corr * comp["O"]
                    if entry.data["oxide_type"] == "hydroxide":
                        ox_corr = self.oxide_correction["oxide"]
                        correction += ox_corr * comp["O"]

                elif hasattr(entry, "structure"):
                    ox_type, nbonds = oxide_type(entry.structure, 1.05, return_nbonds=True)
                    if ox_type in self.oxide_correction:
                        correction += self.oxide_correction[ox_type] * nbonds
                    elif ox_type == "hydroxide":
                        correction += self.oxide_correction["oxide"] * comp["O"]
                else:
                    warnings.warn(
                        "No structure or oxide_type parameter present. Note "
                        "that peroxide/superoxide corrections are not as "
                        "reliable and relies only on detection of special"
                        "formulas, e.g., Li2O2."
                    )
                    rform = entry.composition.reduced_formula
                    if rform in UCorrection.common_peroxides:
                        correction += self.oxide_correction["peroxide"] * comp["O"]
                    elif rform in UCorrection.common_superoxides:
                        correction += self.oxide_correction["superoxide"] * comp["O"]
                    elif rform in UCorrection.ozonides:
                        correction += self.oxide_correction["ozonide"] * comp["O"]
                    elif Element("O") in comp.elements and len(comp.elements) > 1:
                        correction += self.oxide_correction["oxide"] * comp["O"]
            else:
                correction += self.oxide_correction["oxide"] * comp["O"]

        return correction

    def __str__(self):
        return "{} Anion Correction".format(self.name)


@cached_class
class AqueousCorrection(Correction):
    """
    This class implements aqueous phase compound corrections for elements
    and H2O.

    Used only by MITAqueousCompatibility.
    """

    def __init__(self, config_file, error_file=None):
        """
        Args:
            config_file: Path to the selected compatibility.yaml config file.
            error_file: Path to the selected compatibilityErrors.yaml config file.
        """
        c = loadfn(config_file)
        self.cpd_energies = c["AqueousCompoundEnergies"]
        # there will either be a CompositionCorrections OR an OxideCorrections key,
        # but not both, depending on the compatibility scheme we are using.
        # MITCompatibility only uses OxideCorrections, and hence self.comp_correction is none.
        self.comp_correction = c.get("CompositionCorrections", defaultdict(float))
        self.oxide_correction = c.get("OxideCorrections", defaultdict(float))
        self.name = c["Name"]
        if error_file:
            e = loadfn(error_file)
            self.cpd_errors = e.get("AqueousCompoundEnergies", defaultdict(float))
        else:
            self.cpd_errors = defaultdict(float)

    def get_correction(self, entry) -> ufloat:
        """
        :param entry: A ComputedEntry/ComputedStructureEntry
        :return: Correction, Uncertainty.
        """
        from pymatgen.analysis.pourbaix_diagram import MU_H2O

        comp = entry.composition
        rform = comp.reduced_formula
        cpdenergies = self.cpd_energies

        # only correct GGA or GGA+U entries
        if entry.parameters.get("run_type", None) not in ["GGA", "GGA+U"]:
            return ufloat(0.0, 0.0)

        correction = ufloat(0.0, 0.0)

        if rform in cpdenergies:
            if rform in ["H2", "H2O"]:
                corr = cpdenergies[rform] * comp.num_atoms - entry.uncorrected_energy - entry.correction
                err = self.cpd_errors[rform] * comp.num_atoms

                correction += ufloat(corr, err)
            else:
                corr = cpdenergies[rform] * comp.num_atoms
                err = self.cpd_errors[rform] * comp.num_atoms

                correction += ufloat(corr, err)
        if not rform == "H2O":
            # if the composition contains water molecules (e.g. FeO.nH2O),
            # correct the gibbs free energy such that the waters are assigned energy=MU_H2O
            # in other words, we assume that the DFT energy of such a compound is really
            # a superposition of the "real" solid DFT energy (FeO in this case) and the free
            # energy of some water molecules
            # e.g. that E_FeO.nH2O = E_FeO + n * g_H2O
            # so, to get the most accurate gibbs free energy, we want to replace
            # g_FeO.nH2O = E_FeO.nH2O + dE_Fe + (n+1) * dE_O + 2n dE_H
            # with
            # g_FeO = E_FeO.nH2O + dE_Fe + dE_O + n g_H2O
            # where E is DFT energy, dE is an energy correction, and g is gibbs free energy
            # This means we have to 1) remove energy corrections associated with H and O in water
            # and then 2) remove the free energy of the water molecules

            nH2O = int(min(comp["H"] / 2.0, comp["O"]))  # only count whole water molecules
            if nH2O > 0:
                # first, remove any H or O corrections already applied to H2O in the
                # formation energy so that we don't double count them
                # No. of H atoms not in a water
                correction -= ufloat((comp["H"] - nH2O / 2) * self.comp_correction["H"], 0.0)
                # No. of O atoms not in a water
                correction -= ufloat(
                    (comp["O"] - nH2O) * (self.comp_correction["oxide"] + self.oxide_correction["oxide"]),
                    0.0,
                )
                # next, add MU_H2O for each water molecule present
                correction += ufloat(-1 * MU_H2O * nH2O, 0.0)
                # correction += 0.5 * 2.46 * nH2O  # this is the old way this correction was calculated

        return correction

    def __str__(self):
        return "{} Aqueous Correction".format(self.name)


@cached_class
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
    """

    common_peroxides = [
        "Li2O2",
        "Na2O2",
        "K2O2",
        "Cs2O2",
        "Rb2O2",
        "BeO2",
        "MgO2",
        "CaO2",
        "SrO2",
        "BaO2",
    ]
    common_superoxides = ["LiO2", "NaO2", "KO2", "RbO2", "CsO2"]
    ozonides = ["LiO3", "NaO3", "KO3", "NaO5"]

    def __init__(self, config_file, input_set, compat_type, error_file=None):
        """
        Args:
            config_file: Path to the selected compatibility.yaml config file.
            input_set: InputSet object (to check for the +U settings)
            compat_type: Two options, GGA or Advanced.  GGA means all GGA+U
                entries are excluded.  Advanced means mixing scheme is
                implemented to make entries compatible with each other,
                but entries which are supposed to be done in GGA+U will have the
                equivalent GGA entries excluded. For example, Fe oxides should
                have a U value under the Advanced scheme. A GGA Fe oxide run
                will therefore be excluded under the scheme.
            error_file: Path to the selected compatibilityErrors.yaml config file.
        """
        if compat_type not in ["GGA", "Advanced"]:
            raise CompatibilityError("Invalid compat_type {}".format(compat_type))

        c = loadfn(config_file)

        self.input_set = input_set
        if compat_type == "Advanced":
            self.u_settings = self.input_set.CONFIG["INCAR"]["LDAUU"]
            self.u_corrections = c["Advanced"]["UCorrections"]
        else:
            self.u_settings = {}
            self.u_corrections = {}

        self.name = c["Name"]
        self.compat_type = compat_type

        if error_file:
            e = loadfn(error_file)
            self.u_errors = e["Advanced"]["UCorrections"]
        else:
            self.u_errors = {}

    def get_correction(self, entry) -> ufloat:
        """
        :param entry: A ComputedEntry/ComputedStructureEntry
        :return: Correction, Uncertainty.
        """
        if entry.parameters.get("run_type") not in ["GGA", "GGA+U"]:
            raise CompatibilityError(
                "Entry {} has invalid run type {}. Discarding.".format(entry.entry_id, entry.parameters.get("run_type"))
            )

        calc_u = entry.parameters.get("hubbards", None)
        calc_u = defaultdict(int) if calc_u is None else calc_u
        comp = entry.composition

        elements = sorted([el for el in comp.elements if comp[el] > 0], key=lambda el: el.X)
        most_electroneg = elements[-1].symbol
        correction = ufloat(0.0, 0.0)

        # only correct GGA or GGA+U entries
        if entry.parameters.get("run_type", None) not in ["GGA", "GGA+U"]:
            return ufloat(0.0, 0.0)

        ucorr = self.u_corrections.get(most_electroneg, {})
        usettings = self.u_settings.get(most_electroneg, {})
        uerrors = self.u_errors.get(most_electroneg, defaultdict(float))

        for el in comp.elements:
            sym = el.symbol
            # Check for bad U values
            if calc_u.get(sym, 0) != usettings.get(sym, 0):
                raise CompatibilityError("Invalid U value of %s on %s" % (calc_u.get(sym, 0), sym))
            if sym in ucorr:
                correction += ufloat(ucorr[sym], uerrors[sym]) * comp[el]

        return correction

    def __str__(self):
        return "{} {} Correction".format(self.name, self.compat_type)


class Compatibility(MSONable, metaclass=abc.ABCMeta):
    """
    Abstract Compatibility class, not intended for direct use.
    Compatibility classes are used to correct the energies of an entry or a set
    of entries. All Compatibility classes must implement .get_adjustments method.
    """

    @abc.abstractmethod
    def get_adjustments(self, entry: ComputedEntry):
        """
        Get the energy adjustments for a ComputedEntry.

        This method must generate a list of EnergyAdjustment objects
        of the appropriate type (constant, composition-based, or temperature-based)
        to be applied to the ComputedEntry, and must raise a CompatibilityError
        if the entry is not compatible.

        Args:
            entry: A ComputedEntry object.

        Returns:
            [EnergyAdjustment]: A list of EnergyAdjustment to be applied to the
                Entry.

        Raises:
            CompatibilityError if the entry is not compatible
        """
        return

    def process_entry(self, entry):
        """
        Process a single entry with the chosen Corrections. Note
        that this method will change the data of the original entry.

        Args:
            entry: A ComputedEntry object.
        Returns:
            An adjusted entry if entry is compatible, otherwise None is
            returned.
        """
        if self.process_entries(entry):
            return self.process_entries(entry)[0]
        return None

    def process_entries(self, entries: Union[ComputedEntry, list], clean: bool = True, verbose: bool = False):
        """
        Process a sequence of entries with the chosen Compatibility scheme. Note
        that this method will change the data of the original entries.

        Args:
            entries: ComputedEntry or [ComputedEntry]
            clean: bool, whether to remove any previously-applied energy adjustments.
                If True, all EnergyAdjustment are removed prior to processing the Entry.
                Default is True.
            verbose: bool, whether to display progress bar for processing multiple entries.
                Default is False.

        Returns:
            A list of adjusted entries.  Entries in the original list which
            are not compatible are excluded.
        """
        # convert input arg to a list if not already
        if isinstance(entries, ComputedEntry):
            entries = [entries]

        processed_entry_list = []

        for entry in PBar(entries, disable=(not verbose)):
            ignore_entry = False
            # if clean is True, remove all previous adjustments from the entry
            if clean:
                entry.energy_adjustments = []

            # get the energy adjustments
            try:
                adjustments = self.get_adjustments(entry)
            except CompatibilityError:
                ignore_entry = True
                continue

            for ea in adjustments:
                # Has this correction already been applied?
                if (ea.name, ea.cls, ea.value) in [(ea.name, ea.cls, ea.value) for ea in entry.energy_adjustments]:
                    # we already applied this exact correction. Do nothing.
                    pass
                elif (ea.name, ea.cls) in [(ea.name, ea.cls) for ea in entry.energy_adjustments]:
                    # we already applied a correction with the same name
                    # but a different value. Something is wrong.
                    ignore_entry = True
                    warnings.warn(
                        "Entry {} already has an energy adjustment called {}, but its "
                        "value differs from the value of {:.3f} calculated here. This "
                        "Entry will be discarded.".format(entry.entry_id, ea.name, ea.value)
                    )
                else:
                    # Add the correction to the energy_adjustments list
                    entry.energy_adjustments.append(ea)

            if not ignore_entry:
                processed_entry_list.append(entry)

        return processed_entry_list

    @staticmethod
    def explain(entry):
        """
        Prints an explanation of the energy adjustments applied by the
        Compatibility class. Inspired by the "explain" methods in many database
        methodologies.

        Args:
            entry: A ComputedEntry.
        """
        print(
            "The uncorrected energy of {} is {:.3f} eV ({:.3f} eV/atom).".format(
                entry.composition,
                entry.uncorrected_energy,
                entry.uncorrected_energy / entry.composition.num_atoms,
            )
        )

        if len(entry.energy_adjustments) > 0:
            print("The following energy adjustments have been applied to this entry:")
            for e in entry.energy_adjustments:
                print(
                    "\t\t{}: {:.3f} eV ({:.3f} eV/atom)".format(e.name, e.value, e.value / entry.composition.num_atoms)
                )
        elif entry.correction == 0:
            print("No energy adjustments have been applied to this entry.")

        print(
            "The final energy after adjustments is {:.3f} eV ({:.3f} eV/atom).".format(
                entry.energy, entry.energy_per_atom
            )
        )


class CorrectionsList(Compatibility):
    """
    The CorrectionsList class combines a list of corrections to be applied to
    an entry or a set of entries. Note that some of the Corrections have
    interdependencies. For example, PotcarCorrection must always be used
    before any other compatibility. Also, AnionCorrection("MP") must be used
    with PotcarCorrection("MP") (similarly with "MIT"). Typically,
    you should use the specific MaterialsProjectCompatibility and
    MITCompatibility subclasses instead.
    """

    def __init__(self, corrections: Sequence):
        """
        Args:
            corrections: List of corrections to apply.
        """
        self.corrections = corrections
        super().__init__()

    def get_adjustments(self, entry):
        """
        Get the list of energy adjustments to be applied to an entry.
        """
        adjustment_list = []
        corrections, uncertainties = self.get_corrections_dict(entry)

        for k, v in corrections.items():
            if v != 0 and uncertainties[k] == 0:
                uncertainty = np.nan
            else:
                uncertainty = uncertainties[k]
            adjustment_list.append(
                ConstantEnergyAdjustment(
                    v,
                    uncertainty=uncertainty,
                    name=k,
                    cls=self.as_dict(),
                )
            )

        return adjustment_list

    def get_corrections_dict(self, entry):
        """
        Returns the corrections applied to a particular entry.

        Args:
            entry: A ComputedEntry object.

        Returns:
            ({correction_name: value})
        """
        corrections = {}
        uncertainties = {}
        for c in self.corrections:
            val = c.get_correction(entry)
            if val != 0:
                corrections[str(c)] = val.nominal_value
                uncertainties[str(c)] = val.std_dev
        return corrections, uncertainties

    def get_explanation_dict(self, entry):
        """
        Provides an explanation dict of the corrections that are being applied
        for a given compatibility scheme. Inspired by the "explain" methods
        in many database methodologies.

        Args:
            entry: A ComputedEntry.

        Returns:
            (dict) of the form
            {"Compatibility": "string",
            "Uncorrected_energy": float,
            "Corrected_energy": float,
            "correction_uncertainty:" float,
            "Corrections": [{"Name of Correction": {
            "Value": float, "Explanation": "string", "Uncertainty": float}]}
        """
        centry = self.process_entry(entry)
        if centry is None:
            uncorrected_energy = entry.uncorrected_energy
            corrected_energy = None
            correction_uncertainty = None
        else:
            uncorrected_energy = centry.uncorrected_energy
            corrected_energy = centry.energy
            correction_uncertainty = centry.correction_uncertainty
        d = {
            "compatibility": self.__class__.__name__,
            "uncorrected_energy": uncorrected_energy,
            "corrected_energy": corrected_energy,
            "correction_uncertainty": correction_uncertainty,
        }
        corrections = []
        corr_dict, uncer_dict = self.get_corrections_dict(entry)
        for c in self.corrections:
            if corr_dict.get(str(c), 0) != 0 and uncer_dict.get(str(c), 0) == 0:
                uncer = np.nan
            else:
                uncer = uncer_dict.get(str(c), 0)
            cd = {
                "name": str(c),
                "description": c.__doc__.split("Args")[0].strip(),
                "value": corr_dict.get(str(c), 0),
                "uncertainty": uncer,
            }
            corrections.append(cd)
        d["corrections"] = corrections
        return d

    def explain(self, entry):
        """
        Prints an explanation of the corrections that are being applied for a
        given compatibility scheme. Inspired by the "explain" methods in many
        database methodologies.

        Args:
            entry: A ComputedEntry.
        """
        d = self.get_explanation_dict(entry)
        print("The uncorrected value of the energy of %s is %f eV" % (entry.composition, d["uncorrected_energy"]))
        print("The following corrections / screening are applied for %s:\n" % d["compatibility"])
        for c in d["corrections"]:
            print("%s correction: %s\n" % (c["name"], c["description"]))
            print("For the entry, this correction has the value %f eV." % c["value"])
            if c["uncertainty"] != 0 or c["value"] == 0:
                print("This correction has an uncertainty value of %f eV." % c["uncertainty"])
            else:
                print("This correction does not have uncertainty data available")
            print("-" * 30)

        print("The final energy after corrections is %f" % d["corrected_energy"])


class MaterialsProjectCompatibility(CorrectionsList):
    """
    This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. Note that this should only be used for VASP calculations using the
    MaterialsProject parameters (see pymatgen.io.vaspio_set.MPVaspInputSet).
    Using this compatibility scheme on runs with different parameters is not
    valid.
    """

    def __init__(self, compat_type="Advanced", correct_peroxide=True, check_potcar_hash=False):
        """
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
            check_potcar_hash (bool): Use potcar hash to verify potcars are correct.
        """
        self.compat_type = compat_type
        self.correct_peroxide = correct_peroxide
        self.check_potcar_hash = check_potcar_hash
        fp = os.path.join(MODULE_DIR, "MPCompatibility.yaml")
        super().__init__(
            [
                PotcarCorrection(MPRelaxSet, check_hash=check_potcar_hash),
                GasCorrection(fp),
                AnionCorrection(fp, correct_peroxide=correct_peroxide),
                UCorrection(fp, MPRelaxSet, compat_type),
            ]
        )


class MaterialsProject2020Compatibility(Compatibility):
    """
    This class implements the Materials Project 2020 energy correction scheme,
    which incorporates uncertainty quantification and allows for mixing of GGA
    and GGA+U entries (see References).

    Note that this scheme should only be applied to VASP calculations that use the
    Materials Project input set parameters (see pymatgen.io.vasp.sets.MPRelaxSet).
    Using this compatibility scheme on calculations with different parameters is not
    valid.
    """

    def __init__(
        self,
        compat_type="Advanced",
        correct_peroxide=True,
        check_potcar_hash=False,
        config_file=None,
    ):
        """
        Args:
            compat_type: Two options, GGA or Advanced.  GGA means all GGA+U
                entries are excluded. Advanced means the GGA/GGA+U mixing scheme
                of Jain et al. (see References) is implemented. In this case,
                entries which are supposed to be calculated in GGA+U (i.e.,
                transition metal oxides and fluorides) will have the corresponding
                GGA entries excluded. For example, Fe oxides should
                have a U value under the Advanced scheme. An Fe oxide run in GGA
                will therefore be excluded.

                To use the "Advanced" type, Entry.parameters must contain a "hubbards"
                key which is a dict of all non-zero Hubbard U values used in the
                calculation. For example, if you ran a Fe2O3 calculation with
                Materials Project parameters, this would look like
                entry.parameters["hubbards"] = {"Fe": 5.3}. If the "hubbards" key
                is missing, a GGA run is assumed. Entries obtained from the
                MaterialsProject database will automatically have these fields
                populated.

                (Default: "Advanced")
            correct_peroxide: Specify whether peroxide/superoxide/ozonide
                corrections are to be applied or not. If false, all oxygen-containing
                compounds are assigned the 'oxide' correction. (Default: True)
            check_potcar_hash (bool): Use potcar hash to verify POTCAR settings are
                consistent with MPRelaxSet. If False, only the POTCAR symbols will
                be used. (Default: False)
            config_file (Path): Path to the selected compatibility.yaml config file.
                If None, defaults to `MP2020Compatibility.yaml` distributed with
                pymatgen.

        References:
            Wang, A., Kingsbury, R., McDermott, M., Horton, M., Jain. A., Ong, S.P.,
                Dwaraknath, S., Persson, K. A framework for quantifying uncertainty
                in DFT energy corrections. Scientific Reports 11: 15496, 2021.
                https://doi.org/10.1038/s41598-021-94550-5

            Jain, A. et al. Formation enthalpies by mixing GGA and GGA + U calculations.
                Phys. Rev. B - Condens. Matter Mater. Phys. 84, 1–10 (2011).
        """
        if compat_type not in ["GGA", "Advanced"]:
            raise CompatibilityError("Invalid compat_type {}".format(compat_type))

        self.compat_type = compat_type
        self.correct_peroxide = correct_peroxide
        self.check_potcar_hash = check_potcar_hash

        # load corrections and uncertainties
        if config_file:
            if os.path.isfile(config_file):
                self.config_file = config_file
                c = loadfn(self.config_file)
            else:
                raise ValueError(
                    f"Custom MaterialsProject2020Compatibility config_file ({config_file}) does not exist."
                )
        else:
            self.config_file = None
            c = loadfn(os.path.join(MODULE_DIR, "MP2020Compatibility.yaml"))

        self.name = c["Name"]
        self.comp_correction = c["Corrections"].get("CompositionCorrections", defaultdict(float))
        self.comp_errors = c["Uncertainties"].get("CompositionCorrections", defaultdict(float))

        if self.compat_type == "Advanced":
            self.u_settings = MPRelaxSet.CONFIG["INCAR"]["LDAUU"]
            self.u_corrections = c["Corrections"].get("GGAUMixingCorrections", defaultdict(float))
            self.u_errors = c["Uncertainties"].get("GGAUMixingCorrections", defaultdict(float))
        else:
            self.u_settings = {}
            self.u_corrections = {}
            self.u_errors = {}

    def get_adjustments(self, entry: Union[ComputedEntry, ComputedStructureEntry]):
        """
        Get the energy adjustments for a ComputedEntry or ComputedStructureEntry.

        Energy corrections are implemented directly in this method instead of in
        separate AnionCorrection, GasCorrection, or UCorrection classes which
        were used in the legacy correction scheme.

        Args:
            entry: A ComputedEntry or ComputedStructureEntry object.

        Returns:
            [EnergyAdjustment]: A list of EnergyAdjustment to be applied to the
                Entry.

        Raises:
            CompatibilityError if the entry is not compatible
        """
        if entry.parameters.get("run_type") not in ["GGA", "GGA+U"]:
            raise CompatibilityError(
                "Entry {} has invalid run type {}. Must be GGA or GGA+U. Discarding.".format(
                    entry.entry_id, entry.parameters.get("run_type")
                )
            )

        # check the POTCAR symbols
        # this should return ufloat(0, 0) or raise a CompatibilityError or ValueError
        pc = PotcarCorrection(MPRelaxSet, check_hash=self.check_potcar_hash)
        pc.get_correction(entry)

        # apply energy adjustments
        adjustments: List[CompositionEnergyAdjustment] = []

        comp = entry.composition
        rform = comp.reduced_formula
        # sorted list of elements, ordered by electronegativity
        elements = sorted([el for el in comp.elements if comp[el] > 0], key=lambda el: el.X)

        # Skip single elements
        if len(comp) == 1:
            return adjustments

        # Check for sulfide corrections
        if Element("S") in comp:
            sf_type = "sulfide"
            if entry.data.get("sulfide_type"):
                sf_type = entry.data["sulfide_type"]
            elif hasattr(entry, "structure"):
                sf_type = sulfide_type(entry.structure)

            # use the same correction for polysulfides and sulfides
            if sf_type == "polysulfide":
                sf_type = "sulfide"

            if sf_type == "sulfide":
                adjustments.append(
                    CompositionEnergyAdjustment(
                        self.comp_correction["S"],
                        comp["S"],
                        uncertainty_per_atom=self.comp_errors["S"],
                        name="MP2020 anion correction (S)",
                        cls=self.as_dict(),
                    )
                )

        # Check for oxide, peroxide, superoxide, and ozonide corrections.
        if Element("O") in comp:
            if self.correct_peroxide:
                # determine the oxide_type
                if entry.data.get("oxide_type"):
                    ox_type = entry.data["oxide_type"]
                elif hasattr(entry, "structure"):
                    ox_type, nbonds = oxide_type(entry.structure, 1.05, return_nbonds=True)
                else:
                    warnings.warn(
                        "No structure or oxide_type parameter present. Note "
                        "that peroxide/superoxide corrections are not as "
                        "reliable and relies only on detection of special"
                        "formulas, e.g., Li2O2."
                    )

                    common_peroxides = [
                        "Li2O2",
                        "Na2O2",
                        "K2O2",
                        "Cs2O2",
                        "Rb2O2",
                        "BeO2",
                        "MgO2",
                        "CaO2",
                        "SrO2",
                        "BaO2",
                    ]
                    common_superoxides = ["LiO2", "NaO2", "KO2", "RbO2", "CsO2"]
                    ozonides = ["LiO3", "NaO3", "KO3", "NaO5"]

                    if rform in common_peroxides:
                        ox_type = "peroxide"
                    elif rform in common_superoxides:
                        ox_type = "superoxide"
                    elif rform in ozonides:
                        ox_type = "ozonide"
                    else:
                        ox_type = "oxide"
            else:
                ox_type = "oxide"

            if ox_type == "hydroxide":
                ox_type = "oxide"

            adjustments.append(
                CompositionEnergyAdjustment(
                    self.comp_correction[ox_type],
                    comp["O"],
                    uncertainty_per_atom=self.comp_errors[ox_type],
                    name="MP2020 anion correction ({})".format(ox_type),
                    cls=self.as_dict(),
                )
            )

        # Check for anion corrections
        # only apply anion corrections if the element is an anion
        # first check for a pre-populated oxidation states key
        # the key is expected to comprise a dict corresponding to the first element output by
        # Composition.oxi_state_guesses(), e.g. {'Al': 3.0, 'S': 2.0, 'O': -2.0} for 'Al2SO4'
        if "oxidation_states" not in entry.data.keys():
            # try to guess the oxidation states from composition
            # for performance reasons, fail if the composition is too large
            try:
                oxi_states = entry.composition.oxi_state_guesses(max_sites=-20)
            except ValueError:
                oxi_states = []

            if oxi_states == []:
                entry.data["oxidation_states"] = {}
            else:
                entry.data["oxidation_states"] = oxi_states[0]

        if entry.data["oxidation_states"] == {}:
            warnings.warn(
                f"Failed to guess oxidation states for Entry {entry.entry_id} "
                f"({entry.composition.reduced_formula}). Assigning anion correction to "
                "only the most electronegative atom."
            )

        for anion in ["Br", "I", "Se", "Si", "Sb", "Te", "H", "N", "F", "Cl"]:
            if Element(anion) in comp and anion in self.comp_correction:
                apply_correction = False
                # if the oxidation_states key is not populated, only apply the correction if the anion
                # is the most electronegative element
                if entry.data["oxidation_states"].get(anion, 0) < 0:
                    apply_correction = True
                else:
                    most_electroneg = elements[-1].symbol
                    if anion == most_electroneg:
                        apply_correction = True

                if apply_correction:
                    adjustments.append(
                        CompositionEnergyAdjustment(
                            self.comp_correction[anion],
                            comp[anion],
                            uncertainty_per_atom=self.comp_errors[anion],
                            name="MP2020 anion correction ({})".format(anion),
                            cls=self.as_dict(),
                        )
                    )

        # GGA / GGA+U mixing scheme corrections
        calc_u = entry.parameters.get("hubbards", None)
        calc_u = defaultdict(int) if calc_u is None else calc_u
        most_electroneg = elements[-1].symbol
        ucorr = self.u_corrections.get(most_electroneg, defaultdict(float))
        usettings = self.u_settings.get(most_electroneg, defaultdict(float))
        uerrors = self.u_errors.get(most_electroneg, defaultdict(float))

        for el in comp.elements:
            sym = el.symbol
            # Check for bad U values
            if calc_u.get(sym, 0) != usettings.get(sym, 0):
                raise CompatibilityError("Invalid U value of {:.1f} on {}".format(calc_u.get(sym, 0), sym))
            if sym in ucorr:
                adjustments.append(
                    CompositionEnergyAdjustment(
                        ucorr[sym],
                        comp[el],
                        uncertainty_per_atom=uerrors[sym],
                        name="MP2020 GGA/GGA+U mixing correction ({})".format(sym),
                        cls=self.as_dict(),
                    )
                )

        return adjustments


class MaterialsProjectScanCompatibility2020(Compatibility):
    """
    This class implements the SCAN/GGA/GGA+U mixing scheme, which allows mixing of
    entries. Note that this should only be used for VASP calculations using the
    MaterialsProject parameters (see pymatgen.io.vaspio_set.MPVaspInputSet).
    Using this compatibility scheme on runs with different parameters is not
    valid.
    """

    def __init__(
        self,
        structure_matcher: StructureMatcher = StructureMatcher(),
        gga_compat: Optional[Compatibility] = MaterialsProject2020Compatibility(),
        scan_compat: Optional[Compatibility] = None,
        energy_scale: str = "GGA",
    ):
        """
        Args:
            structure_matcher (StructureMatcher): StructureMatcher object used to determine
                whether SCAN calculations exist for GGA reference states.
            gga_compat (Compatibility or None): Compatibility class used to pre-process GGA entries.
                Defaults to MaterialsProjectCompatibility2020.
            scan_compat (Compatibility or None): Compatibility class used to pre-process GGA entries.
                Defaults to None.
            energy_scale (str): Energy scale to correct entries to. Either "SCAN"
                (default) or "GGA".
        """
        self.name = "MP SCAN mixing"
        self.structure_matcher = structure_matcher
        self.gga_compat = gga_compat
        self.scan_compat = scan_compat
        self.energy_scale = energy_scale

    def process_entries(self, entries: Union[ComputedEntry, list], clean: bool = False, verbose: bool = True):
        """
        Args:
            verbose (bool): Print detailed output when processing entries. Default True.
        """
        # State variable to pass to get_adjustments
        gga_entries = []
        scan_entries = []
        gga_ground_states = []
        pd_gga = None
        has_scan_ground_states = []
        has_scan_hull_entries = []

        # We can't operate on single entries in this scheme
        if len(entries) == 1:
            warnings.warn("{} cannot process single entries. Supply a list of entries."
                          .format(self.__class__.__name__)
                          )

        # separate GGA and SCAN entries
        for e in entries:
            # if clean is True, remove all previous adjustments from the entry
            if clean:
                for ea in e.energy_adjustments:
                    e.energy_adjustments.remove(ea)

            if not isinstance(e, ComputedStructureEntry):
                warnings.warn("Entry {} is not a ComputedStructureEntry and will be"
                              "ignored. The SCAN/GGA(+U) mixing scheme requires structures for"
                              " all entries".format(e.entry_id)
                              )
                continue

            if not e.parameters.get("run_type"):
                warnings.warn("Entry {} is missing parameters.run_type! This field"
                              "is required. This entry will be ignored.".format(e.entry_id)
                              )
                continue

            if e.parameters.get("run_type") == "GGA":
                gga_entries.append(e)
            elif e.parameters.get("run_type") == "GGA+U":
                gga_entries.append(e)
            elif e.parameters.get("run_type") == "R2SCAN":
                scan_entries.append(e)
            else:
                warnings.warn("Invalid run_type {} for entry {}. This entry "
                              "will be ignored".format(e.parameters.get("run_type"),
                                                       e.entry_id)
                              )
                continue

        if verbose:
            print("Processing {} valid GGA(+U) and {} valid SCAN entries".format(len(gga_entries),
                                                                                 len(scan_entries)))

        # preprocess entries with any corrections
        if self.gga_compat:
            gga_entries = self.gga_compat.process_entries(gga_entries)
            if verbose:
                print("Processed {} compatible GGA(+U) entries with {}".format(len(gga_entries),
                                                                               self.gga_compat.__class__.__name__))

        if self.scan_compat:
            scan_entries = self.scan_compat.process_entries(scan_entries)
            if verbose:
                print("Processed {} compatible SCAN entries with {}".format(len(scan_entries),
                                                                            self.scan_compat.__class__.__name__))

        # Make deep copies of the entries at this point, because later we'll need their original, corrected
        # energies prior to mixing functionals
        gga_entries = EntrySet([copy.deepcopy(e) for e in gga_entries])
        scan_entries = EntrySet([copy.deepcopy(e) for e in scan_entries])

        # make sure both sets of entries belong to the same chemical system
        # assuming there are any gga entries at all
        if len(gga_entries.chemsys) > 0:
            chemsys = gga_entries.chemsys
            if not scan_entries.chemsys <= gga_entries.chemsys:
                warnings.warn("SCAN entries chemical system {} is larger than GGA entries"
                              " chemical system {}. Entries outside the GGA chemical system"
                              " will be discarded".format(scan_entries.chemsys,
                                                          gga_entries.chemsys)
                              )
            # If GGA entries are present, construct a gga phase diagram and
            # identify the GGA ground states
            try:
                pd_gga = PhaseDiagram(gga_entries)
                gga_hull_entries = pd_gga.stable_entries
            except ValueError:
                pd_gga = None
                gga_hull_entries = []
                warnings.warn("GGA entries do not form a complete PhaseDiagram!")
        else:
            chemsys = scan_entries.chemsys
            gga_ground_states = []
            gga_hull_entries = []

        # bool lists indicating whether there is a SCAN entry for each GGA ground state
        # and each GGA elemental reference
        gga_ground_states = gga_entries
        gga_ground_states.remove_non_ground_states()

        for e in gga_ground_states:
            has_scan = False
            for e2 in scan_entries:
                if self.structure_matcher.fit(e.structure, e2.structure):
                    has_scan = True
                    break

            has_scan_ground_states.append(has_scan)
            if e in gga_hull_entries:
                has_scan_hull_entries.append(has_scan)

        if verbose:
            print("Entries belong to the {} chemical system".format(chemsys))
            print("Entries contain SCAN calculations for {} of {} GGA(+U) hull "
                  "entries.".format(sum(has_scan_hull_entries),
                                    len(gga_hull_entries),
                                    )
                  )

        processed_entry_list = []

        for entry in entries:
            ignore_entry = False
            # get the energy adjustments
            try:
                adjustments = self.get_adjustments(entry,
                                                   gga_entries,
                                                   scan_entries,
                                                   gga_ground_states,
                                                   pd_gga,
                                                   has_scan_ground_states,
                                                   has_scan_hull_entries
                                                   )
            except CompatibilityError as exc:
                ignore_entry = True
                print(exc)
                continue

            for ea in adjustments:
                # Has this correction already been applied?
                if (ea.name, ea.cls, ea.value) in [(ea.name, ea.cls, ea.value) for ea in entry.energy_adjustments]:
                    # we already applied this exact correction. Do nothing.
                    pass
                elif (ea.name, ea.cls) in [(ea.name, ea.cls) for ea in entry.energy_adjustments]:
                    # we already applied a correction with the same name
                    # but a different value. Something is wrong.
                    ignore_entry = True
                    warnings.warn("Entry {} already has an energy adjustment called {}, but its "
                                  "value differs from the value of {:.3f} calculated here. This "
                                  "Entry will be discarded."
                                  .format(entry.entry_id,
                                          ea.name,
                                          ea.value
                                          )
                                  )
                else:
                    # Add the correction to the energy_adjustments list
                    entry.energy_adjustments.append(ea)

            if not ignore_entry:
                processed_entry_list.append(entry)

        return processed_entry_list

    def get_adjustments(self,
                        entry: ComputedStructureEntry,
                        gga_entries: EntrySet,
                        scan_entries: EntrySet,
                        gga_ground_states: List,
                        pd_gga: Optional[PhaseDiagram],
                        has_scan_ground_states: List,
                        has_scan_hull_entries: List
                        ):
        """
        Returns the corrections applied to a particular entry. Note that get_adjustments is not
        intended to be called directly in the SCAN mixing scheme. Call process_entries instead,
        and it will pass the required arguments to get_adjustments.

        Args:
            entry: A ComputedEntry object.
            gga_entries: EntrySet of GGA(+U) ComputedStructureEntry
            scan_entries: EntrySet of SCAN ComputedStructureEntry
            gga_ground_states: List of ground-state GGA(+U) ComputedStructureEntry. Note that this should include
                ground states for both stable and unstable compositions (i.e., compositions not on the hull)
            pd_gga: Phase diagram constructed from the GGA(+U) ground states
            has_scan_ground_states: List of boolean indicating whether scan_entries contain structures
                that match each structure in gga_ground_states.
            has_scan_hull_entries: List of boolean indicating whether scan_entries contain structures that
                match each structure on the stable energy hull in pd_gga.

        Returns:
            [EnergyAdjustment]: Energy adjustments to be applied to entry.

        Raises:
            CompatibilityError if the SCAN/GGA(+U) mixing scheme cannot be applied to the entry.
        """
        adjustments = []

        run_type = entry.parameters.get("run_type")

        if run_type not in ["GGA", "GGA+U", "R2SCAN"]:
            raise CompatibilityError("Invalid run type {} for entry {}. Must be GGA, GGA+U, "
                                     "or R2SCAN".format(run_type, entry.entry_id)
                                     )

        # The correction value depends on how many of the GGA reference states
        # are present as SCAN calculations
        if all(has_scan_hull_entries) and has_scan_hull_entries != [] or all(has_scan_ground_states):
            # If all SCAN reference are present,
            if run_type == 'R2SCAN':
                # For SCAN entries, there is no correction
                return adjustments

            # Discard GGA entries whose structures already exist in SCAN.
            present = False
            for e in scan_entries:
                if self.structure_matcher.fit(entry.structure, e.structure):
                    present = True
                    raise CompatibilityError(
                            "Discarding GGA entry {} for {} whose structure already "
                            "exists in SCAN".format(entry.entry_id, entry.composition.formula)
                            )

            if not present:
                # If a GGA is not present in SCAN, correct its energy to give the same
                # e_above_hull on the SCAN hull that it would have on the GGA hull
                try:
                    pd_scan = PhaseDiagram(scan_entries)
                    comp = entry.composition
                    scan_hull_energy = pd_scan.get_hull_energy_per_atom(comp)
                except ValueError:
                    raise CompatibilityError("Cannot adjust the energy of GGA(+U) entry {} because the "
                                                "SCAN entries given do not form a complete phase diagram. "
                                                "Discarding this entry.".format(entry.entry_id)
                                                )

                gga_hull_energy = pd_gga.get_hull_energy_per_atom(comp)
                correction = (scan_hull_energy - gga_hull_energy) * entry.composition.num_atoms

                adjustments.append(ConstantEnergyAdjustment(correction, 0.0,
                                                            name="MP GGA/SCAN mixing adjustment",
                                                            cls=self.as_dict(),
                                                            description="Place GGA energy onto the SCAN hull"
                                                            ))
                return adjustments

        elif not any(has_scan_ground_states):
            # Discard SCAN entries if there are no SCAN reference states
            if run_type == 'R2SCAN':
                raise CompatibilityError(
                              "Discarding SCAN entry {} for {} because there are no "
                              "SCAN reference structures".format(entry.entry_id, entry.composition.formula)
                              )
            return adjustments

        # Discard SCAN entries if there are fewer than 2 SCAN entries
        elif len(scan_entries) < 2:
            if run_type == 'R2SCAN':
                raise CompatibilityError(
                            "Discarding SCAN {} entry for {} because there are fewer "
                            "than 2 SCAN entries".format(entry.entry_id, entry.composition.formula)
                            )
            return adjustments

        # if there are some SCAN reference states, we may be able to replace
        # an unstable GGA polymorph with a SCAN one
        # in this scenario, energies remain on the "GGA scale"
        else:
            for ref, present in zip(gga_ground_states, has_scan_ground_states):
                if entry.composition.reduced_composition == ref.composition.reduced_composition:
                    if present:
                        # correct the energy of GGA entries to maintain the same polymorph
                        # energy difference above the reference state
                        if run_type == 'R2SCAN':
                            # find the SCAN reference
                            if self.structure_matcher.fit(ref.structure, entry.structure):
                                # this is the SCAN reference energy
                                # correct the energy to match that of the GGA reference
                                correction = (ref.energy_per_atom - entry.energy_per_atom) * entry.composition.num_atoms
                                adjustments.append(
                                            ConstantEnergyAdjustment(correction, 0.0,
                                                                     name="MP GGA/SCAN mixing adjustment",
                                                                     cls=self.as_dict(),
                                                                     description="Place SCAN energy above hull onto "
                                                                     "the GGA hull"
                                                                     ))
                                return adjustments

                            # search the other SCAN entries to find the reference structure
                            for e in scan_entries:
                                if self.structure_matcher.fit(ref.structure, e.structure):
                                    # correct the SCAN energy to the GGA hull using the references
                                    correction = (ref.energy_per_atom - e.energy_per_atom) * \
                                        entry.composition.num_atoms

                                    adjustments.append(
                                        ConstantEnergyAdjustment(correction, 0.0,
                                                                    name="MP GGA/SCAN mixing adjustment",
                                                                    cls=self.as_dict(),
                                                                    description="Place SCAN energy above hull onto "
                                                                    "the GGA hull"
                                                                    ))
                                    return adjustments
                        else:
                            # is there  a SCAN entry for this GGA structure?
                            for e in scan_entries:
                                if self.structure_matcher.fit(entry.structure, e.structure):
                                    # if yes, discard the GGA entry
                                    raise CompatibilityError(
                                                    "Discarding GGA entry {} for {} whose structure already "
                                                    "exists in SCAN".format(entry.entry_id, entry.composition.formula)
                                                    )
                            else:
                                return adjustments
                    else:
                        # There's no SCAN reference state here, so nothing we can do
                        if run_type == 'R2SCAN':
                            raise CompatibilityError(
                                "Discarding SCAN entry {} for {} for which there is no "
                                "SCAN reference structure".format(entry.entry_id, entry.composition.formula)
                                )
                        return adjustments


class MITCompatibility(CorrectionsList):
    """
    This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. Note that this should only be used for VASP calculations using the
    MIT parameters (see pymatgen.io.vaspio_set MITVaspInputSet). Using
    this compatibility scheme on runs with different parameters is not valid.
    """

    def __init__(self, compat_type="Advanced", correct_peroxide=True, check_potcar_hash=False):
        """
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
            check_potcar_hash (bool): Use potcar hash to verify potcars are correct.
        """
        self.compat_type = compat_type
        self.correct_peroxide = correct_peroxide
        self.check_potcar_hash = check_potcar_hash
        fp = os.path.join(MODULE_DIR, "MITCompatibility.yaml")
        super().__init__(
            [
                PotcarCorrection(MITRelaxSet, check_hash=check_potcar_hash),
                GasCorrection(fp),
                AnionCorrection(fp, correct_peroxide=correct_peroxide),
                UCorrection(fp, MITRelaxSet, compat_type),
            ]
        )


class MITAqueousCompatibility(CorrectionsList):
    """
    This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. Note that this should only be used for VASP calculations using the
    MIT parameters (see pymatgen.io.vaspio_set MITVaspInputSet). Using
    this compatibility scheme on runs with different parameters is not valid.
    """

    def __init__(self, compat_type="Advanced", correct_peroxide=True, check_potcar_hash=False):
        """
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
            check_potcar_hash (bool): Use potcar hash to verify potcars are correct.
        """
        self.compat_type = compat_type
        self.correct_peroxide = correct_peroxide
        self.check_potcar_hash = check_potcar_hash
        fp = os.path.join(MODULE_DIR, "MITCompatibility.yaml")
        super().__init__(
            [
                PotcarCorrection(MITRelaxSet, check_hash=check_potcar_hash),
                GasCorrection(fp),
                AnionCorrection(fp, correct_peroxide=correct_peroxide),
                UCorrection(fp, MITRelaxSet, compat_type),
                AqueousCorrection(fp),
            ]
        )


class MaterialsProjectAqueousCompatibility(Compatibility):
    """
    This class implements the Aqueous energy referencing scheme for constructing
    Pourbaix diagrams from DFT energies, as described in Persson et al.

    This scheme applies various energy adjustments to convert DFT energies into
    Gibbs free energies of formation at 298 K and to guarantee that the experimental
    formation free energy of H2O is reproduced. Briefly, the steps are:

        1. Beginning with the DFT energy of O2, adjust the energy of H2 so that
           the experimental reaction energy of -2.458 eV/H2O is reproduced.
        2. Add entropy to the DFT energy of any compounds that are liquid or
           gaseous at room temperature
        3. Adjust the energy of H2O for consistency with the adjusted H2 energy.
        4. Adjust the DFT energies of solid hydrate compounds (compounds that
           contain water, e.g. FeO.nH2O) such that the energies of the embedded
           H2O molecules are equal to the experimental free energy

    The above energy adjustments are computed dynamically based on the input
    Entries.

    References:
        K.A. Persson, B. Waldwick, P. Lazic, G. Ceder, Prediction of solid-aqueous
        equilibria: Scheme to combine first-principles calculations of solids with
        experimental aqueous states, Phys. Rev. B - Condens. Matter Mater. Phys.
        85 (2012) 1–12. doi:10.1103/PhysRevB.85.235438.
    """

    def __init__(
        self,
        solid_compat: Optional[Union[Compatibility, Type[Compatibility]]] = MaterialsProject2020Compatibility,
        o2_energy: Optional[float] = None,
        h2o_energy: Optional[float] = None,
        h2o_adjustments: Optional[float] = None,
    ):
        """
        Initialize the MaterialsProjectAqueousCompatibility class.

        Note that this class requires as inputs the ground-state DFT energies of O2 and H2O, plus the value of any
        energy adjustments applied to an H2O molecule. If these parameters are not provided in __init__, they can
        be automatically populated by including ComputedEntry for the ground state of O2 and H2O in a list of entries
        passed to process_entries. process_entries will fail if one or the other is not provided.

        Args:
            solid_compat: Compatiblity scheme used to pre-process solid DFT energies prior to applying aqueous
                energy adjustments. May be passed as a class (e.g. MaterialsProject2020Compatibility) or an instance
                (e.g., MaterialsProject2020Compatibility()). If None, solid DFT energies are used as-is.
                Default: MaterialsProject2020Compatibility
            o2_energy: The ground-state DFT energy of oxygen gas, including any adjustments or corrections, in eV/atom.
                If not set, this value will be determined from any O2 entries passed to process_entries.
                Default: None
            h2o_energy: The ground-state DFT energy of water, including any adjstments or corrections, in eV/atom.
                If not set, this value will be determined from any H2O entries passed to process_entries.
                Default: None
            h2o_adjustments: Total energy adjustments applied to one water molecule, in eV/atom.
                If not set, this value will be determined from any H2O entries passed to process_entries.
                Default: None
        """
        self.solid_compat = None
        # check whether solid_compat has been instantiated
        if solid_compat is None:
            self.solid_compat = None
        elif isinstance(solid_compat, type) and issubclass(solid_compat, Compatibility):
            self.solid_compat = solid_compat()
        elif issubclass(type(solid_compat), Compatibility):
            self.solid_compat = solid_compat
        else:
            raise ValueError("Expected a Compatability class, instance of a Compatability or None")

        self.o2_energy = o2_energy
        self.h2o_energy = h2o_energy
        self.h2o_adjustments = h2o_adjustments

        if not all([self.o2_energy, self.h2o_energy, self.h2o_adjustments]):
            warnings.warn(
                "You did not provide the required O2 and H2O energies. {} "
                "needs these energies in order to compute the appropriate energy adjustments. It will try "
                "to determine the values from ComputedEntry for O2 and H2O passed to process_entries, but "
                "will fail if these entries are not provided.".format(type(self).__name__)
            )

        # Standard state entropy of molecular-like compounds at 298K (-T delta S)
        # from Kubaschewski Tables (eV/atom)
        self.cpd_entropies = {
            "O2": 0.316731,
            "N2": 0.295729,
            "F2": 0.313025,
            "Cl2": 0.344373,
            "Br": 0.235039,
            "Hg": 0.234421,
            "H2O": 0.071963,  # 0.215891 eV/H2O
        }
        self.name = "MP Aqueous free energy adjustment"
        super().__init__()

    def get_adjustments(self, entry: ComputedEntry):
        """
        Returns the corrections applied to a particular entry.

        Args:
            entry: A ComputedEntry object.

        Returns:
            [EnergyAdjustment]: Energy adjustments to be applied to entry.

        Raises:
            CompatibilityError if the required O2 and H2O energies have not been provided to
            MaterialsProjectAqueousCompatibility during init or in the list of entries passed to process_entries.
        """
        adjustments = []
        if self.o2_energy is None or self.h2o_energy is None or self.h2o_adjustments is None:
            raise CompatibilityError(
                "You did not provide the required O2 and H2O energies. "
                "{} needs these energies in order to compute "
                "the appropriate energy adjustments. Either specify the energies as arguments "
                "to {}.__init__ or run process_entries on a list that includes ComputedEntry for "
                "the ground state of O2 and H2O.".format(type(self).__name__, type(self).__name__)
            )

        # compute the free energies of H2 and H2O (eV/atom) to guarantee that the
        # formationfree energy of H2O is equal to -2.4583 eV/H2O from experiments
        # (MU_H2O from pourbaix module)

        # Free energy of H2 in eV/atom, fitted using Eq. 40 of Persson et al. PRB 2012 85(23)
        # for this calculation ONLY, we need the (corrected) DFT energy of water
        self.h2_energy = round(
            0.5
            * (
                3 * (self.h2o_energy - self.cpd_entropies["H2O"]) - (self.o2_energy - self.cpd_entropies["O2"]) - MU_H2O
            ),
            6,
        )

        # Free energy of H2O, fitted for consistency with the O2 and H2 energies.
        self.fit_h2o_energy = round(
            (2 * self.h2_energy + (self.o2_energy - self.cpd_entropies["O2"]) + MU_H2O) / 3,
            6,
        )

        comp = entry.composition
        rform = comp.reduced_formula

        # pin the energy of all H2 entries to h2_energy
        if rform == "H2":
            adjustments.append(
                ConstantEnergyAdjustment(
                    self.h2_energy * comp.num_atoms - entry.energy,
                    uncertainty=np.nan,
                    name="MP Aqueous H2 / H2O referencing",
                    cls=self.as_dict(),
                    description="Adjusts the H2 and H2O energy to reproduce the experimental "
                    "Gibbs formation free energy of H2O, based on the DFT energy "
                    "of Oxygen",
                )
            )

        # pin the energy of all H2O entries to fit_h2o_energy
        elif rform == "H2O":
            adjustments.append(
                ConstantEnergyAdjustment(
                    self.fit_h2o_energy * comp.num_atoms - entry.energy,
                    uncertainty=np.nan,
                    name="MP Aqueous H2 / H2O referencing",
                    cls=self.as_dict(),
                    description="Adjusts the H2 and H2O energy to reproduce the experimental "
                    "Gibbs formation free energy of H2O, based on the DFT energy "
                    "of Oxygen",
                )
            )

        # add minus T delta S to the DFT energy (enthalpy) of compounds that are
        # molecular-like at room temperature
        elif rform in self.cpd_entropies and rform != "H2O":
            adjustments.append(
                TemperatureEnergyAdjustment(
                    -1 * self.cpd_entropies[rform] / 298,
                    298,
                    comp.num_atoms,
                    uncertainty_per_deg=np.nan,
                    name="Compound entropy at room temperature",
                    cls=self.as_dict(),
                    description="Adds the entropy (T delta S) to energies of compounds that "
                    "are gaseous or liquid at standard state",
                )
            )

        # TODO - detection of embedded water molecules is not very sophisticated
        # Should be replaced with some kind of actual structure detection

        # For any compound except water, check to see if it is a hydrate (contains)
        # H2O in its structure. If so, adjust the energy to remove MU_H2O ev per
        # embedded water molecule.
        # in other words, we assume that the DFT energy of such a compound is really
        # a superposition of the "real" solid DFT energy (FeO in this case) and the free
        # energy of some water molecules
        # e.g. that E_FeO.nH2O = E_FeO + n * g_H2O
        # so, to get the most accurate gibbs free energy, we want to replace
        # g_FeO.nH2O = E_FeO.nH2O + dE_Fe + (n+1) * dE_O + 2n dE_H
        # with
        # g_FeO = E_FeO.nH2O + dE_Fe + dE_O + n g_H2O
        # where E is DFT energy, dE is an energy correction, and g is gibbs free energy
        # This means we have to 1) remove energy corrections associated with H and O in water
        # and then 2) remove the free energy of the water molecules
        if not rform == "H2O":
            # count the number of whole water molecules in the composition
            nH2O = int(min(comp["H"] / 2.0, comp["O"]))
            if nH2O > 0:
                # first, remove any H or O corrections already applied to H2O in the
                # formation energy so that we don't double count them
                # next, remove MU_H2O for each water molecule present
                hydrate_adjustment = -1 * (self.h2o_adjustments * 3 + MU_H2O)

                adjustments.append(
                    CompositionEnergyAdjustment(
                        hydrate_adjustment,
                        nH2O,
                        uncertainty_per_atom=np.nan,
                        name="MP Aqueous hydrate",
                        cls=self.as_dict(),
                        description="Adjust the energy of solid hydrate compounds (compounds "
                        "containing H2O molecules in their structure) so that the "
                        "free energies of embedded H2O molecules match the experimental"
                        " value enforced by the MP Aqueous energy referencing scheme.",
                    )
                )

        return adjustments

    def process_entries(self, entries: Union[ComputedEntry, list], clean: bool = False, verbose: bool = False):
        """
        Process a sequence of entries with the chosen Compatibility scheme.

        Args:
            entries: ComputedEntry or [ComputedEntry]
            clean: bool, whether to remove any previously-applied energy adjustments.
                If True, all EnergyAdjustment are removed prior to processing the Entry.
                Default is False.
            verbose: bool, whether to display progress bar for processing multiple entries.
                Default is False.

        Returns:
            A list of adjusted entries.  Entries in the original list which
            are not compatible are excluded.
        """
        # convert input arg to a list if not already
        if isinstance(entries, ComputedEntry):
            entries = [entries]

        # pre-process entries with the given solid compatibility class
        if self.solid_compat:
            entries = self.solid_compat.process_entries(entries, clean=True)

        # extract the DFT energies of oxygen and water from the list of entries, if present
        if not self.o2_energy:
            o2_entries = [e for e in entries if e.composition.reduced_formula == "O2"]
            if o2_entries:
                self.o2_energy = min(e.energy_per_atom for e in o2_entries)

        if not self.h2o_energy and not self.h2o_adjustments:
            h2o_entries = [e for e in entries if e.composition.reduced_formula == "H2O"]
            if h2o_entries:
                h2o_entries = sorted(h2o_entries, key=lambda e: e.energy_per_atom)
                self.h2o_energy = h2o_entries[0].energy_per_atom
                self.h2o_adjustments = h2o_entries[0].correction / h2o_entries[0].composition.num_atoms

        return super().process_entries(entries, clean=clean, verbose=verbose)
