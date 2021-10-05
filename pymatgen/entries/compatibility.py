# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
This module implements Compatibility corrections for mixing runs of different
functionals.
"""

import abc
import os
import warnings
from collections import defaultdict
from itertools import groupby
from typing import Optional, Sequence, Union, List, Type
import pandas as pd

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
    EnergyAdjustment,
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
    def get_adjustments(self, entry: Union[ComputedEntry, ComputedStructureEntry]) -> List[EnergyAdjustment]:
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

    def process_entries(
        self, entries: Union[ComputedEntry, ComputedStructureEntry, list], clean: bool = True, verbose: bool = False
    ):
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

    def get_adjustments(self, entry):
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


class MaterialsProjectDFTMixingScheme(Compatibility):
    """
    This class implements the Materials Project mixing scheme, which allows mixing of
    energies from different DFT functionals. Note that this should only be used for
    VASP calculations using the MaterialsProject parameters (e.g. MPRelaxSet or
    MPScanRelaxSet). Using this compatibility scheme on runs with different parameters
    may lead to unexpected results.

    This is the scheme used by the Materials Project to generate Phase Diagrams containing
    a mixture of GGA(+U) and r2SCAN calculations. However in principle it can be used to
    mix energies from any two functionals.
    """

    def __init__(
        self,
        structure_matcher: StructureMatcher = StructureMatcher(),
        run_type_1: str = "GGA(+U)",
        run_type_2: str = "R2SCAN",
        compat_1: Optional[Compatibility] = MaterialsProject2020Compatibility(),
        compat_2: Optional[Compatibility] = None,
        fuzzy_diatomic_matching: bool = True,
    ):
        """
        Instantiate the mixing scheme. The init method creates a generator class that
        contains relevant settings (e.g., StrutureMatcher instance, Compatibility settings
        for each functional) for processing groups of entries.

        Args:
            structure_matcher (StructureMatcher): StructureMatcher object used to determine
                whether calculations from different functionals describe the same material.
            run_type_1: The first DFT run_type. All energies are corrected to the energy scale
                of this run_type. Typically this is the majority or run type or
                the "base case" onto which the other calculations are referenced. Valid choices
                are any run_type recognized by Vasprun.run_type, such as "LDA", "GGA", "GGA+U",
                "PBEsol", "SCAN", or "R2SCAN".

                Note that the special string "GGA(+U)" (default) will treat both GGA and GGA+U
                calculations as a single type. This option exists because GGA/GGA+U mixing is
                already handled by MaterialsProject2020Compatibility.

                The class will ignore any entries that have a run_type different than run_type_1
                or run_type_2.
            run_type_2: The second DFT run_type. Typically this is the run_type that is 'preferred'
                but has fewer calculations. If run_type_1 and run_type_2 calculations exist for all
                materials, run_type_2 energies will be used (hence the 'preferred' status). The class
                will ignore any entries that have a run_type different than run_type_1 or run_type_2.
            compat_1: Compatibility class used to pre-process entries of run_type_1.
                Defaults to MaterialsProjectCompatibility2020.
            compat_2: Compatibility class used to pre-process entries of run_type_2.
                Defaults to None.
            fuzzy_diatomic_matching: Whether to use less strict structure matching logic for
                diatomic elements O2, N2, F2, H2, and Cl2. Outputs of DFT relaxations using
                different functionals frequently fail to struture match for these elements
                even though they come from the same original material. Fuzzy structure matching
                considers the materials equivalent if the formula, number of sites, and
                space group are all identical. If there are multiple materials of run_type_2
                that satisfy these criteria, the one with lowest energy is considered to
                match.
        """
        self.name = "MP SCAN mixing"
        self.structure_matcher = structure_matcher
        self.run_type_1 = run_type_1
        self.run_type_2 = run_type_2
        if self.run_type_1 == "GGA(+U)":
            self.valid_rtypes_1 = ["GGA", "GGA+U"]
        else:
            self.valid_rtypes_1 = [self.run_type_1]

        if self.run_type_2 == "GGA(+U)":
            self.valid_rtypes_2 = ["GGA", "GGA+U"]
        else:
            self.valid_rtypes_2 = [self.run_type_2]

        self.compat_1 = compat_1
        self.compat_2 = compat_2
        self.fuzzy_diatomic_matching = fuzzy_diatomic_matching

    def process_entries(self, entries: Union[ComputedStructureEntry, list], clean: bool = True, verbose: bool = True):
        """
        Process a sequence of entries with the DFT mixing scheme. Note
        that this method will change the data of the original entries.

        Args:
            entries: ComputedEntry or [ComputedEntry]. Pass all entries as a single list, even if they are
                computed with different functionals or require different preprocessing.
            clean: bool, whether to remove any previously-applied energy adjustments.
                If True, all EnergyAdjustment are removed prior to processing the Entry.
                Default is True.
            verbose: bool, whether to display progress bar for processing multiple entries.
                Default is True.

        Returns:
            A list of adjusted entries.  Entries in the original list which
            are not compatible are excluded.
        """
        # We can't operate on single entries in this scheme
        if len(entries) == 1:
            warnings.warn("{} cannot process single entries. Supply a list of entries.".format(self.__class__.__name__))

        filtered_entries = []

        for entry in entries:
            # if clean is True, remove all previous adjustments from the entry
            if clean:
                for ea in entry.energy_adjustments:
                    entry.energy_adjustments.remove(ea)

            if not isinstance(entry, ComputedStructureEntry):
                warnings.warn(
                    "Entry {} is not a ComputedStructureEntry and will be"
                    "ignored. The DFT mixing scheme requires structures for"
                    " all entries".format(entry.entry_id)
                )
                continue

            if not entry.parameters.get("run_type"):
                warnings.warn(
                    "Entry {} is missing parameters.run_type! This field"
                    "is required. This entry will be ignored.".format(entry.entry_id)
                )
                continue

            if entry.parameters.get("run_type") not in self.valid_rtypes_1 + self.valid_rtypes_2:
                warnings.warn(
                    "Invalid run_type {} for entry {}. This entry "
                    "will be ignored".format(entry.parameters.get("run_type"), entry.entry_id)
                )
                continue

            filtered_entries.append(entry)

        processed_entry_list = []
        mixing_scheme_state_data = self._generate_mixing_scheme_state_data(filtered_entries, verbose)

        for entry in entries:
            ignore_entry = False
            # get the energy adjustments
            try:
                adjustments = self.get_adjustments(entry, filtered_entries, mixing_scheme_state_data)
            except CompatibilityError as exc:
                if verbose:
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

    def _populate_df_row(self, struct_group, comp, sg, n, pd_type_1, all_entries):
        """
        helper function to populate a row of the mixing state DataFrame, given
        a list of matched structures
        """

        # consider
        # size 1 with run_type _1
        # size 1 with run_type_2
        # size 2 with only run_type_2 is impossible because of remove_non_ground_states
        # size 2 with only run_type_2
        # size 2 with one of each
        # size >2 with one 1 and multiple 2

        if len(struct_group) > 2:
            # discard all but the lowest energy entry of run_type_2
            energies = []
            for i, s in enumerate(struct_group):
                if all_entries[s.index].parameters["run_type"] in self.valid_rtypes_2:
                    energies.append((all_entries[s.index].energy_per_atom, i))

            # sort energies, lowest energy first, remove all entries except the lowest
            energies = sorted(energies)
            for tup in energies[1:]:
                del struct_group[tup[1]]

        first_entry = all_entries[struct_group[0].index]
        second_entry = all_entries[struct_group[1].index] if len(struct_group) > 1 else None

        # generate info for the DataFrame
        idx_1, idx_2 = None, None
        rt1, rt2 = None, None
        stable_1 = None
        ground_state_1, ground_state_2 = False, False

        # get the respective hull energies at this composition, if available
        hull_energy_1 = np.nan
        if pd_type_1:
            hull_energy_1 = pd_type_1.get_hull_energy_per_atom(comp)

        if first_entry.parameters["run_type"] in self.valid_rtypes_1:
            idx_1 = struct_group[0].index
            idx_2 = struct_group[1].index if len(struct_group) > 1 else None
            rt1 = first_entry.parameters["run_type"]
            rt2 = second_entry.parameters["run_type"] if second_entry else None
            # are the entries the lowest energy at this composition?
            ground_state_1 = all_entries[idx_1].energy_per_atom
            ground_state_2 = all_entries[idx_2].energy_per_atom if len(struct_group) > 1 else None
            # are they stable?
            if pd_type_1:
                stable_1 = first_entry in pd_type_1.stable_entries

        elif first_entry.parameters["run_type"] in self.valid_rtypes_2:
            idx_1 = struct_group[1].index if len(struct_group) > 1 else None
            idx_2 = struct_group[0].index
            rt2 = first_entry.parameters["run_type"]
            rt1 = second_entry.parameters["run_type"] if second_entry else None
            # are the entries the lowest energy at this composition?
            ground_state_1 = all_entries[idx_1].energy_per_atom if len(struct_group) > 1 else None
            ground_state_2 = all_entries[idx_2].energy_per_atom
            # are they stable?
            if pd_type_1:
                stable_1 = second_entry in pd_type_1.stable_entries
        else:
            raise ValueError("Problem - entries in group have strange run_type")

        return [
            comp.reduced_composition,
            sg,
            n,
            rt1,
            rt2,
            ground_state_1,
            ground_state_2,
            stable_1,
            hull_energy_1,
            None,
        ]

    def _generate_mixing_scheme_state_data(self, entries, verbose):
        """
        Generate internal state data to be passed to get_adjustments.

        Args:
            entries: The list of ComputedStructureEntry to process. This list will automatically
                be filtered based on run_type_1 and run_type_2, processed accroding to compat_1
                and compat_2, and structure matched using structure_matcher.

        Returns:
            DataFrame: A pandas DataFrame that contains information associating structures from
                different functionals with specific materials and establishing how many run_type_1
                ground states have been computed with run_type_2.
        """
        # Utility functions

        # else:
        #     print([s.index for s in struct_group])
        #     warnings.warn("Not Implemented! More than 2 entries have the same structure.")

        # Actual processing begins here

        # separate by run_type
        entries_type_1 = [e for e in entries if e.parameters["run_type"] in self.valid_rtypes_1]
        entries_type_2 = [e for e in entries if e.parameters["run_type"] in self.valid_rtypes_2]

        # Discard entries if there are fewer than 2 for a specific run_type
        # it isn't possible to implement the mixing scheme in that case
        # this code is placed here for performance reasons, to bypass all the Structure Matching
        if len(entries_type_1) < 2:
            warnings.warn(
                f"Discarding {self.run_type_1} entries because there are fewer "
                f"than 2 of them. {self.run_type_2} energies will not be modified."
            )
            return None

        if len(entries_type_2) < 2:
            warnings.warn(
                f"Discarding {self.run_type_2} entries because there are fewer "
                f"than 2 of them. {self.run_type_1} energies will not be modified."
            )
            return None

        if verbose:
            print(
                f"Processing {len(entries_type_1)} {self.run_type_1} and {len(entries_type_2)} "
                f"{self.run_type_2} entries"
            )

        # preprocess entries with any corrections
        # make an EntrySet to enable some useful methods like .chemsys and .is_ground_state
        # Make deep copies of the entries at this point, because later we'll need their original, corrected
        # energies prior to mixing functionals
        if self.compat_1:
            entries_type_1 = self.compat_1.process_entries(entries_type_1)
            if verbose:
                print(
                    f"Processed {len(entries_type_1)} compatible {self.run_type_1} entries with "
                    f"{self.compat_1.__class__.__name__}"
                )
        entries_type_1 = EntrySet(entries_type_1)

        if self.compat_2:
            entries_type_2 = self.compat_2.process_entries(entries_type_2)
            if verbose:
                print(
                    f"Processed {len(entries_type_2)} compatible {self.run_type_2} entries with "
                    f"{self.compat_2.__class__.__name__}"
                )
        entries_type_2 = EntrySet(entries_type_2)

        # Discard any entries that are not ground states from run_type_1
        # This ensures that we have at most one entry from run_type_1 at each composition
        # We must retain all run_type_2 polymorphs in case the one that matches the
        # run_type_1 ground state is not the ground state in run_type_2
        entries_type_1.remove_non_ground_states()

        # make sure both sets of entries belong to the same chemical system
        # assuming there are any gga entries at all
        if len(entries_type_1.chemsys) > 0:
            chemsys = entries_type_1.chemsys
            if not entries_type_2.chemsys <= entries_type_1.chemsys:
                warnings.warn(
                    f"{self.run_type_2} entries chemical system {entries_type_2.chemsys} is larger than "
                    f"{self.run_type_1} entries chemical system {entries_type_1.chemsys}. Entries outside the "
                    f"{self.run_type_1} chemical system will be discarded"
                )
        else:
            # if only run_type_2 entries are present, then they define the chemsys
            chemsys = entries_type_2.chemsys

        if verbose:
            print("Entries belong to the {} chemical system".format(chemsys))

        # Verify that one of the sets of entries forms a complete PhaseDiagram
        pd_type_1 = None
        try:
            pd_type_1 = PhaseDiagram(entries_type_1)
        except ValueError:
            warnings.warn(f"{self.run_type_1} entries do not form a complete PhaseDiagram.")

        try:
            PhaseDiagram(entries_type_2)
        except ValueError:
            warnings.warn(f"{self.run_type_2} entries do not form a complete PhaseDiagram.")

        # Objective: loop through all the entries, group them by structure matching (or fuzzy structure matching
        # where relevant). For each group, put a row in a pandas DataFrame with the entry_id of the run_type_1 entry,
        # the run_type_2 entry, whether or not that entry is a ground state (not necessarily on the hull), and its
        # e_above_hull
        all_entries = list(entries_type_1) + list(entries_type_2)
        row_list = []
        columns = [
            "composition",
            "spacegroup",
            "num_sites",
            "run_type_1",
            "run_type_2",
            "ground_state_energy_1",
            "ground_state_energy_2",
            "is_stable_1",
            "hull_energy_1",
            "hull_energy_2",
        ]

        def _get_sg(struc) -> int:
            """helper function to get spacegroup with a loose tolerance"""
            try:
                return struc.get_space_group_info(symprec=0.1)[1]
            except Exception:
                return -1

        # loop through all structures
        # this logic follows emmet.builders.vasp.materials.MaterialsBuilder.filter_and_group_tasks
        structures = []
        for idx, entry in enumerate(all_entries):
            s = entry.structure
            s.index: int = idx  # type: ignore
            structures.append(s)

        # First group by composition, then by spacegroup number, then by structure matching
        for comp, compgroup in groupby(sorted(structures, key=lambda s: s.composition), key=lambda s: s.composition):

            # group by spacegroup, then by number of sites (for diatmics) or by structure matching
            for sg, pregroup in groupby(sorted(compgroup, key=_get_sg), key=_get_sg):
                if comp.reduced_formula in ["O2", "H2", "Cl2", "F2", "N2"] and self.fuzzy_diatomic_matching:
                    # group by number of sites
                    for n, sitegroup in groupby(sorted(pregroup, key=lambda s: s.num_sites), key=lambda s: s.num_sites):
                        row_list.append(self._populate_df_row(list(sitegroup), comp, sg, n, pd_type_1, all_entries))
                else:
                    for group in self.structure_matcher.group_structures(list(pregroup)):
                        n = group[0].num_sites
                        # StructureMatcher.group_structures returns a list of lists,
                        # so each group should be a list containing matched structures
                        row_list.append(self._populate_df_row(group, comp, sg, n, pd_type_1, all_entries))

        mixing_state_data = pd.DataFrame(row_list, columns=columns)

        if verbose:
            # how many stable entries from run_type_1 do we have in run_type_2?
            stable_df = mixing_state_data[mixing_state_data["is_stable_1"]]
            hull_entries_2 = len(stable_df["ground_state_energy_2"].notna())
            print(
                f"Entries contain {self.run_type_2} calculations for {hull_entries_2} of {len(stable_df)} "
                f"{self.run_type_1} hull entries."
            )

        df = pd.DataFrame(row_list, columns=columns)
        # do some basic quality checks
        # assert all(rt in self.valid_rtypes_1 for rt in df[df["run_type_1"].notna()]["run_type_1"]), "Problem with the
        # DataFrame!"
        assert all(
            rt in self.valid_rtypes_2 for rt in df[df["run_type_2"].notna()]["run_type_2"]
        ), "Problem with the DataFrame!"
        # Any run_type_1 entry that is stable should have a ground state energy
        assert all(df[df["is_stable_1"]]["ground_state_energy_1"].notna()), "Problem with the DataFrame!"
        # If all run_type_1 stable entries are also present in run_type_2, the hull_energy_2 column should be populated
        if all(df[df["is_stable_1"]]["ground_state_energy_2"].notna()):
            assert all(df["hull_energy_2"].notna())

        return df

    def get_adjustments(self, entry, entries: List = None, mixing_scheme_state_data: pd.DataFrame = None):
        """
        Returns the corrections applied to a particular entry. Note that get_adjustments is not
        intended to be called directly in the SCAN mixing scheme. Call process_entries instead,
        and it will pass the required arguments to get_adjustments.

        Args:
            entry: A ComputedEntry object. The entry must be a member of the list of entries
                used to instantiate the class.
            entries: The list of ComputedEntry used to instantiate the class.
            mixing_scheme_state_data: A DataFrame containing information about which Entries
                correspond to the same materials, which are stable on the phase diagrams of
                the respective run_types, etc. Can be generated from a list of entries using
                MaterialsProjectDFTMixingScheme._generate_mixing_scheme_state_data.

        Returns:
            [EnergyAdjustment]: Energy adjustments to be applied to entry.

        Raises:
            CompatibilityError if the DFT mixing scheme cannot be applied to the entry.
        """
        adjustments: List[ConstantEnergyAdjustment] = []
        run_type = entry.parameters.get("run_type")

        if mixing_scheme_state_data is None:
            raise CompatibilityError(
                "You did not provide `mixing_scheme_state_data`. A DataFrame of mixing scheme state data is required."
            )

        # # determine the position of the entry in the original list
        # index = None
        # for i, e in enumerate(entries):
        #     if e == entry:
        #         index = i
        # if not index:
        #     raise CompatibilityError(
        #         f"Entry {entry.entry_id} is not a member of the list used to instantiate this Compatibility class,
        # and hence its energy cannot be adjusted."
        #     )

        if run_type not in self.valid_rtypes_1 + self.valid_rtypes_2:
            raise CompatibilityError(
                f"Invalid run type {run_type} for entry {entry.entry_id}. Must be one of "
                f"{self.valid_rtypes_1 + self.valid_rtypes_2}."
            )

        # The correction value depends on how many of the run_type_1 stable entries are present as run_type_2
        # calculations construct some views of the data to identify ground and stable states
        run_1_stable = mixing_scheme_state_data[
            mixing_scheme_state_data["is_stable_1"] is True
        ]  # pylint: disable=C0121
        run_1_ground_state = mixing_scheme_state_data[mixing_scheme_state_data["ground_state_energy_1"].notna()]

        # First case, ALL stable entries are present in both run types
        # in that case, prefer run_type_2 energies

        # this means all run_type_1 stable states are present in run_type_2
        # therefore, we don't need to correct run_type_2
        # instead, we discard run_type_1 entries that already exist in run_type_2
        # and correct any other run_type_1 energies to the run_type_2 scale
        if all(run_1_stable["ground_state_energy_2"].notna()):  # pylint: disable=R1705
            if run_type in self.valid_rtypes_2:  # pylint: disable=R1705
                # For run_type_2 entries, there is no correction
                return adjustments

            # Discard GGA entries whose structures already exist in SCAN.
            # TODO - replace the energy check with some kind of hash to identify entries
            elif run_type in self.valid_rtypes_1 and entry.energy_per_atom in list(
                run_1_stable["ground_state_energy_1"]
            ):
                # discard entry if its energy matches one in the DataFrame,
                # because it already exists in run_type_2
                raise CompatibilityError(
                    f"Discarding {run_type} entry {entry.entry_id} for {entry.composition.formula} "
                    f"whose structure already exists in {self.run_type_2}"
                )

            # If a GGA is not present in SCAN, correct its energy to give the same
            # e_above_hull on the SCAN hull that it would have on the GGA hull
            df_slice = mixing_scheme_state_data[
                mixing_scheme_state_data["composition"] == entry.composition.reduced_composition
            ]
            hull_energy_1 = df_slice["hull_energy_1"].iloc[0]
            hull_energy_2 = df_slice["hull_energy_2"].iloc[0]

            if not hull_energy_2:
                raise CompatibilityError(
                    f"Cannot adjust the energy of {self.run_type_1} entry {entry.entry_id} because the "
                    f"{self.run_type_2} entries given do not form a complete phase diagram. "
                    "Discarding this entry."
                )

            correction = (hull_energy_2 - hull_energy_1) * entry.composition.num_atoms

            adjustments.append(
                ConstantEnergyAdjustment(
                    correction,
                    0.0,
                    name=f"MP {self.run_type_1}/{self.run_type_2} mixing adjustment",
                    cls=self.as_dict(),
                    description=f"Place {self.run_type_1} energy onto the {self.run_type_2} hull",
                )
            )
            return adjustments

        # Second case, there are no run_type_2 energies available for any run_type_1
        # ground states. There's no way to use the run_type_2 energies in this case.
        elif sum(run_1_ground_state["ground_state_energy_2"].notna()) == 0:
            # Discard SCAN entries if there are no SCAN reference states
            if run_type in self.valid_rtypes_2:
                raise CompatibilityError(
                    f"Discarding {self.run_type_2} entry {entry.entry_id} for {entry.composition.formula} "
                    f"because there are no {self.run_type_2} reference structures"
                )
            return adjustments

        # Third case, there are run_type_2 energies available for at least some run_type_1
        # reference states. Here, we can correct run_type_2 energies onto the run_type_1 scale
        # by replacing run_type_1 e_above_hull with a run_type_2 e_above_hull for unstable polymorphs
        else:
            if run_type in self.valid_rtypes_1:  # pylint: disable=R1705
                # For run_type_1 entries, there is no correction
                return adjustments

            elif run_type in self.valid_rtypes_2 and any(
                mixing_scheme_state_data[
                    mixing_scheme_state_data["composition"] == entry.composition.reduced_composition
                ]["ground_state_energy_2"].notna()
            ):
                # We have a run_type_2 reference energy for this composition, so we can correct

                # TODO - replace the energy check with some kind of hash to identify entries
                if entry.energy_per_atom in list(run_1_ground_state["ground_state_energy_2"]):
                    # If this run_type_2 entry IS the reference state, discard it
                    raise CompatibilityError(
                        f"Discarding {run_type} entry {entry.entry_id} for {entry.composition.formula} "
                        f"because it is a {self.run_type_1} ground state."
                    )

                # e_above_hull on the run_type_2 hull that it would have on the run_type_1 hull
                df_slice = mixing_scheme_state_data[
                    mixing_scheme_state_data["composition"] == entry.composition.reduced_composition
                ]
                e_above_hull = entry.energy_per_atom - df_slice["ground_state_energy_2"].iloc[0]
                hull_energy_1 = df_slice["hull_energy_1"].iloc[0]
                # preserve the e_above_hull between the two run types
                correction = (hull_energy_1 + e_above_hull - entry.energy_per_atom) * entry.composition.num_atoms

                adjustments.append(
                    ConstantEnergyAdjustment(
                        correction,
                        0.0,
                        name=f"MP {self.run_type_1}/{self.run_type_2} mixing adjustment",
                        cls=self.as_dict(),
                        description=f"Place {self.run_type_2} energy onto the {self.run_type_1} hull",
                    )
                )
                return adjustments

            else:
                # there is no reference energy available at this composition. Discard.
                raise CompatibilityError(
                    f"Discarding {run_type} entry {entry.entry_id} for {entry.composition.formula} "
                    f"because there is no {self.run_type_2} reference energy available."
                )


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

    def get_adjustments(self, entry):
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
