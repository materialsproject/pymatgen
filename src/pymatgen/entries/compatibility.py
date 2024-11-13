"""This module implements Compatibility corrections for mixing runs of different
functionals.
"""

from __future__ import annotations

import abc
import copy
import os
import warnings
from collections import defaultdict
from typing import TYPE_CHECKING, TypeAlias, cast

import numpy as np
from joblib import Parallel, delayed
from monty.design_patterns import cached_class
from monty.json import MSONable
from monty.serialization import loadfn
from tqdm import tqdm
from uncertainties import ufloat

from pymatgen.analysis.structure_analyzer import oxide_type, sulfide_type
from pymatgen.core import SETTINGS, Composition, Element
from pymatgen.entries.computed_entries import (
    CompositionEnergyAdjustment,
    ComputedEntry,
    ComputedStructureEntry,
    ConstantEnergyAdjustment,
    EnergyAdjustment,
    TemperatureEnergyAdjustment,
)
from pymatgen.io.vasp.sets import MITRelaxSet, MPRelaxSet, VaspInputSet
from pymatgen.util.due import Doi, due
from pymatgen.util.joblib import set_python_warnings, tqdm_joblib

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any, Literal

    from pymatgen.util.typing import CompositionLike, PathLike


__author__ = "Amanda Wang, Ryan Kingsbury, Shyue Ping Ong, Anubhav Jain, Stephen Dacek, Sai Jayaraman"
__copyright__ = "Copyright 2012-2020, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "April 2020"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
MU_H2O = -2.4583  # Free energy of formation of water, eV/H2O, used by MaterialsProjectAqueousCompatibility
MP2020_COMPAT_CONFIG = loadfn(f"{MODULE_DIR}/MP2020Compatibility.yaml")
MP_COMPAT_CONFIG = loadfn(f"{MODULE_DIR}/MPCompatibility.yaml")

# This was compiled by cross-referencing structures in Materials Project from exp_compounds.json.gz
# used in the fitting of the MP2020 correction scheme, and applying the BVAnalyzer algorithm to
# determine oxidation state. O and S are not included since these are treated separately.
MP2020_ANION_OXIDATION_STATE_RANGES: dict[str, tuple[int, int]] = {
    "Br": (-1, -1),
    "Cl": (-1, -1),
    "F": (-1, -1),
    "H": (-1, -1),
    "I": (-1, -1),
    "N": (-3, -2),
    "Sb": (-3, -2),
    "Se": (-2, -1),
    "Si": (-4, -1),
    "Te": (-2, -1),
}

# Ping @janosh @rkingsbury on GitHub if this fails
if (
    MP2020_COMPAT_CONFIG["Corrections"]["GGAUMixingCorrections"]["O"]
    != MP2020_COMPAT_CONFIG["Corrections"]["GGAUMixingCorrections"]["F"]
):
    raise RuntimeError("MP2020Compatibility.yaml expected to have the same Hubbard U corrections for O and F")

AnyComputedEntry: TypeAlias = ComputedEntry | ComputedStructureEntry


class CompatibilityError(Exception):
    """Exception class for Compatibility. Raised by attempting correction
    on incompatible calculation.
    """


class Correction(abc.ABC):
    """A Correction class is a pre-defined scheme for correction a computed
    entry based on the type and chemistry of the structure and the
    calculation parameters. All Correction classes must implement a
    correct_entry method.
    """

    @abc.abstractmethod
    def get_correction(self, entry: AnyComputedEntry) -> EnergyAdjustment:
        """Get correction and uncertainty for a single entry.

        Args:
            entry: A ComputedEntry object.

        Returns:
            The energy correction to be applied and the uncertainty of the correction.

        Raises:
            CompatibilityError if entry is not compatible.
        """
        raise NotImplementedError

    def correct_entry(self, entry):
        """Corrects a single entry.

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

        # if there are no error values available for the corrections applied,
        # set correction uncertainty to not a number
        uncertainty = np.nan if updated_corr.nominal_value != 0 and updated_corr.std_dev == 0 else updated_corr.std_dev

        entry.energy_adjustments.append(ConstantEnergyAdjustment(updated_corr.nominal_value, uncertainty))

        return entry


class PotcarCorrection(Correction):
    """Check that POTCARs are valid within a pre-defined input set. This
    ensures that calculations performed using different InputSets are not
    compared against each other.

    Entry.parameters must contain a "potcar_symbols" key that is a list of
    all POTCARs used in the run. Again, using the example of an Fe2O3 run
    using Materials Project parameters, this would look like
    entry.parameters["potcar_symbols"] = ['PAW_PBE Fe_pv 06Sep2000',
    'PAW_PBE O 08Apr2002'].
    """

    def __init__(
        self,
        input_set: type[VaspInputSet],
        check_potcar: bool = True,
        check_hash: bool = False,
    ) -> None:
        """
        Args:
            input_set (InputSet): object used to generate the runs (used to check
                for correct potcar symbols).
            check_potcar (bool): If False, bypass the POTCAR check altogether. Defaults to True.
                Can also be disabled globally by running `pmg config --add PMG_POTCAR_CHECKS false`.
            check_hash (bool): If True, uses the potcar hash to check for valid
                potcars. If false, uses the potcar symbol (less reliable). Defaults to False.

        Raises:
            ValueError: if check_potcar=True and entry does not contain "potcar_symbols" key.
        """
        potcar_settings = input_set.CONFIG["POTCAR"]
        if isinstance(list(potcar_settings.values())[-1], dict):
            self.valid_potcars = {
                key: dct.get("hash" if check_hash else "symbol") for key, dct in potcar_settings.items()
            }
        else:
            if check_hash:
                raise ValueError("Cannot check hashes of potcars, since hashes are not included in the entry.")
            self.valid_potcars = potcar_settings

        self.input_set = input_set
        self.check_hash = check_hash
        self.check_potcar = check_potcar

    def __str__(self) -> str:
        return f"{self.input_set.__name__} Potcar Correction"

    def get_correction(self, entry: AnyComputedEntry) -> ufloat:
        """
        Args:
            entry (AnyComputedEntry): ComputedEntry or ComputedStructureEntry.

        Raises:
            ValueError: If entry does not contain "potcar_symbols" key.
            CompatibilityError: If entry has wrong potcar hash/symbols.

        Returns:
            ufloat: 0.0 +/- 0.0 (from uncertainties package)
        """
        if SETTINGS.get("PMG_POTCAR_CHECKS") is False or not self.check_potcar:
            return ufloat(0.0, 0.0)

        potcar_spec = entry.parameters.get("potcar_spec")
        if self.check_hash:
            if potcar_spec:
                psp_settings = {dct.get("hash") for dct in potcar_spec if dct}
            else:
                raise ValueError("Cannot check hash without potcar_spec field")
        elif potcar_spec:
            psp_settings = {dct.get("titel").split()[1] for dct in potcar_spec if dct}
        else:
            psp_settings = {sym.split()[1] for sym in entry.parameters["potcar_symbols"] if sym}

        expected_psp = {self.valid_potcars.get(el.symbol) for el in entry.elements}
        if expected_psp != psp_settings:
            raise CompatibilityError(f"Incompatible POTCAR {psp_settings}, expected {expected_psp}")
        return ufloat(0.0, 0.0)


@cached_class
class GasCorrection(Correction):
    """Correct gas energies to obtain the right formation energies. Note that
    this depends on calculations being run within the same input set.
    Used by legacy MaterialsProjectCompatibility and MITCompatibility.
    """

    def __init__(self, config_file: PathLike) -> None:
        """
        Args:
            config_file (PathLike): Path to the selected compatibility.yaml config file.
        """
        config = loadfn(config_file)
        self.name = config["Name"]
        self.cpd_energies = config["Advanced"]["CompoundEnergies"]

    def __str__(self) -> str:
        return f"{self.name} Gas Correction"

    def get_correction(self, entry: ComputedEntry | ComputedStructureEntry) -> ufloat:
        """
        Args:
            entry: A ComputedEntry/ComputedStructureEntry.

        Returns:
            ufloat: Correction.
        """
        comp = entry.composition
        # Set error to 0 because old MPCompatibility doesn't have errors
        correction = ufloat(0.0, 0.0)

        rform = entry.reduced_formula
        if rform in self.cpd_energies:
            correction += self.cpd_energies[rform] * comp.num_atoms - entry.uncorrected_energy

        return correction


@cached_class
class AnionCorrection(Correction):
    """Correct anion energies to obtain the right formation energies. Note that
    this depends on calculations being run within the same input set.

    Used by legacy MaterialsProjectCompatibility and MITCompatibility.
    """

    def __init__(self, config_file: PathLike, correct_peroxide: bool = True) -> None:
        """
        Args:
            config_file (PathLike): Path to the selected compatibility.yaml config file.
            correct_peroxide (bool): Whether peroxide/superoxide/ozonide
                corrections are to be applied or not.
        """
        config = loadfn(config_file)
        self.oxide_correction = config["OxideCorrections"]
        self.sulfide_correction = config.get("SulfideCorrections", defaultdict(float))
        self.name = config["Name"]
        self.correct_peroxide = correct_peroxide

    def __str__(self) -> str:
        return f"{self.name} Anion Correction"

    def get_correction(self, entry: ComputedEntry | ComputedStructureEntry) -> ufloat:
        """
        Args:
            entry: A ComputedEntry/ComputedStructureEntry.

        Returns:
            ufloat: Correction.
        """
        comp = entry.composition
        if len(comp) == 1:  # Skip element entry
            return ufloat(0.0, 0.0)

        correction = ufloat(0.0, 0.0)

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
                    ox_type, n_bonds = cast(
                        tuple[str, int],
                        oxide_type(entry.structure, 1.05, return_nbonds=True),
                    )
                    if ox_type in self.oxide_correction:
                        correction += self.oxide_correction[ox_type] * n_bonds
                    elif ox_type == "hydroxide":
                        correction += self.oxide_correction["oxide"] * comp["O"]
                else:
                    warnings.warn(
                        "No structure or oxide_type parameter present. Note that peroxide/superoxide corrections "
                        "are not as reliable and relies only on detection of special formulas, e.g. Li2O2."
                    )
                    rform = entry.reduced_formula
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


@cached_class
class AqueousCorrection(Correction):
    """This class implements aqueous phase compound corrections for elements
    and H2O.

    Used only by MITAqueousCompatibility.
    """

    def __init__(self, config_file: PathLike, error_file: PathLike | None = None) -> None:
        """
        Args:
            config_file (PathLike): Path to the selected compatibility.yaml config file.
            error_file (PathLike): Path to the selected compatibilityErrors.yaml config file.
        """
        config = loadfn(config_file)
        self.cpd_energies = config["AqueousCompoundEnergies"]
        # there will either be a CompositionCorrections OR an OxideCorrections key,
        # but not both, depending on the compatibility scheme we are using.
        # MITCompatibility only uses OxideCorrections, and hence self.comp_correction is none.
        self.comp_correction = config.get("CompositionCorrections", defaultdict(float))
        self.oxide_correction = config.get("OxideCorrections", defaultdict(float))
        self.name = config["Name"]
        if error_file:
            e = loadfn(error_file)
            self.cpd_errors = e.get("AqueousCompoundEnergies", defaultdict(float))
        else:
            self.cpd_errors = defaultdict(float)

    def __str__(self) -> str:
        return f"{self.name} Aqueous Correction"

    def get_correction(self, entry: ComputedEntry | ComputedStructureEntry) -> ufloat:
        """
        Args:
            entry: A ComputedEntry/ComputedStructureEntry.

        Returns:
            ufloat: Correction with uncertainty.
        """
        comp = entry.composition
        rform = comp.reduced_formula
        cpd_energies = self.cpd_energies

        correction = ufloat(0.0, 0.0)

        if rform in cpd_energies:
            if rform in {"H2", "H2O"}:
                corr = cpd_energies[rform] * comp.num_atoms - entry.uncorrected_energy - entry.correction
            else:
                corr = cpd_energies[rform] * comp.num_atoms

            err = self.cpd_errors[rform] * comp.num_atoms
            correction += ufloat(corr, err)

        if rform != "H2O":
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


@cached_class
class UCorrection(Correction):
    """This class implements the GGA/GGA+U mixing scheme, which allows mixing of
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

    common_peroxides = (
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
    )
    common_superoxides = ("LiO2", "NaO2", "KO2", "RbO2", "CsO2")
    ozonides = ("LiO3", "NaO3", "KO3", "NaO5")

    def __init__(
        self,
        config_file: PathLike,
        input_set,
        compat_type: Literal["GGA", "Advanced"],
        error_file: PathLike | None = None,
    ) -> None:
        """
        Args:
            config_file (PathLike): Path to the selected compatibility.yaml config file.
            input_set: InputSet object to check for the +U settings.
            compat_type ("GGA" | "Advanced"): "GGA" means all GGA+U
                entries are excluded. "Advanced" means mixing scheme is
                implemented to make entries compatible with each other,
                but entries which are supposed to be done in GGA+U will have the
                equivalent GGA entries excluded. For example, Fe oxides should
                have a U value under the Advanced scheme. A GGA Fe oxide run
                will therefore be excluded under the scheme.
            error_file (PathLike): Path to the selected compatibilityErrors.yaml config file.
        """
        if compat_type not in {"GGA", "Advanced"}:
            raise CompatibilityError(f"Invalid {compat_type=}")

        config = loadfn(config_file)

        self.input_set = input_set
        if compat_type == "Advanced":
            self.u_settings = self.input_set.CONFIG["INCAR"]["LDAUU"]
            self.u_corrections = config["Advanced"]["UCorrections"]
        else:
            self.u_settings = {}
            self.u_corrections = {}

        self.name = config["Name"]
        self.compat_type = compat_type

        if error_file:
            err = loadfn(error_file)
            self.u_errors = err["Advanced"]["UCorrections"]
        else:
            self.u_errors = {}

    def __str__(self) -> str:
        return f"{self.name} {self.compat_type} Correction"

    def get_correction(self, entry: ComputedEntry | ComputedStructureEntry) -> ufloat:
        """
        Args:
            entry: A ComputedEntry/ComputedStructureEntry.

        Returns:
            ufloat: Correction with Uncertainty.
        """
        calc_u = entry.parameters.get("hubbards") or defaultdict(int)
        comp = entry.composition

        elements = sorted((el for el in comp.elements if comp[el] > 0), key=lambda el: el.X)
        most_electro_neg = elements[-1].symbol
        correction = ufloat(0.0, 0.0)

        u_corr = self.u_corrections.get(most_electro_neg, {})
        u_settings = self.u_settings.get(most_electro_neg, {})
        u_errors = self.u_errors.get(most_electro_neg, defaultdict(float))

        for el in comp.elements:
            sym = el.symbol
            # Check for bad U values
            if calc_u.get(sym, 0) != u_settings.get(sym, 0):
                raise CompatibilityError(f"Invalid U value of {calc_u.get(sym, 0)} on {sym}")
            if sym in u_corr:
                correction += ufloat(u_corr[sym], u_errors[sym]) * comp[el]

        return correction


class Compatibility(MSONable, abc.ABC):
    """Abstract Compatibility class, not intended for direct use.
    Compatibility classes are used to correct the energies of an entry or a set
    of entries. All Compatibility classes must implement get_adjustments() method.
    """

    @abc.abstractmethod
    def get_adjustments(self, entry: AnyComputedEntry) -> list[EnergyAdjustment]:
        """Get the energy adjustments for a ComputedEntry.

        This method must generate a list of EnergyAdjustment objects
        of the appropriate type (constant, composition-based, or temperature-based)
        to be applied to the ComputedEntry, and must raise a CompatibilityError
        if the entry is not compatible.

        Args:
            entry: A ComputedEntry object.

        Returns:
            list[EnergyAdjustment]: A list of EnergyAdjustment to be applied to the
                Entry.

        Raises:
            CompatibilityError if the entry is not compatible
        """
        raise NotImplementedError

    def process_entry(
        self,
        entry: ComputedEntry,
        inplace: bool = True,
        **kwargs,
    ) -> ComputedEntry | None:
        """Process a single entry with the chosen Corrections.
        Note that this method may change the original entry.

        Args:
            entry (ComputedEntry): A ComputedEntry object.
            inplace (bool): Whether to adjust the entry in place. Defaults to True.
            **kwargs: Will be passed to process_entries().

        Returns:
            An adjusted entry if entry is compatible, else None.
        """
        if not inplace:
            entry = copy.deepcopy(entry)

        _entry: tuple[ComputedEntry, bool] | None = self._process_entry_inplace(entry, **kwargs)

        return _entry[0] if _entry is not None else None

    def _process_entry_inplace(
        self,
        entry: AnyComputedEntry,
        clean: bool = True,
        on_error: Literal["ignore", "warn", "raise"] = "ignore",
    ) -> tuple[ComputedEntry, bool] | None:
        """Process a single entry with the chosen Corrections.
        Note that this method will change the original entry.

        Args:
            entry (AnyComputedEntry): An AnyComputedEntry object.
            clean (bool): Whether to remove any previously-applied energy adjustments.
                If True, all EnergyAdjustment are removed prior to processing the Entry.
                Defaults to True.
            on_error ("ignore" | "warn" | "raise"): What to do when get_adjustments(entry)
                raises CompatibilityError. Defaults to "ignore".

        Returns:
            tuple[AnyComputedEntry, ignore_entry (bool)] if entry is compatible, else None.
        """
        ignore_entry: bool = False
        # If clean, remove all previous adjustments from the entry
        if clean:
            entry.energy_adjustments = []

        # Get the energy adjustments
        try:
            adjustments: list[EnergyAdjustment] = self.get_adjustments(entry)

        except CompatibilityError as exc:
            if on_error == "raise":
                raise
            if on_error == "warn":
                warnings.warn(str(exc))
            return None

        for e_adj in adjustments:
            # Check if this correction already been applied
            if (e_adj.name, e_adj.cls, e_adj.value) in [
                (ea2.name, ea2.cls, ea2.value) for ea2 in entry.energy_adjustments
            ]:
                # We already applied this exact correction. Do nothing.
                pass

            elif (e_adj.name, e_adj.cls) in [(ea2.name, ea2.cls) for ea2 in entry.energy_adjustments]:
                # We already applied a correction with the same name
                # but a different value. Something is wrong.
                ignore_entry = True
                warnings.warn(
                    f"Entry {entry.entry_id} already has an energy adjustment called {e_adj.name}, but its "
                    f"value differs from the value of {e_adj.value:.3f} calculated here. This "
                    "Entry will be discarded."
                )

            else:
                # Add the correction to the energy_adjustments list
                entry.energy_adjustments.append(e_adj)

        return entry, ignore_entry

    def process_entries(
        self,
        entries: AnyComputedEntry | list[AnyComputedEntry],
        clean: bool = True,
        verbose: bool = False,
        inplace: bool = True,
        n_workers: int = 1,
        on_error: Literal["ignore", "warn", "raise"] = "ignore",
    ) -> list[AnyComputedEntry]:
        """Process a sequence of entries with the chosen Compatibility scheme.

        Warning: This method changes entries in place! All changes can be undone and original entries
        restored by setting entry.energy_adjustments = [].

        Args:
            entries (AnyComputedEntry | list[AnyComputedEntry]): A sequence of
                Computed(Structure)Entry objects.
            clean (bool): Whether to remove any previously-applied energy adjustments.
                If True, all EnergyAdjustment are removed prior to processing the Entry.
                Defaults to True.
            verbose (bool): Whether to display progress bar for processing multiple entries.
                Defaults to False.
            inplace (bool): Whether to adjust input entries in place. Defaults to True.
            n_workers (int): Number of workers to use for parallel processing. Defaults to 1.
            on_error ('ignore' | 'warn' | 'raise'): What to do when get_adjustments(entry)
                raises CompatibilityError. Defaults to 'ignore'.

        Returns:
            list[AnyComputedEntry]: Adjusted entries. Entries in the original list incompatible with
                chosen correction scheme are excluded from the returned list.
        """
        # if single entry convert to list
        if isinstance(entries, ComputedEntry):  # True for ComputedStructureEntry too
            entries = [entries]

        processed_entry_list: list[AnyComputedEntry] = []

        # if inplace = False, process entries on a copy
        if not inplace:
            entries = copy.deepcopy(entries)

        if n_workers == 1:
            for entry in tqdm(entries, disable=not verbose):
                result = self._process_entry_inplace(entry, clean, on_error)
                if result is None:
                    continue
                entry, ignore_entry = result
                if not ignore_entry:
                    processed_entry_list.append(entry)
        elif not inplace:
            # set python warnings to ignore otherwise warnings will be printed multiple times
            with (
                tqdm_joblib(tqdm(total=len(entries), disable=not verbose)),
                set_python_warnings("ignore"),
            ):
                results = Parallel(n_jobs=n_workers)(
                    delayed(self._process_entry_inplace)(entry, clean, on_error) for entry in entries
                )
            for result in results:
                if result is None:
                    continue
                entry, ignore_entry = result
                if not ignore_entry:
                    processed_entry_list.append(entry)
        else:
            raise ValueError("Parallel processing is not possible with for 'inplace=True'")

        return processed_entry_list

    @staticmethod
    def explain(entry: ComputedEntry) -> None:
        """Print an explanation of the energy adjustments applied by the
        Compatibility class. Inspired by the "explain" methods in many database
        methodologies.

        Args:
            entry: A ComputedEntry.
        """
        print(
            f"The uncorrected energy of {entry.composition} is {entry.uncorrected_energy:.3f} eV "
            f"({entry.uncorrected_energy / entry.composition.num_atoms:.3f} eV/atom)."
        )

        if len(entry.energy_adjustments) > 0:
            print("The following energy adjustments have been applied to this entry:")
            for adj in entry.energy_adjustments:
                print(f"\t\t{adj.name}: {adj.value:.3f} eV ({adj.value / entry.composition.num_atoms:.3f} eV/atom)")
        elif entry.correction == 0:
            print("No energy adjustments have been applied to this entry.")

        print(f"The final energy after adjustments is {entry.energy:.3f} eV ({entry.energy_per_atom:.3f} eV/atom).")


class CorrectionsList(Compatibility):
    """The CorrectionsList class combines a list of corrections to be applied to
    an entry or a set of entries. Note that some of the Corrections have
    interdependencies. For example, PotcarCorrection must always be used
    before any other compatibility. Also, AnionCorrection("MP") must be used
    with PotcarCorrection("MP") (similarly with "MIT"). Typically,
    you should use the specific MaterialsProjectCompatibility and
    MITCompatibility subclasses instead.
    """

    def __init__(
        self,
        corrections: Sequence[Correction],
        run_types: list[Literal["GGA", "GGA+U", "PBE", "PBE+U"]] | None = None,
    ) -> None:
        """
        Args:
            corrections (Sequence[Correction]): Corrections to apply.
            run_types (list[str]): Valid DFT run types for this correction scheme.
                Entries with run type other than those in this list will be excluded
                from the list returned by process_entries. The default value captures
                both GGA and GGA+U run types historically used by the Materials Project.
        """
        if run_types is None:
            run_types = ["GGA", "GGA+U", "PBE", "PBE+U"]
        self.corrections = corrections
        self.run_types = run_types
        super().__init__()

    def get_adjustments(self, entry: AnyComputedEntry) -> list[EnergyAdjustment]:
        """Get the list of energy adjustments to be applied to an entry."""
        adjustment_list = []
        if entry.parameters.get("run_type") not in self.run_types:
            raise CompatibilityError(
                f"Entry {entry.entry_id} has invalid run type {entry.parameters.get('run_type')}. "
                f"Must be GGA or GGA+U. Discarding."
            )

        corrections, uncertainties = self.get_corrections_dict(entry)

        for k, v in corrections.items():
            uncertainty = np.nan if v != 0 and uncertainties[k] == 0 else uncertainties[k]
            adjustment_list.append(ConstantEnergyAdjustment(v, uncertainty=uncertainty, name=k, cls=self.as_dict()))

        return adjustment_list

    def get_corrections_dict(
        self,
        entry: AnyComputedEntry,
    ) -> tuple[dict[str, float], dict[str, float]]:
        """Get the correction values and uncertainties applied to a particular entry.

        Args:
            entry: An AnyComputedEntry object.

        Returns:
            tuple[dict[str, float], dict[str, float]]: Map from correction names to values
                (1st) and uncertainties (2nd).
        """
        corrections = {}
        uncertainties = {}
        for c in self.corrections:
            val = c.get_correction(entry)
            if val != 0:
                corrections[str(c)] = val.nominal_value
                uncertainties[str(c)] = val.std_dev
        return corrections, uncertainties

    def get_explanation_dict(self, entry: ComputedEntry) -> dict[str, Any]:
        """Explain the corrections applied for a given compatibility scheme. Inspired by the
        "explain" methods in many database methodologies.

        Args:
            entry: A ComputedEntry.

        Returns:
            dict[str, str | float | list[dict[str, Union[str, float]]]: of the form
                {"Compatibility": "string",
                "Uncorrected_energy": float,
                "Corrected_energy": float,
                "correction_uncertainty:" float,
                "Corrections": [{"Name of Correction": {
                "Value": float, "Explanation": "string", "Uncertainty": float}]}
        """
        corr_entry = self.process_entry(entry)
        uncorrected_energy = (corr_entry or entry).uncorrected_energy
        corrected_energy = corr_entry.energy if corr_entry else None
        correction_uncertainty = corr_entry.correction_uncertainty if corr_entry else None

        dct = {
            "compatibility": type(self).__name__,
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

            if c.__doc__ is None:
                raise RuntimeError("__doc__ of some correction is None.")
            cd = {
                "name": str(c),
                "description": c.__doc__.split("Args")[0].strip(),
                "value": corr_dict.get(str(c), 0),
                "uncertainty": uncer,
            }
            corrections.append(cd)
        dct["corrections"] = corrections
        return dct

    def explain(self, entry: ComputedEntry) -> None:
        """Print an explanation of the corrections that are being applied for a
        given compatibility scheme. Inspired by the "explain" methods in many
        database methodologies.

        Args:
            entry: A ComputedEntry.
        """
        dct = self.get_explanation_dict(entry)
        print(f"The uncorrected value of the energy of {entry.composition} is {dct['uncorrected_energy']:f} eV")
        print(f"The following corrections / screening are applied for {dct['compatibility']}:\n")
        for corr in dct["corrections"]:
            print(f"{corr['name']} correction: {corr['description']}\n")
            print(f"For the entry, this correction has the value {corr['value']:f} eV.")
            if corr["uncertainty"] != 0 or corr["value"] == 0:
                print(f"This correction has an uncertainty value of {corr['uncertainty']:f} eV.")
            else:
                print("This correction does not have uncertainty data available")
            print("-" * 30)

        print(f"The final energy after corrections is {dct['corrected_energy']:f}")


class MaterialsProjectCompatibility(CorrectionsList):
    """This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. Note that this should only be used for VASP calculations using the
    MaterialsProject parameters (see pymatgen.io.vasp.sets.MPVaspInputSet).
    Using this compatibility scheme on runs with different parameters is not valid.
    """

    def __init__(
        self,
        compat_type: Literal["GGA", "Advanced"] = "Advanced",
        correct_peroxide: bool = True,
        check_potcar_hash: bool = False,
    ) -> None:
        """
        Args:
            compat_type ("GGA" | "Advanced"): "GGA" means all GGA+U
                entries are excluded. "Advanced" means mixing scheme is
                implemented to make entries compatible with each other,
                but entries which are supposed to be done in GGA+U will have the
                equivalent GGA entries excluded. For example, Fe oxides should
                have a U value under the Advanced scheme. A GGA Fe oxide run
                will therefore be excluded under the scheme.
            correct_peroxide: Specify whether peroxide/superoxide/ozonide
                corrections are to be applied or not.
            check_potcar_hash (bool): Use potcar hash to verify potcars are correct.
            silence_deprecation (bool): Silence deprecation warning. Defaults to False.
        """
        warnings.warn(  # added by @janosh on 2023-05-25
            "MaterialsProjectCompatibility is deprecated, Materials Project formation energies "
            "use the newer MaterialsProject2020Compatibility scheme.",
            DeprecationWarning,
        )
        self.compat_type = compat_type
        self.correct_peroxide = correct_peroxide
        self.check_potcar_hash = check_potcar_hash
        file_path = f"{MODULE_DIR}/MPCompatibility.yaml"
        super().__init__(
            [
                PotcarCorrection(MPRelaxSet, check_hash=check_potcar_hash),
                GasCorrection(file_path),
                AnionCorrection(file_path, correct_peroxide=correct_peroxide),
                UCorrection(file_path, MPRelaxSet, compat_type),
            ]
        )


# Note from Ryan Kingsbury (2022-10-14): MaterialsProject2020Compatibility inherits from Compatibility
# instead of CorrectionsList which came before it because CorrectionsList had technical limitations.
# When we did the new scheme (MP2020) we decided to refactor the base Compatibility class to not
# require CorrectionsList.

# This was particularly helpful for the AqueousCorrection class. The new system gives complete
# flexibility to process entries however needed inside the get_adjustments() method, rather than
# having to create a list of separate correction classes.


@cached_class
class MaterialsProject2020Compatibility(Compatibility):
    """This class implements the Materials Project 2020 energy correction scheme, which
    incorporates uncertainty quantification and allows for mixing of GGA and GGA+U entries
    (see References).

    Note that this scheme should only be applied to VASP calculations that use the
    Materials Project input set parameters (see pymatgen.io.vasp.sets.MPRelaxSet). Using
    this compatibility scheme on calculations with different parameters is not valid.

    The option `strict_anions` was added due to a bug. See PR #3803 (May 2024) for
    related discussion. This behavior may change in subsequent versions as a more comprehensive
    fix for this issue may be found.

    Note: While the correction scheme is largely composition-based, the energy corrections
    applied to ComputedEntry and ComputedStructureEntry can differ for O and S-containing
    structures if entry.data['oxidation_states'] is not populated or explicitly set. This
    occurs because pymatgen will use atomic distances to classify O and S anions as
    superoxide/peroxide/oxide and sulfide/polysulfide, resp. when oxidation states are not
    provided. If you want the most accurate corrections possible, supply pre-defined
    oxidation states to entry.data or pass ComputedStructureEntry.
    """

    def __init__(
        self,
        compat_type: str = "Advanced",
        correct_peroxide: bool = True,
        strict_anions: Literal["require_exact", "require_bound", "no_check"] = "require_bound",
        check_potcar: bool = True,
        check_potcar_hash: bool = False,
        config_file: str | None = None,
    ) -> None:
        """
        Args:
            compat_type: Two options, GGA or Advanced. GGA means all GGA+U
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
                populated. Default: "Advanced"
            correct_peroxide: Specify whether peroxide/superoxide/ozonide
                corrections are to be applied or not. If false, all oxygen-containing
                compounds are assigned the 'oxide' correction. Default: True
            strict_anions: only apply the anion corrections to anions. The option
                "require_exact" will only apply anion corrections in cases where the
                anion oxidation state is between the oxidation states used
                in the experimental fitting data. The option "require_bound" will
                define an anion as any species with an oxidation state value of <= -1.
                This prevents the anion correction from being applied to unrealistic
                hypothetical structures containing large proportions of very electronegative
                elements, thus artificially over-stabilizing the compound. Set to "no_check"
                to restore the original behavior described in the associated publication. Default: True
            check_potcar (bool): Check that the POTCARs used in the calculation are consistent
                with the Materials Project parameters. False bypasses this check altogether. Default: True
                Can also be disabled globally by running `pmg config --add PMG_POTCAR_CHECKS false`.
            check_potcar_hash (bool): Use potcar hash to verify POTCAR settings are
                consistent with MPRelaxSet. If False, only the POTCAR symbols will
                be used. Default: False
            config_file (Path): Path to the selected compatibility.yaml config file.
                If None, defaults to `MP2020Compatibility.yaml` distributed with
                pymatgen.

        References:
            Wang, A., Kingsbury, R., McDermott, M., Horton, M., Jain. A., Ong, S.P.,
                Dwaraknath, S., Persson, K. A framework for quantifying uncertainty
                in DFT energy corrections. Scientific Reports 11: 15496, 2021.
                https://doi.org/10.1038/s41598-021-94550-5

            Jain, A. et al. Formation enthalpies by mixing GGA and GGA + U calculations.
                Phys. Rev. B - Condens. Matter Mater. Phys. 84, 1-10 (2011).
        """
        if compat_type not in ["GGA", "Advanced"]:
            raise CompatibilityError(f"Invalid {compat_type=}")

        self.compat_type = compat_type
        self.correct_peroxide = correct_peroxide
        self.strict_anions = strict_anions
        self.check_potcar = check_potcar
        self.check_potcar_hash = check_potcar_hash

        # load corrections and uncertainties
        if config_file:
            if os.path.isfile(config_file):
                self.config_file: str | None = config_file
                config = loadfn(self.config_file)
            else:
                raise ValueError(f"Custom MaterialsProject2020Compatibility {config_file=} does not exist.")
        else:
            self.config_file = None

            config = MP2020_COMPAT_CONFIG

        self.name = config["Name"]
        self.comp_correction = config["Corrections"].get("CompositionCorrections", defaultdict(float))
        self.comp_errors = config["Uncertainties"].get("CompositionCorrections", defaultdict(float))

        if self.compat_type == "Advanced":
            self.u_settings = MPRelaxSet.CONFIG["INCAR"]["LDAUU"]
            self.u_corrections = config["Corrections"].get("GGAUMixingCorrections", defaultdict(float))
            self.u_errors = config["Uncertainties"].get("GGAUMixingCorrections", defaultdict(float))
        else:
            self.u_settings = {}
            self.u_corrections = {}
            self.u_errors = {}

    def get_adjustments(self, entry: AnyComputedEntry) -> list[EnergyAdjustment]:
        """Get the energy adjustments for a ComputedEntry or ComputedStructureEntry.

        Energy corrections are implemented directly in this method instead of in
        separate AnionCorrection, GasCorrection, or UCorrection classes which
        were used in the legacy correction scheme.

        Args:
            entry: A ComputedEntry or ComputedStructureEntry object.

        Returns:
            list[EnergyAdjustment]: A list of EnergyAdjustment to be applied to the Entry.

        Raises:
            CompatibilityError if the entry is not compatible
        """
        if entry.parameters.get("run_type") not in ("GGA", "GGA+U"):
            raise CompatibilityError(
                f"Entry {entry.entry_id} has invalid run type {entry.parameters.get('run_type')}. "
                f"Must be GGA or GGA+U. Discarding."
            )

        # check the POTCAR symbols
        # this should return ufloat(0, 0) or raise a CompatibilityError or ValueError
        if entry.parameters.get("software", "vasp") == "vasp":
            pc = PotcarCorrection(
                MPRelaxSet,
                check_hash=self.check_potcar_hash,
                check_potcar=self.check_potcar,
            )
            pc.get_correction(entry)

        # apply energy adjustments
        adjustments: list[CompositionEnergyAdjustment] = []

        comp = entry.composition
        rform = comp.reduced_formula
        # sorted list of elements, ordered by electronegativity
        sorted_elements = sorted((el for el in comp.elements if comp[el] > 0), key=lambda el: el.X)

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
                    ox_type = oxide_type(entry.structure, 1.05)
                else:
                    warnings.warn(
                        "No structure or oxide_type parameter present. Note that peroxide/superoxide corrections "
                        "are not as reliable and relies only on detection of special formulas, e.g. Li2O2."
                    )

                    common_peroxides = "Li2O2 Na2O2 K2O2 Cs2O2 Rb2O2 BeO2 MgO2 CaO2 SrO2 BaO2".split()
                    common_superoxides = "LiO2 NaO2 KO2 RbO2 CsO2".split()
                    ozonides = "LiO3 NaO3 KO3 NaO5".split()

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
                    name=f"MP2020 anion correction ({ox_type})",
                    cls=self.as_dict(),
                )
            )

        # Check for anion corrections
        # only apply anion corrections if the element is an anion
        # first check for a pre-populated oxidation states key
        # the key is expected to comprise a dict corresponding to the first element output by
        # Composition.oxi_state_guesses(), e.g. {'Al': 3.0, 'S': 2.0, 'O': -2.0} for 'Al2SO4'
        if "oxidation_states" not in entry.data:
            # try to guess the oxidation states from composition
            # for performance reasons, fail if the composition is too large
            try:
                oxi_states = entry.composition.oxi_state_guesses(max_sites=-20)
            except ValueError:
                oxi_states = ({},)

            entry.data["oxidation_states"] = (oxi_states or ({},))[0]

        if entry.data["oxidation_states"] == {}:
            warnings.warn(
                f"Failed to guess oxidation states for Entry {entry.entry_id} "
                f"({entry.reduced_formula}). Assigning anion correction to "
                "only the most electronegative atom."
            )

        for anion in ("Br", "I", "Se", "Si", "Sb", "Te", "H", "N", "F", "Cl"):
            if Element(anion) in comp and anion in self.comp_correction:
                apply_correction = False
                oxidation_state = entry.data["oxidation_states"].get(anion, 0)
                # if the oxidation_states key is not populated, only apply the correction if the anion
                # is the most electronegative element
                if oxidation_state < 0:
                    apply_correction = True
                    if self.strict_anions == "require_bound" and oxidation_state > -1:
                        # This is not an anion. Noting that the rare case of a fractional
                        # oxidation state in range [-1, 0] might be considered an anionic.
                        # This could include suboxides or metal-rich pnictides, chalcogenides etc.
                        # However! these cases are not included in the experimental fitting data
                        # used for the correction scheme, and so there is no information for
                        # whether the corrections are appropriate in this instance, and likely
                        # may.
                        apply_correction = False
                else:
                    most_electroneg = sorted_elements[-1].symbol
                    if anion == most_electroneg:
                        apply_correction = True

                if self.strict_anions == "require_exact":
                    apply_correction = False
                    if (oxi_range := MP2020_ANION_OXIDATION_STATE_RANGES.get(anion)) and (
                        oxi_range[0] <= oxidation_state <= oxi_range[1]
                    ):
                        apply_correction = True

                if apply_correction:
                    adjustments.append(
                        CompositionEnergyAdjustment(
                            self.comp_correction[anion],
                            comp[anion],
                            uncertainty_per_atom=self.comp_errors[anion],
                            name=f"MP2020 anion correction ({anion})",
                            cls=self.as_dict(),
                        )
                    )

        # GGA / GGA+U mixing scheme corrections
        calc_u = entry.parameters.get("hubbards")
        calc_u = defaultdict(int) if calc_u is None else calc_u
        most_electroneg = sorted_elements[-1].symbol
        u_corrections = self.u_corrections.get(most_electroneg, defaultdict(float))
        u_settings = self.u_settings.get(most_electroneg, defaultdict(float))
        u_errors = self.u_errors.get(most_electroneg, defaultdict(float))

        for el in comp.elements:
            symbol = el.symbol
            # Check for bad U values
            expected_u = float(u_settings.get(symbol, 0))
            actual_u = float(calc_u.get(symbol, 0))
            if actual_u != expected_u:
                raise CompatibilityError(
                    f"Invalid U value of {actual_u:.3} on {symbol}, expected {expected_u:.3} for {entry.as_dict()}"
                )
            if symbol in u_corrections:
                adjustments.append(
                    CompositionEnergyAdjustment(
                        u_corrections[symbol],
                        comp[el],
                        uncertainty_per_atom=u_errors[symbol],
                        name=f"MP2020 GGA/GGA+U mixing correction ({symbol})",
                        cls=self.as_dict(),
                    )
                )

        return adjustments


class MITCompatibility(CorrectionsList):
    """This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. Note that this should only be used for VASP calculations using the
    MIT parameters (see pymatgen.io.vasp.sets MITVaspInputSet). Using
    this compatibility scheme on runs with different parameters is not valid.
    """

    def __init__(
        self,
        compat_type: Literal["GGA", "Advanced"] = "Advanced",
        correct_peroxide: bool = True,
        check_potcar_hash: bool = False,
    ) -> None:
        """
        Args:
            compat_type: Two options, GGA or Advanced. GGA means all GGA+U
                entries are excluded. Advanced means mixing scheme is
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
        fp = f"{MODULE_DIR}/MITCompatibility.yaml"
        super().__init__(
            [
                PotcarCorrection(MITRelaxSet, check_hash=check_potcar_hash),
                GasCorrection(fp),
                AnionCorrection(fp, correct_peroxide=correct_peroxide),
                UCorrection(fp, MITRelaxSet, compat_type),
            ]
        )


class MITAqueousCompatibility(CorrectionsList):
    """This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. Note that this should only be used for VASP calculations using the
    MIT parameters (see pymatgen.io.vasp.sets MITVaspInputSet). Using
    this compatibility scheme on runs with different parameters is not valid.
    """

    def __init__(
        self,
        compat_type: Literal["GGA", "Advanced"] = "Advanced",
        correct_peroxide: bool = True,
        check_potcar_hash: bool = False,
    ) -> None:
        """
        Args:
            compat_type: Two options, GGA or Advanced. GGA means all GGA+U
                entries are excluded. Advanced means mixing scheme is
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
        fp = f"{MODULE_DIR}/MITCompatibility.yaml"
        super().__init__(
            [
                PotcarCorrection(MITRelaxSet, check_hash=check_potcar_hash),
                GasCorrection(fp),
                AnionCorrection(fp, correct_peroxide=correct_peroxide),
                UCorrection(fp, MITRelaxSet, compat_type),
                AqueousCorrection(fp),
            ]
        )


@cached_class
@due.dcite(
    Doi(
        "10.1103/PhysRevB.85.235438",
        "Pourbaix scheme to combine calculated and experimental data",
    )
)
class MaterialsProjectAqueousCompatibility(Compatibility):
    """This class implements the Aqueous energy referencing scheme for constructing
    Pourbaix diagrams from DFT energies, as described in Persson et al.

    This scheme applies various energy adjustments to convert DFT energies into
    Gibbs free energies of formation at 298 K and to guarantee that the experimental
    formation free energy of H2O is reproduced. Briefly, the steps are:

        1. Beginning with the DFT energy of O2, adjust the energy of H2 so that
           the experimental reaction energy of -2.458 eV/H2O is reproduced.
        2. Add entropy to the DFT energy of any compounds that are liquid or
           gaseous at room temperature
        3. Adjust the DFT energies of solid hydrate compounds (compounds that
           contain water, e.g. FeO.nH2O) such that the energies of the embedded
           H2O molecules are equal to the experimental free energy

    The above energy adjustments are computed dynamically based on the input
    Entries.

    References:
        K.A. Persson, B. Waldwick, P. Lazic, G. Ceder, Prediction of solid-aqueous
        equilibria: Scheme to combine first-principles calculations of solids with
        experimental aqueous states, Phys. Rev. B - Condens. Matter Mater. Phys.
        85 (2012) 1-12. doi:10.1103/PhysRevB.85.235438.
    """

    def __init__(
        self,
        solid_compat: (Compatibility | type[Compatibility] | None) = MaterialsProject2020Compatibility,
        o2_energy: float | None = None,
        h2o_energy: float | None = None,
        h2o_adjustments: float | None = None,
    ) -> None:
        """Initialize the MaterialsProjectAqueousCompatibility class.

        Note that this class requires as inputs the ground-state DFT energies of O2 and H2O, plus the value of any
        energy adjustments applied to an H2O molecule. If these parameters are not provided in __init__, they can
        be automatically populated by including ComputedEntry for the ground state of O2 and H2O in a list of
        entries passed to process_entries. process_entries will fail if one or the other is not provided.

        Args:
            solid_compat: Compatibility scheme used to pre-process solid DFT energies prior to applying aqueous
                energy adjustments. May be passed as a class (e.g. MaterialsProject2020Compatibility) or an instance
                (e.g., MaterialsProject2020Compatibility()). If None, solid DFT energies are used as-is.
                Default: MaterialsProject2020Compatibility
            o2_energy: The ground-state DFT energy of oxygen gas, including any adjustments or corrections, in eV/atom.
                If not set, this value will be determined from any O2 entries passed to process_entries.
                Default: None
            h2o_energy: The ground-state DFT energy of water, including any adjustments or corrections, in eV/atom.
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
            raise ValueError("Expected a Compatibility class, instance of a Compatibility or None")

        self.o2_energy = o2_energy
        self.h2o_energy = h2o_energy
        self.h2_energy = None
        self.h2o_adjustments = h2o_adjustments

        if not all([self.o2_energy, self.h2o_energy, self.h2o_adjustments]):
            warnings.warn(
                f"You did not provide the required O2 and H2O energies. {type(self).__name__} "
                "needs these energies in order to compute the appropriate energy adjustments. It will try "
                "to determine the values from ComputedEntry for O2 and H2O passed to process_entries, but "
                "will fail if these entries are not provided."
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

    def get_adjustments(self, entry: ComputedEntry) -> list[EnergyAdjustment]:
        """Get the corrections applied to a particular entry.

        Args:
            entry: A ComputedEntry object.

        Returns:
            list[EnergyAdjustment]: Energy adjustments to be applied to entry.

        Raises:
            CompatibilityError if the required O2 and H2O energies have not been provided to
            MaterialsProjectAqueousCompatibility during init or in the list of entries passed to process_entries.
        """
        adjustments = []
        if self.o2_energy is None or self.h2o_energy is None or self.h2o_adjustments is None:
            raise CompatibilityError(
                "You did not provide the required O2 and H2O energies. "
                f"{type(self).__name__} needs these energies in order to compute "
                "the appropriate energy adjustments. Either specify the energies as arguments "
                f"to {type(self).__name__}.__init__ or run process_entries on a list that includes ComputedEntry for "
                "the ground state of O2 and H2O."
            )

        # compute the free energies of H2 and H2O (eV/atom) to guarantee that the
        # formation-free energy of H2O is equal to -2.4583 eV/H2O from experiments
        # (MU_H2O from Pourbaix module)

        # Free energy of H2 in eV/atom, fitted using Eq. 40 of Persson et al. PRB 2012 85(23)
        # https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.235438
        # for this calculation ONLY, we need the (corrected) DFT energy of water
        self.fit_h2_energy = round(
            0.5
            * (
                3 * (self.h2o_energy - self.cpd_entropies["H2O"]) - (self.o2_energy - self.cpd_entropies["O2"]) - MU_H2O
            ),
            6,
        )

        comp = entry.composition
        rform = comp.reduced_formula

        # use fit_h2_energy to adjust the energy of all H2 polymorphs such that
        # the lowest energy polymorph has the correct experimental value
        # if H2O and O2 energies have been set explicitly via kwargs, then
        # all H2 polymorphs will get the same energy.
        if rform == "H2":
            if self.h2_energy is None:
                raise ValueError("H2 energy not set")
            adjustments.append(
                ConstantEnergyAdjustment(
                    (self.fit_h2_energy - self.h2_energy) * comp.num_atoms,
                    uncertainty=np.nan,
                    name="MP Aqueous H2 / H2O referencing",
                    cls=self.as_dict(),
                    description="Adjusts the H2 energy to reproduce the experimental "
                    "Gibbs formation free energy of H2O, based on the DFT energy "
                    "of Oxygen and H2O",
                )
            )

        # add minus T delta S to the DFT energy (enthalpy) of compounds that are
        # molecular-like at room temperature
        if rform in self.cpd_entropies:
            adjustments.append(
                TemperatureEnergyAdjustment(
                    -self.cpd_entropies[rform] / 298,
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

        # For any compound except water, check if it is a hydrate (contains
        # H2O in its structure). If so, adjust the energy to remove MU_H2O eV per
        # embedded water molecule.
        # in other words, we assume that the DFT energy of such a compound is really
        # a superposition of the "real" solid DFT energy (FeO in this case) and the free
        # energy of some water molecules
        # e.g. that E_FeO.nH2O = E_FeO + n * g_H2O
        # so, to get the most accurate Gibbs free energy, we want to replace
        # g_FeO.nH2O = E_FeO.nH2O + dE_Fe + (n+1) * dE_O + 2n dE_H
        # with
        # g_FeO = E_FeO.nH2O + dE_Fe + dE_O + n g_H2O
        # where E is DFT energy, dE is an energy correction, and g is Gibbs free energy
        # of formation
        # This means we have to 1) reverse any energy corrections that have already been
        # applied to H and O in water and then 2) remove the free energy of the water
        # molecules from the hydrated solid energy.
        if rform != "H2O":
            # count the number of whole water molecules in the composition
            rcomp, factor = comp.get_reduced_composition_and_factor()
            nH2O = int(min(rcomp["H"] / 2.0, rcomp["O"])) * factor
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

    def process_entries(
        self,
        entries: list[AnyComputedEntry],
        clean: bool = False,
        verbose: bool = False,
        inplace: bool = True,
        n_workers: int = 1,
        on_error: Literal["ignore", "warn", "raise"] = "ignore",
    ) -> list[AnyComputedEntry]:
        """Process a sequence of entries with the chosen Compatibility scheme.

        Args:
            entries (list[ComputedEntry | ComputedStructureEntry]): Entries to be processed.
            clean (bool): Whether to remove any previously-applied energy adjustments.
                If True, all EnergyAdjustment are removed prior to processing the Entry.
                Default is False.
            verbose (bool): Whether to display progress bar for processing multiple entries.
                Default is False.
            inplace (bool): Whether to modify the entries in place. If False, a copy of the
                entries is made and processed. Default is True.
            n_workers (int): Number of workers to use for parallel processing. Default is 1.
            on_error ('ignore' | 'warn' | 'raise'): What to do when get_adjustments(entry)
                raises CompatibilityError. Defaults to 'ignore'.

        Returns:
            list[AnyComputedEntry]: Adjusted entries. Entries in the original list incompatible with
                chosen correction scheme are excluded from the returned list.
        """
        # Convert input arg to a list if not already
        if isinstance(entries, ComputedEntry):
            entries = [entries]

        # If not inplace, process entries on a copy
        if not inplace:
            entries = copy.deepcopy(entries)

        # Pre-process entries with the given solid compatibility class
        if self.solid_compat:
            entries = self.solid_compat.process_entries(entries, clean=True, inplace=inplace, n_workers=n_workers)

        # when processing single entries, all H2 polymorphs will get assigned the
        # same energy
        if len(entries) == 1 and entries[0].reduced_formula == "H2":
            warnings.warn(
                "Processing single H2 entries will result in the all polymorphs "
                "being assigned the same energy. This should not cause problems "
                "with Pourbaix diagram construction, but may be confusing. "
                "Pass all entries to process_entries() at once in if you want to "
                "preserve H2 polymorph energy differences."
            )

        # extract the DFT energies of oxygen and water from the list of entries, if present
        # do not do this when processing a single entry, as it might lead to unintended
        # results
        if len(entries) > 1:
            if not self.o2_energy and (o2_entries := [e for e in entries if e.reduced_formula == "O2"]):
                self.o2_energy = min(e.energy_per_atom for e in o2_entries)

            if not self.h2o_energy and not self.h2o_adjustments:  # noqa: SIM102
                if h2o_entries := [e for e in entries if e.reduced_formula == "H2O"]:
                    h2o_entries = sorted(h2o_entries, key=lambda e: e.energy_per_atom)
                    self.h2o_energy = h2o_entries[0].energy_per_atom
                    self.h2o_adjustments = h2o_entries[0].correction / h2o_entries[0].composition.num_atoms

        if h2_entries := [e for e in entries if e.reduced_formula == "H2"]:
            h2_entries = sorted(h2_entries, key=lambda e: e.energy_per_atom)
            self.h2_energy = h2_entries[0].energy_per_atom  # type: ignore[assignment]

        return super().process_entries(
            entries,
            clean=clean,
            verbose=verbose,
            inplace=inplace,
            n_workers=n_workers,
            on_error=on_error,
        )


def needs_u_correction(
    comp: CompositionLike,
    u_config: dict[str, dict[str, float]] = MP2020_COMPAT_CONFIG["Corrections"]["GGAUMixingCorrections"],
) -> set[str]:
    """Check if a composition is Hubbard U-corrected in the Materials Project 2020
    GGA/GGA+U mixing scheme.

    Args:
        comp (CompositionLike): The formula/composition to check.
        u_config (dict): The U-correction configuration to use. Default is the
            Materials Project 2020 configuration.

    Returns:
        set[str]: The subset of elements whose combination requires a U-correction. Pass
            return value to bool(ret_val) if you just want True/False.
    """
    elements = set(map(str, Composition(comp).elements))
    has_u_anion = set(u_config) & elements

    u_corrected_cations = set(u_config["O"])
    has_u_cation = u_corrected_cations & elements
    if has_u_cation and has_u_anion:
        return has_u_cation | has_u_anion
    return set()
