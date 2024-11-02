"""This module implements Compatibility corrections for mixing runs of different
functionals.
"""

from __future__ import annotations

import copy
import os
import warnings
from itertools import groupby

import numpy as np
import pandas as pd

from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.entries.compatibility import (
    AnyComputedEntry,
    Compatibility,
    CompatibilityError,
    MaterialsProject2020Compatibility,
)
from pymatgen.entries.computed_entries import ComputedStructureEntry, ConstantEnergyAdjustment
from pymatgen.entries.entry_tools import EntrySet

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

__author__ = "Ryan Kingsbury"
__copyright__ = "Copyright 2019-2021, The Materials Project"
__version__ = "0.1"
__email__ = "RKingsbury@lbl.gov"
__date__ = "October 2021"


class MaterialsProjectDFTMixingScheme(Compatibility):
    """This class implements the Materials Project mixing scheme, which allows mixing of
    energies from different DFT functionals. Note that this should only be used for
    VASP calculations using the MaterialsProject parameters (e.g. MPRelaxSet or
    MPScanRelaxSet). Using this compatibility scheme on runs with different parameters
    may lead to unexpected results.

    This is the scheme used by the Materials Project to generate Phase Diagrams containing
    a mixture of GGA(+U) and R2SCAN calculations. However in principle it can be used to
    mix energies from any two functionals.
    """

    def __init__(
        self,
        structure_matcher: StructureMatcher | None = None,
        run_type_1: str = "GGA(+U)",
        run_type_2: str = "R2SCAN",
        compat_1: (Compatibility | None) = MaterialsProject2020Compatibility(),  # noqa: B008
        compat_2: Compatibility | None = None,
        fuzzy_matching: bool = True,
        check_potcar: bool = True,
    ) -> None:
        """Instantiate the mixing scheme. The init method creates a generator class that
        contains relevant settings (e.g., StructureMatcher instance, Compatibility settings
        for each functional) for processing groups of entries.

        Args:
            structure_matcher (StructureMatcher): StructureMatcher object used to determine
                whether calculations from different functionals describe the same material.
            run_type_1: The first DFT run_type. Typically this is the majority or run type or
                the "base case" onto which the other calculations are referenced. Valid choices
                are any run_type recognized by Vasprun.run_type, such as "LDA", "GGA", "GGA+U",
                "PBEsol", "SCAN", or "R2SCAN". The class will ignore any entries that have a
                run_type different than run_type_1 or run_type_2.

                The list of run_type_1 entries provided to process_entries MUST form a complete
                Phase Diagram in order for the mixing scheme to work. If this condition is not
                satisfied, processing the entries will fail.

                Note that the special string "GGA(+U)" (default) will treat both GGA and GGA+U
                calculations as a single type. This option exists because GGA/GGA+U mixing is
                already handled by MaterialsProject2020Compatibility.
            run_type_2: The second DFT run_type. Typically this is the run_type that is 'preferred'
                but has fewer calculations. If run_type_1 and run_type_2 calculations exist for all
                materials, run_type_2 energies will be used (hence the 'preferred' status). The class
                will ignore any entries that have a run_type different than run_type_1 or run_type_2.
            compat_1: Compatibility class used to pre-process entries of run_type_1.
                Defaults to MaterialsProjectCompatibility2020.
            compat_2: Compatibility class used to pre-process entries of run_type_2.
                Defaults to None.
            fuzzy_matching: Whether to use less strict structure matching logic for
                diatomic elements O2, N2, F2, H2, and Cl2 as well as I and Br. Outputs of DFT
                relaxations using
                different functionals frequently fail to structure match for these elements
                even though they come from the same original material. Fuzzy structure matching
                considers the materials equivalent if the formula, number of sites, and
                space group are all identical. If there are multiple materials of run_type_2
                that satisfy these criteria, the one with lowest energy is considered to
                match.
            check_potcar: Whether to ensure the POTCARs used for the run_type_1 and run_type_2 calculations
                are the same. This is useful for ensuring that the mixing scheme is not used on calculations
                that used different POTCARs, which can lead to unphysical results. Defaults to True.
                Has no effect if neither compat_1 nor compat_2 have a check_potcar attribute.
                Can also be disabled globally by running `pmg config --add PMG_POTCAR_CHECKS false`.
        """
        self.name = "MP DFT mixing scheme"
        self.structure_matcher = structure_matcher or StructureMatcher()
        if run_type_1 == run_type_2:
            raise ValueError(f"run_type_1={run_type_2=}. The mixing scheme is meaningless unless run_types different")
        self.run_type_1 = run_type_1
        self.run_type_2 = run_type_2
        self.valid_rtypes_1 = ["GGA", "GGA+U"] if self.run_type_1 == "GGA(+U)" else [self.run_type_1]
        self.valid_rtypes_2 = ["GGA", "GGA+U"] if self.run_type_2 == "GGA(+U)" else [self.run_type_2]

        self.compat_1 = compat_1
        self.compat_2 = compat_2
        self.fuzzy_matching = fuzzy_matching
        self.check_potcar = check_potcar
        for compat in (self.compat_1, self.compat_2):
            if hasattr(compat, "check_potcar"):
                compat.check_potcar = check_potcar  # type: ignore[union-attr]

    def process_entries(
        self,
        entries: AnyComputedEntry | list[AnyComputedEntry],
        clean: bool = True,
        verbose: bool = False,
        inplace: bool = True,
        mixing_state_data=None,
    ) -> list[AnyComputedEntry]:
        """Process a sequence of entries with the DFT mixing scheme. Note
        that this method will change the data of the original entries.

        Args:
            entries: ComputedEntry or [ComputedEntry]. Pass all entries as a single list, even if they are
                computed with different functionals or require different preprocessing. This list will
                automatically be filtered based on run_type_1 and run_type_2, and processed according to
                compat_1 and compat_2.

                Note that under typical use, when mixing_state_data=None, the entries MUST be
                ComputedStructureEntry. They will be matched using structure_matcher.
            clean (bool): Whether to remove any previously-applied energy adjustments.
                If True, all EnergyAdjustment are removed prior to processing the Entry.
                Default is True.
            verbose (bool): Whether to print verbose error messages about the mixing scheme. Default is False.
            inplace (bool): Whether to adjust input entries in place. Default is True.
            mixing_state_data: A DataFrame containing information about which Entries
                correspond to the same materials, which are stable on the phase diagrams of
                the respective run_types, etc. If None (default), it will be generated from the
                list of entries using MaterialsProjectDFTMixingScheme.get_mixing_state_data.
                This argument is included to facilitate use of the mixing scheme in high-throughput
                databases where an alternative to get_mixing_state_data is desirable for performance
                reasons. In general, it should always be left at the default value (None) to avoid
                inconsistencies between the mixing state data and the properties of the
                ComputedStructureEntry in entries.

        Returns:
            list[AnyComputedEntry]: Adjusted entries. Entries in the original list incompatible with
                chosen correction scheme are excluded from the returned list.
        """
        processed_entry_list: list = []

        # We can't operate on single entries in this scheme
        if len(entries) == 1:
            warnings.warn(f"{type(self).__name__} cannot process single entries. Supply a list of entries.")
            return processed_entry_list

        # if inplace = False, process entries on a copy
        if not inplace:
            entries = copy.deepcopy(entries)

        # if clean is True, remove all previous adjustments from the entry
        # this code must be placed before the next block, because we don't want to remove
        # any corrections added by compat_1 or compat_2.
        if clean:
            for entry in entries:
                for ea in entry.energy_adjustments:
                    entry.energy_adjustments.remove(ea)

        entries_type_1, entries_type_2 = self._filter_and_sort_entries(entries, verbose=verbose)

        if mixing_state_data is None:
            if verbose:
                print("  Generating mixing state data from provided entries.")
            mixing_state_data = self.get_mixing_state_data(entries_type_1 + entries_type_2)

        if verbose:
            # how many stable entries from run_type_1 do we have in run_type_2?
            hull_entries_2 = 0
            stable_df = mixing_state_data[mixing_state_data["is_stable_1"]]
            if len(stable_df) > 0:
                hull_entries_2 = sum(stable_df["energy_2"].notna())
            print(
                f"  Entries contain {self.run_type_2} calculations for {hull_entries_2} of {len(stable_df)} "
                f"{self.run_type_1} hull entries."
            )
            if hull_entries_2 == len(stable_df):
                print(f"  {self.run_type_1} energies will be adjusted to the {self.run_type_2} scale")
            else:
                print(f"  {self.run_type_2} energies will be adjusted to the {self.run_type_1} scale")

            if hull_entries_2 > 0:
                print(
                    f"  The energy above hull for {self.run_type_2} materials at compositions with "
                    f"{self.run_type_2} hull entries will be preserved. For other compositions, "
                    f"Energies of {self.run_type_2} materials will be set equal to those of "
                    f"matching {self.run_type_1} materials"
                )

        # the code below is identical to code inside process_entries in the base
        # Compatibility class, except that an extra kwarg is passed to get_adjustments
        for entry in entries_type_1 + entries_type_2:
            ignore_entry = False
            # get the energy adjustments
            try:
                adjustments = self.get_adjustments(entry, mixing_state_data)
            except CompatibilityError as exc:
                if "WARNING!" in str(exc):
                    warnings.warn(str(exc))
                elif verbose:
                    print(f"  {exc}")
                ignore_entry = True
                continue

            for ea in adjustments:
                # Has this correction already been applied?
                if (ea.name, ea.cls, ea.value) in [(ea2.name, ea2.cls, ea2.value) for ea2 in entry.energy_adjustments]:
                    # we already applied this exact correction. Do nothing.
                    pass
                elif (ea.name, ea.cls) in [(ea2.name, ea2.cls) for ea2 in entry.energy_adjustments]:
                    # we already applied a correction with the same name
                    # but a different value. Something is wrong.
                    ignore_entry = True
                    warnings.warn(
                        f"Entry {entry.entry_id} already has an energy adjustment called {ea.name}, but its "
                        f"value differs from the value of {ea.value:.3f} calculated here. This "
                        "Entry will be discarded."
                    )
                else:
                    # Add the correction to the energy_adjustments list
                    entry.energy_adjustments.append(ea)

            if not ignore_entry:
                processed_entry_list.append(entry)

        if verbose:
            count_type_1 = sum(entry.parameters["run_type"] in self.valid_rtypes_1 for entry in processed_entry_list)
            count_type_2 = sum(entry.parameters["run_type"] in self.valid_rtypes_2 for entry in processed_entry_list)
            print(
                f"\nProcessing complete. Mixed entries contain {count_type_1} {self.run_type_1} and {count_type_2} "
                f"{self.run_type_2} entries.\n"
            )
            self.display_entries(processed_entry_list)

        return processed_entry_list

    def get_adjustments(self, entry, mixing_state_data: pd.DataFrame | None = None):
        """Get the corrections applied to a particular entry. Note that get_adjustments is not
        intended to be called directly in the R2SCAN mixing scheme. Call process_entries instead,
        and it will pass the required arguments to get_adjustments.

        Args:
            entry: A ComputedEntry object. The entry must be a member of the list of entries
                used to create mixing_state_data.
            mixing_state_data: A DataFrame containing information about which Entries
                correspond to the same materials, which are stable on the phase diagrams of
                the respective run_types, etc. Can be generated from a list of entries using
                MaterialsProjectDFTMixingScheme.get_mixing_state_data. This argument is included to
                facilitate use of the mixing scheme in high-throughput databases where an alternative
                to get_mixing_state_data is desirable for performance reasons. In general, it should
                always be left at the default value (None) to avoid inconsistencies between the mixing
                state data and the properties of the ComputedStructureEntry.

        Returns:
            [EnergyAdjustment]: Energy adjustments to be applied to entry.

        Raises:
            CompatibilityError if the DFT mixing scheme cannot be applied to the entry.
        """
        adjustments: list[ConstantEnergyAdjustment] = []
        run_type = entry.parameters.get("run_type")

        if mixing_state_data is None:
            raise CompatibilityError(
                "WARNING! `mixing_state_data` DataFrame is None. No energy adjustments will be applied."
            )

        if not all(mixing_state_data["hull_energy_1"].notna()) and any(mixing_state_data["entry_id_1"].notna()):
            raise CompatibilityError(
                f"WARNING! {self.run_type_1} entries do not form a complete PhaseDiagram."
                " No energy adjustments will be applied."
            )

        if run_type not in self.valid_rtypes_1 + self.valid_rtypes_2:
            raise CompatibilityError(
                f"WARNING! Invalid {run_type=} for entry {entry.entry_id}. Must be one of "
                f"{self.valid_rtypes_1 + self.valid_rtypes_2}. This entry will be ignored."
            )

        # Verify that the entry is included in the mixing state data
        if (entry.entry_id not in mixing_state_data["entry_id_1"].to_numpy()) and (
            entry.entry_id not in mixing_state_data["entry_id_2"].to_numpy()
        ):
            raise CompatibilityError(
                f"WARNING! Discarding {run_type} entry {entry.entry_id} for {entry.formula} "
                f"because it was not found in the mixing state data. This can occur when there are duplicate "
                "structures. In such cases, only the lowest energy entry with that structure appears in the "
                "mixing state data."
            )

        # Verify that the entry's energy has not been modified since mixing state data was generated
        if (entry.energy_per_atom not in mixing_state_data["energy_1"].to_numpy()) and (
            entry.energy_per_atom not in mixing_state_data["energy_2"].to_numpy()
        ):
            raise CompatibilityError(
                f"WARNING! Discarding {run_type} entry {entry.entry_id} for {entry.formula} "
                "because it's energy has been modified since the mixing state data was generated."
            )

        # Compute the energy correction for mixing. The correction value depends on how many of the
        # run_type_1 stable entries are present as run_type_2 calculations

        # First case - ALL run_type_1 stable entries are present in run_type_2
        # In this scenario we construct the hull using run_type_2 energies. We discard any
        # run_type_1 entries that already exist in run_type_2 and correct other run_type_1
        # energies to have the same e_above_hull on the run_type_2 hull as they had on the run_type_1 hull
        if all(mixing_state_data[mixing_state_data["is_stable_1"]]["entry_id_2"].notna()):
            if run_type in self.valid_rtypes_2:
                # For run_type_2 entries, there is no correction
                return adjustments

            # Discard GGA ground states whose structures already exist in R2SCAN.
            df_slice = mixing_state_data[(mixing_state_data["entry_id_1"] == entry.entry_id)]

            if df_slice["entry_id_2"].notna().item():
                # there is a matching run_type_2 entry, so we will discard this entry
                if df_slice["is_stable_1"].item():
                    # this is a GGA ground state.
                    raise CompatibilityError(
                        f"Discarding {run_type} entry {entry.entry_id} for {entry.formula} "
                        f"because it is a {self.run_type_1} ground state that matches a {self.run_type_2} "
                        "material."
                    )

                raise CompatibilityError(
                    f"Discarding {run_type} entry {entry.entry_id} for {entry.formula} "
                    f"because there is a matching {self.run_type_2} material."
                )

            # If a GGA is not present in R2SCAN, correct its energy to give the same
            # e_above_hull on the R2SCAN hull that it would have on the GGA hull
            hull_energy_1 = df_slice["hull_energy_1"].iloc[0]
            hull_energy_2 = df_slice["hull_energy_2"].iloc[0]
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

        # Second case - there are run_type_2 energies available for at least some run_type_1
        # stable entries. Here, we can correct run_type_2 energies at certain compositions
        # to preserve their e_above_hull on the run_type_1 hull
        if any(mixing_state_data[mixing_state_data["is_stable_1"]]["entry_id_2"].notna()):
            if run_type in self.valid_rtypes_1:
                df_slice = mixing_state_data[mixing_state_data["entry_id_1"] == entry.entry_id]

                if df_slice["entry_id_2"].notna().item():
                    # there is a matching run_type_2 entry. We should discard this entry
                    if df_slice["is_stable_1"].item():
                        # this is a GGA ground state.
                        raise CompatibilityError(
                            f"Discarding {run_type} entry {entry.entry_id} for {entry.formula} "
                            f"because it is a {self.run_type_1} ground state that matches a {self.run_type_2} "
                            "material."
                        )

                    raise CompatibilityError(
                        f"Discarding {run_type} entry {entry.entry_id} for {entry.formula} "
                        f"because there is a matching {self.run_type_2} material"
                    )

                # For other run_type_1 entries, there is no correction
                return adjustments

            # for run_type_2, determine whether there is a run_type_2 ground state at this composition
            df_slice = mixing_state_data[mixing_state_data["formula"] == entry.reduced_formula]

            if any(df_slice[df_slice["is_stable_1"]]["entry_id_2"].notna()):
                # there is a run_type_2 entry corresponding to the run_type_1 ground state
                # adjust the run_type_2 energy to preserve the e_above_hull
                gs_energy_type_2 = df_slice[df_slice["is_stable_1"]]["energy_2"].item()
                e_above_hull = entry.energy_per_atom - gs_energy_type_2
                hull_energy_1 = df_slice["hull_energy_1"].iloc[0]
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

            # this composition is not stable in run_type_1. If the run_type_2 entry matches a run_type_1
            # entry, we can adjust the run_type_2 energy to match the run_type_1 energy.
            if any(df_slice[df_slice["entry_id_2"] == entry.entry_id]["entry_id_1"].notna()):
                # adjust the energy of the run_type_2 entry to match that of the run_type_1 entry
                type_1_energy = df_slice[df_slice["entry_id_2"] == entry.entry_id]["energy_1"].iloc[0]
                correction = (type_1_energy - entry.energy_per_atom) * entry.composition.num_atoms
                adjustments.append(
                    ConstantEnergyAdjustment(
                        correction,
                        0.0,
                        name=f"MP {self.run_type_1}/{self.run_type_2} mixing adjustment",
                        cls=self.as_dict(),
                        description=f"Replace {self.run_type_2} energy with {self.run_type_1} energy",
                    )
                )
                return adjustments

            # there is no run_type_1 entry that matches this material, and no ground state. Discard.
            raise CompatibilityError(
                f"Discarding {run_type} entry {entry.entry_id} for {entry.formula} "
                f"because there is no matching {self.run_type_1} entry and no {self.run_type_2} "
                "ground state at this composition."
            )

        # Third case - there are no run_type_2 energies available for any run_type_1
        # ground states. There's no way to use the run_type_2 energies in this case.
        if all(mixing_state_data[mixing_state_data["is_stable_1"]]["entry_id_2"].isna()):
            if run_type in self.valid_rtypes_1:
                # nothing to do for run_type_1, return as is
                return adjustments

            # for run_type_2, discard the entry
            raise CompatibilityError(
                f"Discarding {run_type} entry {entry.entry_id} for {entry.formula} "
                f"because there are no {self.run_type_2} ground states at this composition."
            )

        # this statement is here to guarantee a return or raise
        raise CompatibilityError(
            "WARNING! If you see this Exception it means you have encountered"
            f"an edge case in {type(self).__name__}. Inspect your input carefully and post a bug report."
        )

    def get_mixing_state_data(self, entries: list[ComputedStructureEntry]):
        """Generate internal state data to be passed to get_adjustments.

        Args:
            entries: The list of ComputedStructureEntry to process. It is assumed that the entries have
                already been filtered using _filter_and_sort_entries() to remove any irrelevant run types,
                apply compat_1 and compat_2, and confirm that all have unique entry_id.

        Returns:
            DataFrame: A pandas DataFrame that contains information associating structures from
                different functionals with specific materials and establishing how many run_type_1
                ground states have been computed with run_type_2. The DataFrame contains one row
                for each distinct material (Structure), with the following columns:
                    formula: str the reduced_formula
                    spacegroup: int the spacegroup
                    num_sites: int the number of sites in the Structure
                    entry_id_1: the entry_id of the run_type_1 entry
                    entry_id_2: the entry_id of the run_type_2 entry
                    run_type_1: Optional[str] the run_type_1 value
                    run_type_2: Optional[str] the run_type_2 value
                    energy_1: float or nan the ground state energy in run_type_1 in eV/atom
                    energy_2: float or nan the ground state energy in run_type_2 in eV/atom
                    is_stable_1: bool whether this material is stable on the run_type_1 PhaseDiagram
                    hull_energy_1: float or nan the energy of the run_type_1 hull at this composition in eV/atom
                    hull_energy_2: float or nan the energy of the run_type_1 hull at this composition in eV/atom
            None: Returns None if the supplied ComputedStructureEntry are insufficient for applying
                the mixing scheme.
        """
        filtered_entries = []

        for entry in entries:
            if not isinstance(entry, ComputedStructureEntry):
                warnings.warn(
                    f"Entry {entry.entry_id} is not a ComputedStructureEntry and will be ignored. "
                    "The DFT mixing scheme requires structures for all entries"
                )
                continue

            filtered_entries.append(entry)

        # separate by run_type
        entries_type_1 = [e for e in filtered_entries if e.parameters["run_type"] in self.valid_rtypes_1]
        entries_type_2 = [e for e in filtered_entries if e.parameters["run_type"] in self.valid_rtypes_2]

        # construct PhaseDiagram for each run_type, if possible
        pd_type_1, pd_type_2 = None, None
        try:
            pd_type_1 = PhaseDiagram(entries_type_1)
        except ValueError:
            warnings.warn(f"{self.run_type_1} entries do not form a complete PhaseDiagram.")

        try:
            pd_type_2 = PhaseDiagram(entries_type_2)
        except ValueError:
            warnings.warn(f"{self.run_type_2} entries do not form a complete PhaseDiagram.")

        # Objective: loop through all the entries, group them by structure matching (or fuzzy structure matching
        # where relevant). For each group, put a row in a pandas DataFrame with the composition of the run_type_1 entry,
        # the run_type_2 entry, whether or not that entry is a ground state (not necessarily on the hull), its energy,
        # and the energy of the hull at that composition
        all_entries = list(entries_type_1) + list(entries_type_2)
        row_list = []
        columns = [
            "formula",
            "spacegroup",
            "num_sites",
            "is_stable_1",
            "entry_id_1",
            "entry_id_2",
            "run_type_1",
            "run_type_2",
            "energy_1",
            "energy_2",
            "hull_energy_1",
            "hull_energy_2",
        ]

        def _get_sg(struct) -> int:
            """Helper function to get spacegroup with a loose tolerance."""
            try:
                return struct.get_space_group_info(symprec=0.1)[1]
            except Exception:
                return -1

        # loop through all structures
        # this logic follows emmet.builders.vasp.materials.MaterialsBuilder.filter_and_group_tasks
        structures = []
        for entry in all_entries:
            struct = entry.structure
            struct.entry_id = entry.entry_id
            structures.append(struct)

        # First group by composition, then by spacegroup number, then by structure matching
        for comp, comp_group in groupby(sorted(structures, key=lambda s: s.composition), key=lambda s: s.composition):
            l_comp_group = list(comp_group)
            # group by spacegroup, then by number of sites (for diatmics) or by structure matching
            for sg, pre_group in groupby(sorted(l_comp_group, key=_get_sg), key=_get_sg):
                l_pre_group = list(pre_group)
                if comp.reduced_formula in ["O2", "H2", "Cl2", "F2", "N2", "I", "Br", "H2O"] and self.fuzzy_matching:
                    # group by number of sites
                    for idx, site_group in groupby(sorted(l_pre_group, key=len), key=len):
                        l_sitegroup = list(site_group)
                        row_list.append(
                            self._populate_df_row(
                                l_sitegroup,
                                comp,
                                sg,
                                idx,
                                pd_type_1,
                                pd_type_2,
                                all_entries,
                            )
                        )
                else:
                    for group in self.structure_matcher.group_structures(l_pre_group):
                        group = list(group)
                        idx = len(group[0])
                        # StructureMatcher.group_structures returns a list of lists,
                        # so each group should be a list containing matched structures
                        row_list.append(self._populate_df_row(group, comp, sg, idx, pd_type_1, pd_type_2, all_entries))

        mixing_state_data = pd.DataFrame(row_list, columns=columns)
        return mixing_state_data.sort_values(["formula", "energy_1", "spacegroup", "num_sites"], ignore_index=True)

    def _filter_and_sort_entries(self, entries, verbose=False):
        """Given a single list of entries, separate them by run_type and return two lists, one containing
        only entries of each run_type.
        """
        filtered_entries = []

        for entry in entries:
            entry_id = entry.entry_id
            if not entry.parameters.get("run_type"):
                warnings.warn(
                    f"Entry {entry_id} is missing parameters.run_type! This field"
                    "is required. This entry will be ignored."
                )
                continue

            run_type = entry.parameters.get("run_type")
            if run_type not in [*self.valid_rtypes_1, *self.valid_rtypes_2]:
                warnings.warn(
                    f"Invalid {run_type=} for entry {entry_id}. Must be one of "
                    f"{self.valid_rtypes_1 + self.valid_rtypes_2}. This entry will be ignored."
                )
                continue

            formula = entry.reduced_formula
            if entry_id is None:
                warnings.warn(
                    f"{entry_id=} for {formula=}. Unique entry_ids are required for every ComputedStructureEntry."
                    " This entry will be ignored."
                )
                continue

            filtered_entries.append(entry)

        filtered_entry_ids = {e.entry_id for e in filtered_entries}
        if len(filtered_entry_ids) != len(filtered_entries):
            raise ValueError(
                "The provided ComputedStructureEntry do not all have unique entry_ids."
                " Unique entry_ids are required for every ComputedStructureEntry."
            )

        # separate by run_type
        entries_type_1 = [e for e in filtered_entries if e.parameters["run_type"] in self.valid_rtypes_1]
        entries_type_2 = [e for e in filtered_entries if e.parameters["run_type"] in self.valid_rtypes_2]

        if verbose:
            print(
                f"Processing {len(entries_type_1)} {self.run_type_1} and {len(entries_type_2)} "
                f"{self.run_type_2} entries..."
            )

        # preprocess entries with any corrections
        # make an EntrySet to enable some useful methods like .chemsys and .is_ground_state
        if self.compat_1:
            entries_type_1 = self.compat_1.process_entries(entries_type_1)
            if verbose:
                print(
                    f"  Processed {len(entries_type_1)} compatible {self.run_type_1} entries with "
                    f"{type(self.compat_1).__name__}"
                )
        entries_type_1 = EntrySet(entries_type_1)

        if self.compat_2:
            entries_type_2 = self.compat_2.process_entries(entries_type_2)
            if verbose:
                print(
                    f"  Processed {len(entries_type_2)} compatible {self.run_type_2} entries with "
                    f"{type(self.compat_2).__name__}"
                )
        entries_type_2 = EntrySet(entries_type_2)

        # make sure both sets of entries belong to the same chemical system
        # assuming there are any gga entries at all
        if len(entries_type_1.chemsys) > 0:
            chemsys = entries_type_1.chemsys
            if not entries_type_2.chemsys <= entries_type_1.chemsys:
                warnings.warn(
                    f"  {self.run_type_2} entries chemical system {entries_type_2.chemsys} is larger than "
                    f"{self.run_type_1} entries chemical system {entries_type_1.chemsys}. Entries outside the "
                    f"{self.run_type_1} chemical system will be discarded"
                )
                entries_type_2 = entries_type_2.get_subset_in_chemsys(chemsys)
        else:
            # if only run_type_2 entries are present, then they define the chemsys
            chemsys = entries_type_2.chemsys

        if verbose:
            print(f"  Entries belong to the {chemsys} chemical system")

        return list(entries_type_1), list(entries_type_2)

    def _populate_df_row(self, struct_group, comp, sg, n, pd_type_1, pd_type_2, all_entries):
        """Helper function to populate a row of the mixing state DataFrame, given
        a list of matched structures.
        """
        # within the group of matched structures, keep the lowest energy entry from
        # each run_type
        entries_type_1 = sorted(
            (
                e
                for e in all_entries
                if e.entry_id in [s.entry_id for s in struct_group] and e.parameters["run_type"] in self.valid_rtypes_1
            ),
            key=lambda x: x.energy_per_atom,
        )
        first_entry = entries_type_1[0] if len(entries_type_1) > 0 else None

        entries_type_2 = sorted(
            (
                e
                for e in all_entries
                if e.entry_id in [s.entry_id for s in struct_group] and e.parameters["run_type"] in self.valid_rtypes_2
            ),
            key=lambda x: x.energy_per_atom,
        )
        second_entry = entries_type_2[0] if len(entries_type_2) > 0 else None

        # generate info for the DataFrame
        stable_1 = False

        id1 = first_entry.entry_id if first_entry else None
        id2 = second_entry.entry_id if second_entry else None
        rt1 = first_entry.parameters["run_type"] if first_entry else None
        rt2 = second_entry.parameters["run_type"] if second_entry else None
        # are the entries the lowest energy at this composition?
        energy_1 = first_entry.energy_per_atom if first_entry else np.nan
        energy_2 = second_entry.energy_per_atom if second_entry else np.nan
        # are they stable?
        if pd_type_1:
            stable_1 = first_entry in pd_type_1.stable_entries

        # get the respective hull energies at this composition, if available
        hull_energy_1, hull_energy_2 = np.nan, np.nan
        if pd_type_1:
            hull_energy_1 = pd_type_1.get_hull_energy_per_atom(comp)
        if pd_type_2:
            hull_energy_2 = pd_type_2.get_hull_energy_per_atom(comp)

        return [
            comp.reduced_formula,
            sg,
            n,
            stable_1,
            id1,
            id2,
            rt1,
            rt2,
            energy_1,
            energy_2,
            hull_energy_1,
            hull_energy_2,
        ]

    @staticmethod
    def display_entries(entries):
        """Generate a pretty printout of key properties of a list of ComputedEntry."""
        entries = sorted(entries, key=lambda e: (e.reduced_formula, e.energy_per_atom))
        try:
            pd = PhaseDiagram(entries)
        except ValueError:
            return

        print(
            f"{'entry_id':<12}{'formula':<12}{'spacegroup':<12}{'run_type':<10}{'eV/atom':<8}"
            f"{'corr/atom':<9} {'e_above_hull':<9}"
        )
        for entry in entries:
            print(
                f"{entry.entry_id:<12}{entry.reduced_formula:<12}{entry.structure.get_space_group_info()[0]:<12}"
                f"{entry.parameters['run_type']:<10}{entry.energy_per_atom:<8.3f}"
                f"{entry.correction / entry.composition.num_atoms:<9.3f} {pd.get_e_above_hull(entry):<9.3f}"
            )
        return
