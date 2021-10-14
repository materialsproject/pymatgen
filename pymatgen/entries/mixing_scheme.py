# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
This module implements Compatibility corrections for mixing runs of different
functionals.
"""
# flake8: ignore=E712
import os
import warnings
from itertools import groupby
from typing import Optional, Union, List
import pandas as pd

import numpy as np

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.entries.entry_tools import EntrySet
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.entries.computed_entries import (
    ComputedStructureEntry,
    ConstantEnergyAdjustment,
)
from pymatgen.entries.compatibility import Compatibility, CompatibilityError, MaterialsProject2020Compatibility

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

__author__ = "Ryan Kingsbury"
__copyright__ = "Copyright 2019-2021, The Materials Project"
__version__ = "0.1"
__email__ = "RKingsbury@lbl.gov"
__date__ = "October 2021"


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
                diatomic elements O2, N2, F2, H2, and Cl2 as well as I and Br. Outputs of DFT
                relaxations using
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
            return []

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
                adjustments = self.get_adjustments(entry, mixing_scheme_state_data)
            except CompatibilityError as exc:
                if verbose:
                    print(exc)
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

    def _populate_df_row(self, struct_group, comp, sg, n, pd_type_1, pd_type_2, all_entries):
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
        stable_1 = False
        ground_state_1, ground_state_2 = np.nan, np.nan

        # get the respective hull energies at this composition, if available
        hull_energy_1, hull_energy_2 = np.nan, np.nan
        if pd_type_1:
            hull_energy_1 = pd_type_1.get_hull_energy_per_atom(comp)
        if pd_type_2:
            hull_energy_2 = pd_type_2.get_hull_energy_per_atom(comp)

        if first_entry.parameters["run_type"] in self.valid_rtypes_1:
            idx_1 = struct_group[0].index
            idx_2 = struct_group[1].index if len(struct_group) > 1 else None
            id1 = first_entry.entry_id
            id2 = second_entry.entry_id if second_entry else None
            rt1 = first_entry.parameters["run_type"]
            rt2 = second_entry.parameters["run_type"] if second_entry else None
            # are the entries the lowest energy at this composition?
            ground_state_1 = all_entries[idx_1].energy_per_atom
            ground_state_2 = all_entries[idx_2].energy_per_atom if len(struct_group) > 1 else np.nan
            # are they stable?
            if pd_type_1:
                stable_1 = first_entry in pd_type_1.stable_entries

        elif first_entry.parameters["run_type"] in self.valid_rtypes_2:
            idx_1 = struct_group[1].index if len(struct_group) > 1 else None
            idx_2 = struct_group[0].index
            id2 = first_entry.entry_id
            id1 = second_entry.entry_id if second_entry else None
            rt2 = first_entry.parameters["run_type"]
            rt1 = second_entry.parameters["run_type"] if second_entry else None
            # are the entries the lowest energy at this composition?
            ground_state_1 = all_entries[idx_1].energy_per_atom if len(struct_group) > 1 else np.nan
            ground_state_2 = all_entries[idx_2].energy_per_atom
            # are they stable?
            if pd_type_1:
                stable_1 = second_entry in pd_type_1.stable_entries
        else:
            raise ValueError("Problem - entries in group have strange run_type")

        return [
            comp.reduced_formula,
            sg,
            n,
            id1,
            id2,
            rt1,
            rt2,
            ground_state_1,
            ground_state_2,
            stable_1,
            hull_energy_1,
            hull_energy_2,
        ]

    def _generate_mixing_scheme_state_data(self, entries, verbose=False):
        """
        Generate internal state data to be passed to get_adjustments.

        Args:
            entries: The list of ComputedStructureEntry to process. This list will automatically
                be filtered based on run_type_1 and run_type_2, processed accroding to compat_1
                and compat_2, and structure matched using structure_matcher.

        Returns:
            DataFrame: A pandas DataFrame that contains information associating structures from
                different functionals with specific materials and establishing how many run_type_1
                ground states have been computed with run_type_2. The DataFrame contains one row
                for each distinct material (Structure), with the following columns:
                    composition: str the reduced_formula
                    spacegroup: int the spacegroup
                    num_sites: int the number of sites in the Structure
                    entry_id_1: the entry_id of the run_type_1 entry
                    entry_id_2: the entry_id of the run_type_2 entry
                    run_type_1: Optional[str] the run_type_1 value
                    run_type_2: Optional[str] the run_type_2 value
                    ground_state_energy_1: float or nan the ground state energy in run_type_1 in eV/atom
                    ground_state_energy_2: float or nan the ground state energy in run_type_2 in eV/atom
                    is_stable_1: bool whether this material is stable on the run_type_1 PhaseDiagram
                    hull_energy_1: float or nan the energy of the run_type_1 hull at this composition in eV/atom
                    hull_energy_2: float or nan the energy of the run_type_1 hull at this composition in eV/atom
            None: Returns None if the supplied ComputedStructureEntry are insufficient for applying
                the mixing scheme.
        """
        # if verbose:
        #     warnings.simplefilter("always")

        # separate by run_type
        entries_type_1 = [e for e in entries if e.parameters["run_type"] in self.valid_rtypes_1]
        entries_type_2 = [e for e in entries if e.parameters["run_type"] in self.valid_rtypes_2]

        # Discard entries if there are fewer than 2 for a specific run_type
        # it isn't possible to implement the mixing scheme in that case
        # this code is placed here for performance reasons, to bypass all the Structure Matching

        # We can't do any mixing without at least two of one entry type
        if len(entries_type_1) < 2 and len(entries_type_2) < 2:
            warnings.warn("Not enough entries to apply the mixing scheme. No energy adjustments will be applied.")

        if len(entries_type_1) >= 2:
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
        else:
            warnings.warn(
                f"Discarding {self.run_type_1} entries because there are fewer "
                f"than 2 of them. {self.run_type_2} energies will not be modified."
            )
            entries_type_1 = EntrySet([])

        if len(entries_type_2) >= 2:
            if self.compat_2:
                entries_type_2 = self.compat_2.process_entries(entries_type_2)
                if verbose:
                    print(
                        f"Processed {len(entries_type_2)} compatible {self.run_type_2} entries with "
                        f"{self.compat_2.__class__.__name__}"
                    )
            entries_type_2 = EntrySet(entries_type_2)
        else:
            warnings.warn(
                f"Discarding {self.run_type_2} entries because there are fewer "
                f"than 2 of them. {self.run_type_1} energies will not be modified."
            )
            entries_type_2 = EntrySet([])

        if verbose:
            print(
                f"Processing {len(entries_type_1)} {self.run_type_1} and {len(entries_type_2)} "
                f"{self.run_type_2} entries"
            )

        # Discard any entries that are not ground states from run_type_1
        # This ensures that we have at most one entry from run_type_1 at each composition
        # We must retain all run_type_2 polymorphs in case the one that matches the
        # run_type_1 ground state is not the ground state in run_type_2
        # entries_type_1.remove_non_ground_states()

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
        pd_type_1, pd_type_2 = None, None
        if len(entries_type_1) > 0:
            try:
                pd_type_1 = PhaseDiagram(entries_type_1)
            except ValueError:
                warnings.warn(
                    f"{self.run_type_1} entries do not form a complete PhaseDiagram."
                    " No energy adjustments will be applied."
                )

        if len(entries_type_2) > 0:
            try:
                pd_type_2 = PhaseDiagram(entries_type_2)
            except ValueError:
                warnings.warn(f"{self.run_type_2} entries do not form a complete PhaseDiagram.")

        # Objective: loop through all the entries, group them by structure matching (or fuzzy structure matching
        # where relevant). For each group, put a row in a pandas DataFrame with the composition of the run_type_1 entry,
        # the run_type_2 entry, whether or not that entry is a ground state (not necessarily on the hull), and its
        # e_above_hull
        all_entries = list(entries_type_1) + list(entries_type_2)
        row_list = []
        columns = [
            "composition",
            "spacegroup",
            "num_sites",
            "entry_id_1",
            "entry_id_2",
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
            l_compgroup = list(compgroup)
            # group by spacegroup, then by number of sites (for diatmics) or by structure matching
            for sg, pregroup in groupby(sorted(l_compgroup, key=_get_sg), key=_get_sg):
                l_pregroup = list(pregroup)
                if comp.reduced_formula in ["O2", "H2", "Cl2", "F2", "N2", "I", "Br"] and self.fuzzy_diatomic_matching:
                    # group by number of sites
                    for n, sitegroup in groupby(
                        sorted(l_pregroup, key=lambda s: s.num_sites), key=lambda s: s.num_sites
                    ):
                        l_sitegroup = list(sitegroup)
                        row_list.append(
                            self._populate_df_row(l_sitegroup, comp, sg, n, pd_type_1, pd_type_2, all_entries)
                        )
                else:
                    for group in self.structure_matcher.group_structures(l_pregroup):
                        grp = list(group)
                        n = group[0].num_sites
                        # StructureMatcher.group_structures returns a list of lists,
                        # so each group should be a list containing matched structures
                        row_list.append(self._populate_df_row(grp, comp, sg, n, pd_type_1, pd_type_2, all_entries))

        mixing_state_data = pd.DataFrame(row_list, columns=columns)

        if verbose:
            # how many stable entries from run_type_1 do we have in run_type_2?
            hull_entries_2 = 0
            stable_df = mixing_state_data[mixing_state_data["is_stable_1"]]
            if len(stable_df) > 0:
                hull_entries_2 = sum(stable_df["ground_state_energy_2"].notna())
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
        if len(df[df["is_stable_1"]]) > 0:
            # Any run_type_1 entry that is stable should have a ground state energy
            assert all(df[df["is_stable_1"]]["ground_state_energy_1"].notna()), "Problem with the DataFrame!"
            # If all run_type_1 stable entries are also present in run_type_2, the hull_energy_2
            # column should be populated
            if all(df[df["is_stable_1"]]["ground_state_energy_2"].notna()):
                # assert all(df["hull_energy_2"].notna())
                pass

        return df

    def get_adjustments(self, entry, mixing_scheme_state_data: pd.DataFrame = None):
        """
        Returns the corrections applied to a particular entry. Note that get_adjustments is not
        intended to be called directly in the SCAN mixing scheme. Call process_entries instead,
        and it will pass the required arguments to get_adjustments.

        Args:
            entry: A ComputedEntry object. The entry must be a member of the list of entries
                used to instantiate the class.
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
                "`mixing_scheme_state_data` DataFrame is None. No energy adjustments will be applied."
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
                f"{self.valid_rtypes_1 + self.valid_rtypes_2}. This entry will be ignored."
            )

        # The correction value depends on how many of the run_type_1 stable entries are present as run_type_2
        # calculations

        # if we do not have a complete phase diagram for either run type, there's no way to apply corrections
        type_1_hull = all(mixing_scheme_state_data["hull_energy_1"].notna())
        type_2_hull = all(mixing_scheme_state_data["hull_energy_2"].notna())

        if not type_1_hull and not type_2_hull:
            raise CompatibilityError(
                "Insufficient combination of entries for mixing scheme. No energy adjustments will be applied."
            )

        # TODO - run_1_stable might be empty!
        # run_1_stable = mixing_scheme_state_data[mixing_scheme_state_data["is_stable_1"]]
        run_1_ground_state = mixing_scheme_state_data[mixing_scheme_state_data["ground_state_energy_1"].notna()]

        # First case, ALL stable entries are present in both run types
        # in that case, prefer run_type_2 energies
        # we determine this by testing wheither EITHER all run_type_1 stable entries have run_type_2 ground state
        # energies OR all run_type_2 hull energies are present. In the latter case, we assume that the run_type_2
        # phase diagram encompasses the run_type_1 stable entries.

        # this means all run_type_1 stable states are present in run_type_2
        # therefore, we don't need to correct run_type_2
        # instead, we discard run_type_1 entries that already exist in run_type_2
        # and correct any other run_type_1 energies to the run_type_2 scale
        if all(mixing_scheme_state_data["ground_state_energy_2"].notna()) or all(
            mixing_scheme_state_data["hull_energy_2"].notna()
        ):
            if run_type in self.valid_rtypes_2:  # pylint: disable=R1705
                # For run_type_2 entries, there is no correction
                return adjustments

            # Discard GGA ground states whose structures already exist in SCAN.
            else:
                df_slice = mixing_scheme_state_data[(mixing_scheme_state_data["entry_id_1"] == entry.entry_id)]
                if len(df_slice) == 0:
                    raise CompatibilityError(
                        f"Discarding {run_type} entry {entry.entry_id} for {entry.composition.formula} "
                        f"because it was not included in the mixing state."
                    )

                if any(df_slice["ground_state_energy_1"] == entry.energy_per_atom) and any(
                    df_slice["ground_state_energy_2"].notna()
                ):
                    # discard entry if its energy matches one in the DataFrame,
                    # because it already exists in run_type_2
                    raise CompatibilityError(
                        f"Discarding {run_type} entry {entry.entry_id} for {entry.composition.formula} "
                        f"whose structure already exists in {self.run_type_2}"
                    )

                # If a GGA is not present in SCAN, correct its energy to give the same
                # e_above_hull on the SCAN hull that it would have on the GGA hull
                hull_energy_1 = df_slice["hull_energy_1"].iloc[0]
                hull_energy_2 = df_slice["hull_energy_2"].iloc[0]

                if not np.isfinite(hull_energy_2):
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

        # Third case, there are run_type_2 energies available for at least some run_type_1
        # reference states. Here, we can correct run_type_2 energies onto the run_type_1 scale
        # by replacing run_type_1 e_above_hull with a run_type_2 e_above_hull for unstable polymorphs
        else:
            if run_type in self.valid_rtypes_1:  # pylint: disable=R1705
                # For run_type_1 entries, there is no correction
                return adjustments

            elif run_type in self.valid_rtypes_2 and any(
                mixing_scheme_state_data[mixing_scheme_state_data["composition"] == entry.composition.reduced_formula][
                    "ground_state_energy_1"
                ].notna()
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
                    mixing_scheme_state_data["composition"] == entry.composition.reduced_formula
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
