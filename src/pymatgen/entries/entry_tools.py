"""This module implements functions to perform various useful operations on
entries, such as grouping entries by structure.
"""

from __future__ import annotations

import collections
import csv
import itertools
import json
import logging
import multiprocessing as mp
import re
from collections import defaultdict
from datetime import datetime, timezone
from typing import TYPE_CHECKING

from monty.json import MontyDecoder, MontyEncoder, MSONable

from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.analysis.structure_matcher import SpeciesComparator, StructureMatcher
from pymatgen.core import Composition, Element

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Literal

    from typing_extensions import Self

    from pymatgen.entries import Entry
    from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
    from pymatgen.util.typing import SpeciesLike

logger = logging.getLogger(__name__)


def _get_host(structure, species_to_remove):
    if species_to_remove:
        struct = structure.copy()
        struct.remove_species(species_to_remove)
        return struct
    return structure


def _perform_grouping(args):
    (
        entries_json,
        hosts_json,
        ltol,
        stol,
        angle_tol,
        primitive_cell,
        scale,
        comparator,
        groups,
    ) = args

    entries = json.loads(entries_json, cls=MontyDecoder)
    hosts = json.loads(hosts_json, cls=MontyDecoder)
    unmatched = list(zip(entries, hosts, strict=True))
    while len(unmatched) > 0:
        ref_host = unmatched[0][1]
        logger.info(f"Reference tid = {unmatched[0][0].entry_id}, formula = {ref_host.formula}")
        ref_formula = ref_host.reduced_formula
        logger.info(f"Reference host = {ref_formula}")
        matches = [unmatched[0]]
        for idx in range(1, len(unmatched)):
            test_host = unmatched[idx][1]
            logger.info(f"Testing tid = {unmatched[idx][0].entry_id}, formula = {test_host.formula}")
            test_formula = test_host.reduced_formula
            logger.info(f"Test host = {test_formula}")
            matcher = StructureMatcher(
                ltol=ltol,
                stol=stol,
                angle_tol=angle_tol,
                primitive_cell=primitive_cell,
                scale=scale,
                comparator=comparator,
            )
            if matcher.fit(ref_host, test_host):
                logger.info("Fit found")
                matches.append(unmatched[idx])
        groups.append(json.dumps([m[0] for m in matches], cls=MontyEncoder))
        unmatched = list(filter(lambda x: x not in matches, unmatched))
        logger.info(f"{len(unmatched)} unmatched remaining")


def group_entries_by_structure(
    entries,
    species_to_remove=None,
    ltol=0.2,
    stol=0.4,
    angle_tol=5,
    primitive_cell=True,
    scale=True,
    comparator=None,
    ncpus=None,
):
    """Given a sequence of ComputedStructureEntries, use structure fitter to group
    them by structural similarity.

    Args:
        entries: Sequence of ComputedStructureEntries.
        species_to_remove: Sometimes you want to compare a host framework
            (e.g., in Li-ion battery analysis). This allows you to specify
            species to remove before structural comparison.
        ltol (float): Fractional length tolerance. Default is 0.2.
        stol (float): Site tolerance in Angstrom. Default is 0.4 Angstrom.
        angle_tol (float): Angle tolerance in degrees. Default is 5 degrees.
        primitive_cell (bool): If true: input structures will be reduced to
            primitive cells prior to matching. Defaults to True.
        scale: Input structures are scaled to equivalent volume if true;
            For exact matching, set to False.
        comparator: A comparator object implementing an equals method that
            declares equivalency of sites. Default is SpeciesComparator,
            which implies rigid species mapping.
        ncpus: Number of cpus to use. Use of multiple cpus can greatly improve
            fitting speed. Default of None means serial processing.

    Returns:
        Sequence of sequence of entries by structural similarity. e.g,
        [[ entry1, entry2], [entry3, entry4, entry5]]
    """
    if comparator is None:
        comparator = SpeciesComparator()
    start = datetime.now(tz=timezone.utc)
    logger.info(f"Started at {start}")
    entries_host = [(entry, _get_host(entry.structure, species_to_remove)) for entry in entries]
    if ncpus:
        symm_entries = defaultdict(list)
        for entry, host in entries_host:
            symm_entries[comparator.get_structure_hash(host)].append((entry, host))

        logger.info(f"Using {ncpus} cpus")
        manager = mp.Manager()
        groups = manager.list()
        with mp.Pool(ncpus) as pool:
            # Parallel processing only supports Python primitives and not objects.
            pool.map(
                _perform_grouping,
                [
                    (
                        json.dumps([e[0] for e in eh], cls=MontyEncoder),
                        json.dumps([e[1] for e in eh], cls=MontyEncoder),
                        ltol,
                        stol,
                        angle_tol,
                        primitive_cell,
                        scale,
                        comparator,
                        groups,
                    )
                    for eh in symm_entries.values()
                ],
            )
    else:
        groups = []
        hosts = [host for entry, host in entries_host]
        _perform_grouping(
            (
                json.dumps(entries, cls=MontyEncoder),
                json.dumps(hosts, cls=MontyEncoder),
                ltol,
                stol,
                angle_tol,
                primitive_cell,
                scale,
                comparator,
                groups,
            )
        )
    entry_groups = []
    for g in groups:
        entry_groups.append(json.loads(g, cls=MontyDecoder))
    logger.info(f"Finished at {datetime.now(tz=timezone.utc)}")
    logger.info(f"Took {datetime.now(tz=timezone.utc) - start}")
    return entry_groups


def group_entries_by_composition(entries, sort_by_e_per_atom=True):
    """Given a sequence of Entry-like objects, group them by composition and
        optionally sort by energy above hull.

    Args:
        entries (list): Sequence of Entry-like objects.
        sort_by_e_per_atom (bool): Whether to sort the grouped entries by
            energy per atom (lowest energy first). Default True.

    Returns:
        Sequence of sequence of entries by composition. e.g,
        [[ entry1, entry2], [entry3, entry4, entry5]]
    """
    entry_groups = []
    entries = sorted(entries, key=lambda e: e.reduced_formula)
    for _, g in itertools.groupby(entries, key=lambda e: e.reduced_formula):
        group = list(g)
        if sort_by_e_per_atom:
            group = sorted(group, key=lambda e: e.energy_per_atom)

        entry_groups.append(group)

    return entry_groups


class EntrySet(collections.abc.MutableSet, MSONable):
    """A convenient container for manipulating entries. Allows for generating
    subsets, dumping into files, etc.
    """

    def __init__(self, entries: Iterable[PDEntry | ComputedEntry | ComputedStructureEntry]):
        """
        Args:
            entries: All the entries.
        """
        self.entries = set(entries)

    def __contains__(self, item):
        return item in self.entries

    def __iter__(self):
        return iter(self.entries)

    def __len__(self):
        return len(self.entries)

    def add(self, element):
        """Add an entry.

        Args:
            element: Entry
        """
        self.entries.add(element)

    def discard(self, element):
        """Discard an entry.

        Args:
            element: Entry
        """
        self.entries.discard(element)

    @property
    def chemsys(self) -> set:
        """
        Returns:
            set representing the chemical system, e.g. {"Li", "Fe", "P", "O"}.
        """
        chemsys = set()
        for e in self.entries:
            chemsys.update([el.symbol for el in e.composition])
        return chemsys

    @property
    def ground_states(self) -> set:
        """A set containing only the entries that are ground states, i.e., the lowest energy
        per atom entry at each composition.
        """
        entries = sorted(self.entries, key=lambda e: e.reduced_formula)
        return {
            min(g, key=lambda e: e.energy_per_atom)
            for _, g in itertools.groupby(entries, key=lambda e: e.reduced_formula)
        }

    def remove_non_ground_states(self):
        """Removes all non-ground state entries, i.e., only keep the lowest energy
        per atom entry at each composition.
        """
        self.entries = self.ground_states

    def is_ground_state(self, entry) -> bool:
        """Boolean indicating whether a given Entry is a ground state."""
        return entry in self.ground_states

    def get_subset_in_chemsys(self, chemsys: list[str]):
        """Get an EntrySet containing only the set of entries belonging to
        a particular chemical system (in this definition, it includes all sub
        systems). For example, if the entries are from the
        Li-Fe-P-O system, and chemsys=["Li", "O"], only the Li, O,
        and Li-O entries are returned.

        Args:
            chemsys: Chemical system specified as list of elements. e.g.
                ["Li", "O"]

        Returns:
            EntrySet
        """
        chem_sys = set(chemsys)
        if not chem_sys.issubset(self.chemsys):
            raise ValueError(
                f"{sorted(chem_sys)} is not a subset of {sorted(self.chemsys)}, extra: {chem_sys - self.chemsys}"
            )
        subset = set()
        for e in self.entries:
            elements = [sp.symbol for sp in e.composition]
            if chem_sys.issuperset(elements):
                subset.add(e)
        return EntrySet(subset)

    def as_dict(self) -> dict[Literal["entries"], list[Entry]]:
        """Get MSONable dict."""
        return {"entries": list(self.entries)}

    def to_csv(self, filename: str, latexify_names: bool = False) -> None:
        """Exports PDEntries to a csv.

        Args:
            filename: Filename to write to.
            entries: PDEntries to export.
            latexify_names: Format entry names to be LaTex compatible,
                e.g. Li_{2}O
        """
        els: set[SpeciesLike] = set()
        for entry in self.entries:
            els.update(entry.elements)
        elements = sorted(els, key=lambda a: a.X)
        with open(filename, mode="w", encoding="utf-8") as file:
            writer = csv.writer(
                file,
                delimiter=",",
                quotechar='"',
                quoting=csv.QUOTE_MINIMAL,
            )
            writer.writerow(["Name"] + [el.symbol for el in elements] + ["Energy"])
            for entry in self.entries:
                row: list[str] = [(entry.name if not latexify_names else re.sub(r"([0-9]+)", r"_{\1}", entry.name))]
                row.extend([str(entry.composition[el]) for el in elements])
                row.append(str(entry.energy))
                writer.writerow(row)

    @classmethod
    def from_csv(cls, filename: str) -> Self:
        """Imports PDEntries from a csv.

        Args:
            filename: Filename to import from.

        Returns:
            List of Elements, List of PDEntries
        """
        with open(filename, encoding="utf-8") as file:
            reader = csv.reader(file, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL)
            entries = []
            header_read = False
            elements: list[str] = []
            for row in reader:
                if not header_read:
                    elements = row[1 : (len(row) - 1)]
                    header_read = True
                else:
                    name = row[0]
                    energy = float(row[-1])
                    comp = {}
                    for ind in range(1, len(row) - 1):
                        if float(row[ind]) > 0:
                            comp[Element(elements[ind - 1])] = float(row[ind])
                    entries.append(PDEntry(Composition(comp), energy, name))
        return cls(entries)
