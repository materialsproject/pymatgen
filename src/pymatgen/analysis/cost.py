"""
This module is used to estimate the cost of various compounds. Costs are taken
from the a CostDB instance, for example a CSV file via CostDBCSV.
For compounds with no cost listed, a Phase Diagram style convex hull
optimization is performed to determine a set of compositions that can be mixed
to give the desired compound with lowest total cost.
"""

from __future__ import annotations

import abc
import csv
import itertools
import os
from collections import defaultdict
from typing import TYPE_CHECKING

import scipy.constants as const
from monty.design_patterns import singleton

from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.core import Composition, Element
from pymatgen.util.provenance import is_valid_bibtex

if TYPE_CHECKING:
    from pymatgen.util.typing import CompositionLike

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Aug 27, 2013"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


class CostEntry(PDEntry):
    """Extends PDEntry to include a BibTeX reference and include language about cost."""

    def __init__(self, composition, cost, name, reference):
        """
        Args:
            composition (Composition): chemical composition of the entry
            cost (float): per mol, NOT per kg of the full Composition
            name (str): Optional parameter to name the entry. Defaults to the reduced
                chemical formula as in PDEntry.
            reference (str): Reference data as BiBTeX string.
        """
        super().__init__(composition, cost, name)
        if reference and not is_valid_bibtex(reference):
            raise ValueError("Invalid format for cost reference! Should be BibTeX string.")
        self.reference = reference

    def __repr__(self):
        return f"CostEntry : {self.composition} with cost = {self.energy:.4f}"


class CostDB(abc.ABC):
    """
    Abstract class for representing a Cost database.
    Can be extended, e.g. for file-based or REST-based databases.
    """

    @abc.abstractmethod
    def get_entries(self, chemsys):
        """For a given chemical system, return an array of CostEntries.

        Args:
            chemsys (list[SpeciesLike]): Elements defining the chemical system.

        Returns:
            list[CostEntries]
        """


class CostDBCSV(CostDB):
    """Read a CSV file to get costs. Format is formula,cost_per_kg,name,BibTeX."""

    def __init__(self, filename):
        """
        Args:
            filename (str): Filename of cost database.
        """
        # read in data from file
        self._chemsys_entries = defaultdict(list)
        filename = os.path.join(os.path.dirname(__file__), filename)
        with open(filename, encoding="utf-8") as file:
            reader = csv.reader(file, quotechar="|")
            for row in reader:
                comp = Composition(row[0])
                cost_per_mol = float(row[1]) * comp.weight.to("kg") * const.N_A
                pde = CostEntry(comp.formula, cost_per_mol, row[2], row[3])
                chemsys = "-".join(sorted(el.symbol for el in pde.elements))
                self._chemsys_entries[chemsys].append(pde)

    def get_entries(self, chemsys):
        """For a given chemical system, return an array of CostEntries.

        Args:
            chemsys (list[Element]): Elements defining the chemical system.

        Returns:
            array of CostEntries
        """
        chemsys = "-".join(sorted(el.symbol for el in chemsys))
        return self._chemsys_entries[chemsys]


@singleton
class CostDBElements(CostDBCSV):
    """Singleton that provides the cost data for elements."""

    def __init__(self):
        CostDBCSV.__init__(self, f"{MODULE_DIR}/costdb_elements.csv")


class CostAnalyzer:
    """Given a CostDB, figures out the minimum cost solutions via convex hull."""

    def __init__(self, costdb: CostDB) -> None:
        """
        Args:
            costdb (CostDB): Cost database to use.
        """
        self.costdb = costdb

    def get_lowest_decomposition(self, composition):
        """Get the decomposition leading to lowest cost.

        Args:
            composition:
                Composition as a pymatgen.core.structure.Composition

        Returns:
            Decomposition as a dict of {Entry: amount}
        """
        entries = []
        elements = [e.symbol for e in composition.elements]
        for idx in range(len(elements)):
            for combi in itertools.combinations(elements, idx + 1):
                chemsys = [Element(el) for el in combi]
                x = self.costdb.get_entries(chemsys)
                entries.extend(x)
        try:
            pd = PhaseDiagram(entries)
            return pd.get_decomposition(composition)
        except IndexError:
            raise ValueError("Error during PD building; most likely, cost data does not exist!")

    def get_cost_per_mol(self, comp: CompositionLike) -> float:
        """Get best estimate of minimum cost/mol based on known data.

        Args:
            comp (CompositionLike): chemical formula

        Returns:
            float: energy cost/mol
        """
        comp = Composition(comp)
        decomp = self.get_lowest_decomposition(comp)
        return sum(elem.energy_per_atom * val * comp.num_atoms for elem, val in decomp.items())

    def get_cost_per_kg(self, comp):
        """Get best estimate of minimum cost/kg based on known data.

        Args:
            comp (CompositionLike): chemical formula

        Returns:
            float: energy cost/kg
        """
        comp = comp if isinstance(comp, Composition) else Composition(comp)
        return self.get_cost_per_mol(comp) / (comp.weight.to("kg") * const.N_A)
