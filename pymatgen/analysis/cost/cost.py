# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module is used to estimate the cost of various compounds. Costs are taken
from the a CostDB instance, for example a CSV file via CostDBCSV.
For compounds with no cost listed, a Phase Diagram style convex hull
optimization is performed to determine a set of compositions that can be mixed
to give the desired compound with lowest total cost.
"""

from __future__ import division, unicode_literals
import abc
from collections import defaultdict
import csv
import os
import itertools
from monty.design_patterns import singleton
from monty.string import unicode2str
import six
from pymatgen import Composition, Element
from pymatgen.core.physical_constants import AVOGADROS_CONST
from pymatgen.matproj.snl import is_valid_bibtex
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.phasediagram.pdanalyzer import PDAnalyzer
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from io import open


__author__ = 'Anubhav Jain'
__copyright__ = 'Copyright 2013, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'
__date__ = 'Aug 27, 2013'


module_dir = os.path.dirname(os.path.abspath(__file__))


class CostEntry(PDEntry):
    """
    Extends PDEntry to include a BibTeX reference and include language about
    cost
    """

    def __init__(self, composition, cost, name, reference):
        """
        Args:
            composition:
                Composition as a pymatgen.core.structure.Composition
            cost:
                Cost (per mol, NOT per kg) of the full Composition
            name:
                Optional parameter to name the entry. Defaults to the reduced
                chemical formula as in PDEntry.
            reference:
                Reference data as BiBTeX string
        """
        super(CostEntry, self).__init__(composition, cost, name)
        if reference and not is_valid_bibtex(reference):
            raise ValueError(
                "Invalid format for cost reference! Should be BibTeX string.")
        self.reference = reference

    def __repr__(self):
        return "CostEntry : {} with cost = {:.4f}".format(self.composition,
                                                          self.energy)


class CostDB(six.with_metaclass(abc.ABCMeta)):
    """
    Abstract class for representing a Cost database.
    Can be extended, e.g. for file-based or REST-based databases
    """

    @abc.abstractmethod
    def get_entries(self, chemsys):
        """
        For a given chemical system, return an array of CostEntries

        Args:
            chemsys:
                array of Elements defining the chemical system.

        Returns:
            array of CostEntries
        """
        return


class CostDBCSV(CostDB):
    """
    Read a CSV file to get costs
    Format is formula,cost_per_kg,name,BibTeX
    """

    def __init__(self, filename):
        # read in data from file
        self._chemsys_entries = defaultdict(list)
        filename = os.path.join(os.path.dirname(__file__), filename)
        reader = csv.reader(open(filename, "rt"), quotechar=unicode2str("|"))
        for row in reader:
            comp = Composition(row[0])
            cost_per_mol = float(row[1]) * comp.weight.to("kg") * \
                           AVOGADROS_CONST
            pde = CostEntry(comp.formula, cost_per_mol, row[2], row[3])
            chemsys = "-".join(sorted([el.symbol
                                       for el in pde.composition.elements]))
            self._chemsys_entries[chemsys].append(pde)

    def get_entries(self, chemsys):
        chemsys = "-".join(sorted([el.symbol for el in chemsys]))
        return self._chemsys_entries[chemsys]


@singleton
class CostDBElements(CostDBCSV):
    """
    Singleton object that provides the cost data for elements
    """
    def __init__(self):
        CostDBCSV.__init__(
            self, os.path.join(module_dir, "costdb_elements.csv"))


class CostAnalyzer(object):
    """
    Given a CostDB, figures out the minimum cost solutions via convex hull
    """

    def __init__(self, costdb):
        self.costdb = costdb

    def get_lowest_decomposition(self, composition):
        """
        Get the decomposition leading to lowest cost

        Args:
            composition:
                Composition as a pymatgen.core.structure.Composition
        Returns:
            Decomposition as a dict of {Entry: amount}
        """

        entries_list = []
        elements = [e.symbol for e in composition.elements]
        for i in range(len(elements)):
            for combi in itertools.combinations(elements, i + 1):
                chemsys = [Element(e) for e in combi]
                x = self.costdb.get_entries(chemsys)
                entries_list.extend(x)
        try:
            pd = PhaseDiagram(entries_list)
            return PDAnalyzer(pd).get_decomposition(composition)
        except IndexError:
            raise ValueError("Error during PD building; most likely, "
                             "cost data does not exist!")

    def get_cost_per_mol(self, comp):
        """
        Get best estimate of minimum cost/mol based on known data

        Args:
            comp:
                Composition as a pymatgen.core.structure.Composition
        Returns:
            float of cost/mol
        """

        comp = comp if isinstance(comp, Composition) else Composition(comp)
        decomp = self.get_lowest_decomposition(comp)
        return sum(k.energy_per_atom * v * comp.num_atoms for k, v in
                   decomp.items())

    def get_cost_per_kg(self, comp):
        """
        Get best estimate of minimum cost/kg based on known data

        Args:
            comp:
                Composition as a pymatgen.core.structure.Composition
        Returns:
            float of cost/kg
        """
        comp = comp if isinstance(comp, Composition) else Composition(comp)
        return self.get_cost_per_mol(comp) / (
            comp.weight.to("kg") * AVOGADROS_CONST)