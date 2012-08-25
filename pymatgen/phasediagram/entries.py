#!/usr/bin/env python

"""
This module defines PDEntry, which wraps information (composition and energy)
necessary to create phase diagrams.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ = "May 16, 2011"

import re

from pymatgen.core.structure import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.serializers.json_coders import MSONable, PMGJSONDecoder


class PDEntry(MSONable):
    """
    An object encompassing all relevant data for phase diagrams.

    .. attribute:: name

        A name for the entry. This is the string shown in the phase diagrams.
        By default, this is the reduced formula for the composition, but can be
        set to some other string for display purposes.
    """

    def __init__(self, composition, energy, name=None):
        """
        Args:
            comp:
                Composition as a pymatgen.core.structure.Composition
            energy:
                Energy for composition.
            name:
                Optional parameter to name the entry. Defaults to the reduced
                chemical formula.
        """
        self._energy = energy
        self._composition = Composition(composition)
        self.name = name if name else self._composition.reduced_formula

    @property
    def energy(self):
        """
        Returns the final energy.
        """
        return self._energy

    @property
    def energy_per_atom(self):
        """
        Returns the final energy per atom.
        """
        return self.energy / self.composition.num_atoms

    @property
    def composition(self):
        """
        Returns the composition.
        """
        return self._composition

    @property
    def is_element(self):
        """
        True if the entry is an element.
        """
        return self._composition.is_element

    def __repr__(self):
        return "PDEntry : {} with energy = {:.4f}".format(self.composition,
                                                          self.energy)

    def __str__(self):
        return self.__repr__()

    @property
    def to_dict(self):
        d = {}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["composition"] = self._composition.to_dict
        d["energy"] = self._energy
        d["name"] = self.name
        return d

    @staticmethod
    def from_dict(d):
        return PDEntry(Composition(d["composition"]), d["energy"], d["name"])


class GrandPotPDEntry(PDEntry):
    """
    A grand potential pd entry object encompassing all relevant data for phase
    diagrams.  Chemical potentials are given as a element-chemical potential
    dict.
    """
    def __init__(self, entry, chempots, name=None):
        """
        Args:
            entry:
                A PDEntry-like object.
            chempots:
                Chemical potential specification as {Element: float}.
            name:
                Optional parameter to name the entry. Defaults to the reduced
                chemical formula of the original entry.
        """
        comp = entry.composition
        self._original_entry = entry
        self._original_comp = comp
        grandpot = entry.energy - sum([comp[el] * pot
                                       for el, pot in chempots.items()])
        self.chempots = chempots
        new_comp_map = dict()
        for el in comp.elements:
            if el not in chempots:
                new_comp_map[el] = comp[el]
        newcomposition = Composition(new_comp_map)
        super(GrandPotPDEntry, self).__init__(newcomposition, grandpot,
                                              entry.name)
        self.name = name if name else entry.name

    @property
    def original_entry(self):
        """
        Original entry.
        """
        return self._original_entry

    @property
    def is_element(self):
        """
        True if the entry is an element.
        """
        return self._original_comp.is_element

    def __repr__(self):
        chempot_str = " ".join(["mu_%s = %.4f" % (el, mu)
                                for el, mu in self.chempots.items()])
        return "GrandPotPDEntry with original composition " + \
            "{}, energy = {:.4f}, {}".format(self.original_entry.composition,
                                             self.original_entry.energy,
                                             chempot_str)

    def __str__(self):
        return self.__repr__()

    @property
    def to_dict(self):
        d = {}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["entry"] = self._original_entry.to_dict
        d["chempots"] = {el.symbol: u for el, u in self.chempots.items()}
        d["name"] = self.name
        return d

    @staticmethod
    def from_dict(d):
        chempots = {Element(symbol): u for symbol, u in d["chempots"].items()}
        entry = PMGJSONDecoder().process_decoded(d["entry"])
        return GrandPotPDEntry(entry, chempots, d["name"])

    def __getattr__(self, a):
        """
        Delegate attribute to original entry if available.
        """
        if hasattr(self._original_entry, a):
            return getattr(self._original_entry, a)
        raise AttributeError(a)


class PDEntryIO(object):
    """
    Utility class to export and import PDEntry to and from csv files, as well
    as to and from json.
    """

    @staticmethod
    def to_csv(filename, entries, latexify_names=False):
        """
        Exports PDEntries to a csv

        Args:
            filename:
                Filename to write to.
            entries:
                PDEntries to export.
            latexify_names:
                Format entry names to be LaTex compatible, e.g., Li_{2}O
        """
        import csv
        elements = set()
        map(elements.update, [entry.composition.elements for entry in entries])
        elements = sorted(list(elements), key=lambda a: a.X)
        writer = csv.writer(open(filename, "wb"), delimiter=",",
                            quotechar="\"", quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["Name"] + elements + ["Energy"])
        for entry in entries:
            row = [entry.name if not latexify_names
                   else re.sub(r"([0-9]+)", r"_{\1}", entry.name)]
            row.extend([entry.composition[el] for el in elements])
            row.append(entry.energy)
            writer.writerow(row)

    @staticmethod
    def from_csv(filename):
        """
        Imports PDEntries from a csv

        Args:
            filename - Filename to import from.

        Returns:
            List of PDEntries
        """
        import csv
        reader = csv.reader(open(filename, "rb"), delimiter=",",
                            quotechar="\"", quoting=csv.QUOTE_MINIMAL)
        entries = list()
        header_read = False
        for row in reader:
            if not header_read:
                elements = row[1:(len(row) - 1)]
                header_read = True
            else:
                name = row[0]
                energy = float(row[-1])
                comp = dict()
                for ind in range(1, len(row) - 1):
                    if float(row[ind]) > 0:
                        comp[Element(elements[ind - 1])] = float(row[ind])
                entries.append(PDEntry(Composition(comp), energy, name))
        elements = [Element(el) for el in elements]
        return (elements, entries)


class TransformedPDEntry(PDEntry):
    """
    This class repesents a TransformedPDEntry, which allows for a PDEntry to be
    transformed to a different composition coordinate space. It is used in the
    construction of phase diagrams that do not have elements as the terminal
    compositions.
    """

    def __init__(self, comp, original_entry):
        """
        Args:
            comp:
                Transformed composition as a Composition.
            energy:
                Energy for composition.
            original_entry:
                Original entry that this entry arose from.
        """
        PDEntry.__init__(self, comp, original_entry.energy)
        self._original_entry = original_entry
        self.name = original_entry.name

    @property
    def original_entry(self):
        """
        Original PDEntry object.
        """
        return self._original_entry

    def __getattr__(self, a):
        """
        Delegate attribute to original entry if available.
        """
        if hasattr(self._original_entry, a):
            return getattr(self._original_entry, a)
        raise AttributeError(a)

    def __repr__(self):
        output = ["TransformedPDEntry {}".format(self.composition)]
        output.append(" with original composition {}"
                      .format(self.original_entry.composition))
        output.append(", E = {:.4f}".format(self.original_entry.energy))
        return "".join(output)

    def __str__(self):
        return self.__repr__()

    @property
    def to_dict(self):
        d = {}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["entry"] = self._original_entry.to_dict
        d["composition"] = self.composition
        return d

    @staticmethod
    def from_dict(d):
        entry = PMGJSONDecoder().process_decoded(d["entry"])
        return TransformedPDEntry(d["composition"], entry)
