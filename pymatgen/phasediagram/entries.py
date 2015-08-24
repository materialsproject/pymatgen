# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module defines PDEntry, which wraps information (composition and energy)
necessary to create phase diagrams.
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "May 16, 2011"

import re
import csv

from monty.json import MontyDecoder

from io import open
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.serializers.json_coders import PMGSONable
from monty.string import unicode2str


class PDEntry(PMGSONable):
    """
    An object encompassing all relevant data for phase diagrams.

    .. attribute:: name

        A name for the entry. This is the string shown in the phase diagrams.
        By default, this is the reduced formula for the composition, but can be
        set to some other string for display purposes.

    Args:
        comp: Composition as a pymatgen.core.structure.Composition
        energy: Energy for composition.
        name: Optional parameter to name the entry. Defaults to the reduced
            chemical formula.
        attribute: Optional attribute of the entry. This can be used to
            specify that the entry is a newly found compound, or to specify a
            particular label for the entry, or else ... Used for further
            analysis and plotting purposes. An attribute can be anything
            but must be PMGSONable.
    """

    def __init__(self, composition, energy, name=None, attribute=None):
        self.energy = energy
        self.composition = Composition(composition)
        self.name = name if name else self.composition.reduced_formula
        self.attribute = attribute

    @property
    def energy_per_atom(self):
        """
        Returns the final energy per atom.
        """
        return self.energy / self.composition.num_atoms

    @property
    def is_element(self):
        """
        True if the entry is an element.
        """
        return self.composition.is_element

    def __repr__(self):
        return "PDEntry : {} with energy = {:.4f}".format(self.composition,
                                                          self.energy)

    def __str__(self):
        return self.__repr__()

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "composition": self.composition.as_dict(),
                "energy": self.energy,
                "name": self.name,
                "attribute": self.attribute}

    @classmethod
    def from_dict(cls, d):
        return cls(Composition(d["composition"]), d["energy"], d["name"],
                   d["attribute"] if "attribute" in d else None)


class GrandPotPDEntry(PDEntry):
    """
    A grand potential pd entry object encompassing all relevant data for phase
    diagrams.  Chemical potentials are given as a element-chemical potential
    dict.

    Args:
        entry: A PDEntry-like object.
        chempots: Chemical potential specification as {Element: float}.
        name: Optional parameter to name the entry. Defaults to the reduced
            chemical formula of the original entry.
    """
    def __init__(self, entry, chempots, name=None):
        comp = entry.composition
        self.original_entry = entry
        self.original_comp = comp
        grandpot = entry.energy - sum([comp[el] * pot
                                       for el, pot in chempots.items()])
        self.chempots = chempots
        new_comp_map = {el: comp[el] for el in comp.elements
                        if el not in chempots}
        super(GrandPotPDEntry, self).__init__(new_comp_map, grandpot,
                                              entry.name)
        self.name = name if name else entry.name

    @property
    def is_element(self):
        """
        True if the entry is an element.
        """
        return self.original_comp.is_element

    def __repr__(self):
        chempot_str = " ".join(["mu_%s = %.4f" % (el, mu)
                                for el, mu in self.chempots.items()])
        return "GrandPotPDEntry with original composition " + \
            "{}, energy = {:.4f}, {}".format(self.original_entry.composition,
                                             self.original_entry.energy,
                                             chempot_str)

    def __str__(self):
        return self.__repr__()

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "entry": self.original_entry.as_dict(),
                "chempots": {el.symbol: u for el, u in self.chempots.items()},
                "name": self.name}

    @classmethod
    def from_dict(cls, d):
        chempots = {Element(symbol): u for symbol, u in d["chempots"].items()}
        entry = MontyDecoder().process_decoded(d["entry"])
        return cls(entry, chempots, d["name"])

    def __getattr__(self, a):
        """
        Delegate attribute to original entry if available.
        """
        if hasattr(self.original_entry, a):
            return getattr(self.original_entry, a)
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
            filename: Filename to write to.
            entries: PDEntries to export.
            latexify_names: Format entry names to be LaTex compatible,
                e.g., Li_{2}O
        """
        import csv
        elements = set()
        for entry in entries:
            elements.update(entry.composition.elements)
        elements = sorted(list(elements), key=lambda a: a.X)
        writer = csv.writer(open(filename, "wb"), delimiter=unicode2str(","),
                            quotechar=unicode2str("\""),
                            quoting=csv.QUOTE_MINIMAL)
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
        Imports PDEntries from a csv.

        Args:
            filename: Filename to import from.

        Returns:
            List of Elements, List of PDEntries
        """
        with open(filename, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter=unicode2str(","),
                                quotechar=unicode2str("\""),
                                quoting=csv.QUOTE_MINIMAL)
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
        return elements, entries


class TransformedPDEntry(PDEntry):
    """
    This class repesents a TransformedPDEntry, which allows for a PDEntry to be
    transformed to a different composition coordinate space. It is used in the
    construction of phase diagrams that do not have elements as the terminal
    compositions.

    Args:
        comp: Transformed composition as a Composition.
        energy: Energy for composition.
        original_entry: Original entry that this entry arose from.
    """

    def __init__(self, comp, original_entry):
        super(TransformedPDEntry, self).__init__(comp, original_entry.energy)
        self.original_entry = original_entry
        self.name = original_entry.name

    def __getattr__(self, a):
        """
        Delegate attribute to original entry if available.
        """
        if hasattr(self.original_entry, a):
            return getattr(self.original_entry, a)
        raise AttributeError(a)

    def __repr__(self):
        output = ["TransformedPDEntry {}".format(self.composition),
                  " with original composition {}"
                  .format(self.original_entry.composition),
                  ", E = {:.4f}".format(self.original_entry.energy)]
        return "".join(output)

    def __str__(self):
        return self.__repr__()

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "entry": self.original_entry.as_dict(),
                "composition": self.composition}

    @classmethod
    def from_dict(cls, d):
        entry = MontyDecoder().process_decoded(d["entry"])
        return cls(d["composition"], entry)
