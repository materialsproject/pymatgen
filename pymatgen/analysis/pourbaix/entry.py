#!/usr/bin/env python

"""
Module which defines basic entries for each ion and oxide to compute a
Pourbaix diagram
"""

from __future__ import division

__author__ = "Sai Jayaraman"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.0"
__maintainer__ = "Sai Jayaraman"
__email__ = "sjayaram@mit.edu"
__status__ = "Production"
__date__ = "December 10, 2012"

import re

from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Composition
from pymatgen.serializers.json_coders import MSONable
from pymatgen.core.ion import Ion
from pymatgen.phasediagram.entries import PDEntry


class PourbaixEntry(MSONable):
    """
    An object encompassing all data relevant to an ion in a pourbaix diagram.
    Each bulk solid/ion has a free energy g of the form:
    g = g0_ref + 0.0591 log10(concn) - nO mu_H2O + (nH - 2nO) pH + phi (-nH+2nO+q)
    """
    def __init__(self, entry, g=None):
        """
        Args:
            entry:
                Either an Ion object or a Composition object
            energy:
                Energy of entry
        """
        if (entry.__class__.__name__ == "IonEntry"):
            self._entry = entry
            self._concn = 1.0e-6
            self._phase_type = "Ion"
            self._charge = entry.composition.charge
        else:
            self._entry = entry
            self._concn = 1.0
            self._phase_type = "Solid"
            self._charge = 0.0
        self._npH = None
        self._nPhi = None
        self._nH2O = None
        self._nM = None
        if g is not None:
            self._g0 = g
        else:
            self._g0 = entry.energy
        self._calc_coeff_terms()
        self._name = self._entry.composition.reduced_formula

    @property
    def name(self):
        """
        Returns the entry's name
        """
        return self._name

    def set_name(self, string):
        """
        Set name of entry
        args:
            string: Input string
        """
        self._name = string

    @property
    def npH(self):
        """
        Returns value of npH, the coefficient of pH
        """
        return self._npH

    @property
    def nH2O(self):
        """
        Returns coefficient of Mu_H2O
        """
        return self._nH2O

    @property
    def nPhi(self):
        """
        Returns nPhi, the coefficient of Phi
        """
        return self._nPhi

    @property
    def g0(self):
        """
        Return g0 for the entry
        """
        return self._g0

    @property
    def conc(self):
        """
        Return concentration of the entry. Returns 1 if solid.
        """
        return self._concn

    @property
    def phase_type(self):
        """
        Returns whether the entry is a solid/ion
        """
        return self._phase_type

    def g0_add(self, term):
        """
        Add a correction term to g0
        args: 
            term: 
                Correction term to add to g0
        """
        self._g0 += term

    def g0_replace(self, term):
        """
        Replace g0 by a different value
        args: 
            term: 
                New value for g0
        """
        self._g0 = term

    @property
    def to_dict(self):
        """
        Returns dict which contains Pourbaix Entry data.
        Note that the pH, voltage, H2O factors are always calculated when
        constructing a PourbaixEntry object.
        """
        d = {}
        d["module"] = self.__class__.__module__
        d["class"] = self.__class__.__name__
        if (self._entry.__class__.__name__ == "IonEntry"):
            d["entry type"] = "Ion"
        else:
            d["entry type"] = "Solid"
        d["entry"] = self._entry.to_dict
        d["pH factor"] = self._npH
        d["voltage factor"] = self._nPhi
        d["concentration"] = self._concn
        d["H2O factor"] = self._nH2O
        d["free energy"] = self._g0
        return d

    @staticmethod
    def from_dict(d):
        """
        Returns a PourbaixEntry by reading in an Ion
        """
        entry_type = d["entry type"]
        if (entry_type == "Ion"):
            entry = IonEntry.from_dict(d["entry"])
        else:
            entry = PDEntry.from_dict(d["entry"])
        g = d["free energy"]
        return PourbaixEntry(entry, g)

    def _calc_coeff_terms(self):
        """
        Calculates coefficients of pH, V, H2O
        """
        nH = 0
        nO = 0
        for elt in self._entry.composition.elements:
            if elt == (Element("H")):
                nH = self.entry.composition[elt]
            elif elt == (Element("O")):
                nO = self.entry.composition[elt]
            else:
                self._nM = self.entry.composition[elt]
        self._npH = (nH - 2 * nO)
        self._nH2O = nO
        self._nPhi = (nH - 2 * nO - self._charge)

    @property
    def normalization_factor(self):
        """
        Normalize each entry by nM
        """
        fact = 1.0 / self._nM
        return fact

    def normalize(self, factor):
        """
        Normalize all entries by normalization factor
        args:
            factor:
                Normalization factor
        """
        self._npH *= factor
        self._nPhi *= factor
        self._nH2O *= factor
        self._g0 *= factor

    @property
    def charge(self):
        """
        Returns charge of entry
        """
        return self._charge

    @property
    def entry(self):
        """
        Returns IonEntry/PDEntry object
        """
        return self._entry

    def reduced_entry(self):
        """
        Phase Diagram returns solids which are not in their most reduced state.
        This, if called, will divide all relevant value by reduction factor
        """
        reduction_factor = self.entry.composition.\
                            get_reduced_composition_and_factor()[1]
        self._g0 /= reduction_factor
        self._npH /= reduction_factor
        self._nH2O /= reduction_factor
        self._nPhi /= reduction_factor
        self._nM /= reduction_factor

    @property
    def num_atoms(self):
        """
        Return number of atoms in current formula. Useful for normalization
        """
        return self.entry.composition.num_atoms\
            / self.entry.composition.get_reduced_composition_and_factor()[1]

    def set_conc(self, concn):
        """
        Set concentration manually
        args:
            concn:
                Input concentration
        """
        self._concn = concn

    def __repr__(self):
        return "Pourbaix Entry : {} with energy = {:.4f}, npH = {}, nPhi = {},\
             nH2O = {}".format(self._entry.composition, self.g0, self.npH,\
                                self.nPhi, self.nH2O)

    def __str__(self):
        return self.__repr__()


class IonEntry(PDEntry):
    """
    Object similar to PDEntry, but contains an Ion object instead of a Composition
    object.
    .. attribute:: name
        A name for the entry. This is the string shown in the phase diagrams.
        By default, this is the reduced formula for the composition, but can be
        set to some other string for display purposes.
    """

    def __init__(self, ion, energy, name=None):
        """
        Args:
            comp:
                Ion object
            energy:
                Energy for composition.
            name:
                Optional parameter to name the entry. Defaults to the
                chemical formula.
        """
        self._energy = energy
        self._composition = ion
        self.name = name if name else self._composition.reduced_formula

    @staticmethod
    def from_dict(d):
        """
        Returns an IonEntry object from a dict.
        """
        return IonEntry(Ion.from_dict(d["composition"]), d["energy"])

    @property
    def to_dict(self):
        """
        Creates a dict of composition, energy, and ion name
        """
        d = {}
        d["composition"] = self._composition.to_dict
        d["energy"] = self._energy
        return d

    @property
    def energy(self):
        """
        Return final energy
        """
        return self._energy

    @property
    def energy_per_atom(self):
        """
        Return final energy per atom
        """
        return self._energy / self.composition.num_atoms

    @property
    def composition(self):
        """
        Returns the composition
        """
        return self._composition

    def __repr__(self):
        return "PDEntry : {} with energy = {:.4f}".format(self.composition,
                                                          self.energy)

    def __str__(self):
        return self.__repr__()


class PourbaixEntryIO(object):
    """
    Class to import and export Pourbaix entries from a csv file
    """
    @staticmethod
    def to_csv(filename, entries, latexify_names = False):
        """
        Exports Pourbaix entries to a csv

        Args:
            filename:
                Filename to write to.
            entries:
                Entries to export.
            latexify_names:
                Format entry names to be LaTex compatible, e.g., Li_{2}O
        """
        import csv
        elements = set()
        map(elements.update, [entry.entry.composition.elements for entry in\
                               entries])
        elements = sorted(list(elements), key=lambda a: a.X)
        writer = csv.writer(open(filename, "wb"), delimiter=",",
                            quotechar="\"", quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["Name"] + elements + ["Energy"] + ["Entry Type"]\
                         + ["Charge"] + ["Concentration"])
        for entry in entries:
            row = [entry.name if not latexify_names
                   else re.sub(r"([0-9]+)", r"_{\1}", entry.name)]
            if entry.phase_type == "Solid":
                reduction_fac = entry.entry.composition.\
                                get_reduced_composition_and_factor()[1]
            else:
                reduction_fac = 1.0
            row.extend([entry.entry.composition[el] / reduction_fac\
                            for el in elements])
            if entry.phase_type == "Solid":
                reduction_fac = 1.0
            row.append(entry.g0 / reduction_fac)
            row.append(entry.phase_type)
            row.append(entry.charge / reduction_fac)
            row.append(entry.conc)
            writer.writerow(row)

    @staticmethod
    def from_csv(filename):
        """
        Imports PourbaixEntries from a csv

        Args:
            filename - Filename to import from.

        Returns:
            List of Entries
        """
        import csv
        reader = csv.reader(open(filename, "rb"), delimiter=",",
                            quotechar="\"", quoting=csv.QUOTE_MINIMAL)
        entries = list()
        header_read = False
        for row in reader:
            if not header_read:
                elements = row[1:(len(row) - 4)]
                header_read = True
            else:
                name = row[0]
                energy = float(row[-4])
                conc = float(row[-1])
                comp = dict()
                for ind in range(1, len(row) - 4):
                    if float(row[ind]) > 0:
                        comp[Element(elements[ind - 1])] = float(row[ind])
                phase_type = row[-3]
                if (phase_type == "Ion"):
                    PoE = PourbaixEntry(IonEntry(Ion.from_formula(name),\
                                                  energy))
                    PoE.set_conc(conc)
                    PoE.set_name(name)
                    entries.append(PoE)
                else:
                    entries.append(PourbaixEntry(PDEntry(Composition(comp),\
                                                  energy)))
        elements = [Element(el) for el in elements]
        return (elements, entries)
