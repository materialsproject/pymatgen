# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import json

from monty.json import MontyEncoder, MontyDecoder

from pymatgen.core.composition import Composition
from monty.json import MSONable

"""
This module implements equivalents of the basic ComputedEntry objects, which
is the basic entity that can be used to perform many analyses. ComputedEntries
contain calculated information, typically from VASP or other electronic
structure codes. For example, ComputedEntries can be used as inputs for phase
diagram analysis.
"""

__author__ = "Shyue Ping Ong, Anubhav Jain"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Apr 30, 2012"


class ComputedEntry(MSONable):
    """
    An lightweight ComputedEntry object containing key computed data
    for many purposes. Extends a PDEntry so that it can be used for phase
    diagram generation. The difference between a ComputedEntry and a standard
    PDEntry is that it includes additional parameters like a correction and
    run_parameters.

    """

    def __init__(self, composition, energy, correction=0.0, parameters=None,
                 data=None, entry_id=None, attribute=None):
        """
        Initializes a ComputedEntry.

        Args:
            composition (Composition): Composition of the entry. For
                flexibility, this can take the form of all the typical input
                taken by a Composition, including a {symbol: amt} dict,
                a string formula, and others.
            energy (float): Energy of the entry. Usually the final calculated
                energy from VASP or other electronic structure codes.
            correction (float): A correction to be applied to the energy.
                This is used to modify the energy for certain analyses.
                Defaults to 0.0.
            parameters (dict): An optional dict of parameters associated with
                the entry. Defaults to None.
            data (dict): An optional dict of any additional data associated
                with the entry. Defaults to None.
            entry_id (obj): An optional id to uniquely identify the entry.
            attribute: Optional attribute of the entry. This can be used to
                specify that the entry is a newly found compound, or to specify
                a particular label for the entry, or else ... Used for further
                analysis and plotting purposes. An attribute can be anything
                but must be MSONable.
        """
        self.uncorrected_energy = energy
        self.composition = Composition(composition)
        self.correction = correction
        self.parameters = parameters if parameters else {}
        self.data = data if data else {}
        self.entry_id = entry_id
        self.name = self.composition.reduced_formula
        self.attribute = attribute

    @property
    def is_element(self):
        return self.composition.is_element

    @property
    def energy(self):
        """
        Returns the *corrected* energy of the entry.
        """
        return self.uncorrected_energy + self.correction

    @property
    def energy_per_atom(self):
        return self.energy / self.composition.num_atoms

    def __repr__(self):
        output = ["ComputedEntry {} - {}".format(self.entry_id,
                                                 self.composition.formula),
                  "Energy = {:.4f}".format(self.uncorrected_energy),
                  "Correction = {:.4f}".format(self.correction),
                  "Parameters:"]
        for k, v in self.parameters.items():
            output.append("{} = {}".format(k, v))
        output.append("Data:")
        for k, v in self.data.items():
            output.append("{} = {}".format(k, v))
        return "\n".join(output)

    def __str__(self):
        return self.__repr__()

    @classmethod
    def from_dict(cls, d):
        dec = MontyDecoder()
        return cls(d["composition"], d["energy"], d["correction"],
                   parameters={k: dec.process_decoded(v)
                               for k, v in d.get("parameters", {}).items()},
                   data={k: dec.process_decoded(v)
                         for k, v in d.get("data", {}).items()},
                   entry_id=d.get("entry_id", None),
                   attribute=d["attribute"] if "attribute" in d else None)

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "energy": self.uncorrected_energy,
                "composition": self.composition.as_dict(),
                "correction": self.correction,
                "parameters": json.loads(json.dumps(self.parameters,
                                                    cls=MontyEncoder)),
                "data": json.loads(json.dumps(self.data, cls=MontyEncoder)),
                "entry_id": self.entry_id,
                "attribute": self.attribute}


class ComputedStructureEntry(ComputedEntry):
    """
    A heavier version of ComputedEntry which contains a structure as well. The
    structure is needed for some analyses.
    """

    def __init__(self, structure, energy, correction=0.0, parameters=None,
                 data=None, entry_id=None):
        """
        Initializes a ComputedStructureEntry.

        Args:
            structure (Structure): The actual structure of an entry.
            energy (float): Energy of the entry. Usually the final calculated
                energy from VASP or other electronic structure codes.
            correction (float): A correction to be applied to the energy.
                This is used to modify the energy for certain analyses.
                Defaults to 0.0.
            parameters (dict): An optional dict of parameters associated with
                the entry. Defaults to None.
            data (dict): An optional dict of any additional data associated
                with the entry. Defaults to None.
            entry_id (obj): An optional id to uniquely identify the entry.
        """
        super(ComputedStructureEntry, self).__init__(
            structure.composition, energy, correction=correction,
            parameters=parameters, data=data, entry_id=entry_id)
        self.structure = structure

    def __repr__(self):
        output = ["ComputedStructureEntry {} - {}".format(
            self.entry_id, self.composition.formula),
                  "Energy = {:.4f}".format(self.uncorrected_energy),
                  "Correction = {:.4f}".format(self.correction), "Parameters:"]
        for k, v in self.parameters.items():
            output.append("{} = {}".format(k, v))
        output.append("Data:")
        for k, v in self.data.items():
            output.append("{} = {}".format(k, v))
        return "\n".join(output)

    def __str__(self):
        return self.__repr__()

    def as_dict(self):
        d = super(ComputedStructureEntry, self).as_dict()
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["structure"] = self.structure.as_dict()
        return d

    @classmethod
    def from_dict(cls, d):
        dec = MontyDecoder()
        return cls(dec.process_decoded(d["structure"]),
                   d["energy"], d["correction"],
                   parameters={k: dec.process_decoded(v)
                               for k, v in d.get("parameters", {}).items()},
                   data={k: dec.process_decoded(v)
                         for k, v in d.get("data", {}).items()},
                   entry_id=d.get("entry_id", None))
