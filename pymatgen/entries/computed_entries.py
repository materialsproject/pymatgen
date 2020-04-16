# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements equivalents of the basic ComputedEntry objects, which
is the basic entity that can be used to perform many analyses. ComputedEntries
contain calculated information, typically from VASP or other electronic
structure codes. For example, ComputedEntries can be used as inputs for phase
diagram analysis.
"""

import json

from monty.json import MontyEncoder, MontyDecoder

from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.entries import Entry

__author__ = "Shyue Ping Ong, Anubhav Jain"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Apr 30, 2012"


class ComputedEntry(Entry):
    """
    Lightweight Entry object for computed data. Contains facilities 
    for applying corrections to the .energy attribute and for storing 
    calculation parameters.

    """

    def __init__(self,
                 composition: Composition,
                 energy: float,
                 correction: float = 0.0,
                 parameters: dict = None,
                 data: dict = None,
                 entry_id: object = None):
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
        """
        super().__init__(composition, energy)
        self.uncorrected_energy = self._energy
        self.correction = correction
        self.parameters = parameters if parameters else {}
        self.data = data if data else {}
        self.entry_id = entry_id
        self.name = self.composition.reduced_formula

    @property
    def energy(self) -> float:
        """
        :return: the *corrected* energy of the entry.
        """
        return self._energy + self.correction

    def normalize(self, mode: str = "formula_unit") -> None:
        """
        Normalize the entry's composition and energy.

        Args:
            mode: "formula_unit" is the default, which normalizes to
                composition.reduced_formula. The other option is "atom", which
                normalizes such that the composition amounts sum to 1.
        """
        factor = self._normalization_factor(mode)
        self.correction /= factor
        self.uncorrected_energy /= factor
        if 'energy_adjustments' in self.data.keys():
            for k1, v1 in self.data["energy_adjustments"].items():
                for k2, v2 in v1.items():
                    self.data["energy_adjustments"][k1][k2] /= factor
        super().normalize(mode)

    def __repr__(self):
        output = ["ComputedEntry {} - {}".format(self.entry_id,
                                                 self.composition.formula),
                  "Energy = {:.4f}".format(self._energy),
                  "Correction = {:.4f}".format(self.correction),
                  "Parameters:"]
        for k, v in self.parameters.items():
            output.append("{} = {}".format(k, v))
        output.append("Data:")
        for k, v in self.data.items():
            output.append("{} = {}".format(k, v))
        return "\n".join(output)

    @classmethod
    def from_dict(cls, d) -> 'ComputedEntry':
        """
        :param d: Dict representation.
        :return: ComputedEntry
        """
        dec = MontyDecoder()
        return cls(d["composition"], d["energy"], d["correction"],
                   parameters={k: dec.process_decoded(v)
                               for k, v in d.get("parameters", {}).items()},
                   data={k: dec.process_decoded(v)
                         for k, v in d.get("data", {}).items()},
                   entry_id=d.get("entry_id", None))

    def as_dict(self) -> dict:
        """
        :return: MSONable dict.
        """
        return_dict = super().as_dict()
        return_dict.update({"parameters": json.loads(json.dumps(self.parameters, cls=MontyEncoder)),
                            "data": json.loads(json.dumps(self.data, cls=MontyEncoder)),
                            "entry_id": self.entry_id,
                            "correction": self.correction})
        return return_dict


class ComputedStructureEntry(ComputedEntry):
    """
    A heavier version of ComputedEntry which contains a structure as well. The
    structure is needed for some analyses.
    """

    def __init__(self,
                 structure: Structure,
                 energy: float,
                 correction: float = 0.0,
                 parameters: dict = None,
                 data: dict = None,
                 entry_id: object = None):
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
        super().__init__(
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

    def as_dict(self) -> dict:
        """
        :return: MSONAble dict.
        """
        d = super().as_dict()
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["structure"] = self.structure.as_dict()
        return d

    @classmethod
    def from_dict(cls, d) -> 'ComputedStructureEntry':
        """
        :param d: Dict representation.
        :return: ComputedStructureEntry
        """
        dec = MontyDecoder()
        return cls(dec.process_decoded(d["structure"]),
                   d["energy"], d["correction"],
                   parameters={k: dec.process_decoded(v)
                               for k, v in d.get("parameters", {}).items()},
                   data={k: dec.process_decoded(v)
                         for k, v in d.get("data", {}).items()},
                   entry_id=d.get("entry_id", None))
