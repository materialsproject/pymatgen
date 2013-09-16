#!/usr/bin/env python

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

import json

from pymatgen.phasediagram.entries import PDEntry
from pymatgen.core.composition import Composition
from pymatgen.serializers.json_coders import MSONable, PMGJSONDecoder, \
    PMGJSONEncoder


class ComputedEntry(PDEntry, MSONable):
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
        Args:
            composition:
                Composition of the entry. For flexibility, this can take the
                form of all the typical input taken by a Composition, including
                a {symbol: amt} dict, a string formula, and others.
            energy:
                Energy of the entry. Usually the final calculated energy from
                VASP or other electronic structure codes.
            correction:
                A correction to be applied to the energy. This is used to
                modify the energy for certain analyses. Defaults to 0.0.
            parameters:
                An optional dict of parameters associated with the entry.
                Defaults to None.
            data:
                An optional dict of any additional data associated with the
                entry. Defaults to None.
            entry_id:
                An optional id to uniquely identify the entry.
            attribute:
                Optional attribute of the entry. This can be used to specify that
                the entry is a newly found compound, or to specify a particular label for
                the entry, or else ... Used for further analysis and plotting purposes.
                An attribute can be anything but must be MSONable.
        """
        comp = Composition(composition)
        PDEntry.__init__(self, comp, energy, attribute=attribute)
        self.correction = correction
        self.parameters = parameters if parameters else {}
        self.data = data if data else {}
        self.entry_id = entry_id
        self._attribute = attribute

    @property
    def energy(self):
        """
        Returns the *corrected* energy of the entry.
        """
        return self.uncorrected_energy + self.correction

    @property
    def uncorrected_energy(self):
        """
        Returns the *uncorrected* energy of the entry.
        """
        return super(ComputedEntry, self).energy

    def __repr__(self):
        output = ["ComputedEntry {}".format(self.composition.formula),
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

    @classmethod
    def from_dict(cls, d):
        dec = PMGJSONDecoder()
        return cls(d["composition"], d["energy"], d["correction"],
                   dec.process_decoded(d.get("parameters", {})),
                   dec.process_decoded(d.get("data", {})),
                   entry_id=d.get("entry_id", None),
                   attribute=d["attribute"] if "attribute" in d else None)

    @property
    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "energy": self.uncorrected_energy,
                "composition": self.composition.to_dict,
                "correction": self.correction,
                "parameters": json.loads(json.dumps(self.parameters,
                                                    cls=PMGJSONEncoder)),
                "data": json.loads(json.dumps(self.data, cls=PMGJSONEncoder)),
                "entry_id": self.entry_id,
                "attribute": self._attribute}


class ComputedStructureEntry(ComputedEntry):
    """
    A heavier version of ComputedEntry which contains a structure as well. The
    structure is needed for some analyses.
    """

    def __init__(self, structure, energy, correction=0.0, parameters=None,
                 data=None, entry_id=None):
        """
        Args:
            structure:
                The actual structure of an entry.
            energy:
                Energy of the entry. Usually the final calculated energy from
                VASP or other electronic structure codes.
            correction:
                A correction to be applied to the energy. This is used to
                modify the energy for certain analyses. Defaults to 0.0.
            parameters:
                An optional dict of parameters associated with the entry.
                Defaults to None.
            data:
                An optional dict of any additional data associated with the
                entry. Defaults to None.
            entry_id:
                An optional id to uniquely identify the entry.
        """
        ComputedEntry.__init__(self, structure.composition, energy,
                               correction=correction, parameters=parameters,
                               data=data, entry_id=entry_id)
        self.structure = structure

    def __repr__(self):
        output = ["ComputedStructureEntry {}".format(self.composition.formula),
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

    @property
    def to_dict(self):
        d = super(ComputedStructureEntry, self).to_dict
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["structure"] = self.structure.to_dict
        return d

    @classmethod
    def from_dict(cls, d):
        dec = PMGJSONDecoder()
        return cls(dec.process_decoded(d["structure"]),
                   d["energy"], d["correction"],
                   dec.process_decoded(d.get("parameters", {})),
                   dec.process_decoded(d.get("data", {})),
                   entry_id=d.get("entry_id", None))
