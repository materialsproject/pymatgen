#!/usr/bin/env python

'''
This module implements equivalents of the basic ComputedEntry objects, which
is the basic entity that can be used to perform many analyses. ComputedEntries
contain calculated information, typically from VASP or other electronic structure
codes. For example, ComputedEntries can be used as inputs for phase diagram analysis.
'''

__author__ = "Shyue Ping Ong, Anubhav Jain"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ = "Sep 23, 2011M"

import json

from pymatgen.phasediagram.entries import PDEntry
from pymatgen.core.structure import Structure, Composition

class ComputedEntry(PDEntry):
    """
    An lightweight ComputedEntry object containing key computed data
    for many purposes.  Extends a PDEntry so that it can be used for phase diagram
    generation.
    
    Author: Shyue Ping Ong, Anubhav Jain
    """

    def __init__(self, composition, energy, correction = 0.0, parameters = None, data = None):
        """
        Args:
            composition:
                Composition of the entry. For flexibility, this can take the form
                of a {symbol: amt} dictionary as well the standard pymatgen
                Composition object.
            energy:
                Energy of the entry. Usually the final calculated energy from
                VASP or other electronic structure codes.
            correction:
                A correction to be applied to the energy. This is used to modify
                the energy for certain analyses. Defaults to 0.0.
            parameters:
                An optional dict of parameters associated with the entry. Defaults
                to None.
            data:
                An optional dict of any additional data associated with 
                the entry. Defaults to None
        """
        if not isinstance(composition, Composition):
            comp = Composition.from_dict(composition)
        else:
            comp = composition
        super(ComputedEntry, self).__init__(comp, energy)
        self.correction = correction
        self.parameters = parameters if parameters else {}
        self.data = data if data else {}

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
        return "ComputedEntry {} with energy = {:.4f}, correction = {:.4f}".format(self.composition.formula, self.uncorrected_energy, self.correction)

    def __str__(self):
        outputstr = [self.__repr__()]
        outputstr.append("Parameters:")
        for k, v in self.parameters.items():
            outputstr.append("{} = {}".format(k, v))
        outputstr.append("Data:")
        for k, v in self.data.items():
            outputstr.append("{} = {}".format(k, v))
        return "\n".join(outputstr)

    @staticmethod
    def from_dict(d):
        return ComputedEntry(d['composition'], d['energy'], d['correction'],
                             d['parameters'], d['data'])

    @property
    def to_dict(self):
        d = dict()
        d['energy'] = self.uncorrected_energy
        d['composition'] = self.composition.to_dict
        d['correction'] = self.correction
        d['parameters'] = self.parameters
        d['data'] = self.data
        return d

class ComputedStructureEntry(ComputedEntry):
    """
    A heavier version of ComputedEntry which contains a structure as well. The 
    structure is needed for some analyses.
    
    Author: Shyue Ping Ong, Anubhav Jain
    """

    def __init__(self, structure, energy, correction = 0.0, parameters = None, data = None):
        """
        Args:
            structure:
                The actual structure of an entry.
            energy:
                Energy of the entry. Usually the final calculated energy from
                VASP or other electronic structure codes.
            correction:
                A correction to be applied to the energy. This is used to modify
                the energy for certain analyses. Defaults to 0.0.
            parameters:
                An optional dict of parameters associated with the entry. Defaults
                to None.
            data:
                An optional dict of any additional data associated with 
                the entry. Defaults to None
        """
        super(ComputedStructureEntry, self).__init__(structure.composition, energy, correction = correction, parameters = parameters, data = data)
        self.structure = structure

    def __repr__(self):
        return "ComputedStructureEntry {} with energy = {:.4f}, correction = {:.4f}".format(self.composition.formula, self.uncorrected_energy, self.correction)

    def __str__(self):
        outputstr = [self.__repr__()]
        outputstr.append("Parameters:")
        for k, v in self.parameters.items():
            outputstr.append("{} = {}".format(k, v))
        outputstr.append("Data:")
        for k, v in self.data.items():
            outputstr.append("{} = {}".format(k, v))
        return "\n".join(outputstr)

    @property
    def to_dict(self):
        d = super(ComputedStructureEntry, self).to_dict
        d['structure'] = self.structure.to_dict
        return d

    @staticmethod
    def from_dict(d):
        return ComputedStructureEntry(Structure.from_dict(d['structure']), d['energy'],
                             d['correction'], d['parameters'], d['data'])


def computed_entries_to_json(all_entries):
    jsonlist = list()
    for entry in all_entries:
        jsonlist.append(entry.to_dict)
    return json.dumps(jsonlist)

def computed_entries_from_json(jsonstr):
    json_list = json.loads(jsonstr)
    entries = list()
    for d in json_list:
        if 'structure' in d:
            entries.append(ComputedStructureEntry.from_dict(d))
        else:
            entries.append(ComputedEntry.from_dict(d))
    return entries
