from pymatgen.core.structure import Structure
from monty.json import MSONable
import numpy as np


class Trajectory(MSONable):
    def __init__(self, structure, displacements, site_properties):
        self.structure = structure
        self.displacements = displacements
        self.sequential_displacements = None
        self.site_properties = site_properties
        self.__structure = structure

    def __getattr__(self, attr):

        if hasaatr(self.__structure, attr):
            return self.__structure.__getattr__(attr)
        else:
            raise Exception("Neither Trajectory nor structure has attribute: {}".format(attr))

    def change_index(self, index):
        if index > len(self.xyz_trajectory):
            raise Exception
        else:
            coords = self.structure.frac_coords + self.displacements[index]
            self.__structure = Structure(self.structure.lattice, self.structure.species,
                                         coords, site_properties=self.site_properties[index])

    def get_diplacements(self, sequential=True):
        if sequential:
            if not self.sequential_displacements:
                seq_displacements = np.subtract(self.displacements, self.structure.frac_coords)
                self.sequential_displacements = seq_displacements
            return self.sequential_displacements
        else:
            return self.displacements

    @classmethod
    def from_structures(cls, structures):
        """
        Convenience constructor to make a Trajectory from a list of Structures
        """
        displacements = [[] * len(structures)]
        site_properties = [[] * len(structures)]
        for i, structure in enumerate(structures):
            displacements[i] = structure.frac_coords - structures[0].frac_coords
            site_properties[i] = structure.site_properties
        return cls(structures[0], displacements, site_properties)

    def as_dict(self):
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__}
        d["structure"] = self.structure.as_dict()
        d["displacements"] = self.displacements
        d["site_properties"] = self.site_properties
        return d

    @classmethod
    def from_dict(cls, d):
        structure = Structure.from_dict(d["structure"])
        return cls(structure, d["displacements"], d["site_properties"])

