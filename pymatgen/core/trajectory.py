from pymatgen.core.structure import Structure
from monty.json import MSONable
from monty.functools import lru_cache
import numpy as np


class Trajectory(MSONable):
    def __init__(self, base_structure, displacements, site_properties):
        self.base_structure = base_structure
        self.displacements = np.array(displacements)
        self.site_properties = site_properties
        self.index = 0
        self.structure = base_structure

    def __getattr__(self, attr):

        if hasattr(self.structure, attr):
            return self.structure.__getattr__(attr)
        else:
            raise Exception("Neither Trajectory nor structure has attribute: {}".format(attr))

    def change_index(self, index):
        if index > len(self.displacements):
            raise Exception
        else:
            coords = self.base_structure.frac_coords + self.displacements[index]
            self.structure = Structure(self.base_structure.lattice, self.base_structure.species,
                                       coords, site_properties=self.site_properties[index])
            self.index = index

    def change_next(self):
        if self.index + 1 < len(self.displacements):
            self.change_index(self.index+1)
        else:
            raise Exception

    def change_previous(self):
        if self.index > 0:
            self.change_index(self.index+1)
        else:
            raise Exception

    def combine(self, trajectory):
        if trajectory.base_structure.lattice != self.base_structure.lattice:
            raise Exception("Lattices are incompatible")
        if trajectory.base_structure.species != self.base_structure.species:
            raise Exception("Elements are not consistent between both trajectories")

        _coords = np.add(trajectory.displacements, trajectory.base_structure.frac_coords)
        _displacements = np.subtract(_coords, self.base_structure.frac_coords)
        self.displacements = np.concatenate((self.displacements, _displacements), axis=0)
        self.site_properties.extend(trajectory.site_properties)

    @lru_cache()
    def as_structures(self):
        structures = [0]*len(self.displacements)
        for i in range(len(self.displacements)):
            self.change_index(i)
            structures[i] = self.structure.copy()
        return structures

    @property
    @lru_cache()
    def sequential_displacements(self, skip=1):
        seq_displacements = np.subtract(self.displacements[::skip],
                                        np.roll(self.displacements[::skip], 1, axis=0))
        return seq_displacements

    @classmethod
    def from_structures(cls, structures):
        """
        Convenience constructor to make a Trajectory from a list of Structures
        """
        displacements = [0]*len(structures)
        site_properties = [0]*len(structures)
        for i, structure in enumerate(structures):
            displacements[i] = structure.frac_coords - structures[0].frac_coords
            site_properties[i] = structure.site_properties
        return cls(structures[0], displacements, site_properties)

    @classmethod
    def from_ionic_steps(cls, ionic_steps_dict):
        """
        Convenience constructor to make a Trajectory from a list of Structures
        :param ionic_steps_dict:
        :return:
        """

        structure = Structure.from_dict(ionic_steps_dict[0]["structure"])
        positions = [0] * len(ionic_steps_dict)
        for i, step in enumerate(ionic_steps_dict):
            _step = [atom['abc'] for atom in step["structure"]["sites"]]
            positions[i] = _step

        displacements = [[]] * len(ionic_steps_dict)
        site_properties = [{}] * len(ionic_steps_dict)
        for i, step in enumerate(positions):
            displacements[i] = np.subtract(positions[i], positions[0])
            # TODO: Add site properties from ionic_steps

        return cls(structure, displacements, site_properties)

    def as_dict(self):
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__}
        d["structure"] = self.structure.as_dict()
        d["displacements"] = self.displacements.tolist()
        d["site_properties"] = self.site_properties
        return d

    @classmethod
    def from_dict(cls, d):
        structure = Structure.from_dict(d["structure"])
        return cls(structure, np.array(d["displacements"]), d["site_properties"])

