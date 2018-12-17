from pymatgen.core.structure import Structure
from monty.json import MSONable
import numpy as np


class Trajectory(MSONable):
    """
    This class stores MD trajectories or any series of structures (with the same lattice, number of atoms, etc.).
    The base structure is stored along with a list of displacements and other site properties. Structures for different
    timesteps are pulled by fast-forwarding the 'current structure' to the desired time step using the displacements.
    Stored displacements correspond to the positional differences between the 0th and nth frame.

    This is preferred to a list of structures, which is memory inefficient.

    This class is a wrapper class for pymatgen.core.structure.Structure.
    """
    def __init__(self, base_structure, displacements, site_properties):
        self.base_structure = base_structure
        self.displacements = np.array(displacements)
        self.site_properties = site_properties
        self.index = 0
        self.structure = base_structure
        self.length = len(displacements)

    def __getattr__(self, attr):
        if hasattr(self.structure, attr):
            return getattr(self.structure, attr)
        else:
            raise Exception("Neither Trajectory nor structure has attribute: {}".format(attr))

    def change_index(self, index):
        """
        Morphs the structure to the specified timestep

        Args:
            index: Index of the desired timestep
        """
        if index > len(self.displacements):
            raise Exception
        else:
            coords = self.base_structure.frac_coords + self.displacements[index]
            self.structure = Structure(self.base_structure.lattice, self.base_structure.species,
                                       coords, site_properties=self.site_properties[index])
            self.index = index

    def change_next(self):
        """
        Increment the structure to the next timestep.
        """
        if self.index + 1 < len(self.displacements):
            self.change_index(self.index+1)
        else:
            raise Exception("At the end of the trajectory, no more steps to increment")

    def change_previous(self):
        """
        Decrement the structure to the previous timestep.
        """
        if self.index > 0:
            self.change_index(self.index+1)
        else:
            raise Exception("At the start of the trajectory, no more steps to decrement")

    def combine(self, trajectory):
        """
        Extends self with the given trajectory.

        Trajectory must have same lattice and contain the same number and type of species.
        """
        if trajectory.base_structure.lattice != self.base_structure.lattice:
            raise Exception("Lattices are incompatible")
        if trajectory.base_structure.species != self.base_structure.species:
            raise Exception("Elements are not consistent between both trajectories")

        _coords = np.add(trajectory.displacements, trajectory.base_structure.frac_coords)
        _displacements = np.subtract(_coords, self.base_structure.frac_coords)
        self.displacements = np.concatenate((self.displacements, _displacements), axis=0)
        self.site_properties.extend(trajectory.site_properties)
        self.length = len(_displacements)

    def as_structures(self):
        """
        For compatibility with existing codes which use [Structure] inputs..

        :return: A list of structures
        """
        structures = [0]*len(self.displacements)
        for i in range(len(self.displacements)):
            self.change_index(i)
            structures[i] = self.structure.copy()
        return structures

    def sequential_displacements(self, skip=1, wrapped=False, cartesian=True):
        """
        Returns the frame(n) to frame(n+1) displacements. Useful for summing to obtain MSD's

        :param skip: (int) Number of time steps to skip between each returned displacement
        :param wrapped: (bool) Specifies whether the displacements returned will be wrapped or unwrapped
        :param cartesian: (bool) Specifies whether returned coordinates will be cartesian(True) or fractional(False)
        :return:
        """
        seq_displacements = np.subtract(self.displacements[::skip],
                                        np.roll(self.displacements[::skip], 1, axis=0))
        seq_displacements[0] = np.zeros(np.shape(seq_displacements[0]))

        if wrapped:
            seq_displacements = [np.subtract(item, np.round(item)) for item in seq_displacements]

        if cartesian:
            seq_displacements = np.multiply(seq_displacements, self.structure.lattice.abc)

        return seq_displacements

    def positions(self, skip=1):
        """
        :param skip: (int) Number of time steps to skip between each returned displacement
        :return:
        """
        positions = np.add(self.structure.frac_coords[::skip], self.displacements[::skip])
        return positions

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
        Convenience constructor to make a Trajectory from an the ionic steps dictionary from the parsed vasprun.xml
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
        """
        Dict representation of Trajectory.

        Returns:
            JSON serializable dict representation.
        """
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__}
        d["structure"] = self.base_structure.as_dict()
        d["displacements"] = self.displacements.tolist()
        d["site_properties"] = self.site_properties
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Reconstitute a Trajectory object from a dict representation of Structure
        created using as_dict().

        Args:
            d (dict): Dict representation of structure.

        Returns:
            Structure object
        """
        structure = Structure.from_dict(d["structure"])
        return cls(structure, np.array(d["displacements"]), d["site_properties"])

