from pymatgen.core.structure import Structure, Lattice
from monty.json import MSONable
import numpy as np
import warnings


class Trajectory(MSONable):
    def __init__(self, frac_coords, lattice, species, time_step=2, site_properties=None, constant_lattice=True,
                 coords_are_displacement=False, base_positions=None):
        # To support from_dict and as_dict
        if type(frac_coords) == list:
            frac_coords = np.array(frac_coords)

        self.frac_coords = frac_coords
        if coords_are_displacement:
            if base_positions == None:
                warnings.warn(
                    "Without providing an array of starting positions, the positions for each time step will not be available")
            self.base_positions = base_positions
        else:
            self.base_positions = frac_coords[0]
        self.coords_are_displacement = coords_are_displacement

        self.lattice = lattice
        self.constant_lattice = constant_lattice
        self.species = species
        self.site_properties = site_properties
        self.time_step = time_step

    def get_structure(self, i):
        return self[i]

    def to_positions(self):
        if self.coords_are_displacement:
            cumulative_displacements = np.cumsum(self.frac_coords, axis=0)
            positions = self.base_positions + cumulative_displacements
            self.frac_coords = positions
            self.coords_are_displacement = False
        return

    def to_displacements(self):
        if not self.coords_are_displacement:
            displacements = np.subtract(self.frac_coords, np.roll(self.frac_coords, 1, axis=0))
            displacements[0] = np.zeros(np.shape(self.frac_coords[0]))
            # Deal with PBC
            displacements = [np.subtract(item, np.round(item)) for item in displacements]

            self.frac_coords = displacements
            self.coords_are_displacement = True
        return

    def extend(self, trajectory):
        if self.time_step != trajectory.time_step:
            warnings.warn('Trajectory not extended: Time steps of trajectories is incompatible')
            return
            # raise Exception('Time steps of trajectories is incompatible')

        if self.species == trajectory.species:
            warnings.warn('Trajectory not extended: species in trajectory do not match')
            return
            # raise Exception('Species in trajectory do not match')

        if self.coords_are_displacement:
            self.to_positions()

        if trajectory.coords_are_displacement:
            trajectory.to_positions()

        self.frac_coords = np.concatenate((self.frac_coords, trajectory), axis=0)
        self.lattice = combine_attribute(self.lattice, trajectory.lattice, self.shape[0], trajectory.size[0])
        self.site_properties = combine_attribute(self.site_properties, trajectory.site_properties,
                                                self.frac_coords.shape[0], trajectory.frac_coords.shape[0])
        return

    def __iter__(self):
        for i in range(self.frac_coords.shape[0]):
            yield self[i]

    def __len__(self):
        return np.shape(self.frac_coords)[0]

    def __getitem__(self, frames):
        """
        Gets a subset of the trajectory if a slice is given, if an int is given, return a structure
        :param frames:
        :return:
        """
        if type(frames) == int and frames < self.frac_coords.shape[0]:
            lattice = self.lattice if np.shape(self.lattice) == (3, 3) else self.lattice[frames]
            site_properties = self.site_properties[frames] if self.site_properties else None
            return Structure(Lattice(lattice), self.species, self.frac_coords[frames], site_properties=site_properties,
                             to_unit_cell=True)

        if type(frames) == slice:
            frames = np.arange(frames.start, frames.stop, frames.step)
        elif type(frames) not in [list, np.ndarray]:
            try:
                frames = np.asarray(frames)
            except:
                raise Exception('Given accessor is not of type int, slice, tuple, list, or array')

        if type(frames) in [list, np.ndarray] and (np.asarray([frames]) < self.frac_coords.shape[0]).all():
            if np.shape(self.lattice) == (3, 3):
                lattice = self.lattice
            else:
                lattice = self.lattice[frames, :]
            return Trajectory(self.frac_coords[frames, :], lattice, self.species, self.time_step,
                              self.site_properties)
        else:
            warnings.warn('Some or all selected frames exceed trajectory length')
        return

    def copy(self):
        return Trajectory(self.frac_coords, self.lattice, self.species, self.time_step, self.site_properties,
                          self.constant_lattice, self.coords_are_displacement, self.base_positions)

    @classmethod
    def from_structures(cls, structures, const_lattice=False, **kwargs):
        """
        Assumes no atoms removed from simulation
        """
        frac_coords = [structure.frac_coords for structure in structures]
        if const_lattice:
            lattice = structures[0].lattice.matrix
        else:
            lattice = [structure.lattice.matrix for structure in structures]
        site_properties = [structure.site_properties for structure in structures]
        return cls(frac_coords, lattice, species=structures[0].species, site_properties=site_properties, **kwargs)


def combine_attribute(attr_1, attr_2, len_1, len_2):
    if type(attr_1) == list or type(attr_2) == list:
        attribute = np.concatenate((attr_1, attr_2), axis=0)
    else:
        if type(attr_1) != list and type(attr_2) != list and attr_1 == attr_2:
            attribute = attr_1
        else:
            attribute = [attr_1.copy()] * len_1 if type(attr_1) != list else attr_1.copy()
            attribute.extend([attr_2.copy()] * len_2 if type(attr_2 != list) else attr_2.copy())
    return attribute
