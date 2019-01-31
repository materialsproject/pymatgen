# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from pymatgen.core.structure import Structure, Lattice
from pymatgen.io.vasp.outputs import Xdatcar, Vasprun
from monty.json import MSONable
from fnmatch import fnmatch
import numpy as np
import warnings
import os

"""
This module provides classes used to define a MD trajectory.
"""

__author__ = "Eric Sivonxay, Shyam Dwaraknath"
__version__ = "0.0"
__date__ = "Jan 25, 2019"


class Trajectory(MSONable):
    """
    Trajectory object that stores structural information related to a MD simulation.
    Provides basic functions such as slicing trajectory or obtaining displacements.
    """
    def __init__(self, lattice, species, frac_coords, time_step=2, site_properties=None, constant_lattice=True,
                 coords_are_displacement=False, base_positions=None):
        """
        Create a trajectory object

        Args:
            lattice: The lattice as any 2D array. Each row should correspond to a lattice
                vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
                lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
            species: List of species on each site. Can take in flexible input,
                including:

                i.  A sequence of element / specie specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g., (3, 56, ...) or actual Element or Specie objects.

                ii. List of dict of elements/species and occupancies, e.g.,
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            frac_coords (MxNx3 array): list of fractional coordinates of
                each species
            time_step (int, float): Timestep of simulation in femtoseconds. Defaults to 2fs.
            site_properties (list): Properties associated with the sites as a list of
                dicts of sequences, e.g., [{"magmom":[5,5,5,5]}, {"magmom":[5,5,5,5]}]. The sequences
                have to be the same length as the atomic species and fractional_coords. Number of supplied
                dicts should match number of frames in trajectory
                Defaults to None for no properties.
            constant_lattice (bool): Whether the lattice changes during the simulation, such as in an NPT MD simulation.
            coords_are_displacement (bool): Whether supplied coordinates are given in displacements (True) or
                positions (False)
            base_positions (Nx3 array): The starting positions of all atoms in trajectory. Used to reconstruct positions
                when converting from displacements to positions. Only needs to be specified if
                coords_are_displacement=True. Defaults to first index of frac_coords if coords_are_displacement=False.
        """
        # To support from_dict and as_dict
        if isinstance(frac_coords, list):
            frac_coords = np.array(frac_coords)

        if isinstance(lattice, list):
            lattice = np.array(lattice)

        self.frac_coords = frac_coords
        if coords_are_displacement:
            if base_positions is None:
                warnings.warn(
                    "Without providing an array of starting positions, the positions for each time step will not be available")
            self.base_positions = base_positions
        else:
            self.base_positions = frac_coords[0]
        self.coords_are_displacement = coords_are_displacement

        if not constant_lattice and np.shape(lattice) == (3, 3):
            self.lattice = [lattice for i in range(self.frac_coords.shape[0])]
        else:
            self.lattice = lattice
        self.constant_lattice = constant_lattice
        self.species = species
        self.site_properties = site_properties
        self.time_step = time_step

    def get_structure(self, i):
        """
        Returns structure at specified index
        Args:
            i (int): Index of structure

        Returns:
            (Structure) pymatgen structure object
        """
        return self[i]

    def to_positions(self):
        """
        Converts fractional coordinates of trajectory into positions
        """
        if self.coords_are_displacement:
            cumulative_displacements = np.cumsum(self.frac_coords, axis=0)
            positions = self.base_positions + cumulative_displacements
            self.frac_coords = positions
            self.coords_are_displacement = False
        return

    def to_displacements(self):
        """
        Converts position coordinates of trajectory into displacements between consecutive frames
        """
        if not self.coords_are_displacement:
            displacements = np.subtract(self.frac_coords, np.roll(self.frac_coords, 1, axis=0))
            displacements[0] = np.zeros(np.shape(self.frac_coords[0]))
            # Deal with PBC
            displacements = [np.subtract(item, np.round(item)) for item in displacements]

            self.frac_coords = displacements
            self.coords_are_displacement = True
        return

    def extend(self, trajectory):
        """
        Concatenate another trajectory

        Args:
            trajectory (Trajectory): Trajectory to add

        """
        if self.time_step != trajectory.time_step:
            raise ValueError('Trajectory not extended: Time steps of trajectories is incompatible')

        if len(self.species) != len(trajectory.species) and self.species != trajectory.species:
            raise ValueError('Trajectory not extended: species in trajectory do not match')

        self.to_positions()
        trajectory.to_positions()

        self.frac_coords = np.concatenate((self.frac_coords, trajectory.frac_coords), axis=0)
        self.lattice, self.constant_lattice = self._combine_attribute(self.lattice, trajectory.lattice,
                                                                      self.frac_coords.shape[0],
                                                                      trajectory.frac_coords.shape[0])
        self.site_properties = self._combine_attribute(self.site_properties, trajectory.site_properties,
                                                       self.frac_coords.shape[0], trajectory.frac_coords.shape[0])

    def __iter__(self):
        for i in range(self.frac_coords.shape[0]):
            yield self[i]

    def __len__(self):
        return np.shape(self.frac_coords)[0]

    def __getitem__(self, frames):
        """
        Gets a subset of the trajectory if a slice is given, if an int is given, return a structure

        Args:
            frames (int, slice): int or slice of trajectory to return

        Return:
            (Trajectory, Structure) Subset of trajectory
        """
        if isinstance(frames, int) and frames < self.frac_coords.shape[0]:
            lattice = self.lattice if self.constant_lattice else self.lattice[frames]
            site_properties = self.site_properties[frames] if self.site_properties else None
            return Structure(Lattice(lattice), self.species, self.frac_coords[frames], site_properties=site_properties,
                             to_unit_cell=True)

        if isinstance(frames, slice):
            frames = np.arange(frames.start, frames.stop, frames.step)
        elif not (isinstance(frames, list) or isinstance(frames, np.ndarray)):
            try:
                frames = np.asarray(frames)
            except:
                raise Exception('Given accessor is not of type int, slice, tuple, list, or array')

        if (isinstance(frames, list) or isinstance(frames, np.ndarray)) and \
                (np.asarray([frames]) < self.frac_coords.shape[0]).all():
            if self.constant_lattice:
                lattice = self.lattice
            else:
                lattice = self.lattice[frames, :]
            return Trajectory(lattice, self.species, self.frac_coords[frames, :], self.time_step,
                              self.site_properties)
        else:
            warnings.warn('Some or all selected frames exceed trajectory length')
        return

    def copy(self):
        return Trajectory(self.lattice, self.species, self.frac_coords, self.time_step, self.site_properties,
                          self.constant_lattice, self.coords_are_displacement, self.base_positions)

    @classmethod
    def from_structures(cls, structures, constant_lattice=True, **kwargs):
        """
        Convenience constructor to obtain trajectory from a list of structures.
        Note: Assumes no atoms removed during simulation

        Args:
            structures (list): list of pymatgen Structure objects.
            constant_lattice (bool): Whether the lattice changes during the simulation, such as in an NPT MD
                simulation. True results in

        Returns:
            (Trajectory)

        """
        frac_coords = [structure.frac_coords for structure in structures]
        if constant_lattice:
            lattice = structures[0].lattice.matrix
        else:
            lattice = [structure.lattice.matrix for structure in structures]
        site_properties = [structure.site_properties for structure in structures]
        return cls(lattice, structures[0].species, frac_coords, site_properties=site_properties,
                   constant_lattice=constant_lattice, **kwargs)

    @classmethod
    def from_file(cls, filename, constant_lattice=True, **kwargs):
        """
        Convenience constructor to obtain trajectory from XDATCAR or vasprun.xml file

        Args:
            filename (str): The filename to read from.
            constant_lattice (bool): Whether the lattice changes during the simulation, such as in an NPT MD
                simulation. True results in

        Returns:
            (Trajectory)
        """
        # TODO: Support other filetypes

        fname = os.path.basename(filename)
        if fnmatch(fname, "*XDATCAR*"):
            structures = Xdatcar(filename).structures
        elif fnmatch(fname, "vasprun*.xml*"):
            structures = Vasprun(filename).structures
        else:
            raise ValueError("Unsupported file")

        return cls.from_structures(structures, constant_lattice=constant_lattice, **kwargs)

    def as_dict(self):
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "species": self.species, "time_step": self.time_step,
             "site_properties": self.site_properties,
             "constant_lattice": self.constant_lattice,
             "coords_are_displacement": self.coords_are_displacement,
             "base_positions": self.base_positions}
        d["lattice"] = self.lattice.tolist()
        d["frac_coords"] = self.frac_coords.tolist()

        return d

    @staticmethod
    def _combine_attribute(attr_1, attr_2, len_1, len_2):
        """
        Helper function to combine trajectory properties such as site_properties or lattice
        """
        if isinstance(attr_1, list) or isinstance(attr_2, list):
            attribute = np.concatenate((attr_1, attr_2), axis=0)
            attribute_changes = True
        else:
            if isinstance(attr_1, list) and isinstance(attr_2, list) and np.allclose(attr_1, attr_2):
                attribute = attr_1
                attribute_changes = False
            else:
                attribute = [attr_1.copy()] * len_1 if type(attr_1) != list else attr_1.copy()
                attribute.extend([attr_2.copy()] * len_2 if type(attr_2 != list) else attr_2.copy())
                attribute_changes = True
        return attribute, attribute_changes
