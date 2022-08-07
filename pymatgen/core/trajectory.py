# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes used to define a MD trajectory.
"""

from __future__ import annotations

import itertools
import os
import warnings
from fnmatch import fnmatch
from typing import Sequence

import numpy as np
from monty.io import zopen
from monty.json import MSONable

from pymatgen.core.structure import (
    Composition,
    DummySpecies,
    Element,
    Lattice,
    Species,
    Structure,
)
from pymatgen.io.vasp.outputs import Vasprun, Xdatcar

__author__ = "Eric Sivonxay, Shyam Dwaraknath"
__version__ = "0.0"
__date__ = "Jan 25, 2019"


class Trajectory(MSONable):
    """
    Trajectory object that stores structural information related to a MD simulation.
    Provides basic functions such as slicing trajectory or obtaining displacements.
    """

    def __init__(
        self,
        lattice: Sequence[float] | Sequence[Sequence[float]] | np.ndarray | Lattice,
        species: list[str | Element | Species | DummySpecies | Composition],
        frac_coords: list[Sequence[Sequence[float]]] | np.ndarray,
        time_step: float = 2,
        site_properties: dict = None,
        frame_properties: dict = None,
        constant_lattice: bool = True,
        coords_are_displacement: bool = False,
        base_positions: Sequence[Sequence[float]] = None,
    ):
        """
        Create a trajectory object
        Args:
            lattice: The lattice as any 2D array. Each row should correspond to a lattice
                vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
                lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
            species: List of species on each site. Can take in flexible input,
                including:
                i.  A sequence of element / species specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g., (3, 56, ...) or actual Element or Species objects.
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
            frame_properties (dict): Properties of the trajectory such as energy, pressure, etc. each property should
                have a length equal to the trajectory length. eg: {'energy': [#, #, #, #], 'pressure': [0, 0.1, 0 0.02]}
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

        if isinstance(lattice, Lattice):
            lattice = lattice.matrix

        if isinstance(lattice, list):
            lattice = np.array(lattice)

        self.frac_coords = frac_coords
        if coords_are_displacement:
            if base_positions is None:
                warnings.warn(
                    "Without providing an array of starting positions, \
                               the positions for each time step will not be available"
                )
            self.base_positions = base_positions
        else:
            self.base_positions = frac_coords[0]
        self.coords_are_displacement = coords_are_displacement

        if not constant_lattice and np.shape(lattice) == (3, 3):
            self.lattice = [lattice for i in range(np.shape(self.frac_coords)[0])]
        else:
            self.lattice = lattice  # type: ignore
        self.constant_lattice = constant_lattice
        self.species = species
        self.site_properties = site_properties
        self.frame_properties = frame_properties
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

    def extend(self, trajectory):
        """
        Concatenate another trajectory
        Args:
            trajectory (Trajectory): Trajectory to add
        """
        if self.time_step != trajectory.time_step:
            raise ValueError("Trajectory not extended: Time steps of trajectories is incompatible")

        if len(self.species) != len(trajectory.species) and self.species != trajectory.species:
            raise ValueError("Trajectory not extended: species in trajectory do not match")

        # Ensure both trajectories are in positions before combining
        self.to_positions()
        trajectory.to_positions()

        self.site_properties = self._combine_site_props(
            self.site_properties,
            trajectory.site_properties,
            np.shape(self.frac_coords)[0],
            np.shape(trajectory.frac_coords)[0],
        )
        self.frame_properties = self._combine_frame_props(
            self.frame_properties,
            trajectory.frame_properties,
            np.shape(self.frac_coords)[0],
            np.shape(trajectory.frac_coords)[0],
        )

        self.frac_coords = np.concatenate((self.frac_coords, trajectory.frac_coords), axis=0)
        self.lattice, self.constant_lattice = self._combine_lattice(
            self.lattice,
            trajectory.lattice,
            np.shape(self.frac_coords)[0],
            np.shape(trajectory.frac_coords)[0],
        )

    def __iter__(self):
        for i in range(np.shape(self.frac_coords)[0]):
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
        # If trajectory is in displacement mode, return the displacements at that frame
        if self.coords_are_displacement:
            if isinstance(frames, int):
                if frames >= np.shape(self.frac_coords)[0]:
                    raise ValueError("Selected frame exceeds trajectory length")
                # For integer input, return the displacements at that timestep
                return self.frac_coords[frames]
            if isinstance(frames, slice):
                # For slice input, return a list of the displacements
                start, stop, step = frames.indices(len(self))
                return [self.frac_coords[i] for i in range(start, stop, step)]
            if isinstance(frames, (list, np.ndarray)):
                # For list input, return a list of the displacements
                pruned_frames = [i for i in frames if i < len(self)]  # Get rid of frames that exceed trajectory length
                if len(pruned_frames) < len(frames):
                    warnings.warn("Some or all selected frames exceed trajectory length")
                return [self.frac_coords[i] for i in pruned_frames]
            raise Exception("Given accessor is not of type int, slice, list, or array")

        # If trajectory is in positions mode, return a structure for the given frame or trajectory for the given frames
        if isinstance(frames, int):
            if frames >= np.shape(self.frac_coords)[0]:
                raise ValueError("Selected frame exceeds trajectory length")
            # For integer input, return the structure at that timestep
            lattice = self.lattice if self.constant_lattice else self.lattice[frames]
            site_properties = self.site_properties[frames] if self.site_properties else None
            site_properties = self.site_properties[frames] if self.site_properties else None
            return Structure(
                Lattice(lattice),
                self.species,
                self.frac_coords[frames],
                site_properties=site_properties,
                to_unit_cell=True,
            )
        if isinstance(frames, slice):
            # For slice input, return a trajectory of the sliced time
            start, stop, step = frames.indices(len(self))
            pruned_frames = range(start, stop, step)
            lattice = self.lattice if self.constant_lattice else [self.lattice[i] for i in pruned_frames]
            frac_coords = [self.frac_coords[i] for i in pruned_frames]
            if self.site_properties is not None:
                site_properties = [self.site_properties[i] for i in pruned_frames]
            else:
                site_properties = None
            if self.frame_properties is not None:
                frame_properties = {}
                for key, item in self.frame_properties.items():
                    frame_properties[key] = [item[i] for i in pruned_frames]
            else:
                frame_properties = None
            return Trajectory(
                lattice,
                self.species,
                frac_coords,
                time_step=self.time_step,
                site_properties=site_properties,
                frame_properties=frame_properties,
                constant_lattice=self.constant_lattice,
                coords_are_displacement=False,
                base_positions=self.base_positions,
            )
        if isinstance(frames, (list, np.ndarray)):
            # For list input, return a trajectory of the specified times
            pruned_frames = [i for i in frames if i < len(self)]  # Get rid of frames that exceed trajectory length
            if len(pruned_frames) < len(frames):
                warnings.warn("Some or all selected frames exceed trajectory length")
            lattice = self.lattice if self.constant_lattice else [self.lattice[i] for i in pruned_frames]
            frac_coords = [self.frac_coords[i] for i in pruned_frames]
            if self.site_properties is not None:
                site_properties = [self.site_properties[i] for i in pruned_frames]
            else:
                site_properties = None
            if self.frame_properties is not None:
                frame_properties = {}
                for key, item in self.frame_properties.items():
                    frame_properties[key] = [item[i] for i in pruned_frames]
            else:
                frame_properties = None
            return Trajectory(
                lattice,
                self.species,
                frac_coords,
                time_step=self.time_step,
                site_properties=site_properties,
                frame_properties=frame_properties,
                constant_lattice=self.constant_lattice,
                coords_are_displacement=False,
                base_positions=self.base_positions,
            )
        raise Exception("Given accessor is not of type int, slice, tuple, list, or array")

    def copy(self):
        """
        :return: Copy of Trajectory.
        """
        return Trajectory(
            self.lattice,
            self.species,
            self.frac_coords,
            time_step=self.time_step,
            site_properties=self.site_properties,
            frame_properties=self.frame_properties,
            constant_lattice=self.constant_lattice,
            coords_are_displacement=False,
            base_positions=self.base_positions,
        )

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
        site_properties = {}

        site_properties = [structure.site_properties for structure in structures]
        return cls(
            lattice,
            structures[0].species,
            frac_coords,
            site_properties=site_properties,
            constant_lattice=constant_lattice,
            **kwargs,
        )

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
        """
        :return: MSONAble dict.
        """
        d = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "species": self.species,
            "time_step": self.time_step,
            "site_properties": self.site_properties,
            "frame_properties": self.frame_properties,
            "constant_lattice": self.constant_lattice,
            "coords_are_displacement": self.coords_are_displacement,
            "base_positions": self.base_positions,
        }
        d["lattice"] = self.lattice.tolist()
        d["frac_coords"] = self.frac_coords.tolist()

        return d

    @staticmethod
    def _combine_lattice(attr_1, attr_2, len_1, len_2):
        """
        Helper function to combine trajectory properties such as site_properties or lattice
        """
        if np.shape(attr_1) == (3, 3) and np.shape(attr_2) == (3, 3):
            attribute = attr_1
            attribute_constant = True
        elif np.shape(attr_1) == 3 and np.shape(attr_2) == 3:
            attribute = np.concatenate((attr_1, attr_2), axis=0)
            attribute_constant = False
        else:
            attribute = [attr_1.copy()] * len_1 if isinstance(attr_1, list) else attr_1.copy()
            attribute.extend([attr_2.copy()] * len_2 if isinstance(attr_2, list) else attr_2.copy())
            attribute_constant = False
        return attribute, attribute_constant

    @staticmethod
    def _combine_site_props(attr_1, attr_2, len_1, len_2):
        """
        Helper function to combine site properties of 2 trajectories
        """
        if attr_1 is None and attr_2 is None:
            return None
        if attr_1 is None or attr_2 is None:
            new_site_properties = []
            if attr_1 is None:
                new_site_properties.extend([None for i in range(len_1)])
            elif len(attr_1) == 1:
                new_site_properties.extend([attr_1[0] for i in range(len_1)])
            elif len(attr_1) > 1:
                new_site_properties.extend(attr_1)

            if attr_2 is None:
                new_site_properties.extend([None for i in range(len_2)])
            elif len(attr_2) == 1:
                new_site_properties.extend([attr_2[0] for i in range(len_2)])
            elif len(attr_2) > 1:
                new_site_properties.extend(attr_2)

            return new_site_properties
        if len(attr_1) == 1 and len(attr_2) == 1:
            # If both properties lists are do not change within their respective trajectory
            if attr_1 == attr_2:
                # If both site_properties are the same, only store one
                return attr_1
            new_site_properties = [attr_1[0] for i in range(len_1)]
            new_site_properties.extend([attr_2[0] for i in range(len_2)])
            return new_site_properties
        if len(attr_1) > 1 and len(attr_2) > 1:
            # Both properties have site properties that change within the trajectory, concat both together
            return [*attr_1, *attr_2]

        new_site_properties = []
        if attr_1 is None:
            new_site_properties.extend([None for i in range(len_1)])
        elif len(attr_1) == 1:
            new_site_properties.extend([attr_1[0] for i in range(len_1)])
        elif len(attr_1) > 1:
            new_site_properties.extend(attr_1)

        if attr_2 is None:
            new_site_properties.extend([None for i in range(len_2)])
        elif len(attr_2) == 1:
            new_site_properties.extend([attr_2[0] for i in range(len_2)])
        elif len(attr_2) > 1:
            new_site_properties.extend(attr_2)

        return new_site_properties

    @staticmethod
    def _combine_frame_props(attr_1, attr_2, len_1, len_2):
        """
        Helper function to combine frame properties such as energy or pressure
        """
        if attr_1 is None and attr_2 is None:
            return None

        # Find all common keys
        all_keys = set(attr_1).union(set(attr_2))

        # Initialize dict with the common keys
        new_frame_props = dict(zip(all_keys, [[] for i in all_keys]))
        for key in all_keys:
            if key in attr_1:
                new_frame_props[key].extend(attr_1[key])
            else:
                # If key doesn't exist in the first trajectory, append None for each index
                new_frame_props[key].extend([None for i in range(len_1)])

            if key in attr_2:
                new_frame_props[key].extend(attr_2[key])
            else:
                # If key doesn't exist in the second trajectory, append None for each index
                new_frame_props[key].extend([None for i in range(len_2)])

        return new_frame_props

    def write_Xdatcar(self, filename="XDATCAR", system=None, significant_figures=6):
        """
        Writes Xdatcar to a file. The supported kwargs are the same as those for
        the Xdatcar_from_structs.get_string method and are passed through directly.

        Args:
            filename (str): name of file (It's prudent to end the filename with 'XDATCAR',
                as most visualization and analysis software require this for autodetection)
            system (str): Description of system
            significant_figures (int): Significant figures in the output file
        """

        # Ensure trajectory is in position form
        self.to_positions()

        if system is None:
            system = f"{self[0].composition.reduced_formula}"

        lines = []
        format_str = f"{{:.{significant_figures}f}}"
        syms = [site.specie.symbol for site in self[0]]
        site_symbols = [a[0] for a in itertools.groupby(syms)]
        syms = [site.specie.symbol for site in self[0]]
        natoms = [len(tuple(a[1])) for a in itertools.groupby(syms)]

        for si, frac_coords in enumerate(self.frac_coords):
            # Only print out the info block if
            if si == 0 or not self.constant_lattice:
                lines.extend([system, "1.0"])

                if self.constant_lattice:
                    _lattice = self.lattice
                else:
                    _lattice = self.lattice[si]

                for latt_vec in _lattice:
                    lines.append(f'{" ".join([str(el) for el in latt_vec])}')

                lines.append(" ".join(site_symbols))
                lines.append(" ".join([str(x) for x in natoms]))

            lines.append(f"Direct configuration=     {si + 1}")

            for (frac_coord, specie) in zip(frac_coords, self.species):
                coords = frac_coord
                line = f'{" ".join([format_str.format(c) for c in coords])} {specie}'
                lines.append(line)

        xdatcar_string = "\n".join(lines) + "\n"

        with zopen(filename, "wt") as f:
            f.write(xdatcar_string)
