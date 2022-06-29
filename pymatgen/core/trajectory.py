# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes to define a simulation trajectory, which could come from
either relaxation or molecular dynamics.
"""

from __future__ import annotations

import itertools
import os
import warnings
from fnmatch import fnmatch
from typing import Any, Optional, Sequence

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

__author__ = "Eric Sivonxay, Shyam Dwaraknath, Mingjian Wen"
__version__ = "0.1"
__date__ = "Jun 29, 2022"

Vector3D = tuple[float, float, float]
Matrix3D = tuple[Vector3D, Vector3D, Vector3D]


class Trajectory(MSONable):
    """
    Trajectory object that stores structural information related to an MD simulation.

    Provides basic functions such as slicing trajectory or obtaining displacements.
    """

    def __init__(
        self,
        lattice: Lattice | Matrix3D | list[Lattice] | list[Matrix3D] | np.ndarray,
        species: list[str | Element | Species | DummySpecies | Composition],
        frac_coords: list[list[Vector3D]] | np.ndarray,
        *,
        constant_lattice: bool = True,
        site_properties: Optional[list[dict[str, Sequence[Any]]]] = None,
        frame_properties: Optional[list[dict[str, Any]]] = None,
        time_step: Optional[int | float] = None,
        coords_are_displacement: bool = False,
        base_positions: list[list[Vector3D]] = None,
        **kwargs,
    ):
        """
        Create a trajectory of N species with M frames (steps).

        In below, `N` denotes the number of sites in the structure, and `M` denotes the
        number of frames in the trajectory.

        Args:
            lattice: shape (3, 3) or (M, 3, 3). Lattice of the structures in the
                trajectory; should be used together with `constant_lattice`.
                If `constant_lattice=True`, this should be a single lattice that is
                common for all structures in the trajectory (e.g. in an NVT run).
                If `constant_lattice=False`, this should be a list of lattices,
                each for one structure in the trajectory (e.g. in an NPT run or a
                relaxation that allows changing the cell size).
            species: shape (N,). List of species on each site. Can take in flexible
                input, including:
                i.  A sequence of element / species specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g., (3, 56, ...) or actual Element or Species objects.
                ii. List of dict of elements/species and occupancies, e.g.,
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            frac_coords: shape (M, N, 3). fractional coordinates of the sites.
            constant_lattice: Whether the lattice changes during the simulation.
                Should be used together with `lattice`. See usage there.
            time_step: Timestep of MD simulation in femto-seconds. Should be `None`
                for relaxation trajectory.
            site_properties: Properties associated with the sites. This should be a
                sequence of `M` dicts, with each dict providing the site properties for
                a frame. Each value in a dict should be a sequence of length `N`, giving
                the properties of the `N` sites. For example, for a trajectory with
                `M=2` and `N=4`, the `site_properties` can be:
                [{"magmom":[5,5,5,5]}, {"magmom":[5,5,5,5]}].
            frame_properties: Properties associated with the structure (e.g. total
                energy). This should be a sequence of `M` dicts, with each dict
                providing the properties for a frame. For example, for a trajectory with
                `M=2`, the `frame_properties` can be [{'energy':1.0}, {'energy':2.0}].
            coords_are_displacement: Whether `frac_coords` are given in displacements
                (True) or positions (False). Note, if this is `True`, `frac_coords`
                of a frame (say i) should be relative to the previous frame (i.e.
                i-1), but not relative to the `base_position`.
            base_positions: shape (N, 3). The starting positions of all atoms in the
                trajectory. Used to reconstruct positions when converting from
                displacements to positions. Only needs to be specified if
                `coords_are_displacement=True`. Defaults to the first index of
                `frac_coords` when `coords_are_displacement=False`.
            kwargs: additional properties of the trajectory.
        """

        if isinstance(lattice, Lattice):
            lattice = lattice.matrix
        elif isinstance(lattice, list) and isinstance(lattice[0], Lattice):
            lattice = [x.matrix for x in lattice]
        lattice = np.asarray(lattice)

        if not constant_lattice and lattice.shape == (3, 3):
            self.lattice = np.tile(lattice, (len(frac_coords), 1, 1))
            warnings.warn(
                "Get `constant_lattice=False`, but only get a single `lattice`. "
                "Use this single `lattice` as the lattice for all frames."
            )
        else:
            self.lattice = lattice

        self.constant_lattice = constant_lattice

        if coords_are_displacement:
            if base_positions is None:
                warnings.warn(
                    "Without providing an array of starting positions, the positions "
                    "for each time step will not be available."
                )
            self.base_positions = base_positions
        else:
            self.base_positions = frac_coords[0]
        self.coords_are_displacement = coords_are_displacement

        self.species = species
        self.frac_coords = np.asarray(frac_coords)
        self.site_properties = site_properties
        self.frame_properties = frame_properties
        self.time_step = time_step
        self.kwargs = kwargs

    def get_structure(self, i: int) -> Structure:
        """
        Get structure at specified index.

        Args:
            i: Index of structure.

        Returns:
            A pymatgen Structure object.
        """
        return self[i]

    def to_positions(self):
        """
        Convert displacements between consecutive frames into positions.

        `base_positions` and `frac_coords` should both be in fractional coords or
        absolute coords.

        This is the opposite operation of `to_displacements()`.
        """
        if self.coords_are_displacement:
            cumulative_displacements = np.cumsum(self.frac_coords, axis=0)
            positions = self.base_positions + cumulative_displacements
            self.frac_coords = positions
            self.coords_are_displacement = False

    def to_displacements(self):
        """
        Converts positions of trajectory into displacements between consecutive frames.

        `base_positions` and `frac_coords` should both be in fractional coords. Does
        not work for absolute coords because the atoms are to be wrapped into the
        simulation box.

        This is the opposite operation of `to_positions()`.
        """
        if not self.coords_are_displacement:

            displacements = np.subtract(
                self.frac_coords,
                np.roll(self.frac_coords, 1, axis=0),
            )
            displacements[0] = np.zeros(np.shape(self.frac_coords[0]))

            # Deal with PBC.
            # For example - If in one frame an atom has fractional coordinates of
            # [0, 0, 0.98] and in the next its coordinates are [0, 0, 0.01], this atom
            # will have moved 0.03*c, but if we only subtract the positions, we would
            # get a displacement vector of [0, 0, -0.97]. Therefore, we can correct for
            # this by adding or subtracting 1 from the value
            displacements = [np.subtract(d, np.around(d)) for d in displacements]

            self.frac_coords = displacements
            self.coords_are_displacement = True

    def extend(self, trajectory: Trajectory):
        """
        Append a trajectory to the current one.

        The lattice, coords, and all other properties are combined.

        Args:
            trajectory: Trajectory to append.
        """
        if self.time_step != trajectory.time_step:
            raise ValueError(
                "Cannot extend trajectory. Time steps of the trajectories are "
                f"incompatible: {self.time_step} and {trajectory.time_step}."
            )

        if self.species != trajectory.species:
            raise ValueError(
                "Cannot extend trajectory. Species in the trajectories are "
                f"incompatible: {self.species} and {trajectory.species}."
            )

        # Ensure both trajectories are in positions before combining
        self.to_positions()
        trajectory.to_positions()

        self.frac_coords = np.concatenate((self.frac_coords, trajectory.frac_coords))

        self.site_properties = self._combine_site_props(
            self.site_properties,
            trajectory.site_properties,
            len(self),
            len(trajectory),
        )
        self.frame_properties = self._combine_frame_props(
            self.frame_properties,
            trajectory.frame_properties,
            len(self),
            len(trajectory),
        )

        self.lattice, self.constant_lattice = self._combine_lattice(
            self.lattice,
            trajectory.lattice,
            len(self),
            len(trajectory),
        )

    def __iter__(self):
        """
        Iterator of the trajectory, yielding a pymatgen structure for each frame.
        """
        for i in range(len(self)):
            yield self[i]

    def __len__(self):
        """
        Number of frames in the trajectory.
        """
        return len(self.frac_coords)

    def __getitem__(self, frames: int | slice) -> Structure | Trajectory:
        """
        Get a subset of the trajectory.

        The output depends on the type of the input `frames`. If an int is given, return
        a pymatgen Structure at the specified frame. If a slice is given, return a new
        trajectory with the given subset of frames.

        Args:
            frames: Indices of the trajectory to return.

        Return:
            Subset of trajectory
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

    def copy(self) -> Trajectory:
        """
        Copy of Trajectory.
        """
        return Trajectory(
            self.lattice,
            self.species,
            self.frac_coords,
            constant_lattice=self.constant_lattice,
            site_properties=self.site_properties,
            frame_properties=self.frame_properties,
            time_step=self.time_step,
            coords_are_displacement=False,
            base_positions=self.base_positions,
            **self.kwargs,
        )

    @classmethod
    def from_structures(cls, structures, constant_lattice=True, **kwargs):
        """
        Convenience constructor to obtain trajectory from a list of structures.
        Note: Assumes no atoms removed during simulation

        Args:
            structures (list): list of pymatgen Structure objects.
            constant_lattice (bool): Whether the lattice changes during the simulation,
            such as in an NPT MD simulation. True results in
        Returns:
            (Trajectory)
        """
        frac_coords = [structure.frac_coords for structure in structures]
        if constant_lattice:
            lattice = structures[0].lattice.matrix
        else:
            lattice = [structure.lattice.matrix for structure in structures]

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
            constant_lattice (bool): Whether the lattice changes during the simulation,
            such as in an NPT MD simulation. True results in

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

    def as_dict(self) -> dict:
        """
        :return: MSONAble dict.
        """
        d = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "lattice": self.lattice.tolist(),
            "species": self.species,
            "frac_coords": self.frac_coords.tolist(),
            "constant_lattice": self.constant_lattice,
            "site_properties": self.site_properties,
            "frame_properties": self.frame_properties,
            "time_step": self.time_step,
            "coords_are_displacement": self.coords_are_displacement,
            "base_positions": self.base_positions,
            "kwargs": self.kwargs,
        }

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
        all_keys = set(attr_1.keys()).union(set(attr_2.keys()))

        # Initialize dict with the common keys
        new_frame_props = dict(zip(all_keys, [[] for i in all_keys]))
        for key in all_keys:
            if key in attr_1.keys():
                new_frame_props[key].extend(attr_1[key])
            else:
                # If key doesn't exist in the first trajectory, append None for each index
                new_frame_props[key].extend([None for i in range(len_1)])

            if key in attr_2.keys():
                new_frame_props[key].extend(attr_2[key])
            else:
                # If key doesn't exist in the second trajectory, append None for each index
                new_frame_props[key].extend([None for i in range(len_2)])

        return new_frame_props

    def write_Xdatcar(
        self,
        filename: str = "XDATCAR",
        system: str = None,
        significant_figures: int = 6,
    ):
        """
        Writes Xdatcar to a file. The supported kwargs are the same as those for
        the Xdatcar_from_structs.get_string method and are passed through directly.

        Args:
            filename: name of file (It's prudent to end the filename with 'XDATCAR', as
                most visualization and analysis software require this for autodetection)
            system: Description of system
            significant_figures: Significant figures in the output file
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


def _collate(x: Any, y: Any) -> Any:
    """
    Collate two data containers x and y of the same data structure.

    Support data structure include tuple, list, and dict.

    Based on the default_collate function in PyTorch at:
    https://github.com/pytorch/pytorch/blob/da61ec2a4a3e7149a0cea6538d57f73412cc433b/torch/utils/data/_utils/collate.py#L84

    Args:
        x: the first data container.
        y: the second data container.

    Returns:
        Combined data.
    """
