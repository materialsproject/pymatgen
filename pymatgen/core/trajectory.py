# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes to define a simulation trajectory, which could come from
either relaxation or molecular dynamics.
"""

from __future__ import annotations

import itertools
import warnings
from fnmatch import fnmatch
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple, Union

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

Vector3D = Tuple[float, float, float]
Matrix3D = Tuple[Vector3D, Vector3D, Vector3D]
SitePropsType = Union[List[Dict[Any, Sequence[Any]]], Dict[Any, Sequence[Any]]]


class Trajectory(MSONable):
    """
    Trajectory of a relaxation or molecular dynamics simulation.

    Provides basic functions such as slicing trajectory, combining trajectories, and
    obtaining displacements.
    """

    def __init__(
        self,
        lattice: Lattice | Matrix3D | list[Lattice] | list[Matrix3D] | np.ndarray,
        species: list[str | Element | Species | DummySpecies | Composition],
        frac_coords: list[list[Vector3D]] | np.ndarray,
        *,
        site_properties: SitePropsType | None = None,
        frame_properties: list[dict] | None = None,
        constant_lattice: bool = True,
        time_step: int | float | None = None,
        coords_are_displacement: bool = False,
        base_positions: list[list[Vector3D]] | np.ndarray | None = None,
    ):
        """
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
            site_properties: Properties associated with the sites. This should be a
                list of `M` dicts for a single dict. If a list of dicts, each provides
                the site properties for a frame. Each value in a dict should be a
                sequence of length `N`, giving the properties of the `N` sites.
                For example, for a trajectory with `M=2` and `N=4`, the
                `site_properties` can be: [{"magmom":[5,5,5,5]}, {"magmom":[5,5,5,5]}].
                If a single dict, the site properties in the dict apply to all frames
                in the trajectory. For example, for a trajectory with `M=2` and `N=4`,
                {"magmom":[2,2,2,2]} means that, through the entire trajectory,
                the magmom are kept constant at 2 for all four atoms.
            frame_properties: Properties associated with the structure (e.g. total
                energy). This should be a sequence of `M` dicts, with each dict
                providing the properties for a frame. For example, for a trajectory with
                `M=2`, the `frame_properties` can be [{'energy':1.0}, {'energy':2.0}].
            constant_lattice: Whether the lattice changes during the simulation.
                Should be used together with `lattice`. See usage there.
            time_step: Time step of MD simulation in femto-seconds. Should be `None`
                for relaxation trajectory.
            coords_are_displacement: Whether `frac_coords` are given in displacements
                (True) or positions (False). Note, if this is `True`, `frac_coords`
                of a frame (say i) should be relative to the previous frame (i.e.
                i-1), but not relative to the `base_position`.
            base_positions: shape (N, 3). The starting positions of all atoms in the
                trajectory. Used to reconstruct positions when converting from
                displacements to positions. Only needs to be specified if
                `coords_are_displacement=True`. Defaults to the first index of
                `frac_coords` when `coords_are_displacement=False`.
        """
        if isinstance(lattice, Lattice):
            lattice = lattice.matrix
        elif isinstance(lattice, list) and isinstance(lattice[0], Lattice):
            lattice = [x.matrix for x in lattice]  # type: ignore
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
            self.base_positions = frac_coords[0]  # type: ignore[assignment]
        self.coords_are_displacement = coords_are_displacement

        self.species = species
        self.frac_coords = np.asarray(frac_coords)
        self.time_step = time_step

        self._check_site_props(site_properties)
        self.site_properties = site_properties

        self._check_frame_props(frame_properties)
        self.frame_properties = frame_properties

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
            # this by adding or subtracting 1 from the value.
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

        # Note, this should be after the other self._combine... method calls, since
        # len(self) is used there.
        self.frac_coords = np.concatenate((self.frac_coords, trajectory.frac_coords))

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

    def __getitem__(self, frames: int | slice | list[int]) -> Structure | Trajectory:
        """
        Get a subset of the trajectory.

        The output depends on the type of the input `frames`. If an int is given, return
        a pymatgen Structure at the specified frame. If a list or a slice, return a new
        trajectory with a subset of frames.

        Args:
            frames: Indices of the trajectory to return.

        Return:
            Subset of trajectory
        """
        # Convert to position mode if not ready
        self.to_positions()

        # For integer input, return the structure at that frame
        if isinstance(frames, int):

            if frames >= len(self):
                raise IndexError(f"Frame index {frames} out of range.")

            lattice = self.lattice if self.constant_lattice else self.lattice[frames]

            return Structure(
                Lattice(lattice),
                self.species,
                self.frac_coords[frames],
                site_properties=self._get_site_props(frames),  # type: ignore
                to_unit_cell=True,
            )

        # For slice input, return a trajectory
        if isinstance(frames, (slice, list, np.ndarray)):

            if isinstance(frames, slice):
                start, stop, step = frames.indices(len(self))
                selected = list(range(start, stop, step))
            else:
                # Get rid of frames that exceed trajectory length
                selected = [i for i in frames if i < len(self)]

                if len(selected) < len(frames):
                    bad_frames = [i for i in frames if i > len(self)]
                    raise IndexError(f"Frame index {bad_frames} out of range.")

            lattice = self.lattice if self.constant_lattice else self.lattice[selected]
            frac_coords = self.frac_coords[selected]

            if self.frame_properties is not None:
                frame_properties = [self.frame_properties[i] for i in selected]
            else:
                frame_properties = None

            return Trajectory(
                lattice,
                self.species,
                frac_coords,
                site_properties=self._get_site_props(selected),
                frame_properties=frame_properties,
                constant_lattice=self.constant_lattice,
                time_step=self.time_step,
                coords_are_displacement=False,
                base_positions=self.base_positions,
            )

        supported = [int, slice, list or np.ndarray]
        raise ValueError(f"Expect the type of frames be one of {supported}; {type(frames)}.")

    def write_Xdatcar(
        self,
        filename: str | Path = "XDATCAR",
        system: str | None = None,
        significant_figures: int = 6,
    ):
        """
        Writes to Xdatcar file.

        The supported kwargs are the same as those for the
        Xdatcar_from_structs.get_string method and are passed through directly.

        Args:
            filename: Name of file to write.  It's prudent to end the filename with
                'XDATCAR', as most visualization and analysis software require this
                for autodetection.
            system: Description of system (e.g. 2D MoS2).
            significant_figures: Significant figures in the output file.
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
        n_atoms = [len(tuple(a[1])) for a in itertools.groupby(syms)]

        for si, frac_coords in enumerate(self.frac_coords):
            # Only print out the info block if
            if si == 0 or not self.constant_lattice:
                lines.extend([system, "1.0"])

                if self.constant_lattice:
                    _lattice = self.lattice
                else:
                    _lattice = self.lattice[si]

                for latt_vec in _lattice:
                    lines.append(f'{" ".join(map(str, latt_vec))}')

                lines.append(" ".join(site_symbols))
                lines.append(" ".join(map(str, n_atoms)))

            lines.append(f"Direct configuration=     {si + 1}")

            for (frac_coord, specie) in zip(frac_coords, self.species):
                coords = frac_coord
                line = f'{" ".join(format_str.format(c) for c in coords)} {specie}'
                lines.append(line)

        xdatcar_string = "\n".join(lines) + "\n"

        with zopen(filename, "wt") as f:
            f.write(xdatcar_string)

    def as_dict(self) -> dict:
        """
        Return the trajectory as a MSONAble dict.
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "lattice": self.lattice.tolist(),
            "species": self.species,
            "frac_coords": self.frac_coords.tolist(),
            "site_properties": self.site_properties,
            "frame_properties": self.frame_properties,
            "constant_lattice": self.constant_lattice,
            "time_step": self.time_step,
            "coords_are_displacement": self.coords_are_displacement,
            "base_positions": self.base_positions,
        }

    @classmethod
    def from_structures(
        cls,
        structures: list[Structure],
        constant_lattice: bool = True,
        **kwargs,
    ) -> Trajectory:
        """
        Create trajectory from a list of structures.

        Note: Assumes no atoms removed during simulation.

        Args:
            structures: pymatgen Structure objects.
            constant_lattice: Whether the lattice changes during the simulation,
                such as in an NPT MD simulation.

        Returns:
            A trajectory from the structures.
        """
        if constant_lattice:
            lattice = structures[0].lattice.matrix
        else:
            lattice = np.array([structure.lattice.matrix for structure in structures])

        species = structures[0].species
        frac_coords = [structure.frac_coords for structure in structures]
        site_properties = [structure.site_properties for structure in structures]

        return cls(
            lattice,
            species,  # type: ignore
            frac_coords,
            site_properties=site_properties,  # type: ignore
            constant_lattice=constant_lattice,
            **kwargs,
        )

    @classmethod
    def from_file(
        cls,
        filename: str | Path,
        constant_lattice: bool = True,
        **kwargs,
    ) -> Trajectory:
        """
        Create trajectory from XDATCAR or vasprun.xml file.

        Args:
            filename: Path to the file to read from.
            constant_lattice: Whether the lattice changes during the simulation,
                such as in an NPT MD simulation.

        Returns:
            A trajectory from the file.
        """
        fname = Path(filename).expanduser().resolve().name

        if fnmatch(fname, "*XDATCAR*"):
            structures = Xdatcar(filename).structures
        elif fnmatch(fname, "vasprun*.xml*"):
            structures = Vasprun(filename).structures
        else:
            supported = ("XDATCAR", "vasprun.xml")
            raise ValueError(f"Expect file to be one of {supported}; got {filename}.")

        return cls.from_structures(
            structures,
            constant_lattice=constant_lattice,
            **kwargs,
        )

    @staticmethod
    def _combine_lattice(lat1: np.ndarray, lat2: np.ndarray, len1: int, len2: int) -> tuple[np.ndarray, bool]:
        """
        Helper function to combine trajectory lattice.
        """
        if lat1.ndim == lat2.ndim == 2:
            constant_lat = True
            lat = lat1
        else:
            constant_lat = False
            if lat1.ndim == 2:
                lat1 = np.tile(lat1, (len1, 1, 1))
            if lat2.ndim == 2:
                lat2 = np.tile(lat2, (len2, 1, 1))
            lat = np.concatenate((lat1, lat2))

        return lat, constant_lat

    @staticmethod
    def _combine_site_props(
        prop1: SitePropsType | None, prop2: SitePropsType | None, len1: int, len2: int
    ) -> SitePropsType | None:
        """
        Combine site properties.

        Either one of prop1 or prop2 can be None, dict, or a list of dict. All
        possibilities of combining them are considered.
        """
        # special cases

        if prop1 is None and prop2 is None:
            return None

        if isinstance(prop1, dict) and prop1 == prop2:
            return prop1

        # general case

        assert prop1 is None or isinstance(prop1, (list, dict))
        assert prop2 is None or isinstance(prop2, (list, dict))

        p1_candidates = {
            "NoneType": [None] * len1,
            "dict": [prop1] * len1,
            "list": prop1,
        }
        p2_candidates = {
            "NoneType": [None] * len2,
            "dict": [prop2] * len2,
            "list": prop2,
        }
        p1_selected: list = p1_candidates[type(prop1).__name__]  # type: ignore
        p2_selected: list = p2_candidates[type(prop2).__name__]  # type: ignore

        return p1_selected + p2_selected

    @staticmethod
    def _combine_frame_props(prop1: list[dict] | None, prop2: list[dict] | None, len1: int, len2: int) -> list | None:
        """
        Combine frame properties.
        """
        if prop1 is None and prop2 is None:
            return None
        if prop1 is None:
            return [None] * len1 + list(prop2)  # type: ignore
        if prop2 is None:
            return list(prop1) + [None] * len2  # type: ignore
        return list(prop1) + list(prop2)  # type:ignore

    def _check_site_props(self, site_props: SitePropsType | None):
        """
        Check data shape of site properties.
        """
        if site_props is None:
            return

        if isinstance(site_props, dict):
            site_props = [site_props]
        else:
            assert len(site_props) == len(
                self
            ), f"Size of the site properties {len(site_props)} does not equal to the number of frames {len(self)}."

        num_sites = len(self.frac_coords[0])
        for d in site_props:
            for k, v in d.items():
                assert len(v) == num_sites, (
                    f"Size of site property {k} {len(v)}) does not equal to the "
                    f"number of sites in the structure {num_sites}."
                )

    def _check_frame_props(self, frame_props: list[dict] | None):
        """
        Check data shape of site properties.
        """
        if frame_props is None:
            return

        assert len(frame_props) == len(
            self
        ), f"Size of the frame properties {len(frame_props)} does not equal to the number of frames {len(self)}."

    def _get_site_props(self, frames: int | list[int]) -> SitePropsType | None:
        """
        Slice site properties.
        """
        if self.site_properties is None:
            return None
        if isinstance(self.site_properties, dict):
            return self.site_properties
        if isinstance(self.site_properties, list):
            if isinstance(frames, int):
                return self.site_properties[frames]
            if isinstance(frames, list):
                return [self.site_properties[i] for i in frames]
            raise ValueError("Unexpected frames type.")
        raise ValueError("Unexpected site_properties type.")
