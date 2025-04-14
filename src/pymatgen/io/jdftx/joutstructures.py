"""Module for JOutStructures class.

This module contains the JOutStructures class for storing a series of
JOutStructure.
"""

from __future__ import annotations

import pprint
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from monty.dev import deprecated

if TYPE_CHECKING:
    from pymatgen.io.jdftx.jelstep import JElSteps

import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.core.units import bohr_to_ang
from pymatgen.io.jdftx._output_utils import correct_geom_opt_type, find_first_range_key, is_lowdin_start_line
from pymatgen.io.jdftx.joutstructure import JOutStructure

__author__ = "Ben Rich"

_joss_atrs_from_last_slice = (
    "etype",
    "eopt_type",
    "emin_flag",
    "ecomponents",
    "elecmindata",
    "stress",
    "strain",
    "forces",
    "nstep",
    "e",
    "grad_k",
    "alpha",
    "linmin",
    "nelectrons",
    "abs_magneticmoment",
    "tot_magneticmoment",
    "mu",
    "elec_nstep",
    "elec_e",
    "elec_grad_k",
    "elec_alpha",
    "elec_linmin",
    "charges",
    "magnetic_moments",
    "selective_dynamics",
    "structure",
)


@dataclass
class JOutStructures:
    """
    Class for storing a series of JStructure objects.

    A class for storing a series of JStructure objects.

    Attributes:
        out_slice_start_flag (str): The string that marks the beginning of the portion of an out file slice
            that contains data for a JOutStructures object.
        opt_type (str | None): The type of optimization performed on the structures in the JOutStructures object.
        geom_converged (bool): Whether the geometry of the last structure in the list has converged.
        geom_converged_reason (str | None): The reason the geometry of the last structure in the list has converged.
        elec_converged (bool): Whether the electronic density of the last structure in the list has converged.
        elec_converged_reason (str | None): The reason the electronic density of the last structure in the list has
            converged.
        slices (list[JOutStructure]): A list of JOutStructure objects.
        eopt_type (str | None): The type of electronic optimization performed on the last structure in the list.
        etype (str | None): String representation of total energy-type of system. Commonly "G"
        (grand-canonical potential) for GC calculations, and "F" for canonical (fixed electron count) calculations.
        emin_flag (str | None): The flag for the electronic minimization.
        ecomponents (list[str] | None): The components of the electronic minimization.
        elecmindata (JElSteps): The electronic minimization data.
        stress (np.ndarray | None): The stress tensor.
        strain (np.ndarray | None): The strain tensor.
        nstep (int | None): The number of steps in the optimization.
        e (float | None): The total energy.
        grad_k (float | None): The final norm of the preconditioned gradient for geometric optimization of the most
            recent JDFTx call (evaluated as dot(g, Kg), where g is the gradient and Kg is the preconditioned gradient).
            (written as "|grad|_K" in JDFTx output).
        alpha (float | None): The step size of the final geometric step in the most recent JDFTx call.
        linmin (float | None): The final normalized projection of the geometric step direction onto the gradient for
            the most recent JDFTx call.
        abs_magneticmoment (float | None): The absolute magnetic moment of the most recent JDFTx call.
        tot_magneticmoment (float | None): The total magnetic moment of the most recent JDFTx call.
        mu (float | None): The Fermi energy of the most recent JDFTx call.
        elec_e (float) | None: The final energy of the most recent electronic optimization step.
        elec_nstep (int): The number of electronic optimization steps in the most recent JDFTx call.
        elec_grad_k (float | None): The final norm of the preconditioned gradient for electronic optimization of the
            most recent JDFTx call (evaluated as dot(g, Kg), where g is the gradient and Kg is the preconditioned
            gradient). (written as "|grad|_K" in JDFTx output).
        elec_alpha (float) | None: The step size of the final electronic step in the most recent JDFTx call.
        elec_linmin (float | None): The final normalized projection of the electronic step direction onto the gradient
            for the most recent JDFTx call.
        charges (np.ndarray[float] | None): The most recent Lowdin-charges.
        magnetic_moments (np.ndarray[float] | None): The most recent Lowdin-magnetic moments.
        selective_dynamics (list[int] | None): The selective dynamics flags for the most recent JDFTx call.
        structure (Structure | None): Cleaned pymatgen Structure object of final JOutStructure
    """

    out_slice_start_flag = "-------- Electronic minimization -----------"
    opt_type: str | None = None
    geom_converged: bool = False
    geom_converged_reason: str | None = None
    elec_converged: bool = False
    elec_converged_reason: str | None = None
    _t_s: float | None = None
    slices: list[JOutStructure] = field(default_factory=list, init=True)
    eopt_type: str | None = None
    etype: str | None = None
    emin_flag: str | None = None
    ecomponents: list[str] | None = None
    elecmindata: JElSteps = None
    stress: np.ndarray | None = None
    strain: np.ndarray | None = None
    forces: np.ndarray | None = None
    nstep: int | None = None
    e: float | None = None
    grad_k: float | None = None
    alpha: float | None = None
    linmin: float | None = None
    nelectrons: float | None = None
    abs_magneticmoment: float | None = None
    tot_magneticmoment: float | None = None
    mu: float | None = None
    elec_nstep: int | None = None
    elec_e: float | None = None
    elec_grad_k: float | None = None
    elec_alpha: float | None = None
    elec_linmin: float | None = None
    charges: np.ndarray[float] | None = None
    magnetic_moments: np.ndarray[float] | None = None
    selective_dynamics: list[int] | None = None
    structure: Structure | None = None

    @classmethod
    def _from_out_slice(cls, out_slice: list[str], opt_type: str = "IonicMinimize") -> JOutStructures:
        """
        Return JStructures object.

        Create a JStructures object from a slice of an out file's text
        corresponding to a single JDFTx call.

        Args:
            out_slice (list[str]): A slice of a JDFTx out file (individual call of JDFTx).

        Returns:
            JOutStructures: The created JOutStructures object.
        """
        if opt_type not in ["IonicMinimize", "LatticeMinimize"]:
            opt_type = correct_geom_opt_type(opt_type)
        start_idx = _get_joutstructures_start_idx(out_slice)
        init_struc = _get_init_structure(out_slice[:start_idx])
        slices = _get_joutstructure_list(out_slice[start_idx:], opt_type, init_structure=init_struc)
        return cls(slices=slices)

    def __post_init__(self):
        self.opt_type = self.slices[-1].opt_type
        if self.opt_type is None and len(self) > 1:
            raise Warning("iter type interpreted as single-point calculation, but multiple structures found")
        self.t_s = self._get_t_s()
        for var in _joss_atrs_from_last_slice:
            setattr(self, var, getattr(self.slices[-1], var))
        self._check_convergence()

    def _get_t_s(self) -> float | None:
        """Return time of calculation.

        Returns:
            float: The total time in seconds for the calculation.
        """
        if self._t_s is not None:
            return self._t_s
        if len(self):
            if (self.opt_type in ["single point", None]) and (isinstance(self[-1].elecmindata[-1].t_s, float)):
                self._t_s = self[-1].elecmindata[-1].t_s
            else:
                self._t_s = self[-1].t_s
        return self._t_s

    def _check_convergence(self) -> None:
        """Set convergence flags.

        Check if the geometry and electronic density of last structure in the
        list has converged.
        """
        jst = self.slices[-1]
        if jst.elecmindata.converged:
            self.elec_converged = True
            self.elec_converged_reason = jst.elecmindata.converged_reason
        if jst.geom_converged:
            self.geom_converged = True
            self.geom_converged_reason = jst.geom_converged_reason

    def as_dict(self) -> dict:
        """
        Convert the JOutStructures object to a dictionary.

        Returns:
            dict: A dictionary representation of the JOutStructures object.
        """
        dct = {}
        for fld in self.__dataclass_fields__:
            if fld == "slices":
                dct[fld] = [slc.as_dict() for slc in self.slices]
                continue
            value = getattr(self, fld)
            if hasattr(value, "as_dict"):
                dct[fld] = value.as_dict()
            else:
                dct[fld] = value
        return dct

    @deprecated(as_dict, deadline=(2025, 10, 4))
    def to_dict(self):
        return self.as_dict()

    def __getitem__(self, key: int | str) -> JOutStructure | Any:
        """Return item.

        Args:
            key (int | str): The key of the item.

        Returns:
            JOutStructure | Any: The value of the item.
        """
        val = None
        if type(key) is int:
            val = self._getitem_int(key)
        if type(key) is str:
            val = self._getitem_str(key)
        return val

    def _getitem_int(self, key: int) -> JOutStructure:
        """Return a JOutStructure object.

        Args:
            key (int): The index of the JOutStructure object.

        Returns:
            JOutStructure: The JOutStructure object at the key index.
        """
        return self.slices[key]

    def _getitem_str(self, key: str) -> Any:
        """Return attribute value.

        Args:
            key (str): The name of the attribute.

        Returns:
            Any: The value of the attribute.
        """
        return getattr(self, key)

    def __len__(self) -> int:
        """Return length of JOutStructures object.

        Returns:
            int: The number of geometric optimization steps in the JOutStructures object.
        """
        return len(self.slices)

    def __str__(self) -> str:
        """Return string representation.

        Returns:
            str: A string representation of the JOutStructures object.
        """
        return pprint.pformat(self)


_elec_min_start_flag: str = "-------- Electronic minimization -----------"


def _get_joutstructure_step_bounds(
    out_slice: list[str],
    out_slice_start_flag: str = _elec_min_start_flag,
) -> list[list[int]]:
    """Return list of boundary indices for each structure in out_slice.

    Args:
        out_slice (list[str]): A slice of a JDFTx out file (individual call of JDFTx).
        out_slice_start_flag (str): The string that marks the beginning of the portion of an out file slice
            that contains data for a JOutStructures object.

    Returns:
        list[list[int]]: A list of lists of integers where each sublist contains the start and end
        of an individual optimization step (or SCF cycle if no optimization).
    """
    bounds_list = []
    bounds = None
    end_started = False
    for i, line in enumerate(out_slice):
        if not end_started:
            if out_slice_start_flag in line:
                bounds = [i]
            elif (bounds is not None) and (is_lowdin_start_line(line)):
                end_started = True
        elif not len(line.strip()) and bounds is not None:
            bounds.append(i)
            bounds_list.append(bounds)
            bounds = None
            end_started = False
    # This case is for jdftx calls that were interrupted prematurely
    if bounds is not None:
        bounds.append(len(out_slice) - 1)
        bounds_list.append(bounds)
    return bounds_list


def _get_joutstructures_start_idx(
    out_slice: list[str],
    out_slice_start_flag: str = _elec_min_start_flag,
) -> int | None:
    """Return index of first line of first structure.

    Args:
        out_slice (list[str]): A slice of a JDFTx out file (individual call of JDFTx).
        out_slice_start_flag (str): The string that marks the beginning of the portion of an out file slice
            that contains data for a JOutStructures object.

    Returns:
        int | None: The index of the first line of the first structure in the out_slice.
    """
    for i, line in enumerate(out_slice):
        if out_slice_start_flag in line:
            return i
    return None


def _get_init_structure(pre_out_slice: list[str]) -> Structure | None:
    """
    Return initial structure.

    Return the initial structure from the pre_out_slice, corresponding to all data cut from JOutStructure list
    initialization. This is needed to ensure structural data that is not being updated (and therefore not being
    logged in the out file) is still available.

    Args:
        pre_out_slice (list[str]): A slice of a JDFTx out file (individual call of JDFTx) that
        contains the initial structure information.

    Returns:
        Structure | None: The initial structure if available, otherwise None.
    """
    try:
        lat_mat = _get_initial_lattice(pre_out_slice)
        coords = _get_initial_coords(pre_out_slice)
        species = _get_initial_species(pre_out_slice)
        return Structure(lattice=lat_mat, species=species, coords=coords, coords_are_cartesian=True)
    except AttributeError:
        return None


def _get_initial_lattice(pre_out_slice: list[str]) -> np.ndarray:
    """Return initial lattice.

    Return the initial lattice from the pre_out_slice, corresponding to all data cut from JOutStructure list
    initialization. This is needed to ensure lattice data that is not being updated (and therefore not being
    logged in the out file) is still available.

    Args:
        pre_out_slice (list[str]): A slice of a JDFTx out file (individual call of JDFTx) that
            contains the initial lattice information.

    Returns:
        np.ndarray: The initial lattice matrix.
    """
    lat_lines = find_first_range_key("lattice  ", pre_out_slice)
    if len(lat_lines):
        lat_line = lat_lines[0]
        lat_mat = np.zeros([3, 3])
        for i in range(3):
            line_text = pre_out_slice[lat_line + i + 1].strip().split()
            for j in range(3):
                lat_mat[i, j] = float(line_text[j])
        return lat_mat.T * bohr_to_ang
    raise AttributeError("Lattice not found in pre_out_slice")


def _get_initial_coords(pre_out_slice: list[str]) -> np.ndarray:
    """Return initial coordinates.

    Return the initial coordinates from the pre_out_slice, corresponding to all data cut from JOutStructure list
    initialization. This is needed to ensure coordinate data that is not being updated (and therefore not being
    logged in the out file) is still available.

    Args:
        pre_out_slice (list[str]): A slice of a JDFTx out file (individual call of JDFTx) that
            contains the initial coordinates information.

    Returns:
        np.ndarray: The initial coordinates.
    """
    lines = _get_ion_lines(pre_out_slice)
    coords = np.zeros([len(lines), 3])
    for i, line in enumerate(lines):
        line_text = pre_out_slice[line].strip().split()[2:]
        for j in range(3):
            coords[i, j] = float(line_text[j])
    coords_type_lines = find_first_range_key("coords-type", pre_out_slice)
    if len(coords_type_lines):
        coords_type = pre_out_slice[coords_type_lines[0]].strip().split()[1]
        if coords_type.lower() != "cartesian":
            coords = np.dot(coords, _get_initial_lattice(pre_out_slice))
        else:
            coords *= bohr_to_ang
    return coords


def _get_initial_species(pre_out_slice: list[str]) -> list[str]:
    """Return initial species.

    Return the initial species from the pre_out_slice, corresponding to all data cut from JOutStructure list
    initialization. This is needed to ensure species data that is not being updated (and therefore not being
    logged in the out file) is still available.

    Args:
        pre_out_slice (list[str]): A slice of a JDFTx out file (individual call of JDFTx) that
            contains the initial species information.

    Returns:
        list[str]: The initial species.
    """
    lines = _get_ion_lines(pre_out_slice)
    species_strs = []
    for line in lines:
        species_strs.append(pre_out_slice[line].strip().split()[1])
    return species_strs


def _get_ion_lines(pre_out_slice: list[str]) -> list[int]:
    """Return ion lines.

    Return the ion lines from the pre_out_slice, ensuring that all the ion lines are consecutive.

    Args:
        pre_out_slice (list[str]): A slice of a JDFTx out file (individual call of JDFTx) that
            contains the ion lines information.

    Returns:
        list[int]: The ion lines.
    """
    _lines = find_first_range_key("ion ", pre_out_slice)
    if not len(_lines):
        raise AttributeError("Ion lines not found in pre_out_slice")
    gaps = [_lines[i + 1] - _lines[i] for i in range(len(_lines) - 1)]
    if not all(g == 1 for g in gaps):
        # TODO: Write the fix for this case
        raise AttributeError("Ion lines not consecutive in pre_out_slice")
    return _lines


def _get_joutstructure_list(
    out_slice: list[str],
    opt_type: str,
    init_structure: Structure | None = None,
) -> list[JOutStructure]:
    """Return list of JOutStructure objects.

    Get list of JStructure objects by splitting out_slice into slices and constructing
    a JOutStructure object for each slice. Used in initialization.

    Args:
        out_slice (list[str]): A slice of a JDFTx out file (individual call of JDFTx).
        init_structure (Structure | None): The initial structure if available, otherwise None.

    Returns:
        list[JOutStructure]: The list of JOutStructure objects.
    """
    out_bounds = _get_joutstructure_step_bounds(out_slice)
    joutstructure_list: list[Structure | JOutStructure] = []
    for i, bounds in enumerate(out_bounds):
        if i > 0:
            init_structure = joutstructure_list[-1]
        joutstructure_list.append(
            JOutStructure._from_text_slice(
                out_slice[bounds[0] : bounds[1]],
                init_structure=init_structure,
                opt_type=opt_type,
            )
        )
    return joutstructure_list
