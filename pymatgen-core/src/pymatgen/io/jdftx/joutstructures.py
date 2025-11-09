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
    import numpy as np
    from numpy.typing import NDArray

    from pymatgen.core.structure import Structure
    from pymatgen.io.jdftx.jelstep import JElSteps


from pymatgen.io.jdftx._output_utils import correct_geom_opt_type, is_lowdin_start_line
from pymatgen.io.jdftx.joutstructure import JOutStructure

__author__ = "Ben Rich"

# TODO: Break this up into data also stored in `properties` and `site_properties`
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
    "t_s",
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
        initial_structure (Structure | None): Cleaned pymatgen Structure object of initial JOutStructure
    """

    out_slice_start_flag = "-------- Electronic minimization -----------"
    opt_type: str | None = None
    geom_converged: bool = False
    geom_converged_reason: str | None = None
    elec_converged: bool = False
    elec_converged_reason: str | None = None
    t_s: float | None = None
    slices: list[JOutStructure] = field(default_factory=list, init=True)
    eopt_type: str | None = None
    etype: str | None = None
    emin_flag: str | None = None
    ecomponents: list[str] | None = None
    elecmindata: JElSteps | None = None
    stress: NDArray[np.float64] | None = None
    strain: NDArray[np.float64] | None = None
    forces: NDArray[np.float64] | None = None
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
    charges: NDArray[np.float64] | None = None
    magnetic_moments: NDArray[np.float64] | None = None
    selective_dynamics: list[int] | None = None
    structure: Structure | None = None
    initial_structure: Structure | None = None

    @classmethod
    def _from_out_slice(
        cls, out_slice: list[str], opt_type: str = "IonicMinimize", init_struc: Structure | None = None
    ) -> JOutStructures:
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
            _opt_type = correct_geom_opt_type(opt_type)
        start_idx = _get_joutstructures_start_idx(out_slice)
        slices = _get_joutstructure_list(out_slice[start_idx:], opt_type, init_structure=init_struc)
        return cls(slices=slices)

    def __post_init__(self):
        # Calling "self.slices[-1].opt_type" is safe since the "initial_structure" argument
        # ensure a length of at least 1, and opt_type has a default value
        self.opt_type = self.slices[-1].opt_type
        if self.opt_type is None and len(self) > 1:
            raise Warning("iter type interpreted as single-point calculation, but multiple structures found")
        # TODO: Change this to fetch from `properties` or `site_properties` instead of a slices attributes.
        for var in _joss_atrs_from_last_slice:
            val = None
            for i in range(1, len(self.slices) + 1):
                val = getattr(self.slices[-i], var)
                if val is not None:
                    break
            setattr(self, var, val)
        self.initial_structure = self._get_initial_structure()
        self.t_s = self._get_t_s()
        self._check_convergence()

    def _get_initial_structure(self) -> Structure | None:
        """Return initial structure.

        Returns:
            Structure | None: The initial structure.
        """
        if len(self):
            return self.slices[0].structure
        return None

    def _get_t_s(self) -> float | None:
        """Return time of calculation.

        Return time of calculation. Backup function to reference elecmindata for t_s for single point calculations.

        Returns:
            float: The total time in seconds for the calculation.
        """
        if self.t_s is not None:
            return self.t_s
        if len(self) and ((self.opt_type in ["single point", None]) and (isinstance(self.elecmindata.t_s, float))):
            return self.elecmindata.t_s
        return None

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


def _get_joutstructure_list(
    out_slice: list[str],
    opt_type: str | None,
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
    joutstructure_list: list[JOutStructure] = []
    if not len(out_bounds):
        # This case should only be called if no optimization steps are found to avoid errors down the line.
        # If this is changed to always be the first structure, logic down the line on what is considered
        # a "single point" calculation will be broken.
        joutstructure_list.append(JOutStructure._from_text_slice([], init_structure=init_structure, opt_type=opt_type))
    for i, bounds in enumerate(out_bounds):
        if i > 0:
            init_structure = joutstructure_list[-1]
        joutstructure = None
        # The final out_slice slice is protected by the try/except block, as this slice has a high
        # chance of being empty or malformed.
        try:
            joutstructure = JOutStructure._from_text_slice(
                out_slice[bounds[0] : bounds[1]],
                init_structure=init_structure,
                opt_type=opt_type,
            )
        except (ValueError, IndexError, TypeError, KeyError, AttributeError):
            if not i == len(out_bounds) - 1:
                raise
        if joutstructure is not None:
            joutstructure_list.append(joutstructure)
    return joutstructure_list
