"""Module for JOutStructures class.

This module contains the JOutStructures class for storing a series of
JOutStructure.

@mkhorton - this file is ready to review.
"""

from __future__ import annotations

import inspect
import pprint
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.core.units import bohr_to_ang
from pymatgen.io.jdftx.utils import correct_geom_opt_type, is_lowdin_start_line

if TYPE_CHECKING:
    from pymatgen.io.jdftx.jelstep import JElSteps
from pymatgen.io.jdftx.joutstructure import JOutStructure
from pymatgen.io.jdftx.utils import find_first_range_key

__author__ = "Ben Rich"


@dataclass
class JOutStructures:
    """Class for storing a series of JStructure objects.

    A class for storing a series of JStructure objects.

    Attributes
    ----------
    out_slice_start_flag: str
        The string that marks the beginning of the portion of an out file slice
        that contains data for a JOutStructures object.
    """

    out_slice_start_flag = "-------- Electronic minimization -----------"
    opt_type: str | None = None
    geom_converged: bool = False
    geom_converged_reason: str | None = None
    elec_converged: bool = False
    elec_converged_reason: str | None = None
    _t_s: float | None = None
    slices: list[JOutStructure] = field(default_factory=list)

    @classmethod
    def from_out_slice(cls, out_slice: list[str], opt_type: str = "IonicMinimize") -> JOutStructures:
        """Return JStructures object.

        Create a JStructures object from a slice of an out file's text
        corresponding to a single JDFTx call.

        Parameters
        ----------
        out_slice: list[str]
            A slice of a JDFTx out file (individual call of JDFTx)
        """
        instance = cls()
        if opt_type not in ["IonicMinimize", "LatticeMinimize"]:
            opt_type = correct_geom_opt_type(opt_type)
        instance.opt_type = opt_type
        start_idx = get_joutstructures_start_idx(out_slice)
        init_struc = instance.get_init_structure(out_slice[:start_idx])
        instance.set_joutstructure_list(out_slice[start_idx:], init_structure=init_struc)
        if instance.opt_type is None and len(instance) > 1:
            raise Warning("iter type interpreted as single-point calculation, but multiple structures found")
        instance.check_convergence()
        return instance

    # TODO: Move me to the correct spot
    def get_init_structure(self, pre_out_slice: list[str]) -> Structure | None:
        """Return initial structure.

        Return the initial structure from the pre_out_slice, corresponding to all data cut from JOutStructure list
        initialization. This is needed to ensure structural data that is not being updated (and therefore not being
        logged in the out file) is still available.

        Parameters
        ----------
        pre_out_slice: list[str]
            A slice of a JDFTx out file (individual call of JDFTx) that
            contains the initial structure information.
        """
        try:
            lat_mat = self.get_initial_lattice(pre_out_slice)
            coords = self.get_initial_coords(pre_out_slice)
            species = self.get_initial_species(pre_out_slice)
            return Structure(lattice=lat_mat, species=species, coords=coords)
        except AttributeError:
            return None

    def get_initial_lattice(self, pre_out_slice: list[str]) -> np.ndarray:
        """Return initial lattice.

        Return the initial lattice from the pre_out_slice, corresponding to all data cut from JOutStructure list
        initialization. This is needed to ensure lattice data that is not being updated (and therefore not being
        logged in the out file) is still available.

        Parameters
        ----------
        pre_out_slice: list[str]
            A slice of a JDFTx out file (individual call of JDFTx) that
            contains the initial lattice information.
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

    def get_initial_coords(self, pre_out_slice: list[str]) -> np.ndarray:
        """Return initial coordinates.

        Return the initial coordinates from the pre_out_slice, corresponding to all data cut from JOutStructure list
        initialization. This is needed to ensure coordinate data that is not being updated (and therefore not being
        logged in the out file) is still available.

        Parameters
        ----------
        pre_out_slice: list[str]
            A slice of a JDFTx out file (individual call of JDFTx) that
            contains the initial coordinates information.
        """
        lines = self._get_ion_lines(pre_out_slice)
        coords = np.zeros([len(lines), 3])
        for i, line in enumerate(lines):
            line_text = pre_out_slice[line].strip().split()[2:]
            for j in range(3):
                coords[i, j] = float(line_text[j])
        coords_type_lines = find_first_range_key("coords-type", pre_out_slice)
        if len(coords_type_lines):
            coords_type = pre_out_slice[coords_type_lines[0]].strip().split()[1]
            if coords_type.lower() != "cartesian":
                coords = np.dot(coords, self.get_initial_lattice(pre_out_slice))
        return coords

    def get_initial_species(self, pre_out_slice: list[str]) -> list[str]:
        """Return initial species.

        Return the initial species from the pre_out_slice, corresponding to all data cut from JOutStructure list
        initialization. This is needed to ensure species data that is not being updated (and therefore not being
        logged in the out file) is still available.

        Parameters
        ----------
        pre_out_slice: list[str]
            A slice of a JDFTx out file (individual call of JDFTx) that
            contains the initial species information.
        """
        lines = self._get_ion_lines(pre_out_slice)
        species_strs = []
        for line in lines:
            species_strs.append(pre_out_slice[line].strip().split()[1])
        return species_strs

    def _get_ion_lines(self, pre_out_slice: list[str]) -> list[int]:
        """Return ion lines.

        Return the ion lines from the pre_out_slice, ensuring that all the ion lines are consecutive.

        Parameters
        ----------
        pre_out_slice: list[str]
            A slice of a JDFTx out file (individual call of JDFTx) that
            contains the ion lines information.
        """
        _lines = find_first_range_key("ion ", pre_out_slice)
        if not len(_lines):
            raise AttributeError("Ion lines not found in pre_out_slice")
        gaps = [_lines[i + 1] - _lines[i] for i in range(len(_lines) - 1)]
        if not all(g == 1 for g in gaps):
            # TODO: Write the fix for this case
            raise AttributeError("Ion lines not consecutive in pre_out_slice")
        return _lines

    # TODO: This currently returns the most recent t_s, which is not at all helpful.
    # Correct this to be the total time in seconds for the series of structures.
    @property
    def t_s(self) -> float | None:
        """Return time of calculation.

        Return the total time in seconds for the calculation.

        Returns
        -------
        t_s: float
            The total time in seconds for the calculation
        """
        if self._t_s is not None:
            return self._t_s
        if len(self):
            if (self.opt_type in ["single point", None]) and (isinstance(self[-1].elecmindata[-1].t_s, float)):
                self._t_s = self[-1].elecmindata[-1].t_s
            elif isinstance(self[-1].t_s, float):
                self._t_s = self[-1].t_s
            else:
                raise AttributeError("t_s not set in most recent JOutStructure")
        return self._t_s

    ###########################################################################
    # Properties inherited from most recent JOutStructure
    ###########################################################################

    @property
    def etype(self) -> str | None:
        """
        Return etype from most recent JOutStructure.

        Return etype from most recent JOutStructure, where etype corresponds to the string
        representation of the ensemble potential - (ie "F" for Helmholtz, "G" for Grand-Canonical Potential).
        """
        if len(self.slices):
            return self.slices[-1].etype
        raise AttributeError("Property etype inaccessible due to empty slices class field")

    @property
    def eopt_type(self) -> str | None:
        """
        Return eopt_type from most recent JOutStructure.

        Return eopt_type from most recent JOutStructure, where eopt_type corresponds to the
        JDFTx string representation for the minimization program used to minimize the electron density.
        """
        if len(self.slices):
            return self.slices[-1].eopt_type
        raise AttributeError("Property eopt_type inaccessible due to empty slices class field")

    @property
    def emin_flag(self) -> str | None:
        """
        Return emin_flag from most recent JOutStructure.

        Return emin_flag from most recent JOutStructure, where emin_flag corresponds to the
        flag string used to mark the beginning of a section of the out file containing the
        data to construct a JOutStructure object.
        """
        if len(self.slices):
            return self.slices[-1].emin_flag
        raise AttributeError("Property emin_flag inaccessible due to empty slices class field")

    @property
    def ecomponents(self) -> dict | None:
        """
        Return ecomponents from most recent JOutStructure.

        Return ecomponents from most recent JOutStructure, where ecomponents is a dictionary
        mapping string representation of system energy types to their values in eV.
        """
        if len(self.slices):
            return self.slices[-1].ecomponents
        raise AttributeError("Property ecomponents inaccessible due to empty slices class field")

    @property
    def elecmindata(self) -> JElSteps | None:
        """
        Return elecmindata from most recent JOutStructure.

        Return elecmindata from most recent JOutStructure, where elecmindata is a JElSteps object
        created to hold electronic minimization data on the electronic density for this JOuStructure.
        """
        if len(self.slices):
            return self.slices[-1].elecmindata
        raise AttributeError("Property elecmindata inaccessible due to empty slices class field")

    # TODO: Figure out how JDFTx defines the equilibrium lattice parameters and
    # incorporate into this docustring.
    @property
    def stress(self) -> np.ndarray | None:
        """
        Return stress from most recent JOutStructure.

        Return stress from most recent JOutStructure, where stress is the 3x3 unitless
        stress tensor.
        """
        if len(self.slices):
            return self.slices[-1].stress
        raise AttributeError("Property stress inaccessible due to empty slices class field")

    @property
    def strain(self) -> np.ndarray | None:
        """
        Return strain from most recent JOutStructure.

        Return strain from most recent JOutStructure, where strain is the 3x3 strain
        tensor in units eV/A^3.
        """
        if len(self.slices):
            return self.slices[-1].strain
        raise AttributeError("Property strain inaccessible due to empty slices class field")

    @property
    def nstep(self) -> int | None:
        """
        Return nstep from most recent JOutStructure.

        Return nstep from most recent JOutStructure, where nstep corresponds to the step
        number of the geometric optimization.
        """
        if len(self.slices):
            return self.slices[-1].nstep
        raise AttributeError("Property nstep inaccessible due to empty slices class field")

    @property
    def e(self) -> float | None:
        """
        Return e from most recent JOutStructure.

        Return e from most recent JOutStructure, where e corresponds to the system energy
        of the system's "etype" in eV.
        """
        if len(self.slices):
            return self.slices[-1].e
        raise AttributeError("Property e inaccessible due to empty slices class field")

    @property
    def grad_k(self) -> float | None:
        """
        Return grad_k from most recent JOutStructure.

        Return grad_k from most recent JOutStructure, where grad_k corresponds to the geometric
        gradient along the geometric line minimization.
        """
        if len(self.slices):
            return self.slices[-1].grad_k
        raise AttributeError("Property grad_k inaccessible due to empty slices class field")

    @property
    def alpha(self) -> float | None:
        """
        Return alpha from most recent JOutStructure.

        Return alpha from most recent JOutStructure, where alpha corresponds to the geometric
        step size along the geometric line minimization.
        """
        if len(self.slices):
            return self.slices[-1].alpha
        raise AttributeError("Property alpha inaccessible due to empty slices class field")

    @property
    def linmin(self) -> float | None:
        """
        Return linmin from most recent JOutStructure.

        Return linmin from most recent JOutStructure, where linmin corresponds to the normalized
        projection of the geometric gradient to the step direction within the line minimization.
        """
        if len(self.slices):
            return self.slices[-1].linmin
        raise AttributeError("Property linmin inaccessible due to empty slices class field")

    @property
    def nelectrons(self) -> float | None:
        """
        Return nelectrons from most recent JOutStructure.

        Return nelectrons from most recent JOutStructure, where nelectrons corresponds to the
        number of electrons in the electron density.
        """
        if len(self.slices):
            return self.slices[-1].nelectrons
        raise AttributeError("Property nelectrons inaccessible due to empty slices class field")

    @property
    def abs_magneticmoment(self) -> float | None:
        """
        Return abs_magneticmoment from most recent JOutStructure.

        Return abs_magneticmoment from most recent JOutStructure, where abs_magneticmoment corresponds
        to the absolute magnetic moment of the electron density.
        """
        if len(self.slices):
            return self.slices[-1].abs_magneticmoment
        raise AttributeError("Property abs_magneticmoment inaccessible due to empty slices class field")

    @property
    def tot_magneticmoment(self) -> float | None:
        """
        Return tot_magneticmoment from most recent JOutStructure.

        Return tot_magneticmoment from most recent JOutStructure, where tot_magneticmoment corresponds
        to the total magnetic moment of the electron density.
        """
        if len(self.slices):
            return self.slices[-1].tot_magneticmoment
        raise AttributeError("Property tot_magneticmoment inaccessible due to empty slices class field")

    @property
    def mu(self) -> float | None:
        """
        Return mu from most recent JOutStructure.

        Return mu from most recent JOutStructure, where mu corresponds to the electron chemical potential
        (Fermi level) in eV.
        """
        if len(self.slices):
            return self.slices[-1].mu
        raise AttributeError("Property mu inaccessible due to empty slices class field")

    ###########################################################################
    # Electronic properties inherited from most recent JElSteps with symbol
    # disambiguation.
    ###########################################################################

    @property
    def elec_nstep(self) -> int | None:
        """Return the most recent electronic step number.

        Return the most recent elec_nstep, where elec_nstep corresponds to the SCF
        step number.

        Returns
        -------
        elec_nstep: int
        """
        if len(self.slices):
            return self.slices[-1].elec_nstep
        raise AttributeError("Property elec_nstep inaccessible due to empty slices class field")

    @property
    def elec_e(self) -> float | None:
        """Return the most recent elec_e.

        Return the most recent elec_e, where elec_e corresponds to the system's "etype"
        energy as printed within the SCF log.

        Returns
        -------
        elec_e: float
        """
        if len(self.slices):
            return self.slices[-1].elec_e
        raise AttributeError("Property elec_e inaccessible due to empty slices class field")

    @property
    def elec_grad_k(self) -> float | None:
        """Return the most recent elec_grad_k.

        Return the most recent elec_grad_k, where elec_grad_k corresponds to the electronic
        gradient along the line minimization (equivalent to grad_k for a JElSteps object).

        Returns
        -------
        grad_k: float
        """
        if len(self.slices):
            return self.slices[-1].elec_grad_k
        raise AttributeError("Property grad_k inaccessible due to empty slices class field")

    @property
    def elec_alpha(self) -> float | None:
        """Return the most recent elec_alpha.

        Return the most recent elec_alpha, where elec_alpha corresponds to the step size
        of the electronic optimization (equivalent to alpha for a JElSteps object).

        Returns
        -------
        alpha: float
        """
        if len(self.slices):
            return self.slices[-1].elec_alpha
        raise AttributeError("Property alpha inaccessible due to empty slices class field")

    @property
    def elec_linmin(self) -> float | None:
        """Return the most recent elec_linmin.

        Return the most recent elec_linmin, where elec_linmin corresponds to the normalized
        projection of the electronic gradient on the electronic line minimization direction
        (equivalent to linmin for a JElSteps object).

        Returns
        -------
        linmin: float
        """
        if len(self.slices):
            return self.slices[-1].elec_linmin
        raise AttributeError("Property linmin inaccessible due to empty slices class field")

    def get_joutstructure_list(
        self, out_slice: list[str], init_structure: Structure | None = None
    ) -> list[JOutStructure]:
        """Return list of JOutStructure objects.

        Get list of JStructure objects by splitting out_slice into slices and constructing
        a JOuStructure object for each slice. Used in initialization.

        Parameters
        ----------
        out_slice: list[str]
            A slice of a JDFTx out file (individual call of JDFTx)
        """
        out_bounds = get_joutstructure_step_bounds(out_slice)
        joutstructure_list: list[Structure | JOutStructure] = []
        for i, bounds in enumerate(out_bounds):
            if i > 0:
                init_structure = joutstructure_list[-1]
            joutstructure_list.append(
                JOutStructure.from_text_slice(
                    out_slice[bounds[0] : bounds[1]],
                    init_structure=init_structure,
                    opt_type=self.opt_type,
                )
            )
        return joutstructure_list

    def set_joutstructure_list(self, out_slice: list[str], init_structure: Structure | None = None) -> None:
        """Set list of JOutStructure objects to slices.

        Set the list of JOutStructure objects to the slices attribute.

        Parameters
        ----------
        out_slice: list[str]
            A slice of a JDFTx out file (individual call of JDFTx)
        """
        out_list = self.get_joutstructure_list(out_slice, init_structure=init_structure)
        for jos in out_list:
            self.slices.append(jos)

    def check_convergence(self) -> None:
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

    # This method is likely never going to be called as all (currently existing)
    # attributes of the most recent slice are explicitly defined as a class
    # property. However, it is included to reduce the likelihood of errors
    # upon future changes to downstream code.
    def __getattr__(self, name: str) -> Any:
        """Return attribute value.

        Return the value of an attribute.

        Parameters
        ----------
        name: str
            The name of the attribute

        Returns
        -------
        value
            The value of the attribute
        """
        if name in self.__dict__:
            return self.__dict__[name]

        # Check if the attribute is a property of the class
        for cls in inspect.getmro(self.__class__):
            if name in cls.__dict__ and isinstance(cls.__dict__[name], property):
                return cls.__dict__[name].__get__(self)

        # Check if the attribute is in self.jstrucs
        if hasattr(self.slices[-1], name):
            return getattr(self.slices[-1], name)

        # If the attribute is not found in either, raise an AttributeError
        raise AttributeError(f"{self.__class__.__name__} not found: {name}")

    def __getitem__(self, key: int | str) -> JOutStructure | Any:
        """Return item.

        Return the value of an item given an integer or string as a key (Otherwise
        returns None).

        Parameters
        ----------
        key: int | str
            The key of the item

        Returns
        -------
        val
            The value of the item
        """
        val = None
        if type(key) is int:
            val = self.getitem_int(key)
        if type(key) is str:
            val = self.getitem_str(key)
        return val

    def getitem_int(self, key: int) -> JOutStructure:
        """Return a JOutStructure object.

        Return the JOutStructure object at the key index.

        Parameters
        ----------
        key: int
            The index of the JOutStructure object

        Returns
        -------
        joutstructure: JOutStructure
            The JOutStructure object at the key index
        """
        return self.slices[key]

    def getitem_str(self, key: str) -> Any:
        """Return attribute value.

        Return the value of an attribute.

        Parameters
        ----------
        key: str
            The name of the attribute

        Returns
        -------
        value
            The value of the attribute
        """
        return getattr(self, key)

    def __len__(self) -> int:
        """Return length of JOutStructures object.

        Returns the number of geometric optimization steps in the
        JOutStructures object.

        Returns
        -------
        length: int
            The number of geometric optimization steps in the JOutStructures
            object
        """
        return len(self.slices)

    def __str__(self) -> str:
        """Return string representation.

        Return a string representation of the JOutStructures object.

        Returns
        -------
        str
            A string representation of the JOutStructures object
        """
        return pprint.pformat(self)


elec_min_start_flag: str = "-------- Electronic minimization -----------"


def get_joutstructure_step_bounds(
    out_slice: list[str],
    out_slice_start_flag: str = elec_min_start_flag,
) -> list[list[int]]:
    """Return list of boundary indices for each structure in out_slice.

    Return a list of lists of integers where each sublist contains the start and end
    of an individual optimization step (or SCF cycle if no optimization).

    Parameters
    ----------
    out_slice: list[str]
        A slice of a JDFTx out file (individual call of JDFTx)

    Returns
    -------
    bounds_list: list[list[int, int]]
        A list of lists of integers where each sublist contains the start and end
        of an individual optimization step (or SCF cycle if no optimization)
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


def get_joutstructures_start_idx(
    out_slice: list[str],
    out_slice_start_flag: str = elec_min_start_flag,
) -> int | None:
    """Return index of first line of first structure.

    Return the index of the first line of the first structure in the out_slice.

    Parameters
    ----------
    out_slice: list[str]
        A slice of a JDFTx out file (individual call of JDFTx)

    Returns
    -------
    i: int
        The index of the first line of the first structure in the out_slice
    """
    for i, line in enumerate(out_slice):
        if out_slice_start_flag in line:
            return i
    return None
