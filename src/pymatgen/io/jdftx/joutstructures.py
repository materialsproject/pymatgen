"""Module for JOutStructures class.

This module contains the JOutStructures class for storing a series of
JOutStructure.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from pymatgen.io.jdftx.utils import correct_geom_iter_type, get_joutstructure_step_bounds, get_joutstructures_start_idx

if TYPE_CHECKING:
    import numpy as np

    from pymatgen.io.jdftx.jelstep import JElSteps
from pymatgen.io.jdftx.joutstructure import JOutStructure

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
    # TODO: Rename "iter_type" to "geom_opt_type"
    iter_type: str | None = None
    geom_converged: bool = False
    geom_converged_reason: str | None = None
    elec_converged: bool = False
    elec_converged_reason: str | None = None
    _t_s: float | None = None
    slices: list[JOutStructure] = field(default_factory=list)

    @classmethod
    def from_out_slice(cls, out_slice: list[str], iter_type: str = "IonicMinimize") -> JOutStructures:
        """Return JStructures object.

        Create a JStructures object from a slice of an out file's text
        corresponding to a single JDFTx call.

        Parameters
        ----------
        out_slice: list[str]
            A slice of a JDFTx out file (individual call of JDFTx)
        """
        instance = cls()
        if iter_type not in ["IonicMinimize", "LatticeMinimize"]:
            iter_type = correct_geom_iter_type(iter_type)
        instance.iter_type = iter_type
        start_idx = get_joutstructures_start_idx(out_slice)
        instance.set_joutstructure_list(out_slice[start_idx:])
        if instance.iter_type is None and len(instance) > 1:
            raise Warning("iter type interpreted as single-point calculation, but multiple structures found")
        instance.check_convergence()
        return instance

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
            if (self.iter_type in ["single point", None]) and (isinstance(self[-1].elecmindata[-1].t_s, float)):
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
    def etype(self) -> str:
        """
        Return etype from most recent JOutStructure.

        Return etype from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].etype
        raise AttributeError("Property etype inaccessible due to empty slices class field")

    @property
    def eiter_type(self) -> str:
        """
        Return eiter_type from most recent JOutStructure.

        Return eiter_type from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].eiter_type
        raise AttributeError("Property eiter_type inaccessible due to empty slices class field")

    @property
    def emin_flag(self) -> str:
        """
        Return emin_flag from most recent JOutStructure.

        Return emin_flag from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].emin_flag
        raise AttributeError("Property emin_flag inaccessible due to empty slices class field")

    @property
    def ecomponents(self) -> dict:
        """
        Return ecomponents from most recent JOutStructure.

        Return ecomponents from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].ecomponents
        raise AttributeError("Property ecomponents inaccessible due to empty slices class field")

    @property
    def elecmindata(self) -> JElSteps:
        """
        Return elecmindata from most recent JOutStructure.

        Return elecmindata from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].elecmindata
        raise AttributeError("Property elecmindata inaccessible due to empty slices class field")

    @property
    def stress(self) -> np.ndarray | None:
        """
        Return stress from most recent JOutStructure.

        Return stress from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].stress
        raise AttributeError("Property stress inaccessible due to empty slices class field")

    @property
    def strain(self) -> np.ndarray | None:
        """
        Return strain from most recent JOutStructure.

        Return strain from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].strain
        raise AttributeError("Property strain inaccessible due to empty slices class field")

    @property
    def nstep(self) -> int:
        """
        Return nstep from most recent JOutStructure.

        Return nstep from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].nstep
        raise AttributeError("Property nstep inaccessible due to empty slices class field")

    @property
    def e(self) -> float:
        """
        Return E from most recent JOutStructure.

        Return E from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].e
        raise AttributeError("Property E inaccessible due to empty slices class field")

    @property
    def grad_k(self) -> float | None:
        """
        Return grad_k from most recent JOutStructure.

        Return grad_k from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].grad_k
        raise AttributeError("Property grad_k inaccessible due to empty slices class field")

    @property
    def alpha(self) -> float | None:
        """
        Return alpha from most recent JOutStructure.

        Return alpha from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].alpha
        raise AttributeError("Property alpha inaccessible due to empty slices class field")

    @property
    def linmin(self) -> float | None:
        """
        Return linmin from most recent JOutStructure.

        Return linmin from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].linmin
        raise AttributeError("Property linmin inaccessible due to empty slices class field")

    @property
    def nelectrons(self) -> float:
        """
        Return nelectrons from most recent JOutStructure.

        Return nelectrons from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].nelectrons
        raise AttributeError("Property nelectrons inaccessible due to empty slices class field")

    @property
    def abs_magneticmoment(self) -> float | None:
        """
        Return abs_magneticmoment from most recent JOutStructure.

        Return abs_magneticmoment from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].abs_magneticmoment
        raise AttributeError("Property abs_magneticmoment inaccessible due to empty slices class field")

    @property
    def tot_magneticmoment(self) -> float | None:
        """
        Return tot_magneticmoment from most recent JOutStructure.

        Return tot_magneticmoment from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].tot_magneticmoment
        raise AttributeError("Property tot_magneticmoment inaccessible due to empty slices class field")

    @property
    def mu(self) -> float:
        """
        Return mu from most recent JOutStructure.

        Return mu from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].mu
        raise AttributeError("Property mu inaccessible due to empty slices class field")

    ###########################################################################
    # Electronic properties inherited from most recent JElSteps with symbol
    # disambiguation.
    ###########################################################################

    @property
    def elec_nstep(self) -> int:
        """Return the most recent electronic iteration.

        Return the most recent electronic iteration.

        Returns
        -------
        elec_nstep: int
        """
        if len(self.slices):
            return self.slices[-1].elec_nstep
        raise AttributeError("Property elec_nstep inaccessible due to empty slices class field")

    @property
    def elec_e(self) -> float:
        """Return the most recent electronic energy.

        Return the most recent electronic energy.

        Returns
        -------
        elec_e: float
        """
        if len(self.slices):
            return self.slices[-1].elec_e
        raise AttributeError("Property elec_e inaccessible due to empty slices class field")

    @property
    def elec_grad_k(self) -> float | None:
        """Return the most recent electronic grad_k.

        Return the most recent electronic grad_k.

        Returns
        -------
        grad_k: float
        """
        if len(self.slices):
            return self.slices[-1].grad_k
        raise AttributeError("Property grad_k inaccessible due to empty slices class field")

    @property
    def elec_alpha(self) -> float | None:
        """Return the most recent electronic alpha.

        Return the most recent electronic alpha.

        Returns
        -------
        alpha: float
        """
        if len(self.slices):
            return self.slices[-1].alpha
        raise AttributeError("Property alpha inaccessible due to empty slices class field")

    @property
    def elec_linmin(self) -> float | None:
        """Return the most recent electronic linmin.

        Return the most recent electronic linmin.

        Returns
        -------
        linmin: float
        """
        if len(self.slices):
            return self.slices[-1].linmin
        raise AttributeError("Property linmin inaccessible due to empty slices class field")

    def get_joutstructure_list(self, out_slice: list[str]) -> list[JOutStructure]:
        """Return list of JOutStructure objects.

        Set relevant variables for the JStructures object by parsing the
        out_slice.

        Parameters
        ----------
        out_slice: list[str]
            A slice of a JDFTx out file (individual call of JDFTx)
        """
        out_bounds = get_joutstructure_step_bounds(out_slice)
        return [
            JOutStructure.from_text_slice(out_slice[bounds[0] : bounds[1]], iter_type=self.iter_type)
            for bounds in out_bounds
        ]

    def set_joutstructure_list(self, out_slice: list[str]) -> None:
        """Set list of JOutStructure objects to slices.

        Set the list of JOutStructure objects to the slices attribute.

        Parameters
        ----------
        out_slice: list[str]
            A slice of a JDFTx out file (individual call of JDFTx)
        """
        out_list = self.get_joutstructure_list(out_slice)
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
        if len(self.slices):
            if not hasattr(self.slices[-1], name):
                raise AttributeError(f"{self.__class__.__name__} not found: {name}")
            return getattr(self.slices[-1], name)
        raise AttributeError(f"Property {name} inaccessible due to empty slices class field")

    def __dir__(self) -> list:
        """List attributes.

        Returns a list of attributes for the object, including those from
        self.slices[-1].

        Returns
        -------
        list
            A list of attribute names
        """
        # Get the default attributes
        default_attrs = dir(self)
        # Get the attributes from self.slices[-1] if slices is not empty
        slice_attrs = dir(self.slices[-1]) if self.slices else []
        # Combine and return unique attributes
        return list(set(default_attrs + slice_attrs))

    def __getitem__(self, key: int | str) -> JOutStructure | Any:
        """Return item.

        Return the value of an item.

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
        """Return JOutStructure object.

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
