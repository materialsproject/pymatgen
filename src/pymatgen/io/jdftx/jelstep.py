"""Module for parsing single SCF step from JDFTx.

This module contains the JElStep class for parsing single SCF step from a JDFTx out file.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Any, ClassVar

from pymatgen.core.units import Ha_to_eV
from pymatgen.io.jdftx.utils import gather_JElSteps_line_collections, get_colon_var_t1

__author__ = "Ben Rich"


@dataclass
class JElStep:
    """Electronic minimization data for a single SCF step.

    Class object for storing logged electronic minimization data for a single
    SCF step.

    Attributes
    iter_type: str | None
        The type of electronic minimization step (almost always ElecMinimize)

    etype: str | None
        The type of energy component (G, F, or Etot)


    """

    iter_type: str | None = None
    etype: str | None = None
    nstep: int | None = None
    e: float | None = None
    grad_k: float | None = None
    alpha: float | None = None
    linmin: float | None = None
    t_s: float | None = None
    mu: float | None = None
    nelectrons: float | None = None
    abs_magneticmoment: float | None = None
    tot_magneticmoment: float | None = None
    subspacerotationadjust: float | None = None
    converged: bool = False
    converged_reason: str | None = None

    @classmethod
    def from_lines_collect(cls, lines_collect: list[str], iter_type: str, etype: str) -> JElStep:
        """Return JElStep object.

        Create a JElStep object from a list of lines of text from a JDFTx out
        file corresponding to a single SCF step.

        Parameters
        ----------
        lines_collect: list[str]
            A list of lines of text from a JDFTx out file corresponding to a
            single SCF step
        iter_type: str
            The type of electronic minimization step
        etype: str
            The type of energy component
        """
        instance = cls()
        instance.iter_type = iter_type
        instance.etype = etype
        _iter_flag = f"{iter_type}: Iter: "
        for i, line_text in enumerate(lines_collect):
            if instance.is_iter_line(i, line_text, _iter_flag):
                instance.read_iter_line(line_text)
            elif instance.is_fillings_line(i, line_text):
                instance.read_fillings_line(line_text)
            elif instance.is_subspaceadjust_line(i, line_text):
                instance.read_subspaceadjust_line(line_text)
        return instance

    def is_iter_line(self, i: int, line_text: str, _iter_flag: str) -> bool:
        """Return True if opt iter line.

        Return True if the line_text is the start of a log message for a
        JDFTx optimization step.

        Parameters
        ----------
        i: int
            The index of the line in the text slice
        line_text: str
            A line of text from a JDFTx out file
        _iter_flag:  str
            The flag that indicates the start of a log message for a JDFTx
            optimization step

        Returns
        -------
        is_line: bool
            True if the line_text is the start of a log message for a JDFTx
            optimization step
        """
        return _iter_flag in line_text

    def read_iter_line(self, line_text: str) -> None:
        """Set class variables iter, E, grad_K, alpha, linmin, t_s.

        Parse the lines of text corresponding to the electronic minimization
        data of a JDFTx out file.

        Parameters
        ----------
        line_text: str
            A line of text from a JDFTx out file containing the electronic
            minimization data
        """
        nstep_float = get_colon_var_t1(line_text, "Iter: ")
        if isinstance(nstep_float, float):
            self.nstep = int(nstep_float)
        elif nstep_float is None:
            raise ValueError("Could not find nstep in line_text")
        self.e = get_colon_var_t1(line_text, f"{self.etype}: ") * Ha_to_eV
        self.grad_k = get_colon_var_t1(line_text, "|grad|_K: ")
        self.alpha = get_colon_var_t1(line_text, "alpha: ")
        self.linmin = get_colon_var_t1(line_text, "linmin: ")
        self.t_s = get_colon_var_t1(line_text, "t[s]: ")

    def is_fillings_line(self, i: int, line_text: str) -> bool:
        """Return True if fillings line.

        Return True if the line_text is the start of a log message for a
        JDFTx optimization step.

        Parameters
        ----------
        i: int
            The index of the line in the text slice
        line_text: str
            A line of text from a JDFTx out file

        Returns
        -------
        is_line: bool
            True if the line_text is the start of a log message for a JDFTx
            optimization step
        """
        return "FillingsUpdate" in line_text

    def read_fillings_line(self, fillings_line: str) -> None:
        """Set class variables mu, nelectrons, magneticmoment.

        Parse the lines of text corresponding to the electronic minimization
        data of a JDFTx out file.

        Parameters
        ----------
        fillings_line: str
            A line of text from a JDFTx out file containing the electronic
            minimization data
        """
        if "FillingsUpdate:" in fillings_line:
            self.set_mu(fillings_line)
            self.set_nelectrons(fillings_line)
            if "magneticMoment" in fillings_line:
                self.set_magdata(fillings_line)
        else:
            raise ValueError("FillingsUpdate string not found")

    def is_subspaceadjust_line(self, i: int, line_text: str) -> bool:
        """Return True if subspace adjust line.

        Return True if the line_text is the start of a log message for a
        JDFTx optimization step.

        Parameters
        ----------
        i: int
            The index of the line in the text slice
        line_text: str
            A line of text from a JDFTx out file

        Returns
        -------
        is_line: bool
            True if the line_text is the start of a log message for a JDFTx
            optimization step
        """
        return "SubspaceRotationAdjust" in line_text

    def read_subspaceadjust_line(self, line_text: str) -> None:
        """Set class variable subspaceRotationAdjust.

        Parse the lines of text corresponding to the electronic minimization
        data of a JDFTx out file.

        Parameters
        ----------
        line_text: str
            A line of text from a JDFTx out file containing the electronic
            minimization data
        """
        self.subspacerotationadjust = get_colon_var_t1(line_text, "SubspaceRotationAdjust: set factor to")

    def set_magdata(self, fillings_line: str) -> None:
        """Set class variables abs_magneticMoment, tot_magneticMoment.

        Parse the lines of text corresponding to the electronic minimization
        data of a JDFTx out file.

        Parameters
        ----------
        fillings_line: str
            A line of text from a JDFTx out file containing the electronic
            minimization data
        """
        _fillings_line = fillings_line.split("magneticMoment: [ ")[1].split(" ]")[0].strip()
        self.abs_magneticmoment = get_colon_var_t1(_fillings_line, "Abs: ")
        self.tot_magneticmoment = get_colon_var_t1(_fillings_line, "Tot: ")

    def set_mu(self, fillings_line: str) -> None:
        """Set mu class variable.

        Parse the lines of text corresponding to the electronic minimization
        data of a JDFTx out file.

        Parameters
        ----------
        fillings_line: str
            A line of text from a JDFTx out file containing the electronic
            minimization data
        """
        self.mu = get_colon_var_t1(fillings_line, "mu: ") * Ha_to_eV

    def set_nelectrons(self, fillings_line: str) -> None:
        """Set nelectrons class variable.

        Parse the lines of text corresponding to the electronic minimization
        data of a JDFTx out file.

        Parameters
        ----------
        fillings_line: str
            A line of text from a JDFTx out file containing the electronic
            minimization data
        """
        self.nelectrons = get_colon_var_t1(fillings_line, "nElectrons: ")


class JElSteps:
    """Class object for series of SCF steps.

    Class object for collecting and storing a series of SCF steps done between
    geometric optimization steps.
    """

    iter_type: str | None = None
    etype: str | None = None
    iter_flag: str | None = None
    converged: bool = False
    converged_reason: str | None = None
    slices: list[JElStep] = field(default_factory=list)
    _getatr_ignore: ClassVar[list[str]] = [
        "e",
        "t_s",
        "mu",
        "nelectrons",
        "subspacerotationadjust",
    ]

    @property
    def nstep(self) -> int:
        """Return nstep.

        Return the nstep attribute of the last JElStep object in the slices.

        Returns
        -------
        nstep: int
            The nstep attribute of the last JElStep object in the slices
        """
        if len(self.slices):
            if self.slices[-1].nstep is not None:
                return self.slices[-1].nstep
            warnings.warn("No nstep attribute in JElStep object. Returning number of JElStep objects.", stacklevel=2)
            return len(self.slices) - 1
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def e(self) -> float:
        """Return total electronic energy.

        Return the e attribute of the last JElStep object in the slices.

        Returns
        -------
        e: float
            The e attribute of the last JElStep object in the slices
        """
        if len(self.slices):
            if self.slices[-1].e is not None:
                return self.slices[-1].e
            raise AttributeError("No E attribute in final JElStep object.")
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def grad_k(self) -> float | None:
        """Return most recent grad_k.

        Return the grad_k attribute of the last JElStep object in the slices.

        Returns
        -------
        grad_k: float
            The grad_k attribute of the last JElStep object in the slices
        """
        if len(self.slices):
            return self.slices[-1].grad_k
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def alpha(self) -> float | None:
        """Return most recent alpha.

        Return the alpha attribute of the last JElStep object in the slices.

        Returns
        -------
        alpha: float
            The alpha attribute of the last JElStep object in the slices
        """
        if len(self.slices):
            return self.slices[-1].alpha
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def linmin(self) -> float | None:
        """Return most recent linmin.

        Return the linmin attribute of the last JElStep object in the slices.

        Returns
        -------
        linmin: float
            The linmin attribute of the last JElStep object in the slices
        """
        if len(self.slices):
            return self.slices[-1].linmin
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def t_s(self) -> float:
        """Return most recent t_s.

        Return the t_s attribute of the last JElStep object in the slices.

        Returns
        -------
        t_s: float
            The t_s attribute of the last JElStep object in the slices
        """
        if len(self.slices):
            if self.slices[-1].t_s is not None:
                return self.slices[-1].t_s
            raise AttributeError("No t_s attribute in final JElStep object.")
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def mu(self) -> float:
        """Return most recent mu.

        Return the mu attribute of the last JElStep object in the slices.

        Returns
        -------
        mu: float
            The mu attribute of the last JElStep object in the slices
        """
        if len(self.slices):
            if self.slices[-1].mu is not None:
                return self.slices[-1].mu
            raise AttributeError("No mu attribute in final JElStep object.")
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def nelectrons(self) -> float:
        """Return most recent nelectrons.

        Return the nelectrons attribute of the last JElStep object in the slices.

        Returns
        -------
        nelectrons: float
            The nelectrons attribute of the last JElStep object in the slices
        """
        if len(self.slices):
            if self.slices[-1].nelectrons is not None:
                return self.slices[-1].nelectrons
            raise AttributeError("No nelectrons attribute in final JElStep object.")
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def abs_magneticmoment(self) -> float | None:
        """Return most recent abs_magneticmoment.

        Return the abs_magneticmoment attribute of the last JElStep object in the slices.

        Returns
        -------
        abs_magneticmoment: float
            The abs_magneticmoment attribute of the last JElStep object in the slices
        """
        if len(self.slices):
            return self.slices[-1].abs_magneticmoment
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def tot_magneticmoment(self) -> float | None:
        """Return most recent tot_magneticmoment.

        Return the tot_magneticmoment attribute of the last JElStep object in the slices.

        Returns
        -------
        tot_magneticmoment: float
            The tot_magneticmoment attribute of the last JElStep object in the slices
        """
        if len(self.slices):
            return self.slices[-1].tot_magneticmoment
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def subspacerotationadjust(self) -> float:
        """Return most recent subspacerotationadjust.

        Return the subspacerotationadjust attribute of the last JElStep object in the slices.

        Returns
        -------
        subspacerotationadjust: float
            The subspacerotationadjust attribute of the last JElStep object in the slices
        """
        if len(self.slices):
            if self.slices[-1].subspacerotationadjust is not None:
                return self.slices[-1].subspacerotationadjust
            raise AttributeError("No subspacerotationadjust attribute in final JElStep object.")
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @classmethod
    def from_text_slice(cls, text_slice: list[str], iter_type: str = "ElecMinimize", etype: str = "F") -> JElSteps:
        """Return JElSteps object.

        Create a JElSteps object from a slice of an out file's text
        corresponding to a series of SCF steps.

        Parameters
        ----------
        text_slice : list[str]
            A slice of text from a JDFTx out file corresponding to a series of
            SCF steps
        iter_type: str
            The type of electronic minimization step
        etype: str
            The type of energy component
        """
        line_collections, lines_collect = gather_JElSteps_line_collections(iter_type, text_slice)
        # instance = cls.from_lines_collect(line_collections[-1], iter_type, etype)
        instance = cls()
        instance.iter_flag = f"{iter_type}: Iter:"
        instance.iter_type = iter_type
        instance.etype = etype
        instance.slices = []
        for _lines_collect in line_collections:
            instance.slices.append(JElStep.from_lines_collect(_lines_collect, iter_type, etype))
        if len(lines_collect):
            instance.parse_ending_lines(lines_collect)
            lines_collect = []
        return instance

    def parse_ending_lines(self, ending_lines: list[str]) -> None:
        """Parse ending lines.

        Parses the ending lines of text from a JDFTx out file corresponding to
        a series of SCF steps.

        Parameters
        ----------
        ending_lines: list[str]
            The ending lines of text from a JDFTx out file corresponding to a
            series of SCF steps
        """
        for i, line in enumerate(ending_lines):
            if self.is_converged_line(i, line):
                self.read_converged_line(line)

    def is_converged_line(self, i: int, line_text: str) -> bool:
        """Return True if converged line.

        Return True if the line_text is the start of a log message about
        convergence for a JDFTx optimization step

        Parameters
        ----------
        i: int
            The index of the line in the text slice
        line_text: str
            A line of text from a JDFTx out file

        Returns
        -------
        is_line: bool
            True if the line_text is the start of a log message about
            convergence for a JDFTx optimization step
        """
        return f"{self.iter_type}: Converged" in line_text

    def read_converged_line(self, line_text: str) -> None:
        """Set class variables converged and converged_reason.

        Read the convergence message from a JDFTx optimization step

        Parameters
        ----------
        line_text: str
            A line of text from a JDFTx out file containing a message about
            convergence for a JDFTx optimization step
        """
        self.converged = True
        self.converged_reason = line_text.split("(")[1].split(")")[0].strip()

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
            if name not in self._getatr_ignore:
                if not hasattr(self.slices[-1], name):
                    raise AttributeError(f"{self.__class__.__name__} not found: {name}")
                return getattr(self.slices[-1], name)
            raise AttributeError(f"Property {name} inaccessible due to empty slices class field")
        raise AttributeError(f"Property {name} inaccessible due to empty slices class field")

    def __getitem__(self, key: int | str) -> JElStep | Any:
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

    def getitem_int(self, key: int) -> JElStep:
        """Return JElStep object.

        Return the JElStep object at the key index.

        Parameters
        ----------
        key: int
            The index of the JElStep object

        Returns
        -------
        JElStep: JElStep
            The JElStep object at the key index
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
        """Return length of JElSteps object.

        Returns the number of SCF steps in the JElSteps object.

        Returns
        -------
        length: int
            The number of SCF steps in the JElSteps object
        """
        return len(self.slices)
