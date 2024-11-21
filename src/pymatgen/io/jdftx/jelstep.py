"""Module for parsing single SCF step from JDFTx.

This module contains the JElStep class for parsing single SCF step from a JDFTx out file.

@mkhorton - this file is ready to review.
"""

from __future__ import annotations

import inspect
import pprint
import warnings
from dataclasses import dataclass, field
from typing import Any, ClassVar

from pymatgen.core.units import Ha_to_eV
from pymatgen.io.jdftx._output_utils import get_colon_var_t1

__author__ = "Ben Rich"


@dataclass
class JElStep:
    """Electronic minimization data for a single SCF step.

    Class object for storing logged electronic minimization data for a single
    SCF step.

    Attributes:
        opt_type (str | None): The type of electronic minimization step
            (almost always ElecMinimize).
        etype (str | None): The type of energy component (G, F, or Etot).
        nstep (int | None): The SCF step number.
        e (float | None): The total electronic energy in eV.
        grad_k (float | None): The gradient of the Kohn-Sham energy (along the
            line minimization direction).
        alpha (float | None): The step length.
        linmin (float | None): Normalized line minimization direction / energy
            gradient projection (-1 for perfectly opposite, 1 for perfectly aligned).
        t_s (float | None): Time in seconds for the SCF step.
        mu (float | None): The chemical potential in eV.
        nelectrons (float | None): The number of electrons.
        abs_magneticmoment (float | None): The absolute magnetic moment.
        tot_magneticmoment (float | None): The total magnetic moment.
        subspacerotationadjust (float | None): The subspace rotation adjustment factor.
    """

    opt_type: str | None = None
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
    def _from_lines_collect(cls, lines_collect: list[str], opt_type: str, etype: str) -> JElStep:
        """Return JElStep object.

        Create a JElStep object from a list of lines of text from a JDFTx out
        file corresponding to a single SCF step.

        Args:
        lines_collect (list[str]): A list of lines of text from a JDFTx out file corresponding to a single SCF step.
        opt_type (str): The type of electronic minimization step.
        etype (str): The type of energy component.

        Returns:
            JElStep: The created JElStep object.
        """
        instance = cls()
        instance.opt_type = opt_type
        instance.etype = etype
        _iter_flag = f"{opt_type}: Iter: "
        for i, line_text in enumerate(lines_collect):
            if instance._is_iter_line(i, line_text, _iter_flag):
                instance._read_iter_line(line_text)
            elif instance._is_fillings_line(i, line_text):
                instance._read_fillings_line(line_text)
            elif instance._is_subspaceadjust_line(i, line_text):
                instance._read_subspaceadjust_line(line_text)
        return instance

    def _is_iter_line(self, i: int, line_text: str, _iter_flag: str) -> bool:
        """Return True if opt iter line.

        Return True if the line_text is the start of a log message for a
        JDFTx optimization step.

        Args:
            i (int): The index of the line in the text slice.
            line_text (str): A line of text from a JDFTx out file.
            _iter_flag (str): The flag that indicates the start of a log message for a JDFTx optimization step.

        Returns:
            bool: True if the line_text is the start of a log message for a JDFTx optimization step.
        """
        return _iter_flag in line_text

    def _read_iter_line(self, line_text: str) -> None:
        """Set class variables iter, E, grad_K, alpha, linmin, t_s.

        Parse the lines of text corresponding to the electronic minimization
        data of a JDFTx out file.

        Args:
            line_text (str): A line of text from a JDFTx out file containing the electronic minimization data.
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

    def _is_fillings_line(self, i: int, line_text: str) -> bool:
        """Return True if fillings line.

        Return True if the line_text is the start of a log message for a
        JDFTx optimization step.

        Args:
            i (int): The index of the line in the text slice.
            line_text (str): A line of text from a JDFTx out file.

        Returns:
            bool: True if the line_text is the start of a log message for a JDFTx optimization step.
        """
        return "FillingsUpdate" in line_text

    def _read_fillings_line(self, fillings_line: str) -> None:
        """Set class variables mu, nelectrons, magneticmoment.

        Parse the lines of text corresponding to the electronic minimization
        data of a JDFTx out file.

        Args:
            fillings_line (str): A line of text from a JDFTx out file containing the electronic
            minimization data.
        """
        if "FillingsUpdate:" in fillings_line:
            self._set_mu(fillings_line)
            self._set_nelectrons(fillings_line)
            if "magneticMoment" in fillings_line:
                self._set_magdata(fillings_line)
        else:
            raise ValueError("FillingsUpdate string not found")

    def _is_subspaceadjust_line(self, i: int, line_text: str) -> bool:
        """Return True if the line_text is the start of a log message for a JDFTx optimization step.

        Args:
            i (int): The index of the line in the text slice.
            line_text (str): A line of text from a JDFTx out file.

        Returns:
            bool: True if the line_text is the start of a log message for a JDFTx optimization step.
        """
        return "SubspaceRotationAdjust" in line_text

    def _read_subspaceadjust_line(self, line_text: str) -> None:
        """Set class variable subspaceRotationAdjust.

        Parse the lines of text corresponding to the electronic minimization
        data of a JDFTx out file.

        Args:
            line_text (str): A line of text from a JDFTx out file containing the electronic
            minimization data.
        """
        self.subspacerotationadjust = get_colon_var_t1(line_text, "SubspaceRotationAdjust: set factor to")

    def _set_magdata(self, fillings_line: str) -> None:
        """Set class variables abs_magneticMoment, tot_magneticMoment.

        Parse the lines of text corresponding to the electronic minimization
        data of a JDFTx out file.

        Args:
            fillings_line (str): A line of text from a JDFTx out file containing the electronic
            minimization data.
        """
        _fillings_line = fillings_line.split("magneticMoment: [ ")[1].split(" ]")[0].strip()
        self.abs_magneticmoment = get_colon_var_t1(_fillings_line, "Abs: ")
        self.tot_magneticmoment = get_colon_var_t1(_fillings_line, "Tot: ")

    def _set_mu(self, fillings_line: str) -> None:
        """Set mu class variable.

        Parse the lines of text corresponding to the electronic minimization
        data of a JDFTx out file.

        Args:
            fillings_line (str): A line of text from a JDFTx out file containing the electronic
            minimization data.
        """
        self.mu = get_colon_var_t1(fillings_line, "mu: ") * Ha_to_eV

    def _set_nelectrons(self, fillings_line: str) -> None:
        """Set nelectrons class variable.

        Parse the lines of text corresponding to the electronic minimization
        data of a JDFTx out file.

        Args:
            fillings_line(str): A line of text from a JDFTx out file containing the electronic minimization data
        """
        self.nelectrons = get_colon_var_t1(fillings_line, "nElectrons: ")

    def __str__(self) -> str:
        """
        Return string representation of JElStep object.

        Returns:
            str: String representation of JElStep object.
        """
        return pprint.pformat(self)


@dataclass
class JElSteps:
    """Class object for series of SCF steps.

    Class object for collecting and storing a series of SCF steps done between
    geometric optimization steps.

    Attributes:
        opt_type (str | None): The type of electronic minimization step.
        etype (str | None): The type of energy component.
        iter_flag (str | None): The flag that indicates the start of a log message for a JDFTx optimization step.
        converged (bool): True if the SCF steps converged.
        converged_reason (str | None): The reason for convergence.
        slices (list[JElStep]): A list of JElStep objects.
    """

    opt_type: str | None = None
    etype: str | None = None
    iter_flag: str | None = None
    converged: bool = False
    converged_reason: str | None = None
    slices: list[JElStep] = field(default_factory=list)
    # List of attributes to ignore when getting attributes from the most recent slice specified by _getatr_ignore
    _getatr_ignore: ClassVar[list[str]] = [
        "e",
        "t_s",
        "mu",
        "nelectrons",
        "subspacerotationadjust",
    ]

    @property
    def nstep(self) -> int | None:
        """Return the nstep attribute of the last JElStep object in the slices.

        The nstep attribute signifies the SCF step number.

        Returns:
            int: The nstep attribute of the last JElStep object in the slices, or the number of JElStep objects
            if nstep is None.

        Raises:
            AttributeError: If there are no JElStep objects in the slices.
        """
        if len(self.slices):
            if self.slices[-1].nstep is not None:
                return self.slices[-1].nstep
            warnings.warn("No nstep attribute in JElStep object. Returning number of JElStep objects.", stacklevel=2)
            return len(self.slices) - 1
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def e(self) -> float | None:
        """Return total electronic energy.

        Return the e attribute of the last JElStep object in the slices, where e
        signifies the total electronic energy in eV.

        Returns:
            float: The e attribute of the last JElStep object in the slices.
        """
        if len(self.slices):
            return self.slices[-1].e
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def grad_k(self) -> float | None:
        """Return most recent grad_k.

        Return the grad_k attribute of the last JElStep object in the slices, where
        grad_k signifies the gradient of the Kohn-Sham energy (along line minimization direction).

        Returns:
            float: The grad_k attribute of the last JElStep object in the slices.
        """
        if len(self.slices):
            return self.slices[-1].grad_k
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def alpha(self) -> float | None:
        """Return most recent alpha.

        Return the alpha attribute of the last JElStep object in the slices, where
        alpha signifies the step length in the electronic minimization.

        Returns:
            float: The alpha attribute of the last JElStep object in the slices.
        """
        if len(self.slices):
            return self.slices[-1].alpha
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def linmin(self) -> float | None:
        """Return most recent linmin.

        Return the linmin attribute of the last JElStep object in the slices, where
        linmin signifies the normalized line minimization direction / energy gradient projection.

        Returns:
            float: The linmin attribute of the last JElStep object in the slices.
        """
        if len(self.slices):
            return self.slices[-1].linmin
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def t_s(self) -> float | None:
        """Return most recent t_s.

        Return the t_s attribute of the last JElStep object in the slices, where
        t_s signifies the time in seconds for the SCF step.

        Returns:
            float: The t_s attribute of the last JElStep object in the slices.
        """
        if len(self.slices):
            return self.slices[-1].t_s
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def mu(self) -> float | None:
        """Return most recent mu.

        Return the mu attribute of the last JElStep object in the slices, where
        mu signifies the chemical potential (Fermi level) in eV.

        Returns:
            float: The mu attribute of the last JElStep object in the slices.
        """
        if len(self.slices):
            return self.slices[-1].mu
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def nelectrons(self) -> float | None:
        """Return most recent nelectrons.

        Return the nelectrons attribute of the last JElStep object in the slices, where
        nelectrons signifies the total number of electrons being evaluated in the SCF step.

        Returns:
            float: The nelectrons attribute of the last JElStep object in the slices.
        """
        if len(self.slices):
            return self.slices[-1].nelectrons
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def abs_magneticmoment(self) -> float | None:
        """Return most recent abs_magneticmoment.

        Return the abs_magneticmoment attribute of the last JElStep object in the slices, where
        abs_magneticmoment signifies the absolute magnetic moment of the electron density.

        Returns:
            float: The abs_magneticmoment attribute of the last JElStep object in the slices.
        """
        if len(self.slices):
            return self.slices[-1].abs_magneticmoment
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def tot_magneticmoment(self) -> float | None:
        """
        Return most recent tot_magneticmoment.

        Return the tot_magneticmoment attribute of the last JElStep object in the slices, where
        tot_magneticmoment signifies the total magnetic moment of the electron density.

        Returns:
            float: The tot_magneticmoment attribute of the last JElStep object in the slices.
        """
        if len(self.slices):
            return self.slices[-1].tot_magneticmoment
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @property
    def subspacerotationadjust(self) -> float | None:
        """Return most recent subspacerotationadjust.

        Return the subspacerotationadjust attribute of the last JElStep object in the slices, where
        subspacerotationadjust signifies the amount by which the subspace was rotated in the SCF step.

        Returns:
            float: The subspacerotationadjust attribute of the last JElStep object in the slices.
        """
        if len(self.slices):
            return self.slices[-1].subspacerotationadjust
        raise AttributeError("No JElStep objects in JElSteps object slices class variable.")

    @classmethod
    def _from_text_slice(cls, text_slice: list[str], opt_type: str = "ElecMinimize", etype: str = "F") -> JElSteps:
        """Return JElSteps object.

        Create a JElSteps object from a slice of an out file's text
        corresponding to a series of SCF steps.

        Args:
            text_slice (list[str]): A slice of text from a JDFTx out file corresponding to a series of SCF steps.
            opt_type (str): The type of electronic minimization step.
            etype (str): The type of energy component.
        """
        line_collections, lines_collect = _gather_JElSteps_line_collections(opt_type, text_slice)
        instance = cls()
        instance.iter_flag = f"{opt_type}: Iter:"
        instance.opt_type = opt_type
        instance.etype = etype
        instance.slices = []
        for _lines_collect in line_collections:
            instance.slices.append(JElStep._from_lines_collect(_lines_collect, opt_type, etype))
        if len(lines_collect):
            instance._parse_ending_lines(lines_collect)
            lines_collect = []
        return instance

    def _parse_ending_lines(self, ending_lines: list[str]) -> None:
        """Parse ending lines.

        Parses the ending lines of text from a JDFTx out file corresponding to
        a series of SCF steps.

        Args:
            ending_lines (list[str]): The ending lines of text from a JDFTx out file corresponding to a
            series of SCF steps.
        """
        for i, line in enumerate(ending_lines):
            if self._is_converged_line(i, line):
                self._read_converged_line(line)

    def _is_converged_line(self, i: int, line_text: str) -> bool:
        """Return True if converged line.

        Return True if the line_text is the start of a log message about
        convergence for a JDFTx optimization step.

        Args:
            i (int): The index of the line in the text slice.
            line_text (str): A line of text from a JDFTx out file.

        Returns:
            bool: True if the line_text is the start of a log message about
            convergence for a JDFTx optimization step.
        """
        return f"{self.opt_type}: Converged" in line_text

    def _read_converged_line(self, line_text: str) -> None:
        """Set class variables converged and converged_reason.

        Args:
            line_text (str): A line of text from a JDFTx out file containing a message about
            convergence for a JDFTx optimization step.
        """
        self.converged = True
        self.converged_reason = line_text.split("(")[1].split(")")[0].strip()

    # This method is likely never going to be called as all (currently existing)
    # attributes of the most recent slice are explicitly defined as a class
    # property. However, it is included to reduce the likelihood of errors
    # upon future changes to downstream code.
    def __getattr__(self, name: str) -> Any:
        """Return attribute value.

        Args:
            name (str): The name of the attribute.

        Returns:
            Any: The value of the attribute.

        Raises:
            AttributeError: If the attribute is not found.
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
            val = self._getitem_int(key)
        if type(key) is str:
            val = self._getitem_str(key)
        return val

    def _getitem_int(self, key: int) -> JElStep:
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

    def _getitem_str(self, key: str) -> Any:
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

    def __str__(self) -> str:
        """Return string representation of JElSteps object.

        Returns
        -------
        str: str
            String representation of JElSteps object
        """
        return pprint.pformat(self)


def _gather_JElSteps_line_collections(opt_type: str, text_slice: list[str]) -> tuple[list[list[str]], list[str]]:
    """Gather line collections for JElSteps initialization.

    Gathers list of line lists where each line list initializes a JElStep object,
    and the remaining lines that do not initialize a JElStep object are used
    for initialization unique to the JElSteps object.

    Args:
        opt_type (str): The type of electronic minimization step.
        text_slice (list[str]): A slice of text from a JDFTx out file corresponding to a series of SCF steps.

    Returns:
        tuple: A tuple containing:
            line_collections (list[list[str]]): A list of lists of lines of text from a JDFTx out file
            corresponding to a single SCF step.
            lines_collect (list[str]): A list of lines of text from a JDFTx out file corresponding to a single SCF step.
    """
    lines_collect = []
    line_collections = []
    _iter_flag = f"{opt_type}: Iter:"
    for line_text in text_slice:
        if len(line_text.strip()):
            lines_collect.append(line_text)
            if _iter_flag in line_text:
                line_collections.append(lines_collect)
                lines_collect = []
        else:
            break
    return line_collections, lines_collect
