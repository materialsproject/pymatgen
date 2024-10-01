"""Module for JEiters class object.

This module contains the JEiters class object for parsing a series of SCF steps.

"""

from __future__ import annotations

import warnings
from dataclasses import field
from typing import Any, ClassVar

from pymatgen.io.jdftx.jeiter import JEiter
from pymatgen.io.jdftx.utils import gather_jeiters_line_collections

__author__ = "Ben Rich"


# def gather_jeiters_line_collections(iter_type: str, text_slice: list[str]) -> tuple[list[list[str]], list[str]]:
#     """Gather line collections for JEiters initialization.

#     Gathers list of line lists where each line list initializes a JEiter object,
#     and the remaining lines that do not initialize a JEiter object are used
#     for initialization unique to the JEiters object.

#     Parameters
#     ----------
#     iter_type: str
#         The type of electronic minimization step
#     text_slice: list[str]
#         A slice of text from a JDFTx out file corresponding to a series of
#         SCF steps

#     Returns
#     -------
#     line_collections: list[list[str]]
#         A list of lists of lines of text from a JDFTx out file corresponding to
#         a single SCF step
#     lines_collect: list[str]
#         A list of lines of text from a JDFTx out file corresponding to a single
#         SCF step

#     """
#     lines_collect = []
#     line_collections = []
#     _iter_flag = f"{iter_type}: Iter:"
#     for line_text in text_slice:
#         if len(line_text.strip()):
#             lines_collect.append(line_text)
#             if _iter_flag in line_text:
#                 line_collections.append(lines_collect)
#                 lines_collect = []
#         else:
#             break
#     return line_collections, lines_collect


class JEiters:
    """Class object for series of SCF steps.

    Class object for collecting and storing a series of SCF steps done between
    geometric optimization steps.
    """

    iter_type: str | None = None
    etype: str | None = None
    iter_flag: str | None = None
    converged: bool = False
    converged_reason: str | None = None
    slices: list[JEiter] = field(default_factory=list)
    _getatr_ignore: ClassVar[list[str]] = [
        "e",
        "t_s",
        "mu",
        "nelectrons",
        "subspacerotationadjust",
    ]

    @property
    def niter(self) -> int:
        """Return niter.

        Return the niter attribute of the last JEiter object in the slices.

        Returns
        -------
        niter: int
            The niter attribute of the last JEiter object in the slices
        """
        if len(self.slices):
            if self.slices[-1].niter is not None:
                return self.slices[-1].niter
            warnings.warn("No niter attribute in JEiter object. Returning number of JEiter objects.", stacklevel=2)
            return len(self.slices) - 1
        raise AttributeError("No JEiter objects in JEiters object slices class variable.")

    @property
    def e(self) -> float:
        """Return total electronic energy.

        Return the e attribute of the last JEiter object in the slices.

        Returns
        -------
        e: float
            The e attribute of the last JEiter object in the slices
        """
        if len(self.slices):
            if self.slices[-1].e is not None:
                return self.slices[-1].e
            raise AttributeError("No E attribute in final JEiter object.")
        raise AttributeError("No JEiter objects in JEiters object slices class variable.")

    @property
    def grad_k(self) -> float | None:
        """Return most recent grad_k.

        Return the grad_k attribute of the last JEiter object in the slices.

        Returns
        -------
        grad_k: float
            The grad_k attribute of the last JEiter object in the slices
        """
        if len(self.slices):
            return self.slices[-1].grad_k
        raise AttributeError("No JEiter objects in JEiters object slices class variable.")

    @property
    def alpha(self) -> float | None:
        """Return most recent alpha.

        Return the alpha attribute of the last JEiter object in the slices.

        Returns
        -------
        alpha: float
            The alpha attribute of the last JEiter object in the slices
        """
        if len(self.slices):
            return self.slices[-1].alpha
        raise AttributeError("No JEiter objects in JEiters object slices class variable.")

    @property
    def linmin(self) -> float | None:
        """Return most recent linmin.

        Return the linmin attribute of the last JEiter object in the slices.

        Returns
        -------
        linmin: float
            The linmin attribute of the last JEiter object in the slices
        """
        if len(self.slices):
            return self.slices[-1].linmin
        raise AttributeError("No JEiter objects in JEiters object slices class variable.")

    @property
    def t_s(self) -> float:
        """Return most recent t_s.

        Return the t_s attribute of the last JEiter object in the slices.

        Returns
        -------
        t_s: float
            The t_s attribute of the last JEiter object in the slices
        """
        if len(self.slices):
            if self.slices[-1].t_s is not None:
                return self.slices[-1].t_s
            raise AttributeError("No t_s attribute in final JEiter object.")
        raise AttributeError("No JEiter objects in JEiters object slices class variable.")

    @property
    def mu(self) -> float:
        """Return most recent mu.

        Return the mu attribute of the last JEiter object in the slices.

        Returns
        -------
        mu: float
            The mu attribute of the last JEiter object in the slices
        """
        if len(self.slices):
            if self.slices[-1].mu is not None:
                return self.slices[-1].mu
            raise AttributeError("No mu attribute in final JEiter object.")
        raise AttributeError("No JEiter objects in JEiters object slices class variable.")

    @property
    def nelectrons(self) -> float:
        """Return most recent nelectrons.

        Return the nelectrons attribute of the last JEiter object in the slices.

        Returns
        -------
        nelectrons: float
            The nelectrons attribute of the last JEiter object in the slices
        """
        if len(self.slices):
            if self.slices[-1].nelectrons is not None:
                return self.slices[-1].nelectrons
            raise AttributeError("No nelectrons attribute in final JEiter object.")
        raise AttributeError("No JEiter objects in JEiters object slices class variable.")

    @property
    def abs_magneticmoment(self) -> float | None:
        """Return most recent abs_magneticmoment.

        Return the abs_magneticmoment attribute of the last JEiter object in the slices.

        Returns
        -------
        abs_magneticmoment: float
            The abs_magneticmoment attribute of the last JEiter object in the slices
        """
        if len(self.slices):
            return self.slices[-1].abs_magneticmoment
        raise AttributeError("No JEiter objects in JEiters object slices class variable.")

    @property
    def tot_magneticmoment(self) -> float | None:
        """Return most recent tot_magneticmoment.

        Return the tot_magneticmoment attribute of the last JEiter object in the slices.

        Returns
        -------
        tot_magneticmoment: float
            The tot_magneticmoment attribute of the last JEiter object in the slices
        """
        if len(self.slices):
            return self.slices[-1].tot_magneticmoment
        raise AttributeError("No JEiter objects in JEiters object slices class variable.")

    @property
    def subspacerotationadjust(self) -> float:
        """Return most recent subspacerotationadjust.

        Return the subspacerotationadjust attribute of the last JEiter object in the slices.

        Returns
        -------
        subspacerotationadjust: float
            The subspacerotationadjust attribute of the last JEiter object in the slices
        """
        if len(self.slices):
            if self.slices[-1].subspacerotationadjust is not None:
                return self.slices[-1].subspacerotationadjust
            raise AttributeError("No subspacerotationadjust attribute in final JEiter object.")
        raise AttributeError("No JEiter objects in JEiters object slices class variable.")

    @classmethod
    def from_text_slice(cls, text_slice: list[str], iter_type: str = "ElecMinimize", etype: str = "F") -> JEiters:
        """Return JEiters object.

        Create a JEiters object from a slice of an out file's text
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
        line_collections, lines_collect = gather_jeiters_line_collections(iter_type, text_slice)
        # instance = cls.from_lines_collect(line_collections[-1], iter_type, etype)
        instance = cls()
        instance.iter_flag = f"{iter_type}: Iter:"
        instance.iter_type = iter_type
        instance.etype = etype
        instance.slices = []
        for _lines_collect in line_collections:
            instance.slices.append(JEiter.from_lines_collect(_lines_collect, iter_type, etype))
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

    def __getitem__(self, key: int | str) -> JEiter | Any:
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

    def getitem_int(self, key: int) -> JEiter:
        """Return JEiter object.

        Return the JEiter object at the key index.

        Parameters
        ----------
        key: int
            The index of the JEiter object

        Returns
        -------
        jeiter: JEiter
            The JEiter object at the key index
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
        """Return length of JEiters object.

        Returns the number of SCF steps in the JEiters object.

        Returns
        -------
        length: int
            The number of SCF steps in the JEiters object
        """
        return len(self.slices)
