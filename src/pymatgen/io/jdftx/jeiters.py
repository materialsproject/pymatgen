"""Module for JEiters class object.

This module contains the JEiters class object for parsing a series of SCF steps.

"""

from __future__ import annotations

from dataclasses import field
from typing import Any

from atomate2.jdftx.io.jeiter import JEiter


def gather_line_collections(
    iter_type: str, text_slice: list[str]
) -> tuple[list[list[str]], list[str]]:
    """Gather line collections for JEiters initialization.

    Gathers list of line lists where each line list initializes a JEiter object,
    and the remaining lines that do not initialize a JEiter object are used
    for initialization unique to the JEiters object.

    Parameters
    ----------
    iter_type: str
        The type of electronic minimization step
    text_slice: list[str]
        A slice of text from a JDFTx out file corresponding to a series of
        SCF steps

    Returns
    -------
    line_collections: list[list[str]]
        A list of lists of lines of text from a JDFTx out file corresponding to
        a single SCF step
    lines_collect: list[str]
        A list of lines of text from a JDFTx out file corresponding to a single
        SCF step

    """
    lines_collect = []
    line_collections = []
    _iter_flag = f"{iter_type}: Iter:"
    for line_text in text_slice:
        if len(line_text.strip()):
            lines_collect.append(line_text)
            if _iter_flag in line_text:
                line_collections.append(lines_collect)
                lines_collect = []
        else:
            break
    return line_collections, lines_collect


class JEiters:
    """Class object for series of SCF steps.

    Class object for collecting and storing a series of SCF steps done between
    geometric optimization steps.
    """

    iter_type: str = None
    etype: str = None
    iter_flag: str = None
    converged: bool = False
    converged_reason: str = None
    slices: list[JEiter] = field(default_factory=list)

    @classmethod
    def from_text_slice(
        cls, text_slice: list[str], iter_type: str = "ElecMinimize", etype: str = "F"
    ) -> JEiters:
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
        line_collections, lines_collect = gather_line_collections(iter_type, text_slice)
        # instance = cls.from_lines_collect(line_collections[-1], iter_type, etype)
        instance = cls()
        instance.iter_flag = f"{iter_type}: Iter:"
        instance.iter_type = iter_type
        instance.etype = etype
        instance.slices = []
        for _lines_collect in line_collections:
            instance.slices.append(
                JEiter.from_lines_collect(_lines_collect, iter_type, etype)
            )
        if len(lines_collect):
            instance.parse_ending_lines(lines_collect)
            lines_collect = []
        return instance

    def parse_text_slice(self, text_slice: list[str]) -> None:
        """Parse text slice.

        Parse a slice of text from a JDFTx out file corresponding to a series
        of SCF steps.

        Parameters
        ----------
        text_slice: list[str]
            A slice of text from a JDFTx out file corresponding to a series of
            SCF steps
        """
        lines_collect = []
        _iter_flag = f"{self.iter_type}: Iter:"
        for line_text in text_slice:
            if len(line_text.strip()):
                lines_collect.append(line_text)
                if _iter_flag in line_text:
                    self.slices.append(
                        JEiter.from_lines_collect(
                            lines_collect, self.iter_type, self.etype
                        )
                    )
                    lines_collect = []
            else:
                break
        if len(lines_collect):
            self.parse_ending_lines(lines_collect)
            lines_collect = []

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

    def __getatr__(self, name: str) -> Any:
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
        if not hasattr(self, name):
            if not hasattr(self.slices[-1], name):
                raise AttributeError(f"{self.__class__.__name__} not found: {name}")
            return getattr(self.slices[-1], name)
        return getattr(self, name)

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
