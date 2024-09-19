"""Helper functions for creation of JOutStructure(s).

This module contains helper functions for the creation of JOutStructure(s) from
the output files of JDFTx calculations.
"""

from __future__ import annotations


def _get_colon_var_t1(linetext: str, lkey: str) -> float | None:
    """Return float val from '...lkey: val...' in linetext.

    Read a float from an elec minimization line assuming value appears as
    "... lkey value ...".

    Parameters
    ----------
    linetext: str
        A line of text from a JDFTx out file
    lkey: str
        A string that appears before the float value in linetext
    """
    colon_var = None
    if lkey in linetext:
        colon_var = float(linetext.split(lkey)[1].strip().split(" ")[0])
    return colon_var


def correct_iter_type(iter_type: str | None) -> str | None:
    """Return recognizable iter_type string.

    Correct the iter_type string to match the JDFTx convention.

    Parameters
    ----------
    iter_type:
        The type of optimization step

    Returns
    -------
    iter_type: str | None
        The corrected type of optimization step
    """
    if iter_type is not None:
        if "lattice" in iter_type.lower():
            iter_type = "LatticeMinimize"
        elif "ionic" in iter_type.lower():
            iter_type = "IonicMinimize"
        else:
            iter_type = None
    return iter_type


def is_strain_start_line(line_text: str) -> bool:
    """Return True if the line_text is the start of strain log message.

    Return True if the line_text is the start of a log message for a JDFTx
    optimization step.

    Parameters
    ----------
    line_text: str
        A line of text from a JDFTx out file

    Returns
    -------
    is_line: bool
        True if the line_text is the start of a log message for a JDFTx
        optimization step
    """
    return "# Strain tensor in" in line_text


def is_lattice_start_line(line_text: str) -> bool:
    """Return True if the line_text is the start of lattice log message.

    Return True if the line_text is the start of a log message for a JDFTx
    optimization step.

    Parameters
    ----------
    line_text: str
        A line of text from a JDFTx out file

    Returns
    -------
    is_line: bool
        True if the line_text is the start of a log message for a JDFTx
        optimization step
    """
    return "# Lattice vectors:" in line_text


def is_forces_start_line(line_text: str) -> bool:
    """Return True if the line_text is the start of forces log message.

    Return True if the line_text is the start of a log message for a JDFTx
    optimization step.

    Parameters
    ----------
    line_text: str
        A line of text from a JDFTx out file

    Returns
    -------
    is_line: bool
        True if the line_text is the start of a log message for a JDFTx
        optimization step
    """
    return "# Forces in" in line_text


def is_stress_start_line(line_text: str) -> bool:
    """Return True if the line_text is the start of stress log message.

    Return True if the line_text is the start of a log message for a JDFTx
    optimization step.

    Parameters
    ----------
    line_text: str
        A line of text from a JDFTx out file

    Returns
    -------
    is_line: bool
        True if the line_text is the start of a log message for a JDFTx
        optimization step
    """
    return "# Stress tensor in" in line_text


def is_posns_start_line(line_text: str) -> bool:
    """Return True if the line_text is the start of posns log message.

    Return True if the line_text is the start of a log message for a JDFTx
    optimization step.

    Parameters
    ----------
    line_text: str
        A line of text from a JDFTx out file containing the positions of atoms

    Returns
    -------
        is_line: bool
            True if the line_text is the start of a log message for a JDFTx
            optimization step
    """
    return "# Ionic positions" in line_text


def is_ecomp_start_line(line_text: str) -> bool:
    """Return True if the line_text is the start of ecomp log message.

    Return True if the line_text is the start of a log message for a JDFTx
    optimization step.

    Parameters
    ----------
    line_text: str
        A line of text from a JDFTx out file

    Returns
    -------
    is_line: bool
        True if the line_text is the start of a log message for a JDFTx
        optimization step
    """
    return "# Energy components" in line_text


def is_lowdin_start_line(line_text: str) -> bool:
    """Return True if the line_text is the start of Lowdin log message.

    Return True if the line_text is the start of a Lowdin population analysis
    in a JDFTx out file.

    Parameters
    ----------
    line_text: str
        A line of text from a JDFTx out file

    Returns
    -------
    is_line: bool
        True if the line_text is the start of a Lowdin population analysis in a
        JDFTx out file
    """
    return "#--- Lowdin population analysis ---" in line_text


def is_charges_line(line_text: str) -> bool:
    """Return True if the line_text is start of charges log message.

    Return True if the line_text is a line of text from a JDFTx out file
    corresponding to a Lowdin population analysis.

    Parameters
    ----------
    line_text: str
        A line of text from a JDFTx out file

    Returns
    -------
    is_line: bool
        True if the line_text is a line of text from a JDFTx out file
        corresponding to a Lowdin population
    """
    return "oxidation-state" in line_text


def is_moments_line(line_text: str) -> bool:
    """Return True if the line_text is start of moments log message.

    Return True if the line_text is a line of text from a JDFTx out file
    corresponding to a Lowdin population analysis.

    Parameters
    ----------
    line_text: str
        A line of text from a JDFTx out file

    Returns
    -------
    is_line: bool
        True if the line_text is a line of text from a JDFTx out file
        corresponding to a Lowdin population
    """
    return "magnetic-moments" in line_text
