"""Helper functions for creation of JOutStructure(s).

This module contains helper functions for the creation of JOutStructure(s) from
the output files of JDFTx calculations.
"""

from __future__ import annotations

elec_min_start_flag: str = "-------- Electronic minimization -----------"


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
            # else:
            #     warnmsg = f"Line {i-1} ({out_slice[i-1]}) triggered \
            #         end_started, but following line is empty. Final step_bounds \
            #             may be incorrect. "
            #     warnings.warn(warnmsg, stacklevel=2)
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


def correct_geom_iter_type(iter_type: str | None) -> str | None:
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


def is_magnetic_moments_line(line_text: str) -> bool:
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
