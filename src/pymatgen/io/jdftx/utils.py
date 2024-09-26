"""Module for JDFTx utils."""

from __future__ import annotations

from functools import wraps
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from monty.io import zopen

if TYPE_CHECKING:
    from collections.abc import Callable


def gather_jeiters_line_collections(iter_type: str, text_slice: list[str]) -> tuple[list[list[str]], list[str]]:
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


############################################
# HELPERS FOR JDFTXOUTFILE #
############################################


def check_file_exists(func: Callable) -> Any:
    """Check if file exists.

    Check if file exists (and continue normally) or raise an exception if
    it does not.
    """

    @wraps(func)
    def wrapper(filename: str) -> Any:
        filepath = Path(filename)
        if not filepath.is_file():
            raise OSError(f"'{filename}' file doesn't exist!")
        return func(filename)

    return wrapper


@check_file_exists
def read_file(file_name: str) -> list[str]:
    """
    Read file into a list of str.

    Parameters
    ----------
    filename: Path or str
        name of file to read

    Returns
    -------
    text: list[str]
        list of strings from file
    """
    with zopen(file_name, "r") as f:
        text = f.readlines()
    f.close()
    return text


def read_outfile_slices(file_name: str) -> list[list[str]]:
    """
    Read slice of out file into a list of str.

    Parameters
    ----------
    filename: Path or str
        name of file to read
    out_slice_idx: int
        index of slice to read from file

    Returns
    -------
    texts: list[list[str]]
        list of out file slices (individual calls of JDFTx)
    """
    _text = read_file(file_name)
    start_lines = get_start_lines(_text, add_end=True)
    texts = []
    for i in range(len(start_lines) - 1):
        text = _text[start_lines[i] : start_lines[i + 1]]
        texts.append(text)
    return texts


############################################
# HELPERS FOR JDFTXINFILE #
############################################


def multi_hasattr(varbase: Any, varname: str):
    """Check if object has an attribute (capable of nesting with . splits).

    Check if object has an attribute (capable of nesting with . splits).

    Parameters
    ----------
    varbase
        Object to check.
    varname
        Attribute to check for.

    Returns
    -------
    bool
        Whether the object has the attribute.
    """
    varlist = varname.split(".")
    for i, var in enumerate(varlist):
        if i == len(varlist) - 1:
            return hasattr(varbase, var)
        if hasattr(varbase, var):
            varbase = getattr(varbase, var)
        else:
            return False
    return None


def multi_getattr(varbase: Any, varname: str):
    """Check if object has an attribute (capable of nesting with . splits).

    Check if object has an attribute (capable of nesting with . splits).

    Parameters
    ----------
    varbase
        Object to check.
    varname
        Attribute to check for.

    Returns
    -------
    Any
        Attribute of the object.
    """
    if not multi_hasattr(varbase, varname):
        raise AttributeError(f"{varbase} does not have attribute {varname}")
    varlist = varname.split(".")
    for var in varlist:
        varbase = getattr(varbase, var)
    return varbase


############################################
# HELPERS FOR GENERIC_TAGS #
############################################


def flatten_list(tag: str, list_of_lists: list[Any]) -> list[Any]:
    """Flatten list of lists into a single list, then stop.

    Flatten list of lists into a single list, then stop.

    Parameters
    ----------
    tag : str
        The tag to flatten the list of lists for.
    list_of_lists : list[Any]
        The list of lists to flatten.

    Returns
    -------
    list[Any]
        The flattened list.
    """
    if not isinstance(list_of_lists, list):
        raise TypeError(f"{tag}: You must provide a list to flatten_list()!")
    flist = []
    for v in list_of_lists:
        if isinstance(v, list):
            flist.extend(flatten_list(tag, v))
        else:
            flist.append(v)
    return flist


############################################
# HELPERS FOR JOUTSTRUCTURE(S) #
############################################


def _brkt_list_of_3_to_nparray(line: str) -> np.ndarray:
    """Return 3x1 numpy array.

    Convert a string of the form "[ x y z ]" to a 3x1 numpy array

    Parameters
    ----------
    line: str
        A string of the form "[ x y z ]"
    """
    return np.array([float(x) for x in line.split()[1:-1]])


def _brkt_list_of_3x3_to_nparray(lines: list[str], i_start: int = 0) -> np.ndarray:
    """Return 3x3 numpy array.

    Convert a list of strings of the form "[ x y z ]" to a 3x3 numpy array

    Parameters
    ----------
    lines: list[str]
        A list of strings of the form "[ x y z ]"
    i_start: int
        The index of the first line in lines

    Returns
    -------
    out: np.ndarray
        A 3x3 numpy array
    """
    out = np.zeros([3, 3])
    for i in range(3):
        out[i, :] += _brkt_list_of_3_to_nparray(lines[i + i_start])
    return out


elec_min_start_flag: str = "-------- Electronic minimization -----------"


def get_colon_var_t1(linetext: str, lkey: str) -> float | None:
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


############################################
# HELPERS FOR JDFTXOUTFILESLICE #
############################################


class ClassPrintFormatter:
    """Generic class object print formatter.

    Generic class object print formatter.
    """

    def __str__(self) -> str:
        """Return class object as str for readable format in command line."""
        return (
            str(self.__class__)
            + "\n"
            + "\n".join(str(item) + " = " + str(self.__dict__[item]) for item in sorted(self.__dict__))
        )


def get_start_lines(
    text: list[str],
    start_key: str = "*************** JDFTx",
    add_end: bool = False,
) -> list[int]:
    """Get start line numbers for JDFTx calculations.

    Get the line numbers corresponding to the beginning of separate JDFTx calculations
    (in case of multiple calculations appending the same out file).

    Args:
        text: output of read_file for out file
    """
    start_lines = []
    i = None
    for i, line in enumerate(text):
        if start_key in line:
            start_lines.append(i)
    if add_end and i is not None:
        start_lines.append(i)
    if i is None:
        raise ValueError("Outfile parser fed an empty file.")
    if not len(start_lines):
        raise ValueError("No JDFTx calculations found in file.")
    return start_lines


def find_key_first(key_input: str, tempfile: list[str]) -> int | None:
    """

    Find first instance of key in output file.

    Parameters
    ----------
    key_input: str
        key string to match
    tempfile: List[str]
        output from readlines() function in read_file method
    """
    key_input = str(key_input)
    line = None
    for i in range(len(tempfile)):
        if key_input in tempfile[i]:
            line = i
            break
    return line


def find_key(key_input: str, tempfile: list[str]) -> int | None:
    """

    Find last instance of key in output file.

    Parameters
    ----------
    key_input: str
        key string to match
    tempfile: List[str]
        output from readlines() function in read_file method
    """
    key_input = str(key_input)
    line = None
    lines = find_all_key(key_input, tempfile)
    if len(lines):
        line = lines[-1]
    return line


def find_first_range_key(
    key_input: str,
    tempfile: list[str],
    startline: int = 0,
    endline: int = -1,
    skip_pound: bool = False,
) -> list[int]:
    """

    Find all lines that exactly begin with key_input in a range of lines.

    Parameters
    ----------
    key_input: str
        key string to match
    tempfile: List[str]
        output from readlines() function in read_file method
    startline: int
        line to start searching from
    endline: int
        line to stop searching at
    skip_pound: bool
        whether to skip lines that begin with a pound sign

    Returns
    -------
    L: list[int]
        list of line numbers where key_input occurs

    """
    key_input = str(key_input)
    startlen = len(key_input)
    line_list = []

    if endline == -1:
        endline = len(tempfile)
    for i in range(startline, endline):
        line = tempfile[i]
        if skip_pound:
            for _ in range(10):  # repeat to make sure no really weird formatting
                line = line.lstrip()
                line = line.lstrip("#")
        line = line[0:startlen]
        if line == key_input:
            line_list.append(i)
    return line_list


def key_exists(key_input: str, tempfile: list[str]) -> bool:
    """Check if key_input exists in tempfile.

    Search through tempfile for key_input. Return True if found,
    False otherwise.

    Parameters
    ----------
    key_input: str
        key string to match
    tempfile: List[str]
        output from readlines() function in read_file method

    Returns
    -------
    bool
        True if key_input exists in tempfile, False otherwise
    """
    line = find_key(key_input, tempfile)
    return line is not None


def find_all_key(key_input: str, tempfile: list[str], startline: int = 0) -> list[int]:
    """Find all lines containing key_input.

    Search through tempfile for all lines containing key_input. Returns a list
    of line numbers.

    Parameters
    ----------
    key_input: str
        key string to match
    tempfile: List[str]
        output from readlines() function in read_file method
    startline: int
        line to start searching from

    Returns
    -------
    line_list: list[int]
        list of line numbers where key_input occurs
    """
    return [i for i in range(startline, len(tempfile)) if key_input in tempfile[i]]


def get_pseudo_read_section_bounds(text: list[str]) -> list[list[int]]:
    """Get the boundary line numbers for the pseudopotential read section.

    Get the boundary line numbers for the pseudopotential read section.

    Parameters
    ----------
    text: list[str]
        output of read_file for out file

    Returns
    -------
    section_bounds: list[list[int]]
        list of line numbers for the pseudopotential read sections
    """
    start_lines = find_all_key("Reading pseudopotential file", text)
    section_bounds = []
    for start_line in start_lines:
        bounds = [start_line]
        for i in range(start_line, len(text)):
            if not len(text[i].strip()):
                bounds.append(i)
                break
        section_bounds.append(bounds)
    return section_bounds
