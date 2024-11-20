"""Module for JDFTx IO module utils.

Module for JDFTx IO module utils. Functions kept in this module are here if they are
used by multiple submodules, or if they are anticipated to be used by multiple
submodules in the future.

@mkhorton - this file is ready to review.
"""

from __future__ import annotations

from functools import wraps
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from monty.io import zopen

if TYPE_CHECKING:
    from collections.abc import Callable


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


# Named "t1" in unmet anticipation of multiple ways that a float would be needed
# to be read following the variable string with a colon.
def get_colon_var_t1(linetext: str, lkey: str) -> float | None:
    """Return float val from '...lkey: val...' in linetext.

    Read a float from an elec minimization line assuming value appears as
    "... lkey value ...".

    Parameters
    ----------
    linetext: str
        A line of text from a JDFTx out file
    lkey: str
        A string that appears before the float value in linetext. Must include
        the colon.
    """
    colon_var = None
    if lkey in linetext:
        colon_var = float(linetext.split(lkey)[1].strip().split(" ")[0])
    return colon_var


# This function matches the format of the generic "is_<x>_start_line" functions specific to
# the JOutStructure object initialization, but is not moved to the joutstructure module
# as it is also used in methods for other JDFTx IO modules (e.g. JOutStructures) so it
# is kept here to avoid circular imports.
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


def correct_geom_opt_type(opt_type: str | None) -> str | None:
    """Return recognizable opt_type string.

    Correct the opt_type string to match the JDFTx convention.

    Parameters
    ----------
    opt_type:
        The type of optimization step

    Returns
    -------
    opt_type: str | None
        The corrected type of optimization step
    """
    if opt_type is not None:
        if "lattice" in opt_type.lower():
            opt_type = "LatticeMinimize"
        elif "ionic" in opt_type.lower():
            opt_type = "IonicMinimize"
        else:
            opt_type = None
    return opt_type


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
    """Find first instance of key in output file.

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
    """Find last instance of key in output file.

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
    """Find all lines that exactly begin with key_input in a range of lines.

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
