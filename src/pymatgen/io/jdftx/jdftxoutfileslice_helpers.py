"""Module containing helper functions for parsing JDFTx output files.

This module contains helper functions for parsing JDFTx output files.
"""

from __future__ import annotations


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
    for i, line in enumerate(text):
        if start_key in line:
            start_lines.append(i)
    if add_end:
        start_lines.append(i)
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
    if not line_list:
        line_list = [len(tempfile)]
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
