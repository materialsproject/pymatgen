"""This module provides utility classes for io operations."""

from __future__ import annotations

import os
import re
from typing import TYPE_CHECKING

from monty.io import zopen

if TYPE_CHECKING:
    from collections.abc import Generator

__author__ = "Shyue Ping Ong, Rickard Armiento, Anubhav Jain, G Matteo, Ioannis Petousis"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Sep 23, 2011"


def clean_lines(
    string_list,
    remove_empty_lines=True,
    rstrip_only=False,
) -> Generator[str, None, None]:
    """Strips whitespace, carriage returns and empty lines from a list of strings.

    Args:
        string_list: List of strings
        remove_empty_lines: Set to True to skip lines which are empty after
            stripping.
        rstrip_only: Set to True to strip trailing whitespaces only (i.e.,
            to retain leading whitespaces). Defaults to False.

    Yields:
        list: clean strings with no whitespaces. If rstrip_only == True,
            clean strings with no trailing whitespaces.
    """
    for s in string_list:
        clean_s = s
        if "#" in s:
            ind = s.index("#")
            clean_s = s[:ind]
        clean_s = clean_s.rstrip() if rstrip_only else clean_s.strip()
        if (not remove_empty_lines) or clean_s != "":
            yield clean_s


def micro_pyawk(filename, search, results=None, debug=None, postdebug=None):
    """Small awk-mimicking search routine.

    'file' is file to search through.
    'search' is the "search program", a list of lists/tuples with 3 elements;
    i.e. [[regex, test, run], [regex, test, run], ...]
    'results' is a an object that your search program will have access to for
    storing results.

    Here regex is either as a Regex object, or a string that we compile into a
    Regex. test and run are callable objects.

    This function goes through each line in filename, and if regex matches that
    line *and* test(results,line)==True (or test is None) we execute
    run(results,match), where match is the match object from running
    Regex.match.

    The default results is an empty dictionary. Passing a results object let
    you interact with it in run() and test(). Hence, in many occasions it is
    thus clever to use results=self.

    Author: Rickard Armiento, Ioannis Petousis

    Returns:
        dict[str, Any]: The results dictionary.
    """
    if results is None:
        results = {}

    # Compile regex strings
    for entry in search:
        entry[0] = re.compile(entry[0])

    with zopen(filename, mode="rt") as file:
        for line in file:
            for entry in search:
                match = re.search(entry[0], line)
                if match and (entry[1] is None or entry[1](results, line)):
                    if debug is not None:
                        debug(results, match)
                    entry[2](results, match)
                    if postdebug is not None:
                        postdebug(results, match)

    return results


umask = os.umask(0)
os.umask(umask)
