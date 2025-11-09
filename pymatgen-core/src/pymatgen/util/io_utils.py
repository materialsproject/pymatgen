"""This module provides utility classes for io operations."""

from __future__ import annotations

import re
import warnings
from typing import TYPE_CHECKING

from monty.io import zopen

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator
    from pathlib import Path
    from typing import Any

__author__ = "Shyue Ping Ong, Rickard Armiento, Anubhav Jain, G Matteo, Ioannis Petousis"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Sep 23, 2011"


def clean_lines(
    string_list: list[str],
    remove_empty_lines: bool = True,
    rstrip_only: bool = False,
) -> Iterator[str]:
    """Remove leading and trailing whitespaces from a list of strings.

    Args:
        string_list (list[str]): List of strings.
        remove_empty_lines (bool): Set to True to skip lines which are empty after
            stripping.
        rstrip_only (bool): Set to True to strip trailing whitespaces only (i.e.,
            to retain leading whitespaces). Defaults to False.

    Yields:
        str: clean string with no leading and trailing whitespaces.
    """
    for string in string_list:
        clean_string: str = string
        if "#" in string:
            clean_string = string[: string.index("#")]

        clean_string = clean_string.rstrip() if rstrip_only else clean_string.strip()

        if (not remove_empty_lines) or clean_string != "":
            yield clean_string


def micro_pyawk(
    filename: str | Path,
    search: list[tuple[re.Pattern | str, Callable, Callable]],
    results: Any | None = None,
    debug: Callable | None = None,
    postdebug: Callable | None = None,
) -> Any:
    """Small awk-mimicking search routine.

    This function goes through each line in the file, and if `regex` matches that
    line AND test(results, line) is True (OR test is None) we execute
    run(results, match), where match is the Match object from running
    Pattern.match.

    Args:
        filename (PathLike): The file to search through.
        search (list[tuple[Pattern | str, Callable, Callable]]): The "search program" of
            3 elements, i.e. [(regex, test, run), ...].
            Here `regex` is either a Pattern object, or a string that we compile
            into a Pattern.
        results: An object to store results. Default as an empty dictionary.
            Passing a results object let you interact with it via `run` and `test`.
            Hence, in many occasions it is clever to use the instance itself as results.
        debug (Callable): Debug `run`.
        postdebug (Callable): Post debug `run` after debug `run`.

    Returns:
        Any: The updated `results` object.

    Author: Rickard Armiento, Ioannis Petousis
    """
    # TODO: remove `debug` and `postdebug` after 2025-11-09 if no one is opposing
    if debug is not None:
        warnings.warn("arg debug is scheduled for removal, see PR4160", DeprecationWarning, stacklevel=2)
    if postdebug is not None:
        warnings.warn("arg postdebug is scheduled for removal, see PR4160", DeprecationWarning, stacklevel=2)

    if results is None:
        results = {}

    # Compile regex strings to Patterns
    searches: list[tuple[re.Pattern, Callable, Callable]] = [
        (re.compile(regex), test, run) for regex, test, run in search
    ]

    with zopen(filename, mode="rt", encoding="utf-8") as file:
        for line in file:
            for regex, test, run in searches:
                match = re.search(regex, line)

                if match is not None and (test is None or test(results, line)):
                    if debug is not None:
                        debug(results, match)

                    run(results, match)
                    if postdebug is not None:
                        postdebug(results, match)

    return results
