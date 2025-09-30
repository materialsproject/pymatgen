from __future__ import annotations

import re
from enum import Enum
from typing import Any

from pymatgen.io.lobster.future.constants import LOBSTER_ORBITALS


def natural_sort(list_to_sort: str) -> list[Any]:
    """Sort a list of strings in human order.

    This function sorts strings in a way that humans would expect,
    taking into account numerical values within the strings.

    Args:
        list_to_sort (str): List of strings to sort.

    Returns:
        list[Any]: Sorted list of strings in human order.

    Example:
        >>> natural_sort(["file10", "file2", "file1"])
        ['file1', 'file2', 'file10']
    """
    return [
        int(text) if text.isdigit() else text.lower()
        for text in re.split(r"(\d+)", list_to_sort)
    ]


def parse_orbital_from_text(text: str) -> str | None:
    """Parse the orbital from a text string of an ICOXXLIST file.

    This function extracts the orbital information from a given text string.
    It checks for valid orbital patterns and returns the matched orbital.

    Args:
        text (str): Text string to parse the orbital from.

    Returns:
        str | None: Parsed orbital string if a valid orbital is found,
        otherwise None.

    Example:
        >>> parse_orbital_from_text("1s_2p_x")
        '2p_x'
    """
    parts = text.split("_")

    if len(parts) == 1:
        return None

    for orbital in LOBSTER_ORBITALS:
        if match := re.search(rf"\d+{re.escape(orbital)}", "_".join(parts[-2:])):
            return match.group(0)

    return parts[-1] if re.match(r"\d+[a-z]+", parts[-1]) else None


def make_json_compatible(obj: Any) -> Any:
    """Convert an object to a JSON-compatible format recursively.

    This function ensures that the input object is converted into a format
    that can be serialized into JSON. It handles lists, tuples, dictionaries,
    and enums.

    Args:
        obj (Any): Input object to convert.

    Returns:
        Any: JSON-compatible representation of the input object.

    Example:
        >>> make_json_compatible({"key": Enum("Example", "value")})
        {'key': 'value'}
    """
    if isinstance(obj, (list, tuple)):
        return [make_json_compatible(item) for item in obj]

    if isinstance(obj, dict):
        new_dict = {}

        for k, v in obj.items():
            new_key = str(k.value) if isinstance(k, Enum) else k
            new_dict[new_key] = make_json_compatible(v)

        return new_dict

    return obj
