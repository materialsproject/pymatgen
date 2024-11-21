"""Module for JDFTx IO module input utils.

Module for JDFTx IO module input utils. Functions kept in this module are here if they are
used by multiple submodules, or if they are anticipated to be used by multiple
submodules in the future.
"""

from __future__ import annotations

from typing import Any


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
