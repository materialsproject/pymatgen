"""The util package implements various utilities that are commonly used by various
packages.
"""

from __future__ import annotations

import numpy as np


def is_np_dict_equal(dict1, dict2, /) -> bool:
    """Compare two dict whose value could be np arrays.

    Args:
        dict1 (dict): The first dict.
        dict2 (dict): The second dict.

    Returns:
        bool: Whether these two dicts are equal.
    """
    if dict1.keys() != dict2.keys():
        return False

    return all(np.array_equal(dict1[key], dict2[key]) for key in dict1)
