"""This module provides classes for predicting new structures from existing ones."""

from __future__ import annotations

import warnings

warnings.warn(
    "This module has been moved to the pymatgen.core.structure_prediction.substitutor module. "
    "This stub will be removed v2027.1. ",
    DeprecationWarning,
    stacklevel=2,
)

from pymatgen.core.structure_prediction.substitutor import *  # noqa: E402, F403
