"""
This module provides classes and methods used to describe deformations and
strains, including applying those deformations to structure objects and
generating deformed structure sets for further calculations.
"""

from __future__ import annotations

import warnings

warnings.warn(
    "This module has been moved to the pymatgen.core.elasticity.strain module. This stub will be removed v2027.1. ",
    DeprecationWarning,
    stacklevel=2,
)

from pymatgen.core.elasticity.strain import *  # noqa: E402, F403
