"""
This module implements a EnergyModel abstract class and some basic
implementations. Basically, an EnergyModel is any model that returns an
"energy" for any given structure.
"""

from __future__ import annotations

import warnings

warnings.warn(
    "This module has been moved to the pymatgen.core.energy_models module. This stub will be removed v2027.1. ",
    DeprecationWarning,
    stacklevel=2,
)


from pymatgen.core.energy_models import *  # noqa: E402, F403
