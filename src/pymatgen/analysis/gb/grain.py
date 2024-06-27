from __future__ import annotations

import warnings

from pymatgen.core.interface import *  # noqa: F403

warnings.warn(
    "Grain boundary analysis has been moved to pymatgen.core.interface."
    "This stub is retained for backwards compatibility and will be removed Dec 31 2024.",
    DeprecationWarning,
)
