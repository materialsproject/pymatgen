"""This module implements Compatibility corrections for mixing runs of different
functionals.
"""

from __future__ import annotations

import warnings

from pymatgen.analysis.compatibility.mixing_scheme import *  # noqa: F403

warnings.warn(
    "MaterialsProjectDFTMixingScheme has been moved to the pymatgen.analysis.compatibility.mixing_scheme "
    "module. This stub will be removed v2027.1.",
    DeprecationWarning,
    stacklevel=2,
)
