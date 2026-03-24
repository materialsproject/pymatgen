"""This module calculates corrections for the species listed below, fitted to the experimental and computed
entries given to the CorrectionCalculator constructor.
"""

from __future__ import annotations

import warnings

from pymatgen.analysis.compatibility.correction_calculator import *  # noqa: F403

warnings.warn(
    "All CorrectionCalculators have been moved to the "
    "pymatgen.analysis.compatibility.correction_calculator module. This stub will be removed v2027.1.",
    DeprecationWarning,
    stacklevel=2,
)
