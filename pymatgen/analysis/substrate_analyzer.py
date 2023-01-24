# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes to identify optimal substrates for film growth
"""

from __future__ import annotations

import warnings

from pymatgen.analysis.interfaces import SubstrateAnalyzer, ZSLGenerator  # noqa

__author__ = "Shyam Dwaraknath"
__copyright__ = "Copyright 2016, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Production"
__date__ = "Feb, 2016"


warnings.warn(
    "The substrate_analyzer module is being moved to the interfaces submodule in analysis."
    " These imports will break in Pymatgen 2023",
    category=FutureWarning,
    stacklevel=2,
)
