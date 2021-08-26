# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes to store, generate, and manipulate material interfaces.
"""

import warnings
from pymatgen.core.interface import Interface  # noqa
from pymatgen.analysis.interfaces import CoherentInterfaceBuilder  # noqa

__author__ = "Eric Sivonxay, Shyam Dwaraknath, and Kyle Bystrom"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Kyle Bystrom"
__email__ = "kylebystrom@gmail.com"
__date__ = "5/29/2019"
__status__ = "Prototype"

warnings.warn(
    "The substrate_analyzer module is being moved to the interfaces submodule in analysis."
    " These imports will break in Pymatgen 2023",
    category=FutureWarning,
    stacklevel=2,
)
