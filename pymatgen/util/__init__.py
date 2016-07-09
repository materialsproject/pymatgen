# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

"""
The util package implements various utilities that are commonly used by various
packages.
"""

__author__ = "Shyue"
__date__ = "$Jun 6, 2011 7:30:05 AM$"

try:
  import coord_utils_cython
except ImportError:
  import coord_utils_python as coord_utils_cython  # use pure-python fallback
