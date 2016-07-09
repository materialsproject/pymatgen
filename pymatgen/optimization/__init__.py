# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

try:
  import linear_assignment
except ImportError:
  import linear_assignment_python as linear_assignment  # use pure-python fallback
