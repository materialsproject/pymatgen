# coding: utf-8

from __future__ import division, unicode_literals

from monty.dev import deprecated
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


@deprecated(replacement=SpacegroupAnalyzer)
class SymmetryFinder(SpacegroupAnalyzer):
    pass
