# coding: utf-8

from __future__ import division, unicode_literals

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "5/8/13"


import pymatgen
from monty.dev import deprecated


@deprecated(replacement=pymatgen.symmetry.analyzer.PointGroupAnalyzer)
class PointGroupAnalyzer(pymatgen.symmetry.analyzer.PointGroupAnalyzer):
    pass