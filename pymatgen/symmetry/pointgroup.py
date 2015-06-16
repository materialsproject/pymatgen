# coding: utf-8

from __future__ import division, unicode_literals

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "5/8/13"


import pymatgen
import warnings
from monty.dev import deprecated


warnings.warn("The pointgroup module has been deprecated. All classes and "
              "functions have been moved to pymatgen.symmetry.analyzer. "
              "Please change your code. If you are looking for point group "
              "representations, check out pymatgen.symmetry.groups. This "
              "module will be removed in pmg v3.1.")

@deprecated(replacement=pymatgen.symmetry.analyzer.PointGroupAnalyzer,
            message="This class will be removed in pmg v3.1.")
class PointGroupAnalyzer(pymatgen.symmetry.analyzer.PointGroupAnalyzer):
    pass