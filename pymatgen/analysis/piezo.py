# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals
from __future__ import absolute_import

"""
This module provides classes for the Piezoelectric tensor
"""
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.elasticity.tensors import TensorBase
from pymatgen.analysis.elasticity import voigt_map
import numpy as np
import warnings
from six.moves import range

__author__ = "Shyam Dwaraknath"
__copyright__ = "Copyright 2016, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "Feb, 2016"


class PiezoTensor(TensorBase):
    """
    This class describes the 3x6 piezo tensor in Voigt-notation
    """

    def __new__(cls, input_array, tol=1e-3):
        """
        Create an PiezoTensor object.  The constructor throws an error if
        the shape of the input_matrix argument is not 3x3x3, i. e. in true
        tensor notation. Note that the constructor uses __new__ rather than
        __init__ according to the standard method of subclassing numpy
        ndarrays.

        Args:
            input_matrix (3x3x3 array-like): the 3x6 array-like
                representing the piezo tensor
        """
        obj = TensorBase(input_array).view(cls)
        if obj.shape != (3, 3, 3):
            raise ValueError("Default piezo tensor constructor requires "
                             "input argument to be the true 3x3x3 "
                             "array.  To construct from a 3x6 array, use "
                             "PiezoTensor.from_voigt")
        if not (obj - np.transpose(obj, (0, 2, 1)) < tol).all():
            warnings.warn("Input piezo tensor does "
                          "not satisfy standard symmetries")
        return obj
