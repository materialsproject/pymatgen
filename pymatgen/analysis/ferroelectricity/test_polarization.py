from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import unittest
import os

import numpy as np
from pymatgen.analysis.elasticity.elastic import *
from pymatgen.analysis.elasticity.strain import Strain, IndependentStrain, Deformation
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.util.testing import PymatgenTest
from scipy.misc import central_diff_weights
import warnings
import json
import random
from six.moves import zip


class PolarizationTest(PymatgenTest):
    def setUp(self):
        pass
    def test_properties(self):
        pass