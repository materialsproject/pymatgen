#!/usr/bin/env python

"""
TODO: Change the module doc.
"""

from __future__ import division

__author__ = "shyuepingong"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Beta"
__date__ = "5/2/13"

import unittest
import os

from pymatgen.analysis.diffusion_analyzer import DiffusionAnalyzer,\
    get_conversion_factor
from pymatgen.io.smartio import read_structure

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class FuncTest(unittest.TestCase):

    def test_get_conversion_factor(self):
        filepath = os.path.join(test_dir, 'LiFePO4.cif')
        s = read_structure(filepath)
        self.assertAlmostEqual(41370704.1173,
                               get_conversion_factor(s, "Li", 600), 4)

#TODO: Write unittest for DiffusionAnalyzer

if __name__ == '__main__':
    unittest.main()
