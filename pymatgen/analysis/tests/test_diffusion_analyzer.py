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
import json

import numpy as np

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


class DiffusionAnalyzerTest(unittest.TestCase):

    def test_init(self):
        # Diffusion vasprun.xmls are rather large. We are only going to use a
        # very small preprocessed run for testing. Note that the results are
        # unreliable for short runs.
        with open(os.path.join(test_dir, "DiffusionAnalyzer.json")) as f:
            d = DiffusionAnalyzer.from_dict(json.load(f))
            self.assertAlmostEqual(d.conductivity, 66.5376433577, 7)
            self.assertAlmostEqual(d.diffusivity, 1.04185688594e-06, 7)
            self.assertTrue(np.allclose(
                d.conductivity_components,
                [63.08957759, 36.83734991, 100.08896311]))
            self.assertTrue(np.allclose(
                d.diffusivity_components,
                [9.87866530e-07, 5.76805019e-07, 1.56720873e-06]))
            self.assertAlmostEqual(d.max_framework_displacement, 1.1865683960)


if __name__ == '__main__':
    unittest.main()
