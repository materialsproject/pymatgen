# coding: utf-8

from __future__ import division, unicode_literals

"""
TODO: Change the module doc.
"""


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
    get_conversion_factor, fit_arrhenius
import pymatgen.core.physical_constants as phyc
from pymatgen.io.smartio import read_structure
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class FuncTest(unittest.TestCase):

    def test_get_conversion_factor(self):
        filepath = os.path.join(test_dir, 'LiFePO4.cif')
        s = read_structure(filepath)
        self.assertAlmostEqual(41370704.1173,
                               get_conversion_factor(s, "Li", 600), 4)

    def test_fit_arrhenius(self):
        Ea = 0.5
        k = phyc.k_b / phyc.e
        c = 12
        temps = np.array([300, 1000, 500])
        diffusivities = c * np.exp(-Ea/(k * temps))
        r = fit_arrhenius(temps, diffusivities)
        self.assertAlmostEqual(r[0], Ea)
        self.assertAlmostEqual(r[1], c)


class DiffusionAnalyzerTest(PymatgenTest):

    def test_init(self):
        # Diffusion vasprun.xmls are rather large. We are only going to use a
        # very small preprocessed run for testing. Note that the results are
        # unreliable for short runs.
        with open(os.path.join(test_dir, "DiffusionAnalyzer.json")) as f:
            d = DiffusionAnalyzer.from_dict(json.load(f))

            self.assertAlmostEqual(d.conductivity, 74.1362195972, 7)
            self.assertAlmostEqual(d.diffusivity,  1.16083658794e-06, 7)
            self.assertAlmostEqual(d.conductivity_std_dev, 8.301909069566328, 7)
            self.assertAlmostEqual(d.diffusivity_std_dev, 1.29992598086e-07, 7)
            self.assertArrayAlmostEqual(
                d.conductivity_components,
                [47.8728896, 31.3098319, 143.47106767])
            self.assertArrayAlmostEqual(
                d.diffusivity_components,
                [7.49601236e-07, 4.90254273e-07, 2.24649255e-06])
            self.assertArrayAlmostEqual(
                d.conductivity_components_std_dev,
                [8.16076457,  22.74144339, 20.64816641]
            )
            self.assertArrayAlmostEqual(
                d.diffusivity_components_std_dev,
                [1.27782535e-07, 3.56089098e-07, 3.23312238e-07]
            )

            self.assertArrayAlmostEqual(d.framework_ion_displacement,
                                        [0.65556702781541298,
                                         0.59869611646057463,
                                         0.56390914443090445,
                                         0.61660041929540588,
                                         0.59979115804226046,
                                         0.4374606277579815, 1.1865683960470783,
                                         0.90170643716765908,
                                         0.66448403678537671,
                                         1.0346375380664645, 0.6177630142863979,
                                         0.7952002051914302,
                                         0.73426861230540108,
                                         0.78580479569055772,
                                         0.55707323690656607,
                                         1.0942937746885417,
                                         0.65093723953087879,
                                         1.0876687380413455,
                                         0.70581621847249998,
                                         0.82983063175985849,
                                         0.78139137476213427,
                                         0.73376552320561528,
                                         0.90571616162367463,
                                         0.59790930931869191,
                                         0.68303335869850146,
                                         0.79265008940846282,
                                         0.67651800099886084,
                                         0.85558660329689984,
                                         0.71308709164223705,
                                         0.76210076957907491]
            )

            self.assertAlmostEqual(d.max_framework_displacement, 1.18656839605)
            d = DiffusionAnalyzer.from_dict(d.as_dict())
            self.assertIsInstance(d, DiffusionAnalyzer)

            #Ensure summary dict is json serializable.
            json.dumps(d.get_summary_dict(include_msd_t=True))

            d = DiffusionAnalyzer(d.structure, d.disp, d.specie, d.temperature,
                                  d.time_step, d.step_skip, smoothed=True,
                                  weighted=True)
            self.assertAlmostEqual(d.conductivity, 74.16537220815061, 7)
            self.assertAlmostEqual(d.diffusivity, 1.14606446822e-06, 7)

            d = DiffusionAnalyzer(d.structure, d.disp, d.specie, d.temperature,
                                  d.time_step, d.step_skip, smoothed=False)
            self.assertAlmostEqual(d.conductivity, 27.2047915553, 7)
            self.assertAlmostEqual(d.diffusivity, 4.25976905436e-07, 7)


if __name__ == '__main__':
    unittest.main()
