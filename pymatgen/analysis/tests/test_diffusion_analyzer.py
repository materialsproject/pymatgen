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
import random
import numpy as np

from pymatgen.analysis.diffusion_analyzer import DiffusionAnalyzer,\
    get_conversion_factor, fit_arrhenius
import pymatgen.core.physical_constants as phyc
from pymatgen.core.structure import Structure
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class FuncTest(unittest.TestCase):

    def test_get_conversion_factor(self):
        filepath = os.path.join(test_dir, 'LiFePO4.cif')
        s = Structure.from_file(filepath)
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

            self.assertArrayAlmostEqual(
                d.max_ion_displacements,
                [1.4620659693989553, 1.2787303484445025, 3.419618540097756,
                 2.340104469126246, 2.6080973517594233, 1.3928579365672844,
                 1.3561505956708932, 1.6699242923686253, 1.0352389639563648,
                 1.1662520093955808, 1.2322019205885841, 0.8094210554832534,
                 1.9917808504954169, 1.2684148391206396, 2.392633794162402,
                 2.566313049232671, 1.3175030435622759, 1.4628945430952793,
                 1.0984921286753002, 1.2864482076554093, 0.655567027815413,
                 0.5986961164605746, 0.5639091444309045, 0.6166004192954059,
                 0.5997911580422605, 0.4374606277579815, 1.1865683960470783,
                 0.9017064371676591, 0.6644840367853767, 1.0346375380664645,
                 0.6177630142863979, 0.7952002051914302, 0.7342686123054011,
                 0.7858047956905577, 0.5570732369065661, 1.0942937746885417,
                 0.6509372395308788, 1.0876687380413455, 0.7058162184725,
                 0.8298306317598585, 0.7813913747621343, 0.7337655232056153,
                 0.9057161616236746, 0.5979093093186919, 0.6830333586985015,
                 0.7926500894084628, 0.6765180009988608, 0.8555866032968998,
                 0.713087091642237, 0.7621007695790749])

            self.assertEqual(d.sq_disp_ions.shape, (50, 206))

            self.assertAlmostEqual(d.max_framework_displacement, 1.18656839605)

            ss = list(d.get_drift_corrected_structures())
            self.assertEqual(len(ss), 1000)
            n = random.randint(0, 999)
            self.assertArrayAlmostEqual(
                ss[n].cart_coords - d.structure.cart_coords + d.drift[:, n, :],
                d.disp[:, n, :])

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
