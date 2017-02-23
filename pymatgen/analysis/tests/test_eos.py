# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals, division, print_function

import numpy as np
import unittest

from pymatgen.analysis.eos import EOS, numerical_eos


class EOSTest(unittest.TestCase):

    def setUp(self):
        # Si data from Cormac
        self.volumes = [25.987454833, 26.9045702104, 27.8430241908,
                        28.8029649591, 29.7848370694, 30.7887887064,
                        31.814968055, 32.8638196693, 33.9353435494,
                        35.0299842495, 36.1477417695, 37.2892088485,
                        38.4543854865, 39.6437162376, 40.857201102,
                        42.095136449, 43.3579668329, 44.6456922537,
                        45.9587572656, 47.2973100535, 48.6614988019,
                        50.0517680652, 51.4682660281, 52.9112890601,
                        54.3808371612, 55.8775030703, 57.4014349722,
                        58.9526328669]
        self.energies = [-7.63622156576, -8.16831294894, -8.63871612686,
                         -9.05181213218, -9.41170988374, -9.72238224345,
                         -9.98744832526, -10.210309552, -10.3943401353,
                         -10.5427238068, -10.6584266073, -10.7442240979,
                         -10.8027285713, -10.8363890521, -10.8474912964,
                         -10.838157792, -10.8103477586, -10.7659387815,
                         -10.7066179666, -10.6339907853, -10.5495538639,
                         -10.4546677714, -10.3506386542, -10.2386366017,
                         -10.1197772808, -9.99504030111, -9.86535084973,
                         -9.73155247952]

    def test_run_all_models(self):
        for eos_name in EOS.MODELS:
            eos = EOS(eos_name=eos_name)
            _ = eos.fit(self.volumes, self.energies)

    def test_numerical_eos(self):
        e0, b0, b1, v0, fit_coeffs = numerical_eos(self.volumes, self.energies)
        eos = EOS(eos_name="numerical_eos")
        fit = eos.fit(self.volumes, self.energies)
        self.assertGreater(len(fit_coeffs), 3)
        self.assertEqual(e0, fit.e0)
        self.assertEqual(v0, fit.v0)
        self.assertEqual(b0, fit.b0)
        self.assertEqual(b1, fit.b1)
        self.assertEqual(fit_coeffs, fit.eos_params)
        np.testing.assert_almost_equal(e0, -10.847, decimal=3)
        np.testing.assert_almost_equal(v0, 40.877, decimal=3)
        np.testing.assert_almost_equal(b0, 0.5545, decimal=2)
        np.testing.assert_almost_equal(fit.b0_GPa, 88.8456, decimal=2)
        np.testing.assert_almost_equal(b1, 4.3473, decimal=2)


if __name__ == "__main__":
    unittest.main()
