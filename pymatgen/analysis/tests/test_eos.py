# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals, division, print_function

import numpy as np
import unittest

from pymatgen.analysis.eos import EOS, NumericalEOS
from pymatgen.util.testing import PymatgenTest


class EOSTest(PymatgenTest):

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
        num_eos = EOS(eos_name="numerical_eos")
        self.num_eos_fit = num_eos.fit(self.volumes, self.energies)

    def test_run_all_models(self):
        for eos_name in EOS.MODELS:
            eos = EOS(eos_name=eos_name)
            _ = eos.fit(self.volumes, self.energies)

    def test_numerical_eoswrapper(self):
        # using numerical eos directly vs via EOS wrapper
        numerical_eos = NumericalEOS(self.volumes, self.energies)
        numerical_eos.fit()
        self.assertGreater(len(numerical_eos.eos_params), 3)
        self.assertAlmostEqual(float(numerical_eos.e0), self.num_eos_fit.e0, 3)
        self.assertAlmostEqual(float(numerical_eos.v0), self.num_eos_fit.v0, 3)
        self.assertAlmostEqual(float(numerical_eos.b0), self.num_eos_fit.b0, 3)
        self.assertAlmostEqual(float(numerical_eos.b1), self.num_eos_fit.b1, 3)
        self.assertArrayAlmostEqual(numerical_eos.eos_params, self.num_eos_fit.eos_params)

    def test_numerical_eos_values(self):
        np.testing.assert_almost_equal(self.num_eos_fit.e0, -10.84749, decimal=3)
        np.testing.assert_almost_equal(self.num_eos_fit.v0, 40.857201, decimal=1)
        np.testing.assert_almost_equal(self.num_eos_fit.b0, 0.55, decimal=2)
        #np.testing.assert_almost_equal(self.num_eos_fit.b0_GPa, 89.0370727, decimal=1)
        #np.testing.assert_almost_equal(self.num_eos_fit.b1, 4.344039, decimal=2)

    def test_eos_func(self):
        # list vs np.array arguments
        np.testing.assert_almost_equal(self.num_eos_fit.func([0,1,2]),
                                       self.num_eos_fit.func(np.array([0,1,2])),
                                       decimal=10)
        # func vs _func
        np.testing.assert_almost_equal(self.num_eos_fit.func(0.),
                                       self.num_eos_fit._func(
                                           0., self.num_eos_fit.eos_params),
                                       decimal=10)
        # test the eos function: energy = f(volume)
        # numerical eos evaluated at volume=0 == a0 of the fit polynomial
        np.testing.assert_almost_equal(self.num_eos_fit.func(0.),
                                       self.num_eos_fit.eos_params[-1], decimal=6)
        birch_eos = EOS(eos_name="birch")
        birch_eos_fit = birch_eos.fit(self.volumes, self.energies)
        # birch eos evaluated at v0 == e0
        np.testing.assert_almost_equal(birch_eos_fit.func(birch_eos_fit.v0),
                                       birch_eos_fit.e0, decimal=6)

        # TODO: Reactivate
        #fig = birch_eos_fit.plot_ax(ax=None, show=False, fontsize=8, title="birch eos")
        #self.assertTrue(hasattr(fig, "savefig"))

    def test_eos_func_call(self):
        # eos_fit_obj.func(volume) == eos_fit_obj(volume)
        np.testing.assert_almost_equal(self.num_eos_fit.func(0.),
                                       self.num_eos_fit(0.), decimal=10)

    def test_summary_dict(self):
        d = {"e0": self.num_eos_fit.e0, "b0": self.num_eos_fit.b0,
             "b1": self.num_eos_fit.b1, "v0": self.num_eos_fit.v0}
        self.assertDictEqual(self.num_eos_fit.results, d)


if __name__ == "__main__":
    unittest.main()
