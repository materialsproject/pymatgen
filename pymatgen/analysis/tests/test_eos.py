# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from __future__ import annotations

import unittest

import numpy as np
from pytest import approx

from pymatgen.analysis.eos import EOS, NumericalEOS
from pymatgen.util.testing import PymatgenTest


class EOSTest(PymatgenTest):
    def setUp(self):
        # Si data from Cormac
        self.volumes = [
            25.987454833,
            26.9045702104,
            27.8430241908,
            28.8029649591,
            29.7848370694,
            30.7887887064,
            31.814968055,
            32.8638196693,
            33.9353435494,
            35.0299842495,
            36.1477417695,
            37.2892088485,
            38.4543854865,
            39.6437162376,
            40.857201102,
            42.095136449,
            43.3579668329,
            44.6456922537,
            45.9587572656,
            47.2973100535,
            48.6614988019,
            50.0517680652,
            51.4682660281,
            52.9112890601,
            54.3808371612,
            55.8775030703,
            57.4014349722,
            58.9526328669,
        ]
        self.energies = [
            -7.63622156576,
            -8.16831294894,
            -8.63871612686,
            -9.05181213218,
            -9.41170988374,
            -9.72238224345,
            -9.98744832526,
            -10.210309552,
            -10.3943401353,
            -10.5427238068,
            -10.6584266073,
            -10.7442240979,
            -10.8027285713,
            -10.8363890521,
            -10.8474912964,
            -10.838157792,
            -10.8103477586,
            -10.7659387815,
            -10.7066179666,
            -10.6339907853,
            -10.5495538639,
            -10.4546677714,
            -10.3506386542,
            -10.2386366017,
            -10.1197772808,
            -9.99504030111,
            -9.86535084973,
            -9.73155247952,
        ]
        num_eos = EOS(eos_name="numerical_eos")
        self.num_eos_fit = num_eos.fit(self.volumes, self.energies)

    def test_run_all_models(self):
        # these have been checked for plausibility,
        # but are not benchmarked against independently known values
        test_output = {
            "birch": {
                "b0": 0.5369258244952931,
                "b1": 4.178644231838501,
                "e0": -10.8428039082307,
                "v0": 40.98926572870838,
            },
            "birch_murnaghan": {
                "b0": 0.5369258245417454,
                "b1": 4.178644235500821,
                "e0": -10.842803908240892,
                "v0": 40.98926572528106,
            },
            "deltafactor": {
                "b0": 0.5369258245611414,
                "b1": 4.178644231924639,
                "e0": -10.842803908299294,
                "v0": 40.989265727927936,
            },
            "murnaghan": {
                "b0": 0.5144967693786603,
                "b1": 3.9123862262572264,
                "e0": -10.836794514626673,
                "v0": 41.13757930387086,
            },
            "numerical_eos": {
                "b0": 0.5557257614101998,
                "b1": 4.344039148405489,
                "e0": -10.847490826530702,
                "v0": 40.857200064982536,
            },
            "pourier_tarantola": {
                "b0": 0.5667729960804602,
                "b1": 4.331688936974368,
                "e0": -10.851486685041658,
                "v0": 40.86770643373908,
            },
            "vinet": {
                "b0": 0.5493839425156859,
                "b1": 4.3051929654936885,
                "e0": -10.846160810560756,
                "v0": 40.916875663779784,
            },
        }

        for eos_name in EOS.MODELS:
            eos = EOS(eos_name=eos_name)
            _ = eos.fit(self.volumes, self.energies)
            for param in ("b0", "b1", "e0", "b0"):
                # TODO: solutions only stable to 2 decimal places
                # between different machines, this seems far too low?
                self.assertArrayAlmostEqual(_.results[param], test_output[eos_name][param], decimal=1)

    def test_fitting(self):
        # courtesy of @katherinelatimer2013
        # known correct values for Vinet

        # Mg

        mp153_volumes = [
            16.69182365,
            17.25441763,
            17.82951915,
            30.47573817,
            18.41725977,
            29.65211363,
            28.84346369,
            19.01777055,
            28.04965916,
            19.63120886,
            27.27053682,
            26.5059864,
            20.25769112,
            25.75586879,
            20.89736201,
            25.02003097,
            21.55035204,
            24.29834347,
            22.21681221,
            23.59066888,
            22.89687316,
        ]

        mp153_energies = [
            -1.269884575,
            -1.339411225,
            -1.39879471,
            -1.424480995,
            -1.44884184,
            -1.45297499,
            -1.4796246,
            -1.49033594,
            -1.504198485,
            -1.52397006,
            -1.5264432,
            -1.54609291,
            -1.550269435,
            -1.56284009,
            -1.569937375,
            -1.576420935,
            -1.583470925,
            -1.58647189,
            -1.591436505,
            -1.592563495,
            -1.594347355,
        ]

        mp153_known_energies_vinet = [
            -1.270038831,
            -1.339366487,
            -1.398683238,
            -1.424556061,
            -1.448746649,
            -1.453000456,
            -1.479614511,
            -1.490266797,
            -1.504163502,
            -1.523910268,
            -1.526395734,
            -1.546038792,
            -1.550298657,
            -1.562800797,
            -1.570015274,
            -1.576368392,
            -1.583605186,
            -1.586404575,
            -1.591578378,
            -1.592547954,
            -1.594410995,
        ]

        # C: 4.590843262
        # B: 2.031381599
        mp153_known_e0_vinet = -1.594429229
        mp153_known_v0_vinet = 22.95764159

        eos = EOS(eos_name="vinet")

        fit = eos.fit(mp153_volumes, mp153_energies)

        np.testing.assert_array_almost_equal(fit.func(mp153_volumes), mp153_known_energies_vinet, decimal=5)

        assert mp153_known_e0_vinet == approx(fit.e0, abs=1e-4)
        assert mp153_known_v0_vinet == approx(fit.v0, abs=1e-4)

        # expt. value 35.5, known fit 36.16
        assert fit.b0_GPa == approx(36.16258687442761, abs=1e-4)

        # Si

        mp149_volumes = [
            15.40611854,
            14.90378698,
            16.44439516,
            21.0636307,
            17.52829835,
            16.98058208,
            18.08767363,
            18.65882487,
            19.83693435,
            15.91961152,
            22.33987173,
            21.69548924,
            22.99688883,
            23.66666322,
            20.44414922,
            25.75374305,
            19.24187473,
            24.34931029,
            25.04496106,
            27.21116571,
            26.4757653,
        ]

        mp149_energies = [
            -4.866909695,
            -4.7120965,
            -5.10921253,
            -5.42036228,
            -5.27448405,
            -5.200810795,
            -5.331915665,
            -5.3744186,
            -5.420058145,
            -4.99862686,
            -5.3836163,
            -5.40610838,
            -5.353700425,
            -5.31714654,
            -5.425263555,
            -5.174988295,
            -5.403353105,
            -5.27481447,
            -5.227210275,
            -5.058992615,
            -5.118805775,
        ]

        mp149_known_energies_vinet = [
            -4.866834585,
            -4.711786499,
            -5.109642598,
            -5.420093739,
            -5.274605844,
            -5.201025714,
            -5.331899365,
            -5.374315789,
            -5.419671568,
            -4.998827503,
            -5.383703409,
            -5.406038887,
            -5.353926272,
            -5.317484252,
            -5.424963418,
            -5.175090887,
            -5.403166824,
            -5.275096644,
            -5.227427635,
            -5.058639193,
            -5.118654229,
        ]

        # C: 4.986513158
        # B: 4.964976215
        mp149_known_e0_vinet = -5.424963506
        mp149_known_v0_vinet = 20.44670279

        eos = EOS(eos_name="vinet")

        fit = eos.fit(mp149_volumes, mp149_energies)

        np.testing.assert_array_almost_equal(fit.func(mp149_volumes), mp149_known_energies_vinet, decimal=5)

        assert mp149_known_e0_vinet == approx(fit.e0, abs=1e-4)
        assert mp149_known_v0_vinet == approx(fit.v0, abs=1e-4)

        # expt. value 97.9, known fit 88.39
        assert fit.b0_GPa == approx(88.38629337404822, abs=1e-4)

        # Ti

        mp72_volumes = [
            12.49233296,
            12.91339188,
            13.34380224,
            22.80836212,
            22.19195533,
            13.78367177,
            21.58675559,
            14.23310328,
            20.99266009,
            20.4095592,
            14.69220297,
            19.83736385,
            15.16106697,
            19.2759643,
            15.63980711,
            18.72525771,
            16.12851491,
            18.18514127,
            16.62729878,
            17.65550599,
            17.13626153,
        ]

        mp72_energies = [
            -7.189983803,
            -7.33985647,
            -7.468745423,
            -7.47892835,
            -7.54945107,
            -7.578012237,
            -7.61513166,
            -7.66891898,
            -7.67549721,
            -7.73000681,
            -7.74290386,
            -7.77803379,
            -7.801246383,
            -7.818964483,
            -7.84488189,
            -7.85211192,
            -7.87486651,
            -7.876767777,
            -7.892161533,
            -7.892199957,
            -7.897605303,
        ]

        mp72_known_energies_vinet = [
            -7.189911138,
            -7.339810181,
            -7.468716095,
            -7.478678021,
            -7.549402394,
            -7.578034391,
            -7.615240977,
            -7.669091347,
            -7.675683891,
            -7.730188653,
            -7.74314028,
            -7.778175824,
            -7.801363213,
            -7.819030923,
            -7.844878053,
            -7.852099741,
            -7.874737806,
            -7.876686864,
            -7.891937429,
            -7.892053535,
            -7.897414664,
        ]

        # C: 3.958192998
        # B: 6.326790098
        mp72_known_e0_vinet = -7.897414997
        mp72_known_v0_vinet = 17.13223229

        eos = EOS(eos_name="vinet")

        fit = eos.fit(mp72_volumes, mp72_energies)

        np.testing.assert_array_almost_equal(fit.func(mp72_volumes), mp72_known_energies_vinet, decimal=5)

        assert mp72_known_e0_vinet == approx(fit.e0, abs=1e-4)
        assert mp72_known_v0_vinet == approx(fit.v0, abs=1e-4)

        # expt. value 107.3, known fit 112.63
        assert fit.b0_GPa == approx(112.62927187296167, abs=1e-4)

    def test_numerical_eoswrapper(self):
        # using numerical eos directly vs via EOS wrapper
        numerical_eos = NumericalEOS(self.volumes, self.energies)
        numerical_eos.fit()
        assert len(numerical_eos.eos_params) > 3
        assert float(numerical_eos.e0) == approx(self.num_eos_fit.e0, abs=1e-3)
        assert float(numerical_eos.v0) == approx(self.num_eos_fit.v0, abs=1e-3)
        assert float(numerical_eos.b0) == approx(self.num_eos_fit.b0, abs=1e-3)
        assert float(numerical_eos.b1) == approx(self.num_eos_fit.b1, abs=1e-3)
        self.assertArrayAlmostEqual(numerical_eos.eos_params, self.num_eos_fit.eos_params)

    def test_numerical_eos_values(self):
        np.testing.assert_almost_equal(self.num_eos_fit.e0, -10.84749, decimal=3)
        np.testing.assert_almost_equal(self.num_eos_fit.v0, 40.857201, decimal=1)
        np.testing.assert_almost_equal(self.num_eos_fit.b0, 0.55, decimal=2)
        # TODO: why were these tests commented out?
        # np.testing.assert_almost_equal(self.num_eos_fit.b0_GPa, 89.0370727, decimal=1)
        # np.testing.assert_almost_equal(self.num_eos_fit.b1, 4.344039, decimal=2)

    def test_eos_func(self):
        # list vs np.array arguments
        np.testing.assert_almost_equal(
            self.num_eos_fit.func([0, 1, 2]),
            self.num_eos_fit.func(np.array([0, 1, 2])),
            decimal=10,
        )
        # func vs _func
        np.testing.assert_almost_equal(
            self.num_eos_fit.func(0.0),
            self.num_eos_fit._func(0.0, self.num_eos_fit.eos_params),
            decimal=10,
        )
        # test the eos function: energy = f(volume)
        # numerical eos evaluated at volume=0 == a0 of the fit polynomial
        np.testing.assert_almost_equal(self.num_eos_fit.func(0.0), self.num_eos_fit.eos_params[-1], decimal=6)
        birch_eos = EOS(eos_name="birch")
        birch_eos_fit = birch_eos.fit(self.volumes, self.energies)
        # birch eos evaluated at v0 == e0
        np.testing.assert_almost_equal(birch_eos_fit.func(birch_eos_fit.v0), birch_eos_fit.e0, decimal=6)

        # TODO: Reactivate
        # fig = birch_eos_fit.plot_ax(ax=None, show=False, fontsize=8, title="birch eos")
        # self.assertTrue(hasattr(fig, "savefig"))

    def test_eos_func_call(self):
        # eos_fit_obj.func(volume) == eos_fit_obj(volume)
        np.testing.assert_almost_equal(self.num_eos_fit.func(0.0), self.num_eos_fit(0.0), decimal=10)

    def test_summary_dict(self):
        d = {
            "e0": self.num_eos_fit.e0,
            "b0": self.num_eos_fit.b0,
            "b1": self.num_eos_fit.b1,
            "v0": self.num_eos_fit.v0,
        }
        assert self.num_eos_fit.results == d


if __name__ == "__main__":
    unittest.main()
