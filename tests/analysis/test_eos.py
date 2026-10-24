from __future__ import annotations

import warnings

import matplotlib.pyplot as plt
import numpy as np
import pytest
from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.analysis.eos import (
    EOS,
    Birch,
    BirchMurnaghan,
    EOSBase,
    EOSError,
    Murnaghan,
    NumericalEOS,
    PolynomialEOS,
    PourierTarantola,
    Vinet,
)
from pymatgen.util.testing import MatSciTest


class TestEOS(MatSciTest):
    def setup_method(self):
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
                assert_allclose(_.results[param], test_output[eos_name][param], atol=1e-1)

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

        assert_allclose(fit.func(mp153_volumes), mp153_known_energies_vinet, atol=1e-5)

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
            *[
                -4.866834585,
                -4.711786499,
                -5.109642598,
                -5.420093739,
                -5.274605844,
                -5.201025714,
            ],
            *[
                -5.331899365,
                -5.374315789,
                -5.419671568,
                -4.998827503,
                -5.383703409,
                -5.406038887,
            ],
            *[
                -5.353926272,
                -5.317484252,
                -5.424963418,
                -5.175090887,
                -5.403166824,
                -5.275096644,
            ],
            *[-5.227427635, -5.058639193, -5.118654229],
        ]

        # C: 4.986513158
        # B: 4.964976215
        mp149_known_e0_vinet = -5.424963506
        mp149_known_v0_vinet = 20.44670279

        eos = EOS(eos_name="vinet")

        fit = eos.fit(mp149_volumes, mp149_energies)

        assert_allclose(fit.func(mp149_volumes), mp149_known_energies_vinet, atol=1e-5)

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

        assert_allclose(fit.func(mp72_volumes), mp72_known_energies_vinet, atol=1e-5)

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
        assert_allclose(numerical_eos.eos_params, self.num_eos_fit.eos_params)

    def test_numerical_eos_values(self):
        assert_allclose(self.num_eos_fit.e0, -10.84749, atol=1e-3)
        assert_allclose(self.num_eos_fit.v0, 40.857201, atol=1e-1)
        assert_allclose(self.num_eos_fit.b0, 0.55, atol=1e-2)
        assert_allclose(self.num_eos_fit.b0_GPa, 89.0370727, atol=1)
        assert_allclose(self.num_eos_fit.b1, 4.344039, atol=1)

    def test_eos_func(self):
        # list vs np.array arguments
        assert_allclose(
            self.num_eos_fit.func([0, 1, 2]),
            self.num_eos_fit.func(np.array([0, 1, 2])),
        )
        # func vs _func
        assert_allclose(
            self.num_eos_fit.func(0.0),
            self.num_eos_fit._func(0.0, self.num_eos_fit.eos_params),
        )
        # test the eos function: energy = f(volume)
        # numerical eos evaluated at volume=0 == a0 of the fit polynomial
        assert_allclose(self.num_eos_fit.func(0.0), self.num_eos_fit.eos_params[-1])
        birch_eos = EOS(eos_name="birch")
        birch_eos_fit = birch_eos.fit(self.volumes, self.energies)
        # birch eos evaluated at v0 == e0
        assert_allclose(birch_eos_fit.func(birch_eos_fit.v0), birch_eos_fit.e0)

        fig = birch_eos_fit.plot_ax(ax=None, show=False, fontsize=8, title="birch eos")
        assert callable(fig.savefig)

    def test_eos_func_call(self):
        # eos_fit_obj.func(volume) == eos_fit_obj(volume)
        assert_allclose(self.num_eos_fit.func(0.0), self.num_eos_fit(0.0))

    def test_summary_dict(self):
        assert self.num_eos_fit.results == {
            "e0": self.num_eos_fit.e0,
            "b0": self.num_eos_fit.b0,
            "b1": self.num_eos_fit.b1,
            "v0": self.num_eos_fit.v0,
        }


# Reference parameters (e0 [eV], b0 [eV/Ang^3], b1, v0 [Ang^3]) used to generate synthetic E(V) data.
REF_PARAMS = (-10.0, 0.6, 4.5, 40.0)
REF_RESULTS = dict(zip(("e0", "b0", "b1", "v0"), REF_PARAMS, strict=True))
REF_VOLUMES = np.linspace(34, 47, 21)

# deltafactor is a cubic in V^(-2/3), i.e. exactly a 3rd order Birch-Murnaghan EOS.
EXACT_MODELS = {
    "murnaghan": Murnaghan,
    "birch": Birch,
    "birch_murnaghan": BirchMurnaghan,
    "pourier_tarantola": PourierTarantola,
    "vinet": Vinet,
    "deltafactor": BirchMurnaghan,
}


def _synthetic_energies(model: type[EOSBase], volumes=REF_VOLUMES, params=REF_PARAMS) -> np.ndarray:
    volumes = np.asarray(volumes)
    return model(volumes, volumes)._func(volumes, params)


@pytest.mark.parametrize(("eos_name", "generator"), EXACT_MODELS.items())
def test_fit_recovers_known_parameters(eos_name, generator):
    energies = _synthetic_energies(generator)
    fit = EOS(eos_name=eos_name).fit(REF_VOLUMES, energies)

    assert isinstance(fit, EOS.MODELS[eos_name])
    assert fit.results == approx(REF_RESULTS, rel=1e-6)
    # 1 eV/Ang^3 = 160.21766208 GPa
    assert float(fit.b0_GPa) == approx(REF_PARAMS[1] * 160.21766208, rel=1e-6)
    assert str(fit.b0_GPa.unit) == "GPa"
    assert_allclose(fit(REF_VOLUMES), energies, atol=1e-8)
    assert fit.func(REF_PARAMS[3]) == approx(REF_PARAMS[0], abs=1e-8)


def test_numerical_eos_recovers_known_parameters():
    energies = _synthetic_energies(BirchMurnaghan)
    fit = EOS(eos_name="numerical_eos").fit(REF_VOLUMES, energies)
    assert isinstance(fit, NumericalEOS)
    # the weighted-average polynomial is padded to n_data coefficients
    assert len(fit.eos_params) == len(REF_VOLUMES)
    assert fit.results == approx(REF_RESULTS, rel=1e-3)
    assert_allclose(fit(REF_VOLUMES), energies, atol=1e-6)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"min_ndata_factor": 1},
        {"min_ndata_factor": 6, "max_poly_order_factor": 8},
        {"max_poly_order_factor": 2},
        {"min_poly_order": 3},
    ],
)
def test_numerical_eos_fit_options(kwargs):
    fit = NumericalEOS(REF_VOLUMES, _synthetic_energies(BirchMurnaghan))
    fit.fit(**kwargs)
    assert fit.results == approx(REF_RESULTS, rel=1e-3)


def test_numerical_eos_quadratic_only():
    # 7 points with the default factors: only quadratic polynomials are candidates,
    # for which b1 = -V0 * dB/dV / B0 = -1 exactly.
    volumes = np.linspace(36, 44, 7)
    fit = NumericalEOS(volumes, _synthetic_energies(BirchMurnaghan, volumes))
    fit.fit()
    assert fit.b1 == approx(-1)
    assert fit.v0 == approx(40, abs=0.5)


def test_numerical_eos_does_not_leak_warning_filters():
    filters_before = list(warnings.filters)
    NumericalEOS(REF_VOLUMES, _synthetic_energies(BirchMurnaghan)).fit()
    assert warnings.filters == filters_before


@pytest.mark.parametrize(
    ("volumes", "match"),
    [
        (np.linspace(36, 44, 5), "No polynomial fit"),
        (np.linspace(36, 44, 6), "No polynomial fit"),
        (np.linspace(40, 50, 11), "boundary"),
        (np.linspace(30, 40, 11), "boundary"),
    ],
)
def test_numerical_eos_errors(volumes, match):
    fit = NumericalEOS(volumes, _synthetic_energies(BirchMurnaghan, volumes))
    with pytest.raises(EOSError, match=match):
        fit.fit()


def test_polynomial_eos_fit():
    energies = 0.01 * (REF_VOLUMES - 40.0) ** 2 - 5.0
    fit = PolynomialEOS(REF_VOLUMES, energies)
    fit.fit(order=2)
    assert_allclose(fit.eos_params, [0.01, -0.8, 11.0], atol=1e-10)
    # B0 = V0 * d2E/dV2 = 40 * 0.02; b1 of a parabola is exactly -1
    assert fit.results == approx({"e0": -5.0, "b0": 0.8, "b1": -1.0, "v0": 40.0}, rel=1e-6)


def test_deltafactor_no_minimum():
    x = REF_VOLUMES ** (-2 / 3)
    x_mid = x.mean()
    # E(x) = -(x^3 - 3 x_mid^2 x) has its maximum at x_mid and its minimum at -x_mid < 0
    with pytest.raises(EOSError, match="No minimum could be found"):
        EOS(eos_name="deltafactor").fit(REF_VOLUMES, -(x**3 - 3 * x_mid**2 * x))


def test_initial_guess_minimum_below_volume_range():
    volumes = np.linspace(42, 50, 9)
    with pytest.raises(EOSError, match="minimum volume of a fitted parabola"):
        EOS(eos_name="vinet").fit(volumes, _synthetic_energies(Vinet, volumes))


@pytest.mark.xfail(
    reason="EOSBase._initial_guess checks `not vmin < v0 and v0 < vmax`, which never raises when the "
    "parabola minimum lies above the largest input volume.",
    strict=True,
)
def test_initial_guess_minimum_above_volume_range():
    volumes = np.linspace(30, 38, 9)
    with pytest.raises(EOSError, match="minimum volume of a fitted parabola"):
        EOS(eos_name="vinet").fit(volumes, _synthetic_energies(Vinet, volumes))


def test_leastsq_failure(monkeypatch):
    energies = _synthetic_energies(Vinet)
    # ier=5: maximum number of function evaluations exceeded
    monkeypatch.setattr("pymatgen.analysis.eos.leastsq", lambda *args, **kwargs: (np.array(REF_PARAMS), 5))
    with pytest.raises(EOSError, match="Optimal parameters not found"):
        EOS(eos_name="vinet").fit(REF_VOLUMES, energies)


@pytest.mark.parametrize("attr", ["e0", "b0", "b0_GPa"])
def test_params_before_fit(attr):
    eos = Vinet(REF_VOLUMES, _synthetic_energies(Vinet))
    with pytest.raises(RuntimeError, match="params have not be initialized"):
        getattr(eos, attr)


def test_invalid_eos_name():
    with pytest.raises(EOSError, match="'foo' is not supported"):
        EOS(eos_name="foo")


@pytest.fixture
def vinet_fit():
    return EOS(eos_name="vinet").fit(REF_VOLUMES, _synthetic_energies(Vinet))


def test_plot(vinet_fit):
    ax = vinet_fit.plot(width=6, height=4, dpi=50)
    try:
        assert ax.get_xlabel() == "Volume $\\AA^3$"
        assert ax.get_ylabel() == "Energy (eV)"
        data_line, fit_line = ax.get_lines()
        assert_allclose(data_line.get_xdata(), REF_VOLUMES)
        assert_allclose(data_line.get_ydata(), vinet_fit.energies)
        assert data_line.get_marker() == "o"
        assert fit_line.get_label() == "Vinet fit"
        assert fit_line.get_linestyle() == "--"
        assert fit_line.get_color() == "r"
        v_fit = fit_line.get_xdata()
        assert len(v_fit) == 100
        assert v_fit[0] == approx(34 * 0.99)
        assert v_fit[-1] == approx(47 * 1.01)
        assert_allclose(fit_line.get_ydata(), vinet_fit.func(v_fit))
        text = ax.texts[0].get_text()
        assert "Equation of State: Vinet" in text
        assert "Minimum energy = -10.00 eV" in text
        assert "Minimum or reference volume = 40.00 Ang^3" in text
        assert "Bulk modulus = 0.60 eV/Ang^3 = 96.13 GPa" in text
        assert "Derivative of bulk modulus w.r.t. pressure = 4.50" in text
    finally:
        plt.close(ax.figure)


def test_plot_custom_kwargs(vinet_fit):
    fig, ax = plt.subplots()
    try:
        assert vinet_fit.plot(ax=ax, color="b", label="my fit", text="custom text") is ax
        assert [line.get_color() for line in ax.get_lines()] == ["b", "b"]
        assert ax.get_lines()[1].get_label() == "my fit"
        assert ax.texts[0].get_text() == "custom text"
    finally:
        plt.close(fig)


def test_plot_ax(vinet_fit):
    fig_in, ax = plt.subplots()
    try:
        fig = vinet_fit.plot_ax(ax=ax, fontsize=9, show=False, color="g", label="vinet")
        assert fig is fig_in
        data_line, fit_line = ax.get_lines()
        assert_allclose(data_line.get_xdata(), REF_VOLUMES)
        assert fit_line.get_label() == "vinet"
        assert fit_line.get_color() == "g"
        assert_allclose(fit_line.get_ydata(), vinet_fit.func(fit_line.get_xdata()))
        assert ax.texts[0].get_fontsize() == 9
        assert "Bulk modulus = 0.60 eV/Ang^3" in ax.texts[0].get_text()
    finally:
        plt.close(fig_in)
