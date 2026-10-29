from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pytest import approx
from scipy.constants import physical_constants

from pymatgen.analysis.eos import EOS, BirchMurnaghan, EOSError, Vinet
from pymatgen.analysis.quasiharmonic import QuasiHarmonicDebyeApprox
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure

__author__ = "Kiran Mathew"


class TestQuasiHarmonicDebyeApprox:
    def setup_method(self):
        struct = Structure.from_dict(
            {
                "lattice": {
                    "a": 2.5630200477817295,
                    "alpha": 59.999993839702206,
                    "b": 2.563020442699644,
                    "beta": 59.999988742674944,
                    "c": 2.56301993,
                    "gamma": 60.00000504373715,
                    "matrix": [
                        [2.21964022, 0, 1.28151046],
                        [0.73987974, 2.09269747, 1.28151046],
                        [-0, -0, 2.56301993],
                    ],
                    "volume": 11.905318492097948,
                },
                "sites": [
                    {
                        "abc": [-0, -0, -0],
                        "label": "B",
                        "species": [{"element": "B", "occu": 1}],
                        "xyz": [0, 0, 0],
                    },
                    {
                        "abc": [0.25, 0.25, 0.25],
                        "label": "N",
                        "species": [{"element": "N", "occu": 1}],
                        "xyz": [0.73987999, 0.5231743675, 1.2815102125],
                    },
                ],
            }
        )

        self.energies = [
            -15.76315281,
            -16.11541813,
            -16.41784171,
            -16.47471523,
            -16.63624155,
            -16.6741551,
            -16.78661144,
            -16.88768073,
            -16.92450672,
            -17.04863261,
            -17.06126553,
            -17.15786866,
            -17.19784976,
            -17.25078749,
            -17.30017149,
            -17.32578594,
            -17.3708922,
            -17.38125127,
            -17.41231934,
            -17.41534352,
            -17.42636644,
        ]
        self.volumes = [
            8.678977833994137,
            8.971505437031707,
            9.27052889309282,
            15.845976281427582,
            15.417733609491387,
            9.576127994353376,
            14.997270631725604,
            9.888370962140854,
            14.584523227465766,
            14.179424329180256,
            10.20732378093211,
            13.78189117535765,
            10.533067462993838,
            13.391864274742145,
            10.865663655755416,
            13.009260480347871,
            11.205193091129587,
            12.634015019827533,
            11.551718049704352,
            12.26606042141808,
            11.90531496343142,
        ]
        self.eos = "vinet"
        self.T = 300
        self.qhda = QuasiHarmonicDebyeApprox(
            self.energies,
            self.volumes,
            struct,
            t_min=self.T,
            t_max=self.T,
            eos=self.eos,
        )
        self.opt_vol = 11.957803302392925

    def test_bulk_modulus(self):
        eos = EOS(self.eos)
        eos_fit = eos.fit(self.volumes, self.energies)
        bulk_modulus = float(str(eos_fit.b0_GPa).split()[0])
        bulk_modulus_ans = float(str(self.qhda.bulk_modulus).split()[0])
        assert_allclose(bulk_modulus, bulk_modulus_ans, atol=1e-3)

    def test_optimum_volume(self):
        opt_vol = self.qhda.optimum_volumes[0]
        assert_allclose(opt_vol, self.opt_vol, atol=1e-3)

    def test_debye_temperature(self):
        theta = self.qhda.debye_temperature(self.opt_vol)
        assert_allclose(theta, 2559.675227, atol=1e-3)

    def test_gruneisen_parameter(self):
        gamma = self.qhda.gruneisen_parameter(self.T, self.opt_vol)
        assert_allclose(gamma, 1.670486, atol=1e-3)

    def test_thermal_conductivity(self):
        kappa = self.qhda.thermal_conductivity(self.T, self.opt_vol)
        assert_allclose(kappa, 131.736242, atol=1e-1)

    def test_vibrational_internal_energy(self):
        u = self.qhda.vibrational_internal_energy(self.T, self.opt_vol)
        assert_allclose(u, 0.50102, atol=1e-3)

    def test_vibrational_free_energy(self):
        A = self.qhda.vibrational_free_energy(self.T, self.opt_vol)
        assert_allclose(A, 0.494687, atol=1e-3)


class TestAnharmonicQuasiHarmonicDebyeApprox:
    def setup_method(self):
        struct = Structure.from_str(
            """FCC Al
1.0
2.473329 0.000000 1.427977
0.824443 2.331877 1.427977
0.000000 0.000000 2.855955
Al
1
direct
0.000000 0.000000 0.000000 Al""",
            fmt="POSCAR",
        )

        self.energies = [
            -3.69150886,
            -3.70788383,
            -3.71997361,
            -3.72522301,
            -3.73569569,
            -3.73649743,
            -3.74054982,
        ]
        self.volumes = [
            14.824542034870653,
            18.118887714656875,
            15.373596786943025,
            17.569833126580278,
            15.92265868064787,
            17.02077912220064,
            16.471717630914863,
        ]
        self.eos = "vinet"
        self.T = 500
        self.qhda = QuasiHarmonicDebyeApprox(
            self.energies,
            self.volumes,
            struct,
            t_min=self.T,
            t_max=self.T,
            eos=self.eos,
            anharmonic_contribution=True,
        )
        self.opt_vol = 17.216094889116807

    def test_optimum_volume(self):
        opt_vol = self.qhda.optimum_volumes[0]
        assert_allclose(opt_vol, self.opt_vol, atol=1e-3)

    def test_debye_temperature(self):
        theta = self.qhda.debye_temperature(self.opt_vol)
        assert_allclose(theta, 601.239096, 4)

    def test_gruneisen_parameter(self):
        gamma = self.qhda.gruneisen_parameter(0, self.qhda.ev_eos_fit.v0)
        assert_allclose(gamma, 2.188302, atol=1e-2)

    def test_thermal_conductivity(self):
        kappa = self.qhda.thermal_conductivity(self.T, self.opt_vol)
        assert_allclose(kappa, 21.810997, atol=1e-1)

    def test_vibrational_internal_energy(self):
        u = self.qhda.vibrational_internal_energy(self.T, self.opt_vol)
        assert_allclose(u, 0.13845, atol=1e-3)

    def test_vibrational_free_energy(self):
        A = self.qhda.vibrational_free_energy(self.T, self.opt_vol)
        assert_allclose(A, -0.014620, atol=1e-3)


# Synthetic, Al-like 4-atom fcc cell with E(V) generated from known EOS parameters:
# e0 [eV], b0 [eV/Ang^3] (= 75.3 GPa), b1, v0 [Ang^3].
KB = physical_constants["Boltzmann constant in eV/K"][0]
GPA_TO_EV_ANG3 = 1 / 160.21766208
SYN_PARAMS = (-14.9, 0.47, 4.6, 66.0)
SYN_VOLUMES = np.linspace(58, 76, 13)
SYN_STRUCT = Structure(Lattice.cubic(4.05), ["Al"] * 4, [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
SYN_NATOMS = 4


def _energies(model=Vinet):
    return model(SYN_VOLUMES, SYN_VOLUMES)._func(SYN_VOLUMES, SYN_PARAMS)


def _qhda(energies=None, **kwargs):
    kwargs = {"t_min": 300, "t_max": 300, **kwargs}
    return QuasiHarmonicDebyeApprox(_energies() if energies is None else energies, SYN_VOLUMES, SYN_STRUCT, **kwargs)


@pytest.fixture(scope="module")
def qhda_vinet():
    return _qhda()


def test_mie_gruneisen_and_anharmonic_are_exclusive():
    with pytest.raises(ValueError, match="cannot be used together"):
        _qhda(use_mie_gruneisen=True, anharmonic_contribution=True)


@pytest.mark.parametrize(
    ("y", "expected"),
    [
        (1e-4, 1 - 3e-4 / 8),  # D(y) = 1 - 3y/8 + O(y^2)
        (1.0, 0.6744155640778),
        # beyond y=155 the integral is replaced by its limit pi^4/15
        (155.0, np.pi**4 / 5 / 155.0**3),
        (1000.0, np.pi**4 / 5 / 1000.0**3),
    ],
)
def test_debye_integral(y, expected):
    assert QuasiHarmonicDebyeApprox.debye_integral(y) == approx(expected, rel=1e-6)


def test_debye_integral_continuous_at_cutoff():
    assert QuasiHarmonicDebyeApprox.debye_integral(154.999) == approx(
        QuasiHarmonicDebyeApprox.debye_integral(155.0), rel=1e-4
    )


def test_vibrational_energies_classical_limit(qhda_vinet):
    # T >> theta_D: U_vib -> 3 n kB T and A_vib -> n kB T (3 ln(theta_D / T) - 1)
    temp, vol = 1e5, 66.0
    theta = qhda_vinet.debye_temperature(vol)
    assert qhda_vinet.vibrational_internal_energy(temp, vol) == approx(3 * SYN_NATOMS * KB * temp, rel=1e-5)
    assert qhda_vinet.vibrational_free_energy(temp, vol) == approx(
        SYN_NATOMS * KB * temp * (3 * np.log(theta / temp) - 1), rel=1e-5
    )


def test_vibrational_energies_zero_point_limit(qhda_vinet):
    # T << theta_D (theta_D / T > 155, the asymptotic branch): U_vib = A_vib = 9/8 n kB theta_D
    vol = 66.0
    zero_point = 9 / 8 * SYN_NATOMS * KB * qhda_vinet.debye_temperature(vol)
    assert qhda_vinet.vibrational_internal_energy(1.0, vol) == approx(zero_point, rel=1e-8)
    assert qhda_vinet.vibrational_free_energy(1.0, vol) == approx(zero_point, rel=1e-8)


def test_vibrational_energies_thermodynamic_consistency(qhda_vinet):
    # U = A - T dA/dT
    temp, vol, dt = 300.0, 66.0, 1e-2
    dadt = (qhda_vinet.vibrational_free_energy(temp + dt, vol) - qhda_vinet.vibrational_free_energy(temp - dt, vol)) / (
        2 * dt
    )
    assert qhda_vinet.vibrational_free_energy(temp, vol) - temp * dadt == approx(
        qhda_vinet.vibrational_internal_energy(temp, vol), rel=1e-6
    )


def test_debye_temperature_scaling(qhda_vinet):
    # harmonic Debye temperature at fixed bulk modulus scales as V^(-1/6)
    assert qhda_vinet.debye_temperature(2**6 * 66.0) == approx(qhda_vinet.debye_temperature(66.0) / 2)
    # anharmonic: theta(V) = theta_harmonic(V) * (V0 / V)^gamma(0 K, V0)
    qhda_anh = _qhda(anharmonic_contribution=True)
    v0 = qhda_anh.ev_eos_fit.v0
    gamma = qhda_anh.gruneisen_parameter(0, v0)
    assert qhda_anh.debye_temperature(70.0) == approx(qhda_vinet.debye_temperature(70.0) * (v0 / 70.0) ** gamma)
    assert qhda_anh.optimum_volumes[0] > qhda_vinet.optimum_volumes[0]


def test_thermal_expansion():
    qhda = _qhda(t_min=100, t_max=900, t_step=200)
    assert qhda.temperatures == approx([100, 300, 500, 700, 900])
    assert np.all(np.diff(qhda.optimum_volumes) > 0)
    assert np.all(np.diff(qhda.gibbs_free_energy) < 0)
    # zero-point expansion pushes V_opt above V0 even at low T
    assert qhda.optimum_volumes[0] > SYN_PARAMS[3]


def test_pressure(qhda_vinet):
    qhda = _qhda(pressure=1.0)
    v_0, v_p = qhda_vinet.optimum_volumes[0], qhda.optimum_volumes[0]
    # dV/V ~ -P/B for small P
    assert (v_p - v_0) / v_0 == approx(-1.0 / float(qhda.bulk_modulus), rel=0.1)
    # G includes the PV term
    assert qhda.gibbs_free_energy[0] > qhda_vinet.gibbs_free_energy[0]


def test_failed_gibbs_minimization():
    # at 10 GPa the G(V) minimum is below the sampled volumes and the EOS fit fails
    with pytest.raises(EOSError, match="minimum volume of a fitted parabola"):
        _qhda(pressure=10.0)
    # with several temperatures, failures are skipped and stored as None
    qhda = _qhda(pressure=10.0, t_min=300, t_max=500, t_step=200)
    assert qhda.temperatures == approx([300, 500])
    assert qhda.gibbs_free_energy == [None, None]
    assert qhda.optimum_volumes == [None, None]


@pytest.mark.parametrize(
    ("eos", "generator"),
    [
        ("murnaghan", None),
        ("birch", None),
        ("birch_murnaghan", None),
        ("pourier_tarantola", None),
        ("vinet", None),
        ("deltafactor", BirchMurnaghan),
    ],
)
def test_eos_options(eos, generator):
    generator = generator or EOS.MODELS[eos]
    qhda = _qhda(_energies(generator), eos=eos)
    _e0, b0, b1, v0 = SYN_PARAMS
    assert float(qhda.bulk_modulus) == approx(b0 / GPA_TO_EV_ANG3, rel=1e-6)
    assert qhda.ev_eos_fit.v0 == approx(v0)
    # Slater gamma at V0: -(1/6 + V0/(2 B0) dB/dV) = B1/2 - 1/6
    assert qhda.gruneisen_parameter(300, v0) == approx(b1 / 2 - 1 / 6, rel=1e-5)
    # Slater gamma does not depend on T
    assert qhda.gruneisen_parameter(1000, v0) == qhda.gruneisen_parameter(300, v0)
    # zero-point + thermal expansion at 300 K: V0 < V_opt < 1.01 V0
    assert v0 < qhda.optimum_volumes[0] < 1.01 * v0


def test_gruneisen_numerical_eos():
    energies = _energies(BirchMurnaghan)
    qhda_num = _qhda(energies, eos="numerical_eos")
    qhda_bm = _qhda(energies, eos="birch_murnaghan")
    volumes = [62.0, 64.0, 66.0, 68.0, 70.0]
    assert [qhda_num.gruneisen_parameter(300, v) for v in volumes] == approx(
        [qhda_bm.gruneisen_parameter(300, v) for v in volumes], rel=0.05
    )


def _vinet_pressure(volume):
    # analytic Vinet P(V) = -dE/dV in eV/Ang^3
    _e0, b0, b1, v0 = SYN_PARAMS
    eta = (volume / v0) ** (1 / 3)
    return 3 * b0 * (1 - eta) / eta**2 * np.exp(1.5 * (b1 - 1) * (1 - eta))


@pytest.mark.parametrize("pressure", [0.0, 1.0])
def test_mie_gruneisen(pressure):
    qhda = _qhda(use_mie_gruneisen=True, pressure=pressure)
    temp, vol = 300.0, 67.0
    # gamma = V (P + dE/dV) / U_vib
    expected = vol * (pressure * GPA_TO_EV_ANG3 - _vinet_pressure(vol)) / qhda.vibrational_internal_energy(temp, vol)
    assert qhda.gruneisen_parameter(temp, vol) == approx(expected, rel=1e-5)
    # unlike the Slater gamma, the Mie-Gruneisen gamma is temperature dependent
    assert qhda.gruneisen_parameter(600.0, vol) != approx(qhda.gruneisen_parameter(temp, vol))


def test_thermal_conductivity(qhda_vinet):
    vol = 66.0
    gamma = qhda_vinet.gruneisen_parameter(300, vol)
    theta_a = qhda_vinet.debye_temperature(vol) / SYN_NATOMS ** (1 / 3)
    # Slack-type scaling: kappa ~ theta_a^2 V^(1/3) / (gamma^2 - 0.514 gamma + 0.228)
    ratio = qhda_vinet.thermal_conductivity(300, vol) * (gamma**2 - 0.514 * gamma + 0.228) / theta_a**2 / vol ** (1 / 3)
    assert qhda_vinet.thermal_conductivity(300, vol) == approx(12.3535, rel=1e-4)
    gamma_2 = qhda_vinet.gruneisen_parameter(300, 68.0)
    theta_a_2 = qhda_vinet.debye_temperature(68.0) / SYN_NATOMS ** (1 / 3)
    assert qhda_vinet.thermal_conductivity(300, 68.0) == approx(
        ratio * theta_a_2**2 * 68.0 ** (1 / 3) / (gamma_2**2 - 0.514 * gamma_2 + 0.228)
    )
    # with the Slater gamma kappa(V) does not depend on T
    assert qhda_vinet.thermal_conductivity(900, vol) == qhda_vinet.thermal_conductivity(300, vol)


def test_get_summary_dict():
    qhda = _qhda(t_min=300, t_max=500, t_step=200, pressure=0.5, poisson=0.3)
    summary = qhda.get_summary_dict()
    assert summary["pressure"] == 0.5
    assert summary["poisson"] == 0.3
    assert summary["mass"] == approx(4 * 26.9815385, rel=1e-6)
    assert summary["natoms"] == 4
    assert isinstance(summary["natoms"], int)
    assert summary["bulk_modulus"] == qhda.bulk_modulus
    assert summary["temperatures"] == approx([300, 500])
    assert summary["gibbs_free_energy"] == qhda.gibbs_free_energy
    assert summary["optimum_volumes"] == qhda.optimum_volumes
    vols, temps = qhda.optimum_volumes, qhda.temperatures
    assert summary["debye_temperature"] == approx([qhda.debye_temperature(v) for v in vols])
    assert summary["gruneisen_parameter"] == approx(
        [qhda.gruneisen_parameter(t, v) for t, v in zip(temps, vols, strict=True)]
    )
    assert summary["thermal_conductivity"] == approx(
        [qhda.thermal_conductivity(t, v) for t, v in zip(temps, vols, strict=True)]
    )
