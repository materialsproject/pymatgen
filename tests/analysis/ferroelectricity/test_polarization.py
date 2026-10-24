from __future__ import annotations

import logging

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.analysis.ferroelectricity.polarization import (
    EnergyTrend,
    Polarization,
    calc_ionic,
    get_nearest_site,
    get_total_ionic_dipole,
    zval_dict_from_potcar,
)
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Potcar
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.util.testing import MatSciTest
from tests.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/analysis/ferroelectricity"
bto_folders = ["nonpolar_polarization"]
bto_folders += [f"interpolation_{idx}_polarization" for idx in range(8, 0, -1)]
bto_folders += ["polar_polarization"]

structures = [Structure.from_file(f"{TEST_DIR}/{folder}/POSCAR") for folder in bto_folders]

ions = np.array(
    [
        [-44.363484, -44.363484, -44.79259988],
        [-44.324764, -44.324764, -69.43452043],
        [-44.286055, -44.286055, -69.8792077],
        [-44.247335, -44.247335, -70.32475473],
        [-44.208626, -44.208626, -70.77139856],
        [-44.169906, -44.169906, -71.21889307],
        [-44.131197, -44.131197, -71.66746921],
        [-44.092477, -44.092477, -72.1168782],
        [-44.053768, -44.053768, -72.56736141],
        [-44.015048, -44.015048, -73.01874336],
    ]
)


class TestUtils(MatSciTest):
    def setup_method(self):
        self.potcar = Potcar.from_file(f"{TEST_DIR}/POTCAR")
        self.zval_dict = {"Ba": 10, "Ti": 10, "O": 6}
        self.ions = ions
        self.structures = structures

    def test_zval_dict_from_potcar(self):
        zval_dict = zval_dict_from_potcar(self.potcar)
        assert self.zval_dict == zval_dict

    def test_get_total_ionic_dipole(self):
        p_ion = get_total_ionic_dipole(self.structures[-1], self.zval_dict)
        assert_allclose(p_ion, self.ions[-1].ravel().tolist())


class TestPolarization(MatSciTest):
    def setup_method(self):
        self.p_ions = ions
        self.p_ions_outcar = np.array(
            [
                [0, 0, 43.93437],
                [0, 0, 19.81697],
                [0, 0, 19.76076],
                [0, 0, 19.70306],
                [0, 0, 19.64372],
                [0, 0, -5.06619],
                [0, 0, -5.18997],
                [0, 0, -5.31457],
                [0, 0, -5.44026],
                [0, 0, -5.56684],
            ]
        )
        self.p_elecs = np.array(
            [
                [4.03304, -4.03304, -3.60393],
                [4.02958, 4.02958, -3.77177],
                [4.02611, 4.02611, -3.93397],
                [4.02264, 4.02263, -4.08851],
                [4.01916, 4.01914, 3.99662],
                [4.01567, 4.01565, 3.90327],
                [4.01217, 4.01214, 3.81998],
                [4.00867, 4.00863, 3.74561],
                [4.00517, 4.00512, 3.67949],
                [0.00024, 0.00019, 3.61674],
            ]
        )
        self.same_branch = np.array(
            [
                [9.76948106e-05, -9.76948108e-05, 4.59556390e-05],
                [-1.36325612e-03, -1.36325612e-03, 5.99098550e00],
                [-2.54781559e-03, -2.54781559e-03, 1.18312234e01],
                [-3.74896442e-03, -3.50709575e-03, 1.74695147e01],
                [-4.67728039e-03, -4.19508654e-03, 2.28288548e01],
                [-5.38348125e-03, -4.90281328e-03, 2.79488973e01],
                [-5.82178137e-03, -5.10304293e-03, 3.28220345e01],
                [-6.28132190e-03, -5.32598777e-03, 3.74721262e01],
                [-6.71430111e-03, -5.52382219e-03, 4.19231297e01],
                [-5.69679257e-03, -4.50996078e-03, 4.62887982e01],
            ]
        )
        self.same_branch_all_in_polar = np.array(
            [
                [9.42008214e-05, -9.42008221e-05, 4.66464162e-05],
                [-1.31996355e-03, -1.31996355e-03, 6.07559443e00],
                [-2.47709492e-03, -2.47709492e-03, 1.19773649e01],
                [-3.65986380e-03, -3.42374356e-03, 1.76543900e01],
                [-4.58474364e-03, -4.11208966e-03, 2.30300987e01],
                [-5.29836518e-03, -4.82529686e-03, 2.81459084e01],
                [-5.75282166e-03, -5.04259676e-03, 3.29954877e01],
                [-6.23177258e-03, -5.28397447e-03, 3.76040802e01],
                [-6.68784798e-03, -5.50205933e-03, 4.19969216e01],
                [-5.69679257e-03, -4.50996079e-03, 4.62887981e01],
            ]
        )
        self.quanta = np.array(
            [
                [98.50186747, 98.50186747, 98.50186747],
                [98.09416498, 98.09416498, 98.67403571],
                [97.69065056, 97.69065056, 98.84660662],
                [97.29131054, 97.29131054, 99.01967988],
                [96.89603543, 96.89603543, 99.19315873],
                [96.50481368, 96.50481368, 99.36714337],
                [96.11753848, 96.11753848, 99.54153654],
                [95.7342003, 95.7342003, 99.71643897],
                [95.35469487, 95.35469487, 99.89175289],
                [94.97901455, 94.97901455, 100.06757957],
            ]
        )
        self.structures = structures
        # We do not use the p_ions values from Outcar.
        # We calculate using calc_ionic_from_zval because it is more reliable.
        self.polarization = Polarization(self.p_elecs, self.p_ions, self.structures)
        self.outcars = [Outcar(f"{TEST_DIR}/{folder}/OUTCAR") for folder in bto_folders]
        self.change = np.array([[-5.79448738e-03, -4.41226597e-03, 4.62887522e01]])
        self.change_norm = 46.288752795325244
        self.max_jumps = [
            0.00021336004941047062,
            0.00016254800426403291,
            0.038269946959965086,
        ]
        self.smoothness = [
            0.00017013512377086267,
            0.00013467465540412905,
            0.034856268571937743,
        ]
        self.max_jumps_all_in_polar = [
            0.0002131725432761777,
            0.00016247151626362123,
            0.03857992184016461,
        ]
        self.smoothness_all_in_polar = [
            0.00016974252210569685,
            0.0001343093739679674,
            0.03504919463271141,
        ]
        self.decimal_tol = 5

    def test_from_outcars_and_structures(self):
        polarization = Polarization.from_outcars_and_structures(self.outcars, self.structures)
        p_elecs, p_ions = polarization.get_pelecs_and_pions(convert_to_muC_per_cm2=False)
        assert_allclose(p_elecs[0].ravel().tolist(), self.p_elecs[0].ravel().tolist())
        assert_allclose(p_elecs[-1].ravel().tolist(), self.p_elecs[-1].ravel().tolist())
        assert_allclose(p_ions[0].ravel().tolist(), self.p_ions_outcar[0].ravel().tolist())
        assert_allclose(p_ions[-1].ravel().tolist(), self.p_ions_outcar[-1].ravel().tolist())
        # Test for calc_ionic_from_zval=True
        polarization = Polarization.from_outcars_and_structures(
            self.outcars, self.structures, calc_ionic_from_zval=True
        )
        p_elecs, p_ions = polarization.get_pelecs_and_pions(convert_to_muC_per_cm2=False)
        assert_allclose(p_elecs[0].ravel().tolist(), self.p_elecs[0].ravel().tolist())
        assert_allclose(p_elecs[-1].ravel().tolist(), self.p_elecs[-1].ravel().tolist())
        assert_allclose(p_ions[0].ravel().tolist(), self.p_ions[0].ravel().tolist())
        assert_allclose(p_ions[-1].ravel().tolist(), self.p_ions[-1].ravel().tolist())

    def test_get_same_branch_polarization_data(self):
        same_branch = self.polarization.get_same_branch_polarization_data(
            convert_to_muC_per_cm2=True, all_in_polar=False
        )
        assert_allclose(
            same_branch[0].ravel().tolist(),
            self.same_branch[0].ravel().tolist(),
            atol=1e-7,
        )
        assert_allclose(
            same_branch[1].ravel().tolist(),
            self.same_branch[1].ravel().tolist(),
            atol=1e-7,
        )
        assert_allclose(
            same_branch[3].ravel().tolist(),
            self.same_branch[3].ravel().tolist(),
            atol=1e-7,
        )
        assert_allclose(
            same_branch[-1].ravel().tolist(),
            self.same_branch[-1].ravel().tolist(),
            atol=1e-7,
        )
        # This will differ only slightly
        same_branch = self.polarization.get_same_branch_polarization_data(
            convert_to_muC_per_cm2=True, all_in_polar=True
        )
        assert_allclose(
            same_branch[0].ravel().tolist(),
            self.same_branch_all_in_polar[0].ravel().tolist(),
        )
        assert_allclose(
            same_branch[1].ravel().tolist(),
            self.same_branch_all_in_polar[1].ravel().tolist(),
        )
        assert_allclose(
            same_branch[3].ravel().tolist(),
            self.same_branch_all_in_polar[3].ravel().tolist(),
        )
        assert_allclose(
            same_branch[-1].ravel().tolist(),
            self.same_branch_all_in_polar[-1].ravel().tolist(),
        )

    def test_get_lattice_quanta(self):
        quanta = self.polarization.get_lattice_quanta(convert_to_muC_per_cm2=True, all_in_polar=False)
        assert_allclose(quanta[0].ravel().tolist(), self.quanta[0].ravel().tolist())
        assert_allclose(quanta[-1].ravel().tolist(), self.quanta[-1].ravel().tolist())
        # For all_in_polar=True, quanta should be identical to polar quantum
        quanta = self.polarization.get_lattice_quanta(convert_to_muC_per_cm2=True, all_in_polar=True)
        assert_allclose(quanta[0].ravel().tolist(), self.quanta[-1].ravel().tolist())
        assert_allclose(quanta[-1].ravel().tolist(), self.quanta[-1].ravel().tolist())

    def test_get_polarization_change(self):
        change = self.polarization.get_polarization_change(convert_to_muC_per_cm2=True, all_in_polar=False)
        assert_allclose(change, self.change)
        # Because nonpolar polarization is (0, 0, 0), all_in_polar should have no effect on polarization change
        change = self.polarization.get_polarization_change(convert_to_muC_per_cm2=True, all_in_polar=True)
        # No change up to 5 decimal
        assert_allclose(change, self.change, self.decimal_tol)

    def test_get_polarization_change_norm(self):
        change_norm = self.polarization.get_polarization_change_norm(convert_to_muC_per_cm2=True, all_in_polar=False)
        assert change_norm == approx(self.change_norm)
        # Because nonpolar polarization is (0, 0, 0), all_in_polar should have no effect on polarization change norm
        change = self.polarization.get_polarization_change(convert_to_muC_per_cm2=True, all_in_polar=True)
        # No change up to 5 decimal
        assert_allclose(change, self.change, self.decimal_tol)

    def test_max_spline_jumps(self):
        max_jumps = self.polarization.max_spline_jumps(convert_to_muC_per_cm2=True, all_in_polar=False)
        assert_allclose(self.max_jumps, max_jumps, atol=1e-7)
        # This will differ slightly
        max_jumps = self.polarization.max_spline_jumps(convert_to_muC_per_cm2=True, all_in_polar=True)
        assert_allclose(self.max_jumps_all_in_polar, max_jumps)

    def test_smoothness(self):
        smoothness = self.polarization.smoothness(convert_to_muC_per_cm2=True, all_in_polar=False)
        assert_allclose(self.smoothness, smoothness, atol=1e-7)
        # This will differ slightly
        smoothness = self.polarization.smoothness(convert_to_muC_per_cm2=True, all_in_polar=True)
        assert_allclose(self.smoothness_all_in_polar, smoothness)


class TestEnergyTrend(MatSciTest):
    def setup_method(self):
        self.energies = [
            -7.97738049,
            -7.988621176,
            -7.9793246479999995,
            -7.987973192,
            -7.984676138,
            -7.982700144000001,
            -7.986539788,
            -7.980859048000001,
            -7.978240114,
            -7.977637218,
        ]
        self.energy_trend = EnergyTrend(self.energies)
        self.smoothness = 0.0029874731764648306
        self.max_jump = 0.0058893082867133018

    def test_max_spline_jump(self):
        max_jump = self.energy_trend.max_spline_jump()
        assert max_jump == approx(self.max_jump)

    def test_smoothness(self):
        smoothness = self.energy_trend.smoothness()
        assert smoothness == approx(self.smoothness)

    def test_endpoints_minima(self):
        endpoints = self.energy_trend.endpoints_minima(slope_cutoff=1e-2)
        assert endpoints == {"polar": True, "nonpolar": True}


E_ANG2_TO_MUC_CM2 = 1.6021766e-13 * 1e16  # e / Angstrom^2 -> muC / cm^2
LOGGER_NAME = "pymatgen.analysis.ferroelectricity.polarization"
# Electronic dipoles (e A, Cartesian as given by VASP) along the BaTiO3 path
P_ELECS = np.array(
    [
        [4.03304, -4.03304, -3.60393],
        [4.02958, 4.02958, -3.77177],
        [4.02611, 4.02611, -3.93397],
        [4.02264, 4.02263, -4.08851],
        [4.01916, 4.01914, 3.99662],
        [4.01567, 4.01565, 3.90327],
        [4.01217, 4.01214, 3.81998],
        [4.00867, 4.00863, 3.74561],
        [4.00517, 4.00512, 3.67949],
        [0.00024, 0.00019, 3.61674],
    ]
)


@pytest.mark.parametrize("radius", [None, 5.0])
def test_get_nearest_site(radius):
    struct = Structure(Lattice.cubic(10), ["Na", "Cl"], [[0.9, 0, 0], [0.5, 0.5, 0.5]])
    site, dist = get_nearest_site(struct, [0, 0, 0], struct[0], r=radius)
    # The periodic image at x = -1 A is closer to the origin than the one at x = 9 A
    assert dist == approx(1.0)
    assert_allclose(site.coords, [-1, 0, 0], atol=1e-10)
    assert site.specie.symbol == "Na"


def test_calc_ionic():
    struct = structures[-1]
    ti_site = next(site for site in struct if site.specie.symbol == "Ti")
    p_ion = calc_ionic(ti_site, struct, 10)
    assert_allclose(p_ion, -10 * np.array(struct.lattice.abc) * ti_site.frac_coords)
    assert_allclose(calc_ionic(ti_site, struct, 0), np.zeros(3))


class TestPolarizationExtra:
    """Unit conversions, branch/quantum invariants and failure branches."""

    @pytest.fixture
    def polarization(self):
        return Polarization(P_ELECS, ions, structures)

    def test_mismatched_lengths(self):
        with pytest.raises(ValueError, match="must be equal"):
            Polarization(P_ELECS[:-1], ions, structures)
        with pytest.raises(ValueError, match="must be equal"):
            Polarization(P_ELECS, ions, structures[:-1])

    def test_cartesian_flags(self):
        # The BaTiO3 cells are axis-aligned, so Cartesian and lattice-direction components coincide
        pol = Polarization(P_ELECS, ions, structures, p_elecs_in_cartesian=False, p_ions_in_cartesian=True)
        assert_allclose(pol.p_elecs, P_ELECS)
        assert_allclose(pol.p_ions, ions, atol=1e-10)

        # Non-orthogonal cell: lattice-direction components are projections onto unit a, b, c
        struct = Structure(Lattice.hexagonal(3, 5), ["C"], [[0, 0, 0]])
        cart = np.array([[0.0, 1.0, 2.0]])
        pol = Polarization(cart, cart, [struct], p_elecs_in_cartesian=True, p_ions_in_cartesian=True)
        expected = struct.lattice.get_vector_along_lattice_directions(cart[0])
        assert not np.allclose(expected, cart[0])
        assert_allclose(pol.p_elecs[0], expected)
        assert_allclose(pol.p_ions[0], expected)

    def test_get_pelecs_and_pions_units(self, polarization):
        p_elecs, p_ions = polarization.get_pelecs_and_pions()
        assert_allclose(p_elecs, P_ELECS, atol=1e-10)
        assert_allclose(p_ions, ions)

        p_elecs_muC, p_ions_muC = polarization.get_pelecs_and_pions(convert_to_muC_per_cm2=True)
        volumes = np.array([struct.volume for struct in structures])[:, None]
        # Electron charge is negative: dipole (e A) -> polarization (muC/cm^2) flips sign
        assert_allclose(p_elecs_muC, -E_ANG2_TO_MUC_CM2 * p_elecs / volumes)
        assert_allclose(p_ions_muC, -E_ANG2_TO_MUC_CM2 * p_ions / volumes)
        # Nonpolar BaTiO3 ionic term along c: 44.79 e A / 65.60 A^3
        assert p_ions_muC[0, 2] == approx(1094.00, abs=0.01)

    def test_lattice_quanta(self, polarization):
        abc = np.array([struct.lattice.abc for struct in structures])
        volumes = np.array([struct.volume for struct in structures])[:, None]
        # Dipole quanta are the lattice lengths (e = 1)
        assert_allclose(polarization.get_lattice_quanta(convert_to_muC_per_cm2=False), abc)
        # Polarization quanta are e R_i / Omega
        assert_allclose(
            polarization.get_lattice_quanta(convert_to_muC_per_cm2=True, all_in_polar=False),
            E_ANG2_TO_MUC_CM2 * abc / volumes,
        )
        # Cubic nonpolar cell: e / a^2 with a = 4.033044 A
        assert_allclose(polarization.get_lattice_quanta(all_in_polar=False)[0], 98.5019, atol=1e-4)

    @pytest.mark.parametrize("convert", [False, True])
    def test_same_branch_differs_by_whole_quanta(self, polarization, convert):
        p_elecs, p_ions = polarization.get_pelecs_and_pions(convert_to_muC_per_cm2=convert)
        raw = p_elecs + p_ions
        same_branch = polarization.get_same_branch_polarization_data(convert_to_muC_per_cm2=convert, all_in_polar=False)
        quanta = polarization.get_lattice_quanta(convert_to_muC_per_cm2=convert, all_in_polar=False)
        n_quanta = (same_branch - raw) / quanta
        assert_allclose(n_quanta, np.round(n_quanta), atol=1e-6)
        # The nonpolar image is the one closest to zero ...
        assert np.all(np.abs(same_branch[0]) <= quanta[0] / 2 + 1e-8)
        # ... and consecutive images stay within half a quantum of each other
        assert np.all(np.abs(np.diff(same_branch, axis=0)) <= quanta[1:] / 2)

    def test_polarization_change_in_dipole_units(self, polarization):
        same_branch = polarization.get_same_branch_polarization_data(convert_to_muC_per_cm2=False)
        change = polarization.get_polarization_change(convert_to_muC_per_cm2=False)
        assert change.shape == (1, 3)
        assert_allclose(change[0], same_branch[-1] - same_branch[0])
        # Axis-aligned polar cell: the norm is the Euclidean norm of the components
        norm = polarization.get_polarization_change_norm(convert_to_muC_per_cm2=False)
        assert norm == approx(np.linalg.norm(change))

    def test_too_few_points_for_splines(self):
        # A cubic spline needs more than 3 points: each component degrades to None
        pol = Polarization(P_ELECS[:3], ions[:3], structures[:3])
        assert pol.same_branch_splines() == (None, None, None)
        assert pol.max_spline_jumps() == [None, None, None]
        assert pol.smoothness() == [None, None, None]

    def test_smoothness_spline_error(self, polarization, monkeypatch, caplog):
        def raise_error(**_kwargs):
            raise RuntimeError("spline failed")

        monkeypatch.setattr(polarization, "same_branch_splines", raise_error)
        with caplog.at_level(logging.ERROR, logger=LOGGER_NAME):
            assert polarization.smoothness() is None
        assert "Something went wrong" in caplog.text


class TestEnergyTrendExtra:
    def test_quadratic_trend(self):
        energies = 0.01 * (np.arange(10) - 4.5) ** 2
        trend = EnergyTrend(energies)
        assert_allclose(trend.spline()(range(10)), energies, atol=1e-10)
        assert trend.smoothness() == approx(0, abs=1e-10)
        assert trend.max_spline_jump() == approx(0, abs=1e-10)
        # Endpoint slopes are -/+0.09 per image
        assert trend.endpoints_minima() == {"polar": False, "nonpolar": False}
        assert trend.endpoints_minima(slope_cutoff=0.1) == {"polar": True, "nonpolar": True}

    def test_too_few_points(self, caplog):
        # A quartic spline needs more than 4 points
        trend = EnergyTrend([1.0, 0.5, 1.0])
        with caplog.at_level(logging.ERROR, logger=LOGGER_NAME):
            assert trend.smoothness() is None
            assert trend.endpoints_minima() is None
        assert caplog.text.count("Energy spline failed") == 2
        # Older scipy raises dfitpack.error "(m>k) failed", newer ValueError "m must be > k"
        with pytest.raises(Exception, match=r"m must be > k|\(m>k\) failed"):
            trend.max_spline_jump()
