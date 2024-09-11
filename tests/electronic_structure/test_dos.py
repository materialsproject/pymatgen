from __future__ import annotations

import json
import re
from unittest import TestCase

import numpy as np
import pytest
from monty.io import zopen
from monty.serialization import loadfn
from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.core import Element, Structure
from pymatgen.electronic_structure.core import Orbital, OrbitalType, Spin
from pymatgen.electronic_structure.dos import DOS, CompleteDos, FermiDos, LobsterCompleteDos
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

TEST_DIR = f"{TEST_FILES_DIR}/electronic_structure/dos"


class TestDos(TestCase):
    def setUp(self):
        with open(f"{TEST_DIR}/complete_dos.json") as file:
            self.dos = CompleteDos.from_dict(json.load(file))

    def test_get_gap(self):
        assert self.dos.get_gap() == approx(2.0589, abs=1e-4)
        assert len(self.dos.energies) == 301
        assert self.dos.get_interpolated_gap(tol=0.001, abs_tol=False, spin=None)[0] == approx(
            2.16815942458015, abs=1e-7
        )
        assert self.dos.get_cbm_vbm() == approx((3.8729, 1.8140000000000001))

        assert self.dos.get_interpolated_value(9.9)[Spin.up] == approx(1.744588888888891, abs=1e-7)
        assert self.dos.get_interpolated_value(9.9)[Spin.down] == approx(1.756888888888886, abs=1e-7)
        with pytest.raises(ValueError, match=r"x=1000 is out of range of provided x_values \(-23.7934, 14.8107\)"):
            self.dos.get_interpolated_value(1000)

    def test_get_smeared_densities(self):
        smeared = self.dos.get_smeared_densities(0.2)
        dens = self.dos.densities
        for spin in Spin:
            assert sum(dens[spin]) == approx(sum(smeared[spin]))

    def test_as_dict(self):
        dos_dict = self.dos.as_dict()
        assert isinstance(dos_dict["energies"], list)
        assert isinstance(dos_dict["energies"][0], float)
        assert not isinstance(dos_dict["energies"][0], np.float64)

        assert isinstance(dos_dict["densities"]["1"], list)
        assert isinstance(dos_dict["densities"]["1"][0], float)
        assert not isinstance(dos_dict["densities"]["1"][0], np.float64)


class TestFermiDos(TestCase):
    def setUp(self):
        with open(f"{TEST_DIR}/complete_dos.json") as file:
            self.dos = CompleteDos.from_dict(json.load(file))
        self.dos = FermiDos(self.dos)

    def test_doping_fermi(self):
        T = 300
        fermi0 = self.dos.efermi
        fermi_range = [fermi0 - 0.5, fermi0, fermi0 + 2.0, fermi0 + 2.2]
        dopings = [self.dos.get_doping(fermi_level=fermi_lvl, temperature=T) for fermi_lvl in fermi_range]
        ref_dopings = [3.48077e21, 1.9235e18, -2.6909e16, -4.8723e19]
        for idx, c_ref in enumerate(ref_dopings):
            assert abs(dopings[idx] / c_ref - 1.0) <= 0.01

        calc_fermis = [self.dos.get_fermi(concentration=c, temperature=T) for c in ref_dopings]
        for j, f_ref in enumerate(fermi_range):
            assert calc_fermis[j] == approx(f_ref, abs=1e-4)

        sci_dos = FermiDos(self.dos, bandgap=3.0)
        assert sci_dos.get_gap() == 3.0
        old_cbm, old_vbm = self.dos.get_cbm_vbm()
        old_gap = old_cbm - old_vbm
        new_cbm, new_vbm = sci_dos.get_cbm_vbm()
        assert new_cbm - old_cbm == approx((3.0 - old_gap) / 2.0)
        assert old_vbm - new_vbm == approx((3.0 - old_gap) / 2.0)
        for idx, c_ref in enumerate(ref_dopings):
            if c_ref < 0:
                assert sci_dos.get_fermi(c_ref, temperature=T) - fermi_range[idx] == approx(0.47, abs=1e-2)
            else:
                assert sci_dos.get_fermi(c_ref, temperature=T) - fermi_range[idx] == approx(-0.47, abs=1e-2)

        assert sci_dos.get_fermi_interextrapolated(-1e26, 300) == approx(7.5108, abs=1e-4)
        assert sci_dos.get_fermi_interextrapolated(1e26, 300) == approx(-1.4182, abs=1e-4)
        assert sci_dos.get_fermi_interextrapolated(0.0, 300) == approx(2.9071, abs=1e-4)

    def test_as_dict(self):
        dos_dict = self.dos.as_dict()
        assert isinstance(dos_dict["energies"], list)
        assert isinstance(dos_dict["energies"][0], float)
        assert not isinstance(dos_dict["energies"][0], np.float64)

        assert isinstance(dos_dict["densities"]["1"], list)
        assert isinstance(dos_dict["densities"]["1"][0], float)
        assert not isinstance(dos_dict["densities"]["1"][0], np.float64)


class TestCompleteDos(TestCase):
    def setUp(self):
        with open(f"{TEST_DIR}/complete_dos.json") as file:
            self.dos = CompleteDos.from_dict(json.load(file))
        with zopen(f"{TEST_DIR}/pdag3_complete_dos.json.gz") as file:
            self.dos_pdag3 = CompleteDos.from_dict(json.load(file))

    def test_get_gap(self):
        assert self.dos.get_gap() == approx(2.0589, abs=1e-4), "Wrong gap from dos!"
        assert len(self.dos.energies) == 301
        assert self.dos.get_interpolated_gap(tol=0.001, abs_tol=False, spin=None)[0] == approx(
            2.16815942458015, abs=1e-7
        )
        spd_dos = self.dos.get_spd_dos()
        assert len(spd_dos) == 3
        el_dos = self.dos.get_element_dos()
        assert len(el_dos) == 4
        sum_spd = spd_dos[OrbitalType.s] + spd_dos[OrbitalType.p] + spd_dos[OrbitalType.d]
        sum_element = None
        for pdos in el_dos.values():
            if sum_element is None:
                sum_element = pdos
            else:
                sum_element += pdos

        # The sums of the SPD or the element doses should be the same.
        assert (abs(sum_spd.energies - sum_element.energies) < 0.0001).all()
        assert (abs(sum_spd.densities[Spin.up] - sum_element.densities[Spin.up]) < 0.0001).all()
        assert (abs(sum_spd.densities[Spin.down] - sum_element.densities[Spin.down]) < 0.0001).all()

        site = self.dos.structure[0]
        assert self.dos.get_site_dos(site) is not None
        assert sum(self.dos.get_site_dos(site).get_densities(Spin.up)) == approx(2.0391)
        assert sum(self.dos.get_site_dos(site).get_densities(Spin.down)) == approx(2.0331999999999995)
        assert self.dos.get_site_orbital_dos(site, Orbital.s) is not None
        egt2g = self.dos.get_site_t2g_eg_resolved_dos(site)
        assert sum(egt2g["e_g"].get_densities(Spin.up)) == approx(0.0)
        assert sum(egt2g["t2g"].get_densities(Spin.up)) == approx(0.0)
        egt2g = self.dos.get_site_t2g_eg_resolved_dos(self.dos.structure[4])
        assert sum(egt2g["e_g"].get_densities(Spin.up)) == approx(15.004399999999997)
        assert sum(egt2g["t2g"].get_densities(Spin.up)) == approx(22.910399999999999)
        assert self.dos.get_cbm_vbm() == approx((3.8729, 1.8140000000000001))

        assert self.dos.get_interpolated_value(9.9)[Spin.up] == approx(1.744588888888891, abs=1e-7)
        assert self.dos.get_interpolated_value(9.9)[Spin.down] == approx(1.756888888888886, abs=1e-7)
        with pytest.raises(ValueError, match=r"x=1000 is out of range of provided x_values \(-23.7934, 14.8107\)"):
            self.dos.get_interpolated_value(1000)

    def test_as_from_dict(self):
        dct = self.dos.as_dict()
        dos = CompleteDos.from_dict(dct)
        el_dos = dos.get_element_dos()
        assert len(el_dos) == 4
        spd_dos = dos.get_spd_dos()
        sum_spd = spd_dos[OrbitalType.s] + spd_dos[OrbitalType.p] + spd_dos[OrbitalType.d]
        sum_element = None
        for pdos in el_dos.values():
            if sum_element is None:
                sum_element = pdos
            else:
                sum_element += pdos

        # The sums of the SPD or the element doses should be the same.
        assert (abs(sum_spd.energies - sum_element.energies) < 0.0001).all()

    def test_str(self):
        assert str(self.dos).startswith("Complete DOS for Full Formula (Li1 Fe4 P4 O16)\nReduced Formula: LiFe4(PO4)4")

    def test_as_dict(self):
        dos_dict = self.dos.as_dict()
        assert isinstance(dos_dict["energies"], list)
        assert isinstance(dos_dict["energies"][0], float)
        assert not isinstance(dos_dict["energies"][0], np.float64)

        assert isinstance(dos_dict["densities"]["1"], list)
        assert isinstance(dos_dict["densities"]["1"][0], float)
        assert not isinstance(dos_dict["densities"]["1"][0], np.float64)

    def test_get_band_center(self):
        dos = self.dos_pdag3
        struct = dos.structure
        band_center = dos.get_band_center()
        assert band_center == approx(-3.078841005723767)

        band_center = dos.get_band_center(elements=[Element("Ag"), Element("Pd")])
        assert band_center == approx(-3.078841005723767)

        band_center = dos.get_band_center(elements=[Element("Pd")])
        assert band_center == approx(-1.476449501704171)

        band_center = dos.get_band_center(sites=[s for s in struct if s.species_string == "Pd"])
        assert band_center == approx(-1.476449501704171)

        band_center = dos.get_band_center(sites=[struct[-3]])
        assert band_center == approx(-1.4144921311083436)

        band_center = dos.get_band_center(band=OrbitalType.p)
        assert band_center == approx(0.9430322204760462)

        band_center = dos.get_band_center(elements=[Element("Pd")], band=OrbitalType.p)
        assert band_center == approx(0.7825770239165218)

        band_center = dos.get_band_center(elements=[Element("Pd")], erange=[-4, -2])
        assert band_center == approx(-2.8754000116714065)

    def test_get_upper_band_edge(self):
        dos = self.dos_pdag3
        struct = dos.structure
        band_edge = dos.get_upper_band_edge()
        assert band_edge == approx(-1.01246969)

        band_edge = dos.get_upper_band_edge(elements=[Element("Pd")])
        assert band_edge == approx(-1.01246969)

        band_edge = dos.get_upper_band_edge(sites=[struct[-3]])
        assert band_edge == approx(-1.01246969)

        band_edge = dos.get_upper_band_edge(elements=[Element("Pd")], erange=[-4, 0.5])
        assert band_edge == approx(-1.01246969)

    def test_get_n_moment(self):
        dos = self.dos_pdag3
        moment = dos.get_n_moment(1)
        assert moment == approx(0)

        moment = dos.get_n_moment(1, center=False)
        assert moment == approx(-3.078841005723767)

    def test_band_filling(self):
        dos = self.dos_pdag3
        filling = dos.get_band_filling()
        assert filling == approx(0.9583552024357637)

    def test_band_width(self):
        dos = self.dos_pdag3
        width = dos.get_band_width()
        assert width == approx(1.7831724662185575)

    def test_get_band_skewness(self):
        dos = self.dos_pdag3
        skewness = dos.get_band_skewness()
        assert skewness == approx(1.7422716340493507)

    def test_get_band_kurtosis(self):
        dos = self.dos_pdag3
        kurtosis = dos.get_band_kurtosis()
        assert kurtosis == approx(7.764506941340621)

    def test_get_dos_fp(self):
        # normalize is True
        dos_fp = self.dos.get_dos_fp(fp_type="s", min_e=-10, max_e=0, n_bins=56, normalize=True)
        bin_width = np.diff(dos_fp.energies)[0][0]
        assert max(dos_fp.energies[0]) <= 0
        assert min(dos_fp.energies[0]) >= -10
        assert len(dos_fp.energies[0]) == 56
        assert dos_fp.fp_type == "s"
        assert sum(dos_fp.densities * bin_width) == approx(1)
        # normalize is False
        dos_fp2 = self.dos.get_dos_fp(fp_type="s", min_e=-10, max_e=0, n_bins=56, normalize=False)
        bin_width2 = np.diff(dos_fp2.energies)[0][0]
        assert sum(dos_fp2.densities * bin_width2) == approx(7.279303571428509)
        assert dos_fp2.bin_width == approx(bin_width2)
        # binning is False
        dos_fp = self.dos.get_dos_fp(fp_type="s", min_e=None, max_e=None, n_bins=56, normalize=True, binning=False)
        assert dos_fp.n_bins == len(self.dos.energies)

    def test_get_dos_fp_similarity(self):
        # Tanimoto
        dos_fp = self.dos.get_dos_fp(fp_type="s", min_e=-10, max_e=0, n_bins=56, normalize=True)
        dos_fp2 = self.dos.get_dos_fp(fp_type="tdos", min_e=-10, max_e=0, n_bins=56, normalize=True)
        similarity_index = self.dos.get_dos_fp_similarity(dos_fp, dos_fp2, col=1, metric="tanimoto")
        assert similarity_index == approx(0.3342481451042263)

        dos_fp = self.dos.get_dos_fp(fp_type="s", min_e=-10, max_e=0, n_bins=56, normalize=True)
        dos_fp2 = self.dos.get_dos_fp(fp_type="s", min_e=-10, max_e=0, n_bins=56, normalize=True)
        similarity_index = self.dos.get_dos_fp_similarity(dos_fp, dos_fp2, col=1, metric="tanimoto")
        assert similarity_index == approx(1)

        # Wasserstein
        dos_fp = self.dos.get_dos_fp(fp_type="s", min_e=-10, max_e=0, n_bins=56, normalize=True)
        dos_fp2 = self.dos.get_dos_fp(fp_type="tdos", min_e=-10, max_e=0, n_bins=56, normalize=True)
        similarity_index = self.dos.get_dos_fp_similarity(dos_fp, dos_fp2, col=1, metric="wasserstein")
        assert similarity_index == approx(0.2668440595873588)

    def test_dos_fp_exceptions(self):
        dos_fp = self.dos.get_dos_fp(fp_type="s", min_e=-10, max_e=0, n_bins=56, normalize=True)
        dos_fp2 = self.dos.get_dos_fp(fp_type="tdos", min_e=-10, max_e=0, n_bins=56, normalize=True)
        # test exceptions
        with pytest.raises(
            ValueError,
            match="Cannot compute similarity index. When normalize=True, then please set metric=cosine-sim",
        ):
            self.dos.get_dos_fp_similarity(dos_fp, dos_fp2, col=1, metric="tanimoto", normalize=True)
        with pytest.raises(
            ValueError,
            match="Please recheck fp_type requested, either the orbital "
            "projections unavailable in input DOS or there's a typo in type.",
        ):
            self.dos.get_dos_fp(fp_type="k", min_e=-10, max_e=0, n_bins=56, normalize=True)

        valid_metrics = ("tanimoto", "wasserstein", "cosine-sim")
        metric = "Dot"
        with pytest.raises(ValueError, match=re.escape(f"Invalid {metric=}, choose from {valid_metrics}.")):
            self.dos.get_dos_fp_similarity(dos_fp, dos_fp2, col=1, metric=metric, normalize=False)


class TestDOS(PymatgenTest):
    def setUp(self):
        with open(f"{TEST_DIR}/complete_dos.json") as file:
            dct = json.load(file)
            ys = list(zip(dct["densities"]["1"], dct["densities"]["-1"], strict=True))
            self.dos = DOS(dct["energies"], ys, dct["efermi"])

    def test_get_gap(self):
        assert self.dos.get_gap() == approx(2.0589, abs=1e-4)
        assert len(self.dos.x) == 301
        assert self.dos.get_interpolated_gap(tol=0.001, abs_tol=False, spin=None)[0] == approx(
            2.16815942458015, abs=1e-7
        )
        assert_allclose(self.dos.get_cbm_vbm(), (3.8729, 1.8140000000000001))

        assert self.dos.get_interpolated_value(9.9)[0] == approx(1.744588888888891, abs=1e-7)
        assert self.dos.get_interpolated_value(9.9)[1] == approx(1.756888888888886, abs=1e-7)
        with pytest.raises(ValueError, match=r"x=1000 is out of range of provided x_values \(-23.7934, 14.8107\)"):
            self.dos.get_interpolated_value(1000)

        assert_allclose(self.dos.get_cbm_vbm(spin=Spin.up), (3.8729, 1.2992999999999999))

        assert_allclose(self.dos.get_cbm_vbm(spin=Spin.down), (4.645, 1.8140000000000001))


class TestSpinPolarization(TestCase):
    def test_spin_polarization(self):
        dos_path = f"{TEST_DIR}/dos_spin_polarization_mp-865805.json"
        dos = loadfn(dos_path)
        assert dos.spin_polarization == approx(0.6460514663341762)


class TestLobsterCompleteDos(TestCase):
    def setUp(self):
        with open(f"{TEST_DIR}/LobsterCompleteDos_spin.json") as file:
            data_spin = json.load(file)
        self.LobsterCompleteDOS_spin = LobsterCompleteDos.from_dict(data_spin)

        with open(f"{TEST_DIR}/LobsterCompleteDos_nonspin.json") as file:
            data_nonspin = json.load(file)
        self.LobsterCompleteDOS_nonspin = LobsterCompleteDos.from_dict(data_nonspin)

        with open(f"{TEST_DIR}/structure_KF.json") as file:
            data_structure = json.load(file)
        self.structure = Structure.from_dict(data_structure)

        with open(f"{TEST_DIR}/LobsterCompleteDos_MnO.json") as file:
            data_MnO = json.load(file)
        self.LobsterCompleteDOS_MnO = LobsterCompleteDos.from_dict(data_MnO)

        with open(f"{TEST_DIR}/LobsterCompleteDos_MnO_nonspin.json") as file:
            data_MnO_nonspin = json.load(file)
        self.LobsterCompleteDOS_MnO_nonspin = LobsterCompleteDos.from_dict(data_MnO_nonspin)

        with open(f"{TEST_DIR}/structure_MnO.json") as file:
            data_MnO = json.load(file)
        self.structure_MnO = Structure.from_dict(data_MnO)

    def test_get_site_orbital_dos(self):
        # with spin polarization
        energies_spin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        fermi = 0.0

        pdos_f_2s_up = [0.00000, 0.00159, 0.00000, 0.00011, 0.00000, 0.00069]
        pdos_f_2s_down = [0.00000, 0.00159, 0.00000, 0.00011, 0.00000, 0.00069]
        pdos_f_2py_up = [0.00000, 0.00160, 0.00000, 0.25801, 0.00000, 0.00029]
        pdos_f_2py_down = [0.00000, 0.00161, 0.00000, 0.25819, 0.00000, 0.00029]
        pdos_f_2pz_up = [0.00000, 0.00161, 0.00000, 0.25823, 0.00000, 0.00029]
        pdos_f_2pz_down = [0.00000, 0.00160, 0.00000, 0.25795, 0.00000, 0.00029]
        pdos_f_2px_up = [0.00000, 0.00160, 0.00000, 0.25805, 0.00000, 0.00029]
        pdos_f_2px_down = [0.00000, 0.00161, 0.00000, 0.25814, 0.00000, 0.00029]
        assert (
            self.LobsterCompleteDOS_spin.get_site_orbital_dos(site=self.structure[0], orbital="2s").energies.tolist()
            == energies_spin
        )
        assert self.LobsterCompleteDOS_spin.get_site_orbital_dos(site=self.structure[0], orbital="2s").efermi == approx(
            fermi
        )
        assert (
            self.LobsterCompleteDOS_spin.get_site_orbital_dos(site=self.structure[0], orbital="2s")
            .densities[Spin.up]
            .tolist()
            == pdos_f_2s_up
        )
        assert (
            self.LobsterCompleteDOS_spin.get_site_orbital_dos(site=self.structure[0], orbital="2s")
            .densities[Spin.down]
            .tolist()
            == pdos_f_2s_down
        )
        assert (
            self.LobsterCompleteDOS_spin.get_site_orbital_dos(site=self.structure[0], orbital="2p_z").energies.tolist()
            == energies_spin
        )
        assert self.LobsterCompleteDOS_spin.get_site_orbital_dos(
            site=self.structure[0], orbital="2p_z"
        ).efermi == approx(fermi)

        assert (
            self.LobsterCompleteDOS_spin.get_site_orbital_dos(site=self.structure[0], orbital="2p_y")
            .densities[Spin.up]
            .tolist()
            == pdos_f_2py_up
        )
        assert (
            self.LobsterCompleteDOS_spin.get_site_orbital_dos(site=self.structure[0], orbital="2p_y")
            .densities[Spin.down]
            .tolist()
            == pdos_f_2py_down
        )
        assert (
            self.LobsterCompleteDOS_spin.get_site_orbital_dos(site=self.structure[0], orbital="2p_y").energies.tolist()
            == energies_spin
        )
        assert self.LobsterCompleteDOS_spin.get_site_orbital_dos(
            site=self.structure[0], orbital="2p_y"
        ).efermi == approx(fermi)

        assert (
            self.LobsterCompleteDOS_spin.get_site_orbital_dos(site=self.structure[0], orbital="2p_z")
            .densities[Spin.up]
            .tolist()
            == pdos_f_2pz_up
        )
        assert (
            self.LobsterCompleteDOS_spin.get_site_orbital_dos(site=self.structure[0], orbital="2p_z")
            .densities[Spin.down]
            .tolist()
            == pdos_f_2pz_down
        )
        assert self.LobsterCompleteDOS_spin.get_site_orbital_dos(
            site=self.structure[0], orbital="2p_z"
        ).efermi == approx(fermi)

        assert (
            self.LobsterCompleteDOS_spin.get_site_orbital_dos(site=self.structure[0], orbital="2p_x").energies.tolist()
            == energies_spin
        )
        assert (
            self.LobsterCompleteDOS_spin.get_site_orbital_dos(site=self.structure[0], orbital="2p_x")
            .densities[Spin.up]
            .tolist()
            == pdos_f_2px_up
        )
        assert (
            self.LobsterCompleteDOS_spin.get_site_orbital_dos(site=self.structure[0], orbital="2p_x")
            .densities[Spin.down]
            .tolist()
            == pdos_f_2px_down
        )
        assert self.LobsterCompleteDOS_spin.get_site_orbital_dos(
            site=self.structure[0], orbital="2p_x"
        ).efermi == approx(fermi)

        # without spin polarization
        energies_nonspin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        pdos_f_2s = [0.00000, 0.00320, 0.00000, 0.00017, 0.00000, 0.00060]
        pdos_f_2py = [0.00000, 0.00322, 0.00000, 0.51635, 0.00000, 0.00037]
        pdos_f_2pz = [0.00000, 0.00322, 0.00000, 0.51636, 0.00000, 0.00037]
        pdos_f_2px = [0.00000, 0.00322, 0.00000, 0.51634, 0.00000, 0.00037]

        assert (
            self.LobsterCompleteDOS_nonspin.get_site_orbital_dos(site=self.structure[0], orbital="2s").energies.tolist()
            == energies_nonspin
        )
        assert self.LobsterCompleteDOS_nonspin.get_site_orbital_dos(
            site=self.structure[0], orbital="2s"
        ).efermi == approx(fermi)
        assert (
            self.LobsterCompleteDOS_nonspin.get_site_orbital_dos(site=self.structure[0], orbital="2s")
            .densities[Spin.up]
            .tolist()
            == pdos_f_2s
        )

        assert (
            self.LobsterCompleteDOS_nonspin.get_site_orbital_dos(
                site=self.structure[0], orbital="2p_y"
            ).energies.tolist()
            == energies_nonspin
        )
        assert self.LobsterCompleteDOS_nonspin.get_site_orbital_dos(
            site=self.structure[0], orbital="2p_y"
        ).efermi == approx(fermi)
        assert (
            self.LobsterCompleteDOS_nonspin.get_site_orbital_dos(site=self.structure[0], orbital="2p_y")
            .densities[Spin.up]
            .tolist()
            == pdos_f_2py
        )

        assert (
            self.LobsterCompleteDOS_nonspin.get_site_orbital_dos(
                site=self.structure[0], orbital="2p_z"
            ).energies.tolist()
            == energies_nonspin
        )
        assert self.LobsterCompleteDOS_nonspin.get_site_orbital_dos(
            site=self.structure[0], orbital="2p_z"
        ).efermi == approx(fermi)
        assert (
            self.LobsterCompleteDOS_nonspin.get_site_orbital_dos(site=self.structure[0], orbital="2p_z")
            .densities[Spin.up]
            .tolist()
            == pdos_f_2pz
        )

        assert (
            self.LobsterCompleteDOS_nonspin.get_site_orbital_dos(
                site=self.structure[0], orbital="2p_x"
            ).energies.tolist()
            == energies_nonspin
        )
        assert self.LobsterCompleteDOS_nonspin.get_site_orbital_dos(
            site=self.structure[0], orbital="2p_x"
        ).efermi == approx(fermi)
        assert (
            self.LobsterCompleteDOS_nonspin.get_site_orbital_dos(site=self.structure[0], orbital="2p_x")
            .densities[Spin.up]
            .tolist()
            == pdos_f_2px
        )

    def test_get_site_t2g_eg_resolved_dos(self):
        # with spin polarization
        energies = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        efermi = 0.0
        pdos_mn_3dxy_up = [0.00000, 0.00001, 0.10301, 0.16070, 0.00070, 0.00060]
        pdos_mn_3dxy_down = [0.00000, 0.00000, 0.00380, 0.00996, 0.03012, 0.21890]
        pdos_mn_3dyz_up = [0.00000, 0.00001, 0.10301, 0.16070, 0.00070, 0.00060]
        pdos_mn_3dyz_down = [0.00000, 0.00000, 0.00380, 0.00996, 0.03012, 0.21890]
        pdos_mn_3dz2_up = [0.00000, 0.00001, 0.09608, 0.16941, 0.00028, 0.00028]
        pdos_mn_3dz2_down = [0.00000, 0.00000, 0.00433, 0.00539, 0.06000, 0.19427]
        pdos_mn_3dxz_up = [0.00000, 0.00001, 0.09746, 0.16767, 0.00036, 0.00034]
        pdos_mn_3dxz_down = [0.00000, 0.00000, 0.00422, 0.00630, 0.05402, 0.19919]
        pdos_mn_3dx2_up = [0.00000, 0.00001, 0.09330, 0.17289, 0.00011, 0.00015]
        pdos_mn_3dx2_down = [0.00000, 0.00000, 0.00454, 0.00356, 0.07195, 0.18442]

        pdos_mn_eg_up = (np.array(pdos_mn_3dx2_up) + np.array(pdos_mn_3dz2_up)).tolist()
        pdos_mn_eg_down = (np.array(pdos_mn_3dx2_down) + np.array(pdos_mn_3dz2_down)).tolist()
        pdos_mn_t2g_up = (np.array(pdos_mn_3dxy_up) + np.array(pdos_mn_3dxz_up) + np.array(pdos_mn_3dyz_up)).tolist()
        pdos_mn_t2g_down = (
            np.array(pdos_mn_3dxy_down) + np.array(pdos_mn_3dxz_down) + np.array(pdos_mn_3dyz_down)
        ).tolist()

        for iel, el in enumerate(
            self.LobsterCompleteDOS_MnO.get_site_t2g_eg_resolved_dos(self.structure_MnO[1])["e_g"]
            .densities[Spin.up]
            .tolist()
        ):
            assert el == approx(pdos_mn_eg_up[iel])
        for iel, el in enumerate(
            self.LobsterCompleteDOS_MnO.get_site_t2g_eg_resolved_dos(self.structure_MnO[1])["e_g"]
            .densities[Spin.down]
            .tolist()
        ):
            assert el == approx(pdos_mn_eg_down[iel])
        for iel, el in enumerate(
            self.LobsterCompleteDOS_MnO.get_site_t2g_eg_resolved_dos(self.structure_MnO[1])["t2g"]
            .densities[Spin.up]
            .tolist()
        ):
            assert el == approx(pdos_mn_t2g_up[iel])
        for iel, el in enumerate(
            self.LobsterCompleteDOS_MnO.get_site_t2g_eg_resolved_dos(self.structure_MnO[1])["t2g"]
            .densities[Spin.down]
            .tolist()
        ):
            assert el == approx(pdos_mn_t2g_down[iel])
        assert (
            energies
            == self.LobsterCompleteDOS_MnO.get_site_t2g_eg_resolved_dos(self.structure_MnO[1])["e_g"].energies.tolist()
        )
        assert (
            energies
            == self.LobsterCompleteDOS_MnO.get_site_t2g_eg_resolved_dos(self.structure_MnO[1])["t2g"].energies.tolist()
        )
        assert efermi == self.LobsterCompleteDOS_MnO.get_site_t2g_eg_resolved_dos(self.structure_MnO[1])["e_g"].efermi
        assert efermi == self.LobsterCompleteDOS_MnO.get_site_t2g_eg_resolved_dos(self.structure_MnO[1])["t2g"].efermi

        # without spin polarization
        energies_nonspin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        pdos_mn_3dxy = [0.00000, 0.00000, 0.02032, 0.16094, 0.33659, 0.01291]
        pdos_mn_3dyz = [0.00000, 0.00000, 0.02032, 0.16126, 0.33628, 0.01290]
        pdos_mn_3dz2 = [0.00000, 0.00000, 0.02591, 0.31460, 0.18658, 0.00509]
        pdos_mn_3dxz = [0.00000, 0.00000, 0.02484, 0.28501, 0.21541, 0.00663]
        pdos_mn_3dx2 = [0.00000, 0.00000, 0.02817, 0.37594, 0.12669, 0.00194]

        pdos_mn_eg = (np.array(pdos_mn_3dx2) + np.array(pdos_mn_3dz2)).tolist()
        pdos_mn_t2g = (np.array(pdos_mn_3dxy) + np.array(pdos_mn_3dxz) + np.array(pdos_mn_3dyz)).tolist()

        for iel, el in enumerate(
            self.LobsterCompleteDOS_MnO_nonspin.get_site_t2g_eg_resolved_dos(self.structure_MnO[1])["e_g"]
            .densities[Spin.up]
            .tolist()
        ):
            assert el == approx(pdos_mn_eg[iel])
        for iel, el in enumerate(
            self.LobsterCompleteDOS_MnO_nonspin.get_site_t2g_eg_resolved_dos(self.structure_MnO[1])["t2g"]
            .densities[Spin.up]
            .tolist()
        ):
            assert el == approx(pdos_mn_t2g[iel])
        assert (
            energies_nonspin
            == self.LobsterCompleteDOS_MnO_nonspin.get_site_t2g_eg_resolved_dos(self.structure_MnO[1])[
                "e_g"
            ].energies.tolist()
        )
        assert (
            energies_nonspin
            == self.LobsterCompleteDOS_MnO_nonspin.get_site_t2g_eg_resolved_dos(self.structure_MnO[1])[
                "t2g"
            ].energies.tolist()
        )
        assert (
            efermi
            == self.LobsterCompleteDOS_MnO_nonspin.get_site_t2g_eg_resolved_dos(self.structure_MnO[1])["e_g"].efermi
        )
        assert (
            efermi
            == self.LobsterCompleteDOS_MnO_nonspin.get_site_t2g_eg_resolved_dos(self.structure_MnO[1])["t2g"].efermi
        )

    def test_get_spd_dos(self):
        # with spin polarization
        energies_spin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        fermi = 0.0

        pdos_f_2s_up = [0.00000, 0.00159, 0.00000, 0.00011, 0.00000, 0.00069]
        pdos_f_2s_down = [0.00000, 0.00159, 0.00000, 0.00011, 0.00000, 0.00069]
        pdos_f_2py_up = [0.00000, 0.00160, 0.00000, 0.25801, 0.00000, 0.00029]
        pdos_f_2py_down = [0.00000, 0.00161, 0.00000, 0.25819, 0.00000, 0.00029]
        pdos_f_2pz_up = [0.00000, 0.00161, 0.00000, 0.25823, 0.00000, 0.00029]
        pdos_f_2pz_down = [0.00000, 0.00160, 0.00000, 0.25795, 0.00000, 0.00029]
        pdos_f_2px_up = [0.00000, 0.00160, 0.00000, 0.25805, 0.00000, 0.00029]
        pdos_f_2px_down = [0.00000, 0.00161, 0.00000, 0.25814, 0.00000, 0.00029]

        pdos_k_3s_up = [0.00000, 0.00000, 0.00000, 0.00008, 0.00000, 0.00007]
        pdos_k_3s_down = [0.00000, 0.00000, 0.00000, 0.00008, 0.00000, 0.00007]
        pdos_k_4s_up = [0.00000, 0.00018, 0.00000, 0.02035, 0.00000, 0.02411]
        pdos_k_4s_down = [0.00000, 0.00018, 0.00000, 0.02036, 0.00000, 0.02420]
        pdos_k_3py_up = [0.00000, 0.26447, 0.00000, 0.00172, 0.00000, 0.00000]
        pdos_k_3py_down = [0.00000, 0.26446, 0.00000, 0.00172, 0.00000, 0.00000]
        pdos_k_3pz_up = [0.00000, 0.26446, 0.00000, 0.00172, 0.00000, 0.00000]
        pdos_k_3pz_down = [0.00000, 0.26447, 0.00000, 0.00172, 0.00000, 0.00000]
        pdos_k_3px_up = [0.00000, 0.26447, 0.00000, 0.00172, 0.00000, 0.00000]
        pdos_k_3px_down = [0.00000, 0.26446, 0.00000, 0.00172, 0.00000, 0.00000]

        pdos_s_up = (np.array(pdos_f_2s_up) + np.array(pdos_k_3s_up) + np.array(pdos_k_4s_up)).tolist()
        pdos_s_down = (np.array(pdos_f_2s_down) + np.array(pdos_k_3s_down) + np.array(pdos_k_4s_down)).tolist()
        pdos_p_up = (
            np.array(pdos_f_2py_up)
            + np.array(pdos_f_2pz_up)
            + np.array(pdos_f_2px_up)
            + np.array(pdos_k_3py_up)
            + np.array(pdos_k_3pz_up)
            + np.array(pdos_k_3px_up)
        ).tolist()
        pdos_p_down = (
            np.array(pdos_f_2py_down)
            + np.array(pdos_f_2pz_down)
            + np.array(pdos_f_2px_down)
            + np.array(pdos_k_3py_down)
            + np.array(pdos_k_3pz_down)
            + np.array(pdos_k_3px_down)
        ).tolist()
        assert self.LobsterCompleteDOS_spin.get_spd_dos()[OrbitalType(0)].energies.tolist() == energies_spin
        assert self.LobsterCompleteDOS_spin.get_spd_dos()[OrbitalType(0)].efermi == fermi

        for ilistel, listel in enumerate(
            self.LobsterCompleteDOS_spin.get_spd_dos()[OrbitalType(0)].densities[Spin.up].tolist()
        ):
            assert listel == approx(pdos_s_up[ilistel])
        for ilistel, listel in enumerate(
            self.LobsterCompleteDOS_spin.get_spd_dos()[OrbitalType(0)].densities[Spin.down].tolist()
        ):
            assert listel == approx(pdos_s_down[ilistel])
        for ilistel, listel in enumerate(
            self.LobsterCompleteDOS_spin.get_spd_dos()[OrbitalType(1)].densities[Spin.up].tolist()
        ):
            assert listel == approx(pdos_p_up[ilistel])
        for ilistel, listel in enumerate(
            self.LobsterCompleteDOS_spin.get_spd_dos()[OrbitalType(1)].densities[Spin.down].tolist()
        ):
            assert listel == approx(pdos_p_down[ilistel])

        # without spin polarization
        energies_nonspin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        pdos_f_2s = [0.00000, 0.00320, 0.00000, 0.00017, 0.00000, 0.00060]
        pdos_f_2py = [0.00000, 0.00322, 0.00000, 0.51635, 0.00000, 0.00037]
        pdos_f_2pz = [0.00000, 0.00322, 0.00000, 0.51636, 0.00000, 0.00037]
        pdos_f_2px = [0.00000, 0.00322, 0.00000, 0.51634, 0.00000, 0.00037]

        pdos_k_3s = [0.00000, 0.00000, 0.00000, 0.00005, 0.00000, 0.00004]

        pdos_k_4s = [0.00000, 0.00040, 0.00000, 0.04039, 0.00000, 0.02241]

        pdos_k_3py = [0.00000, 0.52891, 0.00000, 0.00345, 0.00000, 0.00000]
        pdos_k_3pz = [0.00000, 0.52891, 0.00000, 0.00345, 0.00000, 0.00000]
        pdos_k_3px = [0.00000, 0.52891, 0.00000, 0.00345, 0.00000, 0.00000]

        pdos_s = (np.array(pdos_f_2s) + np.array(pdos_k_3s) + np.array(pdos_k_4s)).tolist()
        pdos_p = (
            np.array(pdos_f_2py)
            + np.array(pdos_f_2pz)
            + np.array(pdos_f_2px)
            + np.array(pdos_k_3py)
            + np.array(pdos_k_3pz)
            + np.array(pdos_k_3px)
        ).tolist()
        assert self.LobsterCompleteDOS_nonspin.get_spd_dos()[OrbitalType(0)].energies.tolist() == energies_nonspin

        for ilistel, listel in enumerate(
            self.LobsterCompleteDOS_nonspin.get_spd_dos()[OrbitalType(0)].densities[Spin.up].tolist()
        ):
            assert listel == approx(pdos_s[ilistel])
        for ilistel, listel in enumerate(
            self.LobsterCompleteDOS_nonspin.get_spd_dos()[OrbitalType(1)].densities[Spin.up].tolist()
        ):
            assert listel == approx(pdos_p[ilistel])

    def test_get_element_spd_dos(self):
        # with spin polarization
        energies_spin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        fermi = 0.0

        pdos_f_2s_up = [0.00000, 0.00159, 0.00000, 0.00011, 0.00000, 0.00069]
        pdos_f_2s_down = [0.00000, 0.00159, 0.00000, 0.00011, 0.00000, 0.00069]
        pdos_f_2py_up = [0.00000, 0.00160, 0.00000, 0.25801, 0.00000, 0.00029]
        pdos_f_2py_down = [0.00000, 0.00161, 0.00000, 0.25819, 0.00000, 0.00029]
        pdos_f_2pz_up = [0.00000, 0.00161, 0.00000, 0.25823, 0.00000, 0.00029]
        pdos_f_2pz_down = [0.00000, 0.00160, 0.00000, 0.25795, 0.00000, 0.00029]
        pdos_f_2px_up = [0.00000, 0.00160, 0.00000, 0.25805, 0.00000, 0.00029]
        pdos_f_2px_down = [0.00000, 0.00161, 0.00000, 0.25814, 0.00000, 0.00029]

        assert (
            self.LobsterCompleteDOS_spin.get_element_spd_dos(el=Element("F"))[OrbitalType(0)].energies.tolist()
            == energies_spin
        )

        assert (
            self.LobsterCompleteDOS_spin.get_element_spd_dos(el=Element("F"))[OrbitalType(0)]
            .densities[Spin.up]
            .tolist()
            == pdos_f_2s_up
        )
        assert (
            self.LobsterCompleteDOS_spin.get_element_spd_dos(el=Element("F"))[OrbitalType(0)]
            .densities[Spin.down]
            .tolist()
            == pdos_f_2s_down
        )

        for ilistel, listel in enumerate(
            self.LobsterCompleteDOS_spin.get_element_spd_dos(el=Element("F"))[OrbitalType(1)]
            .densities[Spin.up]
            .tolist()
        ):
            assert listel == approx(
                (np.array(pdos_f_2px_up) + np.array(pdos_f_2py_up) + np.array(pdos_f_2pz_up)).tolist()[ilistel]
            )

        for ilistel, listel in enumerate(
            self.LobsterCompleteDOS_spin.get_element_spd_dos(el=Element("F"))[OrbitalType(1)]
            .densities[Spin.down]
            .tolist()
        ):
            assert listel == approx(
                (np.array(pdos_f_2px_down) + np.array(pdos_f_2py_down) + np.array(pdos_f_2pz_down)).tolist()[ilistel]
            )

        assert self.LobsterCompleteDOS_spin.get_element_spd_dos(el=Element("F"))[OrbitalType(0)].efermi == approx(fermi)

        # without spin polarization
        energies_nonspin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        efermi = 0.0
        pdos_f_2s = [0.00000, 0.00320, 0.00000, 0.00017, 0.00000, 0.00060]
        pdos_f_2py = [0.00000, 0.00322, 0.00000, 0.51635, 0.00000, 0.00037]
        pdos_f_2pz = [0.00000, 0.00322, 0.00000, 0.51636, 0.00000, 0.00037]
        pdos_f_2px = [0.00000, 0.00322, 0.00000, 0.51634, 0.00000, 0.00037]

        assert (
            self.LobsterCompleteDOS_nonspin.get_element_spd_dos(el=Element("F"))[OrbitalType(0)].energies.tolist()
            == energies_nonspin
        )

        assert (
            self.LobsterCompleteDOS_nonspin.get_element_spd_dos(el=Element("F"))[OrbitalType(0)]
            .densities[Spin.up]
            .tolist()
            == pdos_f_2s
        )

        for ilistel, listel in enumerate(
            self.LobsterCompleteDOS_nonspin.get_element_spd_dos(el=Element("F"))[OrbitalType(1)]
            .densities[Spin.up]
            .tolist()
        ):
            assert listel == approx(
                (np.array(pdos_f_2px) + np.array(pdos_f_2py) + np.array(pdos_f_2pz)).tolist()[ilistel]
            )

        assert self.LobsterCompleteDOS_nonspin.get_element_spd_dos(el=Element("F"))[OrbitalType(0)].efermi == approx(
            efermi
        )
