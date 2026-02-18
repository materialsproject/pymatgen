from __future__ import annotations

import shutil

import numpy as np
import pytest
from monty.serialization import loadfn
from pytest import approx

from pymatgen.electronic_structure.core import OrbitalType, Spin
from pymatgen.io.vasp import Vasprun
from pymatgen.util.testing import TEST_FILES_DIR

try:
    from pymatgen.electronic_structure.boltztrap2 import (
        BandstructureLoader,
        BztInterpolator,
        BztPlotter,
        BztTransportProperties,
        VasprunBSLoader,
        VasprunLoader,
    )

except ImportError:
    pytest.skip("No boltztrap2", allow_module_level=True)

TEST_DIR = f"{TEST_FILES_DIR}/electronic_structure/boltztrap2"

VASP_RUN_FILE = f"{TEST_DIR}/vasprun.xml"
VASP_RUN = Vasprun(VASP_RUN_FILE, parse_projected_eigen=True)

VASP_RUN_FILE_SPIN = f"{TEST_DIR}/vasprun_spin.xml"
VASP_RUN_SPIN = Vasprun(VASP_RUN_FILE_SPIN, parse_projected_eigen=True)

BAND_STRUCT = loadfn(f"{TEST_DIR}/PbTe_bandstructure.json")
BAND_STRUCT_SPIN = loadfn(f"{TEST_DIR}/N2_bandstructure.json")

BZT_INTERP_FN = f"{TEST_DIR}/bztInterp.json.gz"
BZT_TRANSP_FN = f"{TEST_DIR}/bztTranspProps.json.gz"


class TestVasprunBSLoader:
    def setup_method(self):
        self.loader = VasprunBSLoader(VASP_RUN)
        assert self.loader is not None
        self.loader = VasprunBSLoader(BAND_STRUCT, VASP_RUN.final_structure)
        assert self.loader is not None
        self.loader = VasprunBSLoader.from_file(VASP_RUN_FILE)
        assert self.loader is not None

        self.loader_sp = VasprunBSLoader(VASP_RUN_SPIN)
        assert self.loader_sp is not None
        self.loader_sp = VasprunBSLoader(BAND_STRUCT_SPIN, VASP_RUN_SPIN.final_structure)
        assert self.loader_sp is not None
        self.loader_sp = VasprunBSLoader.from_file(VASP_RUN_FILE_SPIN)
        assert self.loader_sp is not None

    def test_properties(self):
        assert self.loader.is_spin_polarized is False
        assert self.loader.fermi == approx(0.185266535678, abs=1e-5)
        assert self.loader.structure.lattice.a == approx(4.64303565932548, abs=1e-5)
        assert self.loader.nelect_all == approx(20.0)
        assert self.loader_sp.nelect_all == approx(10.0)

        assert self.loader.ebands_all.shape == (20, 120)
        assert self.loader.ebands_all[10, 100] == approx(0.2708057, abs=1e-5)
        assert len(self.loader.proj_all) == 1
        assert self.loader.proj_all[Spin.up].shape == (120, 20, 2, 9)

        assert self.loader_sp.is_spin_polarized
        assert self.loader_sp.ebands_all.shape == (24, 198)
        assert self.loader_sp.ebands_all[10, 100] == approx(0.2543788, abs=1e-4)
        assert self.loader_sp.ebands_all[22, 100] == approx(0.2494617, abs=1e-4)
        assert len(self.loader_sp.proj_all) == 2
        assert self.loader_sp.proj_all[Spin.down].shape == (198, 12, 2, 9)

    def test_get_volume(self):
        assert self.loader.get_volume() == approx(477.6256714925874, abs=1e-5)


@pytest.mark.filterwarnings("ignore:BandstructureLoader is deprecated:DeprecationWarning")
class TestBandstructureLoader:
    def setup_method(self):
        self.loader = BandstructureLoader(BAND_STRUCT, VASP_RUN.structures[-1])
        assert self.loader is not None

        self.loader_sp = BandstructureLoader(BAND_STRUCT_SPIN, VASP_RUN_SPIN.structures[-1])
        assert self.loader_sp is not None
        assert self.loader_sp.ebands_all.shape == (24, 198)

    def test_properties(self):
        assert self.loader.ebands_all.shape == (20, 120)
        assert self.loader.fermi == approx(0.185266535678, abs=1e-5)
        assert self.loader.structure.lattice.a == approx(4.64303565932548, abs=1e-5)

    def test_get_volume(self):
        assert self.loader.get_volume() == approx(477.6256714925874, abs=1e-5)

    @pytest.mark.xfail(reason="TODO: need someone to fix this")
    def test_set_upper_lower_bands(self):
        min_bnd = min(self.loader_sp_up.ebands.min(), self.loader_sp_dn.ebands.min())
        max_bnd = max(self.loader_sp_up.ebands.max(), self.loader_sp_dn.ebands.max())
        self.loader_sp_up.set_upper_lower_bands(min_bnd, max_bnd)
        self.loader_sp_dn.set_upper_lower_bands(min_bnd, max_bnd)
        assert self.loader_sp_up.ebands.shape == (14, 198)
        assert self.loader_sp_dn.ebands.shape == (14, 198)


@pytest.mark.filterwarnings("ignore:VasprunLoader is deprecated:DeprecationWarning")
class TestVasprunLoader:
    def setup_method(self):
        self.loader = VasprunLoader(VASP_RUN)
        assert self.loader.proj.shape == (120, 20, 2, 9)
        assert self.loader is not None

    def test_properties(self):
        assert self.loader.ebands.shape == (20, 120)
        assert self.loader.fermi == approx(0.185266535678, abs=1e-5)
        assert self.loader.structure.lattice.a == approx(4.64303565932548, abs=1e-5)

    def test_get_volume(self):
        assert self.loader.get_volume() == approx(477.6256714925874, abs=1e-5)

    def test_from_file(self):
        self.loader = VasprunLoader().from_file(VASP_RUN_FILE)
        assert self.loader is not None


class TestBztInterpolator:
    def setup_method(self):
        shutil.copy(BZT_INTERP_FN, ".")

        loader = VasprunBSLoader(VASP_RUN)
        self.bztInterp = BztInterpolator(loader, lpfac=2)
        assert self.bztInterp is not None
        self.bztInterp = BztInterpolator(loader, lpfac=2, save_bztInterp=True)
        assert self.bztInterp is not None
        self.bztInterp = BztInterpolator(loader, load_bztInterp=True)
        assert self.bztInterp is not None

        loader_sp = VasprunBSLoader(VASP_RUN_SPIN)
        self.bztInterp_sp = BztInterpolator(loader_sp, lpfac=2)
        assert self.bztInterp_sp is not None
        self.bztInterp_sp = BztInterpolator(loader_sp, lpfac=2, save_bztInterp=True)
        assert self.bztInterp_sp is not None
        self.bztInterp_sp = BztInterpolator(loader_sp, lpfac=2, load_bztInterp=True)
        assert self.bztInterp_sp is not None

    def test_properties(self):
        assert self.bztInterp.cband.shape == (6, 3, 3, 3, 29791)
        assert self.bztInterp.eband.shape == (6, 29791)
        assert self.bztInterp.coeffs.shape == (6, 322)
        assert self.bztInterp.data.nelect == approx(6.0)
        assert self.bztInterp.data.nelect_all == approx(20.0)
        assert self.bztInterp.data.ebands.shape == (6, 120)

        assert self.bztInterp_sp.cband.shape == (10, 3, 3, 3, 23275)
        assert self.bztInterp_sp.eband.shape == (10, 23275)
        assert self.bztInterp_sp.coeffs.shape == (10, 519)
        assert self.bztInterp_sp.data.nelect == approx(6.0)
        assert self.bztInterp_sp.data.nelect_all == approx(10.0)
        assert self.bztInterp_sp.data.ebands.shape == (10, 198)

    def test_get_band_structure(self):
        sbs = self.bztInterp.get_band_structure()
        assert sbs is not None
        assert sbs.bands[Spin.up].shape == (6, 137)
        kpaths = [["L", "K"]]
        kp_lbl = {"L": np.array([0.5, 0.5, 0.5]), "K": np.array([0.375, 0.375, 0.75])}
        sbs = self.bztInterp.get_band_structure(kpaths, kp_lbl)
        assert sbs is not None
        assert sbs.bands[Spin.up].shape == (6, 20)

        sbs = self.bztInterp_sp.get_band_structure()
        assert sbs is not None
        assert sbs.bands[Spin.up].shape == (6, 143)
        assert sbs.bands[Spin.down].shape == (4, 143)

    def test_tot_dos(self):
        tot_dos = self.bztInterp.get_dos(T=200, npts_mu=100)
        assert tot_dos is not None
        assert len(tot_dos.energies) == 100
        assert tot_dos.densities[Spin.up][0] == approx(1.35371715, abs=1e-5)

        tot_dos = self.bztInterp_sp.get_dos(T=200, npts_mu=100)
        assert tot_dos is not None
        assert len(tot_dos.energies) == 100
        assert tot_dos.densities[Spin.up][75] == approx(88.034456, abs=1e-5)
        assert tot_dos.densities[Spin.down][75] == approx(41.421367, abs=1e-5)

    def test_tot_proj_dos(self):
        tot_proj_dos = self.bztInterp.get_dos(partial_dos=True, T=200, npts_mu=100)
        assert tot_proj_dos is not None
        assert len(tot_proj_dos.get_spd_dos().values()) == 3
        pdos = tot_proj_dos.get_spd_dos()[OrbitalType.s].densities[Spin.up][75]
        assert pdos == approx(2490.169396, abs=1e-5)

        tot_proj_dos = self.bztInterp_sp.get_dos(partial_dos=True, T=200, npts_mu=100)
        assert tot_proj_dos is not None
        assert len(tot_proj_dos.get_spd_dos().values()) == 3
        pdos = tot_proj_dos.get_spd_dos()[OrbitalType.s].densities[Spin.up][75]
        assert pdos == approx(166.4933305, abs=1e-5)
        pdos = tot_proj_dos.get_spd_dos()[OrbitalType.s].densities[Spin.down][75]
        assert pdos == approx(272.194174, abs=1e-5)


class TestBztTransportProperties:
    def setup_method(self):
        shutil.copy(BZT_TRANSP_FN, ".")

        # non spin polarized
        loader = VasprunBSLoader(VASP_RUN)
        bztInterp = BztInterpolator(loader, lpfac=2)

        self.bztTransp = BztTransportProperties(
            bztInterp,
            temp_r=np.arange(300, 600, 100),
            save_bztTranspProps=True,
        )
        assert self.bztTransp is not None

        bztInterp = BztInterpolator(loader, lpfac=2)
        self.bztTransp = BztTransportProperties(
            bztInterp,
            load_bztTranspProps=True,
        )
        assert self.bztTransp is not None

        # spin polarized
        loader_sp = VasprunBSLoader(VASP_RUN_SPIN)
        bztInterp_sp = BztInterpolator(loader_sp, lpfac=2)

        self.bztTransp_sp = BztTransportProperties(
            bztInterp_sp,
            temp_r=np.arange(300, 600, 100),
            save_bztTranspProps=True,
        )
        assert self.bztTransp_sp is not None

        bztInterp_sp = BztInterpolator(loader_sp, lpfac=2)
        self.bztTransp_sp = BztTransportProperties(
            bztInterp_sp,
            load_bztTranspProps=True,
        )
        assert self.bztTransp_sp is not None

    def test_init(self):
        # non spin polarized
        loader = VasprunBSLoader(VASP_RUN)
        bztInterp = BztInterpolator(loader, lpfac=2)
        self.bztTransp = BztTransportProperties(bztInterp, temp_r=np.arange(300, 600, 100))
        assert self.bztTransp is not None

        self.bztTransp = BztTransportProperties(
            bztInterp,
            doping=10.0 ** np.arange(20, 22),
            temp_r=np.arange(300, 600, 100),
        )
        assert self.bztTransp is not None
        assert self.bztTransp.contain_props_doping

        # spin polarized
        loader_sp = VasprunBSLoader(VASP_RUN_SPIN)
        bztInterp_sp = BztInterpolator(loader_sp, lpfac=2)
        self.bztTransp_sp = BztTransportProperties(bztInterp_sp, temp_r=np.arange(300, 600, 100))
        assert self.bztTransp_sp is not None

    def test_properties(self):
        for p in (
            self.bztTransp.Conductivity_mu,
            self.bztTransp.Seebeck_mu,
            self.bztTransp.Kappa_mu,
            self.bztTransp.Effective_mass_mu,
            self.bztTransp.Power_Factor_mu,
        ):
            assert p.shape == (3, 3686, 3, 3)

        for p in (
            self.bztTransp.Carrier_conc_mu,
            self.bztTransp.Hall_carrier_conc_trace_mu,
        ):
            assert p.shape == (3, 3686)

        for p in (
            self.bztTransp_sp.Conductivity_mu,
            self.bztTransp_sp.Seebeck_mu,
            self.bztTransp_sp.Kappa_mu,
            self.bztTransp_sp.Effective_mass_mu,
            self.bztTransp_sp.Power_Factor_mu,
        ):
            assert p.shape == (3, 3252, 3, 3)

        for p in (
            self.bztTransp_sp.Carrier_conc_mu,
            self.bztTransp_sp.Hall_carrier_conc_trace_mu,
        ):
            assert p.shape == (3, 3252)

    def test_compute_properties_doping(self):
        self.bztTransp.compute_properties_doping(doping=10.0 ** np.arange(20, 22))
        for p in (
            self.bztTransp.Conductivity_doping,
            self.bztTransp.Seebeck_doping,
            self.bztTransp.Kappa_doping,
            self.bztTransp.Effective_mass_doping,
            self.bztTransp.Power_Factor_doping,
        ):
            assert p["n"].shape == (3, 2, 3, 3)
            assert self.bztTransp.contain_props_doping

        self.bztTransp_sp.compute_properties_doping(doping=10.0 ** np.arange(20, 22))
        for p in (
            self.bztTransp_sp.Conductivity_doping,
            self.bztTransp_sp.Seebeck_doping,
            self.bztTransp_sp.Kappa_doping,
            self.bztTransp_sp.Effective_mass_doping,
            self.bztTransp_sp.Power_Factor_doping,
        ):
            assert p["n"].shape == (3, 2, 3, 3)
            assert self.bztTransp_sp.contain_props_doping


class TestBztPlotter:
    def test_plot(self):
        loader = VasprunBSLoader(VASP_RUN)
        bztInterp = BztInterpolator(loader, lpfac=2)
        bztTransp = BztTransportProperties(bztInterp, temp_r=np.arange(300, 600, 100))
        self.bztPlotter = BztPlotter(bztTransp, bztInterp)
        assert self.bztPlotter is not None
        fig = self.bztPlotter.plot_props("S", "mu", "temp", temps=[300, 500])
        assert fig is not None

        fig = self.bztPlotter.plot_bands()
        assert fig is not None

        fig = self.bztPlotter.plot_dos()
        assert fig is not None
