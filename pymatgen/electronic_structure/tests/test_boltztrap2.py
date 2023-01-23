# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import annotations

import os
import unittest
import warnings

import numpy as np
from monty.serialization import loadfn
from pytest import approx

from pymatgen.electronic_structure.core import OrbitalType, Spin
from pymatgen.io.vasp import Vasprun
from pymatgen.util.testing import PymatgenTest

try:
    from pymatgen.electronic_structure.boltztrap2 import (
        BandstructureLoader,
        BztInterpolator,
        BztPlotter,
        BztTransportProperties,
        VasprunBSLoader,
        VasprunLoader,
    )

    BOLTZTRAP2_PRESENT = True
except Exception:
    BOLTZTRAP2_PRESENT = False

# BOLTZTRAP2_PRESENT = False

test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "boltztrap2")


vrunfile = os.path.join(test_dir, "vasprun.xml")
vrun = Vasprun(vrunfile, parse_projected_eigen=True)

vrunfile_sp = os.path.join(test_dir, "vasprun_spin.xml")
vrun_sp = Vasprun(vrunfile_sp, parse_projected_eigen=True)
bs = loadfn(os.path.join(test_dir, "PbTe_bandstructure.json"))
bs_sp = loadfn(os.path.join(test_dir, "N2_bandstructure.json"))

bztinterp_fn = os.path.join(test_dir, "bztInterp.json.gz")
bzttransp_fn = os.path.join(test_dir, "bztTranspProps.json.gz")


@unittest.skipIf(not BOLTZTRAP2_PRESENT, "No boltztrap2, skipping tests...")
class VasprunBSLoaderTest(unittest.TestCase):
    def setUp(self):
        self.loader = VasprunBSLoader(vrun)
        assert self.loader is not None
        self.loader = VasprunBSLoader(bs, vrun.final_structure)
        assert self.loader is not None
        self.loader = VasprunBSLoader.from_file(vrunfile)
        assert self.loader is not None

        warnings.simplefilter("ignore")

        self.loader_sp = VasprunBSLoader(vrun_sp)
        assert self.loader_sp is not None
        self.loader_sp = VasprunBSLoader(bs_sp, vrun_sp.final_structure)
        assert self.loader_sp is not None
        self.loader_sp = VasprunBSLoader.from_file(vrunfile_sp)
        assert self.loader_sp is not None

        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_properties(self):
        assert self.loader.is_spin_polarized is False
        assert self.loader.fermi == approx(0.185266535678, abs=1e-5)
        assert self.loader.structure.lattice.a == approx(4.64303565932548, abs=1e-5)
        assert self.loader.nelect_all == 20.0
        assert self.loader_sp.nelect_all == 10.0

        assert self.loader.ebands_all.shape == (20, 120)
        assert self.loader.ebands_all[10, 100] == approx(0.2708057, abs=1e-5)
        assert len(self.loader.proj_all) == 1
        assert self.loader.proj_all[Spin.up].shape == (120, 20, 2, 9)

        assert self.loader_sp.is_spin_polarized is True
        assert self.loader_sp.ebands_all.shape == (24, 198)
        assert self.loader_sp.ebands_all[10, 100] == approx(0.2543788, abs=1e-4)
        assert self.loader_sp.ebands_all[22, 100] == approx(0.2494617, abs=1e-4)
        assert len(self.loader_sp.proj_all) == 2
        assert self.loader_sp.proj_all[Spin.down].shape == (198, 12, 2, 9)

    def test_get_volume(self):
        assert self.loader.get_volume() == approx(477.6256714925874, abs=1e-5)


@unittest.skipIf(not BOLTZTRAP2_PRESENT, "No boltztrap2, skipping tests...")
class BandstructureLoaderTest(unittest.TestCase):
    def setUp(self):
        self.loader = BandstructureLoader(bs, vrun.structures[-1])
        assert self.loader is not None

        self.loader_sp = BandstructureLoader(bs_sp, vrun_sp.structures[-1])
        assert self.loader_sp is not None
        assert self.loader_sp.ebands_all.shape == (24, 198)

        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_properties(self):
        assert self.loader.ebands_all.shape == (20, 120)
        assert self.loader.fermi == approx(0.185266535678, abs=1e-5)
        assert self.loader.structure.lattice.a == approx(4.64303565932548, abs=1e-5)

    def test_get_volume(self):
        assert self.loader.get_volume() == approx(477.6256714925874, abs=1e-5)

    # def test_set_upper_lower_bands(self):
    #     min_bnd = min(self.loader_sp_up.ebands.min(),
    #                   self.loader_sp_dn.ebands.min())
    #     max_bnd = max(self.loader_sp_up.ebands.max(),
    #                   self.loader_sp_dn.ebands.max())
    #     self.loader_sp_up.set_upper_lower_bands(min_bnd, max_bnd)
    #     self.loader_sp_dn.set_upper_lower_bands(min_bnd, max_bnd)
    #     self.assertTupleEqual(self.loader_sp_up.ebands.shape, (14, 198))
    #     self.assertTupleEqual(self.loader_sp_dn.ebands.shape, (14, 198))


@unittest.skipIf(not BOLTZTRAP2_PRESENT, "No boltztrap2, skipping tests...")
class VasprunLoaderTest(unittest.TestCase):
    def setUp(self):
        self.loader = VasprunLoader(vrun)
        assert self.loader.proj.shape == (120, 20, 2, 9)
        assert self.loader is not None
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_properties(self):
        assert self.loader.ebands.shape == (20, 120)
        assert self.loader.fermi == approx(0.185266535678, abs=1e-5)
        assert self.loader.structure.lattice.a == approx(4.64303565932548, abs=1e-5)

    def test_get_volume(self):
        assert self.loader.get_volume() == approx(477.6256714925874, abs=1e-5)

    def test_from_file(self):
        self.loader = VasprunLoader().from_file(vrunfile)
        assert self.loader is not None


@unittest.skipIf(not BOLTZTRAP2_PRESENT, "No boltztrap2, skipping tests...")
class BztInterpolatorTest(unittest.TestCase):
    def setUp(self):
        self.loader = VasprunBSLoader(vrun)
        self.bztInterp = BztInterpolator(self.loader, lpfac=2)
        assert self.bztInterp is not None
        self.bztInterp = BztInterpolator(self.loader, lpfac=2, save_bztInterp=True, fname=bztinterp_fn)
        assert self.bztInterp is not None
        self.bztInterp = BztInterpolator(self.loader, load_bztInterp=True, fname=bztinterp_fn)
        assert self.bztInterp is not None

        warnings.simplefilter("ignore")

        self.loader_sp = VasprunBSLoader(vrun_sp)
        self.bztInterp_sp = BztInterpolator(self.loader_sp, lpfac=2)
        assert self.bztInterp_sp is not None
        self.bztInterp_sp = BztInterpolator(self.loader_sp, lpfac=2, save_bztInterp=True, fname=bztinterp_fn)
        assert self.bztInterp_sp is not None
        self.bztInterp_sp = BztInterpolator(self.loader_sp, lpfac=2, load_bztInterp=True, fname=bztinterp_fn)
        assert self.bztInterp_sp is not None

        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_properties(self):
        assert self.bztInterp.cband.shape == (6, 3, 3, 3, 29791)
        assert self.bztInterp.eband.shape == (6, 29791)
        assert self.bztInterp.coeffs.shape == (6, 322)
        assert self.bztInterp.data.nelect == 6.0
        assert self.bztInterp.data.nelect_all == 20.0
        assert self.bztInterp.data.ebands.shape == (6, 120)

        assert self.bztInterp_sp.cband.shape == (10, 3, 3, 3, 23275)
        assert self.bztInterp_sp.eband.shape == (10, 23275)
        assert self.bztInterp_sp.coeffs.shape == (10, 519)
        assert self.bztInterp_sp.data.nelect == 6.0
        assert self.bztInterp_sp.data.nelect_all == 10.0
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


@unittest.skipIf(not BOLTZTRAP2_PRESENT, "No boltztrap2, skipping tests...")
class BztTransportPropertiesTest(unittest.TestCase):
    def setUp(self):
        loader = VasprunBSLoader(vrun)
        bztInterp = BztInterpolator(loader, lpfac=2)
        self.bztTransp = BztTransportProperties(bztInterp, temp_r=np.arange(300, 600, 100))
        assert self.bztTransp is not None
        warnings.simplefilter("ignore")

        self.bztTransp = BztTransportProperties(
            bztInterp, doping=10.0 ** np.arange(20, 22), temp_r=np.arange(300, 600, 100)
        )
        assert self.bztTransp is not None
        assert self.bztTransp.contain_props_doping is True

        warnings.simplefilter("ignore")

        bztInterp = BztInterpolator(loader, lpfac=2)
        self.bztTransp = BztTransportProperties(
            bztInterp,
            temp_r=np.arange(300, 600, 100),
            save_bztTranspProps=True,
            fname=bzttransp_fn,
        )
        assert self.bztTransp is not None
        warnings.simplefilter("ignore")

        bztInterp = BztInterpolator(loader, lpfac=2)
        self.bztTransp = BztTransportProperties(bztInterp, load_bztTranspProps=True, fname=bzttransp_fn)
        assert self.bztTransp is not None
        warnings.simplefilter("ignore")

        loader_sp = VasprunBSLoader(vrun_sp)
        bztInterp_sp = BztInterpolator(loader_sp, lpfac=2)
        self.bztTransp_sp = BztTransportProperties(bztInterp_sp, temp_r=np.arange(300, 600, 100))
        assert self.bztTransp_sp is not None
        warnings.simplefilter("ignore")

        bztInterp_sp = BztInterpolator(loader_sp, lpfac=2)
        self.bztTransp_sp = BztTransportProperties(
            bztInterp_sp,
            temp_r=np.arange(300, 600, 100),
            save_bztTranspProps=True,
            fname=bzttransp_fn,
        )
        assert self.bztTransp_sp is not None
        warnings.simplefilter("ignore")

        bztInterp_sp = BztInterpolator(loader_sp, lpfac=2)
        self.bztTransp_sp = BztTransportProperties(bztInterp_sp, load_bztTranspProps=True, fname=bzttransp_fn)
        assert self.bztTransp_sp is not None
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_properties(self):
        for p in [
            self.bztTransp.Conductivity_mu,
            self.bztTransp.Seebeck_mu,
            self.bztTransp.Kappa_mu,
            self.bztTransp.Effective_mass_mu,
            self.bztTransp.Power_Factor_mu,
        ]:
            assert p.shape == (3, 3686, 3, 3)

        for p in [
            self.bztTransp.Carrier_conc_mu,
            self.bztTransp.Hall_carrier_conc_trace_mu,
        ]:
            assert p.shape == (3, 3686)

        for p in [
            self.bztTransp_sp.Conductivity_mu,
            self.bztTransp_sp.Seebeck_mu,
            self.bztTransp_sp.Kappa_mu,
            self.bztTransp_sp.Effective_mass_mu,
            self.bztTransp_sp.Power_Factor_mu,
        ]:
            assert p.shape == (3, 3252, 3, 3)

        for p in [
            self.bztTransp_sp.Carrier_conc_mu,
            self.bztTransp_sp.Hall_carrier_conc_trace_mu,
        ]:
            assert p.shape == (3, 3252)

    def test_compute_properties_doping(self):
        self.bztTransp.compute_properties_doping(doping=10.0 ** np.arange(20, 22))
        for p in [
            self.bztTransp.Conductivity_doping,
            self.bztTransp.Seebeck_doping,
            self.bztTransp.Kappa_doping,
            self.bztTransp.Effective_mass_doping,
            self.bztTransp.Power_Factor_doping,
        ]:
            assert p["n"].shape == (3, 2, 3, 3)
            assert self.bztTransp.contain_props_doping is True

        self.bztTransp_sp.compute_properties_doping(doping=10.0 ** np.arange(20, 22))
        for p in [
            self.bztTransp_sp.Conductivity_doping,
            self.bztTransp_sp.Seebeck_doping,
            self.bztTransp_sp.Kappa_doping,
            self.bztTransp_sp.Effective_mass_doping,
            self.bztTransp_sp.Power_Factor_doping,
        ]:
            assert p["n"].shape == (3, 2, 3, 3)
            assert self.bztTransp_sp.contain_props_doping is True


@unittest.skipIf(not BOLTZTRAP2_PRESENT, "No boltztrap2, skipping tests...")
class BztPlotterTest(unittest.TestCase):
    def test_plot(self):
        loader = VasprunBSLoader(vrun)
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


if __name__ == "__main__":
    unittest.main()
