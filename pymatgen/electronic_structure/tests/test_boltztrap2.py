# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest
import os
import warnings

from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.core import Spin, OrbitalType

import numpy as np
from monty.serialization import loadfn

try:
    from pymatgen.electronic_structure.boltztrap2 import VasprunBSLoader, \
        BandstructureLoader, VasprunLoader, BztInterpolator, \
        BztTransportProperties, BztPlotter

    BOLTZTRAP2_PRESENT = True
except Exception:
    BOLTZTRAP2_PRESENT = False

# BOLTZTRAP2_PRESENT = False

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files/boltztrap2/')

vrunfile = os.path.join(test_dir, 'vasprun.xml')
vrun = Vasprun(vrunfile, parse_projected_eigen=True)

vrunfile_sp = os.path.join(test_dir, 'vasprun_spin.xml')
vrun_sp = Vasprun(vrunfile_sp, parse_projected_eigen=True)
bs = loadfn(os.path.join(test_dir, "PbTe_bandstructure.json"))
bs_sp = loadfn(os.path.join(test_dir, "N2_bandstructure.json"))

bztinterp_fn = os.path.join(test_dir, "bztInterp.json.gz")
bzttransp_fn = os.path.join(test_dir, "bztTranspProps.json.gz")


@unittest.skipIf(not BOLTZTRAP2_PRESENT, "No boltztrap2, skipping tests...")
class VasprunBSLoaderTest(unittest.TestCase):
    def setUp(self):
        self.loader = VasprunBSLoader(vrun)
        self.assertIsNotNone(self.loader)
        self.loader = VasprunBSLoader(bs, vrun.final_structure)
        self.assertIsNotNone(self.loader)
        self.loader = VasprunBSLoader.from_file(vrunfile)
        self.assertIsNotNone(self.loader)

        warnings.simplefilter("ignore")

        self.loader_sp = VasprunBSLoader(vrun_sp)
        self.assertIsNotNone(self.loader_sp)
        self.loader_sp = VasprunBSLoader(bs_sp, vrun_sp.final_structure)
        self.assertIsNotNone(self.loader_sp)
        self.loader_sp = VasprunBSLoader.from_file(vrunfile_sp)
        self.assertIsNotNone(self.loader_sp)

        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_properties(self):
        self.assertEqual(self.loader.is_spin_polarized, False)
        self.assertAlmostEqual(self.loader.fermi, 0.185266535678, 5)
        self.assertAlmostEqual(self.loader.structure.lattice.a,
                               4.64303565932548, 5)
        self.assertEqual(self.loader.nelect_all, 20.0)
        self.assertEqual(self.loader_sp.nelect_all, 10.0)

        self.assertTupleEqual(self.loader.ebands_all.shape, (20, 120))
        self.assertAlmostEqual(self.loader.ebands_all[10, 100], 0.2708057, 5)
        self.assertEqual(len(self.loader.proj_all), 1)
        self.assertTupleEqual(self.loader.proj_all[Spin.up].shape,
                              (120, 20, 2, 9))

        self.assertEqual(self.loader_sp.is_spin_polarized, True)
        self.assertTupleEqual(self.loader_sp.ebands_all.shape, (24, 198))
        self.assertAlmostEqual(self.loader_sp.ebands_all[10, 100], 0.2543788,
                               4)
        self.assertAlmostEqual(self.loader_sp.ebands_all[22, 100], 0.2494617,
                               4)
        self.assertEqual(len(self.loader_sp.proj_all), 2)
        self.assertTupleEqual(self.loader_sp.proj_all[Spin.down].shape,
                              (198, 12, 2, 9))

    def test_get_volume(self):
        self.assertAlmostEqual(self.loader.get_volume(), 477.6256714925874, 5)


@unittest.skipIf(not BOLTZTRAP2_PRESENT, "No boltztrap2, skipping tests...")
class BandstructureLoaderTest(unittest.TestCase):
    def setUp(self):
        self.loader = BandstructureLoader(bs, vrun.structures[-1])
        self.assertIsNotNone(self.loader)

        self.loader_sp = BandstructureLoader(bs_sp, vrun_sp.structures[-1])
        self.assertIsNotNone(self.loader_sp)
        self.assertTupleEqual(self.loader_sp.ebands_all.shape, (24, 198))

        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_properties(self):
        self.assertTupleEqual(self.loader.ebands_all.shape, (20, 120))
        self.assertAlmostEqual(self.loader.fermi, 0.185266535678, 5)
        self.assertAlmostEqual(self.loader.structure.lattice.a,
                               4.64303565932548, 5)

    def test_get_volume(self):
        self.assertAlmostEqual(self.loader.get_volume(), 477.6256714925874, 5)

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
        self.assertTupleEqual(self.loader.proj.shape, (120, 20, 2, 9))
        self.assertIsNotNone(self.loader)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_properties(self):
        self.assertTupleEqual(self.loader.ebands.shape, (20, 120))
        self.assertAlmostEqual(self.loader.fermi, 0.185266535678, 5)
        self.assertAlmostEqual(self.loader.structure.lattice.a,
                               4.64303565932548, 5)

    def test_get_volume(self):
        self.assertAlmostEqual(self.loader.get_volume(), 477.6256714925874, 5)

    def test_from_file(self):
        self.loader = VasprunLoader().from_file(vrunfile)
        self.assertIsNotNone(self.loader)


@unittest.skipIf(not BOLTZTRAP2_PRESENT, "No boltztrap2, skipping tests...")
class BztInterpolatorTest(unittest.TestCase):
    def setUp(self):
        self.loader = VasprunBSLoader(vrun)
        self.bztInterp = BztInterpolator(self.loader, lpfac=2)
        self.assertIsNotNone(self.bztInterp)
        self.bztInterp = BztInterpolator(self.loader,
                                         lpfac=2,
                                         save_bztInterp=True,
                                         fname=bztinterp_fn)
        self.assertIsNotNone(self.bztInterp)
        self.bztInterp = BztInterpolator(self.loader,
                                         load_bztInterp=True,
                                         fname=bztinterp_fn)
        self.assertIsNotNone(self.bztInterp)

        warnings.simplefilter("ignore")

        self.loader_sp = VasprunBSLoader(vrun_sp)
        self.bztInterp_sp = BztInterpolator(self.loader_sp, lpfac=2)
        self.assertIsNotNone(self.bztInterp_sp)
        self.bztInterp_sp = BztInterpolator(self.loader_sp,
                                            lpfac=2,
                                            save_bztInterp=True,
                                            fname=bztinterp_fn)
        self.assertIsNotNone(self.bztInterp_sp)
        self.bztInterp_sp = BztInterpolator(self.loader_sp,
                                            lpfac=2,
                                            load_bztInterp=True,
                                            fname=bztinterp_fn)
        self.assertIsNotNone(self.bztInterp_sp)

        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_properties(self):
        self.assertTupleEqual(self.bztInterp.cband.shape, (6, 3, 3, 3, 29791))
        self.assertTupleEqual(self.bztInterp.eband.shape, (6, 29791))
        self.assertTupleEqual(self.bztInterp.coeffs.shape, (6, 322))
        self.assertEqual(self.bztInterp.data.nelect, 6.0)
        self.assertEqual(self.bztInterp.data.nelect_all, 20.0)
        self.assertTupleEqual(self.bztInterp.data.ebands.shape, (6, 120))

        self.assertTupleEqual(self.bztInterp_sp.cband.shape,
                              (10, 3, 3, 3, 23275))
        self.assertTupleEqual(self.bztInterp_sp.eband.shape, (10, 23275))
        self.assertTupleEqual(self.bztInterp_sp.coeffs.shape, (10, 519))
        self.assertEqual(self.bztInterp_sp.data.nelect, 6.0)
        self.assertEqual(self.bztInterp_sp.data.nelect_all, 10.0)
        self.assertTupleEqual(self.bztInterp_sp.data.ebands.shape, (10, 198))

    def test_get_band_structure(self):
        sbs = self.bztInterp.get_band_structure()
        self.assertIsNotNone(sbs)
        self.assertTupleEqual(sbs.bands[Spin.up].shape, (6, 137))
        kpaths = [['L', 'K']]
        kp_lbl = {
            'L': np.array([0.5, 0.5, 0.5]),
            'K': np.array([0.375, 0.375, 0.75])
        }
        sbs = self.bztInterp.get_band_structure(kpaths, kp_lbl)
        self.assertIsNotNone(sbs)
        self.assertTupleEqual(sbs.bands[Spin.up].shape, (6, 20))

        sbs = self.bztInterp_sp.get_band_structure()
        self.assertIsNotNone(sbs)
        self.assertTupleEqual(sbs.bands[Spin.up].shape, (6, 143))
        self.assertTupleEqual(sbs.bands[Spin.down].shape, (4, 143))

    def test_tot_dos(self):
        tot_dos = self.bztInterp.get_dos(T=200, npts_mu=100)
        self.assertIsNotNone(tot_dos)
        self.assertEqual(len(tot_dos.energies), 100)
        self.assertAlmostEqual(tot_dos.densities[Spin.up][0], 1.35371715, 5)

        tot_dos = self.bztInterp_sp.get_dos(T=200, npts_mu=100)
        self.assertIsNotNone(tot_dos)
        self.assertEqual(len(tot_dos.energies), 100)
        self.assertAlmostEqual(tot_dos.densities[Spin.up][75], 88.034456, 5)
        self.assertAlmostEqual(tot_dos.densities[Spin.down][75], 41.421367, 5)

    def test_tot_proj_dos(self):
        tot_proj_dos = self.bztInterp.get_dos(partial_dos=True,
                                              T=200,
                                              npts_mu=100)
        self.assertIsNotNone(tot_proj_dos)
        self.assertEqual(len(tot_proj_dos.get_spd_dos().values()), 3)
        pdos = tot_proj_dos.get_spd_dos()[OrbitalType.s].densities[Spin.up][75]
        self.assertAlmostEqual(pdos, 2490.169396, 5)

        tot_proj_dos = self.bztInterp_sp.get_dos(partial_dos=True,
                                                 T=200,
                                                 npts_mu=100)
        self.assertIsNotNone(tot_proj_dos)
        self.assertEqual(len(tot_proj_dos.get_spd_dos().values()), 3)
        pdos = tot_proj_dos.get_spd_dos()[OrbitalType.s].densities[Spin.up][75]
        self.assertAlmostEqual(pdos, 166.4933305, 5)
        pdos = tot_proj_dos.get_spd_dos()[OrbitalType.s].densities[
            Spin.down][75]
        self.assertAlmostEqual(pdos, 272.194174, 5)


@unittest.skipIf(not BOLTZTRAP2_PRESENT, "No boltztrap2, skipping tests...")
class BztTransportPropertiesTest(unittest.TestCase):
    def setUp(self):
        loader = VasprunBSLoader(vrun)
        bztInterp = BztInterpolator(loader, lpfac=2)
        self.bztTransp = BztTransportProperties(bztInterp,
                                                temp_r=np.arange(300, 600, 100))
        self.assertIsNotNone(self.bztTransp)
        warnings.simplefilter("ignore")

        self.bztTransp = BztTransportProperties(bztInterp,
                                                doping=10.**np.arange(20, 22),
                                                temp_r=np.arange(300, 600, 100))
        self.assertIsNotNone(self.bztTransp)
        self.assertEqual(self.bztTransp.contain_props_doping, True)

        warnings.simplefilter("ignore")

        bztInterp = BztInterpolator(loader, lpfac=2)
        self.bztTransp = BztTransportProperties(bztInterp,
                                                temp_r=np.arange(
                                                    300, 600, 100),
                                                save_bztTranspProps=True,
                                                fname=bzttransp_fn)
        self.assertIsNotNone(self.bztTransp)
        warnings.simplefilter("ignore")

        bztInterp = BztInterpolator(loader, lpfac=2)
        self.bztTransp = BztTransportProperties(bztInterp,
                                                load_bztTranspProps=True,
                                                fname=bzttransp_fn)
        self.assertIsNotNone(self.bztTransp)
        warnings.simplefilter("ignore")

        loader_sp = VasprunBSLoader(vrun_sp)
        bztInterp_sp = BztInterpolator(loader_sp, lpfac=2)
        self.bztTransp_sp = BztTransportProperties(bztInterp_sp,
                                                   temp_r=np.arange(
                                                       300, 600, 100))
        self.assertIsNotNone(self.bztTransp_sp)
        warnings.simplefilter("ignore")

        bztInterp_sp = BztInterpolator(loader_sp, lpfac=2)
        self.bztTransp_sp = BztTransportProperties(bztInterp_sp,
                                                   temp_r=np.arange(
                                                       300, 600, 100),
                                                   save_bztTranspProps=True,
                                                   fname=bzttransp_fn)
        self.assertIsNotNone(self.bztTransp_sp)
        warnings.simplefilter("ignore")

        bztInterp_sp = BztInterpolator(loader_sp, lpfac=2)
        self.bztTransp_sp = BztTransportProperties(bztInterp_sp,
                                                   load_bztTranspProps=True,
                                                   fname=bzttransp_fn)
        self.assertIsNotNone(self.bztTransp_sp)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_properties(self):
        for p in [
                self.bztTransp.Conductivity_mu, self.bztTransp.Seebeck_mu,
                self.bztTransp.Kappa_mu, self.bztTransp.Effective_mass_mu,
                self.bztTransp.Power_Factor_mu
        ]:
            self.assertTupleEqual(p.shape, (3, 3686, 3, 3))

        for p in [
                self.bztTransp.Carrier_conc_mu,
                self.bztTransp.Hall_carrier_conc_trace_mu
        ]:
            self.assertTupleEqual(p.shape, (3, 3686))

        for p in [
                self.bztTransp_sp.Conductivity_mu,
                self.bztTransp_sp.Seebeck_mu, self.bztTransp_sp.Kappa_mu,
                self.bztTransp_sp.Effective_mass_mu,
                self.bztTransp_sp.Power_Factor_mu
        ]:
            self.assertTupleEqual(p.shape, (3, 3252, 3, 3))

        for p in [
                self.bztTransp_sp.Carrier_conc_mu,
                self.bztTransp_sp.Hall_carrier_conc_trace_mu
        ]:
            self.assertTupleEqual(p.shape, (3, 3252))

    def test_compute_properties_doping(self):
        self.bztTransp.compute_properties_doping(doping=10.**np.arange(20, 22))
        for p in [
                self.bztTransp.Conductivity_doping,
                self.bztTransp.Seebeck_doping, self.bztTransp.Kappa_doping,
                self.bztTransp.Effective_mass_doping,
                self.bztTransp.Power_Factor_doping
        ]:
            self.assertTupleEqual(p['n'].shape, (3, 2, 3, 3))
            self.assertEqual(self.bztTransp.contain_props_doping, True)

        self.bztTransp_sp.compute_properties_doping(
            doping=10.**np.arange(20, 22))
        for p in [
                self.bztTransp_sp.Conductivity_doping,
                self.bztTransp_sp.Seebeck_doping,
                self.bztTransp_sp.Kappa_doping,
                self.bztTransp_sp.Effective_mass_doping,
                self.bztTransp_sp.Power_Factor_doping
        ]:
            self.assertTupleEqual(p['n'].shape, (3, 2, 3, 3))
            self.assertEqual(self.bztTransp_sp.contain_props_doping, True)


@unittest.skipIf(not BOLTZTRAP2_PRESENT, "No boltztrap2, skipping tests...")
class BztPlotterTest(unittest.TestCase):
    def test_plot(self):
        loader = VasprunBSLoader(vrun)
        bztInterp = BztInterpolator(loader, lpfac=2)
        bztTransp = BztTransportProperties(bztInterp,
                                           temp_r=np.arange(300, 600, 100))
        self.bztPlotter = BztPlotter(bztTransp, bztInterp)
        self.assertIsNotNone(self.bztPlotter)
        fig = self.bztPlotter.plot_props('S', 'mu', 'temp', temps=[300, 500])
        self.assertIsNotNone(fig)
        fig = self.bztPlotter.plot_bands()
        self.assertIsNotNone(fig)
        fig = self.bztPlotter.plot_dos()
        self.assertIsNotNone(fig)


if __name__ == '__main__':
    unittest.main()
