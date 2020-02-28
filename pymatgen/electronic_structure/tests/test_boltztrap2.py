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
    from pymatgen.electronic_structure.boltztrap2 import BandstructureLoader, \
        VasprunLoader, BztInterpolator, BztTransportProperties, BztPlotter, \
        merge_up_down_doses

    BOLTZTRAP2_PRESENT = True
except Exception:
    BOLTZTRAP2_PRESENT = False

BOLTZTRAP2_PRESENT = False

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files/boltztrap2/')

vrunfile = os.path.join(test_dir, 'vasprun.xml')
vrun = Vasprun(vrunfile, parse_projected_eigen=True)

vrunfile_sp = os.path.join(test_dir, 'vasprun_spin.xml')
vrun_sp = Vasprun(vrunfile_sp, parse_projected_eigen=True)


@unittest.skipIf(not BOLTZTRAP2_PRESENT, "No boltztrap2, skipping tests...")
class BandstructureLoaderTest(unittest.TestCase):

    def setUp(self):
        bs = loadfn(os.path.join(test_dir, "PbTe_bandstructure.json"))
        bs_sp = loadfn(os.path.join(test_dir, "N2_bandstructure.json"))
        self.loader = BandstructureLoader(bs, vrun.structures[-1])
        self.assertIsNotNone(self.loader)

        self.loader_sp_up = BandstructureLoader(bs_sp, vrun_sp.structures[-1], spin=1)
        self.loader_sp_dn = BandstructureLoader(bs_sp, vrun_sp.structures[-1], spin=-1)
        self.assertTupleEqual(self.loader_sp_up.ebands.shape, (12, 198))
        self.assertTupleEqual(self.loader_sp_dn.ebands.shape, (12, 198))
        self.assertIsNotNone(self.loader_sp_dn)
        self.assertIsNotNone(self.loader_sp_up)

        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_properties(self):
        self.assertTupleEqual(self.loader.ebands.shape, (20, 120))
        self.assertAlmostEqual(self.loader.fermi, 0.185266535678, 5)
        self.assertAlmostEqual(self.loader.structure.lattice.a, 4.64303565932548, 5)

    def test_get_volume(self):
        self.assertAlmostEqual(self.loader.get_volume(), 477.6256714925874, 5)

    def test_set_upper_lower_bands(self):
        min_bnd = min(self.loader_sp_up.ebands.min(), self.loader_sp_dn.ebands.min())
        max_bnd = max(self.loader_sp_up.ebands.max(), self.loader_sp_dn.ebands.max())
        self.loader_sp_up.set_upper_lower_bands(min_bnd, max_bnd)
        self.loader_sp_dn.set_upper_lower_bands(min_bnd, max_bnd)
        self.assertTupleEqual(self.loader_sp_up.ebands.shape, (14, 198))
        self.assertTupleEqual(self.loader_sp_dn.ebands.shape, (14, 198))


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
        self.assertAlmostEqual(self.loader.structure.lattice.a, 4.64303565932548, 5)

    def test_get_volume(self):
        self.assertAlmostEqual(self.loader.get_volume(), 477.6256714925874, 5)

    def test_from_file(self):
        self.loader = VasprunLoader().from_file(vrunfile)
        self.assertIsNotNone(self.loader)


@unittest.skipIf(not BOLTZTRAP2_PRESENT, "No boltztrap2, skipping tests...")
class BztInterpolatorTest(unittest.TestCase):

    def setUp(self):
        self.loader = VasprunLoader(vrun)
        self.assertTupleEqual(self.loader.proj.shape, (120, 20, 2, 9))
        self.bztInterp = BztInterpolator(self.loader, lpfac=2)
        self.assertIsNotNone(self.bztInterp)
        warnings.simplefilter("ignore")

        bs_sp = loadfn(os.path.join(test_dir, "N2_bandstructure.json"))
        loader_sp_up = BandstructureLoader(bs_sp, vrun_sp.structures[-1], spin=1)
        loader_sp_dn = BandstructureLoader(bs_sp, vrun_sp.structures[-1], spin=-1)

        min_bnd = min(loader_sp_up.ebands.min(), loader_sp_dn.ebands.min())
        max_bnd = max(loader_sp_up.ebands.max(), loader_sp_dn.ebands.max())
        loader_sp_up.set_upper_lower_bands(min_bnd, max_bnd)
        loader_sp_dn.set_upper_lower_bands(min_bnd, max_bnd)

        self.bztI_up = BztInterpolator(loader_sp_up, lpfac=2, energy_range=np.inf, curvature=False)
        self.bztI_dn = BztInterpolator(loader_sp_dn, lpfac=2, energy_range=np.inf, curvature=False)

    def tearDown(self):
        warnings.simplefilter("default")

    def test_properties(self):
        self.assertTupleEqual(self.bztInterp.cband.shape, (5, 3, 3, 3, 29791))
        self.assertTupleEqual(self.bztInterp.eband.shape, (5, 29791))
        self.assertTupleEqual(self.bztInterp.coeffs.shape, (5, 322))
        self.assertEqual(self.bztInterp.nemax, 12)

    def test_get_band_structure(self):
        sbs = self.bztInterp.get_band_structure()
        self.assertIsNotNone(sbs)
        self.assertTupleEqual(sbs.bands[Spin.up].shape, (5, 137))

    def test_tot_dos(self):
        tot_dos = self.bztInterp.get_dos(T=200, npts_mu=100)
        self.assertIsNotNone(tot_dos)
        self.assertEqual(len(tot_dos.energies), 100)
        self.assertAlmostEqual(tot_dos.densities[Spin.up][0], 1.42859939, 5)

        dos_up = self.bztI_up.get_dos(partial_dos=False, npts_mu=100)
        dos_dn = self.bztI_dn.get_dos(partial_dos=False, npts_mu=100)
        cdos = merge_up_down_doses(dos_up, dos_dn)
        self.assertAlmostEqual(cdos.densities[Spin.down][50], 92.87836778, 5)
        self.assertAlmostEqual(cdos.densities[Spin.up][45], 9.564067, 5)
        self.assertEqual(len(cdos.energies), 100)

    def test_tot_proj_dos(self):
        tot_proj_dos = self.bztInterp.get_dos(partial_dos=True, T=200, npts_mu=100)
        self.assertIsNotNone(tot_proj_dos)
        self.assertEqual(len(tot_proj_dos.get_spd_dos().values()), 3)
        pdos = tot_proj_dos.get_spd_dos()[OrbitalType.s].densities[Spin.up][0]
        self.assertAlmostEqual(pdos, 15.474392020, 5)


@unittest.skipIf(not BOLTZTRAP2_PRESENT, "No boltztrap2, skipping tests...")
class BztTransportPropertiesTest(unittest.TestCase):

    def setUp(self):
        loader = VasprunLoader(vrun)
        bztInterp = BztInterpolator(loader, lpfac=2)
        self.bztTransp = BztTransportProperties(bztInterp, temp_r=np.arange(300, 600, 100))
        self.assertIsNotNone(self.bztTransp)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_properties(self):
        for p in [self.bztTransp.Conductivity_mu, self.bztTransp.Seebeck_mu,
                  self.bztTransp.Kappa_mu, self.bztTransp.Effective_mass_mu,
                  self.bztTransp.Power_Factor_mu]:
            self.assertTupleEqual(p.shape, (3, 3670, 3, 3))

        for p in [self.bztTransp.Carrier_conc_mu, self.bztTransp.Hall_carrier_conc_trace_mu]:
            self.assertTupleEqual(p.shape, (3, 3670))

    def test_compute_properties_doping(self):
        self.bztTransp.compute_properties_doping(doping=10. ** np.arange(20, 22))
        for p in [self.bztTransp.Conductivity_doping, self.bztTransp.Seebeck_doping,
                  self.bztTransp.Kappa_doping, self.bztTransp.Effective_mass_doping,
                  self.bztTransp.Power_Factor_doping]:
            self.assertTupleEqual(p['n'].shape, (3, 2, 3, 3))


@unittest.skipIf(not BOLTZTRAP2_PRESENT, "No boltztrap2, skipping tests...")
class BztPlotterTest(unittest.TestCase):

    def test_plot(self):
        loader = VasprunLoader(vrun)
        bztInterp = BztInterpolator(loader, lpfac=2)
        bztTransp = BztTransportProperties(bztInterp, temp_r=np.arange(300, 600, 100))
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
