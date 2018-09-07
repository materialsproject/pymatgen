# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest
import os
import warnings

from pymatgen.io.vasp import Vasprun

from pymatgen.electronic_structure.boltztrap2 import BandstructureLoader,\
    VasprunLoader,BztInterpolator,BztTransportProperties,BztPlotter

import numpy as np
from monty.serialization import loadfn



test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files/boltztrap2/')

vrunfile=os.path.join(test_dir,'vasprun-PbTe_uniform_bs.xml')
vrun = Vasprun(vrunfile)

class BandstructureLoaderTest(unittest.TestCase):

    def setUp(self):
        bs = loadfn(os.path.join(test_dir, "PbTe_bandstructure.json"))
        
        self.loader = BandstructureLoader(bs,vrun.structures[-1])
        self.assertIsNotNone(self.loader)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.resetwarnings()

    def test_properties(self):
        self.assertTupleEqual(self.loader.ebands.shape,(20,120))
        self.assertAlmostEqual(self.loader.fermi,0.185266535678,5)
        self.assertAlmostEqual(self.loader.structure.lattice.a,4.64303565932548,5)
                               
        
    def test_get_volume(self):
        self.assertAlmostEqual(self.loader.get_volume(),477.6256714925874,5)


class VasprunLoaderTest(unittest.TestCase):

    def setUp(self):
        self.loader = VasprunLoader(vrun)
        self.assertIsNotNone(self.loader)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.resetwarnings()

    def test_properties(self):
        self.assertTupleEqual(self.loader.ebands.shape,(20,120))
        self.assertAlmostEqual(self.loader.fermi,0.185266535678,5)
        self.assertAlmostEqual(self.loader.structure.lattice.a,4.64303565932548,5)
                               
        
    def test_get_volume(self):
        self.assertAlmostEqual(self.loader.get_volume(),477.6256714925874,5)

    def test_from_file(self):
        self.loader = VasprunLoader().from_file(vrunfile)
        self.assertIsNotNone(self.loader)



class BztInterpolatorTest(unittest.TestCase):

    def setUp(self):
        loader = VasprunLoader().from_file(vrunfile)
        self.bztInterp = BztInterpolator(loader)
        self.assertIsNotNone(self.bztInterp)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.resetwarnings()

    def test_properties(self):
        self.assertTupleEqual(self.bztInterp.cband.shape,(5, 3, 3, 3, 148877))
        self.assertTupleEqual(self.bztInterp.eband.shape,(5, 148877))
        self.assertTupleEqual(self.bztInterp.coeffs.shape,(5, 1429))
        self.assertEqual(self.bztInterp.nemax,12)
        
    def test_get_band_structure(self):
        sbs = self.bztInterp.get_band_structure()
        self.assertIsNotNone(sbs)
        self.assertTupleEqual(list(sbs.bands.values())[0].shape,(5,137))

    def test_tot_dos(self):
        tot_dos = self.bztInterp.get_dos(T=200)
        self.assertIsNotNone(tot_dos)
        self.assertEqual(len(tot_dos.energies),10000)

    def test_tot_proj_dos(self):
        tot_proj_dos = self.bztInterp.get_dos(partial_dos=True,T=200)
        self.assertIsNotNone(tot_proj_dos)
        self.assertEqual(len(tot_proj_dos.get_spd_dos().values()),3)


class BztTransportPropertiesTest(unittest.TestCase):

    def setUp(self):
        loader = VasprunLoader().from_file(vrunfile)
        bztInterp = BztInterpolator(loader)
        self.bztTransp = BztTransportProperties(bztInterp,temp_r = np.arange(300,600,100))
        self.assertIsNotNone(self.bztTransp)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.resetwarnings()

    def test_properties(self):
        for p in [self.bztTransp.Conductivity_mu,self. bztTransp.Seebeck_mu,
                  self.bztTransp.Kappa_mu,self.bztTransp.Effective_mass_mu,
                  self.bztTransp.Power_Factor_mu]:
            
            self.assertTupleEqual(p.shape,(3, 3670, 3, 3))
        
        for p in [self.bztTransp.Carrier_conc_mu,self.bztTransp.Hall_carrier_conc_trace_mu]:
            self.assertTupleEqual(p.shape,(3, 3670))

        
    def test_compute_properties_doping(self):
        self.bztTransp.compute_properties_doping(doping=10.**np.arange(16,23))
        for p in [self.bztTransp.Conductivity_doping,self. bztTransp.Seebeck_doping,
                  self.bztTransp.Kappa_doping,self.bztTransp.Effective_mass_doping,
                  self.bztTransp.Power_Factor_doping]:
            
            self.assertTupleEqual(p['n'].shape,(3, 7, 3, 3))
        

class BztPlotterTest(unittest.TestCase):
    
    def setUp(self):
        loader = VasprunLoader().from_file(vrunfile)
        bztInterp = BztInterpolator(loader)
        bztTransp = BztTransportProperties(bztInterp,temp_r = np.arange(300,600,100))
        self.bztPlotter = BztPlotter(bztTransp,bztInterp)
        self.assertIsNotNone(self.bztPlotter)
        
    def tearDown(self):
        warnings.resetwarnings()

    def test_plot_props(self):
        fig = self.bztPlotter.plot_props('S','mu','temp',temps=[300,500])
        self.assertIsNotNone(fig)
        
    def test_plot_bands(self):
        fig = self.bztPlotter.plot_bands()
        self.assertIsNotNone(fig)
    
    def test_plot_dos(self):
        fig = self.bztPlotter.plot_dos()
        self.assertIsNotNone(fig)
        
        
        
if __name__ == '__main__':
    unittest.main()
