#!/usr/bin/env python
from __future__ import division, print_function

import unittest
import numpy as np

from pymatgen.io.abinitio.abiobjects import *
from pymatgen.core.physical_constants import Ha_eV

##########################################################################################

class SpinModeTest(unittest.TestCase):

    def test_base(self):
        polarized = SpinMode.asspinmode("polarized")
        other_polarized = SpinMode.asspinmode("polarized")
        unpolarized = SpinMode.asspinmode("unpolarized")

        polarized.to_abivars()

        self.assertTrue(polarized is other_polarized)
        self.assertTrue(polarized == other_polarized)
        self.assertTrue(polarized != unpolarized)

##########################################################################################

class SmearingTest(unittest.TestCase):
    def test_base(self):
        fd1ev = Smearing.assmearing("fermi_dirac:1 eV")
        print(fd1ev)
        fd1ev.to_abivars()

        self.assertTrue(fd1ev)

        same_fd = Smearing.assmearing("fermi_dirac:"+ str(1.0/Ha_eV))

        self.assertTrue(same_fd == fd1ev)

        nosmear = Smearing.nosmearing()

        self.assertFalse(nosmear)
        self.assertTrue(nosmear != fd1ev)
        new_fd1ev = Smearing.from_dict(fd1ev.to_dict)
        self.assertTrue(new_fd1ev == fd1ev)

##########################################################################################

#class ElectronsTest(unittest.TestCase):
#    default_electrons = Electrons()

##########################################################################################

#class AbiStructureTest(unittest.TestCase):

##########################################################################################

#class KSamplingTest(unittest.TestCase):

##########################################################################################

#class RelaxationTest(unittest.TestCase):

##########################################################################################

class PPModelTest(unittest.TestCase):

    def test_base(self):
        godby = PPModel.asppmodel("godby:12 eV")
        print(godby)
        print(repr(godby))
        godby.to_abivars()
        self.assertTrue(godby)

        same_godby = PPModel.asppmodel("godby:"+ str(12.0/Ha_eV))
        self.assertTrue(same_godby == godby)

        noppm = PPModel.noppmodel()

        self.assertFalse(noppm)
        self.assertTrue(noppm != godby)
        new_godby = PPModel.from_dict(godby.to_dict)
        self.assertTrue(new_godby == godby)

##########################################################################################

if __name__ == '__main__':
    unittest.main()
