#!/usr/bin/env python
from __future__ import division, print_function

import unittest
import numpy as np

from pymatgen.io.abinitio.eos import EOS

def have_scipy():
    try:
        import scipy
        return True
    except ImportError:
        return False

class EOSTestCase(unittest.TestCase):

    def setUp(self):
        self.volumes = np.array([13.72, 14.83, 16.0, 17.23, 18.52])
        self.energies = np.array([-56.29, -56.41, -56.46, -56.46, -56.42])

    @unittest.skipUnless(have_scipy(), "test_fit requires scipy")
    def test_fit(self):
        "Test EOS fit"
        for eos_name in EOS.MODELS:
            eos = EOS(eos_name=eos_name)
            fit = eos.fit(self.volumes, self.energies)
            print(fit)

if __name__ == "__main__":
    unittest.main()
