# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest

from pymatgen.analysis.defects.utils import QModel, eV_to_k, generate_reciprocal_vectors_squared, \
    closestsites
from pymatgen.util.testing import PymatgenTest


class DefectsUtilsTest(PymatgenTest):
    def test_qmodel(self):
        qm = QModel()
        modqm = QModel(beta=2., expnorm=0.5, gamma=0.1)

        #test rho_rec
        self.assertEqual(qm.rho_rec(1.), 0.77880078307140488)
        self.assertEqual(modqm.rho_rec(1.), 0.6814583156907158)

        #test rho_rec_limit0
        self.assertEqual(qm.rho_rec_limit0, -0.25)
        self.assertEqual(modqm.rho_rec_limit0, -0.51)

    def test_eV_to_k(self):
        self.assertAlmostEqual(eV_to_k(1.), 0.9681404248678961)

    def test_generate_reciprocal_vectors_squared(self):
        #test cubic case
        a = 6.
        lattvectors = [[a if i==j else 0. for j in range(3)] for i in range(3)]
        brecip = [1.0966227112321507 for i in range(6)]
        self.assertAlmostEqual(
            list(generate_reciprocal_vectors_squared( lattvectors[0], lattvectors[1], lattvectors[2], 1.3)),
            brecip)

        #test orthorhombic case
        lattconsts = [a, a/2., 3.*a]
        lattvectors = [[lattconsts[i] if i==j else 0. for j in range(3)] for i in range(3)]
        brval = 0.4873878716587337
        brecip = [brval, brval/4., brval/4., brval]
        self.assertAlmostEqual(
            list(generate_reciprocal_vectors_squared( lattvectors[0], lattvectors[1], lattvectors[2], 1.)),
            brecip)

        #test triclinic case
        lattvectors = [[ 1.5, 0.2, 0.3], [0.3, 1.2, .2], [0.5, 0.4, 1.3]]
        brval = 24.28330561545568
        brecip = [brval, brval]
        self.assertAlmostEqual(
            list(generate_reciprocal_vectors_squared( lattvectors[0], lattvectors[1], lattvectors[2], 30.)),
            brecip)

    def test_closest_sites(self):
        struct = PymatgenTest.get_structure("VO2")

        #test O vacancy
        dstruct = struct.copy()
        dstruct.remove_sites([0])
        pos = struct.sites[0].coords
        bsite, dsite = closestsites(struct, dstruct, pos)
        self.assertEqual( bsite[2], 0) #test against index
        self.assertEqual( dsite[2], 4)

        #test V vacancy
        dstruct = struct.copy()
        dstruct.remove_sites([4])
        pos = struct.sites[4].coords
        bsite, dsite = closestsites(struct, dstruct, pos)
        self.assertEqual( bsite[2], 4) #test against index
        self.assertEqual( dsite[2], 1)



if __name__ == "__main__":
    unittest.main()
