#!/usr/bin/env python


__author__ = 'waroquiers'

import unittest
import os
import shutil
from pymatgen.util.testing import PymatgenTest

from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import AngleNbSetWeight


class FakeNbSet(object):
    pass

class DummyStructureEnvironments(object):
    pass


class StrategyWeightsTest(PymatgenTest):

    def test_angle_weight(self):
        fake_nb_set = FakeNbSet()
        dummy_se = DummyStructureEnvironments()
        # Angles for 5
        fake_nb_set.angles = [1.8595833644514066, 2.622518848090717, 3.08570351705799,
                              2.2695472184920042, 2.2695338778592387]
        angle_weight = AngleNbSetWeight(aa=1.0)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        self.assertAlmostEqual(aw, 0.9634354419021528, delta=1e-8)

        angle_weight = AngleNbSetWeight(aa=2.0)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        self.assertAlmostEqual(aw, 0.92820785071319645, delta=1e-8)

        angle_weight = AngleNbSetWeight(aa=0.5)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        self.assertAlmostEqual(aw, 0.98154747307613843, delta=1e-8)

        fake_nb_set.angles = [1.8595833644514066, 2.622518848090717, 3.08570351705799,
                              2.2695472184920042, 2.2695338778592387]

        angle_weight = AngleNbSetWeight(aa=1.0)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        self.assertAlmostEqual(aw, 0.9634354419021528, delta=1e-8)

        angle_weight = AngleNbSetWeight(aa=2.0)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        self.assertAlmostEqual(aw, 0.92820785071319645, delta=1e-8)

        angle_weight = AngleNbSetWeight(aa=0.5)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        self.assertAlmostEqual(aw, 0.98154747307613843, delta=1e-8)

if __name__ == "__main__":
    unittest.main()