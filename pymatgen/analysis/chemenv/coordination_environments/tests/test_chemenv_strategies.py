#!/usr/bin/env python


__author__ = 'waroquiers'

import unittest
from pymatgen.util.testing import PymatgenTest

from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import DistanceCutoffFloat
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import AngleCutoffFloat
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import CSMFloat
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import AdditionalConditionInt
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy


class StrategyOptionsTest(PymatgenTest):

    def test_options(self):
        # DistanceCutoffFloat
        with self.assertRaises(ValueError) as cm:
            DistanceCutoffFloat(0.5)
        self.assertEqual(str(cm.exception), 'Distance cutoff should be between 1.0 and +infinity')
        dc1 = DistanceCutoffFloat(1.2)
        dc1_dict = dc1.as_dict()
        dc2 = DistanceCutoffFloat.from_dict(dc1_dict)
        self.assertEqual(dc1, dc2)

        # AngleCutoffFloat
        with self.assertRaises(ValueError) as cm:
            AngleCutoffFloat(1.2)
        self.assertEqual(str(cm.exception), 'Angle cutoff should be between 0.0 and 1.0')
        ac1 = AngleCutoffFloat(0.3)
        ac1_dict = ac1.as_dict()
        ac2 = AngleCutoffFloat.from_dict(ac1_dict)
        self.assertEqual(ac1, ac2)

        # CSMFloat
        with self.assertRaises(ValueError) as cm:
            CSMFloat(100.1)
        self.assertEqual(str(cm.exception), 'Continuous symmetry measure limits should be between 0.0 and 100.0')
        csm1 = CSMFloat(0.458)
        csm1_dict = csm1.as_dict()
        csm2 = CSMFloat.from_dict(csm1_dict)
        self.assertEqual(csm1, csm2)

        # AdditionalConditions
        with self.assertRaises(ValueError) as cm:
            AdditionalConditionInt(5)
        self.assertEqual(str(cm.exception), 'Additional condition 5 is not allowed')
        with self.assertRaises(ValueError) as cm:
            AdditionalConditionInt(0.458)
        self.assertEqual(str(cm.exception), 'Additional condition 0.458 is not an integer')
        acd1 = AdditionalConditionInt(3)
        acd1_dict = acd1.as_dict()
        acd2 = AdditionalConditionInt.from_dict(acd1_dict)
        self.assertEqual(acd1, acd2)

    def test_strategies(self):
        simplest_strategy = SimplestChemenvStrategy()
        self.assertTrue(simplest_strategy.uniquely_determines_coordination_environments)
        self.assertAlmostEqual(simplest_strategy.continuous_symmetry_measure_cutoff, 10.0)
        self.assertAlmostEqual(simplest_strategy.distance_cutoff, 1.4)
        self.assertAlmostEqual(simplest_strategy.angle_cutoff, 0.3)

        simplest_strategy = SimplestChemenvStrategy(distance_cutoff=1.3, angle_cutoff=0.45,
                                                    continuous_symmetry_measure_cutoff=8.26)
        self.assertAlmostEqual(simplest_strategy.continuous_symmetry_measure_cutoff, 8.26)
        self.assertAlmostEqual(simplest_strategy.distance_cutoff, 1.3)
        self.assertAlmostEqual(simplest_strategy.angle_cutoff, 0.45)

        simplest_strategy.set_option('distance_cutoff', 1.5)
        self.assertAlmostEqual(simplest_strategy.distance_cutoff, 1.5)

        with self.assertRaises(ValueError) as cm:
            simplest_strategy.set_option('distance_cutoff', 0.5)
        self.assertEqual(str(cm.exception), 'Distance cutoff should be between 1.0 and +infinity')

        simplest_strategy.set_option('angle_cutoff', 0.2)
        self.assertAlmostEqual(simplest_strategy.angle_cutoff, 0.2)

        with self.assertRaises(ValueError) as cm:
            simplest_strategy.set_option('angle_cutoff', 1.5)
        self.assertEqual(str(cm.exception), 'Angle cutoff should be between 0.0 and 1.0')

        simplest_strategy.setup_options({'distance_cutoff': 1.4,
                                         'additional_condition': 3,
                                         'continuous_symmetry_measure_cutoff': 8.5})
        self.assertAlmostEqual(simplest_strategy.distance_cutoff, 1.4)
        self.assertAlmostEqual(simplest_strategy.continuous_symmetry_measure_cutoff, 8.5)
        self.assertEqual(simplest_strategy.additional_condition, 3)

        with self.assertRaises(ValueError) as cm:
            simplest_strategy.setup_options({'continuous_symmetry_measure_cutoff': -0.1})
        self.assertEqual(str(cm.exception), 'Continuous symmetry measure limits should be between 0.0 and 100.0')

        with self.assertRaises(ValueError) as cm:
            simplest_strategy.setup_options({'continuous_symmetry_measure_cutoff': 100.1})
        self.assertEqual(str(cm.exception), 'Continuous symmetry measure limits should be between 0.0 and 100.0')



if __name__ == "__main__":
    unittest.main()