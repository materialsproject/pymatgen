# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import os
import unittest

from monty.serialization import loadfn

from pymatgen.analysis.defects.thermodynamics import DefectPhaseDiagram
from pymatgen.util.testing import PymatgenTest


class DefectsThermodynamicsTest(PymatgenTest):
    def setUp(self):

        self.entries = list(loadfn(os.path.join(os.path.dirname(__file__), "GaAs_test_defentries.json")).values())
        for entry in self.entries:
            entry.parameters.update( {'vbm': 2.6682})

    def test_good_test_data(self):
        self.assertEqual(len(self.entries), 48)

    def test_suggest_charges(self):
        pd = DefectPhaseDiagram(self.entries, 2.6682, 2.0)
        suggested_charges = pd.suggest_charges()
        for k in [
                "Int_As_mult1", "Int_Ga_mult1", "Sub_As_on_Ga_mult4", "Sub_Ga_on_As_mult4", "Vac_As_mult4",
                "Vac_Ga_mult4"
        ]:
            self.assertTrue(
                len(suggested_charges[k]) > 0, "Could not find any suggested charges for {} with band_gap of {}".format(
                    k, pd.band_gap))

        pd = DefectPhaseDiagram(self.entries, 2.6682, 1.0)
        suggested_charges = pd.suggest_charges()
        for k in ["Vac_As_mult4", "Vac_Ga_mult4"]:
            self.assertTrue(
                len(suggested_charges[k]) > 0, "Could not find any suggested charges for {} with band_gap of {}".format(
                    k, pd.band_gap))

    def test_entries(self):
        pd = DefectPhaseDiagram(self.entries, 2.6682, 2.0)

        all_stable_entries = pd.all_stable_entries

        self.assertEqual(len(pd.defect_types), 6)
        self.assertEqual(len(all_stable_entries), sum([len(v) for v in pd.stable_charges.values()]))


if __name__ == "__main__":
    unittest.main()
