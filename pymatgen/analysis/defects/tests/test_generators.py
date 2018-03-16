# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.sites import Site
from pymatgen.analysis.defects.generators import VacancyGenerator, InterstitialGenerator


class VacancyGeneratorTest(PymatgenTest):
    def test_vacancy_gen(self):
        struc = PymatgenTest.get_structure("VO2")
        vac_gen = VacancyGenerator(struc)

        vacs = list(vac_gen)
        self.assertEqual(len(vacs), 2)

        multiplicities = {str(v.defect_site.specie): v.multiplicity for v in vacs}
        self.assertEqual(multiplicities, {"O": 4, "V": 2})

    def test_vacancy_gen_charges(self):
        # Ensure correct BV charges are assigned
        struc = PymatgenTest.get_structure("VO2")
        vac_gen = VacancyGenerator(struc, include_bv_charge=True)
        for vac in vac_gen:
            if str(vac.defect_site.specie) == "V":
                self.assertEqual(vac.charge, -4)
            if str(vac.defect_site.specie) == "O":
                self.assertEqual(vac.charge, 2)


class InterstitialGeneratorTest(PymatgenTest):
    def test_int_gen(self):
        struc = PymatgenTest.get_structure("VO2")
        int_gen = InterstitialGenerator(struc, "Li")

        ints = list(int_gen)
        self.assertEqual(len(ints), 2)

        multiplicities = [i.multiplicity for i in ints]
        self.assertEqual(multiplicities, [1, 1])

        self.assertEqual(str(ints[0].specie), "Li")
        self.assertEqual(str(ints[1].specie), "Li")

        self.assertArrayAlmostEqual(ints[0].defect_site.coords, (0.9106, 0.3078, 0.3078))
        self.assertArrayAlmostEqual(ints[1].defect_site.coords, (1.5177, 0.7183, 0.7183))


if __name__ == "__main__":
    unittest.main()
