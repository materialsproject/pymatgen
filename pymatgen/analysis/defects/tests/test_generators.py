# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest

from pymatgen.util.testing import PymatgenTest
from pymatgen.analysis.defects.generators import VacancyGenerator, InterstitialGenerator, VoronoiInterstitialGenerator


class VacancyGeneratorTest(PymatgenTest):
    def test_vacancy_gen(self):
        struc = PymatgenTest.get_structure("VO2")
        vac_gen = VacancyGenerator(struc)

        vacs = list(vac_gen)
        self.assertEqual(len(vacs), 2)

        multiplicities = {str(v.site.specie): v.multiplicity for v in vacs}
        self.assertEqual(multiplicities, {"O": 4, "V": 2})

    def test_vacancy_gen_charges(self):
        # Ensure correct BV charges are assigned
        struc = PymatgenTest.get_structure("VO2")
        vac_gen = VacancyGenerator(struc, include_bv_charge=True)
        for vac in vac_gen:
            if str(vac.site.specie) == "V":
                self.assertEqual(vac.charge, -4)
            if str(vac.site.specie) == "O":
                self.assertEqual(vac.charge, 2)


class InterstitialGeneratorTest(PymatgenTest):
    def test_int_gen(self):
        struc = PymatgenTest.get_structure("VO2")
        int_gen = InterstitialGenerator(struc, "Li")

        ints = list(int_gen)
        self.assertEqual(len(ints), 2)

        multiplicities = [i.multiplicity for i in ints]
        self.assertEqual(multiplicities, [4, 2])

        self.assertEqual(str(ints[0].site.specie), "Li")
        self.assertEqual(str(ints[1].site.specie), "Li")

        self.assertArrayAlmostEqual(ints[0].site.coords, (0.9106, 0.3078, 0.3078), decimal=4)
        self.assertArrayAlmostEqual(ints[1].site.coords, (1.5177, 0.7183, 0.7183), decimal=4)


class VoronoiInterstitialGeneratorTest(PymatgenTest):
    def test_int_gen(self):
        struc = PymatgenTest.get_structure("VO2")
        int_gen = VoronoiInterstitialGenerator(struc, "Li")

        ints = list(int_gen)
        self.assertEqual(len(ints), 4)

        multiplicities = [i.multiplicity for i in ints]
        self.assertEqual(multiplicities, [8, 8, 4, 4])

        self.assertEqual(str(ints[0].site.specie), "Li")
        self.assertEqual(str(ints[1].site.specie), "Li")
        self.assertEqual(str(ints[2].site.specie), "Li")
        self.assertEqual(str(ints[3].site.specie), "Li")

        self.assertArrayAlmostEqual(ints[0].site.coords, (1.5177146, 2.6784354, 3.9481299))
        self.assertArrayAlmostEqual(ints[1].site.coords, (1.7357692, 3.8392513, 3.8392513))
        self.assertArrayAlmostEqual(ints[2].site.coords, (1.5177146, 3.7168193, 3.7168193))
        self.assertArrayAlmostEqual(ints[3].site.coords, (2.2765713, 2.2575138, 4.5150233))


if __name__ == "__main__":
    unittest.main()
