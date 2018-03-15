# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.sites import Site
from pymatgen.analysis.defects.generators import VacancyGenerator
from pymatgen.analysis.defects.core import Vacancy


class VacancyGeneratorTest(PymatgenTest):
    def test_vacancy_gen(self):
        struc = PymatgenTest.get_structure("VO2")
        vac_gen = VacancyGenerator(struc)

        vacs = list(vac_gen)
        self.assertEqual(len(vacs), 2)

        multiplicities = {str(v.defect_site.specie): v.multiplicity for v in vacs}
        self.assertEqual(multiplicities, {"O": 4, "V": 2})

        vac_gen = VacancyGenerator(struc)
        for vac in vac_gen:
            self.assertIsInstance(vac, Vacancy)

        # Ensure correct BV charges are assigned
        vac_gen = VacancyGenerator(struc, include_bv_charge=True)
        for vac in vac_gen:
            if str(vac.defect_site.specie) == "V":
                self.assertEqual(vac.charge, -4)
            if str(vac.defect_site.specie) == "O":
                self.assertEqual(vac.charge, 2)


if __name__ == "__main__":
    unittest.main()
