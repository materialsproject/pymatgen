# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest

from pymatgen.core.sites import Site
from pymatgen.analysis.defects.core import Vacancy, Interstitial, Substitution
from pymatgen.util.testing import PymatgenTest


class DefectsCoreTest(PymatgenTest):
    def test_vacancy(self):
        struc = PymatgenTest.get_structure("VO2")
        V_index = struc.indices_from_symbol("V")[0]
        vac = Vacancy(struc, struc[V_index])

        # test generation and super cell
        vac_struc = vac.generate_defect_structure(1)
        self.assertEqual(vac_struc.composition.as_dict(), {"V": 1, "O": 4})

        vac_struc = vac.generate_defect_structure(2)
        self.assertEqual(vac_struc.composition.as_dict(), {"V": 15, "O": 32})

        vac_struc = vac.generate_defect_structure(3)
        self.assertEqual(vac_struc.composition.as_dict(), {"V": 53, "O": 108})

        # test charge
        vac = Vacancy(struc, struc[V_index])
        vac_struc = vac.generate_defect_structure(1)
        self.assertEqual(vac_struc.charge, 0.0)

        vac = Vacancy(struc, struc[V_index], charge=1.0)
        vac_struc = vac.generate_defect_structure(1)
        self.assertEqual(vac_struc.charge, 1.0)

        vac = Vacancy(struc, struc[V_index], charge=-1.0)
        vac_struc = vac.generate_defect_structure(1)
        self.assertEqual(vac_struc.charge, -1.0)

        # test multiplicity
        vac = Vacancy(struc, struc[V_index])
        self.assertEqual(vac.multiplicity, 1.0)

        vac = Vacancy(struc, struc[V_index], multiplicity=4.0)
        self.assertEqual(vac.multiplicity, 4.0)

    def test_interstitial(self):
        struc = PymatgenTest.get_structure("VO2")
        V_index = struc.indices_from_symbol("V")[0]

        int_site = Site("V", struc[V_index].coords)
        interstitial = Interstitial(struc, int_site)

        # test generation and super cell
        int_struc = interstitial.generate_defect_structure(1)
        self.assertEqual(int_struc.composition.as_dict(), {"V": 3, "O": 4})

        int_struc = interstitial.generate_defect_structure(2)
        self.assertEqual(int_struc.composition.as_dict(), {"V": 17, "O": 32})

        int_struc = interstitial.generate_defect_structure(3)
        self.assertEqual(int_struc.composition.as_dict(), {"V": 55, "O": 108})

        # test charge
        interstitial = Interstitial(struc, int_site)
        int_struc = interstitial.generate_defect_structure(1)
        self.assertEqual(int_struc.charge, 0.0)

        interstitial = Interstitial(struc, int_site, charge=1.0)
        int_struc = interstitial.generate_defect_structure(1)
        self.assertEqual(int_struc.charge, 1.0)

        interstitial = Interstitial(struc, int_site, charge=-1.0)
        int_struc = interstitial.generate_defect_structure(1)
        self.assertEqual(int_struc.charge, -1.0)

        # test multiplicity
        interstitial = Interstitial(struc, int_site)
        self.assertEqual(interstitial.multiplicity, 1.0)

        interstitial = Interstitial(struc, int_site, multiplicity=4.0)
        self.assertEqual(interstitial.multiplicity, 4.0)


if __name__ == "__main__":
    unittest.main()
