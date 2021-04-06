# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest

from pymatgen.analysis.defects.core import Vacancy
from pymatgen.transformations.defect_transformations import DefectTransformation
from pymatgen.util.testing import PymatgenTest


class DefectTransformationTest(PymatgenTest):
    def test_apply_transformation(self):
        struc = PymatgenTest.get_structure("VO2")
        vac = Vacancy(struc, struc[0], charge=1)

        def_transform = DefectTransformation([2, 2, 2], vac)
        trans_structure = def_transform.apply_transformation(struc)
        self.assertEqual(len(trans_structure), 47)

        # confirm that transformation doesnt work for bulk structures
        # which are slightly different than those used for defect object
        # scaled volume
        scaled_struc = struc.copy()
        scaled_struc.scale_lattice(1.1 * struc.volume)
        self.assertRaises(ValueError, def_transform.apply_transformation, scaled_struc)
        # slightly different atomic positions
        pert_struc = struc.copy()
        pert_struc.perturb(0.1)
        self.assertRaises(ValueError, def_transform.apply_transformation, pert_struc)


if __name__ == "__main__":
    unittest.main()
