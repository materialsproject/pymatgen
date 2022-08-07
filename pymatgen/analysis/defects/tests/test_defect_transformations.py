# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest

from pymatgen.analysis.defects.core import Vacancy
from pymatgen.analysis.defects.defect_transformations import DefectTransformation
from pymatgen.util.testing import PymatgenTest


class DefectTransformationTest(PymatgenTest):
    def test_apply_transformation(self):
        struct = PymatgenTest.get_structure("VO2")
        vac = Vacancy(struct, struct[0], charge=1)

        def_transform = DefectTransformation([2, 2, 2], vac)
        trans_structure = def_transform.apply_transformation(struct)
        self.assertEqual(len(trans_structure), 47)

        # confirm that transformation doesn't work for bulk structures
        # which are slightly different than those used for defect object
        # scaled volume
        scaled_struct = struct.copy()
        scaled_struct.scale_lattice(1.1 * struct.volume)
        self.assertRaises(ValueError, def_transform.apply_transformation, scaled_struct)
        # slightly different atomic positions
        pert_struc = struct.copy()
        pert_struc.perturb(0.1)
        self.assertRaises(ValueError, def_transform.apply_transformation, pert_struc)


if __name__ == "__main__":
    unittest.main()
