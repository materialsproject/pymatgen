#!/usr/bin/env python

import unittest
import os
import itertools

from pymatgen.analysis.symmetry_fitter import SymmetryFitter
from pymatgen import __file__
from pymatgen.transformations.site_transformations import \
    RemoveSitesTransformation
from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen.symmetry.finder import SymmetryFinder

test_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                        'test_files')


class SymmetryFitterTest(unittest.TestCase):

    def test_init(self):
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR.LiFePO4'),
                             check_for_POTCAR=False)
        struct = p.structure
        structs = []
        for i, j in itertools.combinations(xrange(4, 8), 2):
            trans = RemoveSitesTransformation([i, j])
            structs.append(trans.apply_transformation(struct))
        sg = SymmetryFinder(struct, 0.1).get_spacegroup()
        fitter = SymmetryFitter(structs, sg, 0.1)

        self.assertEqual(len(fitter.unique_groups), 3)

        structs = []
        for i in xrange(4, 8):
            trans = RemoveSitesTransformation([i])
            structs.append(trans.apply_transformation(struct))
        fitter = SymmetryFitter(structs, sg, 0.1)
        self.assertEqual(len(fitter.unique_groups), 1)

if __name__ == '__main__':
    unittest.main()
