"""
Unit tests for defect transformations
"""

from __future__ import division

__author__ = "Bharat Medasani"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "0.1"
__maintainier__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__date__ = "Jul 1 2014"

import unittest

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.transformations.defect_transformations import \
    VacancyTransformation, SubstitutionDefectTransformation, \
    AntisiteDefectTransformation, InterstitialTransformation

try:
    import zeo
except ImportError:
    zeo = None

class VacancyTransformationTest(unittest.TestCase):

    def test_apply_transformation(self):
        t = VacancyTransformation([2,2,2])
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        structures = t.apply_transformation(struct,return_ranked_list=100)
        self.assertEqual(len(structures),2)
        for struct in structures:
            self.assertIn(struct['structure'].composition.formula,
                          ["Li16 O16", "Li15 O16", "Li16 O15"])


class SubstitutionDefectTransformationTest(unittest.TestCase):

    def test_apply_transformation(self):
        t = SubstitutionDefectTransformation({"Li+":"Na+","O2-":"S2-"},[2,2,2])
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        scs = t.apply_transformation(struct,return_ranked_list=100)
        self.assertEqual(len(scs),2)
        for sc in scs:
            self.assertIn(sc['structure'].composition.formula,
                          ["Li16 O16", "Na1 Li15 O16", "Li16 S1 O15"])


class AntisiteDefectTransformationTest(unittest.TestCase):

    def test_apply_transformation(self):
        t = AntisiteDefectTransformation([2,2,2])
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        scs = t.apply_transformation(struct,return_ranked_list=100)
        self.assertEqual(len(scs),2)
        for sc in scs:
            self.assertIn(sc['structure'].composition.formula,
                          ["Li16 O16", "Li15 O17", "Li17 O15"])

@unittest.skipIf(not zeo, "zeo not present.")
class InterstitialTransformationTest(unittest.TestCase):

    def test_apply_transformation(self):
        t = InterstitialTransformation("Na+",[2,2,2])
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        scs = t.apply_transformation(struct,return_ranked_list=100000)
        #self.assertEqual(len(scs),3)
        for sc in scs:
            #print sc.composition.formula
            self.assertIn(sc['structure'].composition.formula,
                          ["Li16 O16", "Na1 Li16 O16", "Li16 Na1 O16"])


if __name__ == '__main__':
    unittest.main()
