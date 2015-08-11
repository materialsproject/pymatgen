# coding: utf-8

from __future__ import unicode_literals

import unittest
import sys

from pymatgen.analysis.defects.point_defects import *
from pymatgen.analysis.defects.charge_defects_analyzer import *
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.analysis.bond_valence import BVAnalyzer
from monty.os.path import which
from pymatgen.io.cifio import CifParser

sxdefectalign_present = which('sxdefectalign')
test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')


class ParsedChargeDefectTest(unittest.TestCase):
    """
    This class is just an aggregate for the data. No tests are necessary.
    """
    pass


@unittest.skipIf(not sxdefectalign_present, "sxdefectalign not present.")
class CorrectionTest(unittest.TestCase):
    """
    Tests for the correction obtained Freysoldt method
    """
    def setUp(self):
        pass

    def test_get_correction(self):
        pass


class ChargeDefectsAnalyzerTest(unittest.TestCase):
    def setUp(self):
        pass

    def test_as_dict(self):
        pass

    def test_from_dict(self):
        pass

    def test_add_parsed_defect(self):
        pass

    def test_change_charge_correction(self):
        pass

    def test_correct_bg_simple(self):
        pass

    def test_correct_bg(self):
        pass

    def test_get_defects_concentration(self):
        pass

    def test_get_eq_ef(self):
        pass

    def test_get_non_eq_ef(self):
        pass


if __name__ == "__main__":
    unittest.main()
