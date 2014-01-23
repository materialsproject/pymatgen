import os
from unittest import TestCase

__author__ = 'xiaohuiqu'


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "molecules", "structural_change")

class TestMoleculeStructureComparator(TestCase):
    def test_are_equal(self):
        pass

    def test_get_bonds(self):
        print "hello"