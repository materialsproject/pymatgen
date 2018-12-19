# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
Created on Apr 25, 2012
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Apr 25, 2012"

# import unittest
#
# from pymatgen.vis.structure_vtk import StructureVis
# from pymatgen.core.structure import Structure
# from mock import MagicMock
#
#
# class StructureVisTest(unittest.TestCase):
#
#     def setUp(self):
#         self.vis = StructureVis()
#
#     def test_set_structure(self):
#         coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
#         structure = Structure([[3.8401979337, 0.00, 0.00],
#                                [1.9200989668, 3.3257101909, 0.00],
#                                [0.00, -2.2171384943, 3.1355090603]],
#                               ["Si"] * 2, coords)
#         self.vis.add_site = MagicMock(return_value=3)
#         self.vis.set_structure(structure)
#         self.assertEqual(self.vis.add_site.call_count, 2)
#
#
# if __name__ == "__main__":
#     unittest.main()
