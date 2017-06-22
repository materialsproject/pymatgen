# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import unittest

from pymatgen.ext.cod import get_structure


"""
Created on Jun 9, 2012
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jun 9, 2012"


class CODTest(unittest.TestCase):

    def test_get_structure(self):
        s = get_structure(2002926)
        self.assertEqual(s.formula, "Be8 H64 N16 F32")


if __name__ == "__main__":
    unittest.main()