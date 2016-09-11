#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import unittest

from pymatgen.analysis.path_finder import NEBPathfinder, ChgcarPotential
from pymatgen.io.vasp import Poscar, Chgcar, Element

__author__ = 'Ziqin (Shaun) Rong'
__version__ = '0.1'
__maintainer__ = 'Ziqin (Shaun) Rong'
__email__ = 'rongzq08@gmail.com'


class PathFinderTest(unittest.TestCase):
    """
    Uses Li migration in LiFePO4
    """
    def setUp(self):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        test_file_dir = os.path.join(module_dir, "..", "..", "..", "test_files",
                                     "path_finder")
        self.start_s = Poscar.from_file(os.path.join(test_file_dir, 'LFP_POSCAR_s')).structure
        self.end_s = Poscar.from_file(os.path.join(test_file_dir, 'LFP_POSCAR_e')).structure
        self.chg = Chgcar.from_file(os.path.join(test_file_dir, 'LFP_CHGCAR'))
        moving_cation_specie = Element('Li')
        self.relax_sites = []
        for site_i, site in enumerate(self.start_s.sites):
            if site.specie == moving_cation_specie:
                self.relax_sites.append(site_i)
        self.pf = NEBPathfinder(self.start_s, self.end_s, relax_sites=self.relax_sites,
                                v=ChgcarPotential(self.chg).get_v(), n_images=(8 * 3))
        self.images = []
        for i, image in enumerate(self.pf.images):
            if i % 3 == 0:
                self.images.append(image)

    def test_image_num(self):
        self.assertEqual(len(self.images), 9)

if __file__ == '__main__':
    unittest.main()
