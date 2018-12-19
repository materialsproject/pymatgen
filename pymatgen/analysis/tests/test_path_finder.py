#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
    def test_image_num(self):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        test_file_dir = os.path.join(module_dir, "..", "..", "..", "test_files",
                                     "path_finder")
        start_s = Poscar.from_file(os.path.join(test_file_dir, 'LFP_POSCAR_s')).structure
        end_s = Poscar.from_file(os.path.join(test_file_dir, 'LFP_POSCAR_e')).structure
        chg = Chgcar.from_file(os.path.join(test_file_dir, 'LFP_CHGCAR.gz'))
        moving_cation_specie = Element('Li')
        relax_sites = []
        for site_i, site in enumerate(start_s.sites):
            if site.specie == moving_cation_specie:
                relax_sites.append(site_i)
        pf = NEBPathfinder(start_s, end_s, relax_sites=relax_sites,
                           v=ChgcarPotential(chg).get_v(), n_images=(8 * 3))
        images = []
        for i, image in enumerate(pf.images):
            if i % 3 == 0:
                images.append(image)
        self.assertEqual(len(images), 9)


if __name__ == '__main__':
    unittest.main()
