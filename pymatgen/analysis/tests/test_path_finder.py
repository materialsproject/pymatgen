from __future__ import annotations

import os
import unittest

from numpy import mean

from pymatgen.analysis.path_finder import ChgcarPotential, NEBPathfinder
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp import Chgcar, Poscar
from pymatgen.util.testing import PymatgenTest

__author__ = "Ziqin (Shaun) Rong"
__version__ = "0.1"
__maintainer__ = "Ziqin (Shaun) Rong"
__email__ = "rongzq08@gmail.com"


class PathFinderTest(unittest.TestCase):
    """
    Uses Li migration in LiFePO4
    """

    def test_image_num(self):
        os.path.dirname(os.path.abspath(__file__))
        test_file_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "path_finder")
        start_s = Poscar.from_file(os.path.join(test_file_dir, "LFP_POSCAR_s")).structure
        end_s = Poscar.from_file(os.path.join(test_file_dir, "LFP_POSCAR_e")).structure
        chg = Chgcar.from_file(os.path.join(test_file_dir, "LFP_CHGCAR.gz"))
        moving_cation_specie = Element("Li")
        relax_sites = []
        for site_i, site in enumerate(start_s.sites):
            if site.specie == moving_cation_specie:
                relax_sites.append(site_i)
        pf = NEBPathfinder(
            start_s,
            end_s,
            relax_sites=relax_sites,
            v=ChgcarPotential(chg).get_v(),
            n_images=(8 * 3),
        )
        images = []
        for i, image in enumerate(pf.images):
            if i % 3 == 0:
                images.append(image)
        self.assertEqual(len(images), 9)

        moving_site = relax_sites[0]
        dists = [s1.sites[moving_site].distance(s2.sites[moving_site]) for s1, s2 in zip(pf.images[:-1], pf.images[1:])]
        # check that all the small distances are about equal
        self.assertTrue(abs(min(dists) - max(dists)) / mean(dists) < 0.02)


if __name__ == "__main__":
    unittest.main()
