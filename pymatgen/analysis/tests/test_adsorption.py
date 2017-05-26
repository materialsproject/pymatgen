from __future__ import absolute_import
from __future__ import print_function

import unittest
import os

import numpy as np
from pymatgen.util.testing import PymatgenTest
from pymatgen.util.coord_utils import in_coord_list, in_coord_list_pbc
from pymatgen.core.surface import generate_all_slabs
from pymatgen.analysis.adsorption import *
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen import Structure, Lattice, Molecule
import json
from six.moves import zip

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')


class AdsorbateSiteFinderTest(PymatgenTest):
    def setUp(self):
        self.structure = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3.5),
                                                ["Ni"], [[0, 0, 0]])
        slabs = generate_all_slabs(self.structure, max_index=2,
                                   min_slab_size=6.0, min_vacuum_size=15.0,
                                   max_normal_search=1, center_slab=True)
        self.slab_dict = {''.join([str(i) for i in slab.miller_index]):
                          slab for slab in slabs}
        self.asf_211 = AdsorbateSiteFinder(self.slab_dict["211"])
        self.asf_100 = AdsorbateSiteFinder(self.slab_dict["100"])
        self.asf_111 = AdsorbateSiteFinder(self.slab_dict["111"])
        self.asf_110 = AdsorbateSiteFinder(self.slab_dict["110"])

    def test_init(self):
        asf_100 = AdsorbateSiteFinder(self.slab_dict["100"])
        asf_111 = AdsorbateSiteFinder(self.slab_dict["111"])

    def test_from_bulk_and_miller(self):
        # Standard site finding
        asf = AdsorbateSiteFinder.from_bulk_and_miller(self.structure, (1, 1, 1))
        sites = asf.find_adsorption_sites()
        self.assertEqual(len(sites['hollow']), 2)
        self.assertEqual(len(sites['bridge']), 1)
        self.assertEqual(len(sites['ontop']), 1)
        self.assertEqual(len(sites['all']), 4)
        asf = AdsorbateSiteFinder.from_bulk_and_miller(self.structure, (1, 0, 0))
        sites = asf.find_adsorption_sites()
        self.assertEqual(len(sites['all']), 3)
        self.assertEqual(len(sites['bridge']), 2)
        asf = AdsorbateSiteFinder.from_bulk_and_miller(self.structure, (1, 1, 0),
                                                       undercoord_threshold=0.1)
        self.assertEqual(len(asf.surface_sites), 1)
        # Subsurface site finding
        asf = AdsorbateSiteFinder.from_bulk_and_miller(self.structure, (1, 1, 1))
        sites = asf.find_adsorption_sites(positions=["ontop", "subsurface", "bridge"])
        self.assertEqual(len(sites['all']), 4)
        self.assertEqual(len(sites['subsurface']), 3)

    def test_find_adsorption_sites(self):
        sites = self.asf_100.find_adsorption_sites()
        self.assertEqual(len(sites['all']), 3)
        self.assertEqual(len(sites['hollow']), 0)
        self.assertEqual(len(sites['bridge']), 2)
        self.assertEqual(len(sites['ontop']), 1)
        sites = self.asf_111.find_adsorption_sites()
        self.assertEqual(len(sites['all']), 4)
        sites = self.asf_110.find_adsorption_sites()
        self.assertEqual(len(sites['all']), 4)
        sites = self.asf_211.find_adsorption_sites()

    def test_generate_adsorption_structures(self):
        co = Molecule("CO", [[0, 0, 0], [0, 0, 1.23]])
        structures = self.asf_111.generate_adsorption_structures(co, repeat=[2, 2, 1])
        self.assertEqual(len(structures), 4)
        sites = self.asf_111.find_adsorption_sites()
        # Check repeat functionality
        self.assertEqual(len([site for site in structures[0] if 
                              site.properties['surface_properties'] != 'adsorbate']),
                         4*len(self.asf_111.slab))
        for n, structure in enumerate(structures):
            self.assertArrayAlmostEqual(structure[-2].coords, sites['all'][n])
        find_args = {"positions":["hollow"]}
        structures_hollow = self.asf_111.\
                generate_adsorption_structures(co, find_args=find_args)
        self.assertEqual(len(structures_hollow), len(sites['hollow']))
        for n, structure in enumerate(structures_hollow):
            self.assertTrue(in_coord_list(sites['hollow'], structure[-2].coords))

    def test_functions(self):
        slab = self.slab_dict["111"]
        rot = get_rot(slab)
        reoriented = reorient_z(slab)
        self.assertArrayAlmostEqual(slab.frac_coords[0],
                                    cart_to_frac(slab.lattice, 
                                                 slab.cart_coords[0]))
        self.assertArrayAlmostEqual(slab.cart_coords[0],
                                    frac_to_cart(slab.lattice,
                                                 slab.frac_coords[0]))

if __name__ == '__main__':
    unittest.main()
