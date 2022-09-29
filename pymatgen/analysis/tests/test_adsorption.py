import unittest

import numpy as np

from pymatgen.analysis.adsorption import (
    AdsorbateSiteFinder,
    generate_all_slabs,
    get_rot,
    reorient_z,
)
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Molecule, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord import in_coord_list
from pymatgen.util.testing import PymatgenTest


class AdsorbateSiteFinderTest(PymatgenTest):
    def setUp(self):
        self.structure = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3.5), ["Ni"], [[0, 0, 0]])
        lattice = Lattice.cubic(3.010)
        frac_coords = [
            [0.00000, 0.00000, 0.00000],
            [0.00000, 0.50000, 0.50000],
            [0.50000, 0.00000, 0.50000],
            [0.50000, 0.50000, 0.00000],
            [0.50000, 0.00000, 0.00000],
            [0.50000, 0.50000, 0.50000],
            [0.00000, 0.00000, 0.50000],
            [0.00000, 0.50000, 0.00000],
        ]
        species = ["Mg", "Mg", "Mg", "Mg", "O", "O", "O", "O"]
        self.MgO = Structure(lattice, species, frac_coords)

        slabs = generate_all_slabs(
            self.structure,
            max_index=2,
            min_slab_size=6.0,
            min_vacuum_size=15.0,
            max_normal_search=1,
            center_slab=True,
        )
        self.slab_dict = {"".join([str(i) for i in slab.miller_index]): slab for slab in slabs}
        self.asf_211 = AdsorbateSiteFinder(self.slab_dict["211"])
        self.asf_100 = AdsorbateSiteFinder(self.slab_dict["100"])
        self.asf_111 = AdsorbateSiteFinder(self.slab_dict["111"])
        self.asf_110 = AdsorbateSiteFinder(self.slab_dict["110"])
        self.asf_struct = AdsorbateSiteFinder(Structure.from_sites(self.slab_dict["111"].sites))

    def test_init(self):
        AdsorbateSiteFinder(self.slab_dict["100"])
        AdsorbateSiteFinder(self.slab_dict["111"])

    def test_from_bulk_and_miller(self):
        # Standard site finding
        asf = AdsorbateSiteFinder.from_bulk_and_miller(self.structure, (1, 1, 1))
        sites = asf.find_adsorption_sites()
        self.assertEqual(len(sites["hollow"]), 2)
        self.assertEqual(len(sites["bridge"]), 1)
        self.assertEqual(len(sites["ontop"]), 1)
        self.assertEqual(len(sites["all"]), 4)
        asf = AdsorbateSiteFinder.from_bulk_and_miller(self.structure, (1, 0, 0))
        sites = asf.find_adsorption_sites()
        self.assertEqual(len(sites["all"]), 3)
        self.assertEqual(len(sites["bridge"]), 2)
        asf = AdsorbateSiteFinder.from_bulk_and_miller(self.structure, (1, 1, 0), undercoord_threshold=0.1)
        self.assertEqual(len(asf.surface_sites), 1)
        # Subsurface site finding
        asf = AdsorbateSiteFinder.from_bulk_and_miller(self.structure, (1, 1, 1))
        sites = asf.find_adsorption_sites(positions=["ontop", "subsurface", "bridge"])
        self.assertEqual(len(sites["all"]), 5)
        self.assertEqual(len(sites["subsurface"]), 3)

    def test_find_adsorption_sites(self):
        sites = self.asf_100.find_adsorption_sites()
        self.assertEqual(len(sites["all"]), 3)
        self.assertEqual(len(sites["hollow"]), 0)
        self.assertEqual(len(sites["bridge"]), 2)
        self.assertEqual(len(sites["ontop"]), 1)
        sites = self.asf_111.find_adsorption_sites()
        self.assertEqual(len(sites["all"]), 4)
        sites = self.asf_110.find_adsorption_sites()
        self.assertEqual(len(sites["all"]), 4)
        sites = self.asf_211.find_adsorption_sites()
        # Test on structure
        sites = self.asf_struct.find_adsorption_sites()

    def test_generate_adsorption_structures(self):
        co = Molecule("CO", [[0, 0, 0], [0, 0, 1.23]])
        structures = self.asf_111.generate_adsorption_structures(co, repeat=[2, 2, 1])
        self.assertEqual(len(structures), 4)
        sites = self.asf_111.find_adsorption_sites()
        # Check repeat functionality
        self.assertEqual(
            len([site for site in structures[0] if site.properties["surface_properties"] != "adsorbate"]),
            4 * len(self.asf_111.slab),
        )
        for n, structure in enumerate(structures):
            self.assertArrayAlmostEqual(structure[-2].coords, sites["all"][n])
        find_args = {"positions": ["hollow"]}
        structures_hollow = self.asf_111.generate_adsorption_structures(co, find_args=find_args)
        self.assertEqual(len(structures_hollow), len(sites["hollow"]))
        for structure in structures_hollow:
            self.assertTrue(in_coord_list(sites["hollow"], structure[-2].coords, 1e-4))
        # Check molecule not changed after rotation when added to surface
        co = Molecule("CO", [[1.0, -0.5, 3], [0.8, 0.46, 3.75]])
        structures = self.asf_211.generate_adsorption_structures(co)
        self.assertEqual(co, Molecule("CO", [[1.0, -0.5, 3], [0.8, 0.46, 3.75]]))
        # Check translation
        sites = self.asf_211.find_adsorption_sites()
        ads_site_coords = sites["all"][0]
        c_site = structures[0].sites[-2]
        self.assertEqual(str(c_site.specie), "C")
        self.assertArrayAlmostEqual(c_site.coords, sites["all"][0])
        # Check no translation
        structures = self.asf_111.generate_adsorption_structures(co, translate=False)
        self.assertEqual(co, Molecule("CO", [[1.0, -0.5, 3], [0.8, 0.46, 3.75]]))
        sites = self.asf_111.find_adsorption_sites()
        ads_site_coords = sites["all"][0]
        c_site = structures[0].sites[-2]
        self.assertArrayAlmostEqual(c_site.coords, ads_site_coords + np.array([1.0, -0.5, 3]))

    def test_adsorb_both_surfaces(self):

        # Test out for monatomic adsorption
        o = Molecule("O", [[0, 0, 0]])
        adslabs = self.asf_100.adsorb_both_surfaces(o)
        adslabs_one = self.asf_100.generate_adsorption_structures(o)
        self.assertEqual(len(adslabs), len(adslabs_one))
        for adslab in adslabs:
            sg = SpacegroupAnalyzer(adslab)
            sites = sorted(adslab, key=lambda site: site.frac_coords[2])
            self.assertTrue(sites[0].species_string == "O")
            self.assertTrue(sites[-1].species_string == "O")
            self.assertTrue(sg.is_laue())

        # Test out for molecular adsorption
        oh = Molecule(["O", "H"], [[0, 0, 0], [0, 0, 1]])
        adslabs = self.asf_100.adsorb_both_surfaces(oh)
        adslabs_one = self.asf_100.generate_adsorption_structures(oh)
        self.assertEqual(len(adslabs), len(adslabs_one))
        for adslab in adslabs:
            sg = SpacegroupAnalyzer(adslab)
            sites = sorted(adslab, key=lambda site: site.frac_coords[2])
            self.assertTrue(sites[0].species_string in ["O", "H"])
            self.assertTrue(sites[-1].species_string in ["O", "H"])
            self.assertTrue(sg.is_laue())

    def test_generate_substitution_structures(self):

        # Test this for a low miller index halite structure
        slabs = generate_all_slabs(self.MgO, 1, 10, 10, center_slab=True, max_normal_search=1)
        for slab in slabs:
            adsgen = AdsorbateSiteFinder(slab)

            adslabs = adsgen.generate_substitution_structures("Ni")
            # There should be 2 configs (sub O and sub
            # Mg) for (110) and (100), 1 for (111)
            if tuple(slab.miller_index) != (1, 1, 1):
                self.assertEqual(len(adslabs), 2)
            else:
                self.assertEqual(len(adslabs), 1)

            # Test out whether it can correctly dope both
            # sides. Avoid (111) because it is not symmetric
            if tuple(slab.miller_index) != (1, 1, 1):
                adslabs = adsgen.generate_substitution_structures("Ni", sub_both_sides=True, target_species=["Mg"])
                # Test if default parameters dope the surface site
                for i, site in enumerate(adslabs[0]):
                    if adsgen.slab[i].surface_properties == "surface" and site.species_string == "Mg":
                        print(
                            adslabs[0][i].surface_properties,
                            adsgen.slab[i].surface_properties,
                        )
                        self.assertTrue(adslabs[0][i].surface_properties == "substitute")

                self.assertTrue(adslabs[0].is_symmetric())
                # Correctly dope the target species
                self.assertEqual(
                    adslabs[0].composition.as_dict()["Mg"],
                    slab.composition.as_dict()["Mg"] - 2,
                )
                # There should be one config (sub Mg)
                self.assertEqual(len(adslabs), 1)

    def test_functions(self):
        slab = self.slab_dict["111"]
        get_rot(slab)
        reorient_z(slab)


if __name__ == "__main__":
    unittest.main()
