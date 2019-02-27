# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from pymatgen.core import Specie, Element, Lattice, Structure
from pymatgen.io.cif import CifParser
from pymatgen.analysis.magnetism import *

import os
import unittest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')


class CollinearMagneticStructureAnalyzerTest(unittest.TestCase):

    def setUp(self):
        parser = CifParser(os.path.join(test_dir, 'Fe.cif'))
        self.Fe = parser.get_structures()[0]

        parser = CifParser(os.path.join(test_dir, 'LiFePO4.cif'))
        self.LiFePO4 = parser.get_structures()[0]

        parser = CifParser(os.path.join(test_dir, 'Fe3O4.cif'))
        self.Fe3O4 = parser.get_structures()[0]

        parser = CifParser(os.path.join(test_dir, 'magnetic.ncl.example.GdB4.mcif'))
        self.GdB4 = parser.get_structures()[0]

        parser = CifParser(os.path.join(test_dir, 'magnetic.example.NiO.mcif'))
        self.NiO_expt = parser.get_structures()[0]

        latt = Lattice.cubic(4.17)
        species = ["Ni", "O"]
        coords = [[0, 0, 0],
                  [0.5, 0.5, 0.5]]
        self.NiO = Structure.from_spacegroup(225, latt, species, coords)

        latt = Lattice([[2.085, 2.085, 0.0],
                        [0.0, -2.085, -2.085],
                        [-2.085, 2.085, -4.17]])
        species = ["Ni", "Ni", "O", "O"]
        coords = [[0.5, 0, 0.5],
                  [0, 0, 0],
                  [0.25, 0.5, 0.25],
                  [0.75, 0.5, 0.75]]
        self.NiO_AFM_111 = Structure(latt, species, coords,
                                     site_properties={'magmom': [-5, 5, 0, 0]})

        latt = Lattice([[2.085, 2.085, 0],
                        [0, 0, -4.17],
                        [-2.085, 2.085, 0]])
        species = ["Ni", "Ni", "O", "O"]
        coords = [[0.5, 0.5, 0.5],
                  [0, 0, 0],
                  [0, 0.5, 0],
                  [0.5, 0, 0.5]]
        self.NiO_AFM_001 = Structure(latt, species, coords,
                                     site_properties={'magmom': [-5, 5, 0, 0]})

        latt = Lattice([[2.085, 2.085, 0],
                        [0, 0, -4.17],
                        [-2.085, 2.085, 0]])
        species = ["Ni", "Ni", "O", "O"]
        coords = [[0.5, 0.5, 0.5],
                  [0, 0, 0],
                  [0, 0.5, 0],
                  [0.5, 0, 0.5]]
        self.NiO_AFM_001_opposite = Structure(latt, species, coords,
                                              site_properties={'magmom': [5, -5, 0, 0]})

        latt = Lattice([[2.085, 2.085, 0],
                        [0, 0, -4.17],
                        [-2.085, 2.085, 0]])
        species = ["Ni", "Ni", "O", "O"]
        coords = [[0.5, 0.5, 0.5],
                  [0, 0, 0],
                  [0, 0.5, 0],
                  [0.5, 0, 0.5]]
        self.NiO_unphysical = Structure(latt, species, coords,
                                        site_properties={'magmom': [-3, 0, 0, 0]})

        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_get_representations(self):

        # tests to convert between storing magnetic moment information
        # on site_properties or on Specie 'spin' property

        # test we store magnetic moments on site properties
        self.Fe.add_site_property('magmom', [5])
        msa = CollinearMagneticStructureAnalyzer(self.Fe)
        self.assertEqual(msa.structure.site_properties['magmom'][0], 5)

        # and that we can retrieve a spin representaiton
        Fe_spin = msa.get_structure_with_spin()
        self.assertFalse('magmom' in Fe_spin.site_properties)
        self.assertEqual(Fe_spin[0].specie.spin, 5)

        # test we can remove magnetic moment information
        Fe_none = msa.get_nonmagnetic_structure()
        self.assertFalse('magmom' in Fe_spin.site_properties)

        # test with disorder on magnetic site
        self.Fe[0] = {Specie('Fe', oxidation_state=0, properties={'spin': 5}): 0.5, 'Ni': 0.5}
        self.assertRaises(NotImplementedError, CollinearMagneticStructureAnalyzer, self.Fe)

    def test_matches(self):

        self.assertTrue(self.NiO.matches(self.NiO_AFM_111))
        self.assertTrue(self.NiO.matches(self.NiO_AFM_001))

        # MSA adds magmoms to Structure, so not equal
        msa = CollinearMagneticStructureAnalyzer(self.NiO,
                                                 overwrite_magmom_mode="replace_all")
        self.assertFalse(msa.matches_ordering(self.NiO))
        self.assertFalse(msa.matches_ordering(self.NiO_AFM_111))
        self.assertFalse(msa.matches_ordering(self.NiO_AFM_001))

        msa = CollinearMagneticStructureAnalyzer(self.NiO_AFM_001,
                                                 overwrite_magmom_mode="respect_sign")
        self.assertFalse(msa.matches_ordering(self.NiO))
        self.assertFalse(msa.matches_ordering(self.NiO_AFM_111))
        self.assertTrue(msa.matches_ordering(self.NiO_AFM_001))
        self.assertTrue(msa.matches_ordering(self.NiO_AFM_001_opposite))

        msa = CollinearMagneticStructureAnalyzer(self.NiO_AFM_111,
                                                 overwrite_magmom_mode="respect_sign")
        self.assertFalse(msa.matches_ordering(self.NiO))
        self.assertTrue(msa.matches_ordering(self.NiO_AFM_111))
        self.assertFalse(msa.matches_ordering(self.NiO_AFM_001))
        self.assertFalse(msa.matches_ordering(self.NiO_AFM_001_opposite))

    def test_modes(self):

        mode = "none"
        msa = CollinearMagneticStructureAnalyzer(self.NiO,
                                                 overwrite_magmom_mode=mode)
        magmoms = msa.structure.site_properties['magmom']
        self.assertEqual(magmoms, [0, 0])

        mode = "respect_sign"
        msa = CollinearMagneticStructureAnalyzer(self.NiO_unphysical,
                                                 overwrite_magmom_mode=mode)
        magmoms = msa.structure.site_properties['magmom']
        self.assertEqual(magmoms, [-5, 0, 0, 0])

        mode = "respect_zeros"
        msa = CollinearMagneticStructureAnalyzer(self.NiO_unphysical,
                                                 overwrite_magmom_mode=mode)
        magmoms = msa.structure.site_properties['magmom']
        self.assertEqual(magmoms, [5, 0, 0, 0])

        mode = "replace_all"
        msa = CollinearMagneticStructureAnalyzer(self.NiO_unphysical,
                                                 overwrite_magmom_mode=mode,
                                                 make_primitive=False)
        magmoms = msa.structure.site_properties['magmom']
        self.assertEqual(magmoms, [5, 5, 0, 0])

        mode = "replace_all_if_undefined"
        msa = CollinearMagneticStructureAnalyzer(self.NiO,
                                                 overwrite_magmom_mode=mode)
        magmoms = msa.structure.site_properties['magmom']
        self.assertEqual(magmoms, [5, 0])

        mode = "normalize"
        msa = CollinearMagneticStructureAnalyzer(msa.structure,
                                                 overwrite_magmom_mode='normalize')
        magmoms = msa.structure.site_properties['magmom']
        self.assertEqual(magmoms, [1, 0])

    def test_net_positive(self):
        msa = CollinearMagneticStructureAnalyzer(self.NiO_unphysical)
        magmoms = msa.structure.site_properties['magmom']
        self.assertEqual(magmoms, [3, 0, 0, 0])

    def test_get_ferromagnetic_structure(self):

        msa = CollinearMagneticStructureAnalyzer(self.NiO,
                                                 overwrite_magmom_mode="replace_all_if_undefined")
        s1 = msa.get_ferromagnetic_structure()
        s1_magmoms = [float(m) for m in s1.site_properties['magmom']]
        s1_magmoms_ref = [5.0, 0.0]
        self.assertListEqual(s1_magmoms, s1_magmoms_ref)

        msa2 = CollinearMagneticStructureAnalyzer(self.NiO_AFM_111,
                                                 overwrite_magmom_mode="replace_all_if_undefined")
        s2 = msa.get_ferromagnetic_structure(make_primitive=False)
        s2_magmoms = [float(m) for m in s2.site_properties['magmom']]
        s2_magmoms_ref = [5.0, 0.0]
        self.assertListEqual(s2_magmoms, s2_magmoms_ref)

        s2_prim = msa.get_ferromagnetic_structure(make_primitive=True)
        self.assertTrue(CollinearMagneticStructureAnalyzer(s1).matches_ordering(s2_prim))

    def test_magnetic_properties(self):

        msa = CollinearMagneticStructureAnalyzer(self.GdB4)
        self.assertFalse(msa.is_collinear)

        msa = CollinearMagneticStructureAnalyzer(self.Fe)
        self.assertFalse(msa.is_magnetic)

        self.Fe.add_site_property('magmom', [5])

        msa = CollinearMagneticStructureAnalyzer(self.Fe)
        self.assertTrue(msa.is_magnetic)
        self.assertTrue(msa.is_collinear)
        self.assertEqual(msa.ordering, Ordering.FM)

        msa = CollinearMagneticStructureAnalyzer(self.NiO, make_primitive=False,
                                                 overwrite_magmom_mode="replace_all_if_undefined")
        self.assertEqual(msa.number_of_magnetic_sites, 4)
        self.assertEqual(msa.number_of_unique_magnetic_sites(), 1)
        self.assertEqual(msa.types_of_magnetic_specie, [Element('Ni')])
        self.assertEqual(msa.get_exchange_group_info(), ('Fm-3m', 225))

    def test_str(self):

        msa = CollinearMagneticStructureAnalyzer(self.NiO_AFM_001)

        ref_msa_str = """Structure Summary
Lattice
    abc : 2.948635277547903 4.17 2.948635277547903
 angles : 90.0 90.0 90.0
 volume : 36.2558565
      A : 2.085 2.085 0.0
      B : 0.0 0.0 -4.17
      C : -2.085 2.085 0.0
Magmoms Sites
+5.00   PeriodicSite: Ni (0.0000, 0.0000, 0.0000) [0.0000, 0.0000, 0.0000]
        PeriodicSite: O (0.0000, 0.0000, -2.0850) [0.0000, 0.5000, 0.0000]
        PeriodicSite: O (0.0000, 2.0850, 0.0000) [0.5000, 0.0000, 0.5000]
-5.00   PeriodicSite: Ni (0.0000, 2.0850, -2.0850) [0.5000, 0.5000, 0.5000]"""

        # just compare lines form 'Magmoms Sites',
        # since lattice param string can vary based on machine precision
        self.assertEqual("\n".join(str(msa).split("\n")[-5:-1]),
                         "\n".join(ref_msa_str.split("\n")[-5:-1]))

    def test_round_magmoms(self):

        struct = self.NiO_AFM_001.copy()
        struct.add_site_property('magmom', [-5.0143, -5.02, 0.147, 0.146])

        msa = CollinearMagneticStructureAnalyzer(struct, round_magmoms=0.001, make_primitive=False)
        self.assertTrue(np.allclose(msa.magmoms, [5.0171, 5.0171, -0.1465, -0.1465]))
        self.assertAlmostEqual(msa.magnetic_species_and_magmoms['Ni'], 5.0171)
        self.assertAlmostEqual(msa.magnetic_species_and_magmoms['O'], 0.1465)

        struct.add_site_property('magmom', [-5.0143, 4.5, 0.147, 0.146])
        msa = CollinearMagneticStructureAnalyzer(struct, round_magmoms=0.001, make_primitive=False)
        self.assertTrue(np.allclose(msa.magmoms, [5.0143, -4.5, -0.1465, -0.1465]))
        self.assertAlmostEqual(msa.magnetic_species_and_magmoms['Ni'][0], 4.5)
        self.assertAlmostEqual(msa.magnetic_species_and_magmoms['Ni'][1], 5.0143)
        self.assertAlmostEqual(msa.magnetic_species_and_magmoms['O'], 0.1465)
        

class MagneticDeformationTest(unittest.TestCase):

    def test_magnetic_deformation(self):

        test_structs = loadfn(os.path.join(test_dir, 'magnetic_deformation.json'))
        mag_def = magnetic_deformation(test_structs[0], test_structs[1])

        self.assertEqual(mag_def.type, "NM-FM")
        self.assertAlmostEqual(mag_def.deformation, 5.0130859485170971)

if __name__ == '__main__':
    unittest.main()
