# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest
import os
import pandas as pd
import numpy as np

from pymatgen.core.structure import Molecule
from pymatgen.io.xyz import XYZ
from pymatgen.io.vasp.inputs import Poscar

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class XYZTest(unittest.TestCase):

    def setUp(self):
        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]
        coords2 = [[x + 10.0 for x in atom] for atom in coords]
        self.mol = Molecule(["C", "H", "H", "H", "H"], coords)
        self.multi_mols = [Molecule(["C", "H", "H", "H", "H"], coords)
                           for coords in [coords, coords2]]
        self.xyz = XYZ(self.mol)
        self.multi_xyz = XYZ(self.multi_mols)

    def test_str(self):
        ans = """5
H4 C1
C 0.000000 0.000000 0.000000
H 0.000000 0.000000 1.089000
H 1.026719 0.000000 -0.363000
H -0.513360 -0.889165 -0.363000
H -0.513360 0.889165 -0.363000"""
        self.assertEqual(str(self.xyz), ans)

        mxyz = XYZ(self.multi_mols, coord_precision=3)
        mxyz_text = str(mxyz)
        ans_multi = """5
H4 C1
C 0.000 0.000 0.000
H 0.000 0.000 1.089
H 1.027 0.000 -0.363
H -0.513 -0.889 -0.363
H -0.513 0.889 -0.363
5
H4 C1
C 10.000 10.000 10.000
H 10.000 10.000 11.089
H 11.027 10.000 9.637
H 9.487 9.111 9.637
H 9.487 10.889 9.637"""
        self.assertEqual(mxyz_text, ans_multi)

    def test_from_string(self):
        ans = """5
H4 C1
C 0.000000 0.000000 0.000000
H 0.000000 0.000000 1.089000
H 1.026719 0.000000 -0.363000
H -0.513360 -0.889165 -0.363000
H -0.513360 0.889165 -0.363000"""
        xyz = XYZ.from_string(ans)
        mol = xyz.molecule
        sp = ["C", "H", "H", "H", "H"]
        for i, site in enumerate(mol):
            self.assertEqual(site.species_string, sp[i])
            self.assertEqual(len(site.coords), 3)
            if i == 0:
                self.assertTrue(all([c == 0 for c in site.coords]))

        mol_str = """2
Random
C 2.39132145462 -0.700993488928 -7.22293142224e-06
C 1.16730636786 -1.38166622735 -2.77112970359e-06
"""
        xyz = XYZ.from_string(mol_str)
        mol = xyz.molecule
        self.assertTrue(abs(mol[0].z) < 1e-5)
        self.assertTrue(abs(mol[1].z) < 1e-5)

        mol_str = """2
Random, Alternate Scientific Notation
C 2.39132145462 -0.700993488928 -7.222*^-06
C 1.16730636786 -1.38166622735 -2.771*^-06
"""
        xyz = XYZ.from_string(mol_str)
        mol = xyz.molecule
        self.assertEqual(mol[0].z, -7.222e-06)
        self.assertEqual(mol[1].z, -2.771e-06)

        mol_str = """3
Random
C   0.000000000000E+00  2.232615992397E+01  0.000000000000E+00
C  -2.383225420567E-31  1.116307996198E+01  1.933502166311E+01
C  -4.440892098501D-01 -1.116307996198d+01  1.933502166311E+01
"""
        xyz = XYZ.from_string(mol_str)
        mol = xyz.molecule
        self.assertAlmostEqual(mol[0].x, 0)
        self.assertAlmostEqual(mol[1].y, 11.16307996198)
        self.assertAlmostEqual(mol[2].x, -0.4440892098501)
        self.assertAlmostEqual(mol[2].y, -11.16307996198)
        # self.assertTrue(abs(mol[1].z) < 1e-5)

        mol_str = """    5
C32-C2-1                                                                        
 C     2.70450   1.16090  -0.14630     1     3    23     2
 C     1.61930   1.72490  -0.79330     2     1     5    26
 C     2.34210   1.02670   1.14620     3     1     8     6
 C    -0.68690   2.16170  -0.13790     4     5    18     7
 C     0.67160   2.15830   0.14350     5     4     2     6
 """
        xyz = XYZ.from_string(mol_str)
        mol = xyz.molecule
        self.assertAlmostEqual(mol[0].x, 2.70450)
        self.assertAlmostEqual(mol[1].y, 1.72490)
        self.assertAlmostEqual(mol[2].x, 2.34210)
        self.assertAlmostEqual(mol[3].z, -0.13790)

    def test_from_file(self):
        filepath = os.path.join(test_dir, 'multiple_frame_xyz.xyz')
        mxyz = XYZ.from_file(filepath)
        self.assertEqual(len(mxyz.all_molecules), 302)
        self.assertEqual(list(mxyz.all_molecules[0].cart_coords[0]),
                         [0.20303525080000001, 2.8569761204000002, 0.44737723190000001])
        self.assertEqual(list(mxyz.all_molecules[-1].cart_coords[-1]),
                         [5.5355550720000002, 0.0282305931, -0.30993102189999999])
        self.assertEqual(list(mxyz.molecule.cart_coords[-1]),
                         [5.5355550720000002, 0.0282305931, -0.30993102189999999])

    def test_init_from_structure(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        poscar = Poscar.from_file(filepath)
        struct = poscar.structure
        xyz = XYZ(struct)
        ans = """24
Fe4 P4 O16
Fe 2.277347 4.550379 2.260125
Fe 2.928536 1.516793 4.639870
Fe 7.483231 4.550379 0.119620
Fe 8.134420 1.516793 2.499364
P 0.985089 1.516793 1.990624
P 4.220794 4.550379 4.370369
P 6.190973 1.516793 0.389120
P 9.426677 4.550379 2.768865
O 0.451582 4.550379 3.365614
O 1.006219 1.516793 3.528306
O 1.725331 0.279529 1.358282
O 1.725331 2.754057 1.358282
O 3.480552 3.313115 3.738027
O 3.480552 5.787643 3.738027
O 4.199665 4.550379 1.148562
O 4.754301 1.516793 0.985870
O 5.657466 4.550379 3.773620
O 6.212102 1.516793 3.610928
O 6.931215 0.279529 1.021463
O 6.931215 2.754057 1.021463
O 8.686436 3.313115 3.401208
O 8.686436 5.787643 3.401208
O 9.405548 4.550379 1.231183
O 9.960184 1.516793 1.393875"""
        self.assertEqual(str(xyz), ans)

    def test_as_dataframe(self):
        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]
        test_df = pd.DataFrame(coords, columns=['x', 'y', 'z'])
        test_df.insert(0, "atom", ["C", "H", "H", "H", "H"])
        test_df.index += 1
        coords2 = [[0.000000, 0.000000, 0.000000],
                   [0.000000, 0.000000, 1.089000],
                   [1.026719, 0.000000, 0.363000],
                   [0.513360, 0.889165, 0.363000],
                   [0.513360, 0.889165, 0.363000]]
        test_df2 = pd.DataFrame(coords2, columns=['x', 'y', 'z'])
        test_df2.insert(0, "atom", ["C", "H", "H", "H", "H"])
        test_df2.index += 1
        mol_df = self.xyz.as_dataframe()

        # body tests
        pd.testing.assert_frame_equal(mol_df, test_df)

        # index tests
        np.testing.assert_array_equal(mol_df.columns, test_df.columns)
        np.testing.assert_array_equal(mol_df.index, test_df.index)


if __name__ == "__main__":
    unittest.main()
