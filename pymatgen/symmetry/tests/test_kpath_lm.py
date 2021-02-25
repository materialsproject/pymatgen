# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import os
import unittest

import numpy as np

from pymatgen.analysis.magnetism.analyzer import CollinearMagneticStructureAnalyzer
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.kpath import KPathLatimerMunro
from pymatgen.util.testing import PymatgenTest


test_dir_structs = PymatgenTest.TEST_FILES_DIR


class KPathLatimerMunroTest(PymatgenTest):
    def test_kpath_generation(self):
        triclinic = [1, 2]
        monoclinic = range(3, 16)
        orthorhombic = range(16, 75)
        tetragonal = range(75, 143)
        rhombohedral = range(143, 168)
        hexagonal = range(168, 195)
        cubic = range(195, 231)

        species = ["K", "La", "Ti"]
        coords = [[0.345, 5, 0.77298], [0.1345, 5.1, 0.77298], [0.7, 0.8, 0.9]]
        for i in range(230):
            sg_num = i + 1
            if sg_num in triclinic:
                lattice = Lattice(
                    [
                        [3.0233057319441246, 1, 0],
                        [0, 7.9850357844548681, 1],
                        [0, 1.2, 8.1136762279561818],
                    ]
                )
            elif sg_num in monoclinic:
                lattice = Lattice.monoclinic(2, 9, 1, 99)
            elif sg_num in orthorhombic:
                lattice = Lattice.orthorhombic(2, 9, 1)
            elif sg_num in tetragonal:
                lattice = Lattice.tetragonal(2, 9)
            elif sg_num in rhombohedral:
                lattice = Lattice.hexagonal(2, 95)
            elif sg_num in hexagonal:
                lattice = Lattice.hexagonal(2, 9)
            elif sg_num in cubic:
                lattice = Lattice.cubic(2)

            struct = Structure.from_spacegroup(sg_num, lattice, species, coords)
            kpath = KPathLatimerMunro(struct)  # Throws error if something doesn't work, causing test to fail.

        struct_file_path = os.path.join(test_dir_structs, "AgO_kpath_test.cif")
        struct = Structure.from_file(struct_file_path)
        kpath = KPathLatimerMunro(struct)  # Throws error if something doesn't work, causing test to fail.

    def test_kpath_acentered(self):
        species = ["K", "La", "Ti"]
        coords = [[0.345, 5, 0.77298], [0.1345, 5.1, 0.77298], [0.7, 0.8, 0.9]]
        lattice = Lattice.orthorhombic(2, 9, 1)
        struct = Structure.from_spacegroup(38, lattice, species, coords)
        sga = SpacegroupAnalyzer(struct)
        struct_prim = sga.get_primitive_standard_structure(international_monoclinic=False)
        kpath = KPathLatimerMunro(struct_prim)

        kpoints = kpath._kpath["kpoints"]
        labels = list(kpoints.keys())

        self.assertEqual(
            sorted(labels),
            sorted(["a", "b", "c", "d", "d_{1}", "e", "f", "q", "q_{1}", "Γ"]),
        )

        self.assertAlmostEqual(kpoints["a"][0], 0.0)
        self.assertAlmostEqual(kpoints["a"][1], 0.4999999999999997)
        self.assertAlmostEqual(kpoints["a"][2], 0.0)

        self.assertAlmostEqual(kpoints["f"][0], -0.49999999999999933)
        self.assertAlmostEqual(kpoints["f"][1], 0.4999999999999992)
        self.assertAlmostEqual(kpoints["f"][2], 0.4999999999999999)

        self.assertAlmostEqual(kpoints["c"][0], 0.0)
        self.assertAlmostEqual(kpoints["c"][1], 0.0)
        self.assertAlmostEqual(kpoints["c"][2], 0.5)

        self.assertAlmostEqual(kpoints["b"][0], -0.5000000000000002)
        self.assertAlmostEqual(kpoints["b"][1], 0.500000000000000)
        self.assertAlmostEqual(kpoints["b"][2], 0.0)

        self.assertAlmostEqual(kpoints["Γ"][0], 0)
        self.assertAlmostEqual(kpoints["Γ"][1], 0)
        self.assertAlmostEqual(kpoints["Γ"][2], 0)

        self.assertAlmostEqual(kpoints["e"][0], 0.0)
        self.assertAlmostEqual(kpoints["e"][1], 0.49999999999999956)
        self.assertAlmostEqual(kpoints["e"][2], 0.5000000000000002)

        d = False
        if np.allclose(kpoints["d_{1}"], [0.2530864197530836, 0.25308641975308915, 0.0], atol=1e-5) or np.allclose(
            kpoints["d"], [0.2530864197530836, 0.25308641975308915, 0.0], atol=1e-5
        ):
            d = True

        self.assertTrue(d)

        q = False
        if np.allclose(kpoints["q_{1}"], [0.2530864197530836, 0.25308641975308915, 0.5], atol=1e-5) or np.allclose(
            kpoints["q"], [0.2530864197530836, 0.25308641975308915, 0.5], atol=1e-5
        ):
            q = True

        self.assertTrue(q)

    def test_magnetic_kpath_generation(self):
        struct_file_path = os.path.join(test_dir_structs, "LaMnO3_magnetic.mcif")
        struct = Structure.from_file(struct_file_path)
        mga = CollinearMagneticStructureAnalyzer(struct)
        col_spin_orig = mga.get_structure_with_spin()
        col_spin_orig.add_spin_by_site([0.0] * 20)
        col_spin_sym = SpacegroupAnalyzer(col_spin_orig)
        col_spin_prim = col_spin_sym.get_primitive_standard_structure(international_monoclinic=False)

        magmom_vec_list = [np.zeros(3) for site in col_spin_prim]
        magmom_vec_list[4:8] = [
            np.array([3.87, 3.87, 0.0]),
            np.array([3.87, 3.87, 0.0]),
            np.array([-3.87, -3.87, 0.0]),
            np.array([-3.87, -3.87, 0.0]),
        ]
        col_spin_prim.add_site_property("magmom", magmom_vec_list)

        kpath = KPathLatimerMunro(col_spin_prim, has_magmoms=True)

        kpoints = kpath._kpath["kpoints"]
        labels = list(kpoints.keys())

        self.assertEqual(
            sorted(labels),
            sorted(["a", "b", "c", "d", "d_{1}", "e", "f", "g", "g_{1}", "Γ"]),
        )

        self.assertAlmostEqual(kpoints["e"][0], -0.4999999999999998)
        self.assertAlmostEqual(kpoints["e"][1], 0.0)
        self.assertAlmostEqual(kpoints["e"][2], 0.5000000000000002)

        self.assertAlmostEqual(kpoints["g"][0], -0.4999999999999999)
        self.assertAlmostEqual(kpoints["g"][1], -0.49999999999999994)
        self.assertAlmostEqual(kpoints["g"][2], 0.5000000000000002)

        self.assertAlmostEqual(kpoints["a"][0], -0.4999999999999999)
        self.assertAlmostEqual(kpoints["a"][1], 0.0)
        self.assertAlmostEqual(kpoints["a"][2], 0.0)

        self.assertAlmostEqual(kpoints["g_{1}"][0], 0.4999999999999999)
        self.assertAlmostEqual(kpoints["g_{1}"][1], -0.5)
        self.assertAlmostEqual(kpoints["g_{1}"][2], 0.5000000000000001)

        self.assertAlmostEqual(kpoints["f"][0], 0.0)
        self.assertAlmostEqual(kpoints["f"][1], -0.5)
        self.assertAlmostEqual(kpoints["f"][2], 0.5000000000000002)

        self.assertAlmostEqual(kpoints["c"][0], 0.0)
        self.assertAlmostEqual(kpoints["c"][1], 0.0)
        self.assertAlmostEqual(kpoints["c"][2], 0.5000000000000001)

        self.assertAlmostEqual(kpoints["b"][0], 0.0)
        self.assertAlmostEqual(kpoints["b"][1], -0.5)
        self.assertAlmostEqual(kpoints["b"][2], 0.0)

        self.assertAlmostEqual(kpoints["Γ"][0], 0)
        self.assertAlmostEqual(kpoints["Γ"][1], 0)
        self.assertAlmostEqual(kpoints["Γ"][2], 0)

        d = False
        if np.allclose(kpoints["d_{1}"], [-0.5, -0.5, 0.0], atol=1e-5) or np.allclose(
            kpoints["d"], [-0.5, -0.5, 0.0], atol=1e-5
        ):
            d = True

        self.assertTrue(d)

        g = False
        if np.allclose(kpoints["g_{1}"], [-0.5, -0.5, 0.5], atol=1e-5) or np.allclose(
            kpoints["g"], [-0.5, -0.5, 0.5], atol=1e-5
        ):
            g = True

        self.assertTrue(g)


if __name__ == "__main__":
    unittest.main()
