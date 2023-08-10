from __future__ import annotations

import unittest

import numpy as np
import pytest
from pytest import approx

from pymatgen.core.periodic_table import Species
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Molecule, Structure
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.symmetry.analyzer import PointGroupAnalyzer, SpacegroupAnalyzer, cluster_sites, iterative_symmetrize
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

test_dir_mol = f"{TEST_FILES_DIR}/molecules"


class TestSpacegroupAnalyzer(PymatgenTest):
    def setUp(self):
        p = Poscar.from_file(f"{TEST_FILES_DIR}/POSCAR")
        self.structure = p.structure
        self.sg = SpacegroupAnalyzer(self.structure, 0.001)
        self.disordered_structure = self.get_structure("Li10GeP2S12")
        self.disordered_sg = SpacegroupAnalyzer(self.disordered_structure, 0.001)
        struct = p.structure.copy()
        site = struct[0]
        del struct[0]
        struct.append(site.species, site.frac_coords)
        self.sg3 = SpacegroupAnalyzer(struct, 0.001)
        graphite = self.get_structure("Graphite")
        graphite.add_site_property("magmom", [0.1] * len(graphite))
        self.sg4 = SpacegroupAnalyzer(graphite, 0.001)
        self.structure4 = graphite

    def test_primitive(self):
        struct = Structure.from_spacegroup("Fm-3m", np.eye(3) * 3, ["Cu"], [[0, 0, 0]])
        spga = SpacegroupAnalyzer(struct)
        assert len(struct) == 4
        assert len(spga.find_primitive()) == 1

    def test_is_laue(self):
        struct = Structure.from_spacegroup("Fm-3m", np.eye(3) * 3, ["Cu"], [[0, 0, 0]])
        assert SpacegroupAnalyzer(struct).is_laue()

        assert self.sg.is_laue()

        assert self.disordered_sg.is_laue()

    def test_magnetic(self):
        lfp = PymatgenTest.get_structure("LiFePO4")
        sg = SpacegroupAnalyzer(lfp, 0.1)
        assert sg.get_space_group_symbol() == "Pnma"
        magmoms = [0] * len(lfp)
        magmoms[4] = 1
        magmoms[5] = -1
        magmoms[6] = 1
        magmoms[7] = -1
        lfp.add_site_property("magmom", magmoms)
        sg = SpacegroupAnalyzer(lfp, 0.1)
        assert sg.get_space_group_symbol() == "Pnma"

    def test_get_space_symbol(self):
        assert self.sg.get_space_group_symbol() == "Pnma"
        assert self.disordered_sg.get_space_group_symbol() == "P4_2/nmc"
        assert self.sg3.get_space_group_symbol() == "Pnma"
        assert self.sg4.get_space_group_symbol() == "P6_3/mmc"

    def test_get_space_number(self):
        assert self.sg.get_space_group_number() == 62
        assert self.disordered_sg.get_space_group_number() == 137
        assert self.sg4.get_space_group_number() == 194

    def test_get_hall(self):
        assert self.sg.get_hall() == "-P 2ac 2n"
        assert self.disordered_sg.get_hall() == "P 4n 2n -1n"

    def test_get_pointgroup(self):
        assert self.sg.get_point_group_symbol() == "mmm"
        assert self.disordered_sg.get_point_group_symbol() == "4/mmm"

    def test_get_point_group_operations(self):
        for sg, structure in [(self.sg, self.structure), (self.sg4, self.structure4)]:
            pg_ops = sg.get_point_group_operations()
            frac_symm_ops = sg.get_symmetry_operations()
            symm_ops = sg.get_symmetry_operations(True)
            for fop, op, pgop in zip(frac_symm_ops, symm_ops, pg_ops):
                # translation vector values should all be 0 or 0.5
                t = fop.translation_vector * 2
                assert np.allclose(t - np.round(t), 0)

                assert np.allclose(fop.rotation_matrix, pgop.rotation_matrix)
                for site in structure:
                    new_frac = fop.operate(site.frac_coords)
                    new_cart = op.operate(site.coords)
                    assert np.allclose(structure.lattice.get_fractional_coords(new_cart), new_frac)
                    found = False
                    new_site = PeriodicSite(site.species, new_cart, structure.lattice, coords_are_cartesian=True)
                    for test_site in structure:
                        if new_site.is_periodic_image(test_site, 1e-3):
                            found = True
                            break
                    assert found

                # Make sure this works for any position, not just the atomic
                # ones.
                random_fcoord = np.random.uniform(size=(3))
                random_ccoord = structure.lattice.get_cartesian_coords(random_fcoord)
                new_frac = fop.operate(random_fcoord)
                new_cart = op.operate(random_ccoord)
                assert np.allclose(structure.lattice.get_fractional_coords(new_cart), new_frac)

    def test_get_point_group_operations_uniq(self):
        # https://github.com/materialsproject/pymatgen/pull/2942
        struct = Structure.from_spacegroup(223, np.eye(3), ["Ni"], [[0.0, 0.0, 0.0]])
        pg_ops = SpacegroupAnalyzer(struct).get_point_group_operations(cartesian=True)
        assert len(pg_ops) == 48

    def test_get_symmetry_dataset(self):
        ds = self.sg.get_symmetry_dataset()
        assert ds["international"] == "Pnma"

    def test_init_cell(self):
        # see https://github.com/materialsproject/pymatgen/pull/3179
        li2o = Structure.from_file(f"{TEST_FILES_DIR}/Li2O.cif")

        # test that magmoms are not included in spacegroup analyzer when species have no spin
        # or no magmom site properties are set
        assert all(species.spin is None for species in li2o.types_of_species)
        assert "magmom" not in li2o.site_properties
        sga = SpacegroupAnalyzer(li2o)
        assert len(sga._cell) == 3  # no magmoms should be added!

        # give a magmom to random Li site
        li2o[0].properties["magmom"] = 0
        sga = SpacegroupAnalyzer(li2o)
        assert len(sga._cell) == 4  # magmoms should be added!
        assert sga._cell[3] == 12 * (0,)

        # now set spin for O only
        li2o = Structure.from_file(f"{TEST_FILES_DIR}/Li2O.cif")
        li2o.replace_species({"O2-": Species("O", oxidation_state=-2, spin=1)})
        assert not all(species.spin is None for species in li2o.types_of_species)
        sga = SpacegroupAnalyzer(li2o)
        assert len(sga._cell) == 4  # magmoms should be added!
        assert sga._cell[3] == tuple(
            8
            * [
                0,
            ]
            + 4
            * [
                1,
            ]
        )

    def test_get_symmetry(self):
        # see discussion in https://github.com/materialsproject/pymatgen/pull/2724
        Co8 = Structure.from_file(f"{TEST_FILES_DIR}/Co8.cif")
        symprec = 1e-1

        sga = SpacegroupAnalyzer(Co8, symprec=symprec)
        magmoms = [0] * len(Co8)  # bad magmoms, see https://github.com/materialsproject/pymatgen/pull/2727
        sga._cell = (*sga._cell, magmoms)
        with pytest.raises(
            ValueError,
            match=f"Symmetry detection failed for structure with formula {Co8.formula}. "
            f"Try setting {symprec=} to a different value.",
        ):
            sga._get_symmetry()

    def test_get_crystal_system(self):
        crystal_system = self.sg.get_crystal_system()
        assert crystal_system == "orthorhombic"
        assert self.disordered_sg.get_crystal_system() == "tetragonal"

        orig_spg = self.sg._space_group_data["number"]
        self.sg._space_group_data["number"] = 0
        with pytest.raises(ValueError, match="Received invalid space group 0"):
            self.sg.get_crystal_system()

        self.sg._space_group_data["number"] = orig_spg

    def test_get_refined_structure(self):
        for a in self.sg.get_refined_structure().lattice.angles:
            assert a == 90
        refined = self.disordered_sg.get_refined_structure()
        for a in refined.lattice.angles:
            assert a == 90
        assert refined.lattice.a == refined.lattice.b

        structure = self.get_structure("Li2O")
        structure.add_site_property("magmom", [1.0] * len(structure))
        sg = SpacegroupAnalyzer(structure, 0.01)
        refined_struct = sg.get_refined_structure(keep_site_properties=True)
        assert refined_struct.site_properties["magmom"] == [1.0] * len(refined_struct)

        structure = self.get_structure("Li2O")
        structure.add_site_property("magmom", [1.0] * len(structure))
        sg = SpacegroupAnalyzer(structure, 0.01)
        refined_struct = sg.get_refined_structure(keep_site_properties=False)
        assert refined_struct.site_properties.get("magmom") is None

    def test_symmetrized_structure(self):
        symm_struct = self.sg.get_symmetrized_structure()
        assert symm_struct.lattice.angles == (90, 90, 90)
        assert len(symm_struct.equivalent_sites) == 5

        symm_struct = self.disordered_sg.get_symmetrized_structure()
        assert len(symm_struct.equivalent_sites) == 8
        assert [len(i) for i in symm_struct.equivalent_sites] == [16, 4, 8, 4, 2, 8, 8, 8]
        s1 = symm_struct.equivalent_sites[1][1]
        s2 = symm_struct[symm_struct.equivalent_indices[1][1]]
        assert s1 == s2
        assert self.sg4.get_symmetrized_structure()[0].magmom == 0.1
        assert symm_struct.wyckoff_symbols[0] == "16h"

        # Check copying
        assert symm_struct.copy() == symm_struct
        d = symm_struct.as_dict()
        from pymatgen.symmetry.structure import SymmetrizedStructure

        ss = SymmetrizedStructure.from_dict(d)
        assert ss.wyckoff_symbols[0] == "16h"
        assert str(ss).startswith("SymmetrizedStructure\nFull Formula (Li20.2 Ge2.06 P3.94 S24)\nReduced Formula: ")

    def test_find_primitive(self):
        """F m -3 m Li2O testing of converting to primitive cell."""
        parser = CifParser(f"{TEST_FILES_DIR}/Li2O.cif")
        structure = parser.get_structures(False)[0]
        spga = SpacegroupAnalyzer(structure)
        primitive_structure = spga.find_primitive()
        assert primitive_structure.formula == "Li2 O1"
        assert primitive_structure.site_properties.get("magmom") is None
        # This isn't what is expected. All the angles should be 60
        assert primitive_structure.lattice.alpha == approx(60)
        assert primitive_structure.lattice.beta == approx(60)
        assert primitive_structure.lattice.gamma == approx(60)
        assert primitive_structure.volume == approx(structure.volume / 4.0)

        structure = parser.get_structures(False)[0]
        structure.add_site_property("magmom", [1.0] * len(structure))
        spga = SpacegroupAnalyzer(structure)
        primitive_structure = spga.find_primitive(keep_site_properties=True)
        assert primitive_structure.site_properties["magmom"] == [1.0] * len(primitive_structure)

        structure = parser.get_structures(False)[0]
        structure.add_site_property("magmom", [1.0] * len(structure))
        spga = SpacegroupAnalyzer(structure)
        primitive_structure = spga.find_primitive(keep_site_properties=False)
        assert primitive_structure.site_properties.get("magmom") is None

    def test_get_ir_reciprocal_mesh(self):
        grid = self.sg.get_ir_reciprocal_mesh()
        assert len(grid) == 216
        assert grid[1][0][0] == approx(0.1)
        assert grid[1][0][1] == approx(0.0)
        assert grid[1][0][2] == approx(0.0)
        assert grid[1][1] == approx(2)

    def test_get_ir_reciprocal_mesh_map(self):
        mesh = (6, 6, 6)
        grid = self.sg.get_ir_reciprocal_mesh(mesh=mesh)
        full_grid, mapping = self.sg.get_ir_reciprocal_mesh_map(mesh=mesh)
        assert len(np.unique(mapping)) == len(grid)
        for _, i in enumerate(np.unique(mapping)):
            assert full_grid[i][0] == approx(grid[_][0][0])
            assert full_grid[i][1] == approx(grid[_][0][1])
            assert full_grid[i][2] == approx(grid[_][0][2])

    def test_get_conventional_standard_structure(self):
        parser = CifParser(f"{TEST_FILES_DIR}/bcc_1927.cif")
        structure = parser.get_structures(False)[0]
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = spga.get_conventional_standard_structure()
        assert conv.lattice.angles == (90, 90, 90)
        assert conv.lattice.lengths == approx([9.1980270633769461] * 3)

        parser = CifParser(f"{TEST_FILES_DIR}/btet_1915.cif")
        structure = parser.get_structures(False)[0]
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = spga.get_conventional_standard_structure()
        assert conv.lattice.angles == (90, 90, 90)
        assert conv.lattice.a == approx(5.0615106678044235)
        assert conv.lattice.b == approx(5.0615106678044235)
        assert conv.lattice.c == approx(4.2327080177761687)

        parser = CifParser(f"{TEST_FILES_DIR}/orci_1010.cif")
        structure = parser.get_structures(False)[0]
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = spga.get_conventional_standard_structure()
        assert conv.lattice.angles == (90, 90, 90)
        assert conv.lattice.a == approx(2.9542233922299999)
        assert conv.lattice.b == approx(4.6330325651443296)
        assert conv.lattice.c == approx(5.373703587040775)

        parser = CifParser(f"{TEST_FILES_DIR}/orcc_1003.cif")
        structure = parser.get_structures(False)[0]
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = spga.get_conventional_standard_structure()
        assert conv.lattice.angles == (90, 90, 90)
        assert conv.lattice.a == approx(4.1430033493799998)
        assert conv.lattice.b == approx(31.437979757624728)
        assert conv.lattice.c == approx(3.99648651)

        parser = CifParser(f"{TEST_FILES_DIR}/orac_632475.cif")
        structure = parser.get_structures(False)[0]
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = spga.get_conventional_standard_structure()
        assert conv.lattice.angles == (90, 90, 90)
        assert conv.lattice.lengths == approx([3.1790663399999999, 9.9032878699999998, 3.5372412099999999])

        parser = CifParser(f"{TEST_FILES_DIR}/monoc_1028.cif")
        structure = parser.get_structures(False)[0]
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = spga.get_conventional_standard_structure()
        assert conv.lattice.angles == approx([90, 117.53832420192903, 90])
        assert conv.lattice.lengths == approx([14.033435583000625, 3.96052850731, 6.8743926325200002])
        parser = CifParser(f"{TEST_FILES_DIR}/hex_1170.cif")
        structure = parser.get_structures(False)[0]
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = spga.get_conventional_standard_structure()
        assert conv.lattice.angles == approx([90, 90, 120])
        assert conv.lattice.lengths == approx([3.699919902005897, 3.699919902005897, 6.9779585500000003])

        structure = Structure.from_file(f"{TEST_FILES_DIR}/tric_684654.json")
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = spga.get_conventional_standard_structure()
        assert conv.lattice.alpha == approx(74.09581916308757)
        assert conv.lattice.beta == approx(75.72817279281173)
        assert conv.lattice.gamma == approx(63.63234318667333)
        assert conv.lattice.a == approx(3.741372924048738)
        assert conv.lattice.b == approx(3.9883228679270686)
        assert conv.lattice.c == approx(7.288495840048958)

        structure = Structure.from_file(f"{TEST_FILES_DIR}/tric_684654.json")
        structure.add_site_property("magmom", [1.0] * len(structure))
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = spga.get_conventional_standard_structure(keep_site_properties=True)
        assert conv.site_properties["magmom"] == [1.0] * len(conv)

        structure = Structure.from_file(f"{TEST_FILES_DIR}/tric_684654.json")
        structure.add_site_property("magmom", [1.0] * len(structure))
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        conv = spga.get_conventional_standard_structure(keep_site_properties=False)
        assert conv.site_properties.get("magmom") is None

    def test_get_primitive_standard_structure(self):
        parser = CifParser(f"{TEST_FILES_DIR}/bcc_1927.cif")
        structure = parser.get_structures(False)[0]
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = spga.get_primitive_standard_structure()
        assert prim.lattice.alpha == approx(109.47122063400001)
        assert prim.lattice.beta == approx(109.47122063400001)
        assert prim.lattice.gamma == approx(109.47122063400001)
        assert prim.lattice.a == approx(7.9657251015812145)
        assert prim.lattice.b == approx(7.9657251015812145)
        assert prim.lattice.c == approx(7.9657251015812145)

        parser = CifParser(f"{TEST_FILES_DIR}/btet_1915.cif")
        structure = parser.get_structures(False)[0]
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = spga.get_primitive_standard_structure()
        assert prim.lattice.alpha == approx(105.015053349)
        assert prim.lattice.beta == approx(105.015053349)
        assert prim.lattice.gamma == approx(118.80658411899999)
        assert prim.lattice.a == approx(4.1579321075608791)
        assert prim.lattice.b == approx(4.1579321075608791)
        assert prim.lattice.c == approx(4.1579321075608791)

        parser = CifParser(f"{TEST_FILES_DIR}/orci_1010.cif")
        structure = parser.get_structures(False)[0]
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = spga.get_primitive_standard_structure()
        assert prim.lattice.alpha == approx(134.78923546600001)
        assert prim.lattice.beta == approx(105.856239333)
        assert prim.lattice.gamma == approx(91.276341676000001)
        assert prim.lattice.a == approx(3.8428217771014852)
        assert prim.lattice.b == approx(3.8428217771014852)
        assert prim.lattice.c == approx(3.8428217771014852)

        parser = CifParser(f"{TEST_FILES_DIR}/orcc_1003.cif")
        structure = parser.get_structures(False)[0]
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = spga.get_primitive_standard_structure()
        assert prim.lattice.alpha == approx(90)
        assert prim.lattice.beta == approx(90)
        assert prim.lattice.gamma == approx(164.985257335)
        assert prim.lattice.a == approx(15.854897098324196)
        assert prim.lattice.b == approx(15.854897098324196)
        assert prim.lattice.c == approx(3.99648651)

        parser = CifParser(f"{TEST_FILES_DIR}/orac_632475.cif")
        structure = parser.get_structures(False)[0]
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = spga.get_primitive_standard_structure()
        assert prim.lattice.alpha == approx(90)
        assert prim.lattice.beta == approx(90)
        assert prim.lattice.gamma == approx(144.40557588533386)
        assert prim.lattice.a == approx(5.2005185662155391)
        assert prim.lattice.b == approx(5.2005185662155391)
        assert prim.lattice.c == approx(3.5372412099999999)

        parser = CifParser(f"{TEST_FILES_DIR}/monoc_1028.cif")
        structure = parser.get_structures(False)[0]
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = spga.get_primitive_standard_structure()
        assert prim.lattice.alpha == approx(63.579155761999999)
        assert prim.lattice.beta == approx(116.42084423747779)
        assert prim.lattice.gamma == approx(148.47965136208569)
        assert prim.lattice.a == approx(7.2908007159612325)
        assert prim.lattice.b == approx(7.2908007159612325)
        assert prim.lattice.c == approx(6.8743926325200002)

        parser = CifParser(f"{TEST_FILES_DIR}/hex_1170.cif")
        structure = parser.get_structures(False)[0]
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = spga.get_primitive_standard_structure()
        assert prim.lattice.alpha == approx(90)
        assert prim.lattice.beta == approx(90)
        assert prim.lattice.gamma == approx(120)
        assert prim.lattice.a == approx(3.699919902005897)
        assert prim.lattice.b == approx(3.699919902005897)
        assert prim.lattice.c == approx(6.9779585500000003)

        parser = CifParser(f"{TEST_FILES_DIR}/rhomb_3478_conv.cif")
        structure = parser.get_structures(False)[0]
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = spga.get_primitive_standard_structure()
        assert prim.lattice.alpha == approx(28.049186140546812)
        assert prim.lattice.beta == approx(28.049186140546812)
        assert prim.lattice.gamma == approx(28.049186140546812)
        assert prim.lattice.a == approx(5.9352627428399982)
        assert prim.lattice.b == approx(5.9352627428399982)
        assert prim.lattice.c == approx(5.9352627428399982)

        parser = CifParser(f"{TEST_FILES_DIR}/rhomb_3478_conv.cif")
        structure = parser.get_structures(False)[0]
        structure.add_site_property("magmom", [1.0] * len(structure))
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = spga.get_primitive_standard_structure(keep_site_properties=True)
        assert prim.site_properties["magmom"] == [1.0] * len(prim)

        parser = CifParser(f"{TEST_FILES_DIR}/rhomb_3478_conv.cif")
        structure = parser.get_structures(False)[0]
        structure.add_site_property("magmom", [1.0] * len(structure))
        spga = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim = spga.get_primitive_standard_structure(keep_site_properties=False)
        assert prim.site_properties.get("magmom") is None

    def test_tricky_structure(self):
        # for some reason this structure kills spglib1.9
        # 1.7 can't find symmetry either, but at least doesn't kill python
        struct = Structure.from_file(f"{TEST_FILES_DIR}/POSCAR.tricky_symmetry")
        sa = SpacegroupAnalyzer(struct, 0.1)
        assert sa.get_space_group_symbol() == "I4/mmm"
        assert sa.get_space_group_number() == 139
        assert sa.get_point_group_symbol() == "4/mmm"
        assert sa.get_crystal_system() == "tetragonal"
        assert sa.get_hall() == "-I 4 2"


class TestSpacegroup(unittest.TestCase):
    def setUp(self):
        p = Poscar.from_file(f"{TEST_FILES_DIR}/POSCAR")
        self.structure = p.structure
        self.sg1 = SpacegroupAnalyzer(self.structure, 0.001).get_space_group_operations()

    def test_are_symmetrically_equivalent(self):
        sites1 = [self.structure[i] for i in [0, 1]]
        sites2 = [self.structure[i] for i in [2, 3]]
        assert self.sg1.are_symmetrically_equivalent(sites1, sites2, 1e-3)

        sites1 = [self.structure[i] for i in [0, 1]]
        sites2 = [self.structure[i] for i in [0, 2]]
        assert not self.sg1.are_symmetrically_equivalent(sites1, sites2, 1e-3)


H2O2 = Molecule(
    ["O", "O", "H", "H"],
    [
        [0, 0.727403, -0.050147],
        [0, -0.727403, -0.050147],
        [0.83459, 0.897642, 0.401175],
        [-0.83459, -0.897642, 0.401175],
    ],
)

C2H2F2Br2 = Molecule(
    ["C", "C", "F", "Br", "H", "F", "H", "Br"],
    [
        [-0.752000, 0.001000, -0.141000],
        [0.752000, -0.001000, 0.141000],
        [-1.158000, 0.991000, 0.070000],
        [-1.240000, -0.737000, 0.496000],
        [-0.924000, -0.249000, -1.188000],
        [1.158000, -0.991000, -0.070000],
        [0.924000, 0.249000, 1.188000],
        [1.240000, 0.737000, -0.496000],
    ],
)

H2O = Molecule(
    ["H", "O", "H"],
    [[0, 0.780362, -0.456316], [0, 0, 0.114079], [0, -0.780362, -0.456316]],
)

C2H4 = Molecule(
    ["C", "C", "H", "H", "H", "H"],
    [
        [0.0000, 0.0000, 0.6695],
        [0.0000, 0.0000, -0.6695],
        [0.0000, 0.9289, 1.2321],
        [0.0000, -0.9289, 1.2321],
        [0.0000, 0.9289, -1.2321],
        [0.0000, -0.9289, -1.2321],
    ],
)

NH3 = Molecule(
    ["N", "H", "H", "H"],
    [
        [0.0000, 0.0000, 0.0000],
        [0.0000, -0.9377, -0.3816],
        [0.8121, 0.4689, -0.3816],
        [-0.8121, 0.4689, -0.3816],
    ],
)

BF3 = Molecule(
    ["B", "F", "F", "F"],
    [
        [0.0000, 0.0000, 0.0000],
        [0.0000, -0.9377, 0.00],
        [0.8121, 0.4689, 0],
        [-0.8121, 0.4689, 0],
    ],
)

CH4 = Molecule(
    ["C", "H", "H", "H", "H"],
    [
        [0.000000, 0.000000, 0.000000],
        [0.000000, 0.000000, 1.08],
        [1.026719, 0.000000, -0.363000],
        [-0.513360, -0.889165, -0.363000],
        [-0.513360, 0.889165, -0.363000],
    ],
)

PF6 = Molecule(
    ["P", "F", "F", "F", "F", "F", "F"],
    [[0, 0, 0], [0, 0, 1], [0, 0, -1], [0, 1, 0], [0, -1, 0], [1, 0, 0], [-1, 0, 0]],
)


class TestPointGroupAnalyzer(PymatgenTest):
    def test_spherical(self):
        a = PointGroupAnalyzer(CH4)
        assert a.sch_symbol == "Td"
        assert len(a.get_pointgroup()) == 24
        assert a.get_rotational_symmetry_number() == 12
        a = PointGroupAnalyzer(H2O)
        assert a.get_rotational_symmetry_number() == 2
        a = PointGroupAnalyzer(PF6)
        assert a.sch_symbol == "Oh"
        assert len(a.get_pointgroup()) == 48
        m = Molecule.from_file(f"{test_dir_mol}/c60.xyz")
        a = PointGroupAnalyzer(m)
        assert a.sch_symbol == "Ih"

        cube_species = ["C", "C", "C", "C", "C", "C", "C", "C"]
        cube_coords = [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [1, 1, 0],
            [0, 0, 1],
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 1],
        ]

        m = Molecule(cube_species, cube_coords)
        a = PointGroupAnalyzer(m, 0.1)
        assert a.sch_symbol == "Oh"

    def test_tricky(self):
        m = Molecule.from_file(f"{test_dir_mol}/dh.xyz")
        a = PointGroupAnalyzer(m, 0.1)
        assert a.sch_symbol == "D*h"

    def test_linear(self):
        coords = [
            [0.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, 1.08],
            [0, 0.000000, -1.08],
        ]
        mol = Molecule(["C", "H", "H"], coords)
        a = PointGroupAnalyzer(mol)
        assert a.sch_symbol == "D*h"
        mol = Molecule(["C", "H", "N"], coords)
        a = PointGroupAnalyzer(mol)
        assert a.sch_symbol == "C*v"

    def test_asym_top(self):
        coords = [
            [0.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, 1.08],
            [1.026719, 0.000000, -0.363000],
            [-0.513360, -0.889165, -0.363000],
            [-0.513360, 0.889165, -0.363000],
        ]
        mol = Molecule(["C", "H", "F", "Br", "Cl"], coords)
        a = PointGroupAnalyzer(mol)

        assert a.sch_symbol == "C1"
        assert len(a.get_pointgroup()) == 1
        coords = [
            [0.000000, 0.000000, 1.08],
            [1.026719, 0.000000, -0.363000],
            [-0.513360, -0.889165, -0.363000],
            [-0.513360, 0.889165, -0.363000],
        ]
        cs_mol = Molecule(["H", "F", "Cl", "Cl"], coords)
        a = PointGroupAnalyzer(cs_mol)
        assert a.sch_symbol == "Cs"
        assert len(a.get_pointgroup()) == 2
        a = PointGroupAnalyzer(C2H2F2Br2)
        assert a.sch_symbol == "Ci"
        assert len(a.get_pointgroup()) == 2

    def test_cyclic(self):
        a = PointGroupAnalyzer(H2O2)
        assert a.sch_symbol == "C2"
        assert len(a.get_pointgroup()) == 2
        a = PointGroupAnalyzer(H2O)
        assert a.sch_symbol == "C2v"
        assert len(a.get_pointgroup()) == 4
        a = PointGroupAnalyzer(NH3)
        assert a.sch_symbol == "C3v"
        assert len(a.get_pointgroup()) == 6
        cs2 = Molecule.from_file(f"{test_dir_mol}/Carbon_Disulfide.xyz")
        a = PointGroupAnalyzer(cs2, eigen_tolerance=0.001)
        assert a.sch_symbol == "C2v"

    def test_dihedral(self):
        a = PointGroupAnalyzer(C2H4)
        assert a.sch_symbol == "D2h"
        assert len(a.get_pointgroup()) == 8
        a = PointGroupAnalyzer(BF3)
        assert a.sch_symbol == "D3h"
        assert len(a.get_pointgroup()) == 12
        m = Molecule.from_file(f"{test_dir_mol}/b12h12.xyz")
        a = PointGroupAnalyzer(m)
        assert a.sch_symbol == "Ih"

    def test_symmetrize_molecule1(self):
        np.random.seed(77)
        distortion = np.random.randn(len(C2H4), 3) / 10
        dist_mol = Molecule(C2H4.species, C2H4.cart_coords + distortion)

        eq = iterative_symmetrize(dist_mol, max_n=100, epsilon=1e-7)
        sym_mol, eq_sets, ops = eq["sym_mol"], eq["eq_sets"], eq["sym_ops"]

        assert {0, 1} in eq_sets.values()
        assert {2, 3, 4, 5} in eq_sets.values()

        coords = sym_mol.cart_coords
        for i, eq_set in eq_sets.items():
            for j in eq_set:
                _ = np.dot(ops[i][j], coords[i])
                assert np.allclose(np.dot(ops[i][j], coords[i]), coords[j])

    def test_symmetrize_molecule2(self):
        np.random.seed(77)
        distortion = np.random.randn(len(C2H2F2Br2), 3) / 20
        dist_mol = Molecule(C2H2F2Br2.species, C2H2F2Br2.cart_coords + distortion)
        pa1 = PointGroupAnalyzer(C2H2F2Br2, tolerance=0.1)
        assert pa1.get_pointgroup().sch_symbol == "Ci"
        pa2 = PointGroupAnalyzer(dist_mol, tolerance=0.1)
        assert pa2.get_pointgroup().sch_symbol == "C1"
        eq = iterative_symmetrize(dist_mol, tolerance=0.3)
        pa3 = PointGroupAnalyzer(eq["sym_mol"], tolerance=0.1)
        assert pa3.get_pointgroup().sch_symbol == "Ci"

    def test_get_kpoint_weights(self):
        for name in ["SrTiO3", "LiFePO4", "Graphite"]:
            struct = PymatgenTest.get_structure(name)
            spga = SpacegroupAnalyzer(struct)
            ir_mesh = spga.get_ir_reciprocal_mesh((4, 4, 4))
            weights = [i[1] for i in ir_mesh]
            weights = np.array(weights) / sum(weights)
            for i, w in zip(weights, spga.get_kpoint_weights([i[0] for i in ir_mesh])):
                assert i == approx(w)

        for name in ["SrTiO3", "LiFePO4", "Graphite"]:
            struct = PymatgenTest.get_structure(name)
            spga = SpacegroupAnalyzer(struct)
            ir_mesh = spga.get_ir_reciprocal_mesh((1, 2, 3))
            weights = [i[1] for i in ir_mesh]
            weights = np.array(weights) / sum(weights)
            for i, w in zip(weights, spga.get_kpoint_weights([i[0] for i in ir_mesh])):
                assert i == approx(w)

        vasp_run = Vasprun(f"{TEST_FILES_DIR}/vasprun.xml")
        spga = SpacegroupAnalyzer(vasp_run.final_structure)
        wts = spga.get_kpoint_weights(vasp_run.actual_kpoints)

        for w1, w2 in zip(vasp_run.actual_kpoints_weights, wts):
            assert w1 == approx(w2)

        kpts = [[0, 0, 0], [0.15, 0.15, 0.15], [0.2, 0.2, 0.2]]
        with pytest.raises(
            ValueError, match="Unable to find 1:1 corresponding between input kpoints and irreducible grid"
        ):
            spga.get_kpoint_weights(kpts)


class TestFunc(unittest.TestCase):
    def test_cluster_sites(self):
        o, c = cluster_sites(CH4, 0.1)
        assert o.specie.symbol == "C"
        assert len(c) == 1
        o, c = cluster_sites(C2H2F2Br2.get_centered_molecule(), 0.1)
        assert o is None
        assert len(c) == 4
