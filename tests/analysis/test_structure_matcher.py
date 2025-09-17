from __future__ import annotations

import json

import numpy as np
import pytest
from monty.json import MontyDecoder
from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.analysis.structure_matcher import (
    ElementComparator,
    FrameworkComparator,
    OccupancyComparator,
    OrderDisorderElementComparator,
    StructureMatcher,
    get_linear_assignment_solution,
)
from pymatgen.core import Element, Lattice, Structure, SymmOp
from pymatgen.util.coord import find_in_coord_list_pbc
from pymatgen.util.testing import TEST_FILES_DIR, VASP_IN_DIR, MatSciTest

TEST_DIR = f"{TEST_FILES_DIR}/analysis/structure_matcher"


class TestLinearAssignment:
    def test_square_cost_matrix(self):
        costs_0 = np.array(
            [
                [19, 95, 9, 43, 62, 90, 10, 77, 71, 27],
                [26, 30, 88, 78, 87, 2, 14, 71, 78, 11],
                [48, 70, 26, 82, 32, 16, 36, 26, 42, 79],
                [47, 46, 93, 66, 38, 20, 73, 39, 55, 51],
                [1, 81, 31, 49, 20, 24, 95, 80, 82, 11],
                [81, 48, 35, 54, 35, 55, 27, 87, 96, 7],
                [42, 17, 60, 73, 37, 36, 79, 3, 60, 82],
                [14, 57, 23, 69, 93, 78, 56, 49, 83, 36],
                [11, 37, 24, 70, 62, 35, 64, 18, 99, 20],
                [73, 11, 98, 50, 19, 96, 61, 73, 98, 14],
            ]
        )
        sol, min_cost = get_linear_assignment_solution(costs_0)
        assert min_cost == 194, "Incorrect cost"

        costs_1 = np.array(
            [
                [95, 60, 89, 38, 36, 38, 58, 94, 66, 23],
                [37, 0, 40, 58, 97, 85, 18, 54, 86, 21],
                [9, 74, 11, 45, 65, 64, 27, 88, 24, 26],
                [58, 90, 6, 36, 17, 21, 2, 12, 80, 90],
                [33, 0, 74, 75, 11, 84, 34, 7, 39, 0],
                [17, 61, 94, 68, 27, 41, 33, 86, 59, 2],
                [61, 94, 36, 53, 66, 33, 15, 87, 97, 11],
                [22, 20, 57, 69, 15, 9, 15, 8, 82, 68],
                [40, 0, 13, 61, 67, 40, 29, 25, 72, 44],
                [13, 97, 97, 54, 5, 30, 44, 75, 16, 0],
            ]
        )
        sol, min_cost = get_linear_assignment_solution(costs_1)
        assert min_cost == 125, "Incorrect cost"

        costs_2 = np.array(
            [
                [34, 44, 72, 13, 10, 58, 16, 1, 10, 61],
                [54, 70, 99, 4, 64, 0, 15, 94, 39, 46],
                [49, 21, 80, 68, 96, 58, 24, 87, 79, 67],
                [86, 46, 58, 83, 83, 56, 83, 65, 4, 96],
                [48, 95, 64, 34, 75, 82, 64, 47, 35, 19],
                [11, 49, 6, 57, 80, 26, 47, 63, 75, 75],
                [74, 7, 15, 83, 64, 26, 78, 17, 67, 46],
                [19, 13, 2, 26, 52, 16, 65, 24, 2, 98],
                [36, 7, 93, 93, 11, 39, 94, 26, 46, 69],
                [32, 95, 37, 50, 97, 96, 12, 70, 40, 93],
            ]
        )
        sol, min_cost = get_linear_assignment_solution(costs_2)
        assert min_cost == 110, "Incorrect cost"

    def test_rectangular_cost_matrix(self):
        costs_0 = np.array(
            [
                [19, 95, 9, 43, 62, 90, 10, 77, 71, 27],
                [26, 30, 88, 78, 87, 2, 14, 71, 78, 11],
                [48, 70, 26, 82, 32, 16, 36, 26, 42, 79],
                [47, 46, 93, 66, 38, 20, 73, 39, 55, 51],
                [1, 81, 31, 49, 20, 24, 95, 80, 82, 11],
                [81, 48, 35, 54, 35, 55, 27, 87, 96, 7],
                [42, 17, 60, 73, 37, 36, 79, 3, 60, 82],
                [14, 57, 23, 69, 93, 78, 56, 49, 83, 36],
                [11, 37, 24, 70, 62, 35, 64, 18, 99, 20],
            ]
        )
        sol, min_cost_0 = get_linear_assignment_solution(costs_0)

        costs_1 = np.array(
            [
                [19, 95, 9, 43, 62, 90, 10, 77, 71, 27],
                [26, 30, 88, 78, 87, 2, 14, 71, 78, 11],
                [48, 70, 26, 82, 32, 16, 36, 26, 42, 79],
                [47, 46, 93, 66, 38, 20, 73, 39, 55, 51],
                [1, 81, 31, 49, 20, 24, 95, 80, 82, 11],
                [81, 48, 35, 54, 35, 55, 27, 87, 96, 7],
                [42, 17, 60, 73, 37, 36, 79, 3, 60, 82],
                [14, 57, 23, 69, 93, 78, 56, 49, 83, 36],
                [11, 37, 24, 70, 62, 35, 64, 18, 99, 20],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            ]
        )
        sol, min_cost_1 = get_linear_assignment_solution(costs_1)
        assert list(sol) == [6, 5, 8, 4, 0, 9, 1, 2, 7, 3]

        assert min_cost_0 == min_cost_1

        with pytest.raises(ValueError, match="cost matrix must have at least as many columns as rows"):
            _ = get_linear_assignment_solution(costs_0.T)

    def test_another_case(self):
        costs = np.array(
            [
                [
                    0.03900238875468465,
                    0.003202415721817453,
                    0.20107156847937024,
                    0.0,
                    0.5002116398420846,
                    0.11951326861160616,
                    0.0,
                    0.5469032363997579,
                    0.3243791041219123,
                    0.1119882291981289,
                ],
                [
                    0.6048342640688928,
                    0.3847629088356139,
                    0.0,
                    0.44358269535118944,
                    0.45925670625165016,
                    0.31416882324798145,
                    0.8065128182180494,
                    0.0,
                    0.26153475286065075,
                    0.6862799559241944,
                ],
                [
                    0.5597215814025246,
                    0.15133664165478322,
                    0.0,
                    0.6218101659263295,
                    0.15438455134183793,
                    0.17281467064043232,
                    0.8458127968475472,
                    0.020860721537078075,
                    0.1926886361228456,
                    0.0,
                ],
                [
                    0.0,
                    0.0,
                    0.6351848838666995,
                    0.21261247074659906,
                    0.4811603832432241,
                    0.6663733668270337,
                    0.63970145187428,
                    0.1415815172623256,
                    0.5294574133825874,
                    0.5576702829768786,
                ],
                [
                    0.25052904388309016,
                    0.2309392544588127,
                    0.0656162006684271,
                    0.0248922362001176,
                    0.0,
                    0.2101808638720748,
                    0.6529031699724193,
                    0.1503003886507902,
                    0.375576165698992,
                    0.7368328849560374,
                ],
                [
                    0.0,
                    0.042215873587668984,
                    0.10326920761908365,
                    0.3562551151517992,
                    0.9170343984958856,
                    0.818783531026254,
                    0.7896770426052844,
                    0.0,
                    0.6573135097946438,
                    0.17806189728574429,
                ],
                [
                    0.44992199118890386,
                    0.0,
                    0.38548898339412585,
                    0.6269193883601244,
                    1.0022861602564634,
                    0.0,
                    0.1869765500803764,
                    0.03474156273982543,
                    0.3715310534696664,
                    0.6197122486230232,
                ],
                [
                    0.37939853696836545,
                    0.2421427374018027,
                    0.5586150342727723,
                    0.0,
                    0.7171485794073893,
                    0.8021029235865014,
                    0.11213464903613135,
                    0.6497896761660467,
                    0.3274108706187846,
                    0.0,
                ],
                [
                    0.6674685746225324,
                    0.5347953626128863,
                    0.11461835366075113,
                    0.0,
                    0.8170639855163434,
                    0.7291931505979982,
                    0.3149153087053108,
                    0.1008681103294512,
                    0.0,
                    0.18751172321112997,
                ],
                [
                    0.6985944652913342,
                    0.6139921045056471,
                    0.0,
                    0.4393266955771965,
                    0.0,
                    0.47265399761400695,
                    0.3674241844351025,
                    0.04731761392352629,
                    0.21484886069716147,
                    0.16488710920126137,
                ],
            ]
        )
        sol, min_cost = get_linear_assignment_solution(costs)
        assert min_cost == approx(0)

    def test_small_range(self):
        # can be tricky for the augment step
        costs = np.array(
            [
                [4, 5, 5, 6, 8, 4, 7, 4, 7, 8],
                [5, 6, 6, 6, 7, 6, 6, 5, 6, 7],
                [4, 4, 5, 7, 7, 4, 8, 4, 7, 7],
                [6, 7, 6, 6, 7, 6, 6, 6, 6, 6],
                [4, 4, 4, 6, 6, 4, 7, 4, 7, 7],
                [4, 5, 5, 6, 8, 4, 7, 4, 7, 8],
                [5, 7, 5, 5, 5, 6, 4, 5, 4, 6],
                [8, 9, 8, 4, 5, 9, 4, 8, 4, 4],
                [5, 6, 6, 6, 7, 6, 6, 5, 6, 7],
                [5, 6, 6, 6, 7, 6, 6, 5, 6, 7],
            ]
        )
        assert get_linear_assignment_solution(costs)[1] == 48

    def test_boolean_inputs(self):
        costs_ones = np.ones((135, 135), dtype=bool)
        np.fill_diagonal(costs_ones, val=False)


class TestStructureMatcher(MatSciTest):
    def setup_method(self):
        with open(f"{TEST_FILES_DIR}/entries/TiO2_entries.json", encoding="utf-8") as file:
            entries = json.load(file, cls=MontyDecoder)
        self.struct_list = [ent.structure for ent in entries]
        self.oxi_structs = [
            self.get_structure("Li2O"),
            Structure.from_file(f"{VASP_IN_DIR}/POSCAR_Li2O"),
        ]

    def test_ignore_species(self):
        s1 = Structure.from_file(f"{TEST_FILES_DIR}/cif/LiFePO4.cif")
        s2 = Structure.from_file(f"{VASP_IN_DIR}/POSCAR")
        matcher = StructureMatcher(ignored_species=["Li"], primitive_cell=False, attempt_supercell=True)
        assert matcher.fit(s1, s2)
        assert matcher.fit_anonymous(s1, s2)
        groups = matcher.group_structures([s1, s2])
        assert len(groups) == 1
        s2.make_supercell((2, 1, 1))
        ss1 = matcher.get_s2_like_s1(s2, s1, include_ignored_species=True)
        assert ss1.lattice.a == approx(20.820740000000001)
        assert ss1.reduced_formula == "LiFePO4"

        assert {
            k.symbol: v.symbol for k, v in matcher.get_best_electronegativity_anonymous_mapping(s1, s2).items()
        } == {
            "Fe": "Fe",
            "P": "P",
            "O": "O",
        }

    def test_get_supercell_size(self):
        latt1 = Lattice.cubic(1)
        latt2 = Lattice.cubic(0.9)
        s1 = Structure(latt1, ["Mg", "Cu", "Ag", "Cu", "Ag"], [[0] * 3] * 5)
        s2 = Structure(latt2, ["Cu", "Cu", "Ag"], [[0] * 3] * 3)

        sm = StructureMatcher(supercell_size="volume")
        assert sm._get_supercell_size(s1, s2) == (1, True)
        assert sm._get_supercell_size(s2, s1) == (1, True)

        sm = StructureMatcher(supercell_size="num_sites")
        assert sm._get_supercell_size(s1, s2) == (2, False)
        assert sm._get_supercell_size(s2, s1) == (2, True)

        sm = StructureMatcher(supercell_size="Ag")
        assert sm._get_supercell_size(s1, s2) == (2, False)
        assert sm._get_supercell_size(s2, s1) == (2, True)

        sm = StructureMatcher(supercell_size=["Ag", "Cu"])
        assert sm._get_supercell_size(s1, s2) == (1, True)
        assert sm._get_supercell_size(s2, s1) == (1, True)

        sm = StructureMatcher(supercell_size="invalid")
        with pytest.raises(ValueError, match="Can't parse Element or Species from 'invalid'"):
            sm._get_supercell_size(s1, s2)

    def test_cmp_fstruct(self):
        sm = StructureMatcher()

        s1 = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
        s2 = np.array([[0.11, 0.22, 0.33]])
        frac_tol = np.array([0.02, 0.03, 0.04])
        mask = np.array([[False, False]])
        mask2 = np.array([[True, False]])

        with pytest.raises(ValueError, match=r"len\(s1\)=1 must be larger than len\(s2\)=2"):
            sm._cmp_fstruct(s2, s1, frac_tol, mask.T)
        with pytest.raises(ValueError, match="mask has incorrect shape"):
            sm._cmp_fstruct(s1, s2, frac_tol, mask.T)

        assert sm._cmp_fstruct(s1, s2, frac_tol, mask)
        assert not sm._cmp_fstruct(s1, s2, frac_tol / 2, mask)
        assert not sm._cmp_fstruct(s1, s2, frac_tol, mask2)

    def test_cart_dists(self):
        sm = StructureMatcher()
        lattice = Lattice.orthorhombic(1, 2, 3)

        s1 = np.array([[0.13, 0.25, 0.37], [0.1, 0.2, 0.3]])
        s2 = np.array([[0.11, 0.22, 0.33]])
        s3 = np.array([[0.1, 0.2, 0.3], [0.11, 0.2, 0.3]])
        s4 = np.array([[0.1, 0.2, 0.3], [0.1, 0.6, 0.7]])
        mask = np.array([[False, False]])
        mask2 = np.array([[False, True]])
        mask3 = np.array([[False, False], [False, False]])
        mask4 = np.array([[False, True], [False, True]])

        n1 = (len(s1) / lattice.volume) ** (1 / 3)
        n2 = (len(s2) / lattice.volume) ** (1 / 3)

        with pytest.raises(ValueError, match=r"len\(s1\)=1 must be larger than len\(s2\)=2"):
            sm._cart_dists(s2, s1, lattice, mask.T, n2)
        with pytest.raises(ValueError, match="mask has incorrect shape"):
            sm._cart_dists(s1, s2, lattice, mask.T, n1)

        distances, trac_trans_vec, solution = sm._cart_dists(s1, s2, lattice, mask, n1)
        assert_allclose(distances, [0])
        assert_allclose(trac_trans_vec, [-0.01, -0.02, -0.03])
        assert_allclose(solution, [1])

        # check that masking best value works
        distances, trac_trans_vec, solution = sm._cart_dists(s1, s2, lattice, mask2, n1)
        assert_allclose(distances, [0])
        assert_allclose(trac_trans_vec, [0.02, 0.03, 0.04])
        assert_allclose(solution, [0])

        # check that averaging of translation is done properly
        distances, trac_trans_vec, solution = sm._cart_dists(s1, s3, lattice, mask3, n1)
        assert_allclose(distances, [0.08093341] * 2)
        assert_allclose(trac_trans_vec, [0.01, 0.025, 0.035])
        assert_allclose(solution, [1, 0])

        # check distances are large when mask allows no 'real' mapping
        distances, trac_trans_vec, solution = sm._cart_dists(s1, s4, lattice, mask4, n1)
        assert np.min(distances) > 1e8
        assert np.min(trac_trans_vec) > 1e8

    def test_get_mask(self):
        sm = StructureMatcher(comparator=ElementComparator())
        lattice = Lattice.cubic(1)
        s1 = Structure(lattice, ["Mg", "Cu", "Ag", "Cu"], [[0] * 3] * 4)
        s2 = Structure(lattice, ["Cu", "Cu", "Ag"], [[0] * 3] * 3)

        result = [
            [True, False, True, False],
            [True, False, True, False],
            [True, True, False, True],
        ]
        mask, inds, idx = sm._get_mask(s1, s2, 1, s1_supercell=True)
        assert np.all(mask == result)
        assert idx == 2
        assert inds == [2]

        # test supercell with match
        result = [
            [1, 1, 0, 0, 1, 1, 0, 0],
            [1, 1, 0, 0, 1, 1, 0, 0],
            [1, 1, 1, 1, 0, 0, 1, 1],
        ]
        mask, inds, idx = sm._get_mask(s1, s2, 2, s1_supercell=True)
        assert np.all(mask == result)
        assert idx == 2
        assert list(inds) == [4]

        # test supercell without match
        result = [
            [1, 1, 1, 1, 1, 1],
            [0, 0, 0, 0, 1, 1],
            [1, 1, 1, 1, 0, 0],
            [0, 0, 0, 0, 1, 1],
        ]
        mask, inds, idx = sm._get_mask(s2, s1, 2, s1_supercell=True)
        assert np.all(mask == result)
        assert idx == 0
        assert list(inds) == []

        # test s2_supercell
        result = [
            [1, 1, 1],
            [1, 1, 1],
            [0, 0, 1],
            [0, 0, 1],
            [1, 1, 0],
            [1, 1, 0],
            [0, 0, 1],
            [0, 0, 1],
        ]
        mask, inds, idx = sm._get_mask(s2, s1, 2, s1_supercell=False)
        assert np.all(mask == result)
        assert idx == 0
        assert list(inds) == []

        # test for multiple translation indices
        s1 = Structure(lattice, ["Cu", "Ag", "Cu", "Ag", "Ag"], [[0] * 3] * 5)
        s2 = Structure(lattice, ["Ag", "Cu", "Ag"], [[0] * 3] * 3)
        result = [[1, 0, 1, 0, 0], [0, 1, 0, 1, 1], [1, 0, 1, 0, 0]]
        mask, inds, idx = sm._get_mask(s1, s2, 1, s1_supercell=True)

        assert np.all(mask == result)
        assert idx == 1
        assert list(inds) == [0, 2]

    def test_get_supercells(self):
        sm = StructureMatcher(comparator=ElementComparator())
        lattice = Lattice.cubic(1)
        l2 = Lattice.cubic(0.5)
        s1 = Structure(lattice, ["Mg", "Cu", "Ag", "Cu"], [[0] * 3] * 4)
        s2 = Structure(l2, ["Cu", "Cu", "Ag"], [[0] * 3] * 3)
        scs = list(sm._get_supercells(s1, s2, fu=8, s1_supercell=False))
        for x in scs:
            assert abs(np.linalg.det(x[3])) == approx(8)
            assert len(x[0]) == 4
            assert len(x[1]) == 24
        assert len(scs) == 48

        scs = list(sm._get_supercells(s2, s1, fu=8, s1_supercell=True))
        for x in scs:
            assert abs(np.linalg.det(x[3])) == approx(8)
            assert len(x[0]) == 24
            assert len(x[1]) == 4
        assert len(scs) == 48

    def test_fit(self):
        """
        Take two known matched structures
            1) Ensure match
            2) Ensure match after translation and rotations
            3) Ensure no-match after large site translation
            4) Ensure match after site shuffling.
        """
        sm = StructureMatcher()

        assert sm.fit(self.struct_list[0], self.struct_list[1])

        # Test rotational/translational invariance
        op = SymmOp.from_axis_angle_and_translation([0, 0, 1], 30, translation_vec=[0.4, 0.7, 0.9])
        self.struct_list[1].apply_operation(op)
        assert sm.fit(self.struct_list[0], self.struct_list[1])

        # Test failure under large atomic translation
        self.struct_list[1].translate_sites([0], [0.4, 0.4, 0.2], frac_coords=True)
        assert not sm.fit(self.struct_list[0], self.struct_list[1])

        self.struct_list[1].translate_sites([0], [-0.4, -0.4, -0.2], frac_coords=True)
        # random.shuffle(editor._sites)
        assert sm.fit(self.struct_list[0], self.struct_list[1])
        # Test FrameworkComparator
        sm2 = StructureMatcher(comparator=FrameworkComparator())
        lfp = self.get_structure("LiFePO4")
        nfp = self.get_structure("NaFePO4")
        assert sm2.fit(lfp, nfp)
        assert not sm.fit(lfp, nfp)

        # Test anonymous fit.
        assert sm.fit_anonymous(lfp, nfp)
        assert sm.get_rms_anonymous(lfp, nfp)[0] == approx(0.060895871160262717)

        # Test partial occupancies.
        struct1 = Structure(
            Lattice.cubic(3),
            [{"Fe": 0.5}, {"Fe": 0.5}, {"Fe": 0.5}, {"Fe": 0.5}],
            [[0, 0, 0], [0.25, 0.25, 0.25], [0.5, 0.5, 0.5], [0.75, 0.75, 0.75]],
        )
        struct2 = Structure(
            Lattice.cubic(3),
            [{"Fe": 0.25}, {"Fe": 0.5}, {"Fe": 0.5}, {"Fe": 0.75}],
            [[0, 0, 0], [0.25, 0.25, 0.25], [0.5, 0.5, 0.5], [0.75, 0.75, 0.75]],
        )
        assert not sm.fit(struct1, struct2)
        assert not sm.fit(struct2, struct1)
        struct2 = Structure(
            Lattice.cubic(3),
            [{"Mn": 0.5}, {"Mn": 0.5}, {"Mn": 0.5}, {"Mn": 0.5}],
            [[0, 0, 0], [0.25, 0.25, 0.25], [0.5, 0.5, 0.5], [0.75, 0.75, 0.75]],
        )
        assert sm.fit_anonymous(struct1, struct2)

        assert sm.get_rms_anonymous(struct1, struct2)[0] == approx(0)

        # test symmetric
        sm_coarse = sm = StructureMatcher(comparator=ElementComparator(), ltol=0.6, stol=0.6, angle_tol=6)

        struct1 = Structure.from_file(f"{VASP_IN_DIR}/POSCAR_fit_symm_s1")
        struct2 = Structure.from_file(f"{VASP_IN_DIR}/POSCAR_fit_symm_s2")
        assert sm_coarse.fit(struct1, struct2)
        assert sm_coarse.fit(struct2, struct1) is False
        assert sm_coarse.fit(struct1, struct2, symmetric=True) is False
        assert sm_coarse.fit(struct2, struct1, symmetric=True) is False

    def test_oxi(self):
        """Test oxidation state removal matching."""
        sm = StructureMatcher()
        assert not sm.fit(self.oxi_structs[0], self.oxi_structs[1])
        sm = StructureMatcher(comparator=ElementComparator())
        assert sm.fit(self.oxi_structs[0], self.oxi_structs[1])

    def test_primitive(self):
        """Test primitive cell reduction."""
        sm = StructureMatcher(primitive_cell=True)
        self.struct_list[1].make_supercell([[2, 0, 0], [0, 3, 0], [0, 0, 1]])
        assert sm.fit(self.struct_list[0], self.struct_list[1])

    def test_class(self):
        # Tests entire class as single working unit
        sm = StructureMatcher()
        # Test group_structures and find_indices
        out = sm.group_structures(self.struct_list)
        assert list(map(len, out)) == [4, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1]
        assert sum(map(len, out)) == len(self.struct_list)
        for s in self.struct_list[::2]:
            s.replace_species({"Ti": "Zr", "O": "Ti"})
        out = sm.group_structures(self.struct_list, anonymous=True)
        assert list(map(len, out)) == [4, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1]

    def test_mix(self):
        structures = list(map(self.get_structure, ["Li2O", "Li2O2", "LiFePO4"]))
        structures += [Structure.from_file(f"{VASP_IN_DIR}/{fname}") for fname in ["POSCAR_Li2O", "POSCAR_LiFePO4"]]
        sm = StructureMatcher(comparator=ElementComparator())
        groups = sm.group_structures(structures)
        for group in groups:
            formula = group[0].reduced_formula
            assert len(group) == (2 if formula in {"Li2O", "LiFePO4"} else 1)

    def test_left_handed_lattice(self):
        """Ensure Left handed lattices are accepted."""
        sm = StructureMatcher()
        struct = Structure.from_file(f"{TEST_DIR}/Li3GaPCO7.json")
        assert sm.fit(struct, struct)

    def test_as_dict_and_from_dict(self):
        sm = StructureMatcher(
            ltol=0.1,
            stol=0.2,
            angle_tol=2,
            primitive_cell=False,
            scale=False,
            comparator=FrameworkComparator(),
        )
        dct = sm.as_dict()
        sm2 = StructureMatcher.from_dict(dct)
        assert sm2.as_dict() == dct

    def test_no_scaling(self):
        sm = StructureMatcher(ltol=0.1, stol=0.1, angle_tol=2, scale=False, comparator=ElementComparator())
        assert sm.fit(self.struct_list[0], self.struct_list[1])

        assert sm.get_rms_dist(self.struct_list[0], self.struct_list[1])[0] < 0.0008

    def test_supercell_fit(self):
        sm = StructureMatcher(attempt_supercell=False)
        s1 = Structure.from_file(f"{TEST_DIR}/Al3F9.json")
        s2 = Structure.from_file(f"{TEST_DIR}/Al3F9_distorted.json")

        assert not sm.fit(s1, s2)

        sm = StructureMatcher(attempt_supercell=True)

        assert sm.fit(s1, s2)
        assert sm.fit(s2, s1)

    def test_get_lattices(self):
        sm = StructureMatcher(
            ltol=0.2,
            stol=0.3,
            angle_tol=5,
            primitive_cell=True,
            scale=True,
            attempt_supercell=False,
        )
        l1 = Lattice.from_parameters(1, 2.1, 1.9, 90, 89, 91)
        l2 = Lattice.from_parameters(1.1, 2, 2, 89, 91, 90)
        s1 = Structure(l1, [], [])
        s2 = Structure(l2, [], [])

        lattices = list(sm._get_lattices(s=s1, target_lattice=s2.lattice))
        assert len(lattices) == 16

        l3 = Lattice.from_parameters(1.1, 2, 20, 89, 91, 90)
        s3 = Structure(l3, [], [])

        lattices = list(sm._get_lattices(s=s1, target_lattice=s3.lattice))
        assert len(lattices) == 0

    def test_find_match1(self):
        sm = StructureMatcher(
            ltol=0.2,
            stol=0.3,
            angle_tol=5,
            primitive_cell=True,
            scale=True,
            attempt_supercell=False,
        )
        lattice = Lattice.orthorhombic(1, 2, 3)
        s1 = Structure(lattice, ["Si", "Si", "Ag"], [[0, 0, 0.1], [0, 0, 0.2], [0.7, 0.4, 0.5]])
        s2 = Structure(
            lattice,
            ["Si", "Si", "Ag"],
            [[0, 0.1, 0], [0, 0.1, -0.95], [0.7, 0.5, 0.375]],
        )

        s1, s2, fu, s1_supercell = sm._preprocess(s1, s2, niggli=False)
        assert s1_supercell is True
        match = sm._strict_match(s1, s2, fu, s1_supercell=True, use_rms=True, break_on_match=False)
        scale_matrix = match[2]
        s2.make_supercell(scale_matrix)
        fc = s2.frac_coords + match[3]
        fc -= np.round(fc)
        assert np.sum(fc) == approx(0.9)
        assert np.sum(fc[:, :2]) == approx(0.1)
        cart_dist = np.sum(match[1] * (lattice.volume / 3) ** (1 / 3))
        assert cart_dist == approx(0.15)

    def test_find_match2(self):
        sm = StructureMatcher(
            ltol=0.2,
            stol=0.3,
            angle_tol=5,
            primitive_cell=True,
            scale=True,
            attempt_supercell=False,
        )
        lattice = Lattice.orthorhombic(1, 2, 3)
        s1 = Structure(lattice, ["Si", "Si"], [[0, 0, 0.1], [0, 0, 0.2]])
        s2 = Structure(lattice, ["Si", "Si"], [[0, 0.1, 0], [0, 0.1, -0.95]])

        s1, s2, fu, _s1_supercell = sm._preprocess(s1, s2, niggli=False)

        match = sm._strict_match(s1, s2, fu, s1_supercell=False, use_rms=True, break_on_match=False)
        scale_matrix = match[2]
        s2.make_supercell(scale_matrix)
        s2.translate_sites(range(len(s2)), match[3])

        assert np.sum(s2.frac_coords) % 1 == approx(0.3)
        assert np.sum(s2.frac_coords[:, :2]) % 1 == approx(0)

    def test_supercell_subsets(self):
        sm = StructureMatcher(
            ltol=0.2,
            stol=0.3,
            angle_tol=5,
            primitive_cell=False,
            scale=True,
            attempt_supercell=True,
            allow_subset=True,
            supercell_size="volume",
        )
        sm_no_s = StructureMatcher(
            ltol=0.2,
            stol=0.3,
            angle_tol=5,
            primitive_cell=False,
            scale=True,
            attempt_supercell=True,
            allow_subset=False,
            supercell_size="volume",
        )
        lattice = Lattice.orthorhombic(1, 2, 3)
        s1 = Structure(lattice, ["Ag", "Si", "Si"], [[0.7, 0.4, 0.5], [0, 0, 0.1], [0, 0, 0.2]])
        s1.make_supercell([2, 1, 1])
        s2 = Structure(
            lattice,
            ["Si", "Si", "Ag"],
            [[0, 0.1, -0.95], [0, 0.1, 0], [-0.7, 0.5, 0.375]],
        )

        shuffle = [0, 2, 1, 3, 4, 5]
        s1 = Structure.from_sites([s1[i] for i in shuffle])

        # test when s1 is exact supercell of s2
        result = sm.get_s2_like_s1(s1, s2)
        for a, b in zip(s1, result, strict=True):
            assert a.distance(b) < 0.08
            assert a.species == b.species

        assert sm.fit(s1, s2)
        assert sm.fit(s2, s1)
        assert sm_no_s.fit(s1, s2)
        assert sm_no_s.fit(s2, s1)

        rms = (0.048604032430991401, 0.059527539448807391)
        assert_allclose(sm.get_rms_dist(s1, s2), rms)
        assert_allclose(sm.get_rms_dist(s2, s1), rms)

        # test when the supercell is a subset of s2
        subset_supercell = s1.copy()
        del subset_supercell[0]
        result = sm.get_s2_like_s1(subset_supercell, s2)
        assert len(result) == 6
        for a, b in zip(subset_supercell, result, strict=False):
            assert a.distance(b) < 0.08
            assert a.species == b.species

        assert sm.fit(subset_supercell, s2)
        assert sm.fit(s2, subset_supercell)
        assert not sm_no_s.fit(subset_supercell, s2)
        assert not sm_no_s.fit(s2, subset_supercell)

        rms = (0.053243049896333279, 0.059527539448807336)
        assert_allclose(sm.get_rms_dist(subset_supercell, s2), rms)
        assert_allclose(sm.get_rms_dist(s2, subset_supercell), rms)

        # test when s2 (once made a supercell) is a subset of s1
        s2_missing_site = s2.copy()
        del s2_missing_site[1]
        result = sm.get_s2_like_s1(s1, s2_missing_site)
        for a, b in zip((s1[i] for i in (0, 2, 4, 5)), result, strict=True):
            assert a.distance(b) < 0.08
            assert a.species == b.species

        assert sm.fit(s1, s2_missing_site)
        assert sm.fit(s2_missing_site, s1)
        assert not sm_no_s.fit(s1, s2_missing_site)
        assert not sm_no_s.fit(s2_missing_site, s1)

        rms = (0.029763769724403633, 0.029763769724403987)
        assert_allclose(sm.get_rms_dist(s1, s2_missing_site), rms)
        assert_allclose(sm.get_rms_dist(s2_missing_site, s1), rms)

    def test_get_s2_large_s2(self):
        sm = StructureMatcher(
            ltol=0.2,
            stol=0.3,
            angle_tol=5,
            primitive_cell=False,
            scale=False,
            attempt_supercell=True,
            allow_subset=False,
            supercell_size="volume",
        )

        lattice = Lattice.orthorhombic(1, 2, 3)
        s1 = Structure(lattice, ["Ag", "Si", "Si"], [[0.7, 0.4, 0.5], [0, 0, 0.1], [0, 0, 0.2]])

        l2 = Lattice.orthorhombic(1.01, 2.01, 3.01)
        s2 = Structure(l2, ["Si", "Si", "Ag"], [[0, 0.1, -0.95], [0, 0.1, 0], [-0.7, 0.5, 0.375]])
        s2.make_supercell([[0, -1, 0], [1, 0, 0], [0, 0, 1]])

        result = sm.get_s2_like_s1(s1, s2)

        for x, y in zip(s1, result, strict=True):
            assert x.distance(y) < 0.08

    def test_get_mapping(self):
        sm = StructureMatcher(
            ltol=0.2,
            stol=0.3,
            angle_tol=5,
            primitive_cell=False,
            scale=True,
            attempt_supercell=False,
            allow_subset=True,
        )
        lattice = Lattice.orthorhombic(1, 2, 3)
        struct1 = Structure(lattice, ["Ag", "Si", "Si"], [[0.7, 0.4, 0.5], [0, 0, 0.1], [0, 0, 0.2]])
        struct1.make_supercell([2, 1, 1])
        struct2 = Structure(
            lattice,
            ["Si", "Si", "Ag"],
            [[0, 0.1, -0.95], [0, 0.1, 0], [-0.7, 0.5, 0.375]],
        )

        shuffle = [2, 0, 1, 3, 5, 4]
        struct1 = Structure.from_sites([struct1[i] for i in shuffle])
        # test the mapping
        struct2.make_supercell([2, 1, 1])
        # equal sizes
        for ii, jj in enumerate(sm.get_mapping(struct1, struct2)):
            assert struct1[jj].species == struct2[ii].species

        del struct1[0]
        # s1 is subset of s2
        for ii, jj in enumerate(sm.get_mapping(struct2, struct1)):
            assert struct1[ii].species == struct2[jj].species
        # s2 is smaller than s1
        del struct2[0]
        del struct2[1]
        with pytest.raises(ValueError, match="subset is larger than superset"):
            sm.get_mapping(struct2, struct1)

    def test_get_supercell_matrix(self):
        sm = StructureMatcher(
            ltol=0.1,
            stol=0.3,
            angle_tol=2,
            primitive_cell=False,
            scale=True,
            attempt_supercell=True,
        )

        lattice = Lattice.orthorhombic(1, 2, 3)

        s1 = Structure(lattice, ["Si", "Si", "Ag"], [[0, 0, 0.1], [0, 0, 0.2], [0.7, 0.4, 0.5]])
        s1.make_supercell([2, 1, 1])
        s2 = Structure(
            lattice,
            ["Si", "Si", "Ag"],
            [[0, 0.1, 0], [0, 0.1, -0.95], [-0.7, 0.5, 0.375]],
        )
        result = sm.get_supercell_matrix(s1, s2)
        assert (result == [[-2, 0, 0], [0, 1, 0], [0, 0, 1]]).all()

        s1 = Structure(lattice, ["Si", "Si", "Ag"], [[0, 0, 0.1], [0, 0, 0.2], [0.7, 0.4, 0.5]])
        s1.make_supercell([[1, -1, 0], [0, 0, -1], [0, 1, 0]])

        s2 = Structure(
            lattice,
            ["Si", "Si", "Ag"],
            [[0, 0.1, 0], [0, 0.1, -0.95], [-0.7, 0.5, 0.375]],
        )
        result = sm.get_supercell_matrix(s1, s2)
        assert (result == [[-1, -1, 0], [0, 0, -1], [0, 1, 0]]).all()

        # test when the supercell is a subset
        sm = StructureMatcher(
            ltol=0.1,
            stol=0.3,
            angle_tol=2,
            primitive_cell=False,
            scale=True,
            attempt_supercell=True,
            allow_subset=True,
        )
        del s1[0]
        result = sm.get_supercell_matrix(s1, s2)
        assert (result == [[-1, -1, 0], [0, 0, -1], [0, 1, 0]]).all()

    def test_subset(self):
        sm = StructureMatcher(
            ltol=0.2,
            stol=0.3,
            angle_tol=5,
            primitive_cell=False,
            scale=True,
            attempt_supercell=False,
            allow_subset=True,
        )
        lattice = Lattice.orthorhombic(10, 20, 30)
        s1 = Structure(lattice, ["Si", "Si", "Ag"], [[0, 0, 0.1], [0, 0, 0.2], [0.7, 0.4, 0.5]])
        s2 = Structure(lattice, ["Si", "Ag"], [[0, 0.1, 0], [-0.7, 0.5, 0.4]])
        result = sm.get_s2_like_s1(s1, s2)

        assert len(find_in_coord_list_pbc(result.frac_coords, [0, 0, 0.1])) == 1
        assert len(find_in_coord_list_pbc(result.frac_coords, [0.7, 0.4, 0.5])) == 1

        # test with fewer species in s2
        s1 = Structure(lattice, ["Si", "Ag", "Si"], [[0, 0, 0.1], [0, 0, 0.2], [0.7, 0.4, 0.5]])
        s2 = Structure(lattice, ["Si", "Si"], [[0, 0.1, 0], [-0.7, 0.5, 0.4]])
        result = sm.get_s2_like_s1(s1, s2)
        mindists = np.min(s1.lattice.get_all_distances(s1.frac_coords, result.frac_coords), axis=0)
        assert np.max(mindists) < 1e-6

        assert len(find_in_coord_list_pbc(result.frac_coords, [0, 0, 0.1])) == 1
        assert len(find_in_coord_list_pbc(result.frac_coords, [0.7, 0.4, 0.5])) == 1

        # test with not enough sites in s1
        # test with fewer species in s2
        s1 = Structure(lattice, ["Si", "Ag", "Cl"], [[0, 0, 0.1], [0, 0, 0.2], [0.7, 0.4, 0.5]])
        s2 = Structure(lattice, ["Si", "Si"], [[0, 0.1, 0], [-0.7, 0.5, 0.4]])
        assert sm.get_s2_like_s1(s1, s2) is None

    def test_out_of_cell_s2_like_s1(self):
        lattice = Lattice.cubic(5)
        s1 = Structure(lattice, ["Si", "Ag", "Si"], [[0, 0, -0.02], [0, 0, 0.001], [0.7, 0.4, 0.5]])
        s2 = Structure(lattice, ["Si", "Ag", "Si"], [[0, 0, 0.98], [0, 0, 0.99], [0.7, 0.4, 0.5]])
        new_s2 = StructureMatcher(primitive_cell=False).get_s2_like_s1(s1, s2)
        dists = np.sum((s1.cart_coords - new_s2.cart_coords) ** 2, axis=-1) ** 0.5
        assert np.max(dists) < 0.1

    def test_disordered_primitive_to_ordered_supercell(self):
        sm_atoms = StructureMatcher(
            ltol=0.2,
            stol=0.3,
            angle_tol=5,
            primitive_cell=False,
            scale=True,
            attempt_supercell=True,
            allow_subset=True,
            supercell_size="num_atoms",
            comparator=OrderDisorderElementComparator(),
        )
        sm_sites = StructureMatcher(
            ltol=0.2,
            stol=0.3,
            angle_tol=5,
            primitive_cell=False,
            scale=True,
            attempt_supercell=True,
            allow_subset=True,
            supercell_size="num_sites",
            comparator=OrderDisorderElementComparator(),
        )
        lp = Lattice.orthorhombic(10, 20, 30)
        pcoords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        ls = Lattice.orthorhombic(20, 20, 30)
        scoords = [[0, 0, 0], [0.75, 0.5, 0.5]]
        prim = Structure(lp, [{"Na": 0.5}, {"Cl": 0.5}], pcoords)
        supercell = Structure(ls, ["Na", "Cl"], scoords)
        supercell.make_supercell([[-1, 1, 0], [0, 1, 1], [1, 0, 0]])

        assert not sm_sites.fit(prim, supercell)
        assert sm_atoms.fit(prim, supercell)

        with pytest.raises(ValueError, match="Struct1 must be the supercell, not the other way around"):
            sm_atoms.get_s2_like_s1(prim, supercell)
        assert len(sm_atoms.get_s2_like_s1(supercell, prim)) == 4

    def test_ordered_primitive_to_disordered_supercell(self):
        sm_atoms = StructureMatcher(
            ltol=0.2,
            stol=0.3,
            angle_tol=5,
            primitive_cell=False,
            scale=True,
            attempt_supercell=True,
            allow_subset=True,
            supercell_size="num_atoms",
            comparator=OrderDisorderElementComparator(),
        )
        sm_sites = StructureMatcher(
            ltol=0.2,
            stol=0.3,
            angle_tol=5,
            primitive_cell=False,
            scale=True,
            attempt_supercell=True,
            allow_subset=True,
            supercell_size="num_sites",
            comparator=OrderDisorderElementComparator(),
        )
        lp = Lattice.orthorhombic(10, 20, 30)
        pcoords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        ls = Lattice.orthorhombic(20, 20, 30)
        scoords = [[0, 0, 0], [0.5, 0, 0], [0.25, 0.5, 0.5], [0.75, 0.5, 0.5]]
        s1 = Structure(lp, ["Na", "Cl"], pcoords)
        s2 = Structure(ls, [{"Na": 0.5}, {"Na": 0.5}, {"Cl": 0.5}, {"Cl": 0.5}], scoords)

        assert sm_sites.fit(s1, s2)
        assert not sm_atoms.fit(s1, s2)

    def test_disordered_to_disordered(self):
        sm_atoms = StructureMatcher(
            ltol=0.2,
            stol=0.3,
            angle_tol=5,
            primitive_cell=False,
            scale=True,
            attempt_supercell=True,
            allow_subset=False,
            comparator=OrderDisorderElementComparator(),
        )
        lp = Lattice.orthorhombic(10, 20, 30)
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        s1 = Structure(lp, [{"Na": 0.5, "Cl": 0.5}, {"Na": 0.5, "Cl": 0.5}], coords)
        s2 = Structure(lp, [{"Na": 0.5, "Cl": 0.5}, {"Na": 0.5, "Br": 0.5}], coords)

        assert not sm_atoms.fit(s1, s2)

    def test_occupancy_comparator(self):
        lp = Lattice.orthorhombic(10, 20, 30)
        pcoords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        s1 = Structure(lp, [{"Na": 0.6, "K": 0.4}, "Cl"], pcoords)
        s2 = Structure(lp, [{"Xa": 0.4, "Xb": 0.6}, "Cl"], pcoords)
        s3 = Structure(lp, [{"Xa": 0.5, "Xb": 0.5}, "Cl"], pcoords)

        sm_sites = StructureMatcher(
            ltol=0.2,
            stol=0.3,
            angle_tol=5,
            primitive_cell=False,
            scale=True,
            attempt_supercell=True,
            allow_subset=True,
            supercell_size="num_sites",
            comparator=OccupancyComparator(),
        )

        assert sm_sites.fit(s1, s2)
        assert not sm_sites.fit(s1, s3)

    def test_electronegativity(self):
        sm = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5)

        s1 = Structure.from_file(f"{TEST_DIR}/Na2Fe2PAsO4S4.json")
        s2 = Structure.from_file(f"{TEST_DIR}/Na2Fe2PNO4Se4.json")
        assert sm.get_best_electronegativity_anonymous_mapping(s1, s2) == {
            Element("S"): Element("Se"),
            Element("As"): Element("N"),
            Element("Fe"): Element("Fe"),
            Element("Na"): Element("Na"),
            Element("P"): Element("P"),
            Element("O"): Element("O"),
        }
        assert len(sm.get_all_anonymous_mappings(s1, s2)) == 2

        # test include_dist
        dists = {Element("N"): 0, Element("P"): 0.0010725064}
        for mapping, d in sm.get_all_anonymous_mappings(s1, s2, include_dist=True):
            assert dists[mapping[Element("As")]] == approx(d)

    def test_rms_vs_minimax(self):
        # This tests that structures with adjusted RMS less than stol, but minimax
        # greater than stol are treated properly
        # stol=0.3 gives exactly an ftol of 0.1 on the c axis
        sm = StructureMatcher(ltol=0.2, stol=0.301, angle_tol=1, primitive_cell=False)
        lattice = Lattice.orthorhombic(1, 2, 12)

        sp = ["Si", "Si", "Al"]
        s1 = Structure(lattice, sp, np.diag((0.5, 0, 0.5)))
        s2 = Structure(lattice, sp, np.diag((0.5, 0, 0.6)))
        assert_allclose(sm.get_rms_dist(s1, s2), (0.32**0.5 / 2, 0.4))

        assert sm.fit(s1, s2) is False
        assert sm.fit_anonymous(s1, s2) is False
        assert sm.get_mapping(s1, s2) is None
