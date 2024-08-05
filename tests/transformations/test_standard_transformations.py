from __future__ import annotations

import functools
import json
import operator
from shutil import which
from unittest import TestCase

import numpy as np
import pytest
from monty.json import MontyDecoder
from pytest import approx

from pymatgen.core import Element, PeriodicSite
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.transformations.standard_transformations import (
    AutoOxiStateDecorationTransformation,
    ChargedCellTransformation,
    ConventionalCellTransformation,
    DeformStructureTransformation,
    DiscretizeOccupanciesTransformation,
    EwaldSummation,
    Fraction,
    OrderDisorderedStructureTransformation,
    OxidationStateDecorationTransformation,
    OxidationStateRemovalTransformation,
    PartialRemoveSpecieTransformation,
    PerturbStructureTransformation,
    PrimitiveCellTransformation,
    RemoveSpeciesTransformation,
    RotationTransformation,
    ScaleToRelaxedTransformation,
    Structure,
    SubstitutionTransformation,
    SupercellTransformation,
)
from pymatgen.util.testing import TEST_FILES_DIR, VASP_IN_DIR

enumlib_present = which("enum.x") and which("makestr.x")


class TestRotationTransformations(TestCase):
    def setUp(self):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = [
            [3.8401979337, 0, 0],
            [1.9200989668, 3.3257101909, 0],
            [0, -2.2171384943, 3.1355090603],
        ]
        self.struct = Structure(lattice, ["Si"] * 2, coords)

    def test_as_from_dict(self):
        trafo = RotationTransformation([0, 1, 0], 30)
        dct = trafo.as_dict()
        assert isinstance(RotationTransformation.from_dict(dct), RotationTransformation)

    def test_rotation_transformation(self):
        trafo = RotationTransformation([0, 1, 0], 30)
        s2 = trafo.apply_transformation(self.struct)
        s1 = trafo.inverse.apply_transformation(s2)
        assert (abs(s1.lattice.matrix - self.struct.lattice.matrix) < 1e-8).all()


class TestRemoveSpeciesTransformation:
    def test_apply_transformation(self):
        trafo = RemoveSpeciesTransformation(["Li+"])
        coords = [[0, 0, 0], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]]
        lattice = [
            [3.8401979337, 0, 0],
            [1.9200989668, 3.3257101909, 0],
            [0, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        struct = trafo.apply_transformation(struct)
        assert struct.formula == "O2"

        dct = trafo.as_dict()
        assert isinstance(RemoveSpeciesTransformation.from_dict(dct), RemoveSpeciesTransformation)


class TestSubstitutionTransformation:
    def test_apply_transformation(self):
        trafo = SubstitutionTransformation({"Li+": "Na+", "O2-": "S2-"})
        coords = [[0, 0, 0], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]]
        lattice = [
            [3.8401979337, 0, 0],
            [1.9200989668, 3.3257101909, 0],
            [0, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        struct_trafo = trafo.apply_transformation(struct)
        assert struct_trafo.formula == "Na2 S2"

    def test_fractional_substitution(self):
        trafo = SubstitutionTransformation({"Li+": "Na+", "O2-": {"S2-": 0.5, "Se2-": 0.5}})
        # test the to and from dict on the nested dictionary
        trafo = SubstitutionTransformation.from_dict(trafo.as_dict())
        coords = [[0, 0, 0], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]]
        lattice = [
            [3.8401979337, 0, 0],
            [1.9200989668, 3.3257101909, 0],
            [0, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        struct_trafo = trafo.apply_transformation(struct)
        assert struct_trafo.formula == "Na2 Se1 S1"


class TestSupercellTransformation(TestCase):
    def setUp(self):
        coords = [[0, 0, 0], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]]
        lattice = [
            [3.8401979337, 0, 0],
            [1.9200989668, 3.3257101909, 0],
            [0, -2.2171384943, 3.1355090603],
        ]
        self.struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)

    def test_apply_transformation(self):
        trafo = SupercellTransformation([[2, 1, 0], [0, 2, 0], [1, 0, 2]])
        struct = trafo.apply_transformation(self.struct)
        assert struct.formula == "Li16 O16"

    def test_from_scaling_factors(self):
        scale_factors = np.random.default_rng().integers(1, 5, 3)
        trafo = SupercellTransformation.from_scaling_factors(*scale_factors)
        struct = trafo.apply_transformation(self.struct)
        assert len(struct) == 4 * functools.reduce(operator.mul, scale_factors)

    def test_from_boundary_distance(self):
        struct_cubic = Structure.from_spacegroup("Pm-3m", 4 * np.eye(3), ["H"], [[0, 0, 0]])

        for struct in (struct_cubic, self.struct):
            for min_dist in range(6, 19, 4):
                trafo = SupercellTransformation.from_boundary_distance(
                    structure=struct, min_boundary_dist=min_dist, allow_rotation=False
                )
                trafo_allow_rotate = SupercellTransformation.from_boundary_distance(
                    structure=struct, min_boundary_dist=min_dist, allow_rotation=True
                )
                struct_super = trafo.apply_transformation(struct.copy())
                struct_super_allow_rotate = trafo_allow_rotate.apply_transformation(struct.copy())

                if min_dist == 9 and struct is struct_cubic:
                    assert len(struct_super_allow_rotate) == 14
                    assert len(struct_super) == 27

                min_expand_allow_rotate_true = np.int8(
                    min_dist / np.array([struct_super_allow_rotate.lattice.d_hkl(plane) for plane in np.eye(3)])
                )
                min_expand_allow_rotate_false = np.int8(
                    min_dist / np.array([struct_super.lattice.d_hkl(plane) for plane in np.eye(3)])
                )

                assert sum(min_expand_allow_rotate_true != 0) == 0
                assert sum(min_expand_allow_rotate_false != 0) == 0
                assert len(struct_super_allow_rotate) <= len(struct_super)

        max_atoms = 10
        with pytest.raises(
            RuntimeError,
            match=f"{max_atoms=} exceeded while trying to solve for supercell. "
            f"You can try lowering min_boundary_dist=6",
        ):
            SupercellTransformation.from_boundary_distance(structure=self.struct, max_atoms=max_atoms)


class TestOxidationStateDecorationTransformation:
    def test_apply_transformation(self):
        trafo = OxidationStateDecorationTransformation({"Li": 1, "O": -2})
        coords = [[0, 0, 0], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]]
        lattice = [
            [3.8401979337, 0, 0],
            [1.9200989668, 3.3257101909, 0],
            [0, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, ["Li", "Li", "O", "O"], coords)
        struct_trafo = trafo.apply_transformation(struct)
        assert struct_trafo[0].species_string == "Li+"
        assert struct_trafo[2].species_string == "O2-"
        dct = trafo.as_dict()
        assert isinstance(OxidationStateDecorationTransformation.from_dict(dct), OxidationStateDecorationTransformation)


class TestAutoOxiStateDecorationTransformation:
    def test_apply_transformation(self):
        trafo = AutoOxiStateDecorationTransformation()
        struct = trafo.apply_transformation(Structure.from_file(f"{VASP_IN_DIR}/POSCAR_LiFePO4"))
        expected_oxi = {"Li": 1, "P": 5, "O": -2, "Fe": 2}
        for site in struct:
            assert site.specie.oxi_state == expected_oxi[site.specie.symbol]

    def test_as_from_dict(self):
        trafo = AutoOxiStateDecorationTransformation()
        dct = trafo.as_dict()
        trafo = AutoOxiStateDecorationTransformation.from_dict(dct)
        assert trafo.analyzer.dist_scale_factor == 1.015


class TestOxidationStateRemovalTransformation:
    def test_apply_transformation(self):
        trafo = OxidationStateRemovalTransformation()
        coords = [[0, 0, 0], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]]
        lattice = [
            [3.8401979337, 0, 0],
            [1.9200989668, 3.3257101909, 0],
            [0, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        struct_trafo = trafo.apply_transformation(struct)
        assert struct_trafo[0].species_string == "Li"
        assert struct_trafo[2].species_string == "O"

        dct = trafo.as_dict()
        assert isinstance(OxidationStateRemovalTransformation.from_dict(dct), OxidationStateRemovalTransformation)


@pytest.mark.skipif(not enumlib_present, reason="enum_lib not present.")
class TestPartialRemoveSpecieTransformation:
    def test_apply_transformation(self):
        trafo = PartialRemoveSpecieTransformation("Li+", 1.0 / 3, 3)
        coords = [[0, 0, 0], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]]
        lattice = [
            [3.8401979337, 0, 0],
            [1.9200989668, 3.3257101909, 0],
            [0, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, ["Li+", "Li+", "Li+", "O2-"], coords)
        assert len(trafo.apply_transformation(struct, 100)) == 2

        dct = trafo.as_dict()
        assert isinstance(PartialRemoveSpecieTransformation.from_dict(dct), PartialRemoveSpecieTransformation)

    def test_apply_transformation_fast(self):
        trafo = PartialRemoveSpecieTransformation("Li+", 0.5)
        coords = [[0, 0, 0], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25], [0.1, 0.1, 0.1], [0.3, 0.75, 0.3]]
        lattice = Lattice([[10, 0, 0], [0, 10, 0], [0, 0, 10]])
        struct = Structure(lattice, ["Li+"] * 6, coords)
        fast_opt_s = trafo.apply_transformation(struct)
        trafo = PartialRemoveSpecieTransformation("Li+", 0.5, PartialRemoveSpecieTransformation.ALGO_COMPLETE)
        slow_opt_s = trafo.apply_transformation(struct)
        assert EwaldSummation(fast_opt_s).total_energy == approx(EwaldSummation(slow_opt_s).total_energy, abs=1e-4)
        assert fast_opt_s == slow_opt_s

    def test_apply_transformations_complete_ranking(self):
        t1 = OxidationStateDecorationTransformation({"Li": 1, "Fe": 2, "P": 5, "O": -2})
        struct = t1.apply_transformation(Structure.from_file(f"{VASP_IN_DIR}/POSCAR_LiFePO4"))
        trafo = PartialRemoveSpecieTransformation("Li+", 0.5, PartialRemoveSpecieTransformation.ALGO_COMPLETE)
        assert len(trafo.apply_transformation(struct, 10)) == 6

    def test_apply_transformations_best_first(self):
        t1 = OxidationStateDecorationTransformation({"Li": 1, "Fe": 2, "P": 5, "O": -2})
        struct = t1.apply_transformation(Structure.from_file(f"{VASP_IN_DIR}/POSCAR_LiFePO4"))
        trafo = PartialRemoveSpecieTransformation("Li+", 0.5, PartialRemoveSpecieTransformation.ALGO_BEST_FIRST)
        assert len(trafo.apply_transformation(struct)) == 26


class TestOrderDisorderedStructureTransformation:
    def test_apply_transformation(self):
        trafo = OrderDisorderedStructureTransformation()
        coords = [[0, 0, 0], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]]
        lattice = [
            [3.8401979337, 0, 0],
            [1.9200989668, 3.3257101909, 0],
            [0, -2.2171384943, 3.1355090603],
        ]

        struct = Structure(
            lattice,
            [
                {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
            ],
            coords,
        )
        output = trafo.apply_transformation(struct, return_ranked_list=50)
        assert len(output) == 12
        assert isinstance(output[0]["structure"], Structure)

        struct = Structure(
            lattice,
            [
                {"Si4+": 0.5},
                {"Si4+": 0.5},
                {"P5+": 0.5, "O2-": 0.5},
                {"P5+": 0.5, "O2-": 0.5},
            ],
            coords,
        )
        output = trafo.apply_transformation(struct, return_ranked_list=50)
        assert isinstance(output, list)
        assert len(output) == 4
        assert trafo.lowest_energy_structure == output[0]["structure"]

        struct = Structure(lattice, [{"Si4+": 0.5}, {"Si4+": 0.5}, {"O2-": 0.5}, {"O2-": 0.5}], coords)
        all_structs = trafo.apply_transformation(struct, 50)
        assert len(all_structs) == 4

        struct = Structure(lattice, [{"Si4+": 0.333}, {"Si4+": 0.333}, {"Si4+": 0.333}, "O2-"], coords)
        all_structs = trafo.apply_transformation(struct, 50)
        assert len(all_structs) == 3

        dct = trafo.as_dict()
        assert isinstance(OrderDisorderedStructureTransformation.from_dict(dct), OrderDisorderedStructureTransformation)

    def test_no_oxidation(self):
        specie = {"Cu1+": 0.5, "Au2+": 0.5}
        cu_au = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3.677), [specie], [[0, 0, 0]])
        trans = OrderDisorderedStructureTransformation()
        ss = trans.apply_transformation(cu_au, return_ranked_list=100)
        assert ss[0]["structure"].composition["Cu+"] == 2
        trans = OrderDisorderedStructureTransformation(no_oxi_states=True)
        ss = trans.apply_transformation(cu_au, return_ranked_list=100)
        assert ss[0]["structure"].composition["Cu+"] == 0
        assert ss[0]["structure"].composition["Cu"] == 2

    def test_symmetrized_structure(self):
        trafo = OrderDisorderedStructureTransformation(symmetrized_structures=True)
        lattice = Lattice.cubic(5)
        coords = [[0.5, 0.5, 0.5], [0.45, 0.45, 0.45], [0.56, 0.56, 0.56], [0.25, 0.75, 0.75], [0.75, 0.25, 0.25]]
        struct = Structure(lattice, [{"Si4+": 1}, *[{"Si4+": 0.5}] * 4], coords)
        test_site = PeriodicSite("Si4+", coords[2], lattice)
        struct = SymmetrizedStructure(struct, "not_real", [0, 1, 1, 2, 2], ["a", "b", "b", "c", "c"])
        output = trafo.apply_transformation(struct)
        assert test_site in output

    def test_too_small_cell(self):
        trafo = OrderDisorderedStructureTransformation()
        coords = [[0.5, 0.5, 0.5]]
        lattice = Lattice([[3.8401979337, 0, 0], [1.9200989668, 3.3257101909, 0], [0, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, [{"X4+": 0.33, "O2-": 0.33, "P5+": 0.33}], coords)
        with pytest.raises(ValueError, match="Occupancy fractions not consistent with size of unit cell"):
            trafo.apply_transformation(struct)

    def test_best_first(self):
        trafo = OrderDisorderedStructureTransformation(algo=2)
        coords = [[0, 0, 0], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]]
        lattice = [
            [3.8401979337, 0, 0],
            [1.9200989668, 3.3257101909, 0],
            [0, -2.2171384943, 3.1355090603],
        ]

        struct = Structure(
            lattice,
            [
                {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
            ],
            coords,
        )
        output = trafo.apply_transformation(struct, return_ranked_list=3)
        assert output[0]["energy"] == approx(-234.57813667648315, abs=1e-4)


class TestPrimitiveCellTransformation:
    def test_apply_transformation(self):
        trafo = PrimitiveCellTransformation()
        coords = [
            [0, 0, 0],
            [0.375, 0.375, 0.375],
            [0.5, 0.5, 0.5],
            [0.875, 0.875, 0.875],
            [0.125, 0.125, 0.125],
            [0.25, 0.25, 0.25],
            [0.625, 0.625, 0.625],
            [0.75, 0.75, 0.75],
        ]

        lattice = [
            [3.8401979337, 0, 0],
            [1.9200989668, 3.3257101909, 0],
            [0, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-", "O2-", "O2-", "O2-"], coords)
        struct = trafo.apply_transformation(struct)
        assert len(struct) == 4

        with open(f"{TEST_FILES_DIR}/transformations/TiO2_super.json") as file:
            struct = json.load(file, cls=MontyDecoder)
            prim = trafo.apply_transformation(struct)
            assert prim.formula == "Ti4 O8"

        dct = trafo.as_dict()
        assert isinstance(PrimitiveCellTransformation.from_dict(dct), PrimitiveCellTransformation)


class TestConventionalCellTransformation:
    def test_apply_transformation(self):
        trafo = ConventionalCellTransformation()
        coords = [[0, 0, 0], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]]
        lattice = [
            [3.8401979337, 0, 0],
            [1.9200989668, 3.3257101909, 0],
            [0, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        conventional_struct = trafo.apply_transformation(struct)
        assert conventional_struct.lattice.alpha == 90


class TestPerturbStructureTransformation:
    def test_apply_transformation(self):
        trafo = PerturbStructureTransformation(0.05)
        coords = [
            [0, 0, 0],
            [0.375, 0.375, 0.375],
            [0.5, 0.5, 0.5],
            [0.875, 0.875, 0.875],
            [0.125, 0.125, 0.125],
            [0.25, 0.25, 0.25],
            [0.625, 0.625, 0.625],
            [0.75, 0.75, 0.75],
        ]

        lattice = [
            [3.8401979337, 0, 0],
            [1.9200989668, 3.3257101909, 0],
            [0, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-", "O2-", "O2-", "O2-"], coords)
        transformed_struct = trafo.apply_transformation(struct)
        for idx, site in enumerate(transformed_struct):
            assert site.distance(struct[idx]) == approx(0.05)

        dct = trafo.as_dict()
        assert isinstance(PerturbStructureTransformation.from_dict(dct), PerturbStructureTransformation)

        t2 = PerturbStructureTransformation(0.05, 0)
        transformed_s2 = t2.apply_transformation(struct)
        for idx, site in enumerate(transformed_s2):
            assert site.distance(struct[idx]) <= 0.05
            assert site.distance(struct[idx]) >= 0

        dct = t2.as_dict()
        assert isinstance(PerturbStructureTransformation.from_dict(dct), PerturbStructureTransformation)


class TestDeformStructureTransformation:
    def test_apply_transformation(self):
        trafo = DeformStructureTransformation([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.05, 1.0]])
        coords = [
            [0, 0, 0],
            [0.375, 0.375, 0.375],
            [0.5, 0.5, 0.5],
            [0.875, 0.875, 0.875],
            [0.125, 0.125, 0.125],
            [0.25, 0.25, 0.25],
            [0.625, 0.625, 0.625],
            [0.75, 0.75, 0.75],
        ]

        lattice = [
            [3.8401979337, 0, 0],
            [1.9200989668, 3.3257101909, 0],
            [0, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-", "O2-", "O2-", "O2-"], coords)
        transformed_s = trafo.apply_transformation(struct)
        assert transformed_s.lattice.a == approx(3.84019793)
        assert transformed_s.lattice.b == approx(3.84379750)
        assert transformed_s.lattice.c == approx(3.75022981)

        dct = json.loads(json.dumps(trafo.as_dict()))
        assert isinstance(DeformStructureTransformation.from_dict(dct), DeformStructureTransformation)


class TestDiscretizeOccupanciesTransformation:
    def test_apply_transformation(self):
        lattice = Lattice.cubic(4)
        struct_orig = Structure(
            lattice,
            [{"Li": 0.19, "Na": 0.19, "K": 0.62}, {"O": 1}],
            [[0, 0, 0], [0.5, 0.5, 0.5]],
        )
        dot = DiscretizeOccupanciesTransformation(max_denominator=5, tol=0.5)
        struct = dot.apply_transformation(struct_orig)
        assert dict(struct[0].species) == {Element("Li"): 0.2, Element("Na"): 0.2, Element("K"): 0.6}

        dot = DiscretizeOccupanciesTransformation(max_denominator=5, tol=0.01)
        with pytest.raises(RuntimeError, match="Cannot discretize structure within tolerance!"):
            dot.apply_transformation(struct_orig)

        struct_orig_2 = Structure(
            lattice,
            [{"Li": 0.5, "Na": 0.25, "K": 0.25}, {"O": 1}],
            [[0, 0, 0], [0.5, 0.5, 0.5]],
        )

        dot = DiscretizeOccupanciesTransformation(max_denominator=9, tol=0.25, fix_denominator=False)

        struct = dot.apply_transformation(struct_orig_2)
        assert dict(struct[0].species) == {
            Element("Li"): Fraction(1 / 2),
            Element("Na"): Fraction(1 / 4),
            Element("K"): Fraction(1 / 4),
        }

        dot = DiscretizeOccupanciesTransformation(max_denominator=9, tol=0.05, fix_denominator=True)
        with pytest.raises(RuntimeError, match="Cannot discretize structure within tolerance"):
            dot.apply_transformation(struct_orig_2)


class TestChargedCellTransformation:
    def test_apply_transformation(self):
        struct_orig = Structure(
            np.eye(3) * 4,
            [{"Li": 0.19, "Na": 0.19, "K": 0.62}, {"O": 1}],
            [[0, 0, 0], [0.5, 0.5, 0.5]],
        )
        cct = ChargedCellTransformation(charge=3)
        struct_trafo = cct.apply_transformation(struct_orig)
        assert struct_trafo.charge == 3


class TestScaleToRelaxedTransformation:
    def test_apply_transformation(self):
        # Test on slab relaxation where volume is fixed
        surf_dir = f"{TEST_FILES_DIR}/surfaces"
        Cu_fin = Structure.from_file(f"{surf_dir}/Cu_slab_fin.cif")
        Cu_init = Structure.from_file(f"{surf_dir}/Cu_slab_init.cif")
        slab_scaling = ScaleToRelaxedTransformation(Cu_init, Cu_fin)
        Au_init = Structure.from_file(f"{surf_dir}/Au_slab_init.cif")
        Au_fin = slab_scaling.apply_transformation(Au_init)
        assert Au_fin.volume == approx(Au_init.volume)

        # Test on gb relaxation
        gb_dir = f"{TEST_FILES_DIR}/core/grain_boundary"
        Be_fin = Structure.from_file(f"{gb_dir}/Be_gb_fin.cif")
        Be_init = Structure.from_file(f"{gb_dir}/Be_gb_init.cif")
        Zn_init = Structure.from_file(f"{gb_dir}/Zn_gb_init.cif")
        gb_scaling = ScaleToRelaxedTransformation(Be_init, Be_fin)
        Zn_fin = gb_scaling.apply_transformation(Zn_init)
        assert all(site.species_string == "Zn" for site in Zn_fin)
        assert (Be_init.lattice.a < Be_fin.lattice.a) == (Zn_init.lattice.a < Zn_fin.lattice.a)
        assert (Be_init.lattice.b < Be_fin.lattice.b) == (Zn_init.lattice.b < Zn_fin.lattice.b)
        assert (Be_init.lattice.c < Be_fin.lattice.c) == (Zn_init.lattice.c < Zn_fin.lattice.c)
        Fe_fin = Structure.from_file(f"{gb_dir}/Fe_gb_fin.cif")
        Fe_init = Structure.from_file(f"{gb_dir}/Fe_gb_init.cif")
        Mo_init = Structure.from_file(f"{gb_dir}/Mo_gb_init.cif")
        gb_scaling = ScaleToRelaxedTransformation(Fe_init, Fe_fin)
        Mo_fin = gb_scaling.apply_transformation(Mo_init)
        assert all(site.species_string == "Mo" for site in Mo_fin)
        assert (Fe_init.lattice.a < Fe_fin.lattice.a) == (Mo_init.lattice.a < Mo_fin.lattice.a)
        assert (Fe_init.lattice.b < Fe_fin.lattice.b) == (Mo_init.lattice.b < Mo_fin.lattice.b)
        assert (Fe_init.lattice.c < Fe_fin.lattice.c) == (Mo_init.lattice.c < Mo_fin.lattice.c)
