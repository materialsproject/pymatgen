# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import annotations

import functools
import json
import os
import random
import unittest
import warnings
from shutil import which

import pytest
from monty.json import MontyDecoder
from pytest import approx

from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.core.sites import PeriodicSite
from pymatgen.io.vasp.inputs import Poscar
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
from pymatgen.util.testing import PymatgenTest

enumlib_present = which("enum.x") and which("makestr.x")


class RotationTransformationsTest(unittest.TestCase):
    def setUp(self):
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ]
        )
        self.struct = Structure(lattice, ["Si"] * 2, coords)

    def test_as_from_dict(self):
        t = RotationTransformation([0, 1, 0], 30, False)
        d = t.as_dict()
        assert isinstance(RotationTransformation.from_dict(d), RotationTransformation)

    def test_rotation_transformation(self):
        t = RotationTransformation([0, 1, 0], 30, False)
        s2 = t.apply_transformation(self.struct)
        s1 = t.inverse.apply_transformation(s2)
        assert (abs(s1.lattice.matrix - self.struct.lattice.matrix) < 1e-8).all()


class RemoveSpeciesTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        t = RemoveSpeciesTransformation(["Li+"])
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ]
        )
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        s = t.apply_transformation(struct)
        assert s.composition.formula == "O2"

        d = t.as_dict()
        assert isinstance(RemoveSpeciesTransformation.from_dict(d), RemoveSpeciesTransformation)


class SubstitutionTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        t = SubstitutionTransformation({"Li+": "Na+", "O2-": "S2-"})
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ]
        )
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        s = t.apply_transformation(struct)
        assert s.composition.formula == "Na2 S2"

    def test_fractional_substitution(self):
        t = SubstitutionTransformation({"Li+": "Na+", "O2-": {"S2-": 0.5, "Se2-": 0.5}})
        # test the to and from dict on the nested dictionary
        t = SubstitutionTransformation.from_dict(t.as_dict())
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ]
        )
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        s = t.apply_transformation(struct)
        assert s.composition.formula == "Na2 Se1 S1"


class SupercellTransformationTest(unittest.TestCase):
    def setUp(self):
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ]
        )
        self.struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)

    def test_apply_transformation(self):
        t = SupercellTransformation([[2, 1, 0], [0, 2, 0], [1, 0, 2]])
        s = t.apply_transformation(self.struct)
        assert s.composition.formula == "Li16 O16"

    def test_from_scaling_factors(self):
        scale_factors = [random.randint(1, 5) for i in range(3)]
        t = SupercellTransformation.from_scaling_factors(*scale_factors)
        s = t.apply_transformation(self.struct)
        assert s.num_sites == 4 * functools.reduce(lambda a, b: a * b, scale_factors)


class OxidationStateDecorationTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        t = OxidationStateDecorationTransformation({"Li": 1, "O": -2})
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ]
        )
        struct = Structure(lattice, ["Li", "Li", "O", "O"], coords)
        s = t.apply_transformation(struct)
        assert s[0].species_string == "Li+"
        assert s[2].species_string == "O2-"
        d = t.as_dict()
        assert isinstance(OxidationStateDecorationTransformation.from_dict(d), OxidationStateDecorationTransformation)


class AutoOxiStateDecorationTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR.LiFePO4"), check_for_POTCAR=False)
        t = AutoOxiStateDecorationTransformation()
        s = t.apply_transformation(p.structure)
        expected_oxi = {"Li": 1, "P": 5, "O": -2, "Fe": 2}
        for site in s:
            assert site.specie.oxi_state == expected_oxi[site.specie.symbol]

    def test_as_from_dict(self):
        t = AutoOxiStateDecorationTransformation()
        d = t.as_dict()
        t = AutoOxiStateDecorationTransformation.from_dict(d)
        assert t.analyzer.dist_scale_factor == 1.015


class OxidationStateRemovalTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        t = OxidationStateRemovalTransformation()
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ]
        )
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        s = t.apply_transformation(struct)
        assert s[0].species_string == "Li"
        assert s[2].species_string == "O"

        d = t.as_dict()
        assert isinstance(OxidationStateRemovalTransformation.from_dict(d), OxidationStateRemovalTransformation)


@unittest.skipIf(not enumlib_present, "enum_lib not present.")
class PartialRemoveSpecieTransformationTest(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_apply_transformation(self):
        t = PartialRemoveSpecieTransformation("Li+", 1.0 / 3, 3)
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ]
        )
        struct = Structure(lattice, ["Li+", "Li+", "Li+", "O2-"], coords)
        assert len(t.apply_transformation(struct, 100)) == 2

        d = t.as_dict()
        assert isinstance(PartialRemoveSpecieTransformation.from_dict(d), PartialRemoveSpecieTransformation)

    def test_apply_transformation_fast(self):
        t = PartialRemoveSpecieTransformation("Li+", 0.5)
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        coords.append([0.1, 0.1, 0.1])
        coords.append([0.3, 0.75, 0.3])
        lattice = Lattice([[10, 0.00, 0.00], [0, 10, 0.00], [0.00, 0, 10]])
        struct = Structure(lattice, ["Li+"] * 6, coords)
        fast_opt_s = t.apply_transformation(struct)
        t = PartialRemoveSpecieTransformation("Li+", 0.5, PartialRemoveSpecieTransformation.ALGO_COMPLETE)
        slow_opt_s = t.apply_transformation(struct)
        assert EwaldSummation(fast_opt_s).total_energy == approx(EwaldSummation(slow_opt_s).total_energy, abs=1e-4)
        assert fast_opt_s == slow_opt_s

    def test_apply_transformations_complete_ranking(self):
        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR.LiFePO4"), check_for_POTCAR=False)
        t1 = OxidationStateDecorationTransformation({"Li": 1, "Fe": 2, "P": 5, "O": -2})
        s = t1.apply_transformation(p.structure)
        t = PartialRemoveSpecieTransformation("Li+", 0.5, PartialRemoveSpecieTransformation.ALGO_COMPLETE)
        assert len(t.apply_transformation(s, 10)) == 6

    def test_apply_transformations_best_first(self):
        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR.LiFePO4"), check_for_POTCAR=False)
        t1 = OxidationStateDecorationTransformation({"Li": 1, "Fe": 2, "P": 5, "O": -2})
        s = t1.apply_transformation(p.structure)
        t = PartialRemoveSpecieTransformation("Li+", 0.5, PartialRemoveSpecieTransformation.ALGO_BEST_FIRST)
        assert len(t.apply_transformation(s)) == 26


class OrderDisorderedStructureTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        t = OrderDisorderedStructureTransformation()
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ]
        )

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
        output = t.apply_transformation(struct, return_ranked_list=50)
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
        output = t.apply_transformation(struct, return_ranked_list=50)
        assert isinstance(output, list)
        assert len(output) == 4
        assert t.lowest_energy_structure == output[0]["structure"]

        struct = Structure(lattice, [{"Si4+": 0.5}, {"Si4+": 0.5}, {"O2-": 0.5}, {"O2-": 0.5}], coords)
        allstructs = t.apply_transformation(struct, 50)
        assert len(allstructs) == 4

        struct = Structure(lattice, [{"Si4+": 0.333}, {"Si4+": 0.333}, {"Si4+": 0.333}, "O2-"], coords)
        allstructs = t.apply_transformation(struct, 50)
        assert len(allstructs) == 3

        d = t.as_dict()
        assert isinstance(OrderDisorderedStructureTransformation.from_dict(d), OrderDisorderedStructureTransformation)

    def test_no_oxidation(self):
        specie = {"Cu1+": 0.5, "Au2+": 0.5}
        cuau = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3.677), [specie], [[0, 0, 0]])
        trans = OrderDisorderedStructureTransformation()
        ss = trans.apply_transformation(cuau, return_ranked_list=100)
        assert ss[0]["structure"].composition["Cu+"] == 2
        trans = OrderDisorderedStructureTransformation(no_oxi_states=True)
        ss = trans.apply_transformation(cuau, return_ranked_list=100)
        assert ss[0]["structure"].composition["Cu+"] == 0
        assert ss[0]["structure"].composition["Cu"] == 2

    def test_symmetrized_structure(self):
        t = OrderDisorderedStructureTransformation(symmetrized_structures=True)
        c = []
        sp = []
        c.append([0.5, 0.5, 0.5])
        sp.append("Si4+")
        c.append([0.45, 0.45, 0.45])
        sp.append({"Si4+": 0.5})
        c.append([0.56, 0.56, 0.56])
        sp.append({"Si4+": 0.5})
        c.append([0.25, 0.75, 0.75])
        sp.append({"Si4+": 0.5})
        c.append([0.75, 0.25, 0.25])
        sp.append({"Si4+": 0.5})
        latt = Lattice.cubic(5)
        struct = Structure(latt, sp, c)
        test_site = PeriodicSite("Si4+", c[2], latt)
        struct = SymmetrizedStructure(struct, "not_real", [0, 1, 1, 2, 2], ["a", "b", "b", "c", "c"])
        output = t.apply_transformation(struct)
        assert test_site in output.sites

    def test_too_small_cell(self):
        t = OrderDisorderedStructureTransformation()
        coords = []
        coords.append([0.5, 0.5, 0.5])
        lattice = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ]
        )
        struct = Structure(lattice, [{"X4+": 0.33, "O2-": 0.33, "P5+": 0.33}], coords)
        with pytest.raises(ValueError):
            t.apply_transformation(struct)

    def test_best_first(self):
        t = OrderDisorderedStructureTransformation(algo=2)
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ]
        )

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
        output = t.apply_transformation(struct, return_ranked_list=3)
        assert output[0]["energy"] == approx(-234.57813667648315, abs=1e-4)


class PrimitiveCellTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        t = PrimitiveCellTransformation()
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.375, 0.375, 0.375])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.875, 0.875, 0.875])
        coords.append([0.125, 0.125, 0.125])
        coords.append([0.25, 0.25, 0.25])
        coords.append([0.625, 0.625, 0.625])
        coords.append([0.75, 0.75, 0.75])

        lattice = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ]
        )
        struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-", "O2-", "O2-", "O2-"], coords)
        s = t.apply_transformation(struct)
        assert len(s) == 4

        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "TiO2_super.json")) as f:
            s = json.load(f, cls=MontyDecoder)
            prim = t.apply_transformation(s)
            assert prim.formula == "Ti4 O8"

        d = t.as_dict()
        assert isinstance(PrimitiveCellTransformation.from_dict(d), PrimitiveCellTransformation)


class ConventionalCellTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        t = ConventionalCellTransformation()
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ]
        )
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        conventional_struct = t.apply_transformation(struct)
        assert conventional_struct.lattice.alpha == 90


class PerturbStructureTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        t = PerturbStructureTransformation(0.05)
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.375, 0.375, 0.375])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.875, 0.875, 0.875])
        coords.append([0.125, 0.125, 0.125])
        coords.append([0.25, 0.25, 0.25])
        coords.append([0.625, 0.625, 0.625])
        coords.append([0.75, 0.75, 0.75])

        lattice = [
            [3.8401979337, 0.00, 0.00],
            [1.9200989668, 3.3257101909, 0.00],
            [0.00, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-", "O2-", "O2-", "O2-"], coords)
        transformed_s = t.apply_transformation(struct)
        for i, site in enumerate(transformed_s):
            assert site.distance(struct[i]) == approx(0.05)

        d = t.as_dict()
        assert isinstance(PerturbStructureTransformation.from_dict(d), PerturbStructureTransformation)

        t2 = PerturbStructureTransformation(0.05, 0)
        transformed_s2 = t2.apply_transformation(struct)
        for i, site in enumerate(transformed_s2):
            assert site.distance(struct[i]) <= 0.05
            assert site.distance(struct[i]) >= 0

        d = t2.as_dict()
        assert isinstance(PerturbStructureTransformation.from_dict(d), PerturbStructureTransformation)


class DeformStructureTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        t = DeformStructureTransformation([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.05, 1.0]])
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.375, 0.375, 0.375])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.875, 0.875, 0.875])
        coords.append([0.125, 0.125, 0.125])
        coords.append([0.25, 0.25, 0.25])
        coords.append([0.625, 0.625, 0.625])
        coords.append([0.75, 0.75, 0.75])

        lattice = [
            [3.8401979337, 0.00, 0.00],
            [1.9200989668, 3.3257101909, 0.00],
            [0.00, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-", "O2-", "O2-", "O2-"], coords)
        transformed_s = t.apply_transformation(struct)
        assert transformed_s.lattice.a == approx(3.84019793)
        assert transformed_s.lattice.b == approx(3.84379750)
        assert transformed_s.lattice.c == approx(3.75022981)

        d = json.loads(json.dumps(t.as_dict()))
        assert isinstance(DeformStructureTransformation.from_dict(d), DeformStructureTransformation)


class DiscretizeOccupanciesTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        latt = Lattice.cubic(4)
        s_orig = Structure(
            latt,
            [{"Li": 0.19, "Na": 0.19, "K": 0.62}, {"O": 1}],
            [[0, 0, 0], [0.5, 0.5, 0.5]],
        )
        dot = DiscretizeOccupanciesTransformation(max_denominator=5, tol=0.5)
        s = dot.apply_transformation(s_orig)
        assert dict(s[0].species) == {Element("Li"): 0.2, Element("Na"): 0.2, Element("K"): 0.6}

        dot = DiscretizeOccupanciesTransformation(max_denominator=5, tol=0.01)
        with pytest.raises(RuntimeError):
            dot.apply_transformation(s_orig)

        s_orig_2 = Structure(
            latt,
            [{"Li": 0.5, "Na": 0.25, "K": 0.25}, {"O": 1}],
            [[0, 0, 0], [0.5, 0.5, 0.5]],
        )

        dot = DiscretizeOccupanciesTransformation(max_denominator=9, tol=0.25, fix_denominator=False)

        s = dot.apply_transformation(s_orig_2)
        assert dict(s[0].species) == {
            Element("Li"): Fraction(1 / 2),
            Element("Na"): Fraction(1 / 4),
            Element("K"): Fraction(1 / 4),
        }

        dot = DiscretizeOccupanciesTransformation(max_denominator=9, tol=0.05, fix_denominator=True)
        with pytest.raises(RuntimeError):
            dot.apply_transformation(s_orig_2)


class ChargedCellTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        lattice = Lattice.cubic(4)
        s_orig = Structure(
            lattice,
            [{"Li": 0.19, "Na": 0.19, "K": 0.62}, {"O": 1}],
            [[0, 0, 0], [0.5, 0.5, 0.5]],
        )
        cct = ChargedCellTransformation(charge=3)
        s = cct.apply_transformation(s_orig)
        assert s.charge == 3


class ScaleToRelaxedTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        # Test on slab relaxation where volume is fixed
        f = os.path.join(PymatgenTest.TEST_FILES_DIR, "surface_tests")
        Cu_fin = Structure.from_file(os.path.join(f, "Cu_slab_fin.cif"))
        Cu_init = Structure.from_file(os.path.join(f, "Cu_slab_init.cif"))
        slab_scaling = ScaleToRelaxedTransformation(Cu_init, Cu_fin)
        Au_init = Structure.from_file(os.path.join(f, "Au_slab_init.cif"))
        Au_fin = slab_scaling.apply_transformation(Au_init)
        assert Au_fin.lattice.volume == approx(Au_init.lattice.volume)

        # Test on gb relaxation
        f = os.path.join(PymatgenTest.TEST_FILES_DIR, "grain_boundary")
        Be_fin = Structure.from_file(os.path.join(f, "Be_gb_fin.cif"))
        Be_init = Structure.from_file(os.path.join(f, "Be_gb_init.cif"))
        Zn_init = Structure.from_file(os.path.join(f, "Zn_gb_init.cif"))
        gb_scaling = ScaleToRelaxedTransformation(Be_init, Be_fin)
        Zn_fin = gb_scaling.apply_transformation(Zn_init)
        assert all(site.species_string == "Zn" for site in Zn_fin)
        assert (Be_init.lattice.a < Be_fin.lattice.a) == (Zn_init.lattice.a < Zn_fin.lattice.a)
        assert (Be_init.lattice.b < Be_fin.lattice.b) == (Zn_init.lattice.b < Zn_fin.lattice.b)
        assert (Be_init.lattice.c < Be_fin.lattice.c) == (Zn_init.lattice.c < Zn_fin.lattice.c)
        Fe_fin = Structure.from_file(os.path.join(f, "Fe_gb_fin.cif"))
        Fe_init = Structure.from_file(os.path.join(f, "Fe_gb_init.cif"))
        Mo_init = Structure.from_file(os.path.join(f, "Mo_gb_init.cif"))
        gb_scaling = ScaleToRelaxedTransformation(Fe_init, Fe_fin)
        Mo_fin = gb_scaling.apply_transformation(Mo_init)
        assert all(site.species_string == "Mo" for site in Mo_fin)
        assert (Fe_init.lattice.a < Fe_fin.lattice.a) == (Mo_init.lattice.a < Mo_fin.lattice.a)
        assert (Fe_init.lattice.b < Fe_fin.lattice.b) == (Mo_init.lattice.b < Mo_fin.lattice.b)
        assert (Fe_init.lattice.c < Fe_fin.lattice.c) == (Mo_init.lattice.c < Mo_fin.lattice.c)


if __name__ == "__main__":
    unittest.main()
