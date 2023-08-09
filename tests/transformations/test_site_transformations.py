from __future__ import annotations

import unittest
from shutil import which

import numpy as np
import pytest

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Molecule, Structure
from pymatgen.transformations.site_transformations import (
    AddSitePropertyTransformation,
    InsertSitesTransformation,
    PartialRemoveSitesTransformation,
    RadialSiteDistortionTransformation,
    RemoveSitesTransformation,
    ReplaceSiteSpeciesTransformation,
    TranslateSitesTransformation,
)
from pymatgen.util.testing import PymatgenTest

enum_cmd = which("enum.x") or which("multienum.x")
makestr_cmd = which("makestr.x") or which("makeStr.x") or which("makeStr.py")
enumlib_present = enum_cmd and makestr_cmd


class TestTranslateSitesTransformation(PymatgenTest):
    def setUp(self):
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
        self.struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-", "O2-", "O2-", "O2-"], coords)

    def test_apply_transformation(self):
        trafo = TranslateSitesTransformation([0, 1], [0.1, 0.2, 0.3])
        s = trafo.apply_transformation(self.struct)
        assert np.allclose(s[0].frac_coords, [0.1, 0.2, 0.3])
        assert np.allclose(s[1].frac_coords, [0.475, 0.575, 0.675])
        inv_t = trafo.inverse
        s = inv_t.apply_transformation(s)
        assert s[0].distance_and_image_from_frac_coords([0, 0, 0])[0] == 0
        assert np.allclose(s[1].frac_coords, [0.375, 0.375, 0.375])

    def test_apply_transformation_site_by_site(self):
        trafo = TranslateSitesTransformation([0, 1], [[0.1, 0.2, 0.3], [-0.075, -0.075, -0.075]])
        s = trafo.apply_transformation(self.struct)
        assert np.allclose(s[0].frac_coords, [0.1, 0.2, 0.3])
        assert np.allclose(s[1].frac_coords, [0.3, 0.3, 0.3])
        inv_t = trafo.inverse
        s = inv_t.apply_transformation(s)
        assert s[0].distance_and_image_from_frac_coords([0, 0, 0])[0] == 0
        assert np.allclose(s[1].frac_coords, [0.375, 0.375, 0.375])

    def test_to_from_dict(self):
        d1 = TranslateSitesTransformation([0], [0.1, 0.2, 0.3]).as_dict()
        d2 = TranslateSitesTransformation([0, 1], [[0.1, 0.2, 0.3], [-0.075, -0.075, -0.075]]).as_dict()
        t1 = TranslateSitesTransformation.from_dict(d1)
        t2 = TranslateSitesTransformation.from_dict(d2)
        s1 = t1.apply_transformation(self.struct)
        s2 = t2.apply_transformation(self.struct)
        assert np.allclose(s1[0].frac_coords, [0.1, 0.2, 0.3])
        assert np.allclose(s2[0].frac_coords, [0.1, 0.2, 0.3])
        assert np.allclose(s2[1].frac_coords, [0.3, 0.3, 0.3])
        str(t1)
        str(t2)


class TestReplaceSiteSpeciesTransformation(unittest.TestCase):
    def setUp(self):
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
        self.struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-", "O2-", "O2-", "O2-"], coords)

    def test_apply_transformation(self):
        trafo = ReplaceSiteSpeciesTransformation({0: "Na"})
        s = trafo.apply_transformation(self.struct)
        assert s.formula == "Na1 Li3 O4"

    def test_to_from_dict(self):
        d = ReplaceSiteSpeciesTransformation({0: "Na"}).as_dict()
        trafo = ReplaceSiteSpeciesTransformation.from_dict(d)
        s = trafo.apply_transformation(self.struct)
        assert s.formula == "Na1 Li3 O4"


class TestRemoveSitesTransformation(unittest.TestCase):
    def setUp(self):
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
        self.struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-", "O2-", "O2-", "O2-"], coords)

    def test_apply_transformation(self):
        trafo = RemoveSitesTransformation(range(2))
        s = trafo.apply_transformation(self.struct)
        assert s.formula == "Li2 O4"

    def test_to_from_dict(self):
        d = RemoveSitesTransformation(range(2)).as_dict()
        trafo = RemoveSitesTransformation.from_dict(d)
        s = trafo.apply_transformation(self.struct)
        assert s.formula == "Li2 O4"


class TestInsertSitesTransformation(unittest.TestCase):
    def setUp(self):
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
        self.struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-", "O2-", "O2-", "O2-"], coords)

    def test_apply_transformation(self):
        trafo = InsertSitesTransformation(["Fe", "Mn"], [[0.0, 0.5, 0], [0.5, 0.2, 0.2]])
        s = trafo.apply_transformation(self.struct)
        assert s.formula == "Li4 Mn1 Fe1 O4"
        trafo = InsertSitesTransformation(["Fe", "Mn"], [[0.001, 0, 0], [0.1, 0.2, 0.2]])

        # Test validate proximity
        with pytest.raises(ValueError, match="New site is too close to an existing site!"):
            trafo.apply_transformation(self.struct)

    def test_to_from_dict(self):
        d = InsertSitesTransformation(["Fe", "Mn"], [[0.5, 0, 0], [0.1, 0.5, 0.2]]).as_dict()
        trafo = InsertSitesTransformation.from_dict(d)
        s = trafo.apply_transformation(self.struct)
        assert s.formula == "Li4 Mn1 Fe1 O4"


class TestPartialRemoveSitesTransformation(unittest.TestCase):
    def setUp(self):
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
        self.struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-", "O2-", "O2-", "O2-"], coords)

    def test_apply_transformation_complete(self):
        trafo = PartialRemoveSitesTransformation(
            [tuple(range(4)), tuple(range(4, 8))],
            [0.5, 0.5],
            PartialRemoveSitesTransformation.ALGO_COMPLETE,
        )
        s = trafo.apply_transformation(self.struct)
        assert s.formula == "Li2 O2"
        s = trafo.apply_transformation(self.struct, 12)
        assert len(s) == 12

    @unittest.skipIf(not enumlib_present, "enum_lib not present.")
    def test_apply_transformation_enumerate(self):
        trafo = PartialRemoveSitesTransformation(
            [tuple(range(4)), tuple(range(4, 8))],
            [0.5, 0.5],
            PartialRemoveSitesTransformation.ALGO_ENUMERATE,
        )
        s = trafo.apply_transformation(self.struct)
        assert s.formula == "Li2 O2"
        s = trafo.apply_transformation(self.struct, 12)
        assert len(s) == 12

    def test_apply_transformation_best_first(self):
        trafo = PartialRemoveSitesTransformation(
            [tuple(range(4)), tuple(range(4, 8))],
            [0.5, 0.5],
            PartialRemoveSitesTransformation.ALGO_BEST_FIRST,
        )
        s = trafo.apply_transformation(self.struct)
        assert s.formula == "Li2 O2"

    def test_apply_transformation_fast(self):
        trafo = PartialRemoveSitesTransformation(
            [tuple(range(4)), tuple(range(4, 8))],
            [0.5, 0.5],
            PartialRemoveSitesTransformation.ALGO_FAST,
        )
        s = trafo.apply_transformation(self.struct)
        assert s.formula == "Li2 O2"
        trafo = PartialRemoveSitesTransformation([tuple(range(8))], [0.5], PartialRemoveSitesTransformation.ALGO_FAST)
        s = trafo.apply_transformation(self.struct)
        assert s.formula == "Li2 O2"

    def test_to_from_dict(self):
        d = PartialRemoveSitesTransformation([tuple(range(4))], [0.5]).as_dict()
        trafo = PartialRemoveSitesTransformation.from_dict(d)
        s = trafo.apply_transformation(self.struct)
        assert s.formula == "Li2 O4"

    def test_str(self):
        d = PartialRemoveSitesTransformation([tuple(range(4))], [0.5]).as_dict()
        assert str(d) is not None


class TestAddSitePropertyTransformation(PymatgenTest):
    def test_apply_transformation(self):
        struct = self.get_structure("Li2O2")
        sd = [[True, True, True] for _ in struct]
        bader = np.random.random(len(struct)).tolist()
        site_props = {"selective_dynamics": sd, "bader": bader}
        trans = AddSitePropertyTransformation(site_props)
        manually_set = struct.copy()
        for prop, value in site_props.items():
            manually_set.add_site_property(prop, value)
        trans_set = trans.apply_transformation(struct)
        for prop in site_props:
            assert np.allclose(trans_set.site_properties[prop], manually_set.site_properties[prop])


class TestRadialSiteDistortionTransformation(PymatgenTest):
    def setUp(self):
        self.molecule = Molecule(
            species=["C", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H"],
            coords=[
                [0, 0, 0],
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
                [-1, 0, 0],
                [0, -1, 0],
                [0, 0, -1],
                [3, 0, 0],
                [0, 3, 0],
                [0, 0, 3],
                [-3, 0, 0],
                [0, -3, 0],
                [0, 0, -3],
            ],
        )
        self.structure = Structure(
            species=["C", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H"],
            coords=[
                [0, 0, 0],
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
                [-1, 0, 0],
                [0, -1, 0],
                [0, 0, -1],
                [3, 0, 0],
                [0, 3, 0],
                [0, 0, 3],
                [-3, 0, 0],
                [0, -3, 0],
                [0, 0, -3],
            ],
            lattice=[[10, 0, 0], [0, 10, 0], [0, 0, 10]],
            coords_are_cartesian=True,
        )

    def test(self):
        trafo = RadialSiteDistortionTransformation(0, 1, nn_only=True)
        struct = trafo.apply_transformation(self.molecule)
        assert np.array_equal(struct[0].coords, [0, 0, 0])
        assert np.array_equal(struct[1].coords, [2, 0, 0])
        assert np.array_equal(struct[2].coords, [0, 2, 0])
        assert np.array_equal(struct[3].coords, [0, 0, 2])
        assert np.array_equal(struct[4].coords, [-2, 0, 0])
        assert np.array_equal(struct[5].coords, [0, -2, 0])
        assert np.array_equal(struct[6].coords, [0, 0, -2])

        trafo = RadialSiteDistortionTransformation(0, 1, nn_only=True)
        struct = trafo.apply_transformation(self.structure)
        for c1, c2 in zip(self.structure[1:7], struct[1:7]):
            assert c1.distance(c2) == 1.0

        assert np.array_equal(struct[0].coords, [0, 0, 0])
        assert np.array_equal(struct[1].coords, [2, 0, 0])
        assert np.array_equal(struct[2].coords, [0, 2, 0])
        assert np.array_equal(struct[3].coords, [0, 0, 2])
        assert np.array_equal(struct[4].coords, [8, 0, 0])
        assert np.array_equal(struct[5].coords, [0, 8, 0])
        assert np.array_equal(struct[6].coords, [0, 0, 8])

    def test_second_nn(self):
        trafo = RadialSiteDistortionTransformation(0, 1, nn_only=False)
        s = trafo.apply_transformation(self.molecule)
        for c1, c2 in zip(self.molecule[7:], s[7:]):
            assert abs(round(sum(c2.coords - c1.coords), 2)) == 0.33
