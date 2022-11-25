# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import json
import os
import random
import unittest
import warnings
from pathlib import Path
from shutil import which

import numpy as np
import pytest
from monty.json import MontyDecoder, MontyEncoder
from monty.tempfile import ScratchDir

from pymatgen.core.composition import Composition
from pymatgen.core.lattice import Lattice
from pymatgen.core.operations import SymmOp
from pymatgen.core.periodic_table import Element, Species
from pymatgen.core.structure import (
    IMolecule,
    IStructure,
    Molecule,
    PeriodicNeighbor,
    Structure,
    StructureError,
)
from pymatgen.electronic_structure.core import Magmom
from pymatgen.util.testing import PymatgenTest

enum_cmd = which("enum.x") or which("multienum.x")
mcsqs_cmd = which("mcsqs")

try:
    import m3gnet
except ImportError:
    m3gnet = None


class NeighborTest(PymatgenTest):
    def test_msonable(self):
        s = PymatgenTest.get_structure("Li2O")
        nn = s.get_neighbors(s[0], r=3)
        self.assertIsInstance(nn[0], PeriodicNeighbor)
        str_ = json.dumps(nn, cls=MontyEncoder)
        nn = json.loads(str_, cls=MontyDecoder)
        self.assertIsInstance(nn[0], PeriodicNeighbor)


class IStructureTest(PymatgenTest):
    def setUp(self):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        self.lattice = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ]
        )
        self.struct = IStructure(self.lattice, ["Si"] * 2, coords)
        assert len(self.struct) == 2, "Wrong number of sites in structure!"
        assert self.struct.is_ordered
        assert self.struct.ntypesp == 1
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.0, 0, 0.0000001])
        with pytest.raises(StructureError):
            IStructure(
                self.lattice,
                ["Si"] * 2,
                coords,
                validate_proximity=True,
            )
        self.propertied_structure = IStructure(self.lattice, ["Si"] * 2, coords, site_properties={"magmom": [5, -5]})

        self.lattice_pbc = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ],
            pbc=(True, True, False),
        )

    @unittest.skipIf(not (mcsqs_cmd and enum_cmd), "enumlib or mcsqs executable not present")
    def test_get_orderings(self):
        ordered = Structure.from_spacegroup("Im-3m", Lattice.cubic(3), ["Fe"], [[0, 0, 0]])
        assert ordered.get_orderings()[0] == ordered
        disordered = Structure.from_spacegroup("Im-3m", Lattice.cubic(3), [Composition("Fe0.5Mn0.5")], [[0, 0, 0]])
        orderings = disordered.get_orderings()
        assert len(orderings) == 1
        super_cell = disordered * 2
        orderings = super_cell.get_orderings()
        assert len(orderings) == 59
        sqs = disordered.get_orderings(mode="sqs", scaling=[2, 2, 2])
        assert sqs[0].formula == "Mn8 Fe8"

        sqs = super_cell.get_orderings(mode="sqs")
        assert sqs[0].formula == "Mn8 Fe8"

    def test_as_dataframe(self):
        df = self.propertied_structure.as_dataframe()
        assert df.attrs["Reduced Formula"] == "Si"
        assert df.shape == (2, 8)

    def test_equal(self):
        struct = self.struct
        assert struct == struct
        assert struct == struct.copy()
        assert not struct == 2 * struct

        assert not struct == "a" * len(struct)  # GH-2584
        assert struct is not None
        assert not struct == 42  # GH-2587

        assert struct == Structure.from_dict(struct.as_dict())

        struct_2 = Structure.from_sites(struct)
        assert struct, struct_2
        struct_2.apply_strain(0.5)
        assert struct != struct_2

    def test_matches(self):
        ss = self.struct * 2
        assert ss.matches(self.struct)

    def test_bad_structure(self):
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        coords.append([0.75, 0.5, 0.75])
        with pytest.raises(StructureError):
            IStructure(
                self.lattice,
                ["Si"] * 3,
                coords,
                validate_proximity=True,
            )
        # these shouldn't raise an error
        IStructure(self.lattice, ["Si"] * 2, coords[:2], True)
        IStructure(self.lattice, ["Si"], coords[:1], True)

    def test_volume_and_density(self):
        assert round(abs(self.struct.volume - 40.04), 2) == 0, "Volume wrong!"
        assert round(abs(self.struct.density - 2.33), 2) == 0, "Incorrect density"

    def test_specie_init(self):
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        s = IStructure(self.lattice, [{Species("O", -2): 1.0}, {Species("Mg", 2): 0.8}], coords)
        assert s.composition.formula == "Mg0.8 O1"

    def test_get_sorted_structure(self):
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        s = IStructure(self.lattice, ["O", "Li"], coords, site_properties={"charge": [-2, 1]})
        sorted_s = s.get_sorted_structure()
        assert sorted_s[0].species == Composition("Li")
        assert sorted_s[1].species == Composition("O")
        assert sorted_s[0].charge == 1
        assert sorted_s[1].charge == -2
        s = IStructure(
            self.lattice,
            ["Se", "C", "Se", "C"],
            [[0] * 3, [0.5] * 3, [0.25] * 3, [0.75] * 3],
        )
        assert [site.specie.symbol for site in s.get_sorted_structure()] == ["C", "C", "Se", "Se"]

    def test_get_space_group_data(self):
        assert self.struct.get_space_group_info() == ("Fd-3m", 227)

    def test_fractional_occupations(self):
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        s = IStructure(self.lattice, [{"O": 1.0}, {"Mg": 0.8}], coords)
        assert s.composition.formula == "Mg0.8 O1"
        assert not s.is_ordered

    def test_get_distance(self):
        assert round(abs(self.struct.get_distance(0, 1) - 2.35), 2) == 0, "Distance calculated wrongly!"
        pt = [0.9, 0.9, 0.8]
        assert (
            round(abs(self.struct[0].distance_from_point(pt) - 1.50332963784), 2) == 0
        ), "Distance calculated wrongly!"

    def test_as_dict(self):
        si = Species("Si", 4)
        mn = Element("Mn")
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        struct = IStructure(self.lattice, [{si: 0.5, mn: 0.5}, {si: 0.5}], coords)
        assert "lattice" in struct.as_dict()
        assert "sites" in struct.as_dict()
        d = self.propertied_structure.as_dict()
        assert d["sites"][0]["properties"]["magmom"] == 5
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        s = IStructure(
            self.lattice,
            [
                {Species("O", -2, properties={"spin": 3}): 1.0},
                {Species("Mg", 2, properties={"spin": 2}): 0.8},
            ],
            coords,
            site_properties={"magmom": [5, -5]},
        )
        d = s.as_dict()
        assert d["sites"][0]["properties"]["magmom"] == 5
        assert d["sites"][0]["species"][0]["properties"]["spin"] == 3

        d = s.as_dict(0)
        assert "volume" not in d["lattice"]
        assert "xyz" not in d["sites"][0]

    def test_from_dict(self):

        d = self.propertied_structure.as_dict()
        s = IStructure.from_dict(d)
        assert s[0].magmom == 5
        d = self.propertied_structure.as_dict(0)
        s2 = IStructure.from_dict(d)
        assert s == s2

        d = {
            "lattice": {
                "a": 3.8401979337,
                "volume": 40.044794644251596,
                "c": 3.8401979337177736,
                "b": 3.840198994344244,
                "matrix": [
                    [3.8401979337, 0.0, 0.0],
                    [1.9200989668, 3.3257101909, 0.0],
                    [0.0, -2.2171384943, 3.1355090603],
                ],
                "alpha": 119.9999908639842,
                "beta": 90.0,
                "gamma": 60.000009137322195,
            },
            "sites": [
                {
                    "properties": {"magmom": 5},
                    "abc": [0.0, 0.0, 0.0],
                    "occu": 1.0,
                    "species": [
                        {
                            "occu": 1.0,
                            "oxidation_state": -2,
                            "properties": {"spin": 3},
                            "element": "O",
                        }
                    ],
                    "label": "O2-",
                    "xyz": [0.0, 0.0, 0.0],
                },
                {
                    "properties": {"magmom": -5},
                    "abc": [0.75, 0.5, 0.75],
                    "occu": 0.8,
                    "species": [
                        {
                            "occu": 0.8,
                            "oxidation_state": 2,
                            "properties": {"spin": 2},
                            "element": "Mg",
                        }
                    ],
                    "label": "Mg2+:0.800",
                    "xyz": [3.8401979336749994, 1.2247250003039056e-06, 2.351631795225],
                },
            ],
        }
        s = IStructure.from_dict(d)
        assert s[0].magmom == 5
        assert s[0].specie.spin == 3
        assert isinstance(s, IStructure)

    def test_site_properties(self):
        site_props = self.propertied_structure.site_properties
        assert site_props["magmom"] == [5, -5]
        assert self.propertied_structure[0].magmom == 5
        assert self.propertied_structure[1].magmom == -5

    def test_copy(self):
        new_struct = self.propertied_structure.copy(site_properties={"charge": [2, 3]})
        assert new_struct[0].magmom == 5
        assert new_struct[1].magmom == -5
        assert new_struct[0].charge == 2
        assert new_struct[1].charge == 3

        coords = []
        coords.append([0, 0, 0])
        coords.append([0.0, 0, 0.0000001])

        structure = IStructure(self.lattice, ["O", "Si"], coords, site_properties={"magmom": [5, -5]})

        new_struct = structure.copy(site_properties={"charge": [2, 3]}, sanitize=True)
        assert new_struct[0].magmom == -5
        assert new_struct[1].magmom == 5
        assert new_struct[0].charge == 3
        assert new_struct[1].charge == 2
        assert round(abs(new_struct.volume - structure.volume), 7) == 0

    def test_interpolate(self):
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        struct = IStructure(self.lattice, ["Si"] * 2, coords)
        coords2 = []
        coords2.append([0, 0, 0])
        coords2.append([0.5, 0.5, 0.5])
        struct2 = IStructure(self.struct.lattice, ["Si"] * 2, coords2)
        int_s = struct.interpolate(struct2, 10)
        for s in int_s:
            assert s is not None, "Interpolation Failed!"
            assert int_s[0].lattice == s.lattice
        self.assertArrayEqual(int_s[1][1].frac_coords, [0.725, 0.5, 0.725])

        # test ximages
        int_s = struct.interpolate(struct2, nimages=np.linspace(0.0, 1.0, 3))
        for s in int_s:
            assert s is not None, "Interpolation Failed!"
            assert int_s[0].lattice == s.lattice
        self.assertArrayEqual(int_s[1][1].frac_coords, [0.625, 0.5, 0.625])

        badlattice = [[1, 0.00, 0.00], [0, 1, 0.00], [0.00, 0, 1]]
        struct2 = IStructure(badlattice, ["Si"] * 2, coords2)
        with pytest.raises(ValueError):
            struct.interpolate(struct2)

        coords2 = []
        coords2.append([0, 0, 0])
        coords2.append([0.5, 0.5, 0.5])
        struct2 = IStructure(self.struct.lattice, ["Si", "Fe"], coords2)
        with pytest.raises(ValueError):
            struct.interpolate(struct2)

        # Test autosort feature.
        s1 = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3), ["Fe"], [[0, 0, 0]])
        s1.pop(0)
        s2 = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3), ["Fe"], [[0, 0, 0]])
        s2.pop(2)
        random.shuffle(s2)

        for s in s1.interpolate(s2, autosort_tol=0.5):
            self.assertArrayAlmostEqual(s1[0].frac_coords, s[0].frac_coords)
            self.assertArrayAlmostEqual(s1[2].frac_coords, s[2].frac_coords)

        # Make sure autosort has no effect on simpler interpolations,
        # and with shuffled sites.
        s1 = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3), ["Fe"], [[0, 0, 0]])
        s2 = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3), ["Fe"], [[0, 0, 0]])
        s2[0] = "Fe", [0.01, 0.01, 0.01]
        random.shuffle(s2)

        for s in s1.interpolate(s2, autosort_tol=0.5):
            self.assertArrayAlmostEqual(s1[1].frac_coords, s[1].frac_coords)
            self.assertArrayAlmostEqual(s1[2].frac_coords, s[2].frac_coords)
            self.assertArrayAlmostEqual(s1[3].frac_coords, s[3].frac_coords)

        # Test non-hexagonal setting.
        lattice = Lattice.rhombohedral(4.0718, 89.459)
        species = [{"S": 1.0}, {"Ni": 1.0}]
        coordinate = [(0.252100, 0.252100, 0.252100), (0.500000, 0.244900, -0.244900)]
        s = Structure.from_spacegroup("R32:R", lattice, species, coordinate)
        assert s.formula == "Ni3 S2"

        # test pbc
        coords = [[0.85, 0.85, 0.85]]
        coords2 = [[0.25, 0.25, 0.25]]
        struct_pbc = IStructure(self.lattice_pbc, ["Si"], coords)
        struct2_pbc = IStructure(self.lattice_pbc, ["Si"], coords2)
        int_s_pbc = struct_pbc.interpolate(struct2_pbc, nimages=2)
        self.assertArrayAlmostEqual(int_s_pbc[1][0].frac_coords, [1.05, 1.05, 0.55])

    def test_interpolate_lattice(self):
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        struct = IStructure(self.lattice, ["Si"] * 2, coords)
        coords2 = []
        coords2.append([0, 0, 0])
        coords2.append([0.5, 0.5, 0.5])
        l2 = Lattice.from_parameters(3, 4, 4, 100, 100, 70)
        struct2 = IStructure(l2, ["Si"] * 2, coords2)
        int_s = struct.interpolate(struct2, 2, interpolate_lattices=True)
        self.assertArrayAlmostEqual(struct.lattice.abc, int_s[0].lattice.abc)
        self.assertArrayAlmostEqual(struct.lattice.angles, int_s[0].lattice.angles)
        self.assertArrayAlmostEqual(struct2.lattice.abc, int_s[2].lattice.abc)
        self.assertArrayAlmostEqual(struct2.lattice.angles, int_s[2].lattice.angles)
        int_angles = [110.3976469, 94.5359731, 64.5165856]
        self.assertArrayAlmostEqual(int_angles, int_s[1].lattice.angles)

        # Assert that volume is monotonic
        assert struct2.lattice.volume >= int_s[1].lattice.volume
        assert int_s[1].lattice.volume >= struct.lattice.volume

    def test_interpolate_lattice_rotation(self):
        l1 = Lattice([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        l2 = Lattice([[-1.01, 0, 0], [0, -1.01, 0], [0, 0, 1]])
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        struct1 = IStructure(l1, ["Si"] * 2, coords)
        struct2 = IStructure(l2, ["Si"] * 2, coords)
        int_s = struct1.interpolate(struct2, 2, interpolate_lattices=True)

        # Assert that volume is monotonic
        assert struct2.lattice.volume >= int_s[1].lattice.volume
        assert int_s[1].lattice.volume >= struct1.lattice.volume

    def test_get_primitive_structure(self):
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]]
        fcc_ag = IStructure(Lattice.cubic(4.09), ["Ag"] * 4, coords)
        assert len(fcc_ag.get_primitive_structure()) == 1
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        bcc_li = IStructure(Lattice.cubic(4.09), ["Li"] * 2, coords)
        bcc_prim = bcc_li.get_primitive_structure()
        assert len(bcc_prim) == 1
        assert round(abs(bcc_prim.lattice.alpha - 109.47122), 3) == 0
        bcc_li = IStructure(Lattice.cubic(4.09), ["Li"] * 2, coords, site_properties={"magmom": [1, -1]})
        bcc_prim = bcc_li.get_primitive_structure()
        assert len(bcc_prim) == 1
        assert round(abs(bcc_prim.lattice.alpha - 109.47122), 3) == 0
        bcc_prim = bcc_li.get_primitive_structure(use_site_props=True)
        assert len(bcc_prim) == 2
        assert round(abs(bcc_prim.lattice.alpha - 90), 3) == 0

        coords = [[0] * 3, [0.5] * 3, [0.25] * 3, [0.26] * 3]
        s = IStructure(Lattice.cubic(4.09), ["Ag"] * 4, coords)
        assert len(s.get_primitive_structure()) == 4

    def test_primitive_cell_site_merging(self):
        l = Lattice.cubic(10)
        coords = [[0, 0, 0], [0, 0, 0.5], [0, 0, 0.26], [0, 0, 0.74]]
        sp = ["Ag", "Ag", "Be", "Be"]
        s = Structure(l, sp, coords)
        dm = s.get_primitive_structure().distance_matrix
        self.assertArrayAlmostEqual(dm, [[0, 2.5], [2.5, 0]])

    def test_primitive_on_large_supercell(self):
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]]
        fcc_ag = Structure(Lattice.cubic(4.09), ["Ag"] * 4, coords)
        fcc_ag.make_supercell([2, 2, 2])
        fcc_ag_prim = fcc_ag.get_primitive_structure()
        assert len(fcc_ag_prim) == 1
        assert round(abs(fcc_ag_prim.volume - 17.10448225), 7) == 0

    def test_primitive_positions(self):
        coords = [[0, 0, 0], [0.3, 0.35, 0.45]]
        s = Structure(Lattice.from_parameters(1, 2, 3, 50, 66, 88), ["Ag"] * 2, coords)

        c = [[2, 0, 0], [1, 3, 0], [1, 1, 1]]

        for sc_matrix in [c]:
            sc = s.copy()
            sc.make_supercell(sc_matrix)
            prim = sc.get_primitive_structure(0.01)

            assert len(prim) == 2
            assert round(abs(prim.distance_matrix[0, 1] - 1.0203432356739286), 7) == 0

    def test_primitive_structure_volume_check(self):
        l = Lattice.tetragonal(10, 30)
        coords = [
            [0.5, 0.8, 0],
            [0.5, 0.2, 0],
            [0.5, 0.8, 0.333],
            [0.5, 0.5, 0.333],
            [0.5, 0.5, 0.666],
            [0.5, 0.2, 0.666],
        ]
        s = IStructure(l, ["Ag"] * 6, coords)
        sprim = s.get_primitive_structure(tolerance=0.1)
        assert len(sprim) == 6

    def test_get_miller_index(self):
        """Test for get miller index convenience method"""
        struct = Structure(
            [2.319, -4.01662582, 0.0, 2.319, 4.01662582, 0.0, 0.0, 0.0, 7.252],
            ["Sn", "Sn", "Sn"],
            [
                [2.319, 1.33887527, 6.3455],
                [1.1595, 0.66943764, 4.5325],
                [1.1595, 0.66943764, 0.9065],
            ],
            coords_are_cartesian=True,
        )
        hkl = struct.get_miller_index_from_site_indexes([0, 1, 2])
        assert hkl == (2, -1, 0)

    def test_get_all_neighbors_and_get_neighbors(self):
        s = self.struct
        nn = s.get_neighbors_in_shell(s[0].frac_coords, 2, 4, include_index=True, include_image=True)
        assert len(nn) == 47
        r = random.uniform(3, 6)
        all_nn = s.get_all_neighbors(r, True, True)
        for idx, site in enumerate(s):
            assert 4 == len(all_nn[idx][0])
            assert len(all_nn[idx]) == len(s.get_neighbors(site, r))

        for site, nns in zip(s, all_nn):
            for nn in nns:
                assert nn[0].is_periodic_image(s[nn[2]])
                d = sum((site.coords - nn[0].coords) ** 2) ** 0.5
                assert round(abs(d - nn[1]), 7) == 0

        s = Structure(Lattice.cubic(1), ["Li"], [[0, 0, 0]])
        s.make_supercell([2, 2, 2])
        assert sum(map(len, s.get_all_neighbors(3))) == 976

        all_nn = s.get_all_neighbors(0.05)
        assert [len(nn) for nn in all_nn] == [0] * len(s)

        # the following test is from issue #2226
        poscar = """POSCAR
 1.0000000000000000
     6.9082208665474800    0.0000000000000005    0.0000000000000011
    -0.0000000000000008   14.1745610988433715    0.0000000000000004
     0.0000000000000005   -0.0000000000000019   20.0189293405157045
C
 14
Direct
  1.2615805559557376  1.4846778919841908  0.3162565598606580
 -0.7615805559557360 -0.9846778919841882 -0.1837434401393425
  1.2615805559557380 -0.4846778919841886 -0.1837434401393425
 -0.7615805559557358  0.9846778919841882  0.3162565598606581
 -0.2615805559557363 -0.4846778919841892  0.6837434401393432
  0.7615805559557360  0.9846778919841881  0.1837434401393430
 -0.2653510291469959  0.5275483828607898  0.1211255106369795
  0.7653510291469972 -0.0275483828607906  0.6211255106369804
  1.2653510291469956  0.5275483828607898  0.3788744893630209
 -0.7653510291469972 -0.0275483828607906 -0.1211255106369797
  1.2653510291469956  0.4724516171392097 -0.1211255106369793
 -0.7653510291469972  0.0275483828607905  0.3788744893630209
 -0.2653510291469959  0.4724516171392097  0.6211255106369801
 -0.2705230397846415  1.4621722452479102  0.0625618775773844
"""
        s = Structure.from_str(poscar, fmt="poscar")
        site0 = s.sites[1]
        site1 = s.sites[9]
        neigh_sites = s.get_neighbors(site0, 2.0)
        assert len(neigh_sites) == 1
        neigh_sites = s.get_neighbors(site1, 2.0)
        assert len(neigh_sites) == 1

    def test_get_neighbor_list(self):
        s = self.struct
        c_indices1, c_indices2, c_offsets, c_distances = s.get_neighbor_list(3)
        p_indices1, p_indices2, p_offsets, p_distances = s._get_neighbor_list_py(3)
        self.assertArrayAlmostEqual(sorted(c_distances), sorted(p_distances))

    # @unittest.skipIf(not os.environ.get("CI"), "Only run this in CI tests.")
    # def test_get_all_neighbors_crosscheck_old(self):
    #     warnings.simplefilter("ignore")
    #     for i in range(100):
    #         alpha, beta = np.random.rand(2) * 90
    #         a, b, c = 3 + np.random.rand(3) * 5
    #         species = ["H"] * 5
    #         frac_coords = np.random.rand(5, 3)
    #         try:
    #             latt = Lattice.from_parameters(a, b, c, alpha, beta, 90)
    #             s = Structure.from_spacegroup("P1", latt,
    #                                           species, frac_coords)
    #             for nn_new, nn_old in zip(s.get_all_neighbors(4),
    #                                       s.get_all_neighbors_old(4)):
    #                 sites1 = [i[0] for i in nn_new]
    #                 sites2 = [i[0] for i in nn_old]
    #                 self.assertEqual(set(sites1), set(sites2))
    #             break
    #         except Exception as ex:
    #             pass
    #     else:
    #         raise ValueError("No valid structure tested.")
    #
    #     from pymatgen.electronic_structure.core import Spin
    #     d = {'@module': 'pymatgen.core.structure', '@class': 'Structure', 'charge': None, 'lattice': {
    #         'matrix': [[0.0, 0.0, 5.5333], [5.7461, 0.0, 3.518471486290303e-16],
    #                    [-4.692662837312786e-16, 7.6637, 4.692662837312786e-16]], 'a': 5.5333, 'b': 5.7461,
    #                    'c': 7.6637,
    #         'alpha': 90.0, 'beta': 90.0, 'gamma': 90.0, 'volume': 243.66653780778103}, 'sites': [
    #         {'species': [{'element': 'Mn', 'oxidation_state': 0, 'properties': {'spin': Spin.down}, 'occu': 1}],
    #          'abc': [0.0, 0.5, 0.5], 'xyz': [2.8730499999999997, 3.83185, 4.1055671618015446e-16],
    #          'label': 'Mn0+,spin=-1',
    #          'properties': {}},
    #         {'species': [{'element': 'Mn', 'oxidation_state': None, 'occu': 1.0}],
    #          'abc': [1.232595164407831e-32, 0.5, 0.5],
    #          'xyz': [2.8730499999999997, 3.83185, 4.105567161801545e-16], 'label': 'Mn', 'properties': {}}]}
    #     struct = Structure.from_dict(d)
    #     self.assertEqual(set([i[0] for i in struct.get_neighbors(struct[0], 0.05)]),
    #                      set([i[0] for i in struct.get_neighbors_old(struct[0], 0.05)]))
    #
    #     warnings.simplefilter("default")

    def test_get_symmetric_neighbor_list(self):
        # tetragonal group with all bonds related by symmetry
        s = Structure.from_spacegroup(100, [[1, 0, 0], [0, 1, 0], [0, 0, 2]], ["Fe"], [[0.0, 0.0, 0.0]])
        c_indices, p_indices, offsets, distances, s_indices, symops = s.get_symmetric_neighbor_list(0.8, sg=100)
        assert len(np.unique(s_indices)) == 1
        assert s_indices[0] == 0
        assert (~np.isnan(s_indices)).all()
        assert (symops[0].affine_matrix == np.eye(4)).all()
        # now more complicated example with bonds of same length but with different symmetry
        s2 = Structure.from_spacegroup(198, [[8.908, 0, 0], [0, 8.908, 0], [0, 0, 8.908]], ["Cu"], [[0.0, 0.0, 0.0]])
        c_indices2, p_indices2, offsets2, distances2, s_indices2, symops2 = s2.get_symmetric_neighbor_list(7, sg=198)
        assert len(np.unique(s_indices2)) == 2
        assert len(s_indices2) == 48
        assert len(s_indices2[s_indices2 == 0]) == len(s_indices2[s_indices2 == 1])
        assert s_indices2[0] == 0
        assert s_indices2[24] == 1
        assert np.isclose(distances2[0], distances2[24])
        assert (symops2[0].affine_matrix == np.eye(4)).all()
        assert (symops2[24].affine_matrix == np.eye(4)).all()
        from_a2 = s2[c_indices2[0]].frac_coords
        to_a2 = s2[p_indices2[0]].frac_coords
        r_a2 = offsets2[0]
        from_b2 = s2[c_indices2[1]].frac_coords
        to_b2 = s2[p_indices2[1]].frac_coords
        r_b2 = offsets2[1]
        assert symops2[1].are_symmetrically_related_vectors(from_a2, to_a2, r_a2, from_b2, to_b2, r_b2)
        assert symops2[1].are_symmetrically_related_vectors(from_b2, to_b2, r_b2, from_a2, to_a2, r_a2)
        c_indices3, p_indices3, offsets3, distances3, s_indices3, symops3 = s2.get_symmetric_neighbor_list(
            7, sg=198, unique=True
        )
        assert (np.sort(np.array([c_indices3, p_indices3]).flatten()) == np.sort(c_indices2)).all()
        assert (np.sort(np.array([c_indices3, p_indices3]).flatten()) == np.sort(p_indices2)).all()

    def test_get_all_neighbors_outside_cell(self):
        s = Structure(
            Lattice.cubic(2),
            ["Li", "Li", "Li", "Si"],
            [[3.1] * 3, [0.11] * 3, [-1.91] * 3, [0.5] * 3],
        )
        all_nn = s.get_all_neighbors(0.2, True)
        for site, nns in zip(s, all_nn):
            for nn in nns:
                assert nn[0].is_periodic_image(s[nn[2]])
                d = sum((site.coords - nn[0].coords) ** 2) ** 0.5
                assert round(abs(d - nn[1]), 7) == 0
        assert list(map(len, all_nn)) == [2, 2, 2, 0]

    def test_get_all_neighbors_small_cutoff(self):
        s = Structure(
            Lattice.cubic(2),
            ["Li", "Li", "Li", "Si"],
            [[3.1] * 3, [0.11] * 3, [-1.91] * 3, [0.5] * 3],
        )
        all_nn = s.get_all_neighbors(1e-5, True)
        assert len(all_nn) == len(s)
        assert [] == all_nn[0]

        all_nn = s.get_all_neighbors(0, True)
        assert len(all_nn) == len(s)
        assert [] == all_nn[0]

    def test_coincide_sites(self):
        s = Structure(
            Lattice.cubic(5),
            ["Li", "Li", "Li"],
            [[0.1, 0.1, 0.1], [0.1, 0.1, 0.1], [3, 3, 3]],
            coords_are_cartesian=True,
        )
        all_nn = s.get_all_neighbors(1e-5, True)
        assert [len(i) for i in all_nn] == [0, 0, 0]

    def test_get_all_neighbors_equal(self):
        with pytest.warns(FutureWarning, match="get_all_neighbors_old is deprecated"):
            s = Structure(
                Lattice.cubic(2),
                ["Li", "Li", "Li", "Si"],
                [[3.1] * 3, [0.11] * 3, [-1.91] * 3, [0.5] * 3],
            )
            nn_traditional = s.get_all_neighbors_old(4, include_index=True, include_image=True, include_site=True)
            nn_cell_lists = s.get_all_neighbors(4, include_index=True, include_image=True)

            for i in range(4):
                assert len(nn_traditional[i]) == len(nn_cell_lists[i])
                assert (
                    np.linalg.norm(
                        np.array(sorted(j[1] for j in nn_traditional[i]))
                        - np.array(sorted(j[1] for j in nn_cell_lists[i]))
                    )
                    < 1e-3
                )

    def test_get_dist_matrix(self):
        ans = [[0.0, 2.3516318], [2.3516318, 0.0]]
        self.assertArrayAlmostEqual(self.struct.distance_matrix, ans)

    def test_to_from_file_string(self):
        with ScratchDir("."):
            for fmt in ["cif", "json", "poscar", "cssr"]:
                s = self.struct.to(fmt=fmt)
                assert s is not None
                ss = IStructure.from_str(s, fmt=fmt)
                self.assertArrayAlmostEqual(ss.lattice.parameters, self.struct.lattice.parameters, decimal=5)
                self.assertArrayAlmostEqual(ss.frac_coords, self.struct.frac_coords)
                assert isinstance(ss, IStructure)

            assert "Fd-3m" in self.struct.to(fmt="CIF", symprec=0.1)

            self.struct.to(filename="POSCAR.testing")
            assert os.path.exists("POSCAR.testing")

            self.struct.to(filename="Si_testing.yaml")
            assert os.path.exists("Si_testing.yaml")
            s = Structure.from_file("Si_testing.yaml")
            assert s == self.struct
            # Test Path support.
            s = Structure.from_file(Path("Si_testing.yaml"))
            assert s == self.struct

            # Test .yml extension works too.
            os.replace("Si_testing.yaml", "Si_testing.yml")
            s = Structure.from_file("Si_testing.yml")
            assert s == self.struct

            with pytest.raises(ValueError):
                self.struct.to(filename="whatever")
            with pytest.raises(ValueError):
                self.struct.to(fmt="badformat")

            self.struct.to(filename="POSCAR.testing.gz")
            s = Structure.from_file("POSCAR.testing.gz")
            assert s == self.struct

    def test_pbc(self):
        self.assertArrayEqual(self.struct.pbc, (True, True, True))
        assert self.struct.is_3d_periodic
        struct_pbc = Structure(self.lattice_pbc, ["Si"] * 2, self.struct.frac_coords)
        self.assertArrayEqual(struct_pbc.pbc, (True, True, False))
        assert not struct_pbc.is_3d_periodic


class StructureTest(PymatgenTest):
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
        self.structure = Structure(lattice, ["Si", "Si"], coords)

    def test_mutable_sequence_methods(self):
        s = self.structure
        s[0] = "Fe"
        assert s.formula == "Fe1 Si1"
        s[0] = "Fe", [0.5, 0.5, 0.5]
        assert s.formula == "Fe1 Si1"
        self.assertArrayAlmostEqual(s[0].frac_coords, [0.5, 0.5, 0.5])
        s.reverse()
        assert s[0].specie == Element("Si")
        self.assertArrayAlmostEqual(s[0].frac_coords, [0.75, 0.5, 0.75])
        s[0] = {"Mn": 0.5}
        assert s.formula == "Mn0.5 Fe1"
        del s[1]
        assert s.formula == "Mn0.5"
        s[0] = "Fe", [0.9, 0.9, 0.9], {"magmom": 5}
        assert s.formula == "Fe1"
        assert s[0].magmom == 5

        # Test atomic replacement.
        s["Fe"] = "Mn"
        assert s.formula == "Mn1"

        # Test slice replacement.
        s = PymatgenTest.get_structure("Li2O")
        s[0:2] = "S"
        assert s.formula == "Li1 S2"

    def test_non_hash(self):
        with pytest.raises(TypeError):
            dict([(self.structure, 1)])

    def test_sort(self):
        s = self.structure
        s[0] = "F"
        s.sort()
        assert s[0].species_string == "Si"
        assert s[1].species_string == "F"
        s.sort(key=lambda site: site.species_string)
        assert s[0].species_string == "F"
        assert s[1].species_string == "Si"
        s.sort(key=lambda site: site.species_string, reverse=True)
        assert s[0].species_string == "Si"
        assert s[1].species_string == "F"

    def test_append_insert_remove_replace_substitute(self):
        s = self.structure
        s.insert(1, "O", [0.5, 0.5, 0.5])
        assert s.formula == "Si2 O1"
        assert s.ntypesp == 2
        assert s.symbol_set == ("O", "Si")
        assert s.indices_from_symbol("Si") == (0, 2)
        assert s.indices_from_symbol("O") == (1,)
        del s[2]
        assert s.formula == "Si1 O1"
        assert s.indices_from_symbol("Si") == (0,)
        assert s.indices_from_symbol("O") == (1,)
        s.append("N", [0.25, 0.25, 0.25])
        assert s.formula == "Si1 N1 O1"
        assert s.ntypesp == 3
        assert s.symbol_set == ("N", "O", "Si")
        assert s.indices_from_symbol("Si") == (0,)
        assert s.indices_from_symbol("O") == (1,)
        assert s.indices_from_symbol("N") == (2,)
        s[0] = "Ge"
        assert s.formula == "Ge1 N1 O1"
        assert s.symbol_set == ("Ge", "N", "O")
        s.replace_species({"Ge": "Si"})
        assert s.formula == "Si1 N1 O1"
        assert s.ntypesp == 3

        s.replace_species({"Si": {"Ge": 0.5, "Si": 0.5}})
        assert s.formula == "Si0.5 Ge0.5 N1 O1"
        # this should change the .5Si .5Ge sites to .75Si .25Ge
        s.replace_species({"Ge": {"Ge": 0.5, "Si": 0.5}})
        assert s.formula == "Si0.75 Ge0.25 N1 O1"

        assert s.ntypesp == 4

        s.replace_species({"Ge": "Si"})
        s.substitute(1, "hydroxyl")
        assert s.formula == "Si1 H1 N1 O1"
        assert s.symbol_set == ("H", "N", "O", "Si")
        # Distance between O and H
        assert round(abs(s.get_distance(2, 3) - 0.96), 7) == 0
        # Distance between Si and H
        assert round(abs(s.get_distance(0, 3) - 2.09840889), 7) == 0

        s.remove_species(["H"])
        assert s.formula == "Si1 N1 O1"

        s.remove_sites([1, 2])
        assert s.formula == "Si1"

    def test_add_remove_site_property(self):
        s = self.structure
        s.add_site_property("charge", [4.1, -5])
        assert s[0].charge == 4.1
        assert s[1].charge == -5
        s.add_site_property("magmom", [3, 2])
        assert s[0].charge == 4.1
        assert s[0].magmom == 3
        s.remove_site_property("magmom")
        with pytest.raises(AttributeError):
            s[0].magmom

    def test_propertied_structure(self):
        # Make sure that site properties are set to None for missing values.
        s = self.structure
        s.add_site_property("charge", [4.1, -5])
        s.append("Li", [0.3, 0.3, 0.3])
        assert len(s.site_properties["charge"]) == 3

    def test_perturb(self):
        d = 0.1
        pre_perturbation_sites = self.structure.copy()
        self.structure.perturb(distance=d)
        post_perturbation_sites = self.structure.sites

        for i, x in enumerate(pre_perturbation_sites):
            assert round(abs(x.distance(post_perturbation_sites[i]) - d), 3) == 0, "Bad perturbation distance"

        structure2 = pre_perturbation_sites.copy()
        structure2.perturb(distance=d, min_distance=0)
        post_perturbation_sites2 = structure2.sites

        for i, x in enumerate(pre_perturbation_sites):
            assert x.distance(post_perturbation_sites2[i]) <= d
            assert x.distance(post_perturbation_sites2[i]) >= 0

    def test_add_oxidation_states(self):
        oxidation_states = {"Si": -4}
        self.structure.add_oxidation_state_by_element(oxidation_states)
        for site in self.structure:
            for k in site.species:
                assert k.oxi_state == oxidation_states[k.symbol], "Wrong oxidation state assigned!"
        oxidation_states = {"Fe": 2}
        with pytest.raises(ValueError):
            self.structure.add_oxidation_state_by_element(oxidation_states)
        self.structure.add_oxidation_state_by_site([2, -4])
        assert self.structure[0].specie.oxi_state == 2
        with pytest.raises(ValueError):
            self.structure.add_oxidation_state_by_site([1])

    def test_remove_oxidation_states(self):
        co_elem = Element("Co")
        o_elem = Element("O")
        co_specie = Species("Co", 2)
        o_specie = Species("O", -2)
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = Lattice.cubic(10)
        s_elem = Structure(lattice, [co_elem, o_elem], coords)
        s_specie = Structure(lattice, [co_specie, o_specie], coords)
        s_specie.remove_oxidation_states()
        assert s_elem == s_specie, "Oxidation state remover failed"

    def test_add_oxidation_states_by_guess(self):
        s = PymatgenTest.get_structure("Li2O")
        s.add_oxidation_state_by_guess()
        for i in s:
            assert i.specie in [Species("Li", 1), Species("O", -2)]

    def test_add_remove_spin_states(self):

        latt = Lattice.cubic(4.17)
        species = ["Ni", "O"]
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        nio = Structure.from_spacegroup(225, latt, species, coords)

        # should do nothing, but not fail
        nio.remove_spin()

        spins = {"Ni": 5}
        nio.add_spin_by_element(spins)
        assert nio[0].specie.spin == 5, "Failed to add spin states"

        nio.remove_spin()
        with pytest.raises(AttributeError):
            nio[0].specie.spin

        spins = [5, -5, -5, 5, 0, 0, 0, 0]  # AFM on (001)
        nio.add_spin_by_site(spins)
        assert nio[1].specie.spin == -5, "Failed to add spin states"

    def test_apply_operation(self):
        op = SymmOp.from_axis_angle_and_translation([0, 0, 1], 90)
        s = self.structure.copy()
        s.apply_operation(op)
        self.assertArrayAlmostEqual(
            s.lattice.matrix,
            [
                [0.000000, 3.840198, 0.000000],
                [-3.325710, 1.920099, 0.000000],
                [2.217138, -0.000000, 3.135509],
            ],
            5,
        )

        op = SymmOp([[1, 1, 0, 0.5], [1, 0, 0, 0.5], [0, 0, 1, 0.5], [0, 0, 0, 1]])
        s = self.structure.copy()
        s.apply_operation(op, fractional=True)
        self.assertArrayAlmostEqual(
            s.lattice.matrix,
            [
                [5.760297, 3.325710, 0.000000],
                [3.840198, 0.000000, 0.000000],
                [0.000000, -2.217138, 3.135509],
            ],
            5,
        )

    def test_apply_strain(self):
        s = self.structure
        initial_coord = s[1].coords
        s.apply_strain(0.01)
        assert pytest.approx(s.lattice.abc) == (3.8785999130369997, 3.878600984287687, 3.8785999130549516)
        self.assertArrayAlmostEqual(s[1].coords, initial_coord * 1.01)
        a1, b1, c1 = s.lattice.abc
        s.apply_strain([0.1, 0.2, 0.3])
        a2, b2, c2 = s.lattice.abc
        assert round(abs(a2 / a1 - 1.1), 7) == 0
        assert round(abs(b2 / b1 - 1.2), 7) == 0
        assert round(abs(c2 / c1 - 1.3), 7) == 0

    def test_scale_lattice(self):
        initial_coord = self.structure[1].coords
        self.structure.scale_lattice(self.structure.volume * 1.01**3)
        self.assertArrayAlmostEqual(
            self.structure.lattice.abc,
            (3.8785999130369997, 3.878600984287687, 3.8785999130549516),
        )
        self.assertArrayAlmostEqual(self.structure[1].coords, initial_coord * 1.01)

    def test_translate_sites(self):
        self.structure.translate_sites([0, 1], [0.5, 0.5, 0.5], frac_coords=True)
        self.assertArrayAlmostEqual(self.structure.frac_coords[0], [0.5, 0.5, 0.5])

        self.structure.translate_sites([0], [0.5, 0.5, 0.5], frac_coords=False)
        self.assertArrayAlmostEqual(self.structure.cart_coords[0], [3.38014845, 1.05428585, 2.06775453])

        self.structure.translate_sites([0], [0.5, 0.5, 0.5], frac_coords=True, to_unit_cell=False)
        self.assertArrayAlmostEqual(self.structure.frac_coords[0], [1.00187517, 1.25665291, 1.15946374])

        lattice_pbc = Lattice(self.structure.lattice.matrix, pbc=(True, True, False))
        struct_pbc = Structure(lattice_pbc, ["Si"], [[0.75, 0.75, 0.75]])
        struct_pbc.translate_sites([0], [0.5, 0.5, 0.5], frac_coords=True, to_unit_cell=True)
        self.assertArrayAlmostEqual(struct_pbc.frac_coords[0], [0.25, 0.25, 1.25])

    def test_rotate_sites(self):
        self.structure.rotate_sites(
            indices=[1],
            theta=2.0 * np.pi / 3.0,
            anchor=self.structure.sites[0].coords,
            to_unit_cell=False,
        )
        self.assertArrayAlmostEqual(self.structure.frac_coords[1], [-1.25, 1.5, 0.75], decimal=6)
        self.structure.rotate_sites(
            indices=[1],
            theta=2.0 * np.pi / 3.0,
            anchor=self.structure.sites[0].coords,
            to_unit_cell=True,
        )
        self.assertArrayAlmostEqual(self.structure.frac_coords[1], [0.75, 0.5, 0.75], decimal=6)

    def test_mul(self):
        self.structure *= [2, 1, 1]
        assert self.structure.formula == "Si4"
        s = [2, 1, 1] * self.structure
        assert s.formula == "Si8"
        assert isinstance(s, Structure)
        s = self.structure * [[1, 0, 0], [2, 1, 0], [0, 0, 2]]
        assert s.formula == "Si8"
        self.assertArrayAlmostEqual(s.lattice.abc, [7.6803959, 17.5979979, 7.6803959])

    def test_make_supercell(self):
        self.structure.make_supercell([2, 1, 1])
        assert self.structure.formula == "Si4"
        self.structure.make_supercell([[1, 0, 0], [2, 1, 0], [0, 0, 1]])
        assert self.structure.formula == "Si4"
        self.structure.make_supercell(2)
        assert self.structure.formula == "Si32"
        self.assertArrayAlmostEqual(self.structure.lattice.abc, [15.360792, 35.195996, 7.680396], 5)

    def test_disordered_supercell_primitive_cell(self):
        l = Lattice.cubic(2)
        f = [[0.5, 0.5, 0.5]]
        sp = [{"Si": 0.54738}]
        s = Structure(l, sp, f)
        # this supercell often breaks things
        s.make_supercell([[0, -1, 1], [-1, 1, 0], [1, 1, 1]])
        assert len(s.get_primitive_structure()) == 1

    def test_another_supercell(self):
        # this is included b/c for some reason the old algo was failing on it
        s = self.structure.copy()
        s.make_supercell([[0, 2, 2], [2, 0, 2], [2, 2, 0]])
        assert s.formula == "Si32"
        s = self.structure.copy()
        s.make_supercell([[0, 2, 0], [1, 0, 0], [0, 0, 1]])
        assert s.formula == "Si4"

    def test_to_from_dict(self):
        d = self.structure.as_dict()
        s2 = Structure.from_dict(d)
        assert isinstance(s2, Structure)

    def test_default_dict_attrs(self):
        d = self.structure.as_dict()
        assert d["charge"] == 0

    def test_to_from_abivars(self):
        """Test as_dict, from_dict with fmt == abivars."""
        d = self.structure.as_dict(fmt="abivars")
        s2 = Structure.from_dict(d, fmt="abivars")
        assert s2 == self.structure
        assert isinstance(s2, Structure)

    def test_to_from_file_string(self):
        with ScratchDir("."):
            for fmt in ["cif", "json", "poscar", "cssr", "yaml", "xsf", "res"]:
                s = self.structure.to(fmt=fmt)
                assert s is not None
                ss = Structure.from_str(s, fmt=fmt)
                self.assertArrayAlmostEqual(ss.lattice.parameters, self.structure.lattice.parameters, decimal=5)
                self.assertArrayAlmostEqual(ss.frac_coords, self.structure.frac_coords)
                assert isinstance(ss, Structure)

            self.structure.to(filename="POSCAR.testing")
            assert os.path.exists("POSCAR.testing")

            self.structure.to(filename="structure_testing.json")
            assert Structure.from_file("structure_testing.json") == self.structure

    def test_from_spacegroup(self):
        s1 = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3), ["Li", "O"], [[0.25, 0.25, 0.25], [0, 0, 0]])
        assert s1.formula == "Li8 O4"
        s2 = Structure.from_spacegroup(225, Lattice.cubic(3), ["Li", "O"], [[0.25, 0.25, 0.25], [0, 0, 0]])
        assert s1 == s2

        s2 = Structure.from_spacegroup(
            225,
            Lattice.cubic(3),
            ["Li", "O"],
            [[0.25, 0.25, 0.25], [0, 0, 0]],
            site_properties={"charge": [1, -2]},
        )
        assert sum(s2.site_properties["charge"]) == 0

        s = Structure.from_spacegroup("Pm-3m", Lattice.cubic(3), ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        assert s.formula == "Cs1 Cl1"

        with pytest.raises(ValueError):
            Structure.from_spacegroup(
                "Pm-3m",
                Lattice.tetragonal(1, 3),
                ["Cs", "Cl"],
                [[0, 0, 0], [0.5, 0.5, 0.5]],
            )

        with pytest.raises(ValueError):
            Structure.from_spacegroup(
                "Pm-3m",
                Lattice.cubic(3),
                ["Cs"],
                [[0, 0, 0], [0.5, 0.5, 0.5]],
            )
        from fractions import Fraction

        s = Structure.from_spacegroup(139, np.eye(3), ["H"], [[Fraction(1, 2), Fraction(1, 4), Fraction(0)]])
        assert len(s) == 8

    def test_from_magnetic_spacegroup(self):

        # AFM MnF
        s1 = Structure.from_magnetic_spacegroup(
            "P4_2'/mnm'",
            Lattice.tetragonal(4.87, 3.30),
            ["Mn", "F"],
            [[0, 0, 0], [0.30, 0.30, 0.00]],
            {"magmom": [4, 0]},
        )

        assert s1.formula == "Mn2 F4"
        assert sum(map(float, s1.site_properties["magmom"])) == 0
        assert max(map(float, s1.site_properties["magmom"])) == 4
        assert min(map(float, s1.site_properties["magmom"])) == -4

        # AFM LaMnO3, ordered on (001) planes
        s2 = Structure.from_magnetic_spacegroup(
            "Pn'ma'",
            Lattice.orthorhombic(5.75, 7.66, 5.53),
            ["La", "Mn", "O", "O"],
            [
                [0.05, 0.25, 0.99],
                [0.00, 0.00, 0.50],
                [0.48, 0.25, 0.08],
                [0.31, 0.04, 0.72],
            ],
            {"magmom": [0, Magmom([4, 0, 0]), 0, 0]},
        )

        assert s2.formula == "La4 Mn4 O12"
        assert sum(map(float, s2.site_properties["magmom"])) == 0
        assert max(map(float, s2.site_properties["magmom"])) == 4
        assert min(map(float, s2.site_properties["magmom"])) == -4

    def test_merge_sites(self):
        species = [
            {"Ag": 0.5},
            {"Cl": 0.25},
            {"Cl": 0.1},
            {"Ag": 0.5},
            {"F": 0.15},
            {"F": 0.1},
        ]
        coords = [
            [0, 0, 0],
            [0.5, 0.5, 0.5],
            [0.5, 0.5, 0.5],
            [0, 0, 0],
            [0.5, 0.5, 1.501],
            [0.5, 0.5, 1.501],
        ]
        s = Structure(Lattice.cubic(1), species, coords)
        s.merge_sites(mode="s")
        assert s[0].specie.symbol == "Ag"
        assert s[1].species == Composition({"Cl": 0.35, "F": 0.25})
        self.assertArrayAlmostEqual(s[1].frac_coords, [0.5, 0.5, 0.5005])

        # Test for TaS2 with spacegroup 166 in 160 setting.
        l = Lattice.hexagonal(3.374351, 20.308941)
        species = ["Ta", "S", "S"]
        coords = [
            [0.000000, 0.000000, 0.944333],
            [0.333333, 0.666667, 0.353424],
            [0.666667, 0.333333, 0.535243],
        ]
        tas2 = Structure.from_spacegroup(160, l, species, coords)
        assert len(tas2) == 13
        tas2.merge_sites(mode="d")
        assert len(tas2) == 9

        l = Lattice.hexagonal(3.587776, 19.622793)
        species = ["Na", "V", "S", "S"]
        coords = [
            [0.333333, 0.666667, 0.165000],
            [0.000000, 0.000000, 0.998333],
            [0.333333, 0.666667, 0.399394],
            [0.666667, 0.333333, 0.597273],
        ]
        navs2 = Structure.from_spacegroup(160, l, species, coords)
        assert len(navs2) == 18
        navs2.merge_sites(mode="d")
        assert len(navs2) == 12

        # Test that we can average the site properties that are floats
        l = Lattice.hexagonal(3.587776, 19.622793)
        species = ["Na", "V", "S", "S"]
        coords = [
            [0.333333, 0.666667, 0.165000],
            [0.000000, 0.000000, 0.998333],
            [0.333333, 0.666667, 0.399394],
            [0.666667, 0.333333, 0.597273],
        ]
        site_props = {"prop1": [3.0, 5.0, 7.0, 11.0]}
        navs2 = Structure.from_spacegroup(160, l, species, coords, site_properties=site_props)
        navs2.insert(0, "Na", coords[0], properties={"prop1": 100.0})
        navs2.merge_sites(mode="a")
        assert len(navs2) == 12
        assert 51.5 in [itr.properties["prop1"] for itr in navs2.sites]

    def test_properties(self):
        assert self.structure.num_sites == len(self.structure)
        self.structure.make_supercell(2)
        self.structure[1] = "C"
        sites = list(self.structure.group_by_types())
        assert sites[-1].specie.symbol == "C"
        self.structure.add_oxidation_state_by_element({"Si": 4, "C": 2})
        assert self.structure.charge == 62

    def test_set_item(self):
        s = self.structure.copy()
        s[0] = "C"
        assert s.formula == "Si1 C1"
        s[(0, 1)] = "Ge"
        assert s.formula == "Ge2"
        s[0:2] = "Sn"
        assert s.formula == "Sn2"

        s = self.structure.copy()
        s["Si"] = "C"
        assert s.formula == "C2"
        s["C"] = "C0.25Si0.5"
        assert s.formula == "Si1 C0.5"
        s["C"] = "C0.25Si0.5"
        assert s.formula == "Si1.25 C0.125"

    def test_init_error(self):
        with pytest.raises(StructureError):
            Structure(
                Lattice.cubic(3),
                ["Si"],
                [[0, 0, 0], [0.5, 0.5, 0.5]],
            )

    def test_from_sites(self):
        self.structure.add_site_property("hello", [1, 2])
        s = Structure.from_sites(self.structure, to_unit_cell=True)
        assert s.site_properties["hello"][1] == 2

    def test_charge(self):
        s = Structure.from_sites(self.structure)
        assert s.charge == 0, "Initial Structure not defaulting to behavior in SiteCollection"
        s.add_oxidation_state_by_site([1, 1])
        assert s.charge == 2, "Initial Structure not defaulting to behavior in SiteCollection"
        s = Structure.from_sites(s, charge=1)
        assert s.charge == 1, "Overall charge not being stored in separate property"
        s = s.copy()
        assert s.charge == 1, "Overall charge not being copied properly with no sanitization"
        s = s.copy(sanitize=True)
        assert s.charge == 1, "Overall charge not being copied properly with sanitization"
        super_cell = s * 3
        assert super_cell.charge == 27, "Overall charge is not being properly multiplied in IStructure __mul__"
        assert "Overall Charge: +1" in str(s), "String representation not adding charge"
        sorted_s = super_cell.get_sorted_structure()
        assert sorted_s.charge == 27, "Overall charge is not properly copied during structure sorting"
        super_cell.set_charge(25)
        assert super_cell.charge == 25, "Set charge not properly modifying _charge"

    def test_vesta_lattice_matrix(self):
        silica_zeolite = Molecule.from_file(self.TEST_FILES_DIR / "CON_vesta.xyz")

        s_vesta = Structure(
            lattice=Lattice.from_parameters(22.6840, 13.3730, 12.5530, 90, 69.479, 90, True),
            species=silica_zeolite.species,
            coords=silica_zeolite.cart_coords,
            coords_are_cartesian=True,
            to_unit_cell=True,
        )

        s_vesta = s_vesta.get_primitive_structure()
        s_vesta.merge_sites(0.01, "delete")
        assert s_vesta.formula == "Si56 O112"

        broken_s = Structure(
            lattice=Lattice.from_parameters(22.6840, 13.3730, 12.5530, 90, 69.479, 90),
            species=silica_zeolite.species,
            coords=silica_zeolite.cart_coords,
            coords_are_cartesian=True,
            to_unit_cell=True,
        )

        broken_s.merge_sites(0.01, "delete")
        assert broken_s.formula == "Si56 O134"

    def test_extract_cluster(self):
        coords = [
            [0.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, 1.089000],
            [1.026719, 0.000000, -0.363000],
            [-0.513360, -0.889165, -0.363000],
            [-0.513360, 0.889165, -0.363000],
        ]
        ch4 = ["C", "H", "H", "H", "H"]

        species = []
        allcoords = []
        for vec in ([0, 0, 0], [4, 0, 0], [0, 4, 0], [4, 4, 0]):
            species.extend(ch4)
            for c in coords:
                allcoords.append(np.array(c) + vec)

        structure = Structure(Lattice.cubic(10), species, allcoords, coords_are_cartesian=True)

        for site in structure:
            if site.specie.symbol == "C":
                cluster = Molecule.from_sites(structure.extract_cluster([site]))
                assert cluster.formula == "H4 C1"

    @unittest.skipIf(m3gnet is None, "Relaxation test requires m3gnet.")
    def test_relax(self):
        structure = self.get_structure("Si")
        relaxed = structure.relax()
        assert round(abs(relaxed.lattice.a - 3.849563), 4) == 0

    def test_from_prototype(self):
        for pt in ["bcc", "fcc", "hcp", "diamond"]:
            s = Structure.from_prototype(pt, ["C"], a=3, c=4)
            assert isinstance(s, Structure)

        with pytest.raises(ValueError):
            Structure.from_prototype("hcp", ["C"], a=3)

        s = Structure.from_prototype("rocksalt", ["Li", "Cl"], a=2.56)
        assert (
            str(s)
            == """Full Formula (Li4 Cl4)
Reduced Formula: LiCl
abc   :   2.560000   2.560000   2.560000
angles:  90.000000  90.000000  90.000000
pbc   :       True       True       True
Sites (8)
  #  SP      a    b    c
---  ----  ---  ---  ---
  0  Li    0    0    0
  1  Li    0.5  0.5  0
  2  Li    0.5  0    0.5
  3  Li    0    0.5  0.5
  4  Cl    0.5  0    0.5
  5  Cl    0    0.5  0.5
  6  Cl    0.5  0.5  0
  7  Cl    0    0    0"""
        )
        for pt in ("cscl", "fluorite", "antifluorite", "zincblende"):
            s = Structure.from_prototype(pt, ["Cs", "Cl"], a=5)
            assert s.lattice.is_orthogonal


class IMoleculeTest(PymatgenTest):
    def setUp(self):
        coords = [
            [0.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, 1.089000],
            [1.026719, 0.000000, -0.363000],
            [-0.513360, -0.889165, -0.363000],
            [-0.513360, 0.889165, -0.363000],
        ]
        self.coords = coords
        self.mol = Molecule(["C", "H", "H", "H", "H"], coords)

    def test_set_item(self):
        s = self.mol.copy()
        s[0] = "Si"
        assert s.formula == "Si1 H4"
        s[(0, 1)] = "Ge"
        assert s.formula == "Ge2 H3"
        s[0:2] = "Sn"
        assert s.formula == "Sn2 H3"

        s = self.mol.copy()
        s["H"] = "F"
        assert s.formula == "C1 F4"
        s["C"] = "C0.25Si0.5"
        assert s.formula == "Si0.5 C0.25 F4"
        s["C"] = "C0.25Si0.5"
        assert s.formula == "Si0.625 C0.0625 F4"

    def test_bad_molecule(self):
        coords = [
            [0.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, 1.089000],
            [1.026719, 0.000000, -0.363000],
            [-0.513360, -0.889165, -0.363000],
            [-0.513360, 0.889165, -0.363000],
            [-0.513360, 0.889165, -0.36301],
        ]
        with pytest.raises(StructureError):
            Molecule(
                ["C", "H", "H", "H", "H", "H"],
                coords,
                validate_proximity=True,
            )

    def test_get_angle_dihedral(self):
        assert round(abs(self.mol.get_angle(1, 0, 2) - 109.47122144618737), 7) == 0
        assert round(abs(self.mol.get_angle(3, 1, 2) - 60.00001388659683), 7) == 0
        assert round(abs(self.mol.get_dihedral(0, 1, 2, 3) - -35.26438851071765), 7) == 0

        coords = []
        coords.append([0, 0, 0])
        coords.append([0, 0, 1])
        coords.append([0, 1, 1])
        coords.append([1, 1, 1])
        self.mol2 = Molecule(["C", "O", "N", "S"], coords)
        assert round(abs(self.mol2.get_dihedral(0, 1, 2, 3) - -90), 7) == 0

    def test_get_covalent_bonds(self):
        assert len(self.mol.get_covalent_bonds()) == 4

    def test_properties(self):
        assert len(self.mol) == 5
        assert self.mol.is_ordered
        assert self.mol.formula == "H4 C1"

    def test_repr_str(self):
        ans = """Full Formula (H4 C1)
Reduced Formula: H4C
Charge = 0.0, Spin Mult = 1
Sites (5)
0 C     0.000000     0.000000     0.000000
1 H     0.000000     0.000000     1.089000
2 H     1.026719     0.000000    -0.363000
3 H    -0.513360    -0.889165    -0.363000
4 H    -0.513360     0.889165    -0.363000"""
        assert str(self.mol) == ans
        ans = """Molecule Summary
Site: C (0.0000, 0.0000, 0.0000)
Site: H (0.0000, 0.0000, 1.0890)
Site: H (1.0267, 0.0000, -0.3630)
Site: H (-0.5134, -0.8892, -0.3630)
Site: H (-0.5134, 0.8892, -0.3630)"""
        assert repr(self.mol) == ans

    def test_site_properties(self):
        propertied_mol = Molecule(
            ["C", "H", "H", "H", "H"],
            self.coords,
            site_properties={"magmom": [0.5, -0.5, 1, 2, 3]},
        )
        assert propertied_mol[0].magmom == 0.5
        assert propertied_mol[1].magmom == -0.5

    def test_get_boxed_structure(self):
        s = self.mol.get_boxed_structure(9, 9, 9)
        # C atom should be in center of box.
        self.assertArrayAlmostEqual(s[4].frac_coords, [0.50000001, 0.5, 0.5])
        self.assertArrayAlmostEqual(s[1].frac_coords, [0.6140799, 0.5, 0.45966667])
        with pytest.raises(ValueError):
            self.mol.get_boxed_structure(1, 1, 1)
        s2 = self.mol.get_boxed_structure(5, 5, 5, (2, 3, 4))
        assert len(s2) == 24 * 5
        assert s2.lattice.abc == (10, 15, 20)

        # Test offset option
        s3 = self.mol.get_boxed_structure(9, 9, 9, offset=[0.5, 0.5, 0.5])
        self.assertArrayAlmostEqual(s3[4].coords, [5, 5, 5])
        # Test no_cross option
        with pytest.raises(ValueError):
            self.mol.get_boxed_structure(
                5,
                5,
                5,
                offset=[10, 10, 10],
                no_cross=True,
            )

        # Test reorder option
        no_reorder = self.mol.get_boxed_structure(10, 10, 10, reorder=False)
        assert str(s3[0].specie) == "H"
        assert str(no_reorder[0].specie) == "C"

    def test_get_distance(self):
        assert round(abs(self.mol.get_distance(0, 1) - 1.089), 7) == 0

    def test_get_neighbors(self):
        nn = self.mol.get_neighbors(self.mol[0], 1)
        assert len(nn) == 0
        nn = self.mol.get_neighbors(self.mol[0], 2)
        assert len(nn) == 4

    def test_get_neighbors_in_shell(self):
        nn = self.mol.get_neighbors_in_shell([0, 0, 0], 0, 1)
        assert len(nn) == 1
        nn = self.mol.get_neighbors_in_shell([0, 0, 0], 1, 0.9)
        assert len(nn) == 4
        nn = self.mol.get_neighbors_in_shell([0, 0, 0], 1, 0.9)
        assert len(nn) == 4
        nn = self.mol.get_neighbors_in_shell([0, 0, 0], 2, 0.1)
        assert len(nn) == 0

    def test_get_dist_matrix(self):
        ans = [
            [0.0, 1.089, 1.08899995636, 1.08900040717, 1.08900040717],
            [1.089, 0.0, 1.77832952654, 1.7783298026, 1.7783298026],
            [1.08899995636, 1.77832952654, 0.0, 1.77833003783, 1.77833003783],
            [1.08900040717, 1.7783298026, 1.77833003783, 0.0, 1.77833],
            [1.08900040717, 1.7783298026, 1.77833003783, 1.77833, 0.0],
        ]
        self.assertArrayAlmostEqual(self.mol.distance_matrix, ans)

    def test_break_bond(self):
        (mol1, mol2) = self.mol.break_bond(0, 1)
        assert mol1.formula == "H3 C1"
        assert mol2.formula == "H1"

    def test_prop(self):
        assert self.mol.charge == 0
        assert self.mol.spin_multiplicity == 1
        assert self.mol.nelectrons == 10
        self.assertArrayAlmostEqual(self.mol.center_of_mass, [0, 0, 0])
        with pytest.raises(ValueError):
            Molecule(
                ["C", "H", "H", "H", "H"],
                self.coords,
                charge=1,
                spin_multiplicity=1,
            )
        mol = Molecule(["C", "H", "H", "H", "H"], self.coords, charge=1)
        assert mol.spin_multiplicity == 2
        assert mol.nelectrons == 9

        # Triplet O2
        mol = IMolecule(["O"] * 2, [[0, 0, 0], [0, 0, 1.2]], spin_multiplicity=3)
        assert mol.spin_multiplicity == 3

    def test_no_spin_check(self):
        coords = [
            [0.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, 1.089000],
            [1.026719, 0.000000, -0.363000],
            [-0.513360, -0.889165, -0.363000],
        ]
        with pytest.raises(ValueError):
            mol = IMolecule(["C", "H", "H", "H"], coords, charge=0, spin_multiplicity=1)
        mol = IMolecule(["C", "H", "H", "H"], coords, charge=0, spin_multiplicity=1, charge_spin_check=False)
        assert mol.spin_multiplicity == 1
        assert mol.charge == 0

    def test_equal(self):
        mol = IMolecule(["C", "H", "H", "H", "H"], self.coords, charge=1)
        assert mol != self.mol

    def test_get_centered_molecule(self):
        mol = IMolecule(["O"] * 2, [[0, 0, 0], [0, 0, 1.2]], spin_multiplicity=3)
        centered = mol.get_centered_molecule()
        self.assertArrayAlmostEqual(centered.center_of_mass, [0, 0, 0])

    def test_to_from_dict(self):
        d = self.mol.as_dict()
        mol2 = IMolecule.from_dict(d)
        assert isinstance(mol2, IMolecule)
        propertied_mol = Molecule(
            ["C", "H", "H", "H", "H"],
            self.coords,
            charge=1,
            site_properties={"magmom": [0.5, -0.5, 1, 2, 3]},
        )
        d = propertied_mol.as_dict()
        assert d["sites"][0]["properties"]["magmom"] == 0.5
        mol = Molecule.from_dict(d)
        assert propertied_mol == mol
        assert mol[0].magmom == 0.5
        assert mol.formula == "H4 C1"
        assert mol.charge == 1

    def test_default_dict_attrs(self):
        d = self.mol.as_dict()
        assert d["charge"] == 0
        assert d["spin_multiplicity"] == 1

    def test_to_from_file_string(self):
        for fmt in ["xyz", "json", "g03", "yaml"]:
            s = self.mol.to(fmt=fmt)
            assert s is not None
            m = IMolecule.from_str(s, fmt=fmt)
            assert m == self.mol
            assert isinstance(m, IMolecule)

        self.mol.to(filename="CH4_testing.xyz")
        assert os.path.exists("CH4_testing.xyz")
        os.remove("CH4_testing.xyz")
        self.mol.to(filename="CH4_testing.yaml")
        assert os.path.exists("CH4_testing.yaml")
        mol = Molecule.from_file("CH4_testing.yaml")
        assert self.mol == mol
        os.remove("CH4_testing.yaml")


class MoleculeTest(PymatgenTest):
    def setUp(self):
        coords = [
            [0.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, 1.089000],
            [1.026719, 0.000000, -0.363000],
            [-0.513360, -0.889165, -0.363000],
            [-0.513360, 0.889165, -0.363000],
        ]
        self.mol = Molecule(["C", "H", "H", "H", "H"], coords)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_mutable_sequence_methods(self):
        s = self.mol
        s[1] = ("F", [0.5, 0.5, 0.5])
        assert s.formula == "H3 C1 F1"
        self.assertArrayAlmostEqual(s[1].coords, [0.5, 0.5, 0.5])
        s.reverse()
        assert s[0].specie == Element("H")
        self.assertArrayAlmostEqual(s[0].coords, [-0.513360, 0.889165, -0.363000])
        del s[1]
        assert s.formula == "H2 C1 F1"
        s[3] = "N", [0, 0, 0], {"charge": 4}
        assert s.formula == "H2 N1 F1"
        assert s[3].charge == 4

    def test_insert_remove_append(self):
        mol = self.mol
        mol.insert(1, "O", [0.5, 0.5, 0.5])
        assert mol.formula == "H4 C1 O1"
        del mol[2]
        assert mol.formula == "H3 C1 O1"
        mol.set_charge_and_spin(0)
        assert mol.spin_multiplicity == 2
        mol.append("N", [1, 1, 1])
        assert mol.formula == "H3 C1 N1 O1"
        with pytest.raises(TypeError):
            dict([(mol, 1)])
        mol.remove_sites([0, 1])
        assert mol.formula == "H3 N1"

    def test_translate_sites(self):
        self.mol.translate_sites([0, 1], [0.5, 0.5, 0.5])
        self.assertArrayEqual(self.mol.cart_coords[0], [0.5, 0.5, 0.5])

    def test_rotate_sites(self):
        self.mol.rotate_sites(theta=np.radians(30))
        self.assertArrayAlmostEqual(self.mol.cart_coords[2], [0.889164737, 0.513359500, -0.363000000])

    def test_replace(self):
        self.mol[0] = "Ge"
        assert self.mol.formula == "Ge1 H4"

        self.mol.replace_species({Element("Ge"): {Element("Ge"): 0.5, Element("Si"): 0.5}})
        assert self.mol.formula == "Si0.5 Ge0.5 H4"

        # this should change the .5Si .5Ge sites to .75Si .25Ge
        self.mol.replace_species({Element("Ge"): {Element("Ge"): 0.5, Element("Si"): 0.5}})
        assert self.mol.formula == "Si0.75 Ge0.25 H4"

        d = 0.1
        pre_perturbation_sites = self.mol.sites[:]
        self.mol.perturb(distance=d)
        post_perturbation_sites = self.mol.sites

        for i, x in enumerate(pre_perturbation_sites):
            assert round(abs(x.distance(post_perturbation_sites[i]) - d), 3) == 0, "Bad perturbation distance"

    def test_add_site_property(self):
        self.mol.add_site_property("charge", [4.1, -2, -2, -2, -2])
        assert self.mol[0].charge == 4.1
        assert self.mol[1].charge == -2

        self.mol.add_site_property("magmom", [3, 2, 2, 2, 2])
        assert self.mol[0].charge == 4.1
        assert self.mol[0].magmom == 3
        self.mol.remove_site_property("magmom")
        with pytest.raises(AttributeError):
            self.mol[0].magmom

    def test_to_from_dict(self):
        self.mol.append("X", [2, 0, 0])
        d = self.mol.as_dict()
        mol2 = Molecule.from_dict(d)
        assert isinstance(mol2, Molecule)
        self.assertMSONable(self.mol)

    def test_apply_operation(self):
        op = SymmOp.from_axis_angle_and_translation([0, 0, 1], 90)
        self.mol.apply_operation(op)
        self.assertArrayAlmostEqual(self.mol[2].coords, [0.000000, 1.026719, -0.363000])

    def test_substitute(self):
        coords = [
            [0.000000, 0.000000, 1.08],
            [0.000000, 0.000000, 0.000000],
            [1.026719, 0.000000, -0.363000],
            [-0.513360, -0.889165, -0.363000],
            [-0.513360, 0.889165, -0.363000],
        ]
        sub = Molecule(["X", "C", "H", "H", "H"], coords)
        self.mol.substitute(1, sub)
        assert round(abs(self.mol.get_distance(0, 4) - 1.54), 7) == 0
        f = Molecule(["X", "F"], [[0, 0, 0], [0, 0, 1.11]])
        self.mol.substitute(2, f)
        assert round(abs(self.mol.get_distance(0, 7) - 1.35), 7) == 0
        oh = Molecule(
            ["X", "O", "H"],
            [[0, 0.780362, -0.456316], [0, 0, 0.114079], [0, -0.780362, -0.456316]],
        )
        self.mol.substitute(1, oh)
        assert round(abs(self.mol.get_distance(0, 7) - 1.43), 7) == 0
        self.mol.substitute(3, "methyl")
        assert self.mol.formula == "H7 C3 O1 F1"
        coords = [
            [0.00000, 1.40272, 0.00000],
            [0.00000, 2.49029, 0.00000],
            [-1.21479, 0.70136, 0.00000],
            [-2.15666, 1.24515, 0.00000],
            [-1.21479, -0.70136, 0.00000],
            [-2.15666, -1.24515, 0.00000],
            [0.00000, -1.40272, 0.00000],
            [0.00000, -2.49029, 0.00000],
            [1.21479, -0.70136, 0.00000],
            [2.15666, -1.24515, 0.00000],
            [1.21479, 0.70136, 0.00000],
            [2.15666, 1.24515, 0.00000],
        ]
        benzene = Molecule(["C", "H", "C", "H", "C", "H", "C", "H", "C", "H", "C", "H"], coords)
        benzene.substitute(1, sub)
        assert benzene.formula == "H8 C7"
        # Carbon attached should be in plane.
        assert round(abs(benzene[11].coords[2] - 0), 7) == 0
        benzene[14] = "Br"
        benzene.substitute(13, sub)
        assert benzene.formula == "H9 C8 Br1"

    def test_to_from_file_string(self):
        for fmt in ["xyz", "json", "g03"]:
            s = self.mol.to(fmt=fmt)
            assert s is not None
            m = Molecule.from_str(s, fmt=fmt)
            assert m == self.mol
            assert isinstance(m, Molecule)

        self.mol.to(filename="CH4_testing.xyz")
        assert os.path.exists("CH4_testing.xyz")
        os.remove("CH4_testing.xyz")

    def test_extract_cluster(self):
        species = self.mol.species * 2
        coords = list(self.mol.cart_coords) + list(self.mol.cart_coords + [10, 0, 0])
        mol = Molecule(species, coords)
        cluster = Molecule.from_sites(mol.extract_cluster([mol[0]]))
        assert mol.formula == "H8 C2"
        assert cluster.formula == "H4 C1"

    def test_no_spin_check(self):
        coords = [
            [0.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, 1.089000],
            [1.026719, 0.000000, -0.363000],
            [-0.513360, -0.889165, -0.363000],
        ]
        with pytest.raises(ValueError):
            mol = Molecule(["C", "H", "H", "H"], coords, charge=0, spin_multiplicity=1)
        mol_valid = Molecule(["C", "H", "H", "H"], coords, charge=0, spin_multiplicity=2)
        with pytest.raises(ValueError):
            mol_valid.set_charge_and_spin(0, 1)
        mol = Molecule(["C", "H", "H", "H"], coords, charge=0, spin_multiplicity=1, charge_spin_check=False)
        assert mol.spin_multiplicity == 1
        assert mol.charge == 0
        mol.set_charge_and_spin(0, 3)
        assert mol.charge == 0
        assert mol.spin_multiplicity == 3


if __name__ == "__main__":
    unittest.main()
