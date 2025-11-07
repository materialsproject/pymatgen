from __future__ import annotations

import itertools
import json
import math
import os
from fractions import Fraction
from pathlib import Path
from shutil import which
from unittest import mock

import numpy as np
import pytest
import requests
import urllib3
from monty.json import MontyDecoder, MontyEncoder
from numpy.testing import assert_allclose, assert_array_equal
from pytest import approx

from pymatgen.core import SETTINGS, Composition, Element, Lattice, Species
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import (
    IMolecule,
    IStructure,
    Molecule,
    Neighbor,
    PeriodicNeighbor,
    Structure,
    StructureError,
)
from pymatgen.electronic_structure.core import Magmom
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.testing import TEST_FILES_DIR, VASP_IN_DIR, MatSciTest

try:
    from ase.atoms import Atoms
    from ase.build import molecule
    from ase.calculators.calculator import Calculator
    from ase.calculators.emt import EMT
except ImportError:
    Atoms = Calculator = EMT = None


ENUM_CMD = which("enum.x") or which("multienum.x")
MCSQS_CMD = which("mcsqs")


class TestNeighbor:
    def test_msonable(self):
        struct = MatSciTest.get_structure("Li2O")
        nn = struct.get_neighbors(struct[0], r=3)
        assert isinstance(nn[0], PeriodicNeighbor)
        str_ = json.dumps(nn, cls=MontyEncoder)
        nn = json.loads(str_, cls=MontyDecoder)
        assert isinstance(nn[0], PeriodicNeighbor)

    def test_neighbor_labels(self):
        comp = Composition("C")
        for label in (None, "", "str label", ("tuple", "label")):
            neighbor = Neighbor(comp, (0, 0, 0), label=label)
            assert neighbor.label == label if label is not None else "C"

            p_neighbor = PeriodicNeighbor(comp, (0, 0, 0), (10, 10, 10), label=label)
            assert p_neighbor.label == label if label is not None else "C"


class TestIStructure(MatSciTest):
    def setup_method(self):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        self.lattice = Lattice(
            [
                [3.8401979337, 0, 0],
                [1.9200989668, 3.3257101909, 0],
                [0, -2.2171384943, 3.1355090603],
            ]
        )
        self.struct = IStructure(self.lattice, ["Si"] * 2, coords)
        assert len(self.struct) == 2, "Wrong number of sites in structure!"
        assert self.struct.is_ordered
        assert self.struct.n_elems == 1
        coords = [[0, 0, 0], [0.0, 0, 0.0000001]]
        with pytest.raises(
            StructureError,
            match=f"sites are less than {self.struct.DISTANCE_TOLERANCE} Angstrom apart",
        ):
            IStructure(self.lattice, ["Si"] * 2, coords, validate_proximity=True)
        self.propertied_structure = IStructure(
            self.lattice,
            ["Si"] * 2,
            coords,
            site_properties={"magmom": [5, -5]},
            properties={"test_property": "test"},
        )
        self.labeled_structure = IStructure(self.lattice, ["Si"] * 2, coords, labels=["Si1", "Si2"])

        self.lattice_pbc = Lattice(
            [
                [3.8401979337, 0, 0],
                [1.9200989668, 3.3257101909, 0],
                [0, -2.2171384943, 3.1355090603],
            ],
            pbc=(True, True, False),
        )
        self.V2O3 = IStructure.from_file(f"{TEST_FILES_DIR}/cif/V2O3.cif")

    @pytest.mark.skipif(not (MCSQS_CMD and ENUM_CMD), reason="enumlib or mcsqs executable not present")
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
        df_struct = self.propertied_structure.as_dataframe()
        assert df_struct.attrs["Reduced Formula"] == self.propertied_structure.reduced_formula
        assert df_struct.shape == (2, 8)
        assert list(df_struct) == ["Species", *"abcxyz", "magmom"]
        assert list(df_struct["magmom"]) == [5, -5]
        assert list(map(str, df_struct["Species"])) == ["Si1", "Si1"]

    def test_equal(self):
        struct = self.struct
        assert struct == struct  # noqa: PLR0124
        assert struct == struct.copy()
        assert struct != 2 * struct

        assert struct != "a" * len(struct)  # GH-2584
        assert struct is not None
        assert struct != 42  # GH-2587

        assert struct == Structure.from_dict(struct.as_dict())

        struct_2 = Structure.from_sites(struct)
        assert struct, struct_2
        struct_2.apply_strain(0.5)
        assert struct != struct_2

    def test_matches(self):
        supercell = self.struct * 2
        assert supercell.matches(self.struct)

    def test_bad_structure(self):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75], [0.75, 0.5, 0.75]]
        with pytest.raises(
            StructureError,
            match=f"sites are less than {Structure.DISTANCE_TOLERANCE} Angstrom apart",
        ):
            IStructure(self.lattice, ["Si"] * 3, coords, validate_proximity=True)
        # these shouldn't raise an error
        IStructure(self.lattice, ["Si"] * 2, coords[:2], validate_proximity=True)
        IStructure(self.lattice, ["Si"], coords[:1], validate_proximity=True)

    def test_volume(self):
        assert self.struct.volume == approx(40.04, abs=1e-2), "Volume wrong!"

    def test_density(self):
        assert self.struct.density == approx(2.33, abs=1e-2), "Incorrect density"

    def test_formula(self):
        assert self.struct.formula == "Si2"
        assert self.labeled_structure.formula == "Si2"
        assert self.propertied_structure.formula == "Si2"
        assert self.V2O3.formula == "V4 O6"

    def test_alphabetical_formula(self):
        assert self.struct.alphabetical_formula == "Si2"
        assert self.labeled_structure.alphabetical_formula == "Si2"
        assert self.propertied_structure.alphabetical_formula == "Si2"
        assert self.V2O3.alphabetical_formula == "O6 V4"

    def test_reduced_formula(self):
        assert self.struct.reduced_formula == "Si"
        assert self.labeled_structure.reduced_formula == "Si"
        assert self.propertied_structure.reduced_formula == "Si"
        assert self.V2O3.reduced_formula == "V2O3"

    def test_elements(self):
        assert self.struct.elements == [Element("Si")]
        assert self.propertied_structure.elements == [Element("Si")]

    def test_chemical_system(self):
        assert self.struct.chemical_system == "Si"
        assert self.propertied_structure.chemical_system == "Si"
        assert self.V2O3.chemical_system == "O-V"

    def test_chemical_system_set(self):
        assert self.struct.chemical_system_set == {"Si"}
        assert self.propertied_structure.chemical_system_set == {"Si"}
        assert self.V2O3.chemical_system_set == {"O", "V"}

    def test_specie_init(self):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        struct = IStructure(self.lattice, [{Species("O", -2): 1.0}, {Species("Mg", 2): 0.8}], coords)
        assert struct.formula == "Mg0.8 O1"

    def test_get_sorted_structure(self):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        struct = IStructure(self.lattice, ["O", "Li"], coords, site_properties={"charge": [-2, 1]})
        sorted_s = struct.get_sorted_structure()
        assert sorted_s[0].species == Composition("Li")
        assert sorted_s[1].species == Composition("O")
        assert sorted_s[0].charge == 1
        assert sorted_s[1].charge == -2
        struct = IStructure(
            self.lattice,
            ["Se", "C", "Se", "C"],
            [[0] * 3, [0.5] * 3, [0.25] * 3, [0.75] * 3],
        )
        assert [site.specie.symbol for site in struct.get_sorted_structure()] == [
            "C",
            "C",
            "Se",
            "Se",
        ]

    def test_get_space_group_data(self):
        assert self.struct.get_space_group_info() == ("Fd-3m", 227)

    def test_fractional_occupations(self):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        struct = IStructure(self.lattice, [{"O": 1.0}, {"Mg": 0.8}], coords)
        assert struct.formula == "Mg0.8 O1"
        assert not struct.is_ordered

    def test_labeled_structure(self):
        assert self.labeled_structure.labels == ["Si1", "Si2"]
        assert self.struct.labels == ["Si", "Si"]

    def test_get_distance(self):
        assert self.struct.get_distance(0, 1) == approx(2.35, abs=1e-2), "Distance calculated wrongly!"
        pt = [0.9, 0.9, 0.8]
        assert self.struct[0].distance_from_point(pt) == approx(1.50332963784, abs=1e-2), "Distance calculated wrongly!"

    def test_as_dict(self):
        si = Species("Si", 4)
        mn = Element("Mn")
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        struct = IStructure(self.lattice, [{si: 0.5, mn: 0.5}, {si: 0.5}], coords)
        assert "lattice" in struct.as_dict()
        assert "sites" in struct.as_dict()
        dct = self.propertied_structure.as_dict()
        assert dct["sites"][0]["properties"]["magmom"] == 5
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        struct = IStructure(
            self.lattice,
            [{Species("O", -2, spin=3): 1.0}, {Species("Mg", 2, spin=2): 0.8}],
            coords,
            site_properties={"magmom": [5, -5]},
            properties={"general_property": "test"},
        )
        dct = struct.as_dict()
        assert dct["sites"][0]["properties"]["magmom"] == 5
        assert dct["sites"][0]["species"][0]["spin"] == 3
        assert dct["properties"]["general_property"] == "test"

        dct = struct.as_dict(0)
        assert "volume" not in dct["lattice"]
        assert "xyz" not in dct["sites"][0]
        assert "label" not in dct["sites"][0]

    def test_from_dict(self):
        dct = self.propertied_structure.as_dict()
        struct = IStructure.from_dict(dct)
        assert struct[0].magmom == 5
        dct = self.propertied_structure.as_dict(0)
        s2 = IStructure.from_dict(dct)
        assert struct == s2

        dct = {
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
                            "spin": 3,
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
                            "spin": 2,
                            "element": "Mg",
                        }
                    ],
                    "label": "Mg2+:0.800",
                    "xyz": [3.8401979336749994, 1.2247250003039056e-06, 2.351631795225],
                },
            ],
            "properties": {"test_property": "test"},
        }
        struct = IStructure.from_dict(dct)
        assert struct[0].magmom == 5
        assert struct[0].specie.spin == 3
        assert struct.properties["test_property"] == "test"
        assert isinstance(struct, IStructure)

    def test_site_properties(self):
        site_props = self.propertied_structure.site_properties
        assert site_props["magmom"] == [5, -5]
        assert self.propertied_structure[0].magmom == 5
        assert self.propertied_structure[1].magmom == -5

    def test_properties_dict(self):
        assert self.propertied_structure.properties == {"test_property": "test"}

    def test_copy(self):
        new_struct = self.propertied_structure.copy(
            site_properties={"charge": [2, 3]}, properties={"another_prop": "test"}
        )
        assert new_struct[0].magmom == 5
        assert new_struct[1].magmom == -5
        assert new_struct[0].charge == 2
        assert new_struct[1].charge == 3
        assert new_struct.properties["another_prop"] == "test"
        assert new_struct.properties["test_property"] == "test"

        coords = [[0, 0, 0], [0.0, 0, 1e-07]]

        struct = IStructure(
            self.lattice,
            ["O", "Si"],
            coords,
            site_properties={"magmom": [5, -5]},
            properties={"test_property": "test"},
        )

        new_struct = struct.copy(
            site_properties={"charge": [2, 3]},
            sanitize=True,
            properties={"another_prop": "test"},
        )
        assert new_struct[0].magmom == -5
        assert new_struct[1].magmom == 5
        assert new_struct[0].charge == 3
        assert new_struct[1].charge == 2
        assert new_struct.volume == approx(struct.volume)
        assert new_struct.properties["another_prop"] == "test"
        assert new_struct.properties["test_property"] == "test"

    def test_interpolate(self):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        struct = IStructure(self.lattice, ["Si"] * 2, coords)
        coords2 = []
        coords2.extend(([0, 0, 0], [0.5, 0.5, 0.5]))
        struct2 = IStructure(self.struct.lattice, ["Si"] * 2, coords2)
        interpolated_structs = struct.interpolate(struct2, 10)
        for inter_struct in interpolated_structs:
            assert inter_struct is not None, "Interpolation Failed!"
            assert interpolated_structs[0].lattice == inter_struct.lattice
        assert_allclose(interpolated_structs[1][1].frac_coords, [0.725, 0.5, 0.725])

        # test ximages
        interpolated_structs = struct.interpolate(struct2, nimages=np.linspace(0.0, 1.0, 3))
        for inter_struct in interpolated_structs:
            assert inter_struct is not None, "Interpolation Failed!"
            assert interpolated_structs[0].lattice == inter_struct.lattice
        assert_allclose(interpolated_structs[1][1].frac_coords, [0.625, 0.5, 0.625])

        bad_lattice = np.eye(3)
        struct2 = IStructure(bad_lattice, ["Si"] * 2, coords2)
        with pytest.raises(ValueError, match="Structures with different lattices"):
            struct.interpolate(struct2)

        coords2 = []
        coords2.append([0, 0, 0])
        coords2.append([0.5, 0.5, 0.5])
        struct2 = IStructure(self.struct.lattice, ["Si", "Fe"], coords2)
        with pytest.raises(ValueError, match="Different species"):
            struct.interpolate(struct2)

        # Test autosort feature.
        s1 = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3), ["Fe"], [[0, 0, 0]])
        s1.pop(0)
        s2 = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3), ["Fe"], [[0, 0, 0]])
        s2.pop(2)
        rng = np.random.default_rng()
        rng.shuffle(s2)

        for struct in s1.interpolate(s2, autosort_tol=0.5):
            assert_allclose(s1[0].frac_coords, struct[0].frac_coords)
            assert_allclose(s1[2].frac_coords, struct[2].frac_coords)

        # Make sure autosort has no effect on simpler interpolations,
        # and with shuffled sites.
        s1 = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3), ["Fe"], [[0, 0, 0]])
        s2 = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3), ["Fe"], [[0, 0, 0]])
        s2[0] = "Fe", [0.01, 0.01, 0.01]
        rng.shuffle(s2)

        for struct in s1.interpolate(s2, autosort_tol=0.5):
            assert_allclose(s1[1].frac_coords, struct[1].frac_coords)
            assert_allclose(s1[2].frac_coords, struct[2].frac_coords)
            assert_allclose(s1[3].frac_coords, struct[3].frac_coords)

        # Test non-hexagonal setting.
        lattice = Lattice.rhombohedral(4.0718, 89.459)
        species = [{"S": 1.0}, {"Ni": 1.0}]
        coordinate = [(0.252100, 0.252100, 0.252100), (0.500000, 0.244900, -0.244900)]
        struct = Structure.from_spacegroup("R32:R", lattice, species, coordinate)
        assert struct.formula == "Ni3 S2"

        # test pbc
        coords = [[0.85, 0.85, 0.85]]
        coords2 = [[0.25, 0.25, 0.25]]
        struct_pbc = IStructure(self.lattice_pbc, ["Si"], coords)
        struct2_pbc = IStructure(self.lattice_pbc, ["Si"], coords2)
        int_s_pbc = struct_pbc.interpolate(struct2_pbc, nimages=2)
        assert_allclose(int_s_pbc[1][0].frac_coords, [1.05, 1.05, 0.55])

        # Test end_amplitude =/= 1
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        struct = IStructure(self.lattice, ["Si"] * 2, coords)
        coords2 = []
        coords2.extend(([0, 0, 0], [0.5, 0.5, 0.5]))
        struct2 = IStructure(self.struct.lattice, ["Si"] * 2, coords2)
        # testing large positive values
        interpolated_structs = struct.interpolate(struct2, 20, end_amplitude=2)
        for inter_struct in interpolated_structs:
            assert inter_struct is not None, "Interpolation Failed!"
            assert interpolated_structs[0].lattice == inter_struct.lattice
        assert_allclose(interpolated_structs[0][1].frac_coords, [0.75, 0.5, 0.75])
        assert_allclose(interpolated_structs[10][1].frac_coords, [0.5, 0.5, 0.5])
        assert_allclose(interpolated_structs[20][1].frac_coords, [0.25, 0.5, 0.25])
        # testing large negative values
        interpolated_structs = struct.interpolate(struct2, 20, end_amplitude=-2)
        for inter_struct in interpolated_structs:
            assert inter_struct is not None, "Interpolation Failed!"
            assert interpolated_structs[0].lattice == inter_struct.lattice
        assert_allclose(interpolated_structs[0][1].frac_coords, [0.75, 0.5, 0.75])
        assert_allclose(interpolated_structs[10][1].frac_coords, [1.0, 0.5, 1.0])
        assert_allclose(interpolated_structs[20][1].frac_coords, [1.25, 0.5, 1.25])
        # testing partial interpolation
        interpolated_structs = struct.interpolate(struct2, 5, end_amplitude=-0.5)
        for inter_struct in interpolated_structs:
            assert inter_struct is not None, "Interpolation Failed!"
            assert interpolated_structs[0].lattice == inter_struct.lattice
        assert_allclose(interpolated_structs[0][1].frac_coords, [0.75, 0.5, 0.75])
        assert_allclose(interpolated_structs[5][1].frac_coords, [0.875, 0.5, 0.875])
        # testing end_amplitude=0
        interpolated_structs = struct.interpolate(struct2, 5, end_amplitude=0)
        for inter_struct in interpolated_structs:
            assert inter_struct is not None, "Interpolation Failed!"
            assert interpolated_structs[0].lattice == inter_struct.lattice
            assert_allclose(inter_struct[1].frac_coords, [0.75, 0.5, 0.75])

    def test_interpolate_lattice(self):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        struct = IStructure(self.lattice, ["Si"] * 2, coords)
        coords2 = [[0, 0, 0], [0.5, 0.5, 0.5]]
        l2 = Lattice.from_parameters(3, 4, 4, 100, 100, 70)
        struct2 = IStructure(l2, ["Si"] * 2, coords2)
        int_s = struct.interpolate(struct2, 2, interpolate_lattices=True)
        assert_allclose(struct.lattice.abc, int_s[0].lattice.abc)
        assert_allclose(struct.lattice.angles, int_s[0].lattice.angles)
        assert_allclose(struct2.lattice.abc, int_s[2].lattice.abc)
        assert_allclose(struct2.lattice.angles, int_s[2].lattice.angles)
        int_angles = [110.3976469, 94.5359731, 64.5165856]
        assert_allclose(int_angles, int_s[1].lattice.angles)
        # Assert that volume is monotonic
        assert struct2.volume >= int_s[1].volume
        assert int_s[1].volume >= struct.volume

        # Repeat for end_amplitude = 0.5
        int_s = struct.interpolate(struct2, 2, interpolate_lattices=True, end_amplitude=0.5)
        assert_allclose(struct.lattice.abc, int_s[0].lattice.abc)
        assert_allclose(struct.lattice.angles, int_s[0].lattice.angles)
        assert_allclose(int_angles, int_s[2].lattice.angles)
        # Assert that volume is monotonic
        assert struct2.volume >= int_s[1].volume
        assert int_s[1].volume >= struct.volume

        # Repeat for end_amplitude = 2
        int_s = struct.interpolate(struct2, 4, interpolate_lattices=True, end_amplitude=2)
        assert_allclose(struct.lattice.abc, int_s[0].lattice.abc)
        assert_allclose(struct.lattice.angles, int_s[0].lattice.angles)
        assert_allclose(int_angles, int_s[1].lattice.angles)
        # Assert that volume is monotonic
        assert struct2.volume >= int_s[1].volume
        assert int_s[1].volume >= struct.volume

        # Repeat for end_amplitude = -1
        int_s = struct.interpolate(struct2, 2, interpolate_lattices=True, end_amplitude=-1)
        assert_allclose(struct.lattice.abc, int_s[0].lattice.abc)
        assert_allclose(struct.lattice.angles, int_s[0].lattice.angles)
        int_angles = [127.72010946461334, 86.27613506707404, 56.52554566317311]
        assert_allclose(int_angles, int_s[1].lattice.angles)
        # Assert that volume is monotonic (should be shrinking for negative end_amplitude)
        assert int_s[1].volume <= struct.volume
        # Assert that coordinate shift is reversed
        assert_allclose(int_s[1][1].frac_coords, [0.875, 0.5, 0.875])
        assert_allclose(int_s[2][1].frac_coords, [1.0, 0.5, 1.0])

    def test_interpolate_lattice_rotation(self):
        l1 = Lattice(np.eye(3))
        l2 = Lattice(np.diag((-1.01, -1.01, 1)))
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        struct1 = IStructure(l1, ["Si"] * 2, coords)
        struct2 = IStructure(l2, ["Si"] * 2, coords)

        # Test positive end_amplitudes
        for end_amplitude in [0, 0.5, 1, 2]:
            int_s = struct1.interpolate(struct2, 2, interpolate_lattices=True, end_amplitude=end_amplitude)
            # Assert that volume is monotonic
            assert struct2.volume >= int_s[1].volume
            assert int_s[1].volume >= struct1.volume

        # Test negative end_amplitudes
        for end_amplitude in [-2, -0.5, 0]:
            int_s = struct1.interpolate(struct2, 2, interpolate_lattices=True, end_amplitude=end_amplitude)
            # Assert that volume is monotonic
            assert struct2.volume >= int_s[1].volume
            assert int_s[1].volume <= struct1.volume

    def test_get_primitive_structure(self):
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]]
        fcc_ag = IStructure(Lattice.cubic(4.09), ["Ag"] * 4, coords)
        assert len(fcc_ag.get_primitive_structure()) == 1
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        bcc_li = IStructure(Lattice.cubic(4.09), ["Li"] * 2, coords)
        bcc_prim = bcc_li.get_primitive_structure()
        assert len(bcc_prim) == 1
        assert bcc_prim.lattice.alpha == approx(109.47122)
        bcc_li = IStructure(Lattice.cubic(4.09), ["Li"] * 2, coords, site_properties={"magmom": [1, -1]})
        bcc_prim = bcc_li.get_primitive_structure()
        assert len(bcc_prim) == 1
        assert bcc_prim.lattice.alpha == approx(109.47122)
        bcc_prim = bcc_li.get_primitive_structure(use_site_props=True)
        assert len(bcc_prim) == 2
        assert bcc_prim.lattice.alpha == approx(90)

        coords = [[0] * 3, [0.5] * 3, [0.25] * 3, [0.26] * 3]
        struct = IStructure(Lattice.cubic(4.09), ["Ag"] * 4, coords)
        assert len(struct.get_primitive_structure()) == 4

    def test_primitive_cell_site_merging(self):
        lattice = Lattice.cubic(10)
        coords = [[0, 0, 0], [0, 0, 0.5], [0, 0, 0.26], [0, 0, 0.74]]
        sp = ["Ag", "Ag", "Be", "Be"]
        struct = Structure(lattice, sp, coords)
        dm = struct.get_primitive_structure().distance_matrix
        assert_allclose(dm, [[0, 2.5], [2.5, 0]])

    def test_primitive_on_large_supercell(self):
        coords = [[0, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]]
        fcc_ag = Structure(Lattice.cubic(4.09), ["Ag"] * 4, coords)
        fcc_ag.make_supercell([2, 2, 2])
        fcc_ag_prim = fcc_ag.get_primitive_structure()
        assert len(fcc_ag_prim) == 1
        assert fcc_ag_prim.volume == approx(17.10448225)

    def test_primitive_with_constrained_lattice(self):
        struct = Structure.from_file(f"{TEST_FILES_DIR}/core/structure/Fe310.json.gz")
        constraints = {"a": 2.83133, "b": 4.69523, "gamma": 107.54840}
        prim_struct = struct.get_primitive_structure(constrain_latt=constraints)
        assert {key: getattr(prim_struct.lattice, key) for key in constraints} == approx(constraints)

    def test_primitive_with_similar_constraints(self):
        struct = Structure.from_file(f"{TEST_FILES_DIR}/core/structure/Fe310.json.gz")
        constraints = {"a": 2.83133, "b": 4.69523, "gamma": 107.54840}
        for perturb in (0, 1e-5, -1e-5):
            copy = constraints.copy()
            copy["a"] += perturb
            prim_struct = struct.get_primitive_structure(constrain_latt=constraints)
            copy_struct = struct.get_primitive_structure(constrain_latt=copy)
            assert prim_struct == copy_struct

    def test_primitive_positions(self):
        coords = [[0, 0, 0], [0.3, 0.35, 0.45]]
        struct = Structure(Lattice.from_parameters(1, 2, 3, 50, 66, 88), ["Ag"] * 2, coords)

        c = [[2, 0, 0], [1, 3, 0], [1, 1, 1]]

        for sc_matrix in [c]:
            sc = struct.copy()
            sc.make_supercell(sc_matrix)
            prim = sc.get_primitive_structure(0.01)

            assert len(prim) == 2
            assert prim.distance_matrix[0, 1] == approx(1.0203432356739286)

    def test_primitive_structure_volume_check(self):
        lattice = Lattice.tetragonal(10, 30)
        coords = [
            [0.5, 0.8, 0],
            [0.5, 0.2, 0],
            [0.5, 0.8, 0.333],
            [0.5, 0.5, 0.333],
            [0.5, 0.5, 0.666],
            [0.5, 0.2, 0.666],
        ]
        struct = IStructure(lattice, ["Ag"] * 6, coords)
        primitive = struct.get_primitive_structure(tolerance=0.1)
        assert len(primitive) == 6

    def test_get_miller_index(self):
        """Test for get miller index convenience method."""
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
        struct = self.struct
        nn = struct.get_neighbors_in_shell(struct[0].frac_coords, 2, 4, include_index=True, include_image=True)
        assert len(nn) == 47
        rand_radius = np.random.default_rng().uniform(3, 6)
        all_nn = struct.get_all_neighbors(rand_radius, include_index=True, include_image=True)
        for idx, site in enumerate(struct):
            assert len(all_nn[idx][0]) == 4
            assert len(all_nn[idx]) == len(struct.get_neighbors(site, rand_radius))

        for site, nns in zip(struct, all_nn, strict=True):
            for nn in nns:
                assert nn[0].is_periodic_image(struct[nn[2]])
                dist = sum((site.coords - nn[0].coords) ** 2) ** 0.5
                assert dist == approx(nn[1])

        struct = Structure(Lattice.cubic(1), ["Li"], [[0, 0, 0]])
        struct.make_supercell([2, 2, 2])
        assert sum(map(len, struct.get_all_neighbors(3))) == 976

        all_nn = struct.get_all_neighbors(0.05)
        assert [len(nn) for nn in all_nn] == [0] * len(struct)

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
        struct = Structure.from_str(poscar, fmt="poscar")
        site0 = struct[1]
        site1 = struct[9]
        neigh_sites = struct.get_neighbors(site0, 2.0)
        assert len(neigh_sites) == 1
        neigh_sites = struct.get_neighbors(site1, 2.0)
        assert len(neigh_sites) == 1

    def test_get_neighbor_list(self):
        # test mutable structure after applying strain (which makes lattice.matrix array no longer contiguous)
        # https://github.com/materialsproject/pymatgen/pull/3108
        mutable_struct = Structure.from_sites(self.struct)
        mutable_struct.apply_strain(0.01)
        for struct in (self.struct, mutable_struct):
            cy_indices1, cy_indices2, cy_offsets, cy_distances = struct.get_neighbor_list(3)
            py_indices1, py_indices2, py_offsets, py_distances = struct._get_neighbor_list_py(3)
            assert_allclose(cy_distances, py_distances)
            assert_allclose(cy_indices1, py_indices1)
            assert_allclose(cy_indices2, py_indices2)
            assert len(cy_offsets) == len(py_offsets)

    @pytest.mark.xfail(reason="TODO: need someone to fix this")
    @pytest.mark.skipif(not os.getenv("CI"), reason="Only run this in CI tests")
    def test_get_all_neighbors_crosscheck_old(self):
        rng = np.random.default_rng()
        for i in range(100):
            alpha, beta = rng.random(2) * 90
            a, b, c = 3 + rng.random(3) * 5
            species = ["H"] * 5
            frac_coords = rng.random(5, 3)
            try:
                lattice = Lattice.from_parameters(a, b, c, alpha, beta, 90)
                struct = Structure.from_spacegroup("P1", lattice, species, frac_coords)
                for nn_new, nn_old in zip(
                    struct.get_all_neighbors(4),
                    struct.get_all_neighbors_old(4),
                    strict=False,
                ):
                    sites1 = [i[0] for i in nn_new]
                    sites2 = [i[0] for i in nn_old]
                    assert set(sites1) == set(sites2)
                break
            except Exception:
                pass
        else:
            raise ValueError("No valid structure tested.")

        from pymatgen.electronic_structure.core import Spin

        d = {
            "@module": "pymatgen.core.structure",
            "@class": "Structure",
            "charge": None,
            "lattice": {
                "matrix": [
                    [0.0, 0.0, 5.5333],
                    [5.7461, 0.0, 3.518471486290303e-16],
                    [-4.692662837312786e-16, 7.6637, 4.692662837312786e-16],
                ],
                "a": 5.5333,
                "b": 5.7461,
                "c": 7.6637,
                "alpha": 90.0,
                "beta": 90.0,
                "gamma": 90.0,
                "volume": 243.66653780778103,
            },
            "sites": [
                {
                    "species": [
                        {
                            "element": "Mn",
                            "oxidation_state": 0,
                            "spin": Spin.down,
                            "occu": 1,
                        }
                    ],
                    "abc": [0.0, 0.5, 0.5],
                    "xyz": [2.8730499999999997, 3.83185, 4.1055671618015446e-16],
                    "label": "Mn0+,spin=-1",
                    "properties": {},
                },
                {
                    "species": [{"element": "Mn", "oxidation_state": None, "occu": 1.0}],
                    "abc": [1.232595164407831e-32, 0.5, 0.5],
                    "xyz": [2.8730499999999997, 3.83185, 4.105567161801545e-16],
                    "label": "Mn",
                    "properties": {},
                },
            ],
        }
        struct = Structure.from_dict(d)
        assert {i[0] for i in struct.get_neighbors(struct[0], 0.05)} == {
            i[0] for i in struct.get_neighbors_old(struct[0], 0.05)
        }

    def test_get_symmetric_neighbor_list(self):
        # tetragonal group with all bonds related by symmetry
        struct = Structure.from_spacegroup(100, [[1, 0, 0], [0, 1, 0], [0, 0, 2]], ["Fe"], [[0.0, 0.0, 0.0]])
        c_indices, p_indices, offsets, distances, s_indices, sym_ops = struct.get_symmetric_neighbor_list(0.8, sg=100)
        assert len(c_indices) == len(p_indices) == len(offsets) == len(distances) == 8
        assert_array_equal(c_indices, [0, 1, 1, 1, 0, 0, 0, 0])
        assert len(np.unique(s_indices)) == 1
        assert s_indices[0] == 0
        assert all(~np.isnan(s_indices))
        assert (sym_ops[0].affine_matrix == np.eye(4)).all()
        # now more complicated example with bonds of same length but with different symmetry
        s2 = Structure.from_spacegroup(
            198,
            [[8.908, 0, 0], [0, 8.908, 0], [0, 0, 8.908]],
            ["Cu"],
            [[0.0, 0.0, 0.0]],
        )
        c_indices2, p_indices2, offsets2, distances2, s_indices2, sym_ops2 = s2.get_symmetric_neighbor_list(7, sg=198)
        assert len(np.unique(s_indices2)) == 2
        assert len(s_indices2) == 48
        assert len(s_indices2[s_indices2 == 0]) == len(s_indices2[s_indices2 == 1])
        assert s_indices2[0] == 0
        assert s_indices2[24] == 1
        assert np.isclose(distances2[0], distances2[24])
        assert (sym_ops2[0].affine_matrix == np.eye(4)).all()
        assert (sym_ops2[24].affine_matrix == np.eye(4)).all()
        from_a2 = s2[c_indices2[0]].frac_coords
        to_a2 = s2[p_indices2[0]].frac_coords
        r_a2 = offsets2[0]
        from_b2 = s2[c_indices2[1]].frac_coords
        to_b2 = s2[p_indices2[1]].frac_coords
        r_b2 = offsets2[1]
        assert sym_ops2[1].are_symmetrically_related_vectors(from_a2, to_a2, r_a2, from_b2, to_b2, r_b2)
        assert sym_ops2[1].are_symmetrically_related_vectors(from_b2, to_b2, r_b2, from_a2, to_a2, r_a2)
        c_indices3, p_indices3, offsets3, distances3, s_indices3, _sym_ops3 = s2.get_symmetric_neighbor_list(
            7, sg=198, unique=True
        )
        assert all(np.sort(np.array([c_indices3, p_indices3]).flatten()) == np.sort(c_indices2))
        assert all(np.sort(np.array([c_indices3, p_indices3]).flatten()) == np.sort(p_indices2))
        assert len(offsets3) == len(distances3) == len(s_indices3) == 24

    def test_get_all_neighbors_outside_cell(self):
        struct = Structure(
            Lattice.cubic(2),
            ["Li", "Li", "Li", "Si"],
            [[3.1] * 3, [0.11] * 3, [-1.91] * 3, [0.5] * 3],
        )
        all_nn = struct.get_all_neighbors(0.2, include_index=True)
        for site, nns in zip(struct, all_nn, strict=True):
            for nn in nns:
                assert nn[0].is_periodic_image(struct[nn[2]])
                d = sum((site.coords - nn[0].coords) ** 2) ** 0.5
                assert d == approx(nn[1])
        assert list(map(len, all_nn)) == [2, 2, 2, 0]

    def test_get_all_neighbors_small_cutoff(self):
        struct = Structure(
            Lattice.cubic(2),
            ["Li", "Li", "Li", "Si"],
            [[3.1] * 3, [0.11] * 3, [-1.91] * 3, [0.5] * 3],
        )
        all_nn = struct.get_all_neighbors(1e-5, include_index=True)
        assert len(all_nn) == len(struct)
        assert all_nn[0] == []

        all_nn = struct.get_all_neighbors(0, include_index=True)
        assert len(all_nn) == len(struct)
        assert all_nn[0] == []

    def test_coincide_sites(self):
        struct = Structure(
            Lattice.cubic(5),
            ["Li", "Li", "Li"],
            [[0.1, 0.1, 0.1], [0.1, 0.1, 0.1], [3, 3, 3]],
            coords_are_cartesian=True,
        )
        all_nn = struct.get_all_neighbors(1e-5, include_index=True)
        assert [len(lst) for lst in all_nn] == [0, 0, 0]

    def test_get_all_neighbors_equal(self):
        struct = Structure(
            Lattice.cubic(2),
            ["Li", "Li", "Li", "Si"],
            [[3.1] * 3, [0.11] * 3, [-1.91] * 3, [0.5] * 3],
        )
        with pytest.warns(FutureWarning, match="get_all_neighbors_old is deprecated"):
            nn_traditional = struct.get_all_neighbors_old(4, include_index=True, include_image=True, include_site=True)

        nn_cell_lists = struct.get_all_neighbors(4, include_index=True, include_image=True)

        for idx in range(4):
            assert len(nn_traditional[idx]) == len(nn_cell_lists[idx])
            norm = np.linalg.norm(
                np.array(sorted(j[1] for j in nn_traditional[idx])) - np.array(sorted(j[1] for j in nn_cell_lists[idx]))
            )
            assert norm < 1e-3

    def test_get_dist_matrix(self):
        assert_allclose(self.struct.distance_matrix, [[0.0, 2.3516318], [2.3516318, 0.0]])

    def test_to_from_file_and_string(self):
        for fmt in ("cif", "json", "poscar", "cssr", "pwmat"):
            struct = self.struct.to(fmt=fmt)
            assert struct is not None
            ss = IStructure.from_str(struct, fmt=fmt)
            assert_allclose(ss.lattice.parameters, self.struct.lattice.parameters, atol=1e-5)
            assert_allclose(ss.frac_coords, self.struct.frac_coords)
            assert isinstance(ss, IStructure)

        assert "Fd-3m" in self.struct.to(fmt="CIF", symprec=0.1)

        poscar_path = f"{self.tmp_path}/POSCAR.testing"
        poscar_str = self.struct.to(filename=poscar_path)
        with open(poscar_path, encoding="utf-8") as file:
            assert file.read() == poscar_str
        assert Structure.from_file(poscar_path) == self.struct

        yaml_path = f"{self.tmp_path}/Si_testing.yaml"
        yaml_str = self.struct.to(filename=yaml_path)
        with open(yaml_path, encoding="utf-8") as file:
            assert file.read() == yaml_str
        assert Structure.from_file(yaml_path) == self.struct

        # Test Path support
        struct = Structure.from_file(Path(yaml_path))
        assert struct == self.struct

        # Test .yml extension works too.
        yml_path = yaml_path.replace(".yaml", ".yml")
        os.replace(yaml_path, yml_path)
        assert Structure.from_file(yml_path) == self.struct

        atom_config_path = f"{self.tmp_path}/atom-test.config"
        self.struct.to(filename=atom_config_path)
        assert Structure.from_file(atom_config_path) == self.struct

        with pytest.raises(
            ValueError,
            match="Format not specified and could not infer from filename='whatever'",
        ):
            self.struct.to(filename="whatever")
        with pytest.raises(ValueError, match="Invalid fmt='badformat'"):
            self.struct.to(fmt="badformat")

        # Default as JSON (no exception expected)
        assert self.struct.to() == self.struct.to(fmt="json")

        self.struct.to(filename=(gz_json_path := "POSCAR.testing.gz"))
        struct = Structure.from_file(gz_json_path)
        assert struct == self.struct

        # test CIF file with unicode error
        # https://github.com/materialsproject/pymatgen/issues/2947
        struct = Structure.from_file(f"{TEST_FILES_DIR}/io/cif/mcif/bad-unicode-gh-2947.mcif")
        assert struct.formula == "Ni32 O32"

        # make sure CIfParser.parse_structures() and Structure.from_file() are consistent
        # i.e. uses same merge_tol for site merging, same primitive=False, etc.
        assert struct == CifParser(f"{TEST_FILES_DIR}/io/cif/mcif/bad-unicode-gh-2947.mcif").parse_structures()[0]

        # https://github.com/materialsproject/pymatgen/issues/3551
        json_path = Path("test-with-path.json")
        self.struct.to(filename=json_path)
        assert os.path.isfile(json_path)

    def test_to_file_alias(self):
        out_path = f"{self.tmp_path}/POSCAR"
        assert self.struct.to(out_path) == self.struct.to_file(out_path)
        assert os.path.isfile(out_path)

    def test_to_from_ase_atoms(self):
        pytest.importorskip("ase")

        atoms = self.struct.to_ase_atoms()
        assert isinstance(atoms, Atoms)
        assert len(atoms) == len(self.struct)

        structure = AseAtomsAdaptor.get_structure(atoms)
        assert structure == self.struct

        # Ensure PBC is `bool` type (not `np.bool_`) and JSON serializable
        assert "pymatgen.core.structure" in structure.to(fmt="json")

        assert IStructure.from_ase_atoms(atoms) == self.struct
        assert type(IStructure.from_ase_atoms(atoms)) is IStructure

    def test_pbc(self):
        assert self.struct.pbc == (True, True, True)
        assert self.struct.is_3d_periodic
        struct_pbc = Structure(self.lattice_pbc, ["Si"] * 2, self.struct.frac_coords)
        assert struct_pbc.pbc == (True, True, False)
        assert not struct_pbc.is_3d_periodic

    def test_sites_setter(self):
        struct = self.struct.copy()
        new_sites = struct.sites[::-1]  # reverse order of sites
        struct.sites = new_sites
        assert struct.sites == new_sites

    def test_get_symmetry_dataset(self):
        """Test getting symmetry dataset from structure using different backends."""
        # Test spglib backend
        for backend in ("spglib", "moyopy"):
            pytest.importorskip(
                backend,
            )
            dataset = self.struct.get_symmetry_dataset(backend=backend)
            assert isinstance(dataset, dict)
            assert dataset["number"] == 227
            assert dataset["international"] == "Fd-3m"
            assert len(dataset["orbits"]) == 2
            pytest.approx(dataset["std_origin_shift"], [0, 0, 0])

        dataset = self.struct.get_symmetry_dataset(backend="spglib", return_raw_dataset=True)

        assert dataset.number == 227  # Fd-3m space group
        assert dataset.international == "Fd-3m"
        assert len(dataset.rotations) > 0
        assert len(dataset.translations) > 0

        # Test moyopy backend if available
        moyopy = pytest.importorskip("moyopy")
        dataset = self.struct.get_symmetry_dataset(backend="moyopy", return_raw_dataset=True)
        assert isinstance(dataset, moyopy.MoyoDataset)
        assert dataset.prim_std_cell.numbers == [14, 14]  # Si atomic number is 14

        # Test import error
        with (
            mock.patch.dict("sys.modules", {"moyopy": None}),
            pytest.raises(ImportError, match="moyopy is not installed. Run pip install moyopy."),
        ):
            self.struct.get_symmetry_dataset(backend="moyopy")

        # Test invalid backend
        with pytest.raises(ValueError, match="Invalid backend='42'"):
            self.struct.get_symmetry_dataset(backend="42")


class TestStructure(MatSciTest):
    def setup_method(self):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = Lattice(
            [
                [3.8401979337, 0, 0],
                [1.9200989668, 3.3257101909, 0],
                [0, -2.2171384943, 3.1355090603],
            ]
        )
        self.struct = Structure(lattice, ["Si", "Si"], coords)
        self.struct.properties["foo"] = "bar"
        self.cu_structure = Structure(lattice, ["Cu", "Cu"], coords)
        self.disordered = Structure.from_spacegroup("Im-3m", Lattice.cubic(3), [Composition("Fe0.5Mn0.5")], [[0, 0, 0]])
        self.labeled_structure = Structure(lattice, ["Si", "Si"], coords, labels=["Si1", "Si2"])

    def test_calc_property(self):
        pytest.importorskip("matcalc")
        d = self.struct.calc_property("elasticity")
        assert "bulk_modulus_vrh" in d

    @pytest.mark.skipif(
        SETTINGS.get("PMG_MAPI_KEY", "") == "",
        reason="PMG_MAPI_KEY environment variable not set or MP API is down. This is also the case in a PR.",
    )
    def test_from_id(self):
        s = Structure.from_id("mp-1143")
        assert isinstance(s, Structure)
        assert s.reduced_formula == "Al2O3"

        try:
            website_down = requests.get("https://www.crystallography.net", timeout=60).status_code != 200
        except (
            requests.exceptions.ConnectionError,
            urllib3.exceptions.ConnectTimeoutError,
        ):
            website_down = True
        if not website_down:
            s = Structure.from_id("1101077", source="COD")
            assert s.reduced_formula == "LiV2O4"

    def test_mutable_sequence_methods(self):
        struct = self.struct
        struct[0] = "Fe"
        assert struct.formula == "Fe1 Si1"
        struct[0] = "Fe", [0.5, 0.5, 0.5]
        assert struct.formula == "Fe1 Si1"
        assert_allclose(struct[0].frac_coords, [0.5, 0.5, 0.5])
        struct.reverse()
        assert struct[0].specie == Element("Si")
        assert_allclose(struct[0].frac_coords, [0.75, 0.5, 0.75])
        struct[0] = {"Mn": 0.5}
        assert struct.formula == "Mn0.5 Fe1"
        del struct[1]
        assert struct.formula == "Mn0.5"
        struct[0] = "Fe", [0.9, 0.9, 0.9], {"magmom": 5}
        assert struct.formula == "Fe1"
        assert struct[0].magmom == 5

        # Test atomic replacement.
        struct["Fe"] = "Mn"
        assert struct.formula == "Mn1"

        # Test slice replacement.
        struct = MatSciTest.get_structure("Li2O")
        struct[:2] = "S"
        assert struct.formula == "Li1 S2"

    def test_not_hashable(self):
        with pytest.raises(TypeError, match="unhashable type: 'Structure'"):
            _ = {self.struct: 1}

    def test_sort(self):
        self.struct[0] = "F"
        returned = self.struct.sort()
        assert returned is self.struct
        assert self.struct[0].species_string == "Si"
        assert self.struct[1].species_string == "F"
        self.struct.sort(key=lambda site: site.species_string)
        assert self.struct[0].species_string == "F"
        assert self.struct[1].species_string == "Si"
        self.struct.sort(key=lambda site: site.species_string, reverse=True)
        assert self.struct[0].species_string == "Si"
        assert self.struct[1].species_string == "F"

    def test_replace(self):
        assert self.struct.formula == "Si2"
        struct = self.struct.replace(0, "O")
        assert struct is self.struct
        assert struct.formula == "Si1 O1"
        assert_allclose(struct[0].frac_coords, [0, 0, 0])
        struct.replace(0, "O", coords=[0.25, 0.25, 0.25])
        assert struct.formula == "Si1 O1"
        assert_allclose(struct[0].frac_coords, [0.25, 0.25, 0.25])
        struct.replace(0, "O", properties={"magmom": 1})
        assert struct.formula == "Si1 O1"
        assert struct[0].magmom == 1
        struct.replace(0, "O", properties={"magmom": 2}, coords=[0.9, 0.9, 0.9])
        assert struct.formula == "Si1 O1"
        assert struct[0].magmom == 2
        assert_allclose(struct[0].frac_coords, [0.9, 0.9, 0.9])

    def test_replace_species(self):
        assert self.struct.formula == "Si2"
        struct = self.struct.replace_species({"Si": "Na"})
        assert struct is self.struct
        assert struct.formula == "Na2"

        # test replacement with a dictionary
        struct.replace_species({"Na": {"K": 0.75, "P": 0.25}})
        assert struct.formula == "K1.5 P0.5"

        # test replacement with species not present in structure
        with pytest.warns(
            UserWarning,
            match="Some species to be substituted are not present in structure.",
        ):
            struct.replace_species({"Li": "Na"})

        # test in_place=False
        new_struct = struct.replace_species({"K": "Li"}, in_place=False)
        assert struct.formula == "K1.5 P0.5"
        assert new_struct.formula == "Li1.5 P0.5"

    def test_replace_species_labels(self):
        """https://github.com/materialsproject/pymatgen/issues/3658"""
        struct = self.labeled_structure
        new1 = struct.replace_species({"Si": "Ge"}, in_place=False)
        assert new1.labels == ["Ge", "Ge"]

        replacement = {"Si": 0.5, "Ge": 0.5}
        label = ", ".join(f"{key}:{val:.3}" for key, val in replacement.items())
        new2 = struct.replace_species({"Si": replacement}, in_place=False)
        assert new2.labels == [label] * len(struct)

    def test_relabel_sites(self):
        struct = self.struct.copy()
        assert self.struct.labels == ["Si", "Si"]
        relabeled = struct.relabel_sites()
        assert relabeled is struct
        assert relabeled.labels == struct.labels == ["Si_1", "Si_2"]

        struct.replace_species({"Si": "Li"})
        assert struct.labels == ["Li", "Li"]
        struct.relabel_sites()
        assert struct.labels == ["Li_1", "Li_2"]

        li_si = self.struct.copy().replace(0, "Li")
        assert li_si.labels == ["Li", "Si"]
        li_si.relabel_sites(ignore_uniq=True)  # check no-op for unique labels
        assert li_si.labels == ["Li", "Si"]
        li_si.relabel_sites()  # check no-op for unique labels
        assert li_si.labels == ["Li_1", "Si_1"]

    def test_append_insert_remove_replace_substitute(self):
        struct = self.struct
        struct.insert(1, "O", [0.5, 0.5, 0.5])
        assert struct.formula == "Si2 O1"
        assert struct.n_elems == 2
        assert struct.symbol_set == ("O", "Si")
        assert struct.indices_from_symbol("Si") == (0, 2)
        assert struct.indices_from_symbol("O") == (1,)
        del struct[2]
        assert struct.formula == "Si1 O1"
        assert struct.indices_from_symbol("Si") == (0,)
        assert struct.indices_from_symbol("O") == (1,)
        struct.append("N", [0.25, 0.25, 0.25])
        assert struct.formula == "Si1 N1 O1"
        assert struct.n_elems == 3
        assert struct.symbol_set == ("N", "O", "Si")
        assert struct.indices_from_symbol("Si") == (0,)
        assert struct.indices_from_symbol("O") == (1,)
        assert struct.indices_from_symbol("N") == (2,)
        struct[0] = "Ge"
        assert struct.formula == "Ge1 N1 O1"
        assert struct.symbol_set == ("Ge", "N", "O")
        struct.replace_species({"Ge": "Si"})
        assert struct.formula == "Si1 N1 O1"
        assert struct.n_elems == 3

        struct.replace_species({"Si": {"Ge": 0.5, "Si": 0.5}})
        assert struct.formula == "Si0.5 Ge0.5 N1 O1"
        # this should change the .5Si .5Ge sites to .75Si .25Ge
        struct.replace_species({"Ge": {"Ge": 0.5, "Si": 0.5}})
        assert struct.formula == "Si0.75 Ge0.25 N1 O1"

        assert struct.n_elems == 4

        struct.replace_species({"Ge": "Si"})
        substituted = struct.substitute(1, "hydroxyl")
        assert substituted is struct
        assert struct.formula == "Si1 H1 N1 O1"
        assert struct.symbol_set == ("H", "N", "O", "Si")
        with pytest.raises(
            ValueError,
            match="Can't find functional group 'OH' in list. Provide explicit coordinates instead",
        ):
            substituted = struct.substitute(2, "OH")
        # Distance between O and H
        assert struct.get_distance(2, 3) == approx(0.96)
        # Distance between Si and H
        assert struct.get_distance(0, 3) == approx(2.09840889)

        h_removed = struct.remove_species(["H"])
        assert h_removed is struct
        assert struct.formula == "Si1 N1 O1"

        sites_removed = struct.remove_sites([1, 2])
        assert sites_removed is struct
        assert struct.formula == "Si1"

    def test_add_remove_site_property(self):
        struct = self.struct
        returned = struct.add_site_property("charge", [4.1, -5])
        assert returned is struct
        assert struct[0].charge == approx(4.1)
        assert struct[1].charge == approx(-5)
        struct.add_site_property("magmom", [3, 2])
        assert struct[0].charge == approx(4.1)
        assert struct[0].magmom == 3
        returned = struct.remove_site_property("magmom")
        assert returned is struct
        with pytest.raises(AttributeError, match="attr='magmom' not found on PeriodicSite"):
            _ = struct[0].magmom

    def test_site_properties(self):
        # Make sure that site properties are set to None for missing values.
        self.struct.add_site_property("charge", [4.1, -5])
        self.struct.append("Li", [0.3, 0.3, 0.3])
        assert len(self.struct.site_properties["charge"]) == 3

        props = {"test_property": 42}
        struct_with_props = Structure(
            lattice=Lattice.cubic(3),
            species=("Fe", "Fe"),
            coords=((0, 0, 0), (0.5, 0.5, 0.5)),
            site_properties={"magmom": [5, -5]},
            properties=props,
        )

        dct = struct_with_props.as_dict()
        struct = Structure.from_dict(dct)
        assert struct.properties == props
        assert dct == struct.as_dict()

        json_str = struct_with_props.to(fmt="json")
        assert '"test_property":42' in json_str
        struct = Structure.from_str(json_str, fmt="json")
        assert struct.properties == props
        assert dct == struct.as_dict()

    def test_selective_dynamics(self):
        """Ensure selective dynamics as numpy arrays can be JSON serialized."""
        struct = self.get_structure("Li2O")
        struct.add_site_property(
            "selective_dynamics", np.array([[True, True, True], [False, False, False], [True, True, True]])
        )

        orjson_str = struct.to(fmt="json")

        # Also test round trip
        orjson_struct = Structure.from_str(orjson_str, fmt="json")
        assert struct == orjson_struct

        with pytest.raises(TypeError, match="Object of type ndarray is not JSON serializable"):
            # Use a dummy kwarg (default value) to force `json.dumps`
            struct.to(fmt="json", ensure_ascii=True)

    def test_perturb(self):
        struct = self.get_structure("Li2O")
        struct_orig = struct.copy()
        struct.perturb(0.1)
        # Ensure all sites were perturbed by a distance of at most 0.1 Angstroms
        for site, site_orig in zip(struct, struct_orig, strict=True):
            cart_dist = site.distance(site_orig)
            # allow 1e-6 to account for numerical precision
            assert cart_dist <= 0.1 + 1e-6, f"Distance {cart_dist} > 0.1"

        # Check that the perturbation does not result in the same translation
        vecs = [site.coords - site_orig.coords for site, site_orig in zip(struct, struct_orig, strict=True)]
        compare_vecs = [np.allclose(v1, v2) for v1, v2 in itertools.pairwise(vecs)]
        assert not all(compare_vecs)

        # Test that same seed gives same perturbation
        s1 = self.get_structure("Li2O")
        s2 = self.get_structure("Li2O")
        s1.perturb(0.1, seed=42)
        s2.perturb(0.1, seed=42)
        for site1, site2 in zip(s1, s2, strict=True):
            assert site1.distance(site2) < 1e-7  # should be exactly equal up to numerical precision

        # Test that different seeds give different perturbations
        s3 = self.get_structure("Li2O")
        s3.perturb(0.1, seed=100)
        any_different = False
        for site1, site3 in zip(s1, s3, strict=True):
            if site1.distance(site3) > 1e-7:
                any_different = True
                break
        assert any_different, "Different seeds should give different perturbations"

        # Test min_distance
        s4 = self.get_structure("Li2O")
        s4.perturb(0.1, min_distance=0.05, seed=42)
        any_different = False
        for site1, site4 in zip(s1, s4, strict=True):
            if site1.distance(site4) > 1e-7:
                any_different = True
                break
        assert any_different, "Using min_distance should give different perturbations"

    def test_add_oxidation_state_by_element(self):
        oxidation_states = {"Si": -4}
        returned = self.struct.add_oxidation_state_by_element(oxidation_states)
        assert returned is self.struct
        for site in self.struct:
            for specie in site.species:
                assert specie.oxi_state == oxidation_states[specie.symbol], "Wrong oxidation state assigned!"
        oxidation_states = {"Fe": 2}
        with pytest.raises(
            ValueError,
            match="Oxidation states not specified for all elements, missing={'Si'}",
        ):
            self.struct.add_oxidation_state_by_element(oxidation_states)

    def test_add_oxidation_states_by_site(self):
        returned = self.struct.add_oxidation_state_by_site([2, -4])
        assert returned is self.struct
        assert self.struct[0].specie.oxi_state == 2
        with pytest.raises(
            ValueError,
            match="Oxidation states of all sites must be specified, expected 2 values, got 1",
        ):
            self.struct.add_oxidation_state_by_site([1])

    def test_remove_oxidation_states(self):
        co_elem = Element("Co")
        o_elem = Element("O")
        co_specie = Species("Co", 2)
        o_specie = Species("O", -2)
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = Lattice.cubic(10)
        struct_elem = Structure(lattice, [co_elem, o_elem], coords)
        struct_specie = Structure(lattice, [co_specie, o_specie], coords)
        returned = struct_specie.remove_oxidation_states()
        assert returned is struct_specie
        assert struct_elem == struct_specie, "Oxidation state remover failed"

    def test_add_oxidation_state_by_guess(self):
        struct = MatSciTest.get_structure("Li2O")
        returned = struct.add_oxidation_state_by_guess()
        assert returned is struct
        expected = [Species("Li", 1), Species("O", -2)]
        for site in struct:
            assert site.specie in expected

    def test_add_remove_spin_states(self):
        lattice = Lattice.cubic(4.17)
        species = ["Ni", "O"]
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        nio = Structure.from_spacegroup(225, lattice, species, coords)

        # should do nothing, but not fail
        nio1 = nio.remove_spin()
        assert nio1 is nio

        spins = {"Ni": 5}
        nio2 = nio.add_spin_by_element(spins)
        assert nio2 is nio
        assert nio[0].specie.spin == 5, "Failed to add spin states"

        nio3 = nio.remove_spin()
        assert nio3 is nio
        assert nio[0].specie.spin is None

        spins = [5, -5, -5, 5, 0, 0, 0, 0]  # AFM on (001)
        nio4 = nio.add_spin_by_site(spins)
        assert nio4 is nio
        assert nio[1].specie.spin == -5, "Failed to add spin states"

        # Test lengths mismatch
        with pytest.raises(ValueError, match="Spins for all sites must be specified"):
            nio.add_spin_by_site([0, 1, 2])

    def test_apply_operation(self):
        op = SymmOp.from_axis_angle_and_translation([0, 0, 1], 90)
        struct = self.struct.copy()
        spg_info = struct.get_space_group_info()
        returned = struct.apply_operation(op)
        assert returned is struct
        assert_allclose(
            struct.lattice.matrix,
            [[0, 3.840198, 0], [-3.325710, 1.920099, 0], [2.217138, -0, 3.135509]],
            atol=1e-6,
        )
        assert returned.get_space_group_info() == spg_info

        op = SymmOp([[1, 1, 0, 0.5], [1, 0, 0, 0.5], [0, 0, 1, 0.5], [0, 0, 0, 1]])  # not a SymmOp of this struct
        struct = self.struct.copy()
        struct.apply_operation(op, fractional=True)
        assert_allclose(
            struct.lattice.matrix,
            [[5.760297, 3.325710, 0], [3.840198, 0, 0], [0, -2.217138, 3.135509]],
            5,
        )

        struct = self.struct.copy()
        # actual SymmOp of this struct: (SpacegroupAnalyzer(struct).get_symmetry_operations()[-2])
        op = SymmOp(
            [
                [0.0, 0.0, -1.0, 0.75],
                [-1.0, -1.0, 1.0, 0.5],
                [0.0, -1.0, 1.0, 0.75],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )
        struct.apply_operation(op, fractional=True)
        assert struct.get_space_group_info() == spg_info

        # same SymmOp in Cartesian coordinates:
        op = SymmOp(
            [
                [-0.5, -0.288675028, -0.816496280, 3.84019793],
                [-0.866025723, 0.166666176, 0.471404694, 0],
                [0, -0.942808868, 0.333333824, 2.35163180],
                [0, 0, 0, 1.0],
            ]
        )
        struct = self.struct.copy()
        struct.apply_operation(op, fractional=False)
        assert struct.get_space_group_info() == spg_info

    def test_apply_strain(self):
        struct = self.struct
        initial_coord = struct[1].coords
        returned = struct.apply_strain(0.01)
        assert returned is struct
        assert_allclose(
            struct.lattice.abc,
            (
                3.8785999130369997,
                3.878600984287687,
                3.8785999130549516,
            ),
        )
        assert_allclose(struct[1].coords, initial_coord * 1.01)
        a1, b1, c1 = struct.lattice.abc
        struct.apply_strain([0.1, 0.2, 0.3])
        a2, b2, c2 = struct.lattice.abc
        assert a2 / a1 == approx(1.1)
        assert b2 / b1 == approx(1.2)
        assert c2 / c1 == approx(1.3)

        # test inplace = False
        struct_copy = struct.copy()

        # Applying strain inplace=False
        strained_struct = struct.apply_strain(0.01, inplace=False)
        assert struct == struct_copy  # Should be equal
        assert struct is not strained_struct  # Should return a different object
        assert_allclose(strained_struct[1].coords, [4.407059, -0.169626, 3.118569], atol=1e-5)

    def test_scale_lattice(self):
        initial_coord = self.struct[1].coords
        returned = self.struct.scale_lattice(self.struct.volume * 1.01**3)
        assert returned is self.struct
        assert_allclose(
            self.struct.lattice.abc,
            (3.8785999130369997, 3.878600984287687, 3.8785999130549516),
        )
        assert_allclose(self.struct[1].coords, initial_coord * 1.01)

    def test_translate_sites(self):
        returned = self.struct.translate_sites([0, 1], [0.5, 0.5, 0.5], frac_coords=True)
        assert returned is self.struct
        assert_allclose(self.struct.frac_coords[0], [0.5, 0.5, 0.5])

        self.struct.translate_sites([0], [0.5, 0.5, 0.5], frac_coords=False)
        assert_allclose(self.struct.cart_coords[0], [3.38014845, 1.05428585, 2.06775453])

        self.struct.translate_sites([0], [0.5, 0.5, 0.5], frac_coords=True, to_unit_cell=False)
        assert_allclose(self.struct.frac_coords[0], [1.00187517, 1.25665291, 1.15946374])

        lattice_pbc = Lattice(self.struct.lattice.matrix, pbc=(True, True, False))
        struct_pbc = Structure(lattice_pbc, ["Si"], [[0.75, 0.75, 0.75]])
        struct_pbc.translate_sites([0], [0.5, 0.5, 0.5], frac_coords=True, to_unit_cell=True)
        assert_allclose(struct_pbc.frac_coords[0], [0.25, 0.25, 1.25])

        with pytest.raises(IndexError, match="list index out of range"):
            self.struct.translate_sites([5], [0.5, 0.5, 0.5])

        # test inverse operation leaves structure unchanged
        original_struct = self.struct.copy()
        self.struct.translate_sites([0], [0.5, 0.5, 0.5], frac_coords=True, to_unit_cell=False)
        self.struct.translate_sites([0], [-0.5, -0.5, -0.5], frac_coords=True, to_unit_cell=False)
        assert self.struct == original_struct

    def test_rotate_sites(self):
        returned = self.struct.rotate_sites(
            indices=[1],
            theta=2.0 * np.pi / 3.0,
            anchor=self.struct[0].coords,
            to_unit_cell=False,
        )
        assert returned is self.struct
        assert_allclose(self.struct.frac_coords[1], [-1.25, 1.5, 0.75], atol=1e-6)
        self.struct.rotate_sites(
            indices=[1],
            theta=2.0 * np.pi / 3.0,
            anchor=self.struct[0].coords,
            to_unit_cell=True,
        )
        assert_allclose(self.struct.frac_coords[1], [0.75, 0.5, 0.75], atol=1e-6)

        with pytest.raises(IndexError, match="list index out of range"):
            self.struct.rotate_sites([5], 2.0 * np.pi / 3.0, self.struct[0].coords, to_unit_cell=False)

    def test_mul(self):
        self.struct *= [2, 1, 1]
        assert self.struct.formula == "Si4"
        struct = [2, 1, 1] * self.struct
        assert struct.formula == "Si8"
        assert isinstance(struct, Structure)
        struct = self.struct * [[1, 0, 0], [2, 1, 0], [0, 0, 2]]
        assert struct.formula == "Si8"
        assert_allclose(struct.lattice.abc, [7.6803959, 17.5979979, 7.6803959])

    def test_make_supercell(self):
        supercell = self.struct.make_supercell([2, 1, 1])
        assert supercell.formula == "Si4"
        # test that make_supercell modified the original structure
        assert len(self.struct) == len(supercell)

        supercell.make_supercell([[1, 0, 0], [2, 1, 0], [0, 0, 1]])
        assert supercell.formula == "Si4"
        supercell.make_supercell(2)
        assert supercell.formula == "Si32"
        assert_allclose(supercell.lattice.abc, [15.360792, 35.195996, 7.680396], 5)

        # test in_place=False leaves original structure unchanged
        orig_len = len(self.struct)
        # test that make_supercell casts floats to ints
        supercell = self.struct.make_supercell([2.5, 1, 1], in_place=False)
        assert len(self.struct) == orig_len
        assert len(supercell) == 2 * orig_len

    def test_make_supercell_labeled(self):
        struct = self.labeled_structure.copy()
        struct.make_supercell([1, 1, 2])
        assert set(struct.labels) == {"Si1", "Si2"}

    def test_disordered_supercell_primitive_cell(self):
        lattice = Lattice.cubic(2)
        coords = [[0.5, 0.5, 0.5]]
        sp = [{"Si": 0.54738}]
        struct = Structure(lattice, sp, coords)
        # this supercell often breaks things
        struct.make_supercell([[0, -1, 1], [-1, 1, 0], [1, 1, 1]])
        assert len(struct.get_primitive_structure()) == 1

    def test_another_supercell(self):
        # this is included b/c for some reason the old algo was failing on it
        struct = self.struct.copy()
        struct.make_supercell([[0, 2, 2], [2, 0, 2], [2, 2, 0]])
        assert struct.formula == "Si32"
        struct = self.struct.copy()
        struct.make_supercell([[0, 2, 0], [1, 0, 0], [0, 0, 1]])
        assert struct.formula == "Si4"

    def test_as_from_dict(self):
        dct = self.struct.as_dict()
        s1 = Structure.from_dict(dct)
        assert isinstance(s1, Structure)

    def test_default_dict_attrs(self):
        dct = self.struct.as_dict()
        assert dct["charge"] == 0
        assert dct["properties"] == {"foo": "bar"}

    def test_to_from_abivars(self):
        """Test as_dict, from_dict with fmt == abivars."""
        # Properties are not supported if fmt="abivars" as its not a serialization protocol
        # but a format that allows one to get a dict with the abinit variables defining the structure.
        struct = self.struct.copy()
        struct.properties = {}
        dct = struct.as_dict(fmt="abivars")
        assert "properties" not in dct
        s2 = Structure.from_dict(dct, fmt="abivars")
        assert s2 == struct
        assert isinstance(s2, Structure)

    def test_to_from_file_str(self):
        # to/from string
        for fmt in (
            "cif",
            "json",
            "poscar",
            "cssr",
            "yaml",
            "yml",
            "xsf",
            "res",
            "pwmat",
        ):
            struct = self.struct.to(fmt=fmt)
            assert struct is not None
            ss = Structure.from_str(struct, fmt=fmt)
            assert_allclose(ss.lattice.parameters, self.struct.lattice.parameters, atol=1e-5)
            assert_allclose(ss.frac_coords, self.struct.frac_coords)
            assert isinstance(ss, Structure)

        # to/from file
        self.struct.to(filename=(poscar_path := "POSCAR.testing"))
        assert os.path.isfile(poscar_path)

        for ext in (".json", ".json.gz", ".json.bz2", ".json.xz", ".json.lzma"):
            self.struct.to(filename=f"json-struct{ext}")
            assert os.path.isfile(f"json-struct{ext}")
            assert Structure.from_file(f"json-struct{ext}") == self.struct

        # test Structure.from_file with unsupported file extension (using tmp JSON file with wrong ext)
        Path(filename := f"{self.tmp_path}/bad.extension").write_text(self.struct.to(fmt="json"), encoding="utf-8")
        with pytest.raises(ValueError, match="Unrecognized extension in filename="):
            self.struct.from_file(filename=filename)

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
            labels=["A", "B"],
        )
        assert sum(s2.site_properties["charge"]) == 0
        assert s2.labels == ["A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "B"]

        struct = Structure.from_spacegroup("Pm-3m", Lattice.cubic(3), ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        assert struct.formula == "Cs1 Cl1"

        with pytest.raises(
            ValueError,
            match="Supplied lattice with parameters \\(.+\\) is incompatible with supplied spacegroup Pm-3m",
        ):
            Structure.from_spacegroup(
                "Pm-3m",
                Lattice.tetragonal(1, 3),
                ["Cs", "Cl"],
                [[0, 0, 0], [0.5, 0.5, 0.5]],
            )

        with pytest.raises(
            ValueError,
            match=r"Supplied species and coords lengths \(1 vs 2\) are different",
        ):
            Structure.from_spacegroup(
                "Pm-3m",
                Lattice.cubic(3),
                ["Cs"],
                [[0, 0, 0], [0.5, 0.5, 0.5]],
            )

        struct = Structure.from_spacegroup(139, np.eye(3), ["H"], [[Fraction(1, 2), Fraction(1, 4), Fraction(0)]])
        assert len(struct) == 8

    def test_from_magnetic_spacegroup(self):
        # AFM MnF
        s1 = Structure.from_magnetic_spacegroup(
            "P4_2'/mnm'",
            Lattice.tetragonal(4.87, 3.30),
            ["Mn", "F"],
            [[0, 0, 0], [0.30, 0.30, 0]],
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
                [0, 0, 0.50],
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
        struct = Structure(Lattice.cubic(1), species, coords)
        struct.merge_sites(mode="sum")
        assert struct[0].specie.symbol == "Ag"
        assert struct[1].species == Composition({"Cl": 0.35, "F": 0.25})
        assert_allclose(struct[1].frac_coords, [0.5, 0.5, 0.5005])

        # Test illegal mode
        with pytest.raises(ValueError, match="Illegal mode='illegal', should start with a/d/s"):
            struct.merge_sites(mode="illegal")

        # Test for TaS2 with spacegroup 166 in 160 setting
        lattice = Lattice.hexagonal(3.374351, 20.308941)
        species = ["Ta", "S", "S"]
        coords = [
            [0, 0, 0.944333],
            [0.333333, 0.666667, 0.353424],
            [0.666667, 0.333333, 0.535243],
        ]
        struct_tas2 = Structure.from_spacegroup(160, lattice, species, coords)
        assert len(struct_tas2) == 13
        struct_tas2.merge_sites(mode="delete")
        assert len(struct_tas2) == 9

        lattice = Lattice.hexagonal(3.587776, 19.622793)
        species = ["Na", "V", "S", "S"]
        coords = [
            [0.333333, 0.666667, 0.165000],
            [0, 0, 0.998333],
            [0.333333, 0.666667, 0.399394],
            [0.666667, 0.333333, 0.597273],
        ]
        struct_navs2 = Structure.from_spacegroup(160, lattice, species, coords)
        assert len(struct_navs2) == 18
        struct_navs2.merge_sites(mode="delete")
        assert len(struct_navs2) == 12

        # Test that we can average the site properties that are numerical (float/int)
        lattice = Lattice.hexagonal(3.587776, 19.622793)
        species = ["Na", "V", "S", "S"]
        coords = [
            [0.333333, 0.666667, 0.165000],
            [0, 0, 0.998333],
            [0.333333, 0.666667, 0.399394],
            [0.666667, 0.333333, 0.597273],
        ]
        site_props = {"prop1": [3.0, 5.0, 7.0, 11.0]}
        struct_navs2 = Structure.from_spacegroup(160, lattice, species, coords, site_properties=site_props)
        struct_navs2.insert(0, "Na", coords[0], properties={"prop1": 100})  # int property
        struct_navs2.merge_sites(mode="average")
        assert len(struct_navs2) == 12
        assert any(math.isclose(site.properties["prop1"], 51.5) for site in struct_navs2)

        # Test non-numerical property warning
        struct_navs2.insert(0, "Na", coords[0], properties={"prop1": "hi"})
        with pytest.warns(UserWarning, match="But property is set to None"):
            struct_navs2.merge_sites(mode="average")

        # Test property handling for np.array (selective dynamics)
        poscar_str_0 = """Test POSCAR
1.0
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
1 1
Selective dynamics
direct
0.000000 0.000000 0.000000 T T T Si
0.750000 0.500000 0.750000 F F F O
"""
        poscar_str_1 = """offset a bit
1.0
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
1 1
Selective dynamics
direct
0.100000 0.000000 0.000000 T T T Si
0.750000 0.500000 0.750000 F F F O
"""

        struct_0 = Poscar.from_str(poscar_str_0).structure
        struct_1 = Poscar.from_str(poscar_str_1).structure

        for site in struct_0:
            struct_1.append(site.species, site.frac_coords, properties=site.properties)
        struct_1.merge_sites(mode="average")

        cu = Structure(
            Lattice.from_parameters(2.545584, 2.545584, 2.545584, 60, 60, 60), species=["Cu"], coords=[[0, 0, 0]]
        )
        cu.merge_sites(mode="delete")
        assert len(cu) == 1

    def test_properties(self):
        assert self.struct.num_sites == len(self.struct)
        self.struct.make_supercell(2)
        self.struct[1] = "C"
        sites = list(self.struct.group_by_types())
        assert sites[-1].specie.symbol == "C"
        self.struct.add_oxidation_state_by_element({"Si": 4, "C": 2})
        assert self.struct.charge == 62

    def test_species(self):
        assert {*map(str, self.struct.species)} == {"Si"}
        assert len(self.struct.species) == len(self.struct)

        # https://github.com/materialsproject/pymatgen/issues/3033
        with pytest.raises(AttributeError, match="species property only supports ordered structures!"):
            _ = self.disordered.species

    def test_set_item(self):
        struct = self.struct.copy()
        struct[0] = "C"
        assert struct.formula == "Si1 C1"
        struct[0, 1] = "Ge"
        assert struct.formula == "Ge2"
        struct[:2] = "Sn"
        assert struct.formula == "Sn2"

        struct = self.struct.copy()
        struct["Si"] = "C"
        assert struct.formula == "C2"
        struct["C"] = "C0.25Si0.5"
        assert struct.formula == "Si1 C0.5"
        struct["C"] = "C0.25Si0.5"
        assert struct.formula == "Si1.25 C0.125"

    def test_init_error(self):
        with pytest.raises(StructureError, match=r"len\(species\)=1 != len\(coords\)=2"):
            Structure(Lattice.cubic(3), ["Si"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    def test_from_sites(self):
        self.struct.add_site_property("hello", [1, 2])
        struct = Structure.from_sites(self.struct, to_unit_cell=True)
        assert struct.site_properties["hello"][1] == 2

        with pytest.raises(ValueError, match="You need at least 1 site to construct a Structure"):
            Structure.from_sites([])

    def test_charge(self):
        struct = Structure.from_sites(self.struct)
        assert struct.charge == 0, "Initial Structure not defaulting to behavior in SiteCollection"
        struct.add_oxidation_state_by_site([1, 1])
        assert struct.charge == 2, "Initial Structure not defaulting to behavior in SiteCollection"
        struct = Structure.from_sites(struct, charge=1)
        assert struct.charge == 1, "Overall charge not being stored in separate property"
        struct = struct.copy()
        assert struct.charge == 1, "Overall charge not being copied properly with no sanitization"
        struct = struct.copy(sanitize=True)
        assert struct.charge == 1, "Overall charge not being copied properly with sanitization"
        super_cell = struct * 3
        assert super_cell.charge == 27, "Overall charge is not being properly multiplied in IStructure __mul__"
        assert "Overall Charge: +1" in str(struct), "String representation not adding charge"
        sorted_s = super_cell.get_sorted_structure()
        assert sorted_s.charge == 27, "Overall charge is not properly copied during structure sorting"
        super_cell.set_charge(25)
        assert super_cell.charge == 25, "Set charge not properly modifying _charge"

    def test_vesta_lattice_matrix(self):
        silica_zeolite = Molecule.from_file(f"{TEST_FILES_DIR}/core/structure/CON_vesta.xyz")

        s_vesta = Structure(
            lattice=Lattice.from_parameters(22.6840, 13.3730, 12.5530, 90, 69.479, 90, vesta=True),
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
            [0, 0, 0],
            [0, 0, 1.089000],
            [1.026719, 0, -0.363],
            [-0.513360, -0.889165, -0.363],
            [-0.513360, 0.889165, -0.363],
        ]
        ch4 = ["C", "H", "H", "H", "H"]

        vectors = [[0, 0, 0], [4, 0, 0], [0, 4, 0], [4, 4, 0]]
        species = [atom for vec in vectors for atom in ch4]
        all_coords = [np.array(c) + vec for vec in vectors for c in coords]

        structure = Structure(Lattice.cubic(10), species, all_coords, coords_are_cartesian=True)

        for site in structure:
            if site.specie.symbol == "C":
                cluster = Molecule.from_sites(structure.extract_cluster([site]))
                assert cluster.formula == "H4 C1"

    def test_calculate_ase(self):
        pytest.importorskip("ase")
        struct_copy = self.cu_structure.copy()
        calculator = self.cu_structure.calculate(calculator=EMT(asap_cutoff=True))
        assert calculator.results["energies"] == approx([0.91304948, 0.91304948], abs=1e-5)
        assert calculator.results["free_energy"] == approx(1.8260989595, abs=1e-5)
        assert calculator.results["energy"] == approx(1.82609895)
        assert calculator.parameters == {"asap_cutoff": True}
        assert not hasattr(calculator, "dynamics")
        assert self.cu_structure == struct_copy, "original structure was modified"

    def test_relax_ase(self):
        pytest.importorskip("ase")
        struct_copy = self.cu_structure.copy()
        relaxed = self.cu_structure.relax(calculator=EMT(), relax_cell=False, optimizer="BFGS")
        assert relaxed.lattice == self.cu_structure.lattice
        assert {*relaxed.calc.results} >= {"energy", "energies", "free_energy"}
        assert relaxed.calc.results["energy"] == approx(1.82559661)
        assert relaxed.volume == approx(self.cu_structure.volume)
        assert relaxed.calc.parameters == {"asap_cutoff": False}
        assert relaxed.dynamics["optimizer"] == "BFGS"
        assert self.cu_structure == struct_copy, "original structure was modified"

    def test_relax_ase_return_traj(self):
        pytest.importorskip("ase")
        structure = self.cu_structure
        relaxed, traj = structure.relax(calculator=EMT(), fmax=0.01, return_trajectory=True)
        assert relaxed.lattice != structure.lattice
        assert {*relaxed.calc.results} >= {"energy", "energies", "free_energy"}
        assert relaxed.calc.parameters == {"asap_cutoff": False}
        assert relaxed.dynamics["optimizer"] == "FIRE"
        assert len(traj) == 15
        assert traj[0] != traj[-1]

        assert os.path.isfile("opt.traj")

    def test_relax_ase_opt_kwargs(self):
        pytest.importorskip("ase")
        structure = self.cu_structure
        traj_file = f"{self.tmp_path}/testing.traj"
        relaxed, traj = structure.relax(
            calculator=EMT(),
            fmax=0.01,
            steps=2,
            return_trajectory=True,
            opt_kwargs={"trajectory": traj_file},
        )
        assert relaxed.lattice != structure.lattice
        assert {*relaxed.calc.results} >= {"energy", "energies", "free_energy"}
        assert relaxed.calc.parameters == {"asap_cutoff": False}
        assert relaxed.dynamics["optimizer"] == "FIRE"
        assert len(traj) == 3  # there is an off-by-one in how ASE counts steps
        assert traj[0] != traj[-1]
        assert os.path.isfile(traj_file)

    @pytest.mark.xfail(reason="TODO: #3958 wait for matgl resolve of torch dependency")
    def test_calculate_matgl(self):
        pytest.importorskip("matgl")
        calculator = self.get_structure("Si").calculate()
        assert {*calculator.results} >= {"stress", "energy", "free_energy", "forces"}
        # Check the errors of predicted energy, forces and stress to be within
        # 0.1 eV/atom, 0.2 eV/Å, and 2 GPa.
        # The reference values here are predicted by M3GNet-MP-2021.2.8-DIRECT-PES in matgl.
        assert calculator.results["energy"] / self.get_structure("Si").num_sites == approx(-5.4146976, abs=0.1)
        assert np.linalg.norm(calculator.results["forces"]) == approx(7.8123485e-06, abs=0.2)
        assert np.linalg.norm(calculator.results["stress"]) == approx(1.7861567, abs=2)

    def test_relax_matgl(self):
        matgl = pytest.importorskip("matgl")
        struct = self.get_structure("Si")
        relaxed = struct.relax()
        assert relaxed.lattice.a == approx(3.860516230545545, rel=0.01)  # allow 1% error
        assert isinstance(relaxed.calc, matgl.ext.ase.PESCalculator)
        for key, val in {"type": "optimization", "optimizer": "FIRE"}.items():
            actual = relaxed.dynamics[key]
            assert actual == val, f"expected {key} to be {val}, {actual=}"
        relaxed_m3gnet = struct.relax("m3gnet")
        assert relaxed_m3gnet.lattice.a != relaxed.lattice.a
        assert relaxed.lattice.a == approx(3.8534658090100815, rel=0.01)  # allow 1% error

    def test_relax_m3gnet_fixed_lattice(self):
        matgl = pytest.importorskip("matgl")
        struct = self.get_structure("Si")
        relaxed = struct.relax(relax_cell=False, optimizer="BFGS")
        assert relaxed.lattice == struct.lattice
        assert isinstance(relaxed.calc, matgl.ext.ase.PESCalculator)
        assert relaxed.dynamics["optimizer"] == "BFGS"

    def test_relax_m3gnet_with_traj(self):
        pytest.importorskip("matgl")
        struct = self.get_structure("Si")
        relaxed, _trajectory = struct.relax(return_trajectory=True)
        assert relaxed.lattice.a == approx(3.867626620642243, abs=0.039)
        # assert sorted(trajectory.__dict__) == expected_attrs
        # for key in expected_attrs:
        #     # check for 2 atoms in Structure, 1 relax step in all observed trajectory attributes
        #     assert len(getattr(trajectory, key)) == {"atoms": 2}.get(key, 1)

    def test_from_prototype(self):
        for prototype in ["bcc", "fcc", "hcp", "diamond"]:
            struct = Structure.from_prototype(prototype, ["C"], a=3, c=4)
            assert isinstance(struct, Structure)

        with pytest.raises(ValueError, match="Required parameter 'c' not specified as a kwargs"):
            Structure.from_prototype("hcp", ["C"], a=3)

        struct = Structure.from_prototype("rocksalt", ["Li", "Cl"], a=2.56)
        expected_struct_str = """Full Formula (Li4 Cl4)
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
        assert str(struct) == expected_struct_str
        for prototype in ("cscl", "fluorite", "antifluorite", "zincblende"):
            struct = Structure.from_prototype(prototype, ["Cs", "Cl"], a=5)
            assert struct.lattice.is_orthogonal

    def test_to_primitive(self):
        struct = Structure.from_file(f"{TEST_FILES_DIR}/cif/orci_1010.cif")
        primitive = struct.to_primitive()

        assert struct != primitive
        sga = SpacegroupAnalyzer(struct)
        assert primitive == sga.get_primitive_standard_structure()
        assert struct.formula == "Mn1 B4"
        assert primitive.formula == "Mn1 B4"

    def test_to_conventional(self):
        struct = Structure.from_file(f"{TEST_FILES_DIR}/cif/bcc_1927.cif")
        conventional = struct.to_conventional()

        assert struct != conventional
        sga = SpacegroupAnalyzer(struct)
        assert conventional == sga.get_conventional_standard_structure()
        assert struct.formula == "Dy8 Sb6"
        assert conventional.formula == "Dy16 Sb12"

    def test_to_from_ase_atoms(self):
        pytest.importorskip("ase")

        atoms = self.struct.to_ase_atoms()
        assert isinstance(atoms, Atoms)
        assert len(atoms) == len(self.struct)

        assert AseAtomsAdaptor.get_structure(atoms) == self.struct

        assert Structure.from_ase_atoms(atoms) == self.struct
        assert type(Structure.from_ase_atoms(atoms)) is Structure

        labeled_atoms = self.labeled_structure.to_ase_atoms()
        assert Structure.from_ase_atoms(labeled_atoms) == self.labeled_structure

        with pytest.raises(ValueError, match="ASE Atoms only supports ordered structures"):
            self.disordered.to_ase_atoms()

    def test_struct_with_isotope(self):
        struct = Structure.from_file(f"{VASP_IN_DIR}/POSCAR_LiFePO4")
        struct = struct.replace_species({"Li": "H"})

        struct_deuter = struct.copy()
        struct_deuter.replace_species({"H": "D"})

        assert "Deuterium" not in [el.long_name for el in struct.composition.elements]
        assert "Deuterium" in [el.long_name for el in struct_deuter.composition.elements]
        assert struct_deuter == struct

        # test to make sure no Deuteriums are written to POSCAR
        struct_deuter.to(f"{self.tmp_path}/POSCAR_deuter")
        struct = Structure.from_file(f"{self.tmp_path}/POSCAR_deuter")
        assert "Deuterium" not in [el.long_name for el in struct.composition.elements]


class TestIMolecule(MatSciTest):
    def setup_method(self):
        coords = [
            [0, 0, 0],
            [0, 0, 1.089],
            [1.026719, 0, -0.363],
            [-0.513360, -0.889165, -0.363],
            [-0.513360, 0.889165, -0.363],
        ]
        self.coords = coords
        self.mol = Molecule(["C", "H", "H", "H", "H"], coords)

    def test_set_item(self):
        mol = self.mol.copy()
        mol[0] = "Si"
        assert mol.formula == "Si1 H4"
        mol[0, 1] = "Ge"
        assert mol.formula == "Ge2 H3"
        mol[:2] = "Sn"
        assert mol.formula == "Sn2 H3"

        mol = self.mol.copy()
        mol["H"] = "F"
        assert mol.formula == "C1 F4"
        mol["C"] = "C0.25Si0.5"
        assert mol.formula == "Si0.5 C0.25 F4"
        mol["C"] = "C0.25Si0.5"
        assert mol.formula == "Si0.625 C0.0625 F4"

    def test_bad_molecule(self):
        coords = [
            [0, 0, 0],
            [0, 0, 1.089000],
            [1.026719, 0, -0.363000],
            [-0.513360, -0.889165, -0.363000],
            [-0.513360, 0.889165, -0.363000],
            [-0.513360, 0.889165, -0.36301],
        ]
        with pytest.raises(StructureError, match="sites that are less than 0.01 Angstrom"):
            Molecule(["C", "H", "H", "H", "H", "H"], coords, validate_proximity=True)

    def test_get_angle_dihedral(self):
        assert self.mol.get_angle(1, 0, 2) == approx(109.47122144618737)
        assert self.mol.get_angle(3, 1, 2) == approx(60.00001388659683)
        assert self.mol.get_dihedral(0, 1, 2, 3) == approx(-35.26438851071765)

        coords = [[0, 0, 0], [0, 0, 1], [0, 1, 1], [1, 1, 1]]
        self.mol2 = Molecule(["C", "O", "N", "S"], coords)
        assert self.mol2.get_dihedral(0, 1, 2, 3) == approx(-90)

    def test_get_covalent_bonds(self):
        assert len(self.mol.get_covalent_bonds()) == 4

    def test_properties(self):
        assert len(self.mol) == 5
        assert self.mol.is_ordered
        assert self.mol.formula == "H4 C1"

    def test_properties_dict(self):
        properties = {"test_property": "test"}
        self.mol.properties = properties
        assert self.mol.properties == properties

    def test_repr_str(self):
        expected = """Full Formula (H4 C1)
Reduced Formula: H4C
Charge = 0.0, Spin Mult = 1
Sites (5)
0 C     0.000000     0.000000     0.000000
1 H     0.000000     0.000000     1.089000
2 H     1.026719     0.000000    -0.363000
3 H    -0.513360    -0.889165    -0.363000
4 H    -0.513360     0.889165    -0.363000"""
        assert str(self.mol) == expected
        expected = """Molecule Summary
Site: C (0.0000, 0.0000, 0.0000)
Site: H (0.0000, 0.0000, 1.0890)
Site: H (1.0267, 0.0000, -0.3630)
Site: H (-0.5134, -0.8892, -0.3630)
Site: H (-0.5134, 0.8892, -0.3630)"""
        assert repr(self.mol) == expected

    def test_site_properties(self):
        propertied_mol = Molecule(
            ["C", "H", "H", "H", "H"],
            self.coords,
            site_properties={"magmom": [0.5, -0.5, 1, 2, 3]},
        )
        assert propertied_mol[0].magmom == approx(0.5)
        assert propertied_mol[1].magmom == approx(-0.5)

    def test_chemical_system(self):
        assert self.mol.chemical_system == "C-H"

    def test_chemical_system_set(self):
        assert self.mol.chemical_system_set == {"C", "H"}

    def test_get_boxed_structure(self):
        struct = self.mol.get_boxed_structure(9, 9, 9)
        # C atom should be in center of box.
        assert_allclose(struct[4].frac_coords, [0.50000001, 0.5, 0.5])
        assert_allclose(struct[1].frac_coords, [0.6140799, 0.5, 0.45966667])
        with pytest.raises(ValueError, match="Box is not big enough to contain Molecule"):
            self.mol.get_boxed_structure(1, 1, 1)
        s2 = self.mol.get_boxed_structure(5, 5, 5, (2, 3, 4))
        assert len(s2) == 24 * 5
        assert s2.lattice.abc == (10, 15, 20)

        # Test offset option
        s3 = self.mol.get_boxed_structure(9, 9, 9, offset=[0.5, 0.5, 0.5])
        assert_allclose(s3[4].coords, [5, 5, 5])
        # Test no_cross option
        with pytest.raises(ValueError, match="Molecule crosses boundary of box"):
            self.mol.get_boxed_structure(5, 5, 5, offset=[10, 10, 10], no_cross=True)

        # Test reorder option
        no_reorder = self.mol.get_boxed_structure(10, 10, 10, reorder=False)
        assert str(s3[0].specie) == "H"
        assert str(no_reorder[0].specie) == "C"
        assert_allclose(no_reorder[2].frac_coords, [0.60267191, 0.5, 0.4637])

    def test_get_distance(self):
        assert self.mol.get_distance(0, 1) == approx(1.089)

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
        assert_allclose(
            self.mol.distance_matrix,
            [
                [0.0, 1.089, 1.08899995636, 1.08900040717, 1.08900040717],
                [1.089, 0.0, 1.77832952654, 1.7783298026, 1.7783298026],
                [1.08899995636, 1.77832952654, 0.0, 1.77833003783, 1.77833003783],
                [1.08900040717, 1.7783298026, 1.77833003783, 0.0, 1.77833],
                [1.08900040717, 1.7783298026, 1.77833003783, 1.77833, 0.0],
            ],
        )

    def test_get_zmatrix(self):
        mol = IMolecule(["C", "H", "H", "H", "H"], self.coords)
        z_matrix = """C
            H 1 B1
            H 1 B2 2 A2
            H 1 B3 2 A3 3 D3
            H 1 B4 2 A4 4 D4

            B1=1.089000
            B2=1.089000
            A2=109.471221
            B3=1.089000
            A3=109.471213
            D3=120.000017
            B4=1.089000
            A4=109.471213
            D4=119.999966
        """
        self.assert_str_content_equal(mol.get_zmatrix(), z_matrix)

    def test_break_bond(self):
        mol1, mol2 = self.mol.break_bond(0, 1)
        assert mol1.formula == "H3 C1"
        assert mol2.formula == "H1"

    def test_prop(self):
        assert self.mol.charge == 0
        assert self.mol.spin_multiplicity == 1
        assert self.mol.nelectrons == 10
        assert_allclose(self.mol.center_of_mass, [0, 0, 0], atol=1e-7)
        with pytest.raises(
            ValueError,
            match="Charge of 1 and spin multiplicity of 1 is not possible for this molecule",
        ):
            Molecule(
                ["C", "H", "H", "H", "H"],
                self.coords,
                charge=1,
                spin_multiplicity=1,
            )
        mol = Molecule(["C", "H", "H", "H", "H"], self.coords, charge=1)
        assert mol.spin_multiplicity == 2
        assert mol.nelectrons == 9
        # https://github.com/materialsproject/pymatgen/issues/3265
        # replace species and ensure nelectrons is updated
        mol[0] = "N"
        assert mol.nelectrons == 10

        # Triplet O2
        mol = IMolecule(["O"] * 2, [[0, 0, 0], [0, 0, 1.2]], spin_multiplicity=3)
        assert mol.spin_multiplicity == 3

    def test_no_spin_check(self):
        coords = [
            [0, 0, 0],
            [0, 0, 1.089000],
            [1.026719, 0, -0.363000],
            [-0.513360, -0.889165, -0.363000],
        ]
        with pytest.raises(
            ValueError,
            match="Charge of 0 and spin multiplicity of 1 is not possible for this molecule",
        ):
            mol = IMolecule(["C", "H", "H", "H"], coords, charge=0, spin_multiplicity=1)
        mol = IMolecule(
            ["C", "H", "H", "H"],
            coords,
            charge=0,
            spin_multiplicity=1,
            charge_spin_check=False,
        )
        assert mol.spin_multiplicity == 1
        assert mol.charge == 0

    def test_equal(self):
        mol = IMolecule(["C", "H", "H", "H", "H"], self.coords, charge=1)
        assert mol != self.mol

    def test_get_centered_molecule(self):
        mol = IMolecule(["O"] * 2, [[0, 0, 0], [0, 0, 1.2]], spin_multiplicity=3)
        centered = mol.get_centered_molecule()
        assert_allclose(centered.center_of_mass, [0, 0, 0])

    def test_as_from_dict(self):
        dct = self.mol.as_dict()
        mol2 = IMolecule.from_dict(dct)
        assert isinstance(mol2, IMolecule)
        propertied_mol = Molecule(
            ["C", "H", "H", "H", "H"],
            self.coords,
            charge=1,
            site_properties={"magmom": [0.5, -0.5, 1, 2, 3]},
            properties={"test_properties": "test"},
        )
        dct = propertied_mol.as_dict()
        assert dct["sites"][0]["properties"]["magmom"] == approx(0.5)
        mol = Molecule.from_dict(dct)
        assert propertied_mol == mol
        assert mol[0].magmom == approx(0.5)
        assert mol.formula == "H4 C1"
        assert mol.charge == 1
        assert mol.properties == {"test_properties": "test"}

    def test_default_dict_attrs(self):
        dct = self.mol.as_dict()
        assert dct["charge"] == 0
        assert dct["spin_multiplicity"] == 1

    def test_to_from_file_str(self):
        self.mol.properties["test_prop"] = 42
        for fmt in ("xyz", "json", "g03", "yaml", "yml"):
            mol = self.mol.to(fmt=fmt)
            assert isinstance(mol, str)
            mol = IMolecule.from_str(mol, fmt=fmt)
            if not mol.properties:
                # only fmt="json", "yaml", "yml" preserve properties, for other formats
                # properties are lost and we restore manually to make tests pass
                # TODO (janosh) long-term solution is to make all formats preserve properties
                mol.properties = self.mol.properties
            assert mol == self.mol
            assert isinstance(mol, IMolecule)

        ch4_xyz_str = self.mol.to(filename=f"{self.tmp_path}/CH4_testing.xyz")
        with open("CH4_testing.xyz", encoding="utf-8") as xyz_file:
            assert xyz_file.read() == ch4_xyz_str
        ch4_mol = IMolecule.from_file(f"{self.tmp_path}/CH4_testing.xyz")
        ch4_mol.properties = self.mol.properties
        assert self.mol == ch4_mol
        ch4_yaml_str = self.mol.to(filename=f"{self.tmp_path}/CH4_testing.yaml")

        with open("CH4_testing.yaml", encoding="utf-8") as yaml_file:
            assert yaml_file.read() == ch4_yaml_str
        ch4_mol = Molecule.from_file(f"{self.tmp_path}/CH4_testing.yaml")
        ch4_mol.properties = self.mol.properties
        assert self.mol == ch4_mol

    def test_to_file_alias(self):
        out_path = f"{self.tmp_path}/mol.gjf"
        assert self.mol.to(out_path) == self.mol.to_file(out_path)
        assert os.path.isfile(out_path)

    def test_to_from_ase_atoms(self):
        pytest.importorskip("ase")

        atoms = self.mol.to_ase_atoms()
        assert isinstance(atoms, Atoms)

        assert type(IMolecule.from_ase_atoms(atoms)) is IMolecule


class TestMolecule(MatSciTest):
    def setup_method(self):
        coords = [
            [0, 0, 0],
            [0, 0, 1.089000],
            [1.026719, 0, -0.363000],
            [-0.513360, -0.889165, -0.363000],
            [-0.513360, 0.889165, -0.363000],
        ]
        self.mol = Molecule(["C", "H", "H", "H", "H"], coords)

    def test_mutable_sequence_methods(self):
        mol = self.mol
        mol[1] = ("F", [0.5, 0.5, 0.5])
        assert mol.formula == "H3 C1 F1"
        assert_allclose(mol[1].coords, [0.5, 0.5, 0.5])
        mol.reverse()
        assert mol[0].specie == Element("H")
        assert_allclose(mol[0].coords, [-0.513360, 0.889165, -0.363000])
        del mol[1]
        assert mol.formula == "H2 C1 F1"
        mol[3] = "N", [0, 0, 0], {"charge": 4}
        assert mol.formula == "H2 N1 F1"
        assert mol[3].charge == 4

    def test_insert_remove_append(self):
        mol = self.mol
        returned = mol.insert(1, "O", [0.5, 0.5, 0.5])
        assert returned is mol
        assert mol.formula == "H4 C1 O1"
        del mol[2]
        assert mol.formula == "H3 C1 O1"
        mol.set_charge_and_spin(0)
        assert mol.spin_multiplicity == 2
        returned = mol.append("N", [1, 1, 1])
        assert returned is mol
        assert mol.formula == "H3 C1 N1 O1"
        with pytest.raises(TypeError, match="unhashable type: 'Molecule'"):
            _ = {mol: 1}
        returned = mol.remove_sites([0, 1])
        assert returned is mol
        assert mol.formula == "H3 N1"

    def test_from_sites(self):
        mol = Molecule.from_sites(self.mol)
        assert mol == self.mol

        with pytest.raises(ValueError, match="You need at least 1 site to make a Molecule"):
            Molecule.from_sites([])

    def test_translate_sites(self):
        returned = self.mol.translate_sites([0, 1], translation := (0.5, 0.5, 0.5))
        assert returned is self.mol
        assert tuple(self.mol.cart_coords[0]) == translation

    def test_rotate_sites(self):
        returned = self.mol.rotate_sites(theta=np.radians(30))
        assert returned is self.mol
        assert_allclose(self.mol.cart_coords[2], [0.889164737, 0.513359500, -0.363000000])

    def test_replace_species(self):
        self.mol[0] = "Ge"
        assert self.mol.formula == "Ge1 H4"

        returned = self.mol.replace_species({Element("Ge"): {Element("Ge"): 0.5, Element("Si"): 0.5}})
        assert returned is self.mol
        assert self.mol.formula == "Si0.5 Ge0.5 H4"

        # this should change the .5Si .5Ge sites to .75Si .25Ge
        self.mol.replace_species({Element("Ge"): {Element("Ge"): 0.5, Element("Si"): 0.5}})
        assert self.mol.formula == "Si0.75 Ge0.25 H4"

        dist = 0.1
        pre_perturbation_sites = self.mol.sites[:]
        self.mol.perturb(distance=dist)
        post_perturbation_sites = self.mol.sites

        for idx, site in enumerate(pre_perturbation_sites):
            assert site.distance(post_perturbation_sites[idx]) == approx(dist), "Bad perturbation distance"

    def test_add_remove_site_property(self):
        returned = self.mol.add_site_property("charge", [4.1, -2, -2, -2, -2])
        assert returned is self.mol
        assert self.mol[0].charge == approx(4.1)
        assert self.mol[1].charge == approx(-2)

        self.mol.add_site_property("magmom", [3, 2, 2, 2, 2])
        assert self.mol[0].charge == approx(4.1)
        assert self.mol[0].magmom == 3
        returned = self.mol.remove_site_property("magmom")
        assert returned is self.mol

        # test ValueError when values have wrong length
        with pytest.raises(ValueError, match=r"len\(values\)=2 must equal sites in structure=5"):
            self.mol.add_site_property("charge", [4, 2])

        with pytest.raises(AttributeError, match="attr='magmom' not found on Site"):
            _ = self.mol[0].magmom

    def test_as_from_dict(self):
        self.mol.append("X", [2, 0, 0])
        dct = self.mol.as_dict()
        mol2 = Molecule.from_dict(dct)
        assert isinstance(mol2, Molecule)
        self.assert_msonable(self.mol)

    def test_apply_operation(self):
        op = SymmOp.from_axis_angle_and_translation([0, 0, 1], 90)
        returned = self.mol.apply_operation(op)
        assert returned is self.mol
        assert_allclose(self.mol[2].coords, [0, 1.026719, -0.363000], atol=1e-12)

    def test_substitute(self):
        coords = [
            [0, 0, 1.08],
            [0, 0, 0],
            [1.026719, 0, -0.363000],
            [-0.513360, -0.889165, -0.363000],
            [-0.513360, 0.889165, -0.363000],
        ]
        sub = Molecule(["X", "C", "H", "H", "H"], coords)
        returned = self.mol.substitute(1, sub)
        assert returned is self.mol
        assert self.mol.get_distance(0, 4) == approx(1.54)
        self.mol.substitute(2, Molecule(["X", "F"], [[0, 0, 0], [0, 0, 1.11]]))
        assert self.mol.get_distance(0, 7) == approx(1.35)
        oh = Molecule(
            ["X", "O", "H"],
            [[0, 0.780362, -0.456316], [0, 0, 0.114079], [0, -0.780362, -0.456316]],
        )
        self.mol.substitute(1, oh)
        assert self.mol.get_distance(0, 7) == approx(1.43)
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
        assert benzene[11].coords[2] == approx(0)
        benzene[14] = "Br"
        benzene.substitute(13, sub)
        assert benzene.formula == "H9 C8 Br1"

    def test_to_from_file_str(self):
        for fmt in ["xyz", "json", "g03"]:
            mol = self.mol.to(fmt=fmt)
            assert mol is not None
            mol = Molecule.from_str(mol, fmt=fmt)
            assert mol == self.mol
            assert isinstance(mol, Molecule)

        self.mol.to(filename=f"{self.tmp_path}/CH4_testing.xyz")
        assert os.path.isfile(f"{self.tmp_path}/CH4_testing.xyz")

    def test_extract_cluster(self):
        species = self.mol.species * 2
        coords = [*self.mol.cart_coords, *(self.mol.cart_coords + np.array([10, 0, 0]))]
        mol = Molecule(species, coords)
        cluster = Molecule.from_sites(mol.extract_cluster([mol[0]]))
        assert mol.formula == "H8 C2"
        assert cluster.formula == "H4 C1"

    def test_no_spin_check(self):
        coords = [
            [0, 0, 0],
            [0, 0, 1.089000],
            [1.026719, 0, -0.363000],
            [-0.513360, -0.889165, -0.363000],
        ]
        expected_msg = "Charge of 0 and spin multiplicity of 1 is not possible for this molecule"
        with pytest.raises(ValueError, match=expected_msg):
            Molecule(["C", "H", "H", "H"], coords, charge=0, spin_multiplicity=1)
        mol_valid = Molecule(["C", "H", "H", "H"], coords, charge=0, spin_multiplicity=2)
        with pytest.raises(ValueError, match=expected_msg):
            mol_valid.set_charge_and_spin(0, 1)

    def test_set_charge_and_spin(self):
        mol = Molecule.from_dict(self.mol.as_dict() | dict(charge=0, spin_multiplicity=1, charge_spin_check=False))
        assert mol.spin_multiplicity == 1
        assert mol.charge == 0
        returned = mol.set_charge_and_spin(0, 3)
        assert returned is mol
        assert mol.charge == 0
        assert mol.spin_multiplicity == 3

    def test_calculate_ase_mol(self):
        pytest.importorskip("ase")
        mol_copy = self.mol.copy()
        calculator = self.mol.calculate(calculator=EMT(asap_cutoff=True))
        assert {*calculator.results} >= {"energy", "energies", "free_energy"}
        assert calculator.results["energy"] == approx(1.99570042)
        assert calculator.parameters == {"asap_cutoff": True}
        assert not hasattr(calculator, "dynamics")
        assert mol_copy == self.mol, "Molecule should not have been modified by calculation"

    def test_relax_ase_mol(self):
        pytest.importorskip("ase")
        mol = self.mol
        relaxed, traj = mol.relax(calculator=EMT(), fmax=0.01, optimizer="BFGS", return_trajectory=True)
        assert {*relaxed.calc.results} >= {"energy", "energies", "free_energy"}
        assert relaxed.calc.parameters == {"asap_cutoff": False}
        assert relaxed.dynamics["optimizer"] == "BFGS"
        assert len(traj) == 5
        assert traj[0] != traj[-1]

        assert os.path.isfile("opt.traj")

    def test_relax_ase_mol_return_traj(self):
        pytest.importorskip("ase")
        traj_file = f"{self.tmp_path}/testing.traj"
        relaxed, traj = self.mol.relax(
            calculator=EMT(),
            fmax=0.01,
            steps=2,
            return_trajectory=True,
            opt_kwargs={"trajectory": traj_file},
        )
        assert {*relaxed.calc.results} >= {"energy", "energies", "free_energy"}
        assert relaxed.calc.parameters == {"asap_cutoff": False}
        assert relaxed.dynamics["optimizer"] == "FIRE"
        assert len(traj) == 3  # there is an off-by-one in how ASE counts steps
        assert traj[0] != traj[-1]
        assert os.path.isfile(traj_file)

    def test_calculate_gfnxtb(self):
        pytest.importorskip("tblite")
        mol_copy = self.mol.copy()
        calculator = self.mol.calculate()
        assert isinstance(calculator, Calculator)
        assert not hasattr(calculator, "dynamics")
        assert calculator.results["energy"] == approx(-113.61022434200855)
        assert mol_copy == self.mol, "Molecule should not have been modified by calculation"

    @pytest.mark.xfail(reason="Pytorch and TBLite clash. https://github.com/materialsproject/pymatgen/pull/3060")
    def test_relax_gfnxtb(self):
        pytest.importorskip("tblite")
        mol = self.mol
        relaxed = mol.relax()
        assert hasattr(relaxed, "calc")
        assert hasattr(relaxed, "dynamics")
        assert relaxed.calc.results.get("energy")
        assert relaxed.dynamics == {"type": "optimization", "optimizer": "FIRE"}
        assert relaxed.calc.results["energy"] == approx(-113.61346199239306)

    def test_to_from_ase_atoms(self):
        pytest.importorskip("ase")

        atoms = self.mol.to_ase_atoms()
        assert isinstance(atoms, Atoms)
        assert len(atoms) == len(self.mol)

        assert AseAtomsAdaptor.get_molecule(atoms) == self.mol

        assert type(Molecule.from_ase_atoms(atoms)) is Molecule

        atoms = molecule("CH3")
        with pytest.raises(ValueError, match="not possible for this molecule"):
            Molecule.from_ase_atoms(atoms, spin_multiplicity=1)

        # Make sure kwargs are correctly passed
        atoms = molecule("CH3")
        assert type(Molecule.from_ase_atoms(atoms, charge_spin_check=False)) is Molecule

    def test_perturb(self):
        mol = self.mol
        mol_orig = mol.copy()
        mol.perturb(0.1)
        # Ensure all sites were perturbed by a distance of 0.1 Angstroms
        for site, site_orig in zip(mol, mol_orig, strict=True):
            cart_dist = site.distance(site_orig)
            # allow 1e-6 to account for numerical precision
            assert cart_dist == approx(0.1), f"Distance {cart_dist} > 0.1"

        # Check that the perturbation does not result in the same translation
        vecs = [site.coords - site_orig.coords for site, site_orig in zip(mol, mol_orig, strict=True)]
        compare_vecs = [np.allclose(v1, v2) for v1, v2 in itertools.pairwise(vecs)]
        assert not all(compare_vecs)

        # Test that same seed gives same perturbation
        s1 = self.mol.copy()
        s2 = self.mol.copy()
        s1.perturb(0.1, seed=42)
        s2.perturb(0.1, seed=42)
        for site1, site2 in zip(s1, s2, strict=True):
            assert site1.distance(site2) < 1e-7  # should be exactly equal up to numerical precision

        # Test that different seeds give different perturbations
        s3 = self.mol.copy()
        s3.perturb(0.1, seed=100)
        any_different = False
        for site1, site3 in zip(s1, s3, strict=True):
            if site1.distance(site3) > 1e-7:
                any_different = True
                break
        assert any_different, "Different seeds should give different perturbations"

        # Test min_distance
        s4 = self.mol.copy()
        s4.perturb(0.1, min_distance=0.05, seed=42)
        any_different = False
        for site1, site4 in zip(s1, s4, strict=True):
            if site1.distance(site4) > 1e-7:
                any_different = True
                break
        assert any_different, "Using min_distance should give different perturbations"
