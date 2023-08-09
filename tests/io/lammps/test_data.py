from __future__ import annotations

import gzip
import json
import os
import random
import unittest

import numpy as np
import pandas as pd
import pytest
from monty.json import MontyDecoder, MontyEncoder
from numpy.testing import assert_array_almost_equal
from pytest import approx
from ruamel.yaml import YAML

from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Molecule, Structure
from pymatgen.io.lammps.data import CombinedData, ForceField, LammpsBox, LammpsData, Topology, lattice_2_lmpbox
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

test_dir = f"{TEST_FILES_DIR}/lammps"


class TestLammpsBox(PymatgenTest):
    @classmethod
    def setUpClass(cls):
        cls.peptide = LammpsBox(
            bounds=[
                [36.840194, 64.211560],
                [41.013691, 68.385058],
                [29.768095, 57.139462],
            ]
        )
        cls.quartz = LammpsBox(
            bounds=[[0, 4.913400], [0, 4.255129], [0, 5.405200]],
            tilt=[-2.456700, 0.0, 0.0],
        )

    def test_volume(self):
        o_bounds = np.array(self.peptide.bounds)
        ov = np.prod(o_bounds[:, 1] - o_bounds[:, 0])
        assert self.peptide.volume == ov
        assert self.quartz.volume == approx(113.007331)

    def test_get_string(self):
        peptide = self.peptide.get_string(5)
        peptide_5 = """36.84019 64.21156  xlo xhi
41.01369 68.38506  ylo yhi
29.76809 57.13946  zlo zhi"""
        assert peptide == peptide_5
        quartz = self.quartz.get_string(4)
        quartz_4 = """0.0000 4.9134  xlo xhi
0.0000 4.2551  ylo yhi
0.0000 5.4052  zlo zhi
-2.4567 0.0000 0.0000  xy xz yz"""
        assert quartz == quartz_4

    def test_get_box_shift(self):
        peptide = self.peptide
        assert peptide.get_box_shift([1, 0, 0])[0] == 64.211560 - 36.840194
        assert peptide.get_box_shift([0, 0, -1])[-1] == 29.768095 - 57.139462
        quartz = self.quartz
        assert_array_almost_equal(quartz.get_box_shift([0, 0, 1]), [0, 0, 5.4052], 4)
        assert_array_almost_equal(quartz.get_box_shift([0, 1, -1]), [-2.4567, 4.2551, -5.4052], 4)
        assert_array_almost_equal(quartz.get_box_shift([1, -1, 0]), [4.9134 + 2.4567, -4.2551, 0], 4)

    def test_to_lattice(self):
        peptide = self.peptide.to_lattice()
        assert_array_almost_equal(peptide.abc, [27.371367] * 3)
        assert peptide.is_orthogonal
        quartz = self.quartz.to_lattice()
        assert_array_almost_equal(
            quartz.matrix,
            [[4.913400, 0, 0], [-2.456700, 4.255129, 0], [0, 0, 5.405200]],
        )


class TestLammpsData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.peptide = LammpsData.from_file(filename=f"{test_dir}/data.peptide")
        cls.ethane = LammpsData.from_file(filename=f"{test_dir}/ethane.data")
        cls.quartz = LammpsData.from_file(filename=f"{test_dir}/data.quartz", atom_style="atomic")
        cls.virus = LammpsData.from_file(filename=f"{test_dir}/virus.data", atom_style="angle")
        cls.tatb = LammpsData.from_file(
            filename=f"{test_dir}/tatb.data",
            atom_style="charge",
            sort_id=True,
        )

    def test_structure(self):
        quartz = self.quartz.structure
        assert_array_almost_equal(
            quartz.lattice.matrix,
            [[4.913400, 0, 0], [-2.456700, 4.255129, 0], [0, 0, 5.405200]],
        )
        assert quartz.formula == "Si3 O6"
        assert "molecule-ID" not in self.quartz.atoms.columns

        ethane = self.ethane.structure
        assert_array_almost_equal(ethane.lattice.matrix, np.diag([10.0] * 3))
        l_bounds = np.array(self.ethane.box.bounds)[:, 0]
        coords = self.ethane.atoms[["x", "y", "z"]] - l_bounds
        assert_array_almost_equal(ethane.cart_coords, coords)
        assert_array_almost_equal(ethane.site_properties["charge"], self.ethane.atoms["q"])
        tatb = self.tatb.structure
        frac_coords = tatb.frac_coords[381]
        real_frac_coords = frac_coords - np.floor(frac_coords)
        assert_array_almost_equal(real_frac_coords, [0.01553397, 0.71487872, 0.14134139])

        co = Structure.from_spacegroup(194, Lattice.hexagonal(2.50078, 4.03333), ["Co"], [[1 / 3, 2 / 3, 1 / 4]])
        ld_co = LammpsData.from_structure(co)
        assert ld_co.structure.composition.reduced_formula == "Co"
        ni = Structure.from_spacegroup(225, Lattice.cubic(3.50804), ["Ni"], [[0, 0, 0]])
        ld_ni = LammpsData.from_structure(ni)
        assert ld_ni.structure.composition.reduced_formula == "Ni"

    def test_sort_structure(self):
        struct = Structure(Lattice.cubic(4), ["S", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        lmp = LammpsData.from_structure(struct, is_sort=False)
        lmp.write_file("test1.data")
        lmp2 = LammpsData.from_file("test1.data", atom_style="charge")

        # internally element:type will be {Fe: 1, S: 2},
        # therefore without sorting the atom types in structure
        # will be [2, 1], i.e., (S, Fe)
        assert lmp2.atoms["type"].tolist() == [2, 1]

        # with sorting the atom types in structures will be [1, 2]
        lmp = LammpsData.from_structure(struct, is_sort=True)
        lmp.write_file("test1.data")
        lmp2 = LammpsData.from_file("test1.data", atom_style="charge")
        assert lmp2.atoms["type"].tolist() == [1, 2]

    def test_get_string(self):
        pep = self.peptide.get_string(distance=7, velocity=5, charge=4)
        pep_lines = pep.split("\n")
        pep_kws = [
            "Masses",
            "Pair Coeffs",
            "Bond Coeffs",
            "Angle Coeffs",
            "Dihedral Coeffs",
            "Improper Coeffs",
            "Atoms",
            "Velocities",
            "Bonds",
            "Angles",
            "Dihedrals",
            "Impropers",
        ]
        kw_inds = {line: idx for idx, line in enumerate(pep_lines) if line in pep_kws}
        # section sequence
        assert sorted(kw_inds, key=kw_inds.get) == pep_kws
        # header
        pep_header = "\n".join(pep_lines[: kw_inds["Masses"]])
        pep_header_7 = """Generated by pymatgen.io.lammps.data.LammpsData

2004  atoms
1365  bonds
 786  angles
 207  dihedrals
  12  impropers

  14  atom types
  18  bond types
  31  angle types
  21  dihedral types
   2  improper types

36.8401940 64.2115600  xlo xhi
41.0136910 68.3850580  ylo yhi
29.7680950 57.1394620  zlo zhi
"""
        assert pep_header == pep_header_7
        # int vs float for coeffs
        pep_dihedral_coeff = pep_lines[kw_inds["Dihedral Coeffs"] + 2]
        assert pep_dihedral_coeff == "1   0.200  1  180  1.0"
        # distance and charge
        pep_atom = pep_lines[kw_inds["Atoms"] + 2]
        assert pep_atom == "1       1   1  0.5100 43.9999300 58.5267800 36.7855000  0  0  0"
        # velocity
        pep_velo = pep_lines[kw_inds["Velocities"] + 2]
        assert pep_velo == "1    -0.00067 -0.00282  0.00383"
        # no floats in topology sections
        pep_topos = "\n".join(pep_lines[kw_inds["Bonds"] :])
        assert "." not in pep_topos

        c2h6 = self.ethane.get_string(distance=5, charge=3)
        c2h6_lines = c2h6.split("\n")
        c2h6_kws = [
            "Masses",
            "Pair Coeffs",
            "Bond Coeffs",
            "Angle Coeffs",
            "Dihedral Coeffs",
            "Improper Coeffs",
            "BondBond Coeffs",
            "BondAngle Coeffs",
            "MiddleBondTorsion Coeffs",
            "EndBondTorsion Coeffs",
            "AngleTorsion Coeffs",
            "AngleAngleTorsion Coeffs",
            "BondBond13 Coeffs",
            "AngleAngle Coeffs",
            "Atoms",
            "Bonds",
            "Angles",
            "Dihedrals",
            "Impropers",
        ]
        kw_inds = {line: idx for idx, line in enumerate(c2h6_lines) if line in c2h6_kws}
        # section sequence
        assert sorted(kw_inds, key=kw_inds.get) == c2h6_kws
        # header
        c2h6_header = "\n".join(c2h6_lines[: kw_inds["Masses"]])
        c2h6_header_5 = """Generated by pymatgen.io.lammps.data.LammpsData

 8  atoms
 7  bonds
12  angles
 9  dihedrals
 8  impropers

 2  atom types
 2  bond types
 2  angle types
 1  dihedral types
 2  improper types

0.21455 10.21454  xlo xhi
0.11418 10.11418  ylo yhi
-10.00014 -0.00015  zlo zhi
"""
        assert c2h6_header == c2h6_header_5
        # distance and charge
        c2h6_atom = c2h6_lines[kw_inds["Atoms"] + 2]
        assert c2h6_atom == "1  1  1 -0.080 4.46291 5.14833 -5.00041  0  0  0"
        # no floats in topology sections
        c2h6_topos = "\n".join(c2h6_lines[kw_inds["Bonds"] :])
        assert "." not in c2h6_topos

        quartz = self.quartz.get_string(distance=4)
        quartz_lines = quartz.split("\n")
        quartz_kws = ["Masses", "Atoms"]
        kw_inds = {line: idx for idx, line in enumerate(quartz_lines) if line in quartz_kws}
        # header
        quartz_header = "\n".join(quartz_lines[: kw_inds["Masses"]])
        quartz_header_4 = """Generated by pymatgen.io.lammps.data.LammpsData

9  atoms

2  atom types

0.0000 4.9134  xlo xhi
0.0000 4.2551  ylo yhi
0.0000 5.4052  zlo zhi
-2.4567 0.0000 0.0000  xy xz yz
"""
        assert quartz_header == quartz_header_4
        # distance
        quartz_atom = quartz_lines[kw_inds["Atoms"] + 2]
        assert quartz_atom == "1  1  2.3088  0.0000  3.6035"

        virus = self.virus.get_string()
        virus_lines = virus.split("\n")
        pairij_coeff = virus_lines[virus_lines.index("PairIJ Coeffs") + 5]
        assert pairij_coeff.strip().split() == ["1", "4", "1", "1.000", "1.12250"]

    def test_write_file(self):
        filename1 = "test1.data"
        self.ethane.write_file(filename=filename1)
        c2h6 = LammpsData.from_file(filename1)
        pd.testing.assert_frame_equal(c2h6.masses, self.ethane.masses)
        pd.testing.assert_frame_equal(c2h6.atoms, self.ethane.atoms)
        ff_kw = random.sample(sorted(self.ethane.force_field), 1)[0]
        pd.testing.assert_frame_equal(c2h6.force_field[ff_kw], self.ethane.force_field[ff_kw], ff_kw)
        topo_kw = random.sample(sorted(self.ethane.topology), 1)[0]
        pd.testing.assert_frame_equal(c2h6.topology[topo_kw], self.ethane.topology[topo_kw], topo_kw)
        filename2 = "test2.data"
        self.virus.write_file(filename=filename2)
        v = LammpsData.from_file(filename2, atom_style="angle")
        pd.testing.assert_frame_equal(v.force_field["PairIJ Coeffs"], self.virus.force_field["PairIJ Coeffs"])

    def test_disassemble(self):
        # general tests
        c = LammpsData.from_file(f"{test_dir}/crambin.data")
        _, c_ff, topos = c.disassemble()
        mass_info = [
            ("N1", 14.0067),
            ("H1", 1.00797),
            ("C1", 12.01115),
            ("H2", 1.00797),
            ("C2", 12.01115),
            ("O1", 15.9994),
            ("C3", 12.01115),
            ("O2", 15.9994),
            ("H3", 1.00797),
            ("C4", 12.01115),
            ("N2", 14.0067),
            ("C5", 12.01115),
            ("S1", 32.064),
            ("C6", 12.01115),
            ("N3", 14.0067),
            ("C7", 12.01115),
            ("C8", 12.01115),
            ("C9", 12.01115),
            ("O3", 15.9994),
        ]
        assert c_ff.mass_info == mass_info
        np.testing.assert_array_equal(c_ff.nonbond_coeffs, c.force_field["Pair Coeffs"].values)
        base_kws = ["Bond", "Angle", "Dihedral", "Improper"]
        for kw in base_kws:
            ff_kw = kw + " Coeffs"
            i = random.randint(0, len(c_ff.topo_coeffs[ff_kw]) - 1)
            sample_coeff = c_ff.topo_coeffs[ff_kw][i]
            np.testing.assert_array_equal(sample_coeff["coeffs"], c.force_field[ff_kw].iloc[i].values, ff_kw)
        topo = topos[-1]
        atoms = c.atoms[c.atoms["molecule-ID"] == 46]
        assert_array_almost_equal(topo.sites.cart_coords, atoms[["x", "y", "z"]])
        np.testing.assert_array_equal(topo.charges, atoms["q"])
        atom_labels = [m[0] for m in mass_info]
        assert topo.sites.site_properties["ff_map"] == [atom_labels[i - 1] for i in atoms["type"]]
        shift = min(atoms.index)
        for kw in base_kws:
            ff_kw = kw + " Coeffs"
            ff_coeffs = c_ff.topo_coeffs[ff_kw]
            topo_kw = kw + "s"
            topos_df = c.topology[topo_kw]
            topo_df = topos_df[topos_df["atom1"] >= shift]
            topo_arr = topo_df.drop("type", axis=1)
            np.testing.assert_array_equal(topo.topologies[topo_kw], topo_arr - shift, topo_kw)
            sample_topo = random.sample(list(topo_df.itertuples(False, None)), 1)[0]
            topo_type_idx = sample_topo[0] - 1
            topo_type = tuple(atom_labels[i - 1] for i in atoms.loc[list(sample_topo[1:])]["type"])

            assert topo_type in ff_coeffs[topo_type_idx]["types"], ff_kw
        # test no guessing element and pairij as non-bond coeffs
        v = self.virus
        _, v_ff, _ = v.disassemble(guess_element=False)
        assert v_ff.maps["Atoms"] == {"Qa1": 1, "Qb1": 2, "Qc1": 3, "Qa2": 4}
        pair_ij_coeffs = v.force_field["PairIJ Coeffs"].drop(["id1", "id2"], axis=1)
        np.testing.assert_array_equal(v_ff.nonbond_coeffs, pair_ij_coeffs.values)
        # test class2 ff
        _, e_ff, _ = self.ethane.disassemble()
        e_topo_coeffs = e_ff.topo_coeffs
        for k in ["BondBond Coeffs", "BondAngle Coeffs"]:
            assert k in e_topo_coeffs["Angle Coeffs"][0], k
        for k in [
            "MiddleBondTorsion Coeffs",
            "EndBondTorsion Coeffs",
            "AngleTorsion Coeffs",
            "AngleAngleTorsion Coeffs",
            "BondBond13 Coeffs",
        ]:
            assert k in e_topo_coeffs["Dihedral Coeffs"][0], k
        assert "AngleAngle Coeffs" in e_topo_coeffs["Improper Coeffs"][0]

    def test_from_file(self):
        # general tests
        pep = self.peptide
        # header stats and Nos. of columns
        assert pep.masses.shape == (14, 1)
        assert pep.atoms.shape == (2004, 9)
        assert list(pep.atoms.columns) == ["molecule-ID", "type", "q", "x", "y", "z", "nx", "ny", "nz"]
        topo = pep.topology
        assert topo["Bonds"].shape == (1365, 3)
        assert topo["Angles"].shape == (786, 4)
        assert topo["Dihedrals"].shape == (207, 5)
        assert topo["Impropers"].shape == (12, 5)
        ff = pep.force_field
        assert ff["Pair Coeffs"].shape == (14, 4)
        assert ff["Bond Coeffs"].shape == (18, 2)
        assert ff["Angle Coeffs"].shape == (31, 4)
        assert ff["Dihedral Coeffs"].shape == (21, 4)
        assert ff["Improper Coeffs"].shape == (2, 2)
        # header box
        np.testing.assert_array_equal(
            pep.box.bounds,
            [[36.840194, 64.211560], [41.013691, 68.385058], [29.768095, 57.139462]],
        )
        # body
        assert pep.masses.loc[7, "mass"] == 12.0110
        assert ff["Pair Coeffs"].loc[9, "coeff3"] == 0.152100
        assert ff["Bond Coeffs"].loc[5, "coeff2"] == 1.430000
        assert ff["Angle Coeffs"].loc[21, "coeff2"] == 120.000000
        assert ff["Dihedral Coeffs"].loc[10, "coeff1"] == 0.040000
        assert ff["Improper Coeffs"].loc[2, "coeff1"] == 20.000000
        assert pep.atoms.loc[29, "molecule-ID"] == 1
        assert pep.atoms.loc[29, "type"] == 7
        assert pep.atoms.loc[29, "q"] == -0.020
        assert pep.atoms.loc[29, "x"] == approx(42.96709)
        assert pep.atoms.loc[1808, "molecule-ID"] == 576
        assert pep.atoms.loc[1808, "type"] == 14
        assert pep.atoms.loc[1808, "y"] == approx(58.64352)
        assert pep.atoms.loc[1808, "nx"] == -1
        assert pep.velocities.loc[527, "vz"] == approx(-0.010889)
        assert topo["Bonds"].loc[47, "type"] == 8
        assert topo["Bonds"].loc[47, "atom2"] == 54
        assert topo["Bonds"].loc[953, "atom1"] == 1384
        assert topo["Angles"].loc[105, "type"] == 19
        assert topo["Angles"].loc[105, "atom3"] == 51
        assert topo["Angles"].loc[376, "atom2"] == 772
        assert topo["Dihedrals"].loc[151, "type"] == 14
        assert topo["Dihedrals"].loc[151, "atom4"] == 51
        assert topo["Impropers"].loc[4, "atom4"] == 32
        # class 2 and comments
        ethane = self.ethane
        assert ethane.masses.shape == (2, 1)
        assert ethane.atoms.shape == (8, 9)
        class2 = ethane.force_field
        assert class2["Pair Coeffs"].shape == (2, 2)
        assert class2["Bond Coeffs"].shape == (2, 4)
        assert class2["Angle Coeffs"].shape == (2, 4)
        assert class2["Dihedral Coeffs"].shape == (1, 6)
        assert class2["Improper Coeffs"].shape == (2, 2)
        assert class2["BondBond Coeffs"].loc[2, "coeff3"] == 1.1010
        assert class2["BondAngle Coeffs"].loc[2, "coeff4"] == 1.1010
        assert class2["AngleAngle Coeffs"].loc[2, "coeff6"] == 107.6600
        assert class2["AngleAngle Coeffs"].loc[2, "coeff6"] == 107.6600
        assert class2["AngleAngleTorsion Coeffs"].loc[1, "coeff3"] == 110.7700
        assert class2["EndBondTorsion Coeffs"].loc[1, "coeff8"] == 1.1010
        assert class2["MiddleBondTorsion Coeffs"].loc[1, "coeff4"] == 1.5300
        assert class2["BondBond13 Coeffs"].loc[1, "coeff3"] == 1.1010
        assert class2["AngleTorsion Coeffs"].loc[1, "coeff8"] == 110.7700
        # tilt box and another atom_style
        quartz = self.quartz
        np.testing.assert_array_equal(quartz.box.tilt, [-2.456700, 0.0, 0.0])
        assert list(quartz.atoms.columns) == ["type", "x", "y", "z"]
        assert quartz.atoms.loc[7, "x"] == approx(0.299963)
        # PairIJ Coeffs section
        virus = self.virus
        pairij = virus.force_field["PairIJ Coeffs"]
        assert pairij.loc[7, "id1"] == 3
        assert pairij.loc[7, "id2"] == 3
        assert pairij.loc[7, "coeff2"] == 2.1
        # sort_id
        atom_id = random.randint(1, 384)
        assert self.tatb.atoms.loc[atom_id].name == atom_id

    def test_from_ff_and_topologies(self):
        mass = {}
        mass["H"] = 1.0079401
        mass["O"] = 15.999400
        nonbond_coeffs = [[0.00774378, 0.98], [0.1502629, 3.1169]]
        topo_coeffs = {
            "Bond Coeffs": [{"coeffs": [176.864, 0.9611], "types": [("H", "O")]}],
            "Angle Coeffs": [{"coeffs": [42.1845, 109.4712], "types": [("H", "O", "H")]}],
        }
        ff = ForceField(mass.items(), nonbond_coeffs, topo_coeffs)
        with gzip.open(f"{test_dir}/topologies_ice.json.gz") as f:
            topo_dicts = json.load(f)
        topologies = [Topology.from_dict(d) for d in topo_dicts]
        box = LammpsBox([[-0.75694412, 44.165558], [0.38127473, 47.066074], [0.17900842, 44.193867]])
        ice = LammpsData.from_ff_and_topologies(box=box, ff=ff, topologies=topologies)
        atoms = ice.atoms
        bonds = ice.topology["Bonds"]
        angles = ice.topology["Angles"]
        np.testing.assert_array_equal(atoms.index.values, np.arange(1, len(atoms) + 1))
        np.testing.assert_array_equal(bonds.index.values, np.arange(1, len(bonds) + 1))
        np.testing.assert_array_equal(angles.index.values, np.arange(1, len(angles) + 1))

        i = random.randint(0, len(topologies) - 1)
        sample = topologies[i]
        in_atoms = ice.atoms[ice.atoms["molecule-ID"] == i + 1]
        np.testing.assert_array_equal(in_atoms.index.values, np.arange(3 * i + 1, 3 * i + 4))
        np.testing.assert_array_equal(in_atoms["type"].values, [2, 1, 1])
        np.testing.assert_array_equal(in_atoms["q"].values, sample.charges)
        np.testing.assert_array_equal(in_atoms[["x", "y", "z"]].values, sample.sites.cart_coords)
        broken_topo_coeffs = {
            "Bond Coeffs": [{"coeffs": [176.864, 0.9611], "types": [("H", "O")]}],
            "Angle Coeffs": [{"coeffs": [42.1845, 109.4712], "types": [("H", "H", "H")]}],
        }
        broken_ff = ForceField(mass.items(), nonbond_coeffs, broken_topo_coeffs)
        ld_woangles = LammpsData.from_ff_and_topologies(box=box, ff=broken_ff, topologies=[sample])
        assert "Angles" not in ld_woangles.topology

    def test_from_structure(self):
        latt = Lattice.monoclinic(9.78746, 4.75058, 8.95892, 115.9693)
        structure = Structure.from_spacegroup(
            15,
            latt,
            ["Os", "O", "O"],
            [
                [0, 0.25583, 0.75],
                [0.11146, 0.46611, 0.91631],
                [0.11445, 0.04564, 0.69518],
            ],
        )
        velocities = np.random.randn(20, 3) * 0.1
        structure.add_site_property("velocities", velocities)
        ld = LammpsData.from_structure(structure=structure, ff_elements=["O", "Os", "Na"])
        i = random.randint(0, 19)
        a = latt.matrix[0]
        va = velocities[i].dot(a) / np.linalg.norm(a)
        assert va == approx(ld.velocities.loc[i + 1, "vx"])
        assert velocities[i, 1] == approx(ld.velocities.loc[i + 1, "vy"])
        assert_array_almost_equal(ld.masses["mass"], [22.989769, 190.23, 15.9994])
        np.testing.assert_array_equal(ld.atoms["type"], [2] * 4 + [3] * 16)

    def test_set_charge_atom(self):
        peptide = self.peptide
        charges = {1: 0.8803}
        peptide.set_charge_atom(charges)
        assert peptide.atoms.loc[1, "q"] == 0.8803
        assert peptide.atoms.loc[2, "q"] == -0.270

    def test_set_charge_atom_type(self):
        peptide = self.peptide
        charges = {1: 0.8803}
        peptide.set_charge_atom_type(charges)
        assert peptide.atoms.loc[1, "q"] == 0.8803
        assert peptide.atoms.loc[2, "q"] == -0.270
        peptide.set_charge_atom_type({4: 2.345})
        assert peptide.atoms.loc[4, "q"] == 2.345
        assert peptide.atoms.loc[5, "q"] == 2.345
        assert peptide.atoms.loc[2004, "q"] == 0.4170

    def test_json_dict(self):
        encoded = json.dumps(self.ethane.as_dict(), cls=MontyEncoder)
        c2h6 = json.loads(encoded, cls=MontyDecoder)
        c2h6.masses.index = c2h6.masses.index.map(int)
        c2h6.atoms.index = c2h6.atoms.index.map(int)
        pd.testing.assert_frame_equal(c2h6.masses, self.ethane.masses)
        pd.testing.assert_frame_equal(c2h6.atoms, self.ethane.atoms)
        ff = self.ethane.force_field
        key, target_df = random.sample(sorted(ff.items()), 1)[0]
        c2h6.force_field[key].index = c2h6.force_field[key].index.map(int)
        assert pd.testing.assert_frame_equal(c2h6.force_field[key], target_df, check_dtype=False) is None, key
        topo = self.ethane.topology
        key, target_df = random.sample(sorted(topo.items()), 1)[0]
        c2h6.topology[key].index = c2h6.topology[key].index.map(int)
        assert pd.testing.assert_frame_equal(c2h6.topology[key], target_df) is None, key

    @classmethod
    def tearDownClass(cls):
        tmpfiles = ["test1.data", "test2.data"]
        for t in tmpfiles:
            if os.path.exists(t):
                os.remove(t)


class TestTopology(unittest.TestCase):
    def test_init(self):
        inner_charge = np.random.rand(10) - 0.5
        outer_charge = np.random.rand(10) - 0.5
        inner_velo = np.random.rand(10, 3) - 0.5
        outer_velo = np.random.rand(10, 3) - 0.5
        m = Molecule(
            ["H"] * 10,
            np.random.rand(10, 3) * 100,
            site_properties={
                "ff_map": ["D"] * 10,
                "charge": inner_charge,
                "velocities": inner_velo,
            },
        )
        # q and v from site properties, while type from species_string
        topo = Topology(sites=m)
        assert topo.type_by_sites == ["H"] * 10
        np.testing.assert_array_equal(topo.charges, inner_charge)
        np.testing.assert_array_equal(topo.velocities, inner_velo)
        # q and v from overriding, while type from site property
        topo_override = Topology(sites=m, ff_label="ff_map", charges=outer_charge, velocities=outer_velo)
        assert topo_override.type_by_sites == ["D"] * 10
        np.testing.assert_array_equal(topo_override.charges, outer_charge)
        np.testing.assert_array_equal(topo_override.velocities, outer_velo)
        # test using a list of sites instead of SiteCollection
        topo_from_list = Topology(sites=m.sites)
        assert topo_from_list.type_by_sites == topo.type_by_sites
        np.testing.assert_array_equal(topo_from_list.charges, topo.charges)
        np.testing.assert_array_equal(topo_from_list.velocities, topo.velocities)

    def test_from_bonding(self):
        # He: no bonding topologies
        helium = Molecule(["He"], [[0, 0, 0]])
        topo_he = Topology.from_bonding(molecule=helium)
        assert topo_he.topologies is None
        # H2: 1 bond only
        hydrogen = Molecule(["H"] * 2, [[0, 0, 0], [0, 0, 0.7414]])
        topo_h = Topology.from_bonding(molecule=hydrogen)
        tp_h = topo_h.topologies
        assert tp_h["Bonds"] == [[0, 1]]
        assert "Angles" not in tp_h
        assert "Dihedrals" not in tp_h
        # water: 2 bonds and 1 angle only
        water = Molecule(
            ["O", "H", "H"],
            [
                [0.0000, 0.0000, 0.1173],
                [0.0000, 0.7572, -0.4692],
                [0.0000, -0.7572, -0.4692],
            ],
        )
        topo_water = Topology.from_bonding(molecule=water)
        tp_water = topo_water.topologies
        assert tp_water["Bonds"] == [[0, 1], [0, 2]]
        assert tp_water["Angles"] == [[1, 0, 2]]
        assert "Dihedrals" not in tp_water
        # EtOH
        etoh = Molecule(
            ["C", "C", "O", "H", "H", "H", "H", "H", "H"],
            [
                [1.1879, -0.3829, 0.0000],
                [0.0000, 0.5526, 0.0000],
                [-1.1867, -0.2472, 0.0000],
                [-1.9237, 0.3850, 0.0000],
                [2.0985, 0.2306, 0.0000],
                [1.1184, -1.0093, 0.8869],
                [1.1184, -1.0093, -0.8869],
                [-0.0227, 1.1812, 0.8852],
                [-0.0227, 1.1812, -0.8852],
            ],
        )
        topo_etoh = Topology.from_bonding(molecule=etoh)
        tp_etoh = topo_etoh.topologies
        assert len(tp_etoh["Bonds"]) == 8
        etoh_bonds = [[0, 1], [0, 4], [0, 5], [0, 6], [1, 2], [1, 7], [1, 8], [2, 3]]
        np.testing.assert_array_equal(tp_etoh["Bonds"], etoh_bonds)
        assert len(tp_etoh["Angles"]) == 13
        etoh_angles = [
            [1, 0, 4],
            [1, 0, 5],
            [1, 0, 6],
            [4, 0, 5],
            [4, 0, 6],
            [5, 0, 6],
            [0, 1, 2],
            [0, 1, 7],
            [0, 1, 8],
            [2, 1, 7],
            [2, 1, 8],
            [7, 1, 8],
            [1, 2, 3],
        ]
        np.testing.assert_array_equal(tp_etoh["Angles"], etoh_angles)
        assert len(tp_etoh["Dihedrals"]) == 12
        etoh_dihedrals = [
            [4, 0, 1, 2],
            [4, 0, 1, 7],
            [4, 0, 1, 8],
            [5, 0, 1, 2],
            [5, 0, 1, 7],
            [5, 0, 1, 8],
            [6, 0, 1, 2],
            [6, 0, 1, 7],
            [6, 0, 1, 8],
            [0, 1, 2, 3],
            [7, 1, 2, 3],
            [8, 1, 2, 3],
        ]
        np.testing.assert_array_equal(tp_etoh["Dihedrals"], etoh_dihedrals)
        assert json.dumps(topo_etoh.as_dict()) is not None
        # bond flag to off
        topo_etoh0 = Topology.from_bonding(molecule=etoh, bond=False, angle=True, dihedral=True)
        assert topo_etoh0.topologies is None
        # angle or dihedral flag to off
        topo_etoh1 = Topology.from_bonding(molecule=etoh, angle=False)
        assert "Angles" not in topo_etoh1.topologies
        topo_etoh2 = Topology.from_bonding(molecule=etoh, dihedral=False)
        assert "Dihedrals" not in topo_etoh2.topologies


class TestForceField(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mass_info = [
            ("A", "H"),
            ("B", Element("C")),
            ("C", Element("O")),
            ("D", 1.00794),
        ]
        nonbond_coeffs = [
            [1, 1, 1.1225],
            [1, 1.175, 1.31894],
            [1, 1.55, 1.73988],
            [1, 1, 1.1225],
            [1, 1.35, 4],
            [1, 1.725, 1.93631],
            [1, 1.175, 1.31894],
            [1, 2.1, 4],
            [1, 1.55, 1.73988],
            [1, 1, 1.1225],
        ]
        topo_coeffs = {
            "Bond Coeffs": [
                {"coeffs": [50, 0.659469], "types": [("A", "B"), ("C", "D")]},
                {"coeffs": [50, 0.855906], "types": [("B", "C")]},
            ]
        }
        cls.virus = ForceField(mass_info=mass_info, nonbond_coeffs=nonbond_coeffs, topo_coeffs=topo_coeffs)
        cls.ethane = ForceField.from_file(f"{test_dir}/ff_ethane.yaml")

    def test_init(self):
        v = self.virus
        assert v.mass_info == [("A", 1.00794), ("B", 12.0107), ("C", 15.9994), ("D", 1.00794)]
        assert v.masses.loc[3, "mass"] == 15.9994
        v_ff = v.force_field
        assert isinstance(v_ff, dict)
        assert "Pair Coeffs" not in v_ff
        assert v_ff["PairIJ Coeffs"].iloc[5, 4] == 1.93631
        assert v_ff["Bond Coeffs"].loc[2, "coeff2"] == 0.855906
        v_maps = v.maps
        assert v_maps["Atoms"] == {"A": 1, "B": 2, "C": 3, "D": 4}
        assert v_maps["Bonds"] == {
            ("A", "B"): 1,
            ("C", "D"): 1,
            ("B", "A"): 1,
            ("D", "C"): 1,
            ("B", "C"): 2,
            ("C", "B"): 2,
        }
        e = self.ethane
        assert e.masses.loc[1, "mass"] == 12.01115
        e_ff = e.force_field
        assert isinstance(e_ff, dict)
        assert "PairIJ Coeffs" not in e_ff
        assert e_ff["Pair Coeffs"].loc[1, "coeff2"] == 3.854
        assert e_ff["Bond Coeffs"].loc[2, "coeff4"] == 844.6
        assert e_ff["Angle Coeffs"].loc[2, "coeff4"] == -2.4318
        assert e_ff["Dihedral Coeffs"].loc[1, "coeff1"] == -0.1432
        assert e_ff["Improper Coeffs"].loc[2, "coeff2"] == 0.0
        assert e_ff["BondBond Coeffs"].loc[2, "coeff1"] == 5.3316
        assert e_ff["BondAngle Coeffs"].loc[1, "coeff3"] == 1.53
        assert e_ff["MiddleBondTorsion Coeffs"].loc[1, "coeff1"] == -14.261
        assert e_ff["EndBondTorsion Coeffs"].loc[1, "coeff1"] == 0.213
        assert e_ff["AngleTorsion Coeffs"].loc[1, "coeff3"] == -0.2466
        assert e_ff["AngleAngleTorsion Coeffs"].loc[1, "coeff1"] == -12.564
        assert e_ff["BondBond13 Coeffs"].loc[1, "coeff1"] == 0.0
        assert e_ff["AngleAngle Coeffs"].loc[1, "coeff2"] == -0.4825
        e_maps = e.maps
        assert e_maps["Atoms"] == {"c4": 1, "h1": 2}
        assert e_maps["Bonds"] == {("c4", "c4"): 1, ("c4", "h1"): 2, ("h1", "c4"): 2}
        assert e_maps["Angles"] == {("c4", "c4", "h1"): 1, ("h1", "c4", "c4"): 1, ("h1", "c4", "h1"): 2}
        assert e_maps["Impropers"] == {
            ("c4", "c4", "h1", "h1"): 1,
            ("c4", "h1", "c4", "h1"): 1,
            ("h1", "h1", "c4", "c4"): 1,
            ("h1", "c4", "h1", "c4"): 1,
            ("h1", "c4", "h1", "h1"): 2,
            ("h1", "h1", "c4", "h1"): 2,
        }

    def test_to_file(self):
        filename = "ff_test.yaml"
        v = self.virus
        v.to_file(filename=filename)
        yaml = YAML()
        with open(filename) as f:
            d = yaml.load(f)
        # assert d["mass_info"] == [list(m) for m in v.mass_info]
        assert d["nonbond_coeffs"] == v.nonbond_coeffs

    def test_from_file(self):
        e = self.ethane
        assert e.mass_info == [("c4", 12.01115), ("h1", 1.00797)]
        np.testing.assert_array_equal(e.nonbond_coeffs, [[0.062, 3.854], [0.023, 2.878]])
        e_tc = e.topo_coeffs
        assert "Bond Coeffs" in e_tc
        assert "BondAngle Coeffs" in e_tc["Angle Coeffs"][0]
        assert "BondBond Coeffs" in e_tc["Angle Coeffs"][0]
        assert "AngleAngleTorsion Coeffs" in e_tc["Dihedral Coeffs"][0]
        assert "AngleTorsion Coeffs" in e_tc["Dihedral Coeffs"][0]
        assert "BondBond13 Coeffs" in e_tc["Dihedral Coeffs"][0]
        assert "EndBondTorsion Coeffs" in e_tc["Dihedral Coeffs"][0]
        assert "MiddleBondTorsion Coeffs" in e_tc["Dihedral Coeffs"][0]
        assert "AngleAngle Coeffs" in e_tc["Improper Coeffs"][0]

    def test_from_dict(self):
        d = self.ethane.as_dict()
        json_str = json.dumps(d)
        decoded = ForceField.from_dict(json.loads(json_str))
        assert decoded.mass_info == self.ethane.mass_info
        assert decoded.nonbond_coeffs == self.ethane.nonbond_coeffs
        assert decoded.topo_coeffs == self.ethane.topo_coeffs

    @classmethod
    def tearDownClass(cls):
        if os.path.exists("ff_test.yaml"):
            os.remove("ff_test.yaml")


class TestFunc(unittest.TestCase):
    def test_lattice_2_lmpbox(self):
        matrix = np.diag(np.random.randint(5, 14, size=(3,))) + np.random.rand(3, 3) * 0.2 - 0.1
        init_latt = Lattice(matrix)
        frac_coords = np.random.rand(10, 3)
        init_structure = Structure(init_latt, ["H"] * 10, frac_coords)
        origin = np.random.rand(3) * 10 - 5
        box, symmop = lattice_2_lmpbox(lattice=init_latt, origin=origin)
        boxed_latt = box.to_lattice()
        assert_array_almost_equal(init_latt.abc, boxed_latt.abc)
        assert_array_almost_equal(init_latt.angles, boxed_latt.angles)
        cart_coords = symmop.operate_multi(init_structure.cart_coords) - origin
        boxed_structure = Structure(boxed_latt, ["H"] * 10, cart_coords, coords_are_cartesian=True)
        assert_array_almost_equal(boxed_structure.frac_coords, frac_coords)
        tetra_latt = Lattice.tetragonal(5, 5)
        tetra_box, _ = lattice_2_lmpbox(tetra_latt)
        assert tetra_box.tilt is None
        ortho_latt = Lattice.orthorhombic(5, 5, 5)
        ortho_box, _ = lattice_2_lmpbox(ortho_latt)
        assert ortho_box.tilt is None
        rot_tetra_latt = Lattice([[5, 0, 0], [0, 2, 2], [0, -2, 2]])
        _, rotop = lattice_2_lmpbox(rot_tetra_latt)
        assert_array_almost_equal(
            rotop.rotation_matrix,
            [
                [1, 0, 0],
                [0, 2**0.5 / 2, 2**0.5 / 2],
                [0, -(2**0.5) / 2, 2**0.5 / 2],
            ],
        )


class TestCombinedData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ec = LammpsData.from_file(filename=f"{test_dir}/ec.data.gz")
        cls.fec = LammpsData.from_file(filename=f"{test_dir}/fec.data.gz")
        cls.li = LammpsData.from_file(filename=f"{test_dir}/li.data")
        cls.li_minimal = LammpsData.from_file(filename=f"{test_dir}/li_minimal.data")
        cls.coord = CombinedData.parse_xyz(filename=f"{test_dir}/ec_fec.xyz.gz")
        cls.small_coord = CombinedData.parse_xyz(filename=f"{test_dir}/li_ec.xyz")
        cls.small_coord_2 = CombinedData.parse_xyz(filename=f"{test_dir}/li_ec_2.xyz")
        cls.small_coord_3 = CombinedData.parse_xyz(filename=f"{test_dir}/li_2.xyz")
        cls.ec_fec1 = CombinedData.from_files(
            f"{test_dir}/ec_fec.xyz.gz",
            [1200, 300],
            f"{test_dir}/ec.data.gz",
            f"{test_dir}/fec.data.gz",
        )
        cls.ec_fec2 = CombinedData.from_lammpsdata([cls.ec, cls.fec], ["EC", "FEC"], [1200, 300], cls.coord)
        cls.ec_fec_ld = cls.ec_fec1.as_lammpsdata()
        cls.double_coord = pd.concat([cls.coord, cls.coord], ignore_index=True)
        cls.double_coord.index += 1
        cls.ec_fec3 = CombinedData.from_lammpsdata([cls.ec_fec_ld], ["EC FEC"], [2], cls.double_coord)
        cls.li_ec = CombinedData.from_lammpsdata([cls.li, cls.ec], ["Li", "EC"], [1, 1], cls.small_coord)
        cls.ec_li = CombinedData.from_lammpsdata([cls.ec, cls.li], ["EC", "Li"], [1, 1], cls.small_coord_2)
        cls.li_2 = CombinedData.from_lammpsdata([cls.li], ["Li"], [2], cls.small_coord_3)
        cls.li_2_minimal = CombinedData.from_lammpsdata([cls.li_minimal], ["Li_minimal"], [2], cls.small_coord_3)

    def test_from_files(self):
        # general tests
        ec_fec = self.ec_fec1
        # header stats and Nos. of columns
        assert ec_fec.names == ["cluster1", "cluster2"]
        assert ec_fec.nums == [1200, 300]
        assert ec_fec.masses.shape == (12, 1)
        assert ec_fec.atoms.shape == (15000, 6)
        assert list(ec_fec.atoms.columns) == ["molecule-ID", "type", "q", "x", "y", "z"]
        topo = ec_fec.topology
        assert topo["Bonds"].shape == (15000, 3)
        assert topo["Angles"].shape == (25500, 4)
        assert topo["Dihedrals"].shape == (42000, 5)
        assert topo["Impropers"].shape == (1500, 5)
        ff = ec_fec.force_field
        assert ff["Pair Coeffs"].shape == (12, 2)
        assert ff["Bond Coeffs"].shape == (15, 2)
        assert ff["Angle Coeffs"].shape == (24, 2)
        assert ff["Dihedral Coeffs"].shape == (39, 6)
        assert ff["Improper Coeffs"].shape == (2, 3)
        # header box
        np.testing.assert_array_equal(
            ec_fec.box.bounds,
            [[-0.597365, 54.56835], [-0.597365, 54.56835], [-0.597365, 54.56835]],
        )
        # body
        assert ec_fec.masses.loc[7, "mass"] == 1.008
        assert ff["Pair Coeffs"].loc[9, "coeff2"] == 3.750
        assert ff["Bond Coeffs"].loc[5, "coeff2"] == 1.0900
        assert ff["Angle Coeffs"].loc[24, "coeff2"] == 108.46005
        assert np.isnan(ff["Dihedral Coeffs"].loc[30, "coeff6"])
        assert ff["Improper Coeffs"].loc[2, "coeff1"] == 10.5
        assert ec_fec.atoms.loc[29, "molecule-ID"] == 3
        assert ec_fec.atoms.loc[29, "type"] == 5
        assert ec_fec.atoms.loc[29, "q"] == 0.0755
        assert ec_fec.atoms.loc[29, "x"] == approx(14.442260)
        assert ec_fec.atoms.loc[14958, "molecule-ID"] == 1496
        assert ec_fec.atoms.loc[14958, "type"] == 11
        assert ec_fec.atoms.loc[14958, "y"] == approx(41.010962)
        assert topo["Bonds"].loc[47, "type"] == 5
        assert topo["Bonds"].loc[47, "atom2"] == 47
        assert topo["Bonds"].loc[953, "atom1"] == 951
        assert topo["Angles"].loc[105, "type"] == 2
        assert topo["Angles"].loc[105, "atom3"] == 63
        assert topo["Angles"].loc[14993, "atom2"] == 8815
        assert topo["Dihedrals"].loc[151, "type"] == 4
        assert topo["Dihedrals"].loc[151, "atom4"] == 55
        assert topo["Dihedrals"].loc[41991, "type"] == 30
        assert topo["Dihedrals"].loc[41991, "atom2"] == 14994
        assert topo["Impropers"].loc[4, "atom4"] == 34

    def test_from_lammpsdata(self):
        # general tests
        ec_fec = self.ec_fec2
        # header stats and Nos. of columns
        assert ec_fec.names == ["EC", "FEC"]
        assert ec_fec.nums == [1200, 300]
        assert ec_fec.masses.shape == (12, 1)
        assert ec_fec.atoms.shape == (15000, 6)
        assert list(ec_fec.atoms.columns) == ["molecule-ID", "type", "q", "x", "y", "z"]
        topo = ec_fec.topology
        assert topo["Bonds"].shape == (15000, 3)
        assert topo["Angles"].shape == (25500, 4)
        assert topo["Dihedrals"].shape == (42000, 5)
        assert topo["Impropers"].shape == (1500, 5)
        ff = ec_fec.force_field
        assert ff["Pair Coeffs"].shape == (12, 2)
        assert ff["Bond Coeffs"].shape == (15, 2)
        assert ff["Angle Coeffs"].shape == (24, 2)
        assert ff["Dihedral Coeffs"].shape == (39, 6)
        assert ff["Improper Coeffs"].shape == (2, 3)
        # header box
        np.testing.assert_array_equal(
            ec_fec.box.bounds,
            [[-0.597365, 54.56835], [-0.597365, 54.56835], [-0.597365, 54.56835]],
        )
        # body
        assert ec_fec.masses.loc[7, "mass"] == 1.008
        assert ff["Pair Coeffs"].loc[9, "coeff2"] == 3.750
        assert ff["Bond Coeffs"].loc[5, "coeff2"] == 1.0900
        assert ff["Angle Coeffs"].loc[24, "coeff2"] == 108.46005
        assert np.isnan(ff["Dihedral Coeffs"].loc[30, "coeff6"])
        assert ff["Improper Coeffs"].loc[2, "coeff1"] == 10.5
        assert ec_fec.atoms.loc[29, "molecule-ID"] == 3
        assert ec_fec.atoms.loc[29, "type"] == 5
        assert ec_fec.atoms.loc[29, "q"] == 0.0755
        assert ec_fec.atoms.loc[29, "x"] == approx(14.442260)
        assert ec_fec.atoms.loc[14958, "molecule-ID"] == 1496
        assert ec_fec.atoms.loc[14958, "type"] == 11
        assert ec_fec.atoms.loc[14958, "y"] == approx(41.010962)
        assert topo["Bonds"].loc[47, "type"] == 5
        assert topo["Bonds"].loc[47, "atom2"] == 47
        assert topo["Bonds"].loc[953, "atom1"] == 951
        assert topo["Angles"].loc[105, "type"] == 2
        assert topo["Angles"].loc[105, "atom3"] == 63
        assert topo["Angles"].loc[14993, "atom2"] == 8815
        assert topo["Dihedrals"].loc[151, "type"] == 4
        assert topo["Dihedrals"].loc[151, "atom4"] == 55
        assert topo["Dihedrals"].loc[41991, "type"] == 30
        assert topo["Dihedrals"].loc[41991, "atom2"] == 14994
        assert topo["Impropers"].loc[4, "atom4"] == 34

        # non-destructively use of input (ID number)
        fec = self.fec
        topo = fec.topology
        ff = fec.force_field
        assert ff["Pair Coeffs"].index[0] == 1
        assert ff["Bond Coeffs"].index[0] == 1
        assert ff["Angle Coeffs"].index[0] == 1
        assert ff["Dihedral Coeffs"].index[0], 1
        assert ff["Improper Coeffs"].index[0] == 1
        assert fec.atoms.index[0] == 1
        assert fec.atoms.loc[1, "molecule-ID"] == 1
        assert fec.atoms.loc[1, "type"] == 1
        assert topo["Bonds"].index[0] == 1
        assert topo["Bonds"].loc[1, "type"] == 1
        assert topo["Bonds"].loc[1, "atom1"] == 1
        assert topo["Bonds"].loc[1, "atom2"] == 2
        assert topo["Angles"].index[0] == 1
        assert topo["Angles"].loc[1, "atom1"] == 1
        assert topo["Angles"].loc[1, "atom2"] == 3
        assert topo["Angles"].loc[1, "atom3"] == 4
        assert topo["Dihedrals"].index[0] == 1
        assert topo["Dihedrals"].loc[1, "atom1"] == 1
        assert topo["Dihedrals"].loc[1, "atom2"] == 3
        assert topo["Dihedrals"].loc[1, "atom3"] == 4
        assert topo["Dihedrals"].loc[1, "atom4"] == 5
        assert topo["Impropers"].index[0] == 1
        assert topo["Impropers"].loc[1, "atom1"] == 5
        assert topo["Impropers"].loc[1, "atom2"] == 4
        assert topo["Impropers"].loc[1, "atom3"] == 3
        assert topo["Impropers"].loc[1, "atom4"] == 6

        # tests for data objects with different number of ff kw
        li_ec = self.li_ec
        ec_li = self.ec_li
        assert li_ec.force_field["Pair Coeffs"].loc[6, "coeff2"] == 2.42
        assert ec_li.force_field["Pair Coeffs"].loc[6, "coeff2"] == 2.87
        assert li_ec.force_field["Bond Coeffs"].loc[5, "coeff2"] == 1.09
        assert ec_li.force_field["Bond Coeffs"].loc[5, "coeff2"] == 1.09
        assert li_ec.force_field["Angle Coeffs"].loc[7, "coeff2"] == 107.80
        assert ec_li.force_field["Angle Coeffs"].loc[7, "coeff2"] == 107.80
        assert li_ec.force_field["Dihedral Coeffs"].loc[11, "coeff2"] == 0.156
        assert ec_li.force_field["Dihedral Coeffs"].loc[11, "coeff2"] == 0.156
        assert li_ec.force_field["Improper Coeffs"].loc[1, "coeff1"] == 10.5
        assert ec_li.force_field["Improper Coeffs"].loc[1, "coeff1"] == 10.5

        # tests for combining data with no topo info
        li_2 = self.li_2
        assert li_2.topology is None, "Empty topo info should be none"

        # tests for combining data with no topo and ff info
        li_2_minimal = self.li_2_minimal
        assert li_2_minimal.force_field is None, "Empty ff info should be none"
        assert li_2_minimal.topology is None, "Empty topo info should be none"

    def test_get_string(self):
        # general tests
        ec_fec_lines = self.ec_fec1.get_string().splitlines()
        ec_fec_double_lines = self.ec_fec3.get_string().splitlines()
        # header information
        assert ec_fec_lines[1] == "# 1200 cluster1 + 300 cluster2"
        assert ec_fec_double_lines[1] == "# 2(1500) EC_FEC"
        # data type consistency tests
        assert ec_fec_lines[98] == "1  harmonic 3.200000000 -1 2"
        assert ec_fec_lines[109] == "12  charmm 2.700000000 2 180 0.0"
        assert ec_fec_lines[113] == "16  multi/harmonic 0.382999522 -1.148998570 0.000000000 1.531998090 0.000000000"
        assert ec_fec_lines[141] == "1  10.5 -1  2"
        assert ec_fec_double_lines[98] == "1  harmonic 3.200000000 -1 2"
        assert ec_fec_double_lines[109] == "12  charmm 2.700000000 2 180 0.0"
        assert (
            ec_fec_double_lines[113]
            == "16  multi/harmonic 0.382999522 -1.148998570 0.000000000 1.531998090 0.000000000"
        )
        assert ec_fec_double_lines[30146] == "30000  3000  12 -0.2329  4.630985  7.328547 51.604678"
        assert ec_fec_double_lines[141] == "1  10.5 -1  2"
        assert len(ec_fec_lines) == 99159
        assert len(ec_fec_double_lines) == 198159

    def test_structure(self):
        li_ec_structure = self.li_ec.structure
        assert_array_almost_equal(
            li_ec_structure.lattice.matrix,
            [[38.698274, 0, 0], [0, 38.698274, 0], [0, 0, 38.698274]],
        )
        assert_array_almost_equal(
            li_ec_structure.lattice.angles,
            (90.0, 90.0, 90.0),
        )
        assert li_ec_structure.formula == "Li1 H4 C3 O3"
        lbounds = np.array(self.li_ec.box.bounds)[:, 0]
        coords = self.li_ec.atoms[["x", "y", "z"]] - lbounds
        assert_array_almost_equal(li_ec_structure.cart_coords, coords)
        assert_array_almost_equal(li_ec_structure.site_properties["charge"], self.li_ec.atoms["q"])
        frac_coords = li_ec_structure.frac_coords[0]
        real_frac_coords = frac_coords - np.floor(frac_coords)
        assert_array_almost_equal(real_frac_coords, [0.01292047, 0.01292047, 0.01292047])

    def test_from_ff_and_topologies(self):
        with pytest.raises(AttributeError, match="Unsupported constructor for CombinedData objects"):
            CombinedData.from_ff_and_topologies()

    def test_from_structure(self):
        with pytest.raises(AttributeError, match="Unsupported constructor for CombinedData objects"):
            CombinedData.from_structure()

    def test_disassemble(self):
        # general tests
        ld = self.li
        cd = self.li_2
        _, cd_ff, topos = cd.disassemble()[0]
        mass_info = [
            ("Li1", 6.94),
        ]
        assert cd_ff.mass_info == mass_info
        np.testing.assert_array_equal(cd_ff.nonbond_coeffs, cd.force_field["Pair Coeffs"].values)

        topo = topos[-1]
        atoms = ld.atoms[ld.atoms["molecule-ID"] == 1]
        assert_array_almost_equal(topo.sites.cart_coords, atoms[["x", "y", "z"]])
        np.testing.assert_array_equal(topo.charges, atoms["q"])
        atom_labels = [m[0] for m in mass_info]
        assert topo.sites.site_properties["ff_map"] == [atom_labels[i - 1] for i in atoms["type"]]

        # test no guessing element
        v = self.li_2
        _, v_ff, _ = v.disassemble(guess_element=False)[0]
        assert v_ff.maps["Atoms"] == {"Qa1": 1}

    def test_json_dict(self):
        encoded = json.dumps(self.li_ec.as_dict(), cls=MontyEncoder)
        lic3o3h4 = json.loads(encoded, cls=MontyDecoder)
        assert lic3o3h4.nums == self.li_ec.nums
        assert lic3o3h4.names == self.li_ec.names
        assert lic3o3h4.atom_style == self.li_ec.atom_style
        pd.testing.assert_frame_equal(lic3o3h4.masses, self.li_ec.masses)
        pd.testing.assert_frame_equal(lic3o3h4.atoms, self.li_ec.atoms)
        ff = self.li_ec.force_field
        key, target_df = random.sample(sorted(ff.items()), 1)[0]
        lic3o3h4.force_field[key].index = lic3o3h4.force_field[key].index.map(int)
        assert pd.testing.assert_frame_equal(lic3o3h4.force_field[key], target_df, check_dtype=False) is None, key
        topo = self.li_ec.topology
        key, target_df = random.sample(sorted(topo.items()), 1)[0]
        assert pd.testing.assert_frame_equal(lic3o3h4.topology[key], target_df) is None, key
        lic3o3h4.mols[1].masses.index = lic3o3h4.mols[1].masses.index.map(int)
        lic3o3h4.mols[1].atoms.index = lic3o3h4.mols[1].atoms.index.map(int)
        pd.testing.assert_frame_equal(lic3o3h4.mols[1].masses, self.li_ec.mols[1].masses)
        pd.testing.assert_frame_equal(lic3o3h4.mols[1].atoms, self.li_ec.mols[1].atoms)
        ff_1 = self.li_ec.mols[1].force_field
        key, target_df = random.sample(sorted(ff_1.items()), 1)[0]
        lic3o3h4.mols[1].force_field[key].index = lic3o3h4.mols[1].force_field[key].index.map(int)
        assert (
            pd.testing.assert_frame_equal(lic3o3h4.mols[1].force_field[key], target_df, check_dtype=False) is None
        ), key
        topo_1 = self.li_ec.mols[1].topology
        key, target_df = random.sample(sorted(topo_1.items()), 1)[0]
        lic3o3h4.mols[1].topology[key].index = lic3o3h4.mols[1].topology[key].index.map(int)
        assert pd.testing.assert_frame_equal(lic3o3h4.mols[1].topology[key], target_df) is None, key

    def test_as_lammpsdata(self):
        ec_fec = self.ec_fec_ld
        assert ec_fec.masses.shape == (12, 1)
        assert ec_fec.atoms.shape == (15000, 6)
        assert list(ec_fec.atoms.columns) == ["molecule-ID", "type", "q", "x", "y", "z"]
        topo = ec_fec.topology
        assert topo["Bonds"].shape == (15000, 3)
        assert topo["Angles"].shape == (25500, 4)
        assert topo["Dihedrals"].shape == (42000, 5)
        assert topo["Impropers"].shape == (1500, 5)
        ff = ec_fec.force_field
        assert ff["Pair Coeffs"].shape == (12, 2)
        assert ff["Bond Coeffs"].shape == (15, 2)
        assert ff["Angle Coeffs"].shape == (24, 2)
        assert ff["Dihedral Coeffs"].shape == (39, 6)
        assert ff["Improper Coeffs"].shape == (2, 3)
        # header box
        np.testing.assert_array_equal(
            ec_fec.box.bounds,
            [[-0.597365, 54.56835], [-0.597365, 54.56835], [-0.597365, 54.56835]],
        )
        # body
        assert ec_fec.masses.loc[7, "mass"] == 1.008
        assert ff["Pair Coeffs"].loc[9, "coeff2"] == 3.750
        assert ff["Bond Coeffs"].loc[5, "coeff2"] == 1.0900
        assert ff["Angle Coeffs"].loc[24, "coeff2"] == 108.46005
        assert np.isnan(ff["Dihedral Coeffs"].loc[30, "coeff6"])
        assert ff["Improper Coeffs"].loc[2, "coeff1"] == 10.5
        assert ec_fec.atoms.loc[29, "molecule-ID"] == 3
        assert ec_fec.atoms.loc[29, "type"] == 5
        assert ec_fec.atoms.loc[29, "q"] == 0.0755
        assert ec_fec.atoms.loc[29, "x"] == approx(14.442260)
        assert ec_fec.atoms.loc[14958, "molecule-ID"] == 1496
        assert ec_fec.atoms.loc[14958, "type"] == 11
        assert ec_fec.atoms.loc[14958, "y"] == approx(41.010962)
        assert topo["Bonds"].loc[47, "type"] == 5
        assert topo["Bonds"].loc[47, "atom2"] == 47
        assert topo["Bonds"].loc[953, "atom1"] == 951
        assert topo["Angles"].loc[105, "type"] == 2
        assert topo["Angles"].loc[105, "atom3"] == 63
        assert topo["Angles"].loc[14993, "atom2"] == 8815
        assert topo["Dihedrals"].loc[151, "type"] == 4
        assert topo["Dihedrals"].loc[151, "atom4"] == 55
        assert topo["Dihedrals"].loc[41991, "type"] == 30
        assert topo["Dihedrals"].loc[41991, "atom2"] == 14994
        assert topo["Impropers"].loc[4, "atom4"] == 34
        ec_fec_lines = self.ec_fec_ld.get_string().splitlines()
        # header information
        assert ec_fec_lines[1] == ""
        # data type consistency tests
        assert ec_fec_lines[97] == "1  harmonic 3.200000000 -1 2"
        assert ec_fec_lines[108] == "12  charmm 2.700000000 2 180 0.0"
        assert ec_fec_lines[112] == "16  multi/harmonic 0.382999522 -1.148998570 0.000000000 1.531998090 0.000000000"
        assert ec_fec_lines[140] == "1  10.5 -1  2"
        assert len(ec_fec_lines) == 99159
