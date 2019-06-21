# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
import unittest
import os
import random
from collections import OrderedDict
import gzip
import json

import numpy as np
import pandas as pd
from ruamel.yaml import YAML
from pymatgen import Molecule, Element, Lattice, Structure

from pymatgen.io.lammps.data import LammpsBox, LammpsData, Topology,\
    ForceField, lattice_2_lmpbox, structure_2_lmpdata

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files", "lammps")


class LammpsBoxTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.peptide = LammpsBox(bounds=[[36.840194, 64.211560],
                                        [41.013691, 68.385058],
                                        [29.768095, 57.139462]])
        cls.quartz = LammpsBox(bounds=[[0, 4.913400],
                                       [0, 4.255129],
                                       [0, 5.405200]],
                               tilt=[-2.456700, 0.0, 0.0])

    def test_volume(self):
        obounds = np.array(self.peptide.bounds)
        ov = np.prod(obounds[:, 1] - obounds[:, 0])
        self.assertEqual(self.peptide.volume, ov)
        self.assertAlmostEqual(self.quartz.volume, 113.00733165874873)

    def test_get_string(self):
        peptide = self.peptide.get_string(5)
        peptide_5 = """36.84019 64.21156  xlo xhi
41.01369 68.38506  ylo yhi
29.76809 57.13946  zlo zhi"""
        self.assertEqual(peptide, peptide_5)
        quartz = self.quartz.get_string(4)
        quartz_4 = """0.0000 4.9134  xlo xhi
0.0000 4.2551  ylo yhi
0.0000 5.4052  zlo zhi
-2.4567 0.0000 0.0000  xy xz yz"""
        self.assertEqual(quartz, quartz_4)

    def test_get_box_shift(self):
        peptide = self.peptide
        self.assertEqual(peptide.get_box_shift([1, 0, 0])[0],
                         64.211560 - 36.840194)
        self.assertEqual(peptide.get_box_shift([0, 0, -1])[-1],
                         29.768095 - 57.139462)
        quartz = self.quartz
        np.testing.assert_array_almost_equal(quartz.get_box_shift([0, 0, 1]),
                                      [0, 0, 5.4052], 4)
        np.testing.assert_array_almost_equal(quartz.get_box_shift([0, 1, -1]),
                                             [-2.4567, 4.2551, -5.4052], 4)
        np.testing.assert_array_almost_equal(quartz.get_box_shift([1, -1, 0]),
                                             [4.9134 + 2.4567, -4.2551, 0], 4)

    def test_to_lattice(self):
        peptide = self.peptide.to_lattice()
        np.testing.assert_array_almost_equal(peptide.abc, [27.371367] * 3)
        self.assertTrue(peptide.is_orthogonal)
        quartz = self.quartz.to_lattice()
        np.testing.assert_array_almost_equal(quartz.matrix,
                                             [[4.913400, 0, 0],
                                              [-2.456700, 4.255129, 0],
                                              [0, 0, 5.405200]])


class LammpsDataTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.peptide = LammpsData.\
            from_file(filename=os.path.join(test_dir, "data.peptide"))
        cls.ethane = LammpsData.\
            from_file(filename=os.path.join(test_dir, "ethane.data"))
        cls.quartz = LammpsData.\
            from_file(filename=os.path.join(test_dir, "data.quartz"),
                      atom_style="atomic")
        cls.virus = LammpsData.\
            from_file(filename=os.path.join(test_dir, "virus.data"),
                      atom_style="angle")
        cls.tatb = LammpsData.\
            from_file(filename=os.path.join(test_dir, "tatb.data"),
                      atom_style="charge", sort_id=True)

    def test_structure(self):
        quartz = self.quartz.structure
        np.testing.assert_array_equal(quartz.lattice.matrix,
                                      [[4.913400, 0, 0],
                                       [-2.456700, 4.255129, 0],
                                       [0, 0, 5.405200]])
        self.assertEqual(quartz.formula, "Si3 O6")
        self.assertNotIn("molecule-ID", self.quartz.atoms.columns)

        ethane = self.ethane.structure
        np.testing.assert_array_equal(ethane.lattice.matrix,
                                      np.diag([10.0] * 3))
        lbounds = np.array(self.ethane.box.bounds)[:, 0]
        coords = self.ethane.atoms[["x", "y", "z"]].values - lbounds
        np.testing.assert_array_equal(ethane.cart_coords, coords)
        np.testing.assert_array_equal(ethane.site_properties["charge"],
                                      self.ethane.atoms["q"])
        tatb = self.tatb.structure
        frac_coords = tatb.frac_coords[381]
        real_frac_coords = frac_coords - np.floor(frac_coords)
        np.testing.assert_array_almost_equal(real_frac_coords,
                                             [0.01553397,
                                              0.71487872,
                                              0.14134139])

        co = Structure.from_spacegroup(194,
                                       Lattice.hexagonal(2.50078, 4.03333),
                                       ["Co"], [[1/3, 2/3, 1/4]])
        ld_co = LammpsData.from_structure(co)
        self.assertEqual(ld_co.structure.composition.reduced_formula, "Co")
        ni = Structure.from_spacegroup(225, Lattice.cubic(3.50804),
                                       ["Ni"], [[0, 0, 0]])
        ld_ni = LammpsData.from_structure(ni)
        self.assertEqual(ld_ni.structure.composition.reduced_formula, "Ni")

    def test_get_string(self):
        pep = self.peptide.get_string(distance=7, velocity=5, charge=4)
        pep_lines = pep.split("\n")
        pep_kws = ["Masses", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs",
                   "Dihedral Coeffs", "Improper Coeffs", "Atoms",
                   "Velocities", "Bonds", "Angles", "Dihedrals", "Impropers"]
        kw_inds = {l: i for i, l in enumerate(pep_lines) if l in pep_kws}
        # section sequence
        self.assertListEqual([k for k in sorted(kw_inds, key=kw_inds.get)],
                             pep_kws)
        # header
        pep_header = "\n".join(pep_lines[:kw_inds["Masses"]])
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
        self.assertEqual(pep_header, pep_header_7)
        # int vs float for coeffs
        pep_dihedral_coeff = pep_lines[kw_inds["Dihedral Coeffs"] + 2]
        self.assertEqual(pep_dihedral_coeff, "1   0.200  1  180  1.0")
        # distance and charge
        pep_atom = pep_lines[kw_inds["Atoms"] + 2]
        self.assertEqual(pep_atom, "1       1   1  0.5100 43.9999300 "
                                   "58.5267800 36.7855000  0  0  0")
        # velocity
        pep_velo = pep_lines[kw_inds["Velocities"] + 2]
        self.assertEqual(pep_velo, "1    -0.00067 -0.00282  0.00383")
        # no floats in topology sections
        pep_topos = "\n".join(pep_lines[kw_inds["Bonds"]:])
        self.assertNotIn(".", pep_topos)

        c2h6 = self.ethane.get_string(distance=5, charge=3)
        c2h6_lines = c2h6.split("\n")
        c2h6_kws = ["Masses", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs",
                    "Dihedral Coeffs", "Improper Coeffs", "BondBond Coeffs",
                    "BondAngle Coeffs", "MiddleBondTorsion Coeffs",
                    "EndBondTorsion Coeffs", "AngleTorsion Coeffs",
                    "AngleAngleTorsion Coeffs", "BondBond13 Coeffs",
                    "AngleAngle Coeffs", "Atoms", "Bonds", "Angles",
                    "Dihedrals", "Impropers"]
        kw_inds = {l: i for i, l in enumerate(c2h6_lines) if l in c2h6_kws}
        # section sequence
        self.assertListEqual([k for k in sorted(kw_inds, key=kw_inds.get)],
                             c2h6_kws)
        # header
        c2h6_header = "\n".join(c2h6_lines[:kw_inds["Masses"]])
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
        self.assertEqual(c2h6_header, c2h6_header_5)
        # distance and charge
        c2h6_atom = c2h6_lines[kw_inds["Atoms"] + 2]
        self.assertEqual(c2h6_atom, "1  1  1 -0.080 4.46291 5.14833 -5.00041"
                                    "  0  0  0")
        # no floats in topology sections
        c2h6_topos = "\n".join(c2h6_lines[kw_inds["Bonds"]:])
        self.assertNotIn(".", c2h6_topos)

        quartz = self.quartz.get_string(distance=4)
        quartz_lines = quartz.split("\n")
        quartz_kws = ["Masses", "Atoms"]
        kw_inds = {l: i for i, l in enumerate(quartz_lines) if l in quartz_kws}
        # header
        quartz_header = "\n".join(quartz_lines[:kw_inds["Masses"]])
        quartz_header_4 = """Generated by pymatgen.io.lammps.data.LammpsData

9  atoms

2  atom types

0.0000 4.9134  xlo xhi
0.0000 4.2551  ylo yhi
0.0000 5.4052  zlo zhi
-2.4567 0.0000 0.0000  xy xz yz
"""
        self.assertEqual(quartz_header, quartz_header_4)
        # distance
        quartz_atom = quartz_lines[kw_inds["Atoms"] + 2]
        self.assertEqual(quartz_atom, "1  1  2.3088  0.0000  3.6035")

        virus = self.virus.get_string()
        virus_lines = virus.split("\n")
        pairij_coeff = virus_lines[virus_lines.index("PairIJ Coeffs") + 5]
        self.assertEqual(pairij_coeff.strip(), "1  4  1  1.000  1.12250")

    def test_write_file(self):
        filename1 = "test1.data"
        self.ethane.write_file(filename=filename1)
        c2h6 = LammpsData.from_file(filename1)
        pd.testing.assert_frame_equal(c2h6.masses, self.ethane.masses)
        pd.testing.assert_frame_equal(c2h6.atoms, self.ethane.atoms)
        ff_kw = random.sample(self.ethane.force_field.keys(), 1)[0]
        pd.testing.assert_frame_equal(c2h6.force_field[ff_kw],
                                      self.ethane.force_field[ff_kw], ff_kw)
        topo_kw = random.sample(self.ethane.topology.keys(), 1)[0]
        pd.testing.assert_frame_equal(c2h6.topology[topo_kw],
                                      self.ethane.topology[topo_kw], topo_kw)
        filename2 = "test2.data"
        self.virus.write_file(filename=filename2)
        v = LammpsData.from_file(filename2, atom_style="angle")
        pd.testing.assert_frame_equal(v.force_field["PairIJ Coeffs"],
                                      self.virus.force_field["PairIJ Coeffs"])

    def test_disassemble(self):
        # general tests
        c = LammpsData.from_file(os.path.join(test_dir, "crambin.data"))
        _, c_ff, topos = c.disassemble()
        mass_info = [('N1', 14.0067), ('H1', 1.00797), ('C1', 12.01115),
                     ('H2', 1.00797), ('C2', 12.01115), ('O1', 15.9994),
                     ('C3', 12.01115), ('O2', 15.9994), ('H3', 1.00797),
                     ('C4', 12.01115), ('N2', 14.0067), ('C5', 12.01115),
                     ('S1', 32.064), ('C6', 12.01115), ('N3', 14.0067),
                     ('C7', 12.01115), ('C8', 12.01115), ('C9', 12.01115),
                     ('O3', 15.9994)]
        self.assertListEqual(c_ff.mass_info, mass_info)
        np.testing.assert_array_equal(c_ff.nonbond_coeffs,
                                      c.force_field["Pair Coeffs"].values)
        base_kws = ["Bond", "Angle", "Dihedral", "Improper"]
        for kw in base_kws:
            ff_kw = kw + " Coeffs"
            i = random.randint(0, len(c_ff.topo_coeffs[ff_kw]) - 1)
            sample_coeff = c_ff.topo_coeffs[ff_kw][i]
            np.testing.\
                assert_array_equal(sample_coeff["coeffs"],
                                   c.force_field[ff_kw].iloc[i].values,
                                   ff_kw)
        topo = topos[-1]
        atoms = c.atoms[c.atoms["molecule-ID"] == 46]
        np.testing.assert_array_equal(topo.sites.cart_coords,
                                      atoms[["x", "y", "z"]])
        np.testing.assert_array_equal(topo.charges, atoms["q"])
        atom_labels = [m[0] for m in mass_info]
        self.assertListEqual(topo.sites.site_properties["ff_map"],
                             [atom_labels[i - 1] for i in atoms["type"]])
        shift = min(atoms.index)
        for kw in base_kws:
            ff_kw = kw + " Coeffs"
            ff_coeffs = c_ff.topo_coeffs[ff_kw]
            topo_kw = kw + "s"
            topos_df = c.topology[topo_kw]
            topo_df = topos_df[topos_df["atom1"] >= shift]
            topo_arr = topo_df.drop("type", axis=1).values
            np.testing.assert_array_equal(topo.topologies[topo_kw],
                                          topo_arr - shift, topo_kw)
            sample_topo = random.sample(list(topo_df.itertuples(False, None)),
                                        1)[0]
            topo_type_idx = sample_topo[0] - 1
            topo_type = tuple([atom_labels[i - 1] for i in
                               atoms.loc[sample_topo[1:], "type"]])

            self.assertIn(topo_type, ff_coeffs[topo_type_idx]["types"], ff_kw)
        # test no guessing element and pairij as nonbond coeffs
        v = self.virus
        _, v_ff, _ = v.disassemble(guess_element=False)
        self.assertDictEqual(v_ff.maps["Atoms"],
                             dict(Qa1=1, Qb1=2, Qc1=3, Qa2=4))
        pairij_coeffs = v.force_field["PairIJ Coeffs"].drop(["id1", "id2"],
                                                            axis=1)
        np.testing.assert_array_equal(v_ff.nonbond_coeffs,
                                      pairij_coeffs.values)
        # test class2 ff
        _, e_ff, _ = self.ethane.disassemble()
        e_topo_coeffs = e_ff.topo_coeffs
        for k in ["BondBond Coeffs", "BondAngle Coeffs"]:
            self.assertIn(k, e_topo_coeffs["Angle Coeffs"][0], k)
        for k in ["MiddleBondTorsion Coeffs", "EndBondTorsion Coeffs",
                  "AngleTorsion Coeffs", "AngleAngleTorsion Coeffs",
                  "BondBond13 Coeffs"]:
            self.assertIn(k, e_topo_coeffs["Dihedral Coeffs"][0], k)
        self.assertIn("AngleAngle Coeffs",
                      e_topo_coeffs["Improper Coeffs"][0])

    def test_from_file(self):
        # general tests
        pep = self.peptide
        # header stats and Nos. of columns
        self.assertEqual(pep.masses.shape, (14, 1))
        self.assertEqual(pep.atoms.shape, (2004, 9))
        self.assertListEqual(list(pep.atoms.columns),
                             ["molecule-ID", "type", "q", "x", "y", "z",
                              "nx", "ny", "nz"])
        topo = pep.topology
        self.assertEqual(topo["Bonds"].shape, (1365, 3))
        self.assertEqual(topo["Angles"].shape, (786, 4))
        self.assertEqual(topo["Dihedrals"].shape, (207, 5))
        self.assertEqual(topo["Impropers"].shape, (12, 5))
        ff = pep.force_field
        self.assertEqual(ff["Pair Coeffs"].shape, (14, 4))
        self.assertEqual(ff["Bond Coeffs"].shape, (18, 2))
        self.assertEqual(ff["Angle Coeffs"].shape, (31, 4))
        self.assertEqual(ff["Dihedral Coeffs"].shape, (21, 4))
        self.assertEqual(ff["Improper Coeffs"].shape, (2, 2))
        # header box
        np.testing.assert_array_equal(pep.box.bounds,
                                      [[36.840194, 64.211560],
                                       [41.013691, 68.385058],
                                       [29.768095, 57.139462]])
        # body
        self.assertEqual(pep.masses.at[7, "mass"], 12.0110)
        self.assertEqual(ff["Pair Coeffs"].at[9, "coeff3"], 0.152100)
        self.assertEqual(ff["Bond Coeffs"].at[5, "coeff2"], 1.430000)
        self.assertEqual(ff["Angle Coeffs"].at[21, "coeff2"], 120.000000)
        self.assertEqual(ff["Dihedral Coeffs"].at[10, "coeff1"], 0.040000)
        self.assertEqual(ff["Improper Coeffs"].at[2, "coeff1"], 20.000000)
        self.assertEqual(pep.atoms.at[29, "molecule-ID"], 1)
        self.assertEqual(pep.atoms.at[29, "type"], 7)
        self.assertEqual(pep.atoms.at[29, "q"], -0.020)
        self.assertAlmostEqual(pep.atoms.at[29, "x"], 42.96709)
        self.assertEqual(pep.atoms.at[1808, "molecule-ID"], 576)
        self.assertEqual(pep.atoms.at[1808, "type"], 14)
        self.assertAlmostEqual(pep.atoms.at[1808, "y"], 58.64352)
        self.assertEqual(pep.atoms.at[1808, "nx"], -1)
        self.assertAlmostEqual(pep.velocities.at[527, "vz"], -0.010889)
        self.assertEqual(topo["Bonds"].at[47, "type"], 8)
        self.assertEqual(topo["Bonds"].at[47, "atom2"], 54)
        self.assertEqual(topo["Bonds"].at[953, "atom1"], 1384)
        self.assertEqual(topo["Angles"].at[105, "type"], 19)
        self.assertEqual(topo["Angles"].at[105, "atom3"], 51)
        self.assertEqual(topo["Angles"].at[376, "atom2"], 772)
        self.assertEqual(topo["Dihedrals"].at[151, "type"], 14)
        self.assertEqual(topo["Dihedrals"].at[151, "atom4"], 51)
        self.assertEqual(topo["Impropers"].at[4, "atom4"], 32)
        # class 2 and comments
        ethane = self.ethane
        self.assertEqual(ethane.masses.shape, (2, 1))
        self.assertEqual(ethane.atoms.shape, (8, 9))
        class2 = ethane.force_field
        self.assertEqual(class2["Pair Coeffs"].shape, (2, 2))
        self.assertEqual(class2["Bond Coeffs"].shape, (2, 4))
        self.assertEqual(class2["Angle Coeffs"].shape, (2, 4))
        self.assertEqual(class2["Dihedral Coeffs"].shape, (1, 6))
        self.assertEqual(class2["Improper Coeffs"].shape, (2, 2))
        self.assertEqual(class2["BondBond Coeffs"].at[2, "coeff3"], 1.1010)
        self.assertEqual(class2["BondAngle Coeffs"].at[2, "coeff4"], 1.1010)
        self.assertEqual(class2["AngleAngle Coeffs"].at[2, "coeff6"],
                         107.6600)
        self.assertEqual(class2["AngleAngle Coeffs"].at[2, "coeff6"],
                         107.6600)
        self.assertEqual(class2["AngleAngleTorsion Coeffs"].at[1, "coeff3"],
                         110.7700)
        self.assertEqual(class2["EndBondTorsion Coeffs"].at[1, "coeff8"],
                         1.1010)
        self.assertEqual(class2["MiddleBondTorsion Coeffs"].at[1, "coeff4"],
                         1.5300)
        self.assertEqual(class2["BondBond13 Coeffs"].at[1, "coeff3"], 1.1010)
        self.assertEqual(class2["AngleTorsion Coeffs"].at[1, "coeff8"],
                         110.7700)
        # tilt box and another atom_style
        quartz = self.quartz
        np.testing.assert_array_equal(quartz.box.tilt, [-2.456700, 0.0, 0.0])
        self.assertListEqual(list(quartz.atoms.columns),
                             ["type", "x", "y", "z"])
        self.assertAlmostEqual(quartz.atoms.at[7, "x"], 0.299963)
        # PairIJ Coeffs section
        virus = self.virus
        pairij = virus.force_field["PairIJ Coeffs"]
        self.assertEqual(pairij.at[7, "id1"], 3)
        self.assertEqual(pairij.at[7, "id2"], 3)
        self.assertEqual(pairij.at[7, "coeff2"], 2.1)
        # sort_id
        atom_id = random.randint(1, 384)
        self.assertEqual(self.tatb.atoms.loc[atom_id].name, atom_id)

    def test_from_ff_and_topologies(self):
        mass = OrderedDict()
        mass["H"] = 1.0079401
        mass["O"] = 15.999400
        nonbond_coeffs = [[0.00774378, 0.98], [0.1502629, 3.1169]]
        topo_coeffs = {"Bond Coeffs": [{"coeffs": [176.864, 0.9611],
                                        "types": [("H", "O")]}],
                       "Angle Coeffs": [{"coeffs": [42.1845, 109.4712],
                                         "types": [("H", "O", "H")]}]}
        ff = ForceField(mass.items(), nonbond_coeffs, topo_coeffs)
        with gzip.open(os.path.join(test_dir, "topologies_ice.json.gz")) as f:
            topo_dicts = json.load(f)
        topologies = [Topology.from_dict(d) for d in topo_dicts]
        box = LammpsBox([[-0.75694412, 44.165558],
                         [0.38127473, 47.066074],
                         [0.17900842, 44.193867]])
        ice = LammpsData.from_ff_and_topologies(box=box, ff=ff,
                                                topologies=topologies)
        atoms = ice.atoms
        bonds = ice.topology["Bonds"]
        angles = ice.topology["Angles"]
        np.testing.assert_array_equal(atoms.index.values,
                                      np.arange(1, len(atoms) + 1))
        np.testing.assert_array_equal(bonds.index.values,
                                      np.arange(1, len(bonds) + 1))
        np.testing.assert_array_equal(angles.index.values,
                                      np.arange(1, len(angles) + 1))

        i = random.randint(0, len(topologies) - 1)
        sample = topologies[i]
        in_atoms = ice.atoms[ice.atoms["molecule-ID"] == i + 1]
        np.testing.assert_array_equal(in_atoms.index.values,
                                      np.arange(3 * i + 1, 3 * i + 4))
        np.testing.assert_array_equal(in_atoms["type"].values, [2, 1, 1])
        np.testing.assert_array_equal(in_atoms["q"].values, sample.charges)
        np.testing.assert_array_equal(in_atoms[["x", "y", "z"]].values,
                                      sample.sites.cart_coords)
        broken_topo_coeffs = {"Bond Coeffs": [{"coeffs": [176.864, 0.9611],
                                               "types": [("H", "O")]}],
                              "Angle Coeffs": [{"coeffs": [42.1845, 109.4712],
                                                "types": [("H", "H", "H")]}]}
        broken_ff = ForceField(mass.items(), nonbond_coeffs,
                               broken_topo_coeffs)
        ld_woangles = LammpsData.from_ff_and_topologies(box=box, ff=broken_ff,
                                                        topologies=[sample])
        self.assertNotIn("Angles", ld_woangles.topology)

    def test_from_structure(self):
        latt = Lattice.monoclinic(9.78746, 4.75058, 8.95892, 115.9693)
        structure = Structure.from_spacegroup(15, latt, ["Os", "O", "O"],
                                              [[0, 0.25583, 0.75],
                                               [0.11146, 0.46611, 0.91631],
                                               [0.11445, 0.04564, 0.69518]])
        velocities = np.random.randn(20, 3) * 0.1
        structure.add_site_property("velocities", velocities)
        ld = LammpsData.from_structure(structure=structure,
                                       ff_elements=["O", "Os", "Na"])
        i = random.randint(0, 19)
        a = latt.matrix[0]
        va = velocities[i].dot(a) / np.linalg.norm(a)
        self.assertAlmostEqual(va, ld.velocities.loc[i + 1, "vx"])
        self.assertAlmostEqual(velocities[i, 1],
                               ld.velocities.loc[i + 1, "vy"])
        np.testing.assert_array_almost_equal(ld.masses["mass"],
                                             [22.989769, 190.23, 15.9994])
        np.testing.assert_array_equal(ld.atoms["type"], [2] * 4 + [3] * 16)

    def test_json_dict(self):
        encoded = json.dumps(self.ethane.as_dict())
        decoded = json.loads(encoded)
        c2h6 = LammpsData.from_dict(decoded)
        pd.testing.assert_frame_equal(c2h6.masses, self.ethane.masses)
        pd.testing.assert_frame_equal(c2h6.atoms, self.ethane.atoms)
        ff = self.ethane.force_field
        key, target_df = random.sample(ff.items(), 1)[0]
        self.assertIsNone(
            pd.testing.assert_frame_equal(c2h6.force_field[key], target_df,
                                          check_dtype=False),
            key
        )
        topo = self.ethane.topology
        key, target_df = random.sample(topo.items(), 1)[0]
        self.assertIsNone(
            pd.testing.assert_frame_equal(c2h6.topology[key], target_df),
            key
        )

    @classmethod
    def tearDownClass(cls):
        tmpfiles = ["test1.data", "test2.data"]
        for t in tmpfiles:
            if os.path.exists(t):
                os.remove(t)


class TopologyTest(unittest.TestCase):

    def test_init(self):
        inner_charge = np.random.rand(10) - 0.5
        outer_charge = np.random.rand(10) - 0.5
        inner_velo = np.random.rand(10, 3) - 0.5
        outer_velo = np.random.rand(10, 3) - 0.5
        m = Molecule(["H"] * 10, np.random.rand(10, 3) * 100,
                     site_properties={"ff_map": ["D"] * 10,
                                      "charge": inner_charge,
                                      "velocities": inner_velo})
        # q and v from site properties, while type from species_string
        topo = Topology(sites=m)
        self.assertListEqual(topo.type_by_sites, ["H"] * 10)
        np.testing.assert_array_equal(topo.charges, inner_charge)
        np.testing.assert_array_equal(topo.velocities, inner_velo)
        # q and v from overriding, while type from site property
        topo_override = Topology(sites=m, ff_label="ff_map",
                                 charges=outer_charge,
                                 velocities=outer_velo)
        self.assertListEqual(topo_override.type_by_sites, ["D"] * 10)
        np.testing.assert_array_equal(topo_override.charges, outer_charge)
        np.testing.assert_array_equal(topo_override.velocities, outer_velo)
        # test using a list of sites instead of SiteCollection
        topo_from_list = Topology(sites=m.sites)
        self.assertListEqual(topo_from_list.type_by_sites, topo.type_by_sites)
        np.testing.assert_array_equal(topo_from_list.charges, topo.charges)
        np.testing.assert_array_equal(topo_from_list.velocities,
                                      topo.velocities)

    def test_from_bonding(self):
        # He: no bonding topologies
        helium = Molecule(["He"], [[0, 0, 0]])
        topo_he = Topology.from_bonding(molecule=helium)
        self.assertIsNone(topo_he.topologies)
        # H2: 1 bond only
        hydrogen = Molecule(["H"] * 2, [[0, 0, 0], [0, 0, 0.7414]])
        topo_h = Topology.from_bonding(molecule=hydrogen)
        tp_h = topo_h.topologies
        self.assertListEqual(tp_h["Bonds"], [[0, 1]])
        self.assertNotIn("Angles", tp_h)
        self.assertNotIn("Dihedrals", tp_h)
        # water: 2 bonds and 1 angle only
        water = Molecule(["O", "H", "H"], [[0.0000, 0.0000, 0.1173],
                                           [0.0000, 0.7572, -0.4692],
                                           [0.0000, -0.7572, -0.4692]])
        topo_water = Topology.from_bonding(molecule=water)
        tp_water = topo_water.topologies
        self.assertListEqual(tp_water["Bonds"], [[0, 1], [0, 2]])
        self.assertListEqual(tp_water["Angles"], [[1, 0, 2]])
        self.assertNotIn("Dihedrals", tp_water)
        # EtOH
        etoh = Molecule(["C", "C", "O", "H", "H", "H", "H", "H", "H"],
                        [[1.1879, -0.3829, 0.0000],
                         [0.0000, 0.5526, 0.0000],
                         [-1.1867, -0.2472, 0.0000],
                         [-1.9237, 0.3850, 0.0000],
                         [2.0985, 0.2306, 0.0000],
                         [1.1184, -1.0093, 0.8869],
                         [1.1184, -1.0093, -0.8869],
                         [-0.0227, 1.1812, 0.8852],
                         [-0.0227, 1.1812, -0.8852]])
        topo_etoh = Topology.from_bonding(molecule=etoh)
        tp_etoh = topo_etoh.topologies
        self.assertEqual(len(tp_etoh["Bonds"]), 8)
        etoh_bonds = [[0, 1], [0, 4], [0, 5], [0, 6],
                      [1, 2], [1, 7], [1, 8], [2, 3]]
        np.testing.assert_array_equal(tp_etoh["Bonds"], etoh_bonds)
        self.assertEqual(len(tp_etoh["Angles"]), 13)
        etoh_angles = [[1, 0, 4], [1, 0, 5], [1, 0, 6], [4, 0, 5], [4, 0, 6],
                       [5, 0, 6], [0, 1, 2], [0, 1, 7], [0, 1, 8], [2, 1, 7],
                       [2, 1, 8], [7, 1, 8], [1, 2, 3]]
        np.testing.assert_array_equal(tp_etoh["Angles"], etoh_angles)
        self.assertEqual(len(tp_etoh["Dihedrals"]), 12)
        etoh_dihedrals = [[4, 0, 1, 2], [4, 0, 1, 7], [4, 0, 1, 8],
                          [5, 0, 1, 2], [5, 0, 1, 7], [5, 0, 1, 8],
                          [6, 0, 1, 2], [6, 0, 1, 7], [6, 0, 1, 8],
                          [0, 1, 2, 3], [7, 1, 2, 3], [8, 1, 2, 3]]
        np.testing.assert_array_equal(tp_etoh["Dihedrals"], etoh_dihedrals)
        self.assertIsNotNone(json.dumps(topo_etoh.as_dict()))
        # bond flag to off
        topo_etoh0 = Topology.from_bonding(molecule=etoh, bond=False,
                                           angle=True, dihedral=True)
        self.assertIsNone(topo_etoh0.topologies)
        # angle or dihedral flag to off
        topo_etoh1 = Topology.from_bonding(molecule=etoh, angle=False)
        self.assertNotIn("Angles", topo_etoh1.topologies)
        topo_etoh2 = Topology.from_bonding(molecule=etoh, dihedral=False)
        self.assertNotIn("Dihedrals", topo_etoh2.topologies)


class ForceFieldTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        mass_info = [("A", "H"), ("B", Element("C")),
                     ("C", Element("O")), ("D", 1.00794)]
        nonbond_coeffs = [[1, 1, 1.1225], [1, 1.175, 1.31894],
                          [1, 1.55, 1.73988], [1, 1, 1.1225 ],
                          [1, 1.35, 4], [1, 1.725, 1.93631],
                          [1, 1.175, 1.31894], [1, 2.1, 4],
                          [1, 1.55, 1.73988], [1, 1,  1.1225 ]]
        topo_coeffs = {"Bond Coeffs": [{"coeffs": [50, 0.659469],
                                        "types": [("A", "B"), ("C", "D")]},
                                       {"coeffs": [50, 0.855906],
                                        "types": [("B", "C")]}]}
        cls.virus = ForceField(mass_info=mass_info,
                               nonbond_coeffs=nonbond_coeffs,
                               topo_coeffs=topo_coeffs)
        cls.ethane = ForceField.from_file(os.path.join(test_dir,
                                                       "ff_ethane.yaml"))

    def test_init(self):
        v = self.virus
        self.assertListEqual(v.mass_info, [("A", 1.00794), ("B", 12.0107),
                                           ("C", 15.9994), ("D", 1.00794)])
        self.assertEqual(v.masses.at[3, "mass"], 15.9994)
        v_ff = v.force_field
        self.assertNotIn("Pair Coeffs", v_ff)
        self.assertEqual(v_ff["PairIJ Coeffs"].iat[5, 4], 1.93631)
        self.assertEqual(v_ff["Bond Coeffs"].at[2, "coeff2"], 0.855906)
        v_maps = v.maps
        self.assertDictEqual(v_maps["Atoms"],
                             {"A": 1, "B": 2, "C": 3, "D": 4})
        self.assertDictEqual(v_maps["Bonds"],
                             {("A", "B"): 1, ("C", "D"): 1,
                              ("B", "A"): 1, ("D", "C"): 1,
                              ("B", "C"): 2, ("C", "B"): 2})
        e = self.ethane
        self.assertEqual(e.masses.at[1, "mass"], 12.01115)
        e_ff = e.force_field
        self.assertNotIn("PairIJ Coeffs", e_ff)
        self.assertEqual(e_ff["Pair Coeffs"].at[1, "coeff2"], 3.854)
        self.assertEqual(e_ff["Bond Coeffs"].at[2, "coeff4"], 844.6)
        self.assertEqual(e_ff["Angle Coeffs"].at[2, "coeff4"], -2.4318)
        self.assertEqual(e_ff["Dihedral Coeffs"].at[1, "coeff1"], -0.1432)
        self.assertEqual(e_ff["Improper Coeffs"].at[2, "coeff2"], 0.0)
        self.assertEqual(e_ff["BondBond Coeffs"].at[2, "coeff1"], 5.3316)
        self.assertEqual(e_ff["BondAngle Coeffs"].at[1, "coeff3"], 1.53)
        self.assertEqual(e_ff["MiddleBondTorsion Coeffs"].at[1, "coeff1"],
                         -14.261)
        self.assertEqual(e_ff["EndBondTorsion Coeffs"].at[1, "coeff1"], 0.213)
        self.assertEqual(e_ff["AngleTorsion Coeffs"].at[1, "coeff3"], -0.2466)
        self.assertEqual(e_ff["AngleAngleTorsion Coeffs"].at[1, "coeff1"],
                         -12.564)
        self.assertEqual(e_ff["BondBond13 Coeffs"].at[1, "coeff1"], 0.0)
        self.assertEqual(e_ff["AngleAngle Coeffs"].at[1, "coeff2"], -0.4825)
        e_maps = e.maps
        self.assertDictEqual(e_maps["Atoms"], {"c4": 1, "h1": 2})
        self.assertDictEqual(e_maps["Bonds"],{("c4", "c4"): 1,
                                              ("c4", "h1"): 2,
                                              ("h1", "c4"): 2})
        self.assertDictEqual(e_maps["Angles"], {("c4", "c4", "h1"): 1,
                                                ("h1", "c4", "c4"): 1,
                                                ("h1", "c4", "h1"): 2})
        self.assertEqual(e_maps["Impropers"], {("c4", "c4", "h1", "h1"): 1,
                                               ("c4", "h1", "c4", "h1"): 1,
                                               ("h1", "h1", "c4", "c4"): 1,
                                               ("h1", "c4", "h1", "c4"): 1,
                                               ("h1", "c4", "h1", "h1"): 2,
                                               ("h1", "h1", "c4", "h1"): 2})

    def test_to_file(self):
        filename = "ff_test.yaml"
        v = self.virus
        v.to_file(filename=filename)
        yaml = YAML(typ="safe")
        with open(filename, "r") as f:
            d = yaml.load(f)
        self.assertListEqual(d["mass_info"], [list(m) for m in v.mass_info])
        self.assertListEqual(d["nonbond_coeffs"], v.nonbond_coeffs)

    def test_from_file(self):
        e = self.ethane
        self.assertListEqual(e.mass_info, [("c4", 12.01115), ("h1", 1.00797)])
        np.testing.assert_array_equal(e.nonbond_coeffs, [[0.062, 3.854],
                                                         [0.023, 2.878]])
        e_tc = e.topo_coeffs
        self.assertIn("Bond Coeffs", e_tc)
        self.assertIn("BondAngle Coeffs", e_tc["Angle Coeffs"][0])
        self.assertIn("BondBond Coeffs", e_tc["Angle Coeffs"][0])
        self.assertIn("AngleAngleTorsion Coeffs", e_tc["Dihedral Coeffs"][0])
        self.assertIn("AngleTorsion Coeffs", e_tc["Dihedral Coeffs"][0])
        self.assertIn("BondBond13 Coeffs", e_tc["Dihedral Coeffs"][0])
        self.assertIn("EndBondTorsion Coeffs", e_tc["Dihedral Coeffs"][0])
        self.assertIn("MiddleBondTorsion Coeffs", e_tc["Dihedral Coeffs"][0])
        self.assertIn("AngleAngle Coeffs", e_tc["Improper Coeffs"][0])

    def test_from_dict(self):
        d = self.ethane.as_dict()
        json_str = json.dumps(d)
        decoded = ForceField.from_dict(json.loads(json_str))
        self.assertListEqual(decoded.mass_info, self.ethane.mass_info)
        self.assertListEqual(decoded.nonbond_coeffs,
                             self.ethane.nonbond_coeffs)
        self.assertDictEqual(decoded.topo_coeffs, self.ethane.topo_coeffs)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists("ff_test.yaml"):
            os.remove("ff_test.yaml")


class FuncTest(unittest.TestCase):

    def test_lattice_2_lmpbox(self):
        matrix = np.diag(np.random.randint(5, 14, size=(3,))) \
                 + np.random.rand(3, 3) * 0.2 - 0.1
        init_latt = Lattice(matrix)
        frac_coords = np.random.rand(10, 3)
        init_structure = Structure(init_latt, ["H"] * 10, frac_coords)
        origin = np.random.rand(3) * 10 - 5
        box, symmop = lattice_2_lmpbox(lattice=init_latt, origin=origin)
        boxed_latt = box.to_lattice()
        np.testing.assert_array_almost_equal(init_latt.abc, boxed_latt.abc)
        np.testing.assert_array_almost_equal(init_latt.angles,
                                             boxed_latt.angles)
        cart_coords = symmop.operate_multi(init_structure.cart_coords) \
                      - origin
        boxed_structure = Structure(boxed_latt, ["H"] * 10, cart_coords,
                                    coords_are_cartesian=True)
        np.testing.assert_array_almost_equal(boxed_structure.frac_coords,
                                             frac_coords)
        tetra_latt = Lattice.tetragonal(5, 5)
        tetra_box, _ = lattice_2_lmpbox(tetra_latt)
        self.assertIsNone(tetra_box.tilt)
        ortho_latt = Lattice.orthorhombic(5, 5, 5)
        ortho_box, _ = lattice_2_lmpbox(ortho_latt)
        self.assertIsNone(ortho_box.tilt)
        rot_tetra_latt = Lattice([[5, 0, 0], [0, 2, 2], [0, -2, 2]])
        _, rotop = lattice_2_lmpbox(rot_tetra_latt)
        np.testing.\
            assert_array_almost_equal(rotop.rotation_matrix,
                                      [[1, 0, 0],
                                       [0, 2 ** 0.5 / 2, 2 ** 0.5 / 2],
                                       [0, -2 ** 0.5 / 2, 2 ** 0.5 / 2]])

    @unittest.skip("The function is deprecated")
    def test_structure_2_lmpdata(self):
        matrix = np.diag(np.random.randint(5, 14, size=(3,))) \
                 + np.random.rand(3, 3) * 0.2 - 0.1
        latt = Lattice(matrix)
        frac_coords = np.random.rand(10, 3)
        structure = Structure(latt, ["H"] * 10, frac_coords)
        ld = structure_2_lmpdata(structure=structure)
        box_tilt = [0.0, 0.0, 0.0] if not ld.box_tilt else ld.box_tilt
        box_bounds = np.array(ld.box_bounds)
        np.testing.assert_array_equal(box_bounds[:, 0], np.zeros(3))
        new_matrix = np.diag(box_bounds[:, 1])
        new_matrix[1, 0] = box_tilt[0]
        new_matrix[2, 0] = box_tilt[1]
        new_matrix[2, 1] = box_tilt[2]
        new_latt = Lattice(new_matrix)
        np.testing.assert_array_almost_equal(new_latt.abc, latt.abc)
        np.testing.assert_array_almost_equal(new_latt.angles, latt.angles)
        coords = ld.atoms[["x", "y", "z"]].values
        new_structure = Structure(new_latt, ['H'] * 10, coords,
                                  coords_are_cartesian=True)
        np.testing.assert_array_almost_equal(new_structure.frac_coords,
                                             frac_coords)
        self.assertEqual(len(ld.masses), 1)
        # test additional elements
        ld_elements = structure_2_lmpdata(structure=structure,
                                          ff_elements=["C", "H"])
        self.assertEqual(len(ld_elements.masses), 2)
        np.testing.assert_array_almost_equal(ld_elements.masses["mass"],
                                             [1.00794, 12.01070])


if __name__ == "__main__":
    unittest.main()
