# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import copy
import os
import unittest
import warnings

import pytest

from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.molecule_matcher import MoleculeMatcher
from pymatgen.core.structure import Molecule
from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.io.xyz import XYZ
from pymatgen.util.testing import PymatgenTest

pytest.importorskip("openbabel", reason="OpenBabel not installed")


class BabelMolAdaptorTest(unittest.TestCase):
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

    def test_init(self):
        adaptor = BabelMolAdaptor(self.mol)
        obmol = adaptor.openbabel_mol
        pybel = adaptor.pybel_mol
        self.assertEqual(obmol.NumAtoms(), 5)

        adaptor = BabelMolAdaptor(adaptor.openbabel_mol)
        self.assertEqual(adaptor.pymatgen_mol.formula, "H4 C1")

        adaptor = BabelMolAdaptor(pybel)
        self.assertEqual(adaptor.pymatgen_mol.formula, "H4 C1")

    def test_from_file(self):
        adaptor = BabelMolAdaptor.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules/Ethane_e.pdb"), "pdb")
        mol = adaptor.pymatgen_mol
        self.assertEqual(mol.formula, "H6 C2")

    def test_from_file_return_all_molecules(self):
        adaptors = BabelMolAdaptor.from_file(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "multiple_frame_xyz.xyz"),
            "xyz",
            return_all_molecules=True,
        )
        self.assertEqual(len(adaptors), 302)

    def test_from_molecule_graph(self):
        graph = MoleculeGraph.with_empty_graph(self.mol)
        adaptor = BabelMolAdaptor.from_molecule_graph(graph)
        obmol = adaptor.openbabel_mol
        self.assertEqual(obmol.NumAtoms(), 5)
        mol = adaptor.pymatgen_mol
        self.assertEqual(mol.formula, "H4 C1")

    def test_from_string(self):
        xyz = XYZ(self.mol)
        adaptor = BabelMolAdaptor.from_string(str(xyz), "xyz")
        mol = adaptor.pymatgen_mol
        self.assertEqual(mol.formula, "H4 C1")

    def test_localopt(self):
        self.mol[1] = "H", [0, 0, 1.05]
        adaptor = BabelMolAdaptor(self.mol)
        adaptor.localopt()
        optmol = adaptor.pymatgen_mol
        for site in optmol[1:]:
            self.assertAlmostEqual(site.distance(optmol[0]), 1.09216, 1)

    def test_make3d(self):
        from openbabel import pybel as pb

        mol_0d = pb.readstring("smi", "CCCC").OBMol
        adaptor = BabelMolAdaptor(mol_0d)
        adaptor.make3d()
        self.assertEqual(mol_0d.GetDimension(), 3)

    def add_hydrogen(self):
        from openbabel import pybel as pb

        mol_0d = pb.readstring("smi", "CCCC").OBMol
        self.assertEqual(len(pb.Molecule(mol_0d).atoms), 2)
        adaptor = BabelMolAdaptor(mol_0d)
        adaptor.add_hydrogen()
        self.assertEqual(len(adaptor.pymatgen_mol.sites), 14)

    def test_rotor_search_wrs(self):
        mol = copy.deepcopy(self.mol)
        mol[1] = "H", [0, 0, 1.05]
        adaptor = BabelMolAdaptor(mol)
        rotor_args = (250, 50)
        adaptor.rotor_conformer(*rotor_args, algo="WeightedRotorSearch")
        optmol = adaptor.pymatgen_mol
        for site in optmol[1:]:
            self.assertAlmostEqual(site.distance(optmol[0]), 1.09216, 1)

    def test_rotor_search_srs(self):
        mol = copy.deepcopy(self.mol)
        mol[1] = "H", [0, 0, 1.05]
        adaptor = BabelMolAdaptor(mol)
        adaptor.rotor_conformer(200, algo="SystematicRotorSearch")
        optmol = adaptor.pymatgen_mol
        for site in optmol[1:]:
            self.assertAlmostEqual(site.distance(optmol[0]), 1.09216, 1)

    def test_rotor_search_rrs(self):
        mol = copy.deepcopy(self.mol)
        mol[1] = "H", [0, 0, 1.05]
        adaptor = BabelMolAdaptor(mol)
        adaptor.rotor_conformer(250, 50, algo="RandomRotorSearch")
        optmol = adaptor.pymatgen_mol
        for site in optmol[1:]:
            self.assertAlmostEqual(site.distance(optmol[0]), 1.09216, 1)

    def test_confab_conformers(self):
        from openbabel import pybel as pb

        mol = pb.readstring("smi", "CCCC").OBMol
        adaptor = BabelMolAdaptor(mol)
        adaptor.make3d()
        conformers = adaptor.confab_conformers()
        self.assertEqual(adaptor.openbabel_mol.NumRotors(), 1)
        self.assertGreaterEqual(len(conformers), 1)
        if len(conformers) > 1:
            self.assertNotAlmostEqual(MoleculeMatcher().get_rmsd(conformers[0], conformers[1]), 0)


if __name__ == "__main__":
    unittest.main()
