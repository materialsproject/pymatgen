# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import unittest
import os
import warnings
from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.analysis.functional_groups import FunctionalGroupExtractor

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "functional_groups")

try:
    import openbabel as ob
    import pybel as pb
    import networkx as nx
except ImportError:
    pb = None
    ob = None
    nx = None

__author__ = "Evan Spotte-Smith"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "ewcspottesmith@lbl.gov"
__status__ = "Beta"
__date__ = "July 2018"
__credit__ = "Peiyuan Yu"


@unittest.skipIf(not (pb and ob and nx), "OpenBabel or NetworkX not present. Skipping...")
class FunctionalGroupExtractorTest(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore")

        self.file = os.path.join(test_dir, "func_group_test.mol")
        self.mol = Molecule.from_file(self.file)
        self.strat = OpenBabelNN()
        self.mg = MoleculeGraph.with_local_env_strategy(self.mol, self.strat,
                                                        reorder=False,
                                                        extend_structure=False)
        self.extractor = FunctionalGroupExtractor(self.mg)

    def tearDown(self):
        warnings.simplefilter("default")
        del self.extractor
        del self.mg
        del self.strat
        del self.mol
        del self.file

    def test_init(self):
        # Ensure that instantiation is equivalent for all valid input types
        extractor_str = FunctionalGroupExtractor(self.file)
        extractor_mol = FunctionalGroupExtractor(self.mol)
        extractor_mg = self.extractor

        self.assertEqual(extractor_str.molgraph, extractor_mol.molgraph)
        self.assertEqual(extractor_str.molgraph, extractor_mg.molgraph)
        self.assertEqual(extractor_str.species, extractor_mol.species)
        self.assertEqual(extractor_str.species, extractor_mg.species)

        # Test optimization
        file_no_h = os.path.join(test_dir, "func_group_test_no_h.mol")
        extractor_no_h = FunctionalGroupExtractor(file_no_h, optimize=True)

        self.assertEqual(len(extractor_no_h.molecule), len(extractor_mol.molecule))
        self.assertEqual(extractor_no_h.species, extractor_mol.species)

    def test_get_heteroatoms(self):
        heteroatoms = self.extractor.get_heteroatoms()
        hetero_species = [self.extractor.species[x] for x in heteroatoms]

        self.assertEqual(len(heteroatoms), 3)
        self.assertEqual(sorted(hetero_species), ["N", "O", "O"])

        # Test with limitation
        hetero_no_o = self.extractor.get_heteroatoms(elements=["N"])
        self.assertEqual(len(hetero_no_o), 1)

    def test_get_special_carbon(self):
        special_cs = self.extractor.get_special_carbon()

        self.assertEqual(len(special_cs), 4)

        # Test with limitation
        special_cs_no_o = self.extractor.get_special_carbon(elements=["N"])
        self.assertEqual(len(special_cs_no_o), 2)

    def test_link_marked_atoms(self):
        heteroatoms = self.extractor.get_heteroatoms()
        special_cs = self.extractor.get_special_carbon()

        link = self.extractor.link_marked_atoms(heteroatoms.union(special_cs))

        self.assertEqual(len(link), 1)
        self.assertEqual(len(link[0]), 9)

        # Exclude Oxygen-related functional groups
        heteroatoms_no_o = self.extractor.get_heteroatoms(elements=["N"])
        special_cs_no_o = self.extractor.get_special_carbon(elements=["N"])
        all_marked = heteroatoms_no_o.union(special_cs_no_o)

        link_no_o = self.extractor.link_marked_atoms(all_marked)

        self.assertEqual(len(link_no_o), 2)

    def test_get_basic_functional_groups(self):
        basics = self.extractor.get_basic_functional_groups()

        # Molecule has one methyl group which will be caught.
        self.assertEqual(len(basics), 1)
        self.assertEqual(len(basics[0]), 4)

        basics_no_methyl = self.extractor.get_basic_functional_groups(func_groups=["phenyl"])
        self.assertEqual(len(basics_no_methyl), 0)

    def test_get_all_functional_groups(self):
        heteroatoms = self.extractor.get_heteroatoms()
        special_cs = self.extractor.get_special_carbon()

        link = self.extractor.link_marked_atoms(heteroatoms.union(special_cs))
        basics = self.extractor.get_basic_functional_groups()

        all_func = self.extractor.get_all_functional_groups()

        self.assertEqual(len(all_func), (len(link) + len(basics)))
        self.assertEqual(sorted(all_func), sorted(link + basics))

    def test_categorize_functional_groups(self):
        all_func = self.extractor.get_all_functional_groups()
        categorized = self.extractor.categorize_functional_groups(all_func)

        self.assertTrue("O=C1C=CC(=O)[N]1" in categorized.keys())
        self.assertTrue("[CH3]" in categorized.keys())

        total_count = sum([c["count"] for c in categorized.values()])
        self.assertEqual(total_count, 2)

if __name__ == "__main__":
    unittest.main()
