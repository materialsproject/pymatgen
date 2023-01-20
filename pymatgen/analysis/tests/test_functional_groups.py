# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from __future__ import annotations

import os
import unittest
import warnings

import pytest

from pymatgen.analysis.functional_groups import FunctionalGroupExtractor
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.core.structure import Molecule
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "functional_groups")

pytest.importorskip("openbabel", reason="OpenBabel not installed")
pytest.importorskip("networkx", reason="NetworkX not installed")

__author__ = "Evan Spotte-Smith"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "ewcspottesmith@lbl.gov"
__status__ = "Beta"
__date__ = "July 2018"
__credit__ = "Peiyuan Yu"


class FunctionalGroupExtractorTest(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore")

        self.file = os.path.join(test_dir, "func_group_test.mol")
        self.mol = Molecule.from_file(self.file)
        self.strategy = OpenBabelNN()
        self.mg = MoleculeGraph.with_local_env_strategy(self.mol, self.strategy)
        self.extractor = FunctionalGroupExtractor(self.mg)

    def tearDown(self):
        warnings.simplefilter("default")
        del self.extractor
        del self.mg
        del self.strategy
        del self.mol
        del self.file

    def test_init(self):
        # Ensure that instantiation is equivalent for all valid input types
        extractor_str = FunctionalGroupExtractor(self.file)
        extractor_mol = FunctionalGroupExtractor(self.mol)
        extractor_mg = self.extractor

        assert extractor_str.molgraph == extractor_mol.molgraph
        assert extractor_str.molgraph == extractor_mg.molgraph
        assert extractor_str.species == extractor_mol.species
        assert extractor_str.species == extractor_mg.species

        # Test optimization
        file_no_h = os.path.join(test_dir, "func_group_test_no_h.mol")
        extractor_no_h = FunctionalGroupExtractor(file_no_h, optimize=True)

        assert len(extractor_no_h.molecule) == len(extractor_mol.molecule)
        assert extractor_no_h.species == extractor_mol.species

    def test_get_heteroatoms(self):
        heteroatoms = self.extractor.get_heteroatoms()
        hetero_species = [self.extractor.species[x] for x in heteroatoms]

        assert len(heteroatoms) == 3
        assert sorted(hetero_species) == ["N", "O", "O"]

        # Test with limitation
        hetero_no_o = self.extractor.get_heteroatoms(elements=["N"])
        assert len(hetero_no_o) == 1

    def test_get_special_carbon(self):
        special_cs = self.extractor.get_special_carbon()

        assert len(special_cs) == 4

        # Test with limitation
        special_cs_no_o = self.extractor.get_special_carbon(elements=["N"])
        assert len(special_cs_no_o) == 2

    def test_link_marked_atoms(self):
        heteroatoms = self.extractor.get_heteroatoms()
        special_cs = self.extractor.get_special_carbon()

        link = self.extractor.link_marked_atoms(heteroatoms | special_cs)

        assert len(link) == 1
        assert len(link[0]) == 9

        # Exclude Oxygen-related functional groups
        heteroatoms_no_o = self.extractor.get_heteroatoms(elements=["N"])
        special_cs_no_o = self.extractor.get_special_carbon(elements=["N"])
        all_marked = heteroatoms_no_o | special_cs_no_o

        link_no_o = self.extractor.link_marked_atoms(all_marked)

        assert len(link_no_o) == 2

    def test_get_basic_functional_groups(self):
        basics = self.extractor.get_basic_functional_groups()

        # Molecule has one methyl group which will be caught.
        assert len(basics) == 1
        assert len(basics[0]) == 4

        basics_no_methyl = self.extractor.get_basic_functional_groups(func_groups=["phenyl"])
        assert len(basics_no_methyl) == 0

    def test_get_all_functional_groups(self):
        heteroatoms = self.extractor.get_heteroatoms()
        special_cs = self.extractor.get_special_carbon()

        link = self.extractor.link_marked_atoms(heteroatoms | special_cs)
        basics = self.extractor.get_basic_functional_groups()

        all_func = self.extractor.get_all_functional_groups()

        assert len(all_func) == (len(link) + len(basics))
        assert sorted(all_func) == sorted(link + basics)

    def test_categorize_functional_groups(self):
        all_func = self.extractor.get_all_functional_groups()
        categorized = self.extractor.categorize_functional_groups(all_func)

        assert "O=C1C=CC(=O)[N]1" in categorized
        assert "[CH3]" in categorized

        total_count = sum(c["count"] for c in categorized.values())
        assert total_count == 2


if __name__ == "__main__":
    unittest.main()
