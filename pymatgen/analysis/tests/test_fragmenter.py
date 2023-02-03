from __future__ import annotations

import os
import unittest

import pytest

from pymatgen.analysis.fragmenter import Fragmenter
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.core.structure import Molecule
from pymatgen.util.testing import PymatgenTest

__author__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"


test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "fragmenter_files")


class TestFragmentMolecule(PymatgenTest):
    @classmethod
    def setUpClass(cls):
        cls.pc = Molecule.from_file(os.path.join(test_dir, "PC.xyz"))
        cls.ec = Molecule.from_file(os.path.join(test_dir, "EC.xyz"))
        cls.pos_pc = Molecule.from_file(os.path.join(test_dir, "PC.xyz"))
        cls.pos_pc.set_charge_and_spin(charge=1)
        cls.pc_edges = [
            [5, 10],
            [5, 12],
            [5, 11],
            [5, 3],
            [3, 7],
            [3, 4],
            [3, 0],
            [4, 8],
            [4, 9],
            [4, 1],
            [6, 1],
            [6, 0],
            [6, 2],
        ]
        cls.pc_frag1 = Molecule.from_file(os.path.join(test_dir, "PC_frag1.xyz"))
        cls.pc_frag1_edges = [[0, 2], [4, 2], [2, 1], [1, 3]]
        cls.tfsi = Molecule.from_file(os.path.join(test_dir, "TFSI.xyz"))
        cls.tfsi_edges = (
            [14, 1],
            [1, 4],
            [1, 5],
            [1, 7],
            [7, 11],
            [7, 12],
            [7, 13],
            [14, 0],
            [0, 2],
            [0, 3],
            [0, 6],
            [6, 8],
            [6, 9],
            [6, 10],
        )
        cls.LiEC = Molecule.from_file(os.path.join(test_dir, "LiEC.xyz"))

    def test_edges_given_PC_frag1(self):
        fragmenter = Fragmenter(molecule=self.pc_frag1, edges=self.pc_frag1_edges, depth=0)
        assert fragmenter.total_unique_fragments == 12

    def test_babel_PC_frag1(self):
        pytest.importorskip("openbabel", reason="OpenBabel not installed")
        fragmenter = Fragmenter(molecule=self.pc_frag1, depth=0)
        assert fragmenter.total_unique_fragments == 12

    def test_babel_PC_old_defaults(self):
        pytest.importorskip("openbabel", reason="OpenBabel not installed")
        fragmenter = Fragmenter(molecule=self.pc, open_rings=True)
        assert fragmenter.open_rings is True
        assert fragmenter.opt_steps == 10000
        default_mol_graph = MoleculeGraph.with_local_env_strategy(self.pc, OpenBabelNN())
        assert fragmenter.mol_graph == default_mol_graph
        assert fragmenter.total_unique_fragments == 13

    def test_babel_PC_defaults(self):
        pytest.importorskip("openbabel", reason="OpenBabel not installed")
        fragmenter = Fragmenter(molecule=self.pc)
        assert fragmenter.open_rings is False
        assert fragmenter.opt_steps == 10000
        default_mol_graph = MoleculeGraph.with_local_env_strategy(self.pc, OpenBabelNN())
        assert fragmenter.mol_graph == default_mol_graph
        assert fragmenter.total_unique_fragments == 8

    def test_edges_given_PC_not_defaults(self):
        fragmenter = Fragmenter(
            molecule=self.pc,
            edges=self.pc_edges,
            depth=2,
            open_rings=False,
            opt_steps=0,
        )
        assert fragmenter.open_rings is False
        assert fragmenter.opt_steps == 0
        edges = {(e[0], e[1]): None for e in self.pc_edges}
        default_mol_graph = MoleculeGraph.with_edges(self.pc, edges=edges)
        assert fragmenter.mol_graph == default_mol_graph
        assert fragmenter.total_unique_fragments == 20

    def test_edges_given_TFSI(self):
        fragmenter = Fragmenter(molecule=self.tfsi, edges=self.tfsi_edges, depth=0)
        assert fragmenter.total_unique_fragments == 156

    def test_babel_TFSI(self):
        pytest.importorskip("openbabel", reason="OpenBabel not installed")
        fragmenter = Fragmenter(molecule=self.tfsi, depth=0)
        assert fragmenter.total_unique_fragments == 156

    def test_babel_PC_with_RO_depth_0_vs_depth_10(self):
        pytest.importorskip("openbabel", reason="OpenBabel not installed")
        fragmenter0 = Fragmenter(molecule=self.pc, depth=0, open_rings=True, opt_steps=1000)
        assert fragmenter0.total_unique_fragments == 509

        fragmenter10 = Fragmenter(molecule=self.pc, depth=10, open_rings=True, opt_steps=1000)
        assert fragmenter10.total_unique_fragments == 509

        fragments_by_level = fragmenter10.fragments_by_level
        num_frags_by_level = [13, 51, 95, 115, 105, 75, 39, 14, 2, 0]
        for ii in range(10):
            num_frags = 0
            for key in fragments_by_level[str(ii)]:
                num_frags += len(fragments_by_level[str(ii)][key])
            assert num_frags == num_frags_by_level[ii]

    def test_PC_depth_0_vs_depth_10(self):
        fragmenter0 = Fragmenter(molecule=self.pc, edges=self.pc_edges, depth=0, open_rings=False)
        assert fragmenter0.total_unique_fragments == 295

        fragmenter10 = Fragmenter(molecule=self.pc, edges=self.pc_edges, depth=10, open_rings=False)
        assert fragmenter10.total_unique_fragments == 63

        fragments_by_level = fragmenter10.fragments_by_level
        num_frags_by_level = [8, 12, 15, 14, 9, 4, 1]
        for ii in range(7):
            num_frags = 0
            for key in fragments_by_level[str(ii)]:
                num_frags += len(fragments_by_level[str(ii)][key])
            assert num_frags == num_frags_by_level[ii]

    def test_PC_frag1_then_PC(self):
        frag1 = Fragmenter(molecule=self.pc_frag1, edges=self.pc_frag1_edges, depth=0)
        assert frag1.new_unique_fragments == frag1.total_unique_fragments
        frag2 = Fragmenter(
            molecule=self.pc,
            edges=self.pc_edges,
            depth=0,
            open_rings=False,
            prev_unique_frag_dict=frag1.unique_frag_dict,
        )
        assert frag2.new_unique_fragments == 295 - 12

    def test_PC_then_EC_depth_10(self):
        pytest.importorskip("openbabel", reason="OpenBabel not installed")
        fragPC = Fragmenter(molecule=self.pc, depth=10, open_rings=True)
        fragEC = Fragmenter(
            molecule=self.ec,
            depth=10,
            open_rings=True,
            prev_unique_frag_dict=fragPC.unique_frag_dict,
        )
        assert fragEC.new_unique_fragments == 11
        assert fragEC.total_unique_fragments == 509 + 11


if __name__ == "__main__":
    unittest.main()
