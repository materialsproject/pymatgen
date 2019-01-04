# coding: utf-8


import os
import unittest

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.util.testing import PymatgenTest
from pymatgen.analysis.fragmenter import Fragmenter

try:
    import openbabel as ob
except ImportError:
    ob = None

__author__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files', 'fragmenter_files')


class TestFragmentMolecule(PymatgenTest):
    @classmethod
    def setUpClass(cls):
        cls.pc = Molecule.from_file(
            os.path.join(test_dir, "PC.xyz"))
        cls.pos_pc = Molecule.from_file(
            os.path.join(test_dir, "PC.xyz"))
        cls.pos_pc.set_charge_and_spin(charge=1)
        cls.pc_edges = [[5, 10], [5, 12], [5, 11], [5, 3], [3, 7], [3, 4],
                        [3, 0], [4, 8], [4, 9], [4, 1], [6, 1], [6, 0], [6, 2]]
        cls.pc_frag1 = Molecule.from_file(
            os.path.join(test_dir, "PC_frag1.xyz"))
        cls.pc_frag1_edges = [[0, 2], [4, 2], [2, 1], [1, 3]]
        cls.tfsi = Molecule.from_file(os.path.join(test_dir, "TFSI.xyz"))
        cls.tfsi_edges = [14,1],[1,4],[1,5],[1,7],[7,11],[7,12],[7,13],[14,0],[0,2],[0,3],[0,6],[6,8],[6,9],[6,10]

    def test_edges_given_PC_frag1(self):
        fragmenter = Fragmenter(molecule=self.pc_frag1, edges=self.pc_frag1_edges, depth=0)
        self.assertEqual(len(fragmenter.unique_fragments), 12)

    @unittest.skipIf(not ob, "OpenBabel not present. Skipping...")
    def test_babel_PC_frag1(self):
        fragmenter = Fragmenter(molecule=self.pc_frag1, depth=0)
        self.assertEqual(len(fragmenter.unique_fragments), 12)

    @unittest.skipIf(not ob, "OpenBabel not present. Skipping...")
    def test_babel_PC_defaults(self):
        fragmenter = Fragmenter(molecule=self.pc)
        self.assertEqual(fragmenter.open_rings,True)
        self.assertEqual(fragmenter.opt_steps,10000)
        default_mol_graph = MoleculeGraph.with_local_env_strategy(self.pc, OpenBabelNN(),
                                                                  reorder=False, extend_structure=False)
        self.assertEqual(fragmenter.mol_graph,default_mol_graph)
        self.assertEqual(len(fragmenter.unique_fragments), 13)
        self.assertEqual(len(fragmenter.unique_fragments_from_ring_openings), 5)

    def test_edges_given_PC_not_defaults(self):
        fragmenter = Fragmenter(molecule=self.pc, edges=self.pc_edges, depth=2, open_rings=False, opt_steps=0)
        self.assertEqual(fragmenter.open_rings,False)
        self.assertEqual(fragmenter.opt_steps,0)
        edges = {(e[0], e[1]): None for e in self.pc_edges}
        default_mol_graph = MoleculeGraph.with_edges(self.pc, edges=edges)
        self.assertEqual(fragmenter.mol_graph,default_mol_graph)
        self.assertEqual(len(fragmenter.unique_fragments), 20)
        self.assertEqual(len(fragmenter.unique_fragments_from_ring_openings), 0)

    def test_edges_given_TFSI(self):
        fragmenter = Fragmenter(molecule=self.tfsi, edges=self.tfsi_edges, depth=0)
        self.assertEqual(len(fragmenter.unique_fragments), 156)

    @unittest.skipIf(not ob, "OpenBabel not present. Skipping...")
    def test_babel_TFSI(self):
        fragmenter = Fragmenter(molecule=self.tfsi, depth=0)
        self.assertEqual(len(fragmenter.unique_fragments), 156)

    @unittest.skipIf(not ob, "OpenBabel not present. Skipping...")
    def test_babel_PC_with_RO_depth_0_vs_depth_10(self):
        fragmenter0 = Fragmenter(molecule=self.pc, depth=0, open_rings=True, opt_steps=1000)
        self.assertEqual(len(fragmenter0.unique_fragments), 509)

        fragmenter10 = Fragmenter(molecule=self.pc, depth=10, open_rings=True, opt_steps=1000)
        self.assertEqual(len(fragmenter10.unique_fragments), 509)

        fragments_by_level = fragmenter10.fragments_by_level
        num_frags_by_level = [13,51,95,115,105,75,39,14,2,0]
        for ii in range(10):
            self.assertEqual(len(fragments_by_level[str(ii)]),num_frags_by_level[ii])

        for fragment10 in fragmenter10.unique_fragments:
            found = False
            for fragment0 in fragmenter0.unique_fragments:
                if fragment0.isomorphic_to(fragment10):
                    found = True
                    break
            self.assertEqual(found, True)

    def test_PC_depth_0_vs_depth_10(self):
        fragmenter0 = Fragmenter(molecule=self.pc, edges=self.pc_edges, depth=0, open_rings=False)
        self.assertEqual(len(fragmenter0.unique_fragments), 295)

        fragmenter10 = Fragmenter(molecule=self.pc, edges=self.pc_edges, depth=10, open_rings=False)
        self.assertEqual(len(fragmenter10.unique_fragments), 63)

        fragments_by_level = fragmenter10.fragments_by_level
        num_frags_by_level = [8,12,15,14,9,4,1]
        for ii in range(7):
            self.assertEqual(len(fragments_by_level[str(ii)]),num_frags_by_level[ii])

        for fragment10 in fragmenter10.unique_fragments:
            found = False
            for fragment0 in fragmenter0.unique_fragments:
                if fragment0.isomorphic_to(fragment10):
                    found = True
                    break
            self.assertEqual(found, True)


if __name__ == "__main__":
    unittest.main()
