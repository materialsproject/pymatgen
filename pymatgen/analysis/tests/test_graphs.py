# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import unittest
import os
import copy
import matplotlib

from pymatgen.command_line.critic2_caller import Critic2Output
from pymatgen.core.structure import Molecule, Structure, FunctionalGroups, Site
from pymatgen.analysis.graphs import *
from pymatgen.analysis.local_env import MinimumDistanceNN, MinimumOKeeffeNN

try:
    import openbabel as ob
    import networkx as nx
except ImportError:
    ob = None
    nx = None

__author__ = "Matthew Horton, Evan Spotte-Smith"
__version__ = "0.1"
__maintainer__ = "Matthew Horton"
__email__ = "mkhorton@lbl.gov"
__status__ = "Beta"
__date__ = "August 2017"


class StructureGraphTest(unittest.TestCase):

    def setUp(self):

        self.maxDiff = None

        # trivial example, simple square lattice for testing
        structure = Structure(Lattice.tetragonal(5.0, 50.0), ['H'], [[0, 0, 0]])
        self.square_sg = StructureGraph.with_empty_graph(structure, edge_weight_name="", edge_weight_units="")
        self.square_sg.add_edge(0, 0, from_jimage=(0, 0, 0), to_jimage=(1, 0, 0))
        self.square_sg.add_edge(0, 0, from_jimage=(0, 0, 0), to_jimage=(-1, 0, 0))
        self.square_sg.add_edge(0, 0, from_jimage=(0, 0, 0), to_jimage=(0, 1, 0))
        self.square_sg.add_edge(0, 0, from_jimage=(0, 0, 0), to_jimage=(0, -1, 0))
        #TODO: decorating still fails because the structure graph gives a CN of 8 for this square lattice
        # self.square_sg.decorate_structure_with_ce_info()

        # body-centered square lattice for testing
        structure = Structure(Lattice.tetragonal(5.0, 50.0), ['H', 'He'], [[0, 0, 0], [0.5, 0.5, 0.5]])
        self.bc_square_sg = StructureGraph.with_empty_graph(structure, edge_weight_name="", edge_weight_units="")
        self.bc_square_sg.add_edge(0, 0, from_jimage=(0, 0, 0), to_jimage=(1, 0, 0))
        self.bc_square_sg.add_edge(0, 0, from_jimage=(0, 0, 0), to_jimage=(-1, 0, 0))
        self.bc_square_sg.add_edge(0, 0, from_jimage=(0, 0, 0), to_jimage=(0, 1, 0))
        self.bc_square_sg.add_edge(0, 0, from_jimage=(0, 0, 0), to_jimage=(0, -1, 0))
        self.bc_square_sg.add_edge(0, 1, from_jimage=(0, 0, 0), to_jimage=(0, 0, 0))
        self.bc_square_sg.add_edge(0, 1, from_jimage=(0, 0, 0), to_jimage=(-1, 0, 0))
        self.bc_square_sg.add_edge(0, 1, from_jimage=(0, 0, 0), to_jimage=(-1, -1, 0))
        self.bc_square_sg.add_edge(0, 1, from_jimage=(0, 0, 0), to_jimage=(0, -1, 0))

        # body-centered square lattice for testing
        # directions reversed, should be equivalent to as bc_square
        structure = Structure(Lattice.tetragonal(5.0, 50.0), ['H', 'He'], [[0, 0, 0], [0.5, 0.5, 0.5]])
        self.bc_square_sg_r = StructureGraph.with_empty_graph(structure, edge_weight_name="", edge_weight_units="")
        self.bc_square_sg_r.add_edge(0, 0, from_jimage=(0, 0, 0), to_jimage=(1, 0, 0))
        self.bc_square_sg_r.add_edge(0, 0, from_jimage=(0, 0, 0), to_jimage=(-1, 0, 0))
        self.bc_square_sg_r.add_edge(0, 0, from_jimage=(0, 0, 0), to_jimage=(0, 1, 0))
        self.bc_square_sg_r.add_edge(0, 0, from_jimage=(0, 0, 0), to_jimage=(0, -1, 0))
        self.bc_square_sg_r.add_edge(0, 1, from_jimage=(0, 0, 0), to_jimage=(0, 0, 0))
        self.bc_square_sg_r.add_edge(1, 0, from_jimage=(-1, 0, 0), to_jimage=(0, 0, 0))
        self.bc_square_sg_r.add_edge(1, 0, from_jimage=(-1, -1, 0), to_jimage=(0, 0, 0))
        self.bc_square_sg_r.add_edge(1, 0, from_jimage=(0, -1, 0), to_jimage=(0, 0, 0))

        # MoS2 example, structure graph obtained from critic2
        # (not ground state, from mp-1023924, single layer)
        stdout_file = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                   'test_files/critic2/MoS2_critic2_stdout.txt')
        with open(stdout_file, 'r') as f:
            reference_stdout = f.read()
        self.structure = Structure.from_file(os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                             'test_files/critic2/MoS2.cif'))
        c2o = Critic2Output(self.structure, reference_stdout)
        self.mos2_sg = c2o.structure_graph(edge_weight="bond_length", edge_weight_units="Å")

        latt = Lattice.cubic(4.17)
        species = ["Ni", "O"]
        coords = [[0, 0, 0],
                  [0.5, 0.5, 0.5]]
        self.NiO = Structure.from_spacegroup(225, latt, species, coords).get_primitive_structure()

        # BCC example.
        self.bcc = Structure(Lattice.cubic(5.0), ['He', 'He'], [[0, 0, 0], [0.5, 0.5, 0.5]])


        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.resetwarnings()

    def test_properties(self):

        self.assertEqual(self.mos2_sg.name, "bonds")
        self.assertEqual(self.mos2_sg.edge_weight_name, "bond_length")
        self.assertEqual(self.mos2_sg.edge_weight_unit, "Å")
        self.assertEqual(self.mos2_sg.get_coordination_of_site(0), 6)
        self.assertEqual(len(self.mos2_sg.get_connected_sites(0)), 6)
        self.assertTrue(isinstance(self.mos2_sg.get_connected_sites(0)[0].site, PeriodicSite))
        self.assertEqual(str(self.mos2_sg.get_connected_sites(0)[0].site.specie), 'S')

        # these two graphs should be equivalent
        for n in range(len(self.bc_square_sg)):
            self.assertEqual(self.bc_square_sg.get_coordination_of_site(n),
                             self.bc_square_sg_r.get_coordination_of_site(n))

    def test_auto_image_detection(self):

        sg = StructureGraph.with_empty_graph(self.structure)
        sg.add_edge(0, 0)

        ref_edges = [(0, 0, {'to_jimage': (-1, -1, 0)}),
                     (0, 0, {'to_jimage': (-1, 0, 0)}),
                     (0, 0, {'to_jimage': (0, -1, 0)}),
                     (0, 0, {'to_jimage': (0, 0, 0)}),
                     (0, 0, {'to_jimage': (1, 0, 0)})]

        self.assertEqual(list(sg.graph.edges(data=True)), ref_edges)

    def test_str(self):

        square_sg_str_ref = """Structure Graph
Structure: 
Full Formula (H1)
Reduced Formula: H2
abc   :   5.000000   5.000000  50.000000
angles:  90.000000  90.000000  90.000000
Sites (1)
  #  SP      a    b    c
---  ----  ---  ---  ---
  0  H       0    0    0
Graph: bonds
from    to  to_image    
----  ----  ------------
   0     0  (1, 0, 0)   
   0     0  (-1, 0, 0)  
   0     0  (0, 1, 0)   
   0     0  (0, -1, 0)  
"""

        mos2_sg_str_ref = """Structure Graph
Structure: 
Full Formula (Mo1 S2)
Reduced Formula: MoS2
abc   :   3.190316   3.190315  17.439502
angles:  90.000000  90.000000 120.000006
Sites (3)
  #  SP           a         b         c
---  ----  --------  --------  --------
  0  Mo    0.333333  0.666667  0.213295
  1  S     0.666667  0.333333  0.303027
  2  S     0.666667  0.333333  0.123562
Graph: bonds
from    to  to_image      bond_length (A)
----  ----  ------------  ------------------
   0     1  (-1, 0, 0)    2.417e+00
   0     1  (0, 0, 0)     2.417e+00
   0     1  (0, 1, 0)     2.417e+00
   0     2  (0, 1, 0)     2.417e+00
   0     2  (-1, 0, 0)    2.417e+00
   0     2  (0, 0, 0)     2.417e+00
"""

        # don't care about testing Py 2.7 unicode support,
        # change Å to A
        self.mos2_sg.graph.graph['edge_weight_units'] = "A"
        self.assertEqual(str(self.square_sg), square_sg_str_ref)
        self.assertEqual(str(self.mos2_sg), mos2_sg_str_ref)

    def test_mul(self):

        square_sg_mul = self.square_sg * (2, 1, 1)

        square_sg_mul_ref_str = """Structure Graph
Structure: 
Full Formula (H2)
Reduced Formula: H2
abc   :  10.000000   5.000000  50.000000
angles:  90.000000  90.000000  90.000000
Sites (2)
  #  SP      a    b    c
---  ----  ---  ---  ---
  0  H     0      0    0
  1  H     0.5    0   -0
Graph: bonds
from    to  to_image    
----  ----  ------------
   0     0  (0, 1, 0)   
   0     0  (0, -1, 0)  
   0     1  (0, 0, 0)   
   0     1  (-1, 0, 0)  
   1     1  (0, 1, 0)   
   1     1  (0, -1, 0)  
"""
        square_sg_mul_actual_str = str(square_sg_mul)

        # only testing bonds portion,
        # the c frac_coord of the second H can vary from
        # 0 to -0 depending on machine precision
        square_sg_mul_ref_str = "\n".join(square_sg_mul_ref_str.splitlines()[11:])
        square_sg_mul_actual_str = "\n".join(square_sg_mul_actual_str.splitlines()[11:])

        self.assertEqual(square_sg_mul_actual_str, square_sg_mul_ref_str)

        # test sequential multiplication
        sq_sg_1 = self.square_sg*(2, 2, 1)
        sq_sg_1 = sq_sg_1*(2, 2, 1)
        sq_sg_2 = self.square_sg*(4, 4, 1)
        self.assertEqual(sq_sg_1.graph.number_of_edges(), sq_sg_2.graph.number_of_edges())
        #TODO: the below test still gives 8 != 4
        #self.assertEqual(self.square_sg.get_coordination_of_site(0), 4)

        mos2_sg_mul = self.mos2_sg * (3, 3, 1)
        for idx in mos2_sg_mul.structure.indices_from_symbol("Mo"):
            self.assertEqual(mos2_sg_mul.get_coordination_of_site(idx), 6)

        mos2_sg_premul = StructureGraph.with_local_env_strategy(self.structure*(3, 3, 1),
                                                                MinimumDistanceNN())
        self.assertTrue(mos2_sg_mul == mos2_sg_premul)

        # test 3D Structure

        nio_sg = StructureGraph.with_local_env_strategy(self.NiO, MinimumDistanceNN())
        nio_sg = nio_sg*3

        for n in range(len(nio_sg)):
            self.assertEqual(nio_sg.get_coordination_of_site(n), 6)


    @unittest.skipIf(not (which('neato') and which('fdp')), "graphviz executables not present")
    def test_draw(self):

        # draw MoS2 graph
        self.mos2_sg.draw_graph_to_file('MoS2_single.pdf', image_labels=True,
                                        hide_image_edges=False)
        mos2_sg = self.mos2_sg * (9, 9, 1)
        mos2_sg.draw_graph_to_file('MoS2.pdf', algo='neato')

        # draw MoS2 graph that's been successively multiplied
        mos2_sg_2 = self.mos2_sg * (3, 3, 1)
        mos2_sg_2 = mos2_sg_2 * (3, 3, 1)
        mos2_sg_2.draw_graph_to_file('MoS2_twice_mul.pdf', algo='neato', hide_image_edges=True)

        # draw MoS2 graph that's generated from a pre-multiplied Structure
        mos2_sg_premul = StructureGraph.with_local_env_strategy(self.structure*(3, 3, 1),
                                                                MinimumDistanceNN())
        mos2_sg_premul.draw_graph_to_file('MoS2_premul.pdf', algo='neato', hide_image_edges=True)

        # draw graph for a square lattice
        self.square_sg.draw_graph_to_file('square_single.pdf', hide_image_edges=False)
        square_sg = self.square_sg * (5, 5, 1)
        square_sg.draw_graph_to_file('square.pdf', algo='neato', image_labels=True,
                                     node_labels=False)

        # draw graph for a body-centered square lattice
        self.bc_square_sg.draw_graph_to_file('bc_square_single.pdf', hide_image_edges=False)
        bc_square_sg = self.bc_square_sg * (9, 9, 1)
        bc_square_sg.draw_graph_to_file('bc_square.pdf', algo='neato', image_labels=False)

        # draw graph for a body-centered square lattice defined in an alternative way
        self.bc_square_sg_r.draw_graph_to_file('bc_square_r_single.pdf', hide_image_edges=False)
        bc_square_sg_r = self.bc_square_sg_r * (9, 9, 1)
        bc_square_sg_r.draw_graph_to_file('bc_square_r.pdf', algo='neato', image_labels=False)

        # delete generated test files
        test_files = ('bc_square_r_single.pdf', 'bc_square_r.pdf', 'bc_square_single.pdf',
                      'bc_square.pdf', 'MoS2_premul.pdf', 'MOS2_single.pdf', 'MoS2_twice_mul.pdf',
                      'MoS2.pdf', 'square_single.pdf', 'square.pdf')
        for test_file in test_files:
            os.remove(test_file)

    def test_to_from_dict(self):
        d = self.mos2_sg.as_dict()
        sg = StructureGraph.from_dict(d)
        d2 = sg.as_dict()
        self.assertDictEqual(d, d2)

    def test_from_local_env_and_equality_and_diff(self):
        nn = MinimumDistanceNN()
        sg = StructureGraph.with_local_env_strategy(self.structure, nn)

        self.assertEqual(sg.graph.number_of_edges(), 6)

        nn2 = MinimumOKeeffeNN()
        sg2 = StructureGraph.with_local_env_strategy(self.structure, nn2)

        self.assertTrue(sg == sg2)
        self.assertTrue(sg == self.mos2_sg)

        # TODO: find better test case where graphs are different
        diff = sg.diff(sg2)
        self.assertEqual(diff['dist'], 0)

        self.assertEqual(self.square_sg.get_coordination_of_site(0), 4)

    def test_extract_molecules(self):

        structure_file = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                      'test_files/C26H16BeN2O2S2.cif')

        s = Structure.from_file(structure_file)

        nn = MinimumOKeeffeNN()
        sg = StructureGraph.with_local_env_strategy(s, nn)

        molecules = sg.get_subgraphs_as_molecules()

        self.assertEqual(molecules[0].composition.formula, "Be1 H16 C26 S2 N2 O2")
        self.assertEqual(len(molecules), 1)

        molecules = self.mos2_sg.get_subgraphs_as_molecules()
        self.assertEqual(len(molecules), 0)


class MoleculeGraphTest(unittest.TestCase):

    def setUp(self):

        cyclohexene = Molecule.from_file(os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                                      "test_files/graphs/cyclohexene.xyz"))
        self.cyclohexene = MoleculeGraph.with_empty_graph(cyclohexene,
                                                       edge_weight_name="strength",
                                                       edge_weight_units="")
        self.cyclohexene.add_edge(0, 1, weight=1.0)
        self.cyclohexene.add_edge(1, 2, weight=1.0)
        self.cyclohexene.add_edge(2, 3, weight=2.0)
        self.cyclohexene.add_edge(3, 4, weight=1.0)
        self.cyclohexene.add_edge(4, 5, weight=1.0)
        self.cyclohexene.add_edge(5, 0, weight=1.0)
        self.cyclohexene.add_edge(0, 6, weight=1.0)
        self.cyclohexene.add_edge(0, 7, weight=1.0)
        self.cyclohexene.add_edge(1, 8, weight=1.0)
        self.cyclohexene.add_edge(1, 9, weight=1.0)
        self.cyclohexene.add_edge(2, 10, weight=1.0)
        self.cyclohexene.add_edge(3, 11, weight=1.0)
        self.cyclohexene.add_edge(4, 12, weight=1.0)
        self.cyclohexene.add_edge(4, 13, weight=1.0)
        self.cyclohexene.add_edge(5, 14, weight=1.0)
        self.cyclohexene.add_edge(5, 15, weight=1.0)

        butadiene = Molecule.from_file(os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                                    "test_files/graphs/butadiene.xyz"))
        self.butadiene = MoleculeGraph.with_empty_graph(butadiene,
                                                        edge_weight_name="strength",
                                                        edge_weight_units="")
        self.butadiene.add_edge(0, 1, weight=2.0)
        self.butadiene.add_edge(1, 2, weight=1.0)
        self.butadiene.add_edge(2, 3, weight=2.0)
        self.butadiene.add_edge(0, 4, weight=1.0)
        self.butadiene.add_edge(0, 5, weight=1.0)
        self.butadiene.add_edge(1, 6, weight=1.0)
        self.butadiene.add_edge(2, 7, weight=1.0)
        self.butadiene.add_edge(3, 8, weight=1.0)
        self.butadiene.add_edge(3, 9, weight=1.0)

        ethylene = Molecule.from_file(os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                                   "test_files/graphs/ethylene.xyz"))
        self.ethylene = MoleculeGraph.with_empty_graph(ethylene,
                                                       edge_weight_name="strength",
                                                       edge_weight_units="")
        self.ethylene.add_edge(0, 1, weight=2.0)
        self.ethylene.add_edge(0, 2, weight=1.0)
        self.ethylene.add_edge(0, 3, weight=1.0)
        self.ethylene.add_edge(1, 4, weight=1.0)
        self.ethylene.add_edge(1, 5, weight=1.0)

        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.resetwarnings()
        del self.ethylene
        del self.butadiene
        del self.cyclohexene

    def test_properties(self):
        self.assertEqual(self.cyclohexene.name, "bonds")
        self.assertEqual(self.cyclohexene.edge_weight_name, "strength")
        self.assertEqual(self.cyclohexene.edge_weight_unit, "")
        self.assertEqual(self.cyclohexene.get_coordination_of_site(0), 4)
        self.assertEqual(self.cyclohexene.get_coordination_of_site(2), 3)
        self.assertEqual(self.cyclohexene.get_coordination_of_site(15), 1)
        self.assertEqual(len(self.cyclohexene.get_connected_sites(0)), 4)
        self.assertTrue(isinstance(self.cyclohexene.get_connected_sites(0)[0].site, Site))
        self.assertEqual(str(self.cyclohexene.get_connected_sites(0)[0].site.specie), 'H')

    @unittest.skipIf(not nx, "NetworkX not present. Skipping...")
    def test_set_node_attributes(self):
        self.ethylene.set_node_attributes()

        specie = nx.get_node_attributes(self.ethylene.graph, "specie")
        coords = nx.get_node_attributes(self.ethylene.graph, "coords")

        self.assertEqual(str(specie[0]), str(self.ethylene.molecule[0].specie))
        self.assertEqual(str(specie[0]), "C")
        self.assertEqual(coords[0][0], self.ethylene.molecule[0].coords[0])
        self.assertEqual(coords[0][1], self.ethylene.molecule[0].coords[1])
        self.assertEqual(coords[0][2], self.ethylene.molecule[0].coords[2])

    def test_coordination(self):
        molecule = Molecule(['C', 'C'], [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])

        mg = MoleculeGraph.with_empty_graph(molecule)
        self.assertEqual(mg.get_coordination_of_site(0), 0)

        self.assertEqual(self.cyclohexene.get_coordination_of_site(0), 4)

    def test_edge_editing(self):
        self.cyclohexene.alter_edge(0, 1, new_weight=0.0, new_edge_properties={"foo": "bar"})
        new_edge = self.cyclohexene.graph.get_edge_data(0, 1)[0]
        self.assertEqual(new_edge["weight"], 0.0)
        self.assertEqual(new_edge["foo"], "bar")

        self.cyclohexene.break_edge(0, 1)
        self.assertTrue(self.cyclohexene.graph.get_edge_data(0, 1) is None)

        # Replace the now-broken edge
        self.cyclohexene.add_edge(0, 1, weight=1.0)

    def test_split(self):
        bonds = [(0, 1), (4, 5)]
        alterations = {(2, 3): {"weight": 1.0},
                       (0, 5): {"weight": 2.0},
                       (1, 2): {"weight": 2.0},
                       (3, 4): {"weight": 2.0}
                       }
        reactants = self.cyclohexene.split_molecule_subgraphs(bonds, alterations=alterations)
        self.assertTrue(isinstance(reactants, list))

        reactants = sorted(reactants, key=len)

        self.assertEqual(reactants[0], self.ethylene)
        self.assertEqual(reactants[1], self.butadiene)

    def test_find_rings(self):
        rings = self.cyclohexene.find_rings(including=[0])
        self.assertEqual(sorted(rings[0]),
                         [(0, 5), (1, 0), (2, 1), (3, 2), (4, 3), (5, 4)])
        no_rings = self.butadiene.find_rings()
        self.assertEqual(no_rings, [])

    def test_equivalent_to(self):
        ethylene = Molecule.from_file(os.path.join(os.path.dirname(__file__),
                                                   "..", "..", "..",
                                                   "test_files/graphs/ethylene.xyz"))
        # switch carbons
        ethylene[0], ethylene[1] = ethylene[1], ethylene[0]

        eth_copy = MoleculeGraph.with_empty_graph(ethylene,
                                                  edge_weight_name="strength",
                                                  edge_weight_units="")
        eth_copy.add_edge(0, 1, weight=2.0)
        eth_copy.add_edge(1, 2, weight=1.0)
        eth_copy.add_edge(1, 3, weight=1.0)
        eth_copy.add_edge(0, 4, weight=1.0)
        eth_copy.add_edge(0, 5, weight=1.0)

        self.assertTrue(self.ethylene.equivalent_to(eth_copy))
        self.assertFalse(self.ethylene.equivalent_to(self.butadiene))

    def test_substitute(self):
        molecule = FunctionalGroups["methyl"]
        molgraph = MoleculeGraph.with_empty_graph(molecule,
                                                  edge_weight_name="strength",
                                                  edge_weight_units="")
        molgraph.add_edge(0, 1, weight=1.0)
        molgraph.add_edge(0, 2, weight=1.0)
        molgraph.add_edge(0, 3, weight=1.0)

        eth_mol = copy.deepcopy(self.ethylene)
        eth_str = copy.deepcopy(self.ethylene)
        eth_mol.substitute_group(5, molecule, MinimumDistanceNN)
        eth_str.substitute_group(5, "methyl", MinimumDistanceNN)
        self.assertEqual(eth_mol, eth_str)

        graph_dict = {(0, 1): {"weight": 1.0},
                      (0, 2): {"weight": 1.0},
                      (0, 3): {"weight": 1.0},
                      }
        eth_mg = copy.deepcopy(self.ethylene)
        eth_graph = copy.deepcopy(self.ethylene)

        eth_graph.substitute_group(5, molecule, MinimumDistanceNN, graph_dict=graph_dict)
        eth_mg.substitute_group(5, molgraph, MinimumDistanceNN)
        self.assertEqual(eth_graph.graph.get_edge_data(5, 6)[0]["weight"], 1.0)
        self.assertEqual(eth_mg, eth_graph)

    def test_as_from_dict(self):
        d = self.cyclohexene.as_dict()
        mg = MoleculeGraph.from_dict(d)
        d2 = mg.as_dict()
        self.assertDictEqual(d, d2)

if __name__ == "__main__":
    unittest.main()
