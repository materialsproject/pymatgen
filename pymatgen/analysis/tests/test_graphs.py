# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import unittest
import os
import copy

from monty.serialization import loadfn  # , dumpfn

from pymatgen.command_line.critic2_caller import Critic2Analysis
from pymatgen.core.structure import Molecule, Structure, FunctionalGroups, Site
from pymatgen.analysis.graphs import *
from pymatgen.analysis.local_env import (
    MinimumDistanceNN,
    MinimumOKeeffeNN,
    OpenBabelNN,
    CutOffDictNN,
    VoronoiNN,
    CovalentBondNN
)
from pymatgen.util.testing import PymatgenTest
try:
    from openbabel import openbabel as ob
except ImportError:
    ob = None
try:
    import networkx as nx
    import networkx.algorithms.isomorphism as iso
except ImportError:
    nx = None

__author__ = "Matthew Horton, Evan Spotte-Smith"
__version__ = "0.1"
__maintainer__ = "Matthew Horton"
__email__ = "mkhorton@lbl.gov"
__status__ = "Beta"
__date__ = "August 2017"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class StructureGraphTest(PymatgenTest):
    def setUp(self):

        self.maxDiff = None

        # trivial example, simple square lattice for testing
        structure = Structure(Lattice.tetragonal(5.0, 50.0), ["H"], [[0, 0, 0]])
        self.square_sg = StructureGraph.with_empty_graph(
            structure, edge_weight_name="", edge_weight_units=""
        )
        self.square_sg.add_edge(0, 0, from_jimage=(0, 0, 0), to_jimage=(1, 0, 0))
        self.square_sg.add_edge(0, 0, from_jimage=(0, 0, 0), to_jimage=(-1, 0, 0))
        self.square_sg.add_edge(0, 0, from_jimage=(0, 0, 0), to_jimage=(0, 1, 0))
        self.square_sg.add_edge(0, 0, from_jimage=(0, 0, 0), to_jimage=(0, -1, 0))
        # TODO: decorating still fails because the structure graph gives a CN of 8 for this square lattice
        # self.square_sg.decorate_structure_with_ce_info()

        # body-centered square lattice for testing
        structure = Structure(
            Lattice.tetragonal(5.0, 50.0), ["H", "He"], [[0, 0, 0], [0.5, 0.5, 0.5]]
        )
        self.bc_square_sg = StructureGraph.with_empty_graph(
            structure, edge_weight_name="", edge_weight_units=""
        )
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
        structure = Structure(
            Lattice.tetragonal(5.0, 50.0), ["H", "He"], [[0, 0, 0], [0.5, 0.5, 0.5]]
        )
        self.bc_square_sg_r = StructureGraph.with_empty_graph(
            structure, edge_weight_name="", edge_weight_units=""
        )
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
        stdout_file = os.path.join(
            os.path.dirname(__file__),
            "..",
            "..",
            "..",
            "test_files/critic2/MoS2_critic2_stdout.txt",
        )
        with open(stdout_file, "r") as f:
            reference_stdout = f.read()
        self.structure = Structure.from_file(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "..",
                "..",
                "test_files/critic2/MoS2.cif",
            )
        )
        c2o = Critic2Analysis(self.structure, reference_stdout)
        self.mos2_sg = c2o.structure_graph(
            include_critical_points=False
        )

        latt = Lattice.cubic(4.17)
        species = ["Ni", "O"]
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        self.NiO = Structure.from_spacegroup(
            225, latt, species, coords
        ).get_primitive_structure()

        # BCC example.
        self.bcc = Structure(
            Lattice.cubic(5.0), ["He", "He"], [[0, 0, 0], [0.5, 0.5, 0.5]]
        )

        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_inappropriate_construction(self):
        # Check inappropriate strategy
        with self.assertRaises(ValueError):
            StructureGraph.with_local_env_strategy(self.NiO, CovalentBondNN())

    def test_properties(self):

        self.assertEqual(self.mos2_sg.name, "bonds")
        self.assertEqual(self.mos2_sg.edge_weight_name, "bond_length")
        self.assertEqual(self.mos2_sg.edge_weight_unit, "Å")
        self.assertEqual(self.mos2_sg.get_coordination_of_site(0), 6)
        self.assertEqual(len(self.mos2_sg.get_connected_sites(0)), 6)
        self.assertTrue(
            isinstance(self.mos2_sg.get_connected_sites(0)[0].site, PeriodicSite)
        )
        self.assertEqual(str(self.mos2_sg.get_connected_sites(0)[0].site.specie), "S")
        self.assertAlmostEqual(
            self.mos2_sg.get_connected_sites(0, jimage=(0, 0, 100))[0].site.frac_coords[
                2
            ],
            100.303027,
        )

        # these two graphs should be equivalent
        for n in range(len(self.bc_square_sg)):
            self.assertEqual(
                self.bc_square_sg.get_coordination_of_site(n),
                self.bc_square_sg_r.get_coordination_of_site(n),
            )

        # test we're not getting duplicate connected sites
        # thanks to Jack D. Sundberg for reporting this bug

        # known example where this bug occurred due to edge weights not being
        # bit-for-bit identical in otherwise identical edges
        nacl_lattice = Lattice(
            [
                [3.48543625, 0.0, 2.01231756],
                [1.16181208, 3.28610081, 2.01231756],
                [0.0, 0.0, 4.02463512],
            ]
        )
        nacl = Structure(nacl_lattice, ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])

        nacl_graph = StructureGraph.with_local_env_strategy(
            nacl, CutOffDictNN({("Cl", "Cl"): 5.0})
        )

        self.assertEqual(len(nacl_graph.get_connected_sites(1)), 12)
        self.assertEqual(len(nacl_graph.graph.get_edge_data(1, 1)), 12)

    @unittest.skipIf(not nx, "NetworkX not present. Skipping...")
    def test_set_node_attributes(self):
        self.square_sg.set_node_attributes()

        specie = nx.get_node_attributes(self.square_sg.graph, "specie")
        coords = nx.get_node_attributes(self.square_sg.graph, "coords")
        properties = nx.get_node_attributes(self.square_sg.graph, "properties")

        for i in range(len(self.square_sg.structure)):
            self.assertEqual(str(specie[i]), str(self.square_sg.structure[i].specie))
            self.assertEqual(coords[i][0], self.square_sg.structure[i].coords[0])
            self.assertEqual(coords[i][1], self.square_sg.structure[i].coords[1])
            self.assertEqual(coords[i][2], self.square_sg.structure[i].coords[2])
            self.assertEqual(properties[i], self.square_sg.structure[i].properties)

    def test_edge_editing(self):
        square = copy.deepcopy(self.square_sg)

        square.alter_edge(
            0,
            0,
            to_jimage=(1, 0, 0),
            new_weight=0.0,
            new_edge_properties={"foo": "bar"},
        )
        new_edge = square.graph.get_edge_data(0, 0)[0]
        self.assertEqual(new_edge["weight"], 0.0)
        self.assertEqual(new_edge["foo"], "bar")

        square.break_edge(0, 0, to_jimage=(1, 0, 0))
        self.assertEqual(len(square.graph.get_edge_data(0, 0)), 3)

    def test_insert_remove(self):
        struct_copy = copy.deepcopy(self.square_sg.structure)
        square_copy = copy.deepcopy(self.square_sg)

        # Ensure that insert_node appropriately wraps Structure.insert()
        struct_copy.insert(1, "O", [0.5, 0.5, 0.5])
        square_copy.insert_node(1, "O", [0.5, 0.5, 0.5])
        self.assertEqual(struct_copy, square_copy.structure)

        # Test that removal is also equivalent between Structure and StructureGraph.structure
        struct_copy.remove_sites([1])
        square_copy.remove_nodes([1])
        self.assertEqual(struct_copy, square_copy.structure)

        square_copy.insert_node(
            1,
            "O",
            [0.5, 0.5, 0.5],
            edges=[{"from_index": 1, "to_index": 0, "to_jimage": (0, 0, 0)}],
        )
        self.assertEqual(square_copy.get_coordination_of_site(1), 1)

        # Test that StructureGraph.graph is correctly updated
        square_copy.insert_node(1, "H", [0.5, 0.5, 0.75], edges=[{"from_index": 1,
                                                                  "to_index": 2,
                                                                  "to_jimage": (0, 0, 0)}])
        square_copy.remove_nodes([1])
        self.assertEqual(square_copy.graph.number_of_nodes(), 2)
        self.assertEqual(square_copy.graph.number_of_edges(), 5)

    def test_substitute(self):
        structure = Structure.from_file(
            os.path.join(
                os.path.dirname(__file__), "..", "..", "..", "test_files", "Li2O.cif"
            )
        )
        molecule = FunctionalGroups["methyl"]

        structure_copy = copy.deepcopy(structure)
        structure_copy_graph = copy.deepcopy(structure)

        sg = StructureGraph.with_local_env_strategy(structure, MinimumDistanceNN())
        sg_copy = copy.deepcopy(sg)

        # Ensure that strings and molecules lead to equivalent substitutions
        sg.substitute_group(1, molecule, MinimumDistanceNN)
        sg_copy.substitute_group(1, "methyl", MinimumDistanceNN)
        self.assertEqual(sg, sg_copy)

        # Ensure that the underlying structure has been modified as expected
        structure_copy.substitute(1, "methyl")
        self.assertEqual(structure_copy, sg.structure)

        # Test inclusion of graph dictionary
        graph_dict = {
            (0, 1): {"weight": 0.5},
            (0, 2): {"weight": 0.5},
            (0, 3): {"weight": 0.5},
        }

        sg_with_graph = StructureGraph.with_local_env_strategy(
            structure_copy_graph, MinimumDistanceNN()
        )
        sg_with_graph.substitute_group(
            1, "methyl", MinimumDistanceNN, graph_dict=graph_dict
        )
        edge = sg_with_graph.graph.get_edge_data(11, 13)[0]
        self.assertEqual(edge["weight"], 0.5)

    def test_auto_image_detection(self):

        sg = StructureGraph.with_empty_graph(self.structure)
        sg.add_edge(0, 0)

        ref_edges = [
            (0, 0, {"to_jimage": (-1, -1, 0)}),
            (0, 0, {"to_jimage": (-1, 0, 0)}),
            (0, 0, {"to_jimage": (0, -1, 0)}),
            (0, 0, {"to_jimage": (0, 1, 0)}),
            (0, 0, {"to_jimage": (1, 0, 0)}),
        ]
        self.assertEqual(len(list(sg.graph.edges(data=True))), 6)

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
        self.mos2_sg.graph.graph["edge_weight_units"] = "A"
        self.assertStrContentEqual(str(self.square_sg), square_sg_str_ref)
        self.assertStrContentEqual(str(self.mos2_sg), mos2_sg_str_ref)

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

        self.assertStrContentEqual(square_sg_mul_actual_str, square_sg_mul_ref_str)

        # test sequential multiplication
        sq_sg_1 = self.square_sg * (2, 2, 1)
        sq_sg_1 = sq_sg_1 * (2, 2, 1)
        sq_sg_2 = self.square_sg * (4, 4, 1)
        self.assertEqual(
            sq_sg_1.graph.number_of_edges(), sq_sg_2.graph.number_of_edges()
        )
        # TODO: the below test still gives 8 != 4
        # self.assertEqual(self.square_sg.get_coordination_of_site(0), 4)

        mos2_sg_mul = self.mos2_sg * (3, 3, 1)
        for idx in mos2_sg_mul.structure.indices_from_symbol("Mo"):
            self.assertEqual(mos2_sg_mul.get_coordination_of_site(idx), 6)

        mos2_sg_premul = StructureGraph.with_local_env_strategy(
            self.structure * (3, 3, 1), MinimumDistanceNN()
        )
        self.assertTrue(mos2_sg_mul == mos2_sg_premul)

        # test 3D Structure

        nio_sg = StructureGraph.with_local_env_strategy(self.NiO, MinimumDistanceNN())
        nio_sg = nio_sg * 3

        for n in range(len(nio_sg)):
            self.assertEqual(nio_sg.get_coordination_of_site(n), 6)

    @unittest.skipIf(
        not (which("neato") and which("fdp")), "graphviz executables not present"
    )
    def test_draw(self):

        # draw MoS2 graph
        self.mos2_sg.draw_graph_to_file(
            "MoS2_single.pdf", image_labels=True, hide_image_edges=False
        )
        mos2_sg = self.mos2_sg * (9, 9, 1)
        mos2_sg.draw_graph_to_file("MoS2.pdf", algo="neato")

        # draw MoS2 graph that's been successively multiplied
        mos2_sg_2 = self.mos2_sg * (3, 3, 1)
        mos2_sg_2 = mos2_sg_2 * (3, 3, 1)
        mos2_sg_2.draw_graph_to_file(
            "MoS2_twice_mul.pdf", algo="neato", hide_image_edges=True
        )

        # draw MoS2 graph that's generated from a pre-multiplied Structure
        mos2_sg_premul = StructureGraph.with_local_env_strategy(
            self.structure * (3, 3, 1), MinimumDistanceNN()
        )
        mos2_sg_premul.draw_graph_to_file(
            "MoS2_premul.pdf", algo="neato", hide_image_edges=True
        )

        # draw graph for a square lattice
        self.square_sg.draw_graph_to_file("square_single.pdf", hide_image_edges=False)
        square_sg = self.square_sg * (5, 5, 1)
        square_sg.draw_graph_to_file(
            "square.pdf", algo="neato", image_labels=True, node_labels=False
        )

        # draw graph for a body-centered square lattice
        self.bc_square_sg.draw_graph_to_file(
            "bc_square_single.pdf", hide_image_edges=False
        )
        bc_square_sg = self.bc_square_sg * (9, 9, 1)
        bc_square_sg.draw_graph_to_file(
            "bc_square.pdf", algo="neato", image_labels=False
        )

        # draw graph for a body-centered square lattice defined in an alternative way
        self.bc_square_sg_r.draw_graph_to_file(
            "bc_square_r_single.pdf", hide_image_edges=False
        )
        bc_square_sg_r = self.bc_square_sg_r * (9, 9, 1)
        bc_square_sg_r.draw_graph_to_file(
            "bc_square_r.pdf", algo="neato", image_labels=False
        )

        # delete generated test files
        test_files = (
            "bc_square_r_single.pdf",
            "bc_square_r.pdf",
            "bc_square_single.pdf",
            "bc_square.pdf",
            "MoS2_premul.pdf",
            "MOS2_single.pdf",
            "MoS2_twice_mul.pdf",
            "MoS2.pdf",
            "square_single.pdf",
            "square.pdf",
        )
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
        self.assertEqual(diff["dist"], 0)

        self.assertEqual(self.square_sg.get_coordination_of_site(0), 4)

    def test_from_edges(self):
        edges = {
            (0, 0, (0, 0, 0), (1, 0, 0)): None,
            (0, 0, (0, 0, 0), (-1, 0, 0)): None,
            (0, 0, (0, 0, 0), (0, 1, 0)): None,
            (0, 0, (0, 0, 0), (0, -1, 0)): None,
        }

        structure = Structure(Lattice.tetragonal(5.0, 50.0), ["H"], [[0, 0, 0]])

        sg = StructureGraph.with_edges(structure, edges)

        self.assertEqual(sg, self.square_sg)

    def test_extract_molecules(self):

        structure_file = os.path.join(
            os.path.dirname(__file__),
            "..",
            "..",
            "..",
            "test_files/H6PbCI3N_mp-977013_symmetrized.cif",
        )

        s = Structure.from_file(structure_file)

        nn = MinimumDistanceNN()
        sg = StructureGraph.with_local_env_strategy(s, nn)

        molecules = sg.get_subgraphs_as_molecules()
        self.assertEqual(molecules[0].composition.formula, "H3 C1")
        self.assertEqual(len(molecules), 1)

        molecules = self.mos2_sg.get_subgraphs_as_molecules()
        self.assertEqual(len(molecules), 0)

    def test_types_and_weights_of_connections(self):

        types = self.mos2_sg.types_and_weights_of_connections

        self.assertEqual(len(types["Mo-S"]), 6)
        self.assertAlmostEqual(types["Mo-S"][0], 2.416931678417331)

    def test_weight_statistics(self):

        weight_statistics = self.mos2_sg.weight_statistics

        self.assertEqual(len(weight_statistics["all_weights"]), 6)
        self.assertAlmostEqual(weight_statistics["min"], 2.4169314100201875)
        self.assertAlmostEqual(weight_statistics["variance"], 0)

    def test_types_of_coordination_environments(self):

        types = self.mos2_sg.types_of_coordination_environments()
        self.assertListEqual(types, ["Mo-S(6)", "S-Mo(3)"])

        types_anonymous = self.mos2_sg.types_of_coordination_environments(
            anonymous=True
        )
        self.assertListEqual(types_anonymous, ["A-B(3)", "A-B(6)"])


class MoleculeGraphTest(unittest.TestCase):
    def setUp(self):

        cyclohexene = Molecule.from_file(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "..",
                "..",
                "test_files/graphs/cyclohexene.xyz",
            )
        )
        self.cyclohexene = MoleculeGraph.with_empty_graph(
            cyclohexene, edge_weight_name="strength", edge_weight_units=""
        )
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

        butadiene = Molecule.from_file(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "..",
                "..",
                "test_files/graphs/butadiene.xyz",
            )
        )
        self.butadiene = MoleculeGraph.with_empty_graph(
            butadiene, edge_weight_name="strength", edge_weight_units=""
        )
        self.butadiene.add_edge(0, 1, weight=2.0)
        self.butadiene.add_edge(1, 2, weight=1.0)
        self.butadiene.add_edge(2, 3, weight=2.0)
        self.butadiene.add_edge(0, 4, weight=1.0)
        self.butadiene.add_edge(0, 5, weight=1.0)
        self.butadiene.add_edge(1, 6, weight=1.0)
        self.butadiene.add_edge(2, 7, weight=1.0)
        self.butadiene.add_edge(3, 8, weight=1.0)
        self.butadiene.add_edge(3, 9, weight=1.0)

        ethylene = Molecule.from_file(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "..",
                "..",
                "test_files/graphs/ethylene.xyz",
            )
        )
        self.ethylene = MoleculeGraph.with_empty_graph(
            ethylene, edge_weight_name="strength", edge_weight_units=""
        )
        self.ethylene.add_edge(0, 1, weight=2.0)
        self.ethylene.add_edge(0, 2, weight=1.0)
        self.ethylene.add_edge(0, 3, weight=1.0)
        self.ethylene.add_edge(1, 4, weight=1.0)
        self.ethylene.add_edge(1, 5, weight=1.0)

        self.pc = Molecule.from_file(
            os.path.join(module_dir, "..", "..", "..", "test_files", "graphs", "PC.xyz")
        )
        self.pc_edges = [
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
        self.pc_frag1 = Molecule.from_file(
            os.path.join(
                module_dir, "..", "..", "..", "test_files", "graphs", "PC_frag1.xyz"
            )
        )
        self.pc_frag1_edges = [[0, 2], [4, 2], [2, 1], [1, 3]]
        self.tfsi = Molecule.from_file(
            os.path.join(
                module_dir, "..", "..", "..", "test_files", "graphs", "TFSI.xyz"
            )
        )
        self.tfsi_edges = (
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

        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")
        del self.ethylene
        del self.butadiene
        del self.cyclohexene

    @unittest.skipIf(not ob, "OpenBabel not present. Skipping...")
    def test_construction(self):
        edges_frag = {(e[0], e[1]): {"weight": 1.0} for e in self.pc_frag1_edges}
        mol_graph = MoleculeGraph.with_edges(self.pc_frag1, edges_frag)
        # dumpfn(mol_graph.as_dict(), os.path.join(module_dir,"pc_frag1_mg.json"))
        ref_mol_graph = loadfn(os.path.join(module_dir, "pc_frag1_mg.json"))
        self.assertEqual(mol_graph, ref_mol_graph)
        self.assertEqual(mol_graph.graph.adj, ref_mol_graph.graph.adj)
        for node in mol_graph.graph.nodes:
            self.assertEqual(
                mol_graph.graph.nodes[node]["specie"],
                ref_mol_graph.graph.nodes[node]["specie"],
            )
            for ii in range(3):
                self.assertEqual(
                    mol_graph.graph.nodes[node]["coords"][ii],
                    ref_mol_graph.graph.nodes[node]["coords"][ii],
                )

        edges_pc = {(e[0], e[1]): {"weight": 1.0} for e in self.pc_edges}
        mol_graph = MoleculeGraph.with_edges(self.pc, edges_pc)
        # dumpfn(mol_graph.as_dict(), os.path.join(module_dir,"pc_mg.json"))
        ref_mol_graph = loadfn(os.path.join(module_dir, "pc_mg.json"))
        self.assertEqual(mol_graph, ref_mol_graph)
        self.assertEqual(mol_graph.graph.adj, ref_mol_graph.graph.adj)
        for node in mol_graph.graph:
            self.assertEqual(
                mol_graph.graph.nodes[node]["specie"],
                ref_mol_graph.graph.nodes[node]["specie"],
            )
            for ii in range(3):
                self.assertEqual(
                    mol_graph.graph.nodes[node]["coords"][ii],
                    ref_mol_graph.graph.nodes[node]["coords"][ii],
                )

        mol_graph_edges = MoleculeGraph.with_edges(self.pc, edges=edges_pc)
        mol_graph_strat = MoleculeGraph.with_local_env_strategy(self.pc, OpenBabelNN())

        self.assertTrue(mol_graph_edges.isomorphic_to(mol_graph_strat))

        # Check inappropriate strategy
        with self.assertRaises(ValueError):
            MoleculeGraph.with_local_env_strategy(self.pc, VoronoiNN())

    def test_properties(self):
        self.assertEqual(self.cyclohexene.name, "bonds")
        self.assertEqual(self.cyclohexene.edge_weight_name, "strength")
        self.assertEqual(self.cyclohexene.edge_weight_unit, "")
        self.assertEqual(self.cyclohexene.get_coordination_of_site(0), 4)
        self.assertEqual(self.cyclohexene.get_coordination_of_site(2), 3)
        self.assertEqual(self.cyclohexene.get_coordination_of_site(15), 1)
        self.assertEqual(len(self.cyclohexene.get_connected_sites(0)), 4)
        self.assertTrue(
            isinstance(self.cyclohexene.get_connected_sites(0)[0].site, Site)
        )
        self.assertEqual(
            str(self.cyclohexene.get_connected_sites(0)[0].site.specie), "H"
        )

    @unittest.skipIf(not nx, "NetworkX not present. Skipping...")
    def test_set_node_attributes(self):
        self.ethylene.set_node_attributes()

        specie = nx.get_node_attributes(self.ethylene.graph, "specie")
        coords = nx.get_node_attributes(self.ethylene.graph, "coords")
        properties = nx.get_node_attributes(self.ethylene.graph, "properties")

        for i in range(len(self.ethylene.molecule)):
            self.assertEqual(str(specie[i]), str(self.ethylene.molecule[i].specie))
            self.assertEqual(coords[i][0], self.ethylene.molecule[i].coords[0])
            self.assertEqual(coords[i][1], self.ethylene.molecule[i].coords[1])
            self.assertEqual(coords[i][2], self.ethylene.molecule[i].coords[2])
            self.assertEqual(properties[i], self.ethylene.molecule[i].properties)

    def test_coordination(self):
        molecule = Molecule(["C", "C"], [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])

        mg = MoleculeGraph.with_empty_graph(molecule)
        self.assertEqual(mg.get_coordination_of_site(0), 0)

        self.assertEqual(self.cyclohexene.get_coordination_of_site(0), 4)

    def test_edge_editing(self):
        self.cyclohexene.alter_edge(
            0, 1, new_weight=0.0, new_edge_properties={"foo": "bar"}
        )
        new_edge = self.cyclohexene.graph.get_edge_data(0, 1)[0]
        self.assertEqual(new_edge["weight"], 0.0)
        self.assertEqual(new_edge["foo"], "bar")

        self.cyclohexene.break_edge(0, 1)
        self.assertTrue(self.cyclohexene.graph.get_edge_data(0, 1) is None)

        # Replace the now-broken edge
        self.cyclohexene.add_edge(0, 1, weight=1.0)

    def test_insert_remove(self):
        mol_copy = copy.deepcopy(self.ethylene.molecule)
        eth_copy = copy.deepcopy(self.ethylene)

        # Ensure that insert_node appropriately wraps Molecule.insert()
        mol_copy.insert(1, "O", [0.5, 0.5, 0.5])
        eth_copy.insert_node(1, "O", [0.5, 0.5, 0.5])
        self.assertEqual(mol_copy, eth_copy.molecule)

        # Test that removal is also equivalent between Molecule and MoleculeGraph.molecule
        mol_copy.remove_sites([1])
        eth_copy.remove_nodes([1])
        self.assertEqual(mol_copy, eth_copy.molecule)

        eth_copy.insert_node(
            1,
            "O",
            [0.5, 0.5, 0.5],
            edges=[{"from_index": 1, "to_index": 2}, {"from_index": 1, "to_index": 3}],
        )
        self.assertEqual(eth_copy.get_coordination_of_site(1), 2)

        # Test that MoleculeGraph.graph is correctly updated
        eth_copy.remove_nodes([1, 2])
        self.assertEqual(eth_copy.graph.number_of_nodes(), 5)
        self.assertEqual(eth_copy.graph.number_of_edges(), 2)

    @unittest.skipIf(not nx, "NetworkX not present. Skipping...")
    def test_split(self):
        bonds = [(0, 1), (4, 5)]
        alterations = {
            (2, 3): {"weight": 1.0},
            (0, 5): {"weight": 2.0},
            (1, 2): {"weight": 2.0},
            (3, 4): {"weight": 2.0},
        }
        # Perform retro-Diels-Alder reaction - turn product into reactants
        reactants = self.cyclohexene.split_molecule_subgraphs(
            bonds, allow_reverse=True, alterations=alterations
        )
        self.assertTrue(isinstance(reactants, list))

        reactants = sorted(reactants, key=len)
        # After alterations, reactants sholuld be ethylene and butadiene
        self.assertEqual(reactants[0], self.ethylene)
        self.assertEqual(reactants[1], self.butadiene)

        with self.assertRaises(MolGraphSplitError):
            self.cyclohexene.split_molecule_subgraphs([(0, 1)])

        # Test naive charge redistribution
        hydroxide = Molecule(["O", "H"], [[0, 0, 0], [0.5, 0.5, 0.5]], charge=-1)
        oh_mg = MoleculeGraph.with_empty_graph(hydroxide)

        oh_mg.add_edge(0, 1)

        new_mgs = oh_mg.split_molecule_subgraphs([(0, 1)])
        for mg in new_mgs:
            if str(mg.molecule[0].specie) == "O":
                self.assertEqual(mg.molecule.charge, -1)
            else:
                self.assertEqual(mg.molecule.charge, 0)

        # Trying to test to ensure that remapping of nodes to atoms works
        diff_species = Molecule(
            ["C", "I", "Cl", "Br", "F"],
            [
                [0.8314, -0.2682, -0.9102],
                [1.3076, 1.3425, -2.2038],
                [-0.8429, -0.7410, -1.1554],
                [1.9841, -1.7636, -1.2953],
                [1.0098, 0.1231, 0.3916],
            ],
        )

        diff_spec_mg = MoleculeGraph.with_empty_graph(diff_species)
        diff_spec_mg.add_edge(0, 1)
        diff_spec_mg.add_edge(0, 2)
        diff_spec_mg.add_edge(0, 3)
        diff_spec_mg.add_edge(0, 4)

        for i in range(1, 5):
            bond = (0, i)

            split_mgs = diff_spec_mg.split_molecule_subgraphs([bond])
            for split_mg in split_mgs:
                species = nx.get_node_attributes(split_mg.graph, "specie")

                for j in range(len(split_mg.graph.nodes)):
                    atom = split_mg.molecule[j]
                    self.assertEqual(species[j], str(atom.specie))

    @unittest.skipIf(not nx, "NetworkX not present. Skipping...")
    def test_build_unique_fragments(self):
        edges = {(e[0], e[1]): None for e in self.pc_edges}
        mol_graph = MoleculeGraph.with_edges(self.pc, edges)
        unique_fragment_dict = mol_graph.build_unique_fragments()
        unique_fragments = []
        for key in unique_fragment_dict:
            for fragment in unique_fragment_dict[key]:
                unique_fragments.append(fragment)
        self.assertEqual(len(unique_fragments), 295)
        nm = iso.categorical_node_match("specie", "ERROR")
        for ii in range(295):
            # Test that each fragment is unique
            for jj in range(ii + 1, 295):
                self.assertFalse(
                    nx.is_isomorphic(
                        unique_fragments[ii].graph,
                        unique_fragments[jj].graph,
                        node_match=nm,
                    )
                )

            # Test that each fragment correctly maps between Molecule and graph
            self.assertEqual(
                len(unique_fragments[ii].molecule),
                len(unique_fragments[ii].graph.nodes),
            )
            species = nx.get_node_attributes(unique_fragments[ii].graph, "specie")
            coords = nx.get_node_attributes(unique_fragments[ii].graph, "coords")

            mol = unique_fragments[ii].molecule
            for ss, site in enumerate(mol):
                self.assertEqual(str(species[ss]), str(site.specie))
                self.assertEqual(coords[ss][0], site.coords[0])
                self.assertEqual(coords[ss][1], site.coords[1])
                self.assertEqual(coords[ss][2], site.coords[2])

            # Test that each fragment is connected
            self.assertTrue(nx.is_connected(unique_fragments[ii].graph.to_undirected()))

    def test_find_rings(self):
        rings = self.cyclohexene.find_rings(including=[0])
        self.assertEqual(
            sorted(rings[0]), [(0, 5), (1, 0), (2, 1), (3, 2), (4, 3), (5, 4)]
        )
        no_rings = self.butadiene.find_rings()
        self.assertEqual(no_rings, [])

    def test_isomorphic(self):
        ethylene = Molecule.from_file(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "..",
                "..",
                "test_files/graphs/ethylene.xyz",
            )
        )
        # switch carbons
        ethylene[0], ethylene[1] = ethylene[1], ethylene[0]

        eth_copy = MoleculeGraph.with_edges(
            ethylene,
            {
                (0, 1): {"weight": 2},
                (1, 2): {"weight": 1},
                (1, 3): {"weight": 1},
                (0, 4): {"weight": 1},
                (0, 5): {"weight": 1},
            },
        )
        # If they are equal, they must also be isomorphic
        eth_copy = copy.deepcopy(self.ethylene)
        self.assertTrue(self.ethylene.isomorphic_to(eth_copy))
        self.assertFalse(self.butadiene.isomorphic_to(self.ethylene))

    def test_substitute(self):
        molecule = FunctionalGroups["methyl"]
        molgraph = MoleculeGraph.with_edges(
            molecule,
            {(0, 1): {"weight": 1}, (0, 2): {"weight": 1}, (0, 3): {"weight": 1}},
        )

        eth_mol = copy.deepcopy(self.ethylene)
        eth_str = copy.deepcopy(self.ethylene)
        # Ensure that strings and molecules lead to equivalent substitutions
        eth_mol.substitute_group(5, molecule, MinimumDistanceNN)
        eth_str.substitute_group(5, "methyl", MinimumDistanceNN)
        self.assertEqual(eth_mol, eth_str)

        graph_dict = {
            (0, 1): {"weight": 1.0},
            (0, 2): {"weight": 1.0},
            (0, 3): {"weight": 1.0},
        }
        eth_mg = copy.deepcopy(self.ethylene)
        eth_graph = copy.deepcopy(self.ethylene)

        # Check that MoleculeGraph input is handled properly
        eth_graph.substitute_group(
            5, molecule, MinimumDistanceNN, graph_dict=graph_dict
        )
        eth_mg.substitute_group(5, molgraph, MinimumDistanceNN)
        self.assertEqual(eth_graph.graph.get_edge_data(5, 6)[0]["weight"], 1.0)
        self.assertEqual(eth_mg, eth_graph)

    def test_replace(self):
        eth_copy_sub = copy.deepcopy(self.ethylene)
        eth_copy_repl = copy.deepcopy(self.ethylene)
        # First, perform a substiution as above
        eth_copy_sub.substitute_group(5, "methyl", MinimumDistanceNN)
        eth_copy_repl.replace_group(5, "methyl", MinimumDistanceNN)
        # Test that replacement on a terminal atom is equivalent to substitution
        self.assertEqual(eth_copy_repl.molecule, eth_copy_sub.molecule)
        self.assertEqual(eth_copy_repl, eth_copy_sub)

        # Methyl carbon should have coordination 4
        self.assertEqual(eth_copy_repl.get_coordination_of_site(5), 4)
        # Now swap one functional group for another
        eth_copy_repl.replace_group(5, "amine", MinimumDistanceNN)
        self.assertEqual(
            ["C", "C", "H", "H", "H", "N", "H", "H"],
            [str(s) for s in eth_copy_repl.molecule.species],
        )
        self.assertEqual(len(eth_copy_repl.graph.nodes), 8)
        # Amine nitrogen should have coordination 3
        self.assertEqual(eth_copy_repl.get_coordination_of_site(5), 3)

    def test_as_from_dict(self):
        d = self.cyclohexene.as_dict()
        mg = MoleculeGraph.from_dict(d)
        d2 = mg.as_dict()
        self.assertEqual(str(d), str(d2))


if __name__ == "__main__":
    unittest.main()
