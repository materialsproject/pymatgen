__author__ = "waroquiers"

import numpy as np

from pymatgen.analysis.chemenv.connectivity.environment_nodes import EnvironmentNode
from pymatgen.analysis.chemenv.utils.graph_utils import (
    MultiGraphCycle,
    SimpleGraphCycle,
    get_delta,
)
from pymatgen.util.testing import PymatgenTest


class FakeNode:
    def __init__(self, isite):
        self.isite = isite


class FakeNodeWithEqMethod:
    def __init__(self, isite):
        self.isite = isite

    def __eq__(self, other):
        return self.isite == other.isite

    def __hash__(self):
        return 0


class FakeNodeWithEqLtMethods:
    def __init__(self, isite):
        self.isite = isite

    def __eq__(self, other):
        return self.isite == other.isite

    def __lt__(self, other):
        return self.isite < other.isite

    def __str__(self):
        return f"FakeNode_{self.isite:d}"

    def __hash__(self):
        return 0


class FakeNodeWithEqLtMethodsBis(FakeNodeWithEqLtMethods):
    pass


class FakeNodeWithEqMethodWrongSortable:
    def __init__(self, isite):
        self.isite = isite

    def __eq__(self, other):
        return self.isite == other.isite

    def __hash__(self):
        return 0

    def __lt__(self, other):
        return self.isite % 2 < other.isite % 2


class GraphUtilsTest(PymatgenTest):
    def test_get_delta(self):
        n1 = FakeNode(3)
        n2 = FakeNode(7)
        edge_data = {"start": 3, "end": 7, "delta": [2, 6, 4]}
        self.assertTrue(np.allclose(get_delta(n1, n2, edge_data), [2, 6, 4]))
        edge_data = {"start": 7, "end": 3, "delta": [2, 6, 4]}
        self.assertTrue(np.allclose(get_delta(n1, n2, edge_data), [-2, -6, -4]))
        with self.assertRaisesRegex(
            ValueError,
            "Trying to find a delta between two nodes with an edge that seems not to link these nodes.",
        ):
            edge_data = {"start": 6, "end": 3, "delta": [2, 6, 4]}
            get_delta(n1, n2, edge_data)
        with self.assertRaisesRegex(
            ValueError,
            "Trying to find a delta between two nodes with an edge that seems not to link these nodes.",
        ):
            edge_data = {"start": 7, "end": 2, "delta": [2, 6, 4]}
            get_delta(n1, n2, edge_data)

    def test_simple_graph_cycle(self):
        sg_cycle1 = SimpleGraphCycle([0, 1, 2, 3])

        # Test equality
        sg_cycle2 = SimpleGraphCycle([1, 2, 3, 0])
        self.assertEqual(sg_cycle1, sg_cycle2)
        sg_cycle2 = SimpleGraphCycle([2, 3, 0, 1])
        self.assertEqual(sg_cycle1, sg_cycle2)
        sg_cycle2 = SimpleGraphCycle([3, 0, 1, 2])
        self.assertEqual(sg_cycle1, sg_cycle2)

        # Test reversed cycles
        sg_cycle2 = SimpleGraphCycle([0, 3, 2, 1])
        self.assertEqual(sg_cycle1, sg_cycle2)
        sg_cycle2 = SimpleGraphCycle([3, 2, 1, 0])
        self.assertEqual(sg_cycle1, sg_cycle2)
        sg_cycle2 = SimpleGraphCycle([2, 1, 0, 3])
        self.assertEqual(sg_cycle1, sg_cycle2)
        sg_cycle2 = SimpleGraphCycle([1, 0, 3, 2])
        self.assertEqual(sg_cycle1, sg_cycle2)

        # Test different cycle lengths inequality
        sg_cycle2 = SimpleGraphCycle([0, 1, 2, 3, 4])
        self.assertNotEqual(sg_cycle1, sg_cycle2)

        # Test equality of self-loops
        self.assertEqual(SimpleGraphCycle([0]), SimpleGraphCycle([0]))
        self.assertNotEqual(SimpleGraphCycle([0]), SimpleGraphCycle([4]))

        # Test inequality inversing two nodes
        sg_cycle2 = SimpleGraphCycle([0, 1, 3, 2])
        self.assertNotEqual(sg_cycle1, sg_cycle2)

        # Test inequality with different nodes
        sg_cycle2 = SimpleGraphCycle([4, 1, 2, 3])
        self.assertNotEqual(sg_cycle1, sg_cycle2)

        # Test hashing function
        self.assertEqual(hash(sg_cycle1), 4)
        self.assertEqual(hash(SimpleGraphCycle([0])), 1)
        self.assertEqual(hash(SimpleGraphCycle([0, 1, 3, 6])), 4)
        self.assertEqual(hash(SimpleGraphCycle([0, 1, 2])), 3)

        # Test from_edges function
        #   3-nodes cycle
        edges = [(0, 2), (4, 2), (0, 4)]
        sg_cycle = SimpleGraphCycle.from_edges(edges=edges, edges_are_ordered=False)
        self.assertEqual(sg_cycle, SimpleGraphCycle([4, 0, 2]))

        #   Self-loop cycle
        edges = [(2, 2)]
        sg_cycle = SimpleGraphCycle.from_edges(edges=edges)
        self.assertEqual(sg_cycle, SimpleGraphCycle([2]))

        #   5-nodes cycle
        edges = [(0, 2), (4, 7), (2, 7), (4, 5), (5, 0)]
        sg_cycle = SimpleGraphCycle.from_edges(edges=edges, edges_are_ordered=False)
        self.assertEqual(sg_cycle, SimpleGraphCycle([2, 7, 4, 5, 0]))

        #   two identical 3-nodes cycles
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="SimpleGraphCycle is not valid : Duplicate nodes.",
        ):
            edges = [(0, 2), (4, 2), (0, 4), (0, 2), (4, 2), (0, 4)]
            SimpleGraphCycle.from_edges(edges=edges, edges_are_ordered=False)

        #   two cycles in from_edges
        with self.assertRaisesRegex(ValueError, expected_regex="Could not construct a cycle from edges."):
            edges = [(0, 2), (4, 2), (0, 4), (1, 3), (6, 7), (3, 6), (1, 7)]
            SimpleGraphCycle.from_edges(edges=edges, edges_are_ordered=False)

        with self.assertRaisesRegex(ValueError, expected_regex="Could not construct a cycle from edges."):
            edges = [(0, 2), (4, 6), (2, 7), (4, 5), (5, 0)]
            SimpleGraphCycle.from_edges(edges=edges, edges_are_ordered=False)

        with self.assertRaisesRegex(ValueError, expected_regex="Could not construct a cycle from edges."):
            edges = [(0, 2), (4, 7), (2, 7), (4, 10), (5, 0)]
            SimpleGraphCycle.from_edges(edges=edges, edges_are_ordered=False)

        # Test as_dict from_dict and len method
        sg_cycle = SimpleGraphCycle([0, 1, 2, 3])
        self.assertEqual(sg_cycle, SimpleGraphCycle.from_dict(sg_cycle.as_dict()))
        self.assertEqual(len(sg_cycle), 4)
        sg_cycle = SimpleGraphCycle([4])
        self.assertEqual(sg_cycle, SimpleGraphCycle.from_dict(sg_cycle.as_dict()))
        self.assertEqual(len(sg_cycle), 1)
        sg_cycle = SimpleGraphCycle([4, 2, 6, 7, 9, 3, 15])
        self.assertEqual(sg_cycle, SimpleGraphCycle.from_dict(sg_cycle.as_dict()))
        self.assertEqual(len(sg_cycle), 7)

        # Check validation at instance creation time
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="SimpleGraphCycle is not valid : Duplicate nodes.",
        ):
            SimpleGraphCycle([0, 2, 4, 6, 2])

        # Check the validate method
        # Nodes not sortable
        sgc = SimpleGraphCycle(
            [FakeNodeWithEqMethod(1), FakeNodeWithEqMethod(0), FakeNodeWithEqMethod(2)],
            validate=False,
            ordered=False,
        )
        self.assertFalse(sgc.ordered)
        self.assertEqual(
            sgc.nodes,
            (FakeNodeWithEqMethod(1), FakeNodeWithEqMethod(0), FakeNodeWithEqMethod(2)),
        )
        sgc.validate(check_strict_ordering=False)
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="SimpleGraphCycle is not valid : The nodes are not sortable.",
        ):
            sgc.validate(check_strict_ordering=True)

        # Empty cycle not valid
        sgc = SimpleGraphCycle([], validate=False, ordered=False)
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="SimpleGraphCycle is not valid : Empty cycle is not valid.",
        ):
            sgc.validate()

        # Simple graph cycle with 2 nodes not valid
        sgc = SimpleGraphCycle([1, 2], validate=False, ordered=False)
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="SimpleGraphCycle is not valid : Simple graph cycle with 2 nodes is not valid.",
        ):
            sgc.validate()

        # Simple graph cycle with nodes that cannot be strictly ordered
        sgc = SimpleGraphCycle(
            [
                FakeNodeWithEqMethodWrongSortable(0),
                FakeNodeWithEqMethodWrongSortable(1),
                FakeNodeWithEqMethodWrongSortable(2),
                FakeNodeWithEqMethodWrongSortable(3),
            ],
            validate=False,
            ordered=False,
        )
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="SimpleGraphCycle is not valid : "
            "The list of nodes in the cycle cannot be strictly ordered.",
        ):
            sgc.validate(check_strict_ordering=True)

        # Check the order method
        sgc = SimpleGraphCycle(
            [
                FakeNodeWithEqMethodWrongSortable(0),
                FakeNodeWithEqMethodWrongSortable(1),
                FakeNodeWithEqMethodWrongSortable(2),
                FakeNodeWithEqMethodWrongSortable(3),
            ],
            validate=False,
            ordered=False,
        )
        sgc.order(raise_on_fail=False)
        self.assertFalse(sgc.ordered)
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="SimpleGraphCycle is not valid : "
            "The list of nodes in the cycle cannot be strictly ordered.",
        ):
            sgc.order(raise_on_fail=True)

        sgc = SimpleGraphCycle(
            [
                FakeNodeWithEqMethod(8),
                FakeNodeWithEqMethod(0),
                FakeNodeWithEqMethod(3),
                FakeNodeWithEqMethod(6),
            ],
            validate=False,
            ordered=False,
        )
        sgc.order(raise_on_fail=False)
        self.assertFalse(sgc.ordered)
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="SimpleGraphCycle is not valid : The nodes are not sortable.",
        ):
            sgc.order(raise_on_fail=True)

        sgc = SimpleGraphCycle(
            [
                FakeNodeWithEqLtMethods(8),
                FakeNodeWithEqLtMethods(0),
                FakeNodeWithEqLtMethods(6),
                FakeNodeWithEqLtMethods(3),
            ],
            validate=False,
            ordered=False,
        )
        sgc.order(raise_on_fail=True)
        self.assertTrue(sgc.ordered)
        self.assertEqual(
            sgc.nodes,
            (
                FakeNodeWithEqLtMethods(0),
                FakeNodeWithEqLtMethods(6),
                FakeNodeWithEqLtMethods(3),
                FakeNodeWithEqLtMethods(8),
            ),
        )

        sgc = SimpleGraphCycle(
            [
                FakeNodeWithEqLtMethodsBis(8),
                FakeNodeWithEqLtMethods(0),
                FakeNodeWithEqLtMethods(6),
                FakeNodeWithEqLtMethods(3),
            ],
            validate=False,
            ordered=False,
        )
        sgc.order(raise_on_fail=False)
        self.assertFalse(sgc.ordered)
        self.assertEqual(
            sgc.nodes,
            (
                FakeNodeWithEqLtMethodsBis(8),
                FakeNodeWithEqLtMethods(0),
                FakeNodeWithEqLtMethods(6),
                FakeNodeWithEqLtMethods(3),
            ),
        )
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="Could not order simple graph cycle as the nodes are of different classes.",
        ):
            sgc.order(raise_on_fail=True)

        sgc = SimpleGraphCycle([FakeNodeWithEqLtMethods(85)], validate=False, ordered=False)
        self.assertFalse(sgc.ordered)
        sgc.order()
        self.assertTrue(sgc.ordered)
        self.assertEqual(sgc.nodes, tuple([FakeNodeWithEqLtMethods(85)]))

        sgc = SimpleGraphCycle(
            [
                FakeNodeWithEqLtMethods(8),
                FakeNodeWithEqLtMethods(2),
                FakeNodeWithEqLtMethods(6),
                FakeNodeWithEqLtMethods(3),
                FakeNodeWithEqLtMethods(4),
                FakeNodeWithEqLtMethods(1),
                FakeNodeWithEqLtMethods(64),
                FakeNodeWithEqLtMethods(32),
            ],
            validate=False,
            ordered=False,
        )
        self.assertFalse(sgc.ordered)
        sgc.order()
        self.assertTrue(sgc.ordered)
        self.assertEqual(
            sgc.nodes,
            tuple(
                [
                    FakeNodeWithEqLtMethods(1),
                    FakeNodeWithEqLtMethods(4),
                    FakeNodeWithEqLtMethods(3),
                    FakeNodeWithEqLtMethods(6),
                    FakeNodeWithEqLtMethods(2),
                    FakeNodeWithEqLtMethods(8),
                    FakeNodeWithEqLtMethods(32),
                    FakeNodeWithEqLtMethods(64),
                ]
            ),
        )

        # Test str method
        self.assertEqual(
            str(sgc),
            "Simple cycle with nodes :\n"
            "FakeNode_1\n"
            "FakeNode_4\n"
            "FakeNode_3\n"
            "FakeNode_6\n"
            "FakeNode_2\n"
            "FakeNode_8\n"
            "FakeNode_32\n"
            "FakeNode_64",
        )

    def test_multigraph_cycle(self):
        mg_cycle1 = MultiGraphCycle([2, 4, 3, 5], [1, 0, 2, 0])
        # Check is_valid method
        is_valid, msg = MultiGraphCycle._is_valid(mg_cycle1)
        self.assertTrue(is_valid)
        self.assertEqual(msg, "")
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="MultiGraphCycle is not valid : "
            "Number of nodes different from number of "
            "edge indices.",
        ):
            MultiGraphCycle([0, 2, 4], [0, 0])  # number of nodes is different from number of edge_indices
        with self.assertRaisesRegex(ValueError, expected_regex="MultiGraphCycle is not valid : Duplicate nodes."):
            MultiGraphCycle([0, 2, 4, 3, 2], [0, 0, 0, 0, 0])  # duplicated nodes
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="MultiGraphCycle is not valid : "
            "Cycles with two nodes cannot use the same "
            "edge for the cycle.",
        ):
            MultiGraphCycle([3, 5], [1, 1])  # number of nodes is different from number of edge_indices

        # Testing equality

        #   Test different cycle lengths inequality
        mg_cycle2 = MultiGraphCycle([2, 3, 4, 5, 6], [1, 0, 2, 0, 0])
        self.assertFalse(mg_cycle1 == mg_cycle2)

        #   Test equality
        mg_cycle2 = MultiGraphCycle([2, 4, 3, 5], [1, 0, 2, 0])
        self.assertTrue(mg_cycle1 == mg_cycle2)
        mg_cycle2 = MultiGraphCycle([4, 3, 5, 2], [0, 2, 0, 1])
        self.assertTrue(mg_cycle1 == mg_cycle2)
        mg_cycle2 = MultiGraphCycle([3, 5, 2, 4], [2, 0, 1, 0])
        self.assertTrue(mg_cycle1 == mg_cycle2)
        mg_cycle2 = MultiGraphCycle([5, 2, 4, 3], [0, 1, 0, 2])
        self.assertTrue(mg_cycle1 == mg_cycle2)
        #   Test equality (reversed)
        mg_cycle2 = MultiGraphCycle([2, 5, 3, 4], [0, 2, 0, 1])
        self.assertTrue(mg_cycle1 == mg_cycle2)
        mg_cycle2 = MultiGraphCycle([5, 3, 4, 2], [2, 0, 1, 0])
        self.assertTrue(mg_cycle1 == mg_cycle2)
        mg_cycle2 = MultiGraphCycle([3, 4, 2, 5], [0, 1, 0, 2])
        self.assertTrue(mg_cycle1 == mg_cycle2)
        mg_cycle2 = MultiGraphCycle([4, 2, 5, 3], [1, 0, 2, 0])
        self.assertTrue(mg_cycle1 == mg_cycle2)
        #   Test inequality
        mg_cycle2 = MultiGraphCycle([2, 5, 3, 4], [0, 1, 0, 1])
        self.assertFalse(mg_cycle1 == mg_cycle2)
        mg_cycle2 = MultiGraphCycle([2, 5, 3, 4], [1, 0, 2, 0])
        self.assertFalse(mg_cycle1 == mg_cycle2)
        mg_cycle2 = MultiGraphCycle([3, 5, 2, 4], [1, 0, 2, 0])
        self.assertFalse(mg_cycle1 == mg_cycle2)

        #   Test Self-loop case
        self.assertTrue(MultiGraphCycle([2], [1]) == MultiGraphCycle([2], [1]))
        self.assertFalse(MultiGraphCycle([1], [1]) == MultiGraphCycle([2], [1]))
        self.assertFalse(MultiGraphCycle([2], [1]) == MultiGraphCycle([2], [0]))
        self.assertFalse(MultiGraphCycle([2], [1]) == MultiGraphCycle([1], [1]))
        self.assertFalse(MultiGraphCycle([2], [0]) == MultiGraphCycle([2], [1]))

        #   Test special case with two nodes
        self.assertTrue(MultiGraphCycle([2, 4], [1, 3]) == MultiGraphCycle([2, 4], [1, 3]))
        self.assertTrue(MultiGraphCycle([2, 4], [1, 3]) == MultiGraphCycle([2, 4], [3, 1]))
        self.assertTrue(MultiGraphCycle([2, 4], [1, 3]) == MultiGraphCycle([4, 2], [3, 1]))
        self.assertTrue(MultiGraphCycle([2, 4], [1, 3]) == MultiGraphCycle([4, 2], [1, 3]))
        self.assertFalse(MultiGraphCycle([2, 4], [1, 3]) == MultiGraphCycle([4, 2], [1, 2]))
        self.assertFalse(MultiGraphCycle([2, 4], [1, 3]) == MultiGraphCycle([4, 0], [1, 3]))

        # Test hashing function
        self.assertEqual(hash(mg_cycle1), 4)
        self.assertEqual(hash(MultiGraphCycle([0], [0])), 1)
        self.assertEqual(hash(MultiGraphCycle([0, 3], [0, 1])), 2)
        self.assertEqual(hash(MultiGraphCycle([0, 3, 5, 7, 8], [0, 1, 0, 0, 0])), 5)

        # Test as_dict from_dict and len method
        mg_cycle = MultiGraphCycle([0, 1, 2, 3], [1, 2, 0, 3])
        self.assertTrue(mg_cycle == MultiGraphCycle.from_dict(mg_cycle.as_dict()))
        self.assertEqual(len(mg_cycle), 4)
        mg_cycle = MultiGraphCycle([2, 5, 3, 4, 1, 0], [0, 0, 2, 0, 1, 0])
        self.assertTrue(mg_cycle == MultiGraphCycle.from_dict(mg_cycle.as_dict()))
        self.assertEqual(len(mg_cycle), 6)
        mg_cycle = MultiGraphCycle([8], [1])
        self.assertTrue(mg_cycle == MultiGraphCycle.from_dict(mg_cycle.as_dict()))
        self.assertEqual(len(mg_cycle), 1)

        # Check the validate method
        # Number of nodes and edges do not match
        mgc = MultiGraphCycle(
            nodes=[
                FakeNodeWithEqMethod(1),
                FakeNodeWithEqMethod(0),
                FakeNodeWithEqMethod(2),
            ],
            edge_indices=[0, 0],
            validate=False,
            ordered=False,
        )
        self.assertFalse(mgc.ordered)
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="MultiGraphCycle is not valid : "
            "Number of nodes different from "
            "number of edge indices.",
        ):
            mgc.validate(check_strict_ordering=False)

        # Empty cycle not valid
        mgc = MultiGraphCycle([], edge_indices=[], validate=False, ordered=False)
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="MultiGraphCycle is not valid : Empty cycle is not valid.",
        ):
            mgc.validate()

        # Multi graph cycle with duplicate nodes not valid
        mgc = MultiGraphCycle(
            [FakeNodeWithEqMethod(1), FakeNodeWithEqMethod(0), FakeNodeWithEqMethod(1)],
            edge_indices=[0, 1, 0],
            validate=False,
            ordered=False,
        )
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="MultiGraphCycle is not valid : Duplicate nodes.",
        ):
            mgc.validate()

        # Multi graph cycle with two nodes cannot use the same edge
        mgc = MultiGraphCycle(
            [FakeNodeWithEqMethod(1), FakeNodeWithEqMethod(0)],
            edge_indices=[1, 1],
            validate=False,
            ordered=False,
        )
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="MultiGraphCycle is not valid : "
            "Cycles with two nodes cannot use the same edge for the cycle.",
        ):
            mgc.validate()

        # Nodes not sortable
        mgc = MultiGraphCycle(
            [FakeNodeWithEqMethod(1), FakeNodeWithEqMethod(0), FakeNodeWithEqMethod(2)],
            edge_indices=[0, 0, 0],
            validate=False,
            ordered=False,
        )
        self.assertFalse(mgc.ordered)
        self.assertEqual(
            mgc.nodes,
            (FakeNodeWithEqMethod(1), FakeNodeWithEqMethod(0), FakeNodeWithEqMethod(2)),
        )
        mgc.validate(check_strict_ordering=False)
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="MultiGraphCycle is not valid : The nodes are not sortable.",
        ):
            mgc.validate(check_strict_ordering=True)

        # Multi graph cycle with nodes that cannot be strictly ordered
        mgc = MultiGraphCycle(
            [
                FakeNodeWithEqMethodWrongSortable(0),
                FakeNodeWithEqMethodWrongSortable(1),
                FakeNodeWithEqMethodWrongSortable(2),
                FakeNodeWithEqMethodWrongSortable(3),
            ],
            edge_indices=[0, 0, 0, 0],
            validate=False,
            ordered=False,
        )
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="MultiGraphCycle is not valid : "
            "The list of nodes in the cycle cannot be strictly ordered.",
        ):
            mgc.validate(check_strict_ordering=True)

        # Check the order method
        mgc = MultiGraphCycle(
            [
                FakeNodeWithEqMethodWrongSortable(0),
                FakeNodeWithEqMethodWrongSortable(1),
                FakeNodeWithEqMethodWrongSortable(2),
                FakeNodeWithEqMethodWrongSortable(3),
            ],
            edge_indices=[0, 0, 0, 0],
            validate=False,
            ordered=False,
        )
        mgc.order(raise_on_fail=False)
        self.assertFalse(mgc.ordered)
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="MultiGraphCycle is not valid : "
            "The list of nodes in the cycle cannot be strictly ordered.",
        ):
            mgc.order(raise_on_fail=True)

        mgc = MultiGraphCycle(
            [
                FakeNodeWithEqMethod(8),
                FakeNodeWithEqMethod(0),
                FakeNodeWithEqMethod(3),
                FakeNodeWithEqMethod(6),
            ],
            edge_indices=[0, 0, 0, 0],
            validate=False,
            ordered=False,
        )
        mgc.order(raise_on_fail=False)
        self.assertFalse(mgc.ordered)
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="MultiGraphCycle is not valid : The nodes are not sortable.",
        ):
            mgc.order(raise_on_fail=True)

        mgc = MultiGraphCycle(
            [
                FakeNodeWithEqLtMethods(8),
                FakeNodeWithEqLtMethods(0),
                FakeNodeWithEqLtMethods(6),
                FakeNodeWithEqLtMethods(3),
            ],
            edge_indices=[2, 5, 3, 7],
            validate=False,
            ordered=False,
        )
        mgc.order(raise_on_fail=True)
        self.assertTrue(mgc.ordered)
        self.assertEqual(
            mgc.nodes,
            (
                FakeNodeWithEqLtMethods(0),
                FakeNodeWithEqLtMethods(6),
                FakeNodeWithEqLtMethods(3),
                FakeNodeWithEqLtMethods(8),
            ),
        )
        self.assertEqual(mgc.edge_indices, (5, 3, 7, 2))

        mgc = MultiGraphCycle(
            [
                FakeNodeWithEqLtMethodsBis(8),
                FakeNodeWithEqLtMethods(0),
                FakeNodeWithEqLtMethods(6),
                FakeNodeWithEqLtMethods(3),
            ],
            edge_indices=[2, 5, 3, 7],
            validate=False,
            ordered=False,
        )
        mgc.order(raise_on_fail=False)
        self.assertFalse(mgc.ordered)
        self.assertEqual(
            mgc.nodes,
            (
                FakeNodeWithEqLtMethodsBis(8),
                FakeNodeWithEqLtMethods(0),
                FakeNodeWithEqLtMethods(6),
                FakeNodeWithEqLtMethods(3),
            ),
        )
        self.assertEqual(mgc.edge_indices, (2, 5, 3, 7))
        with self.assertRaisesRegex(
            ValueError,
            expected_regex="Could not order simple graph cycle as the nodes are of different classes.",
        ):
            mgc.order(raise_on_fail=True)

        mgc = MultiGraphCycle(
            [FakeNodeWithEqLtMethods(85)],
            edge_indices=[7],
            validate=False,
            ordered=False,
        )
        self.assertFalse(mgc.ordered)
        mgc.order()
        self.assertTrue(mgc.ordered)
        self.assertEqual(mgc.nodes, tuple([FakeNodeWithEqLtMethods(85)]))
        self.assertEqual(mgc.edge_indices, tuple([7]))

        mgc = MultiGraphCycle(
            [
                FakeNodeWithEqLtMethods(8),
                FakeNodeWithEqLtMethods(2),
                FakeNodeWithEqLtMethods(6),
                FakeNodeWithEqLtMethods(3),
                FakeNodeWithEqLtMethods(4),
                FakeNodeWithEqLtMethods(1),
                FakeNodeWithEqLtMethods(64),
                FakeNodeWithEqLtMethods(32),
            ],
            edge_indices=[2, 0, 4, 1, 0, 3, 5, 2],
            validate=False,
            ordered=False,
        )
        self.assertFalse(mgc.ordered)
        mgc.order()
        self.assertTrue(mgc.ordered)
        self.assertEqual(
            mgc.nodes,
            tuple(
                [
                    FakeNodeWithEqLtMethods(1),
                    FakeNodeWithEqLtMethods(4),
                    FakeNodeWithEqLtMethods(3),
                    FakeNodeWithEqLtMethods(6),
                    FakeNodeWithEqLtMethods(2),
                    FakeNodeWithEqLtMethods(8),
                    FakeNodeWithEqLtMethods(32),
                    FakeNodeWithEqLtMethods(64),
                ]
            ),
        )
        self.assertEqual(mgc.edge_indices, tuple([0, 1, 4, 0, 2, 2, 5, 3]))

        # Testing all cases for a length-4 cycle
        nodes_ref = tuple(FakeNodeWithEqLtMethods(inode) for inode in [0, 1, 2, 3])
        edges_ref = (3, 6, 9, 12)
        for inodes, iedges in [
            ((0, 1, 2, 3), (3, 6, 9, 12)),
            ((1, 2, 3, 0), (6, 9, 12, 3)),
            ((2, 3, 0, 1), (9, 12, 3, 6)),
            ((3, 0, 1, 2), (12, 3, 6, 9)),
            ((3, 2, 1, 0), (9, 6, 3, 12)),
            ((2, 1, 0, 3), (6, 3, 12, 9)),
            ((1, 0, 3, 2), (3, 12, 9, 6)),
            ((0, 3, 2, 1), (12, 9, 6, 3)),
        ]:
            mgc = MultiGraphCycle(
                [FakeNodeWithEqLtMethods(inode) for inode in inodes],
                edge_indices=[iedge for iedge in iedges],
            )
            strnodes = ", ".join([str(i) for i in inodes])
            self.assertEqual(
                mgc.nodes,
                nodes_ref,
                msg=f"Nodes not equal for inodes = ({', '.join([str(i) for i in inodes])})",
            )
            self.assertEqual(
                mgc.edge_indices,
                edges_ref,
                msg=f"Edges not equal for inodes = ({strnodes})",
            )


class EnvironmentNodesGraphUtilsTest(PymatgenTest):
    def test_cycle(self):
        e1 = EnvironmentNode(central_site="Si", i_central_site=0, ce_symbol="T:4")
        e2 = EnvironmentNode(central_site="Si", i_central_site=3, ce_symbol="T:4")
        e3 = EnvironmentNode(central_site="Si", i_central_site=2, ce_symbol="T:4")
        e4 = EnvironmentNode(central_site="Si", i_central_site=5, ce_symbol="T:4")
        e5 = EnvironmentNode(central_site="Si", i_central_site=1, ce_symbol="T:4")

        # Tests of SimpleGraphCycle with EnvironmentNodes
        c1 = SimpleGraphCycle([e2])
        c2 = SimpleGraphCycle([e2])
        self.assertEqual(c1, c2)

        c1 = SimpleGraphCycle([e1])
        c2 = SimpleGraphCycle([e2])
        self.assertNotEqual(c1, c2)

        c1 = SimpleGraphCycle([e1, e2, e3])
        c2 = SimpleGraphCycle([e2, e1, e3])
        self.assertEqual(c1, c2)
        c2 = SimpleGraphCycle([e2, e3, e1])
        self.assertEqual(c1, c2)

        c1 = SimpleGraphCycle([e3, e2, e4, e1, e5])
        c2 = SimpleGraphCycle([e1, e4, e2, e3, e5])
        self.assertEqual(c1, c2)
        c2 = SimpleGraphCycle([e2, e3, e5, e1, e4])
        self.assertEqual(c1, c2)

        c1 = SimpleGraphCycle([e2, e3, e4, e1, e5])
        c2 = SimpleGraphCycle([e2, e3, e5, e1, e4])
        self.assertNotEqual(c1, c2)

        # Tests of MultiGraphCycle with EnvironmentNodes
        c1 = MultiGraphCycle([e1], [2])
        c2 = MultiGraphCycle([e1], [2])
        self.assertEqual(c1, c2)
        c2 = MultiGraphCycle([e1], [1])
        self.assertNotEqual(c1, c2)
        c2 = MultiGraphCycle([e2], [2])
        self.assertNotEqual(c1, c2)

        c1 = MultiGraphCycle([e1, e2], [0, 1])
        c2 = MultiGraphCycle([e1, e2], [1, 0])
        self.assertEqual(c1, c2)
        c2 = MultiGraphCycle([e2, e1], [1, 0])
        self.assertEqual(c1, c2)
        c2 = MultiGraphCycle([e2, e1], [0, 1])
        self.assertEqual(c1, c2)
        c2 = MultiGraphCycle([e2, e1], [2, 1])
        self.assertNotEqual(c1, c2)

        c1 = MultiGraphCycle([e1, e2, e3], [0, 1, 2])
        c2 = MultiGraphCycle([e2, e1, e3], [0, 2, 1])
        self.assertEqual(c1, c2)
        c2 = MultiGraphCycle([e2, e3, e1], [1, 2, 0])
        self.assertEqual(c1, c2)


if __name__ == "__main__":
    import unittest

    unittest.main()
