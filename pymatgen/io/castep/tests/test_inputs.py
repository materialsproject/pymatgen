from pymatgen.io.castep.inputs import Cell
from pymatgen.util.testing import PymatgenTest


class CellTest(PymatgenTest):

    def test_cell(self):

        structure = self.get_structure("Si")
        cell = Cell.from_structure(structure)

        # tests a round trip, embedding structure into the .cell and then extracting it again
        self.assertEqual(cell.structure, structure)

        # tests string method
        known_contents = """%block positions_abs
Si 0.0 0.0 0.0
Si 0.75 0.5 0.75
%endblock positions_abs"""
        self.assertIn(known_contents, str(cell))

        # tests adding an unknown element
        self.assertRaises(KeyError, cell.add_block, "unknown keyword", [])
