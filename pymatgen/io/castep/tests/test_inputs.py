from __future__ import annotations

from pymatgen.io.castep.inputs import Cell
from pymatgen.util.testing import PymatgenTest


class CellTest(PymatgenTest):
    def test_cell(self):
        structure = self.get_structure("Si")
        cell = Cell.from_structure(structure)

        # tests a round trip, embedding structure into the .cell and then extracting it again
        assert cell.structure == structure

        # tests string method
        known_contents = """%block positions_frac
Si 0.0 0.0 0.0
Si 0.75 0.5 0.75
%endblock positions_frac"""
        assert known_contents in str(cell)

        # tests adding an unknown element
        self.assertRaises(KeyError, cell.set_block, "unknown keyword", [])

        cell.set_pseudopotential("C19")
        known_contents = """%block species_pot
C19
%endblock species_pot"""
        assert known_contents in str(cell)
