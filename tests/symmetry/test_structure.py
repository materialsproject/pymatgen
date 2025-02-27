from __future__ import annotations

from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from pymatgen.core import Lattice, Structure
from pymatgen.util.testing import PymatgenTest


class TestSymmetrizedStructure(PymatgenTest):
    def setUp(self):
        self.structure = Structure(
            lattice=Lattice.cubic(3),
            species=("Fe", "Fe"),
            coords=((0, 0, 0), (0.5, 0.5, 0.5)),
        )

        self.symm_structure = SpacegroupAnalyzer(self.structure).get_symmetrized_structure()

    def test_as_dict(self):
        self.assert_msonable(self.symm_structure)

    def test_serialize(self):
        self.symm_structure.to(fmt="json")
