from __future__ import annotations

from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.util.testing import MatSciTest


class TestSymmetrizedStructure(MatSciTest):
    def setup_method(self):
        self.structure = Structure(
            lattice=Lattice.cubic(3),
            species=("Fe", "Fe"),
            coords=((0, 0, 0), (0.5, 0.5, 0.5)),
        )

        self.symm_structure = SpacegroupAnalyzer(self.structure).get_symmetrized_structure()

    def test_str_repr(self):
        assert str(self.symm_structure) == repr(self.symm_structure)
        assert "Reduced Formula: Fe" in str(self.symm_structure)

    def test_dict(self):
        dct = self.symm_structure.as_dict()

        assert isinstance(SymmetrizedStructure.from_dict(dct), SymmetrizedStructure)

    def test_serialize(self):
        assert isinstance(self.symm_structure, SymmetrizedStructure)
        self.assert_msonable(self.symm_structure)

        self.symm_structure.to(fmt="json")

    def test_find_equivalent_sites(self):
        site = self.symm_structure.sites[0]
        assert self.symm_structure.find_equivalent_sites(site) == self.symm_structure.sites
